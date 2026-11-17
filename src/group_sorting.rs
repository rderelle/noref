//! Selection of the variant groups after extraction: one middle length per group,
//! then the strongest groups that share no k-mer with each other.
use hashbrown::HashSet;

use anyhow::Result;

use crate::dna::{Roller, MAX_K};
use crate::group_extraction::{AnchorPair, BlockArena, RawExtractedGroups, RawGroup};

fn overlap_kmer_size(k: usize) -> usize {
    (2 * k + 1).min(MAX_K)
}

#[derive(Clone, Debug)]
enum IsolateMask {
    Inline(u64),
    Inline512([u64; 8]),
    Dynamic(Vec<u64>),
}

impl IsolateMask {
    fn new(n_isolates: usize) -> Self {
        if n_isolates <= 64 {
            Self::Inline(0)
        } else if n_isolates <= 512 {
            Self::Inline512([0; 8])
        } else {
            Self::Dynamic(vec![0; n_isolates.div_ceil(64)])
        }
    }

    #[inline]
    fn insert(&mut self, isolate: usize) {
        match self {
            Self::Inline(mask) => *mask |= 1u64 << isolate,
            Self::Inline512(words) => words[isolate >> 6] |= 1u64 << (isolate & 63),
            Self::Dynamic(words) => words[isolate >> 6] |= 1u64 << (isolate & 63),
        }
    }

    #[inline]
    fn count(&self) -> usize {
        match self {
            Self::Inline(mask) => mask.count_ones() as usize,
            Self::Inline512(words) => words.iter().map(|w| w.count_ones() as usize).sum(),
            Self::Dynamic(words) => words.iter().map(|w| w.count_ones() as usize).sum(),
        }
    }
}

struct LengthSupport {
    len: usize,
    isolates: IsolateMask,
}

struct LengthFilteredGroup {
    key: AnchorPair,
    group: RawGroup,
    middle_len: usize,
    coverage: usize,
}

enum LengthSupportReject {
    Tie,
    LowSupport,
    ShortBlock,
}

pub(crate) struct SortedGroup {
    /// Anchor pair used for deterministic ordering and early overlap checks.
    pub(crate) key: AnchorPair,
    /// Number of isolates supporting the selected middle length.
    pub(crate) coverage: usize,
    /// Length of the variable middle between the two anchors.
    pub(crate) middle_len: usize,
    /// Deduplicated sequences and ordered occurrences for the selected middle length;
    /// isolates with the anchor pair at several genome positions have no occurrence.
    pub(crate) group: RawGroup,
}

fn select_unique_best_supported_length(
    key: AnchorPair,
    mut group: RawGroup,
    n_isolates: usize,
    min_needed: usize,
) -> Result<LengthFilteredGroup, LengthSupportReject> {
    let mut by_len: Vec<LengthSupport> = Vec::with_capacity(4);
    for occurrence in &group.occurrences {
        let block = group.sequences[occurrence.sequence_id as usize];
        let len = block.len();
        if let Some(bucket) = by_len.iter_mut().find(|bucket| bucket.len == len) {
            bucket.isolates.insert(occurrence.isolate as usize);
        } else {
            let mut isolates = IsolateMask::new(n_isolates);
            isolates.insert(occurrence.isolate as usize);
            by_len.push(LengthSupport { len, isolates });
        }
    }

    let mut best = (0usize, 0usize);
    let mut tie = false;
    for bucket in by_len {
        let c = bucket.isolates.count();
        if c > best.1 {
            best = (bucket.len, c);
            tie = false;
        } else if c == best.1 {
            tie = true;
        }
    }

    if tie {
        return Err(LengthSupportReject::Tie);
    }
    if best.1 < min_needed {
        return Err(LengthSupportReject::LowSupport);
    }
    if best.0 < 1 {
        return Err(LengthSupportReject::ShortBlock);
    }

    // Isolates with the anchor pair at several genome positions (whatever the middle
    // length) may hold paralogs: their occurrences are dropped, so the group is
    // missing data for them. Sequences are kept for the overlap filtering.
    let mut isolates: Vec<u16> = group.occurrences.iter().map(|o| o.isolate).collect();
    isolates.sort_unstable();
    let mut multi_copy: Vec<u16> = isolates
        .windows(2)
        .filter(|w| w[0] == w[1])
        .map(|w| w[0])
        .collect();
    multi_copy.dedup();

    let mut remap = vec![u16::MAX; group.sequences.len()];
    let mut sequences = Vec::with_capacity(group.sequences.len());
    for (old_id, block) in group.sequences.into_iter().enumerate() {
        if block.len() == best.0 {
            let new_id = sequences.len() as u16;
            remap[old_id] = new_id;
            sequences.push(block);
        }
    }
    group.sequences = sequences;
    group.occurrences.retain_mut(|occurrence| {
        let new_id = remap[occurrence.sequence_id as usize];
        if new_id == u16::MAX || multi_copy.binary_search(&occurrence.isolate).is_ok() {
            false
        } else {
            occurrence.sequence_id = new_id;
            true
        }
    });

    Ok(LengthFilteredGroup {
        key,
        group,
        middle_len: best.0,
        coverage: best.1,
    })
}

/// Stream the canonical `overlap_k`-mers of the full path `left anchor + middle +
/// right anchor`; stops early and returns true as soon as `f` returns true.
fn any_path_kmer(
    middle: impl Iterator<Item = u8>,
    key: AnchorPair,
    anchor_k: usize,
    overlap_k: usize,
    mut f: impl FnMut(u64) -> bool,
) -> bool {
    let anchor = |code: u64| {
        (0..anchor_k)
            .rev()
            .map(move |i| ((code >> (2 * i)) & 3) as u8)
    };
    let mut roller = Roller::new(overlap_k);
    anchor(key.0)
        .chain(middle)
        .chain(anchor(key.1))
        .any(|code| roller.push(code) && f(roller.canonical()))
}

fn retain_nonoverlapping_groups(
    mut groups: Vec<SortedGroup>,
    arena: &BlockArena,
    overlap_k: usize,
    anchor_k: usize,
) -> Vec<SortedGroup> {
    groups.sort_unstable_by(|a, b| {
        b.coverage
            .cmp(&a.coverage)
            .then_with(|| b.middle_len.cmp(&a.middle_len))
            .then_with(|| a.key.cmp(&b.key))
    });

    // Greedily keep the strongest groups whose paths share no canonical k-mer with
    // an already retained group.
    let mut used: HashSet<u64> = HashSet::new();
    groups.retain(|candidate| {
        let sequences = &candidate.group.sequences;
        if sequences.iter().any(|&b| {
            any_path_kmer(arena.codes(b), candidate.key, anchor_k, overlap_k, |km| {
                used.contains(&km)
            })
        }) {
            return false;
        }
        for &b in sequences {
            any_path_kmer(arena.codes(b), candidate.key, anchor_k, overlap_k, |km| {
                used.insert(km);
                false
            });
        }
        true
    });
    groups
}

/// One row per isolate: the isolate's middle sequence, or `-` (missing) when it has
/// no occurrence. After length filtering, an isolate has at most one occurrence.
pub(crate) fn rows_by_isolate(
    group: &RawGroup,
    arena: &BlockArena,
    n_isolates: usize,
    len: usize,
) -> Vec<Vec<u8>> {
    let mut rows = vec![vec![b'-'; len]; n_isolates];
    for occurrence in &group.occurrences {
        let block = group.sequences[occurrence.sequence_id as usize];
        let row = &mut rows[occurrence.isolate as usize];
        debug_assert!(
            row.iter().all(|&b| b == b'-'),
            "isolate with two occurrences"
        );
        for (cell, base) in row.iter_mut().zip(arena.bases(block)).take(len) {
            *cell = base;
        }
    }
    rows
}

pub(crate) struct SortedGroups {
    /// Isolate names in deterministic row order.
    pub(crate) isolate_names: Vec<String>,
    /// Blocks (middles) referenced by the occurrences of the groups.
    pub(crate) arena: BlockArena,
    /// Variant groups retained after k-mer non-overlap filtering.
    pub(crate) groups: Vec<SortedGroup>,
    /// Per genome: coordinate of the end of each contig in the concatenated genome.
    pub(crate) contig_ends: Vec<Vec<u32>>,
    /// Number of isolates represented by `isolate_names`.
    pub(crate) n_isolates: usize,
}

pub(crate) fn sort_and_deduplicate_groups(raw: RawExtractedGroups) -> Result<SortedGroups> {
    let RawExtractedGroups {
        isolate_names,
        arena,
        groups: extracted,
        contig_ends,
        k,
        min_needed,
        n_isolates: n,
    } = raw;

    // For each anchor pair, keep the one middle length found in the most isolates (if
    // unique and found in enough isolates).
    let mut groups = Vec::with_capacity(extracted.len());
    for (key, group) in extracted {
        if let Ok(group) = select_unique_best_supported_length(key, group, n, min_needed) {
            groups.push(SortedGroup {
                key: group.key,
                coverage: group.coverage,
                middle_len: group.middle_len,
                group: group.group,
            });
        }
    }

    // Greedily keep the strongest groups that share no k-mer, so that one
    // biological signal cannot contribute the same transient encoded k-mer
    // twice. The detector uses the shorter of the minimum logical path length
    // (`2 * k + 1`) and the largest packable k-mer (32). Canonical k-mers make
    // the comparison independent of the orientation in which groups are stored.
    let nonoverlapping = retain_nonoverlapping_groups(groups, &arena, overlap_kmer_size(k), k);
    Ok(SortedGroups {
        isolate_names,
        arena,
        groups: nonoverlapping,
        contig_ends,
        n_isolates: n,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::{encode, pack, revcomp};
    use crate::group_extraction::BlockOccurrence;

    fn group(arena: &mut BlockArena, rows: &[(u16, &[u8])]) -> RawGroup {
        let mut group = RawGroup::default();
        for &(isolate, seq) in rows {
            let sequence_id = group
                .sequences
                .iter()
                .position(|&block| arena.get(block) == seq)
                .unwrap_or_else(|| {
                    group.sequences.push(arena.push(seq).unwrap());
                    group.sequences.len() - 1
                }) as u16;
            group.occurrences.push(BlockOccurrence {
                isolate,
                sequence_id,
                pos: 0,
            });
        }
        group
    }

    fn candidate(
        arena: &mut BlockArena,
        key: AnchorPair,
        coverage: usize,
        rows: &[&[u8]],
    ) -> SortedGroup {
        let rows: Vec<_> = rows
            .iter()
            .enumerate()
            .map(|(s, &r)| (s as u16, r))
            .collect();
        SortedGroup {
            key,
            coverage,
            middle_len: rows[0].1.len(),
            group: group(arena, &rows),
        }
    }

    fn path_kmers(middle: &[u8], key: AnchorPair, anchor_k: usize, k: usize) -> HashSet<u64> {
        let mut out = HashSet::new();
        any_path_kmer(middle.iter().map(|&b| encode(b)), key, anchor_k, k, |km| {
            out.insert(km);
            false
        });
        out
    }

    #[test]
    fn overlap_k_is_capped_at_32() {
        assert_eq!(overlap_kmer_size(7), 15);
        assert_eq!(overlap_kmer_size(15), 31);
        assert_eq!(overlap_kmer_size(16), 32);
        assert_eq!(overlap_kmer_size(31), 32);
    }

    #[test]
    fn path_reconstruction_includes_anchor_boundaries_and_is_strand_independent() {
        let (left, middle, right) = (b"ACG".as_slice(), b"TTA".as_slice(), b"GCA".as_slice());
        let full = [left, middle, right].concat();
        let canon = |s: &[u8]| pack(s).min(pack(&revcomp(s)));
        let fwd = path_kmers(middle, (pack(left), pack(right)), 3, 4);
        assert_eq!(fwd, full.windows(4).map(canon).collect());
        // Same locus stored in the opposite orientation.
        let rev = path_kmers(
            &revcomp(middle),
            (pack(&revcomp(right)), pack(&revcomp(left))),
            3,
            4,
        );
        assert_eq!(fwd, rev);
    }

    #[test]
    fn nonoverlap_keeps_strongest_and_detects_reverse_complement_overlap() {
        let mut arena = BlockArena::new();
        let shared = b"ACGTTGCAAGGCTTAC";
        let shared_rc = revcomp(shared);
        let unique = b"GGGGGGGGGGGGGGGG";
        let candidates = vec![
            candidate(&mut arena, (50, 60), 2, &[&shared_rc]),
            candidate(&mut arena, (10, 20), 3, &[shared]),
            candidate(&mut arena, (30, 40), 1, &[unique]),
        ];
        let kept: Vec<_> = retain_nonoverlapping_groups(candidates, &arena, 15, 0)
            .into_iter()
            .map(|c| c.key)
            .collect();
        assert_eq!(kept, vec![(10, 20), (30, 40)]);
    }

    #[test]
    fn isolates_with_several_copies_become_missing() {
        let mut arena = BlockArena::new();
        // Isolate 1: two identical copies; isolate 2: two different copies;
        // isolate 3: one copy of the selected length and one of another length.
        let g = group(
            &mut arena,
            &[
                (0, b"ACG"),
                (1, b"ACG"),
                (1, b"ACG"),
                (2, b"ACG"),
                (2, b"ATG"),
                (3, b"ACG"),
                (3, b"ACGT"),
                (4, b"ATG"),
            ],
        );
        let selected = select_unique_best_supported_length((1, 2), g, 6, 2)
            .ok()
            .unwrap();
        // Support (and the selected length) still counts every isolate with the length.
        assert_eq!((selected.middle_len, selected.coverage), (3, 5));
        // Sequences are kept for the overlap filtering, even those of dropped isolates.
        assert_eq!(selected.group.sequences.len(), 2);
        assert_eq!(
            rows_by_isolate(&selected.group, &arena, 6, 3),
            vec![
                b"ACG".to_vec(),
                b"---".to_vec(),
                b"---".to_vec(),
                b"---".to_vec(),
                b"ATG".to_vec(),
                b"---".to_vec()
            ]
        );
    }

    #[test]
    fn length_filter_rejects_ties_low_support_and_short_blocks() {
        let mut arena = BlockArena::new();
        let tie = group(
            &mut arena,
            &[(0, b"ACGTA"), (1, b"ACGTA"), (2, b"ACGTAC"), (3, b"ACGTAC")],
        );
        assert!(matches!(
            select_unique_best_supported_length((1, 2), tie, 4, 2),
            Err(LengthSupportReject::Tie)
        ));
        let low = group(&mut arena, &[(0, b"ACGTA")]);
        assert!(matches!(
            select_unique_best_supported_length((1, 2), low, 4, 2),
            Err(LengthSupportReject::LowSupport)
        ));
        let short = group(&mut arena, &[(0, b""), (1, b"")]);
        assert!(matches!(
            select_unique_best_supported_length((1, 2), short, 2, 2),
            Err(LengthSupportReject::ShortBlock)
        ));
        // Support counts distinct isolates, not copies.
        let copies = group(
            &mut arena,
            &[(0, b"AAAA"), (0, b"CCCC"), (1, b"GGGGG"), (2, b"TTTTT")],
        );
        let selected = select_unique_best_supported_length((1, 2), copies, 3, 2)
            .ok()
            .unwrap();
        assert_eq!((selected.middle_len, selected.coverage), (5, 2));
        assert_eq!(selected.group.occurrences.len(), 2);
    }

    #[test]
    fn isolate_mask_variants_count_distinct_isolates() {
        for n in [64, 65, 512, 513] {
            let mut mask = IsolateMask::new(n);
            for isolate in [0, n - 1, n - 1, 1] {
                mask.insert(isolate);
            }
            assert_eq!(mask.count(), 3);
        }
    }

    #[test]
    fn sort_keeps_original_key() {
        let mut arena = BlockArena::new();
        let groups = vec![(
            (31, 37),
            group(&mut arena, &[(0, b"ACGTAC"), (1, b"ACTTAC")]),
        )];
        let sorted = sort_and_deduplicate_groups(RawExtractedGroups {
            isolate_names: vec!["a".into(), "b".into()],
            arena,
            groups,
            contig_ends: vec![vec![6], vec![6]],
            k: 3,
            min_needed: 2,
            n_isolates: 2,
        })
        .unwrap();
        assert_eq!(sorted.groups.len(), 1);
        assert_eq!(sorted.groups[0].key, (31, 37));
    }
}

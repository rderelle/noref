//! Column selection in each variant group and alignment construction.
use anyhow::Result;
use rayon::{prelude::*, ThreadPoolBuilder};

use crate::dna::{encode, INVALID};
use crate::group_extraction::BlockArena;
use crate::group_sorting::{SortedGroup, SortedGroups};

/// Fully assembled result returned to the CLI layer before writing files.
pub(crate) struct AlignmentResult {
    /// Isolate names in the same row order as `concat`.
    pub(crate) isolate_names: Vec<String>,
    /// Concatenated alignment rows, one byte vector per isolate.
    pub(crate) concat: Vec<Vec<u8>>,
    /// Number of variant groups in the alignment.
    pub(crate) n_groups: usize,
    /// Alignment column just after the last column of each variant group.
    pub(crate) group_ends: Vec<u32>,
    /// Genome position of each variant group in each isolate (`group * n + isolate`),
    /// or [`NO_POSITION`] when the isolate lacks the group.
    pub(crate) positions: Vec<u32>,
    /// Per isolate: end of each contig in the genome (contigs placed end to end).
    pub(crate) contig_ends: Vec<Vec<u32>>,
}

/// Marks isolates without a position for a variant group.
pub(crate) const NO_POSITION: u32 = u32::MAX;

struct FilteredGroup {
    /// Flattened per-isolate rows after retained-column projection.
    rows: Vec<u8>,
    /// Number of retained columns in each flattened row.
    kept_len: usize,
    /// Genome position of the group in each isolate (`NO_POSITION` when missing).
    positions: Vec<u32>,
}

/// Number of rows with each base code in column `col` (slot 4 collects `N`, gaps
/// and any other non-ACGT byte), each row counting `weights[row]` times.
fn base_counts(rows: &[Vec<u8>], weights: &[usize], col: usize) -> [usize; 5] {
    let mut counts = [0; 5];
    for (row, &weight) in rows.iter().zip(weights) {
        counts[encode(row[col]) as usize] += weight;
    }
    counts
}

/// Per-column majority base (ties keep the first of A, C, G, T); `None` if empty.
fn majority_consensus(rows: &[Vec<u8>], weights: &[usize], len: usize) -> Vec<Option<u8>> {
    (0..len)
        .map(|c| {
            let counts = base_counts(rows, weights, c);
            let (best, &count) = counts[..4]
                .iter()
                .enumerate()
                .rev()
                .max_by_key(|&(_, count)| count)?;
            (count > 0).then_some(b"ACGT"[best])
        })
        .collect()
}

/// Mask (with `N`) rows having at least `mask` consecutive differences to the
/// consensus; missing positions break runs.
fn mask_rows_with_long_divergence(rows: &mut [Vec<u8>], consensus: &[Option<u8>], mask: usize) {
    if mask == 0 {
        return;
    }
    for row in rows {
        let mut run = 0usize;
        let divergent = row.iter().zip(consensus).any(|(&b, &cons)| {
            match cons {
                Some(cons) if encode(b) != INVALID && b != cons => run += 1,
                _ => run = 0,
            }
            run >= mask
        });
        if divergent {
            row.fill(b'N');
        }
    }
}

fn is_polymorphic(counts: &[usize; 5]) -> bool {
    counts[..4].iter().filter(|&&count| count > 0).count() > 1
}

fn has_minimum_presence(counts: &[usize; 5], min_needed: usize) -> bool {
    counts[..4].iter().sum::<usize>() >= min_needed
}

/// Whether a column holds three or four bases, the two most frequent of which are
/// equally frequent. Other columns with three or four bases are split into one SNP per
/// base other than the most frequent one (see [`crate::snp_clusters`]); for these, the
/// most frequent base is not defined, so they are left out of the alignment.
fn has_tied_majority(counts: &[usize; 5]) -> bool {
    let mut bases = [counts[0], counts[1], counts[2], counts[3]];
    bases.sort_unstable_by(|a, b| b.cmp(a));
    bases[2] > 0 && bases[0] == bases[1]
}

/// Columns of `rows` (each standing for `weights[row]` isolates) to keep: those present
/// in at least `min_needed` isolates and without a tied majority (see
/// [`has_tied_majority`]), trimmed to the first and last polymorphic ones. Divergent
/// rows are masked first (`mask > 0`). Returns `None` when no retained column is
/// polymorphic.
fn select_columns(
    rows: &mut [Vec<u8>],
    weights: &[usize],
    len: usize,
    min_needed: usize,
    mask: usize,
) -> Option<Vec<usize>> {
    if mask > 0 {
        let consensus = majority_consensus(rows, weights, len);
        mask_rows_with_long_divergence(rows, &consensus, mask);
    }
    let counts: Vec<[usize; 5]> = (0..len).map(|c| base_counts(rows, weights, c)).collect();
    let retained: Vec<usize> = (0..len)
        .filter(|&c| has_minimum_presence(&counts[c], min_needed) && !has_tied_majority(&counts[c]))
        .collect();
    let first = retained.iter().position(|&c| is_polymorphic(&counts[c]))?;
    let last = retained.iter().rposition(|&c| is_polymorphic(&counts[c]))?;
    Some(retained[first..=last].to_vec())
}

fn filter_group(
    raw: SortedGroup,
    arena: &BlockArena,
    n: usize,
    min_needed: usize,
    mask: usize,
) -> Option<FilteredGroup> {
    // The row of an isolate is either missing (gaps, which count in no column
    // statistic) or one of the few sequences of the group: columns are selected on
    // these sequences, each weighted by its number of isolates.
    let (len, group) = (raw.middle_len, &raw.group);
    let mut weights = vec![0; group.sequences.len()];
    for o in &group.occurrences {
        weights[o.sequence_id as usize] += 1;
    }
    let mut rows: Vec<Vec<u8>> = group
        .sequences
        .iter()
        .map(|&block| arena.bases(block).take(len).collect())
        .collect();
    let keep_cols = select_columns(&mut rows, &weights, len, min_needed, mask)?;
    let kept_len = keep_cols.len();
    let mut kept_rows = vec![b'-'; n * kept_len];
    let mut positions = vec![NO_POSITION; n];
    for o in &group.occurrences {
        let isolate = o.isolate as usize;
        let row = &rows[o.sequence_id as usize];
        let kept = &mut kept_rows[isolate * kept_len..(isolate + 1) * kept_len];
        for (cell, &c) in kept.iter_mut().zip(&keep_cols) {
            *cell = row[c];
        }
        positions[isolate] = o.pos;
    }
    Some(FilteredGroup {
        rows: kept_rows,
        kept_len,
        positions,
    })
}

/// Build the alignment from the variant groups, keeping the columns with a base in at
/// least `min_needed` isolates.
pub(crate) fn filter_groups(
    sorted: SortedGroups,
    min_needed: usize,
    mask: usize,
    num_threads: usize,
) -> Result<AlignmentResult> {
    let SortedGroups {
        isolate_names,
        arena,
        groups,
        contig_ends,
        n_isolates: n,
    } = sorted;

    // Per-group filtering is independent; `collect` on an indexed parallel iterator
    // keeps the deterministic order of the sorted groups.
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads.max(1))
        .build()?;
    let filtered: Vec<FilteredGroup> = pool.install(|| {
        groups
            .into_par_iter()
            .map(|raw| filter_group(raw, &arena, n, min_needed, mask))
            .collect::<Vec<_>>()
            .into_iter()
            .flatten()
            .collect()
    });
    drop(arena);

    let n_groups = filtered.len();
    let final_len: usize = filtered.iter().map(|group| group.kept_len).sum();
    let mut concat = vec![Vec::with_capacity(final_len); n];
    let mut group_ends = Vec::with_capacity(n_groups);
    let mut positions = Vec::with_capacity(n_groups * n);
    for group in filtered {
        for (row_out, row) in concat
            .iter_mut()
            .zip(group.rows.chunks_exact(group.kept_len))
        {
            row_out.extend_from_slice(row);
        }
        group_ends.push(concat[0].len() as u32);
        positions.extend_from_slice(&group.positions);
    }
    Ok(AlignmentResult {
        isolate_names,
        concat,
        n_groups,
        group_ends,
        positions,
        contig_ends,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::group_extraction::{BlockOccurrence, RawGroup};

    fn rows(rows: &[&[u8]]) -> Vec<Vec<u8>> {
        rows.iter().map(|r| r.to_vec()).collect()
    }

    /// Base counts of column `col`, each row counting once.
    fn counts(rows: &[Vec<u8>], col: usize) -> [usize; 5] {
        base_counts(rows, &vec![1; rows.len()], col)
    }

    /// Kept columns of `rows`, as the projected rows.
    fn project(input: &[&[u8]], min_needed: usize, mask: usize) -> Option<Vec<Vec<u8>>> {
        let mut r = rows(input);
        let ones = vec![1; r.len()];
        let cols = select_columns(&mut r, &ones, input[0].len(), min_needed, mask)?;
        Some(
            r.iter()
                .map(|row| cols.iter().map(|&c| row[c]).collect())
                .collect(),
        )
    }

    #[test]
    fn exact_minimum_presence_is_retained_without_float_rounding() {
        let mut rows = vec![vec![b'A']; 17];
        rows.extend(vec![vec![b'-']; 3]);
        assert!(has_minimum_presence(&counts(&rows, 0), 17));
        assert!(!has_minimum_presence(&counts(&rows, 0), 18));
        rows[0][0] = b'N';
        assert!(!has_minimum_presence(&counts(&rows, 0), 17));
    }

    #[test]
    fn consensus_ignores_missing_and_breaks_ties_in_acgt_order() {
        let r = rows(&[b"AAN", b"CC-", b"C-G", b"-A-"]);
        assert_eq!(
            majority_consensus(&r, &[1; 4], 3),
            vec![Some(b'C'), Some(b'A'), Some(b'G')]
        );
        // Weighted rows count as many identical rows.
        assert_eq!(
            majority_consensus(&r, &[1, 1, 0, 3], 3),
            vec![Some(b'A'), Some(b'A'), None]
        );
        assert_eq!(majority_consensus(&rows(&[b"N"]), &[1], 1), vec![None]);
    }

    #[test]
    fn masking_requires_consecutive_differences() {
        let consensus = vec![Some(b'A'); 5];
        let mut r = rows(&[b"ACCAA", b"ACCCA", b"CC-CC", b"CCNCC"]);
        mask_rows_with_long_divergence(&mut r, &consensus, 3);
        assert_eq!(r, rows(&[b"ACCAA", b"NNNNN", b"CC-CC", b"CCNCC"]));
        let mut r = rows(&[b"CCCCC"]);
        let gapped = [Some(b'A'), Some(b'A'), None, Some(b'A'), Some(b'A')];
        mask_rows_with_long_divergence(&mut r, &gapped, 3);
        assert_eq!(r[0], b"CCCCC");
        mask_rows_with_long_divergence(&mut r, &consensus, 0);
        assert_eq!(r[0], b"CCCCC");
    }

    #[test]
    fn polymorphism_ignores_gaps_and_n() {
        let r = rows(&[b"AA", b"NC", b"-A"]);
        assert!(!is_polymorphic(&counts(&r, 0)));
        assert!(is_polymorphic(&counts(&r, 1)));
    }

    #[test]
    fn columns_are_trimmed_to_first_and_last_polymorphic_columns() {
        assert_eq!(
            project(&[b"AGCTA", b"ACCAA"], 2, 0),
            Some(rows(&[b"GCT", b"CCA"]))
        );
    }

    #[test]
    fn conserved_groups_are_discarded() {
        assert_eq!(project(&[b"ACCG", b"ACCG", b"ACNG", b"A-CG"], 2, 0), None);
    }

    #[test]
    fn missingness_filtering_happens_before_polymorphic_trimming() {
        let kept = project(&[b"ATCAG", b"AGCTG", b"A-CAG", b"A-CTG"], 3, 0);
        assert_eq!(kept, Some(rows(&[b"A", b"T", b"A", b"T"])));
    }

    #[test]
    fn n_does_not_create_polymorphism_and_counts_as_missing() {
        let kept = project(&[b"ACCAG", b"ACCTG", b"ANCAG", b"ACNTG"], 3, 0);
        assert_eq!(kept, Some(rows(&[b"A", b"T", b"A", b"T"])));
    }

    #[test]
    fn multiallelic_columns_with_a_tied_majority_are_left_out() {
        let r = rows(&[b"AAAAA", b"AACCA", b"CCGGC", b"CGTTN", b"G-A-T"]);
        // Column 0 (A A C C G): A and C tied, three bases: tied.
        // Column 1 (A A C G -): A most frequent: not tied.
        // Column 2 (A C G T A): A most frequent, minor bases tied: not tied.
        // Column 3 (A C G T -): four bases, once each: tied.
        // Column 4 (A A C N T): A most frequent: not tied.
        let tied: Vec<bool> = (0..5).map(|c| has_tied_majority(&counts(&r, c))).collect();
        assert_eq!(tied, vec![true, false, false, true, false]);
        // Two bases, even tied, are a biallelic SNP.
        assert!(!has_tied_majority(&counts(
            &rows(&[b"A", b"A", b"C", b"C"]),
            0
        )));
        let mut kept = r.clone();
        assert_eq!(
            select_columns(&mut kept, &[1; 5], 5, 4, 0),
            Some(vec![1, 2, 4])
        );
        // A group whose only polymorphic column is tied is discarded.
        let mut only = rows(&[b"AA", b"AC", b"AT"]);
        assert_eq!(select_columns(&mut only, &[1; 3], 2, 3, 0), None);
    }

    #[test]
    fn weighted_rows_select_the_columns_of_the_repeated_rows() {
        let distinct: &[&[u8]] = &[b"ACGTACGA", b"ACCTAGGA", b"TCGAACGT", b"ANGTACG-"];
        for weights in [[3, 1, 0, 2], [1, 2, 2, 1], [5, 0, 1, 1]] {
            let repeated: Vec<&[u8]> = distinct
                .iter()
                .zip(weights)
                .flat_map(|(&row, w)| std::iter::repeat_n(row, w))
                .collect();
            for (min_needed, mask) in [(2, 0), (4, 2), (5, 3), (6, 0)] {
                let mut r = rows(distinct);
                let found = select_columns(&mut r, &weights, 8, min_needed, mask);
                let mut expected_rows = rows(&repeated);
                let ones = vec![1; repeated.len()];
                let expected = select_columns(&mut expected_rows, &ones, 8, min_needed, mask);
                assert_eq!(
                    found, expected,
                    "weights {weights:?} min {min_needed} mask {mask}"
                );
            }
        }
    }

    #[test]
    fn divergent_rows_are_masked_before_polymorphism_detection() {
        let input: &[&[u8]] = &[b"ACGTACA", b"ACGTACA", b"AGTACGA"];
        assert_eq!(project(input, 3, 5), None);
        assert_eq!(project(input, 3, 0).unwrap()[0].len(), 5);
    }

    fn sorted(seqs: &[&[u8]]) -> SortedGroups {
        let mut arena = BlockArena::new();
        let mut group = RawGroup::default();
        for (isolate, &seq) in seqs.iter().enumerate() {
            group.sequences.push(arena.push(seq).unwrap());
            group.occurrences.push(BlockOccurrence {
                isolate: isolate as u16 + 1,
                sequence_id: isolate as u16,
                pos: 100 + isolate as u32,
            });
        }
        SortedGroups {
            arena,
            isolate_names: (0..=seqs.len())
                .map(|isolate| format!("s{isolate}"))
                .collect(),
            groups: vec![SortedGroup {
                key: (1, 2),
                coverage: seqs.len(),
                middle_len: seqs[0].len(),
                group,
            }],
            contig_ends: Vec::new(),
            n_isolates: seqs.len() + 1,
        }
    }

    #[test]
    fn filter_builds_alignment() {
        // Isolate 0 is absent (`-` row).
        let r = filter_groups(sorted(&[b"AGCTA", b"ACCAA"]), 2, 0, 1).unwrap();
        assert_eq!(r.n_groups, 1);
        assert_eq!(r.concat, rows(&[b"---", b"GCT", b"CCA"]));
        // Too much missing data when all three isolates are needed.
        let r = filter_groups(sorted(&[b"AGCTA", b"ACCAA"]), 3, 0, 1).unwrap();
        assert_eq!(r.n_groups, 0);
        assert!(r.concat.iter().all(Vec::is_empty));
    }
}

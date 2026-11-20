//! Step 1 of recombination detection: cluster SNPs by isolate profile.
//!
//! SNPs come from alignment columns (already passing `-f`) with two or more bases
//! among the isolates with data (`-` and `N` are missing). A column with exactly two
//! bases is one SNP, whose *profile* is the set of isolates carrying one of them. A
//! column with three or four bases is split into one SNP per base other than the most
//! frequent one: its profile is the set of isolates carrying that base, among all the
//! isolates with data. Each such base is taken as one mutation from the most frequent
//! base, on one branch. (Columns whose two most frequent bases are equally frequent are
//! left out of the alignment beforehand, see [`crate::group_filtering`].) A profile
//! and its complement describe the same split of the isolates, so profiles are stored
//! in a canonical orientation.
//! Two SNPs *agree* when their profiles are identical (or complementary) on the
//! isolates where both have data; missing isolates are wildcards.
//!
//! Clustering is one greedy pass, with no parameter:
//! 1. identical SNPs (same isolates with data, same profile) are collapsed into
//!    patterns, each weighted by its number of SNPs;
//! 2. patterns without missing data come first, then patterns with missing data;
//!    within each group, patterns are taken from the most to the least informative
//!    (frequency of the minor allele among the isolates with data: 50% first, an
//!    allele carried by a single isolate last), then by fewer missing isolates;
//! 3. a pattern that agrees with no existing cluster representative becomes the
//!    representative of a new cluster. A pattern that agrees with `n` representatives
//!    is split between them: each of its SNPs counts `1/n` in each of these clusters.
//!    Without missing data, every distinct profile is therefore its own cluster.
//!    Representatives never change, so the order only decides which profiles become
//!    representatives. No SNP is discarded.
//!
//! Finding the agreeing representatives is exact: every representative is compared
//! with the pattern. Representatives are stored isolate by isolate, 64 per word (see
//! [`Representatives`]), so that a pattern is compared with 64 representatives at a
//! time, until all 64 disagree. Patterns are compared with the representatives in
//! parallel, by batches (see [`cluster_patterns`]).
use hashbrown::HashMap;
use rayon::prelude::*;

use crate::dna::encode;

/// Columns scanned together, so that each alignment row is read in contiguous runs.
const CHUNK: usize = 256;

/// Chunks scanned in parallel before their SNPs are collapsed into patterns, which
/// bounds the memory used by per-SNP bitsets.
const CHUNKS_PER_BATCH: usize = 64;

/// Patterns compared in parallel with the cluster representatives (see
/// [`cluster_patterns`]). Affects speed, not results.
const BATCH: usize = 1024;

/// One bit per isolate.
type Bits = Box<[u64]>;

/// Counts printed after each run.
#[derive(Debug, Default, PartialEq, Eq)]
pub(crate) struct SnpStats {
    /// Distinct (isolates with data, profile) combinations.
    pub(crate) patterns: usize,
    /// Patterns with data in every isolate.
    pub(crate) complete_patterns: usize,
    pub(crate) clusters: usize,
}

/// Result of step 1: the cluster(s) of each SNP, for the next step.
pub(crate) struct SnpClusters {
    /// Alignment column and pattern of each SNP, in column order.
    pub(crate) snps: Vec<(u32, u32)>,
    /// Clusters of each pattern, in ascending order: a SNP split between `n` clusters
    /// counts `1/n` in each.
    pub(crate) pattern_clusters: Vec<Vec<u32>>,
    pub(crate) stats: SnpStats,
}

/// One SNP: a biallelic column, or one base of a column with three or four bases.
struct ColumnSnp {
    column: u32,
    /// Isolates with data; `None` when every isolate has data.
    present: Option<Bits>,
    /// Canonical profile (subset of the isolates with data).
    profile: Bits,
}

/// Identical SNPs collapsed together.
struct Pattern {
    present: Bits,
    profile: Bits,
    /// Number of isolates without data.
    missing: usize,
    /// Number of SNPs.
    weight: u64,
    /// Column of the first SNP (for a deterministic order).
    first_column: u32,
}

/// Bitset with the first `n` bits set.
fn all_isolates(n: usize) -> Bits {
    let words = n.div_ceil(64).max(1);
    let mut full: Bits = vec![u64::MAX; words].into();
    if !n.is_multiple_of(64) {
        full[words - 1] = (1u64 << (n % 64)) - 1;
    }
    full
}

fn is_set(bits: &[u64], i: usize) -> bool {
    bits[i / 64] >> (i % 64) & 1 == 1
}

/// Index of the first set bit.
fn first_set(bits: &[u64]) -> Option<usize> {
    bits.iter()
        .position(|&w| w != 0)
        .map(|w| w * 64 + bits[w].trailing_zeros() as usize)
}

/// Orient `profile` (a subset of `present`) so that it excludes the first isolate of
/// `present`; returns false if the profile is empty or covers all of `present`.
fn canonicalize(profile: &mut [u64], present: &[u64]) -> bool {
    let first = first_set(present).expect("a SNP has isolates with data");
    if is_set(profile, first) {
        for (p, &m) in profile.iter_mut().zip(present) {
            *p = m & !*p;
        }
    }
    profile.iter().any(|&w| w != 0)
}

/// Scan columns `start..end` of the alignment; returns their SNPs, in column order. A
/// column with three or four bases gives one SNP per base other than the most frequent
/// one; if the two most frequent bases are equally frequent (never
/// the case after [`crate::group_filtering`]), it gives none.
fn scan_columns(rows: &[Vec<u8>], start: usize, end: usize) -> Vec<ColumnSnp> {
    let (n, width) = (rows.len(), end - start);
    let words = n.div_ceil(64);
    // One bitset per column and base (A, C, G, T).
    let mut masks = vec![0u64; width * 4 * words];
    for (i, row) in rows.iter().enumerate() {
        let (word, bit) = (i / 64, 1u64 << (i % 64));
        for (j, &b) in row[start..end].iter().enumerate() {
            let code = encode(b) as usize;
            if code < 4 {
                masks[(j * 4 + code) * words + word] |= bit;
            }
        }
    }
    let full = all_isolates(n);
    let count = |m: &[u64]| m.iter().map(|w| w.count_ones()).sum::<u32>();
    let mut snps = Vec::new();
    for j in 0..width {
        let column = &masks[j * 4 * words..(j + 1) * 4 * words];
        let bases: Vec<&[u64]> = column
            .chunks_exact(words)
            .filter(|m| m.iter().any(|&w| w != 0))
            .collect();
        if bases.len() < 2 {
            continue;
        }
        let present: Bits = (0..words)
            .map(|w| bases.iter().fold(0, |acc, m| acc | m[w]))
            .collect();
        let complete = present == full;
        let column = (start + j) as u32;
        let mut push = |base: &[u64]| {
            let mut profile: Bits = base.into();
            canonicalize(&mut profile, &present);
            snps.push(ColumnSnp {
                column,
                present: (!complete).then(|| present.clone()),
                profile,
            });
        };
        if bases.len() == 2 {
            push(bases[0]);
            continue;
        }
        // The most frequent base (the first one, in A, C, G, T order, among equals).
        let sizes: Vec<u32> = bases.iter().map(|m| count(m)).collect();
        let top = (0..bases.len())
            .max_by_key(|&b| (sizes[b], usize::MAX - b))
            .unwrap();
        if sizes.iter().filter(|&&s| s == sizes[top]).count() > 1 {
            continue;
        }
        for (b, base) in bases.iter().enumerate() {
            if b != top {
                push(base);
            }
        }
    }
    snps
}

/// Whether two profiles agree on the isolates where both have data.
fn agree(present_a: &[u64], a: &[u64], present_b: &[u64], b: &[u64]) -> bool {
    let (mut same, mut complement) = (true, true);
    for i in 0..a.len() {
        let shared = present_a[i] & present_b[i];
        let (x, y) = (a[i] & shared, b[i] & shared);
        same &= x == y;
        complement &= x == shared & !y;
    }
    same || complement
}

/// Cluster representatives, stored isolate by isolate in blocks of 64: in block `b`,
/// bit `j` of word `i` is set when representative `64 b + j` has data in isolate `i`
/// and its profile includes it, and bit `j` of word `n + i` when it has data in
/// isolate `i` and its profile excludes it.
struct Representatives {
    n: usize,
    len: usize,
    blocks: Vec<u64>,
}

impl Representatives {
    fn new(n: usize) -> Self {
        Self {
            n,
            len: 0,
            blocks: Vec::new(),
        }
    }

    fn push(&mut self, pattern: &Pattern) {
        if self.len.is_multiple_of(64) {
            self.blocks.resize(self.blocks.len() + 2 * self.n, 0);
        }
        let block = self.blocks.len() - 2 * self.n;
        for i in (0..self.n).filter(|&i| is_set(&pattern.present, i)) {
            let word = if is_set(&pattern.profile, i) {
                i
            } else {
                self.n + i
            };
            self.blocks[block + word] |= 1 << (self.len % 64);
        }
        self.len += 1;
    }

    /// Appends to `hits`, in ascending order, the representatives that agree with a
    /// pattern described by `checks` (see [`checks`]).
    fn agreeing(&self, checks: &[(usize, usize)], hits: &mut Vec<u32>) {
        for (b, block) in self.blocks.chunks_exact(2 * self.n).enumerate() {
            let remaining = self.len - 64 * b;
            let all = if remaining >= 64 {
                u64::MAX
            } else {
                (1 << remaining) - 1
            };
            // Representatives without disagreement so far, in either orientation.
            let (mut same, mut complement) = (all, all);
            for &(s, c) in checks {
                same &= !block[s];
                complement &= !block[c];
                if same | complement == 0 {
                    break;
                }
            }
            let mut agreeing = same | complement;
            while agreeing != 0 {
                hits.push(64 * b as u32 + agreeing.trailing_zeros());
                agreeing &= agreeing - 1;
            }
        }
    }
}

/// For each isolate with data in `pattern`, the words of a block of [`Representatives`]
/// marking the representatives that disagree with it there: in the same orientation,
/// and in the complementary one. Isolates in and out of the profile alternate, as each
/// kind quickly rules out one orientation.
fn checks(pattern: &Pattern, n: usize) -> Vec<(usize, usize)> {
    let with_data = |inside: bool| {
        (0..n)
            .filter(move |&i| is_set(&pattern.present, i) && is_set(&pattern.profile, i) == inside)
    };
    let inside: Vec<(usize, usize)> = with_data(true).map(|i| (n + i, i)).collect();
    let outside: Vec<(usize, usize)> = with_data(false).map(|i| (i, n + i)).collect();
    let mut checks = Vec::with_capacity(inside.len() + outside.len());
    for k in 0..inside.len().max(outside.len()) {
        checks.extend(inside.get(k));
        checks.extend(outside.get(k));
    }
    checks
}

/// Processing order: patterns without missing data first, then by decreasing
/// informativeness (minor-allele frequency among the isolates with data), then by
/// fewer missing isolates, then by first column.
fn processing_order(patterns: &[Pattern]) -> Vec<usize> {
    let count = |b: &[u64]| b.iter().map(|w| w.count_ones() as u64).sum::<u64>();
    // (minor-allele count, isolates with data) of each pattern.
    let freq: Vec<(u64, u64)> = patterns
        .iter()
        .map(|p| {
            let (present, ones) = (count(&p.present), count(&p.profile));
            (ones.min(present - ones), present)
        })
        .collect();
    let mut order: Vec<usize> = (0..patterns.len()).collect();
    order.sort_unstable_by(|&a, &b| {
        let (pa, pb) = (&patterns[a], &patterns[b]);
        let ((ma, na), (mb, nb)) = (freq[a], freq[b]);
        (pa.missing > 0)
            .cmp(&(pb.missing > 0))
            // ma / na > mb / nb first, compared exactly.
            .then((mb * na).cmp(&(ma * nb)))
            .then(pa.missing.cmp(&pb.missing))
            .then(pa.first_column.cmp(&pb.first_column))
    });
    order
}

/// Greedy clustering of patterns, taken in `order`. Returns the number of clusters and
/// the clusters of each pattern. Patterns are taken by batches of `batch`: those of a
/// batch are compared in parallel with the representatives created before it, then in
/// turn with those created by the earlier patterns of the batch.
fn cluster_patterns(
    patterns: &[Pattern],
    n: usize,
    order: &[usize],
    batch: usize,
) -> (usize, Vec<Vec<u32>>) {
    let mut representatives = Representatives::new(n);
    // Pattern of each representative.
    let mut founders: Vec<usize> = Vec::new();
    let mut clusters_of: Vec<Vec<u32>> = vec![Vec::new(); patterns.len()];
    for batch in order.chunks(batch) {
        let before = founders.len();
        let found: Vec<Vec<u32>> = batch
            .par_iter()
            .map(|&p| {
                let mut hits = Vec::new();
                representatives.agreeing(&checks(&patterns[p], n), &mut hits);
                hits
            })
            .collect();
        for (&p, mut hits) in batch.iter().zip(found) {
            let a = &patterns[p];
            hits.extend((before..founders.len()).filter_map(|c| {
                let b = &patterns[founders[c]];
                agree(&a.present, &a.profile, &b.present, &b.profile).then_some(c as u32)
            }));
            if hits.is_empty() {
                hits.push(founders.len() as u32);
                founders.push(p);
                representatives.push(a);
            }
            clusters_of[p] = hits;
        }
    }
    (founders.len(), clusters_of)
}

/// Cluster the SNPs of an alignment (one row per isolate).
pub(crate) fn cluster_snps(rows: &[Vec<u8>], num_threads: usize) -> anyhow::Result<SnpClusters> {
    let n = rows.len();
    let len = rows.first().map_or(0, Vec::len);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads.max(1))
        .build()?;
    let mut stats = SnpStats::default();
    let full = all_isolates(n);

    // Scan columns in parallel batches and collapse identical SNPs into patterns.
    let mut pattern_of_key: HashMap<Bits, u32> = HashMap::new();
    let mut patterns: Vec<Pattern> = Vec::new();
    let mut snps: Vec<(u32, u32)> = Vec::new();
    let n_chunks = len.div_ceil(CHUNK);
    for batch in (0..n_chunks).step_by(CHUNKS_PER_BATCH) {
        let scanned: Vec<Vec<ColumnSnp>> = pool.install(|| {
            (batch..(batch + CHUNKS_PER_BATCH).min(n_chunks))
                .into_par_iter()
                .map(|c| scan_columns(rows, c * CHUNK, ((c + 1) * CHUNK).min(len)))
                .collect()
        });
        for chunk in scanned {
            for snp in chunk {
                let present = snp.present.unwrap_or_else(|| full.clone());
                let key: Bits = present.iter().chain(snp.profile.iter()).copied().collect();
                let id = *pattern_of_key.entry(key).or_insert_with(|| {
                    let missing = n - present
                        .iter()
                        .map(|w| w.count_ones() as usize)
                        .sum::<usize>();
                    patterns.push(Pattern {
                        present,
                        profile: snp.profile,
                        missing,
                        weight: 0,
                        first_column: snp.column,
                    });
                    patterns.len() as u32 - 1
                });
                patterns[id as usize].weight += 1;
                snps.push((snp.column, id));
            }
        }
    }
    drop(pattern_of_key);
    stats.patterns = patterns.len();
    stats.complete_patterns = patterns.iter().filter(|p| p.missing == 0).count();

    let order = processing_order(&patterns);
    let (clusters, pattern_clusters) =
        pool.install(|| cluster_patterns(&patterns, n, &order, BATCH));
    stats.clusters = clusters;
    Ok(SnpClusters {
        snps,
        pattern_clusters,
        stats,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alignment(columns: &[&str]) -> Vec<Vec<u8>> {
        // Columns are given top to bottom (one character per isolate).
        let n = columns[0].len();
        (0..n)
            .map(|i| columns.iter().map(|c| c.as_bytes()[i]).collect())
            .collect()
    }

    /// Column, cluster and weight of each SNP in each of its clusters.
    fn entries(result: &SnpClusters) -> Vec<(u32, u32, f64)> {
        let mut entries = Vec::new();
        for &(column, p) in &result.snps {
            let clusters = &result.pattern_clusters[p as usize];
            for &c in clusters {
                entries.push((column, c, 1.0 / clusters.len() as f64));
            }
        }
        entries
    }

    /// Number of SNPs in each cluster (split SNPs count fractionally).
    fn sizes(result: &SnpClusters) -> Vec<f64> {
        let mut sizes = vec![0.0; result.stats.clusters];
        for (_, cluster, w) in entries(result) {
            sizes[cluster as usize] += w;
        }
        sizes
    }

    fn bits(n: usize, members: &[usize]) -> Bits {
        let mut b = vec![0u64; n.div_ceil(64)].into_boxed_slice();
        for &i in members {
            b[i / 64] |= 1 << (i % 64);
        }
        b
    }

    #[test]
    fn canonical_profiles_exclude_the_first_isolate_with_data() {
        let present = bits(5, &[1, 2, 3, 4]);
        let mut p = bits(5, &[1, 3]);
        assert!(canonicalize(&mut p, &present));
        assert_eq!(p, bits(5, &[2, 4]));
        let mut p = bits(5, &[2, 4]);
        assert!(canonicalize(&mut p, &present));
        assert_eq!(p, bits(5, &[2, 4]));
        let mut all = present.clone();
        assert!(!canonicalize(&mut all, &present));
    }

    #[test]
    fn columns_are_classified_and_counted() {
        let rows = alignment(&[
            "AAAA", // constant
            "AACC", // biallelic, complete
            "CCAA", // same split as above (complement)
            "ACGA", // triallelic
            "A-CN", // biallelic with missing data
            "----", // empty
        ]);
        let snps = scan_columns(&rows, 0, 6);
        let columns: Vec<u32> = snps.iter().map(|s| s.column).collect();
        assert_eq!(columns, vec![1, 2, 3, 3, 4]);
        assert_eq!(snps[0].profile, bits(4, &[2, 3]));
        assert_eq!(snps[1].profile, snps[0].profile);
        assert!(snps[0].present.is_none());
        // Triallelic column: A is the most frequent base; one SNP for C, one for G.
        assert_eq!(snps[2].profile, bits(4, &[1]));
        assert_eq!(snps[3].profile, bits(4, &[2]));
        assert!(snps[2].present.is_none() && snps[3].present.is_none());
        assert_eq!(snps[4].present.as_deref(), Some(&bits(4, &[0, 2])[..]));
        assert_eq!(snps[4].profile, bits(4, &[2]));
    }

    #[test]
    fn multiallelic_columns_give_one_snp_per_minor_base() {
        let rows = alignment(&[
            "AAACGT-", // four bases, A most frequent: SNPs for C, G and T
            "CAATT-T", // T most frequent, isolate 5 missing: SNPs for A and C
            "AACCG--", // A and C tied: no SNP
        ]);
        let snps = scan_columns(&rows, 0, 3);
        let found: Vec<(u32, Bits)> = snps.iter().map(|s| (s.column, s.profile.clone())).collect();
        assert_eq!(
            found,
            vec![
                (0, bits(7, &[3])),
                (0, bits(7, &[4])),
                (0, bits(7, &[5])),
                // A, then C: isolate 0 (the first with data) carries C, so the C
                // profile is stored as its complement.
                (1, bits(7, &[1, 2])),
                (1, bits(7, &[1, 2, 3, 4, 6])),
            ]
        );
        // Every SNP of a column covers all the isolates with data.
        let present = |i: usize| snps[i].present.as_deref().map(<[u64]>::to_vec);
        assert_eq!(present(0), Some(bits(7, &[0, 1, 2, 3, 4, 5]).to_vec()));
        assert_eq!(present(3), Some(bits(7, &[0, 1, 2, 3, 4, 6]).to_vec()));
        assert_eq!(present(3), present(4));
    }

    #[test]
    fn minor_bases_join_the_clusters_of_their_branches() {
        // 8 isolates: branches {4, 5} and {6, 7}, two SNPs each; a triallelic column
        // with G on one branch and T on the other.
        let rows = alignment(&["AAAACCAA", "AAAACCAA", "AAAAAACC", "AAAAAACC", "AAAAGGTT"]);
        let result = cluster_snps(&rows, 1).unwrap();
        // Four biallelic SNPs and two from the triallelic column.
        assert_eq!(result.snps.len(), 6);
        let entries = entries(&result);
        let cluster_of = |column: u32| -> Vec<u32> {
            entries
                .iter()
                .filter(|s| s.0 == column)
                .map(|s| s.1)
                .collect()
        };
        assert_eq!(result.stats.clusters, 2);
        assert_eq!(cluster_of(4), vec![cluster_of(0)[0], cluster_of(2)[0]]);
        // Each piece is one whole SNP of its cluster.
        assert!(entries.iter().all(|e| e.2 == 1.0));
        assert_eq!(sizes(&result), vec![3.0, 3.0]);
    }

    #[test]
    fn agreement_ignores_missing_isolates_and_orientation() {
        let all = bits(5, &[0, 1, 2, 3, 4]);
        let no3 = bits(5, &[0, 1, 2, 4]);
        // {3,4} vs {4} with isolate 3 missing: agree; {2,3,4} vs {4}: disagree.
        assert!(agree(&all, &bits(5, &[3, 4]), &no3, &bits(5, &[4])));
        assert!(!agree(&all, &bits(5, &[2, 3, 4]), &no3, &bits(5, &[4])));
        // Complementary orientation.
        assert!(agree(&all, &bits(5, &[3, 4]), &all, &bits(5, &[0, 1, 2])));
    }

    #[test]
    fn greedy_clustering_splits_ambiguous_snps() {
        // 6 isolates.
        let rows = alignment(&[
            "AAAACC", "AAAACC", "GGGGTT", // {4,5}: one complete pattern of 3 SNPs
            "AAACCC", // {3,4,5}: most informative complete pattern -> cluster 0
            "AAA-CC", // isolate 3 missing: agrees with clusters 0 and 1 -> split
            "AAACC-", // isolate 5 missing: {3,4} agrees with cluster 0 only
            "-CACAA", // isolate 0 missing: agrees with cluster 2 only
            "-CAC-A", // most informative incomplete pattern: no agreement -> cluster 2
        ]);
        let result = cluster_snps(&rows, 2).unwrap();
        assert_eq!(
            result.stats,
            SnpStats {
                patterns: 6,
                complete_patterns: 2,
                clusters: 3,
            }
        );
        assert_eq!(
            entries(&result),
            vec![
                (0, 1, 1.0),
                (1, 1, 1.0),
                (2, 1, 1.0),
                (3, 0, 1.0),
                (4, 0, 0.5),
                (4, 1, 0.5),
                (5, 0, 1.0),
                (6, 2, 1.0),
                (7, 2, 1.0),
            ]
        );
        assert_eq!(sizes(&result), vec![2.5, 3.5, 2.0]);
    }

    #[test]
    fn complete_profiles_come_first_then_informativeness() {
        // The incomplete pattern is the most informative but comes after the complete
        // one; among incomplete patterns the 50% split comes before the 1-isolate one.
        let rows = alignment(&["AAAAAC", "-AACCC", "-AAAAC"]);
        let snps = scan_columns(&rows, 0, 3);
        let n = 6;
        let full = bits(n, &(0..n).collect::<Vec<_>>());
        let patterns: Vec<Pattern> = snps
            .into_iter()
            .map(|snp| {
                let present = snp.present.unwrap_or_else(|| full.clone());
                Pattern {
                    missing: (0..n).filter(|&i| !is_set(&present, i)).count(),
                    present,
                    profile: snp.profile,
                    weight: 1,
                    first_column: snp.column,
                }
            })
            .collect();
        assert_eq!(processing_order(&patterns), vec![0, 1, 2]);
    }

    #[test]
    fn clusters_form_without_any_complete_snp() {
        // Every SNP misses one isolate, as in large or divergent datasets.
        let rows = alignment(&["-AACC", "A-ACC", "AA-CC", "AAA-C", "-CCAA", "C-CAA"]);
        let result = cluster_snps(&rows, 1).unwrap();
        assert_eq!(result.stats.complete_patterns, 0);
        assert_eq!(result.stats.clusters, 1);
        assert_eq!(sizes(&result), vec![6.0]);
        assert_eq!(result.snps.len(), 6);
    }

    fn random_alignment(n: usize, columns: usize, splits: usize, seed: u64) -> Vec<Vec<u8>> {
        let mut state = seed;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            state >> 33
        };
        // A few founder splits with random missing data, so that patterns overlap.
        let splits: Vec<Vec<bool>> = (0..splits)
            .map(|_| (0..n).map(|_| next() % 3 == 0).collect())
            .collect();
        let cols: Vec<String> = (0..columns)
            .map(|_| {
                let split = &splits[next() as usize % splits.len()];
                let missing_rate = [0, 2, 10, 40][next() as usize % 4];
                (0..n)
                    .map(|i| match () {
                        _ if missing_rate > 0 && next() % 100 < missing_rate => '-',
                        _ if split[i] => 'A',
                        _ => 'C',
                    })
                    .collect()
            })
            .collect();
        let refs: Vec<&str> = cols.iter().map(String::as_str).collect();
        alignment(&refs)
    }

    /// Greedy clustering by comparing each pattern with every representative.
    fn cluster_patterns_directly(patterns: &[Pattern], order: &[usize]) -> Vec<Vec<u32>> {
        let mut representatives: Vec<&Pattern> = Vec::new();
        let mut clusters_of = vec![Vec::new(); patterns.len()];
        for &p in order {
            let a = &patterns[p];
            let mut hits: Vec<u32> = (0..representatives.len() as u32)
                .filter(|&c| {
                    let b = representatives[c as usize];
                    agree(&a.present, &a.profile, &b.present, &b.profile)
                })
                .collect();
            if hits.is_empty() {
                hits.push(representatives.len() as u32);
                representatives.push(a);
            }
            clusters_of[p] = hits;
        }
        clusters_of
    }

    #[test]
    fn representatives_agree_as_in_direct_comparison() {
        // Up to 3 words per isolate set, and several blocks of 64 representatives.
        for (n, splits, seed) in [
            (9, 12, 1),
            (40, 12, 2),
            (70, 12, 3),
            (70, 300, 4),
            (150, 300, 5),
        ] {
            let rows = random_alignment(n, 2000, splits, seed);
            let full: Bits = bits(n, &(0..n).collect::<Vec<_>>());
            let mut keyed: HashMap<Bits, usize> = HashMap::new();
            let mut patterns = Vec::new();
            let snps = scan_columns(&rows, 0, rows[0].len());
            for snp in snps {
                let present = snp.present.unwrap_or_else(|| full.clone());
                let key: Bits = present.iter().chain(snp.profile.iter()).copied().collect();
                keyed.entry(key).or_insert_with(|| {
                    patterns.push(Pattern {
                        missing: (0..n).filter(|&i| !is_set(&present, i)).count(),
                        present,
                        profile: snp.profile,
                        weight: 1,
                        first_column: snp.column,
                    });
                    patterns.len() - 1
                });
            }
            let order = processing_order(&patterns);
            let direct = cluster_patterns_directly(&patterns, &order);
            // Whole batches, a single batch, and batches of one pattern.
            for batch in [1, 7, 64, BATCH] {
                let (clusters, clusters_of) = cluster_patterns(&patterns, n, &order, batch);
                if splits > 64 {
                    assert!(clusters > 128, "n={n}: {clusters} clusters");
                }
                assert_eq!(clusters_of, direct, "n={n} batch={batch}");
            }
        }
    }

    #[test]
    fn results_do_not_depend_on_thread_count() {
        let rows = random_alignment(130, 900, 12, 7);
        let one = cluster_snps(&rows, 1).unwrap();
        let four = cluster_snps(&rows, 4).unwrap();
        assert_eq!(one.stats, four.stats);
        assert_eq!(one.snps, four.snps);
        assert_eq!(one.pattern_clusters, four.pattern_clusters);
        // Each SNP's weights sum to 1 (all columns are biallelic here).
        let total: f64 = entries(&one).iter().map(|e| e.2).sum();
        let mut columns: Vec<u32> = one.snps.iter().map(|s| s.0).collect();
        columns.dedup();
        assert!((total - columns.len() as f64).abs() < 1e-6);
    }
}

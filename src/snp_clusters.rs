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
//! Finding the agreeing representatives is exact, by one of two equivalent methods
//! chosen per pattern by cost: a profile with `u` missing isolates stands for `2^u`
//! complete profiles (each missing isolate filled with either base), and two profiles
//! agree exactly when they share one of these completions. Representatives never agree
//! with each other, so their completions are disjoint and can be kept in a hash map:
//! a pattern is placed with `2^u` lookups, or, when `2^u` exceeds the number of
//! clusters, by comparing it with every representative. Both methods return the same
//! clusters, so the choice only affects speed.
use hashbrown::HashMap;
use rayon::prelude::*;

use crate::dna::encode;

/// Columns scanned together, so that each alignment row is read in contiguous runs.
const CHUNK: usize = 256;

/// Chunks scanned in parallel before their SNPs are collapsed into patterns, which
/// bounds the memory used by per-SNP bitsets.
const CHUNKS_PER_BATCH: usize = 64;

/// Founders with more missing isolates than this are compared directly instead of
/// having their `2^u` completions indexed. Affects memory and speed, not results.
const MAX_INDEXED_MISSING: usize = 10;

/// One bit per isolate.
type Bits = Box<[u64]>;

/// Counts printed after each run.
#[derive(Debug, Default, PartialEq, Eq)]
pub(crate) struct SnpStats {
    /// Columns with at least two bases.
    pub(crate) polymorphic: usize,
    /// Distinct (isolates with data, profile) combinations.
    pub(crate) patterns: usize,
    /// Patterns with data in every isolate.
    pub(crate) complete_patterns: usize,
    pub(crate) clusters: usize,
}

/// Result of step 1: the cluster(s) of each SNP, for the next step.
pub(crate) struct SnpClusters {
    /// Alignment column and cluster of each SNP membership, in column order: a SNP
    /// split between several clusters appears once per cluster.
    pub(crate) snps: Vec<(u32, u32)>,
    /// Weight of each entry of `snps`: `1/n` for a SNP split between `n` clusters.
    pub(crate) weights: Vec<f64>,
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
    /// Isolates without data.
    missing: Vec<u32>,
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

/// Scan columns `start..end` of the alignment; returns their SNPs, in column order, and
/// the number of polymorphic columns. A column with three or four bases gives one SNP per base other than
/// the most frequent one; if the two most frequent bases are equally frequent (never
/// the case after [`crate::group_filtering`]), it gives none.
fn scan_columns(rows: &[Vec<u8>], start: usize, end: usize) -> (Vec<ColumnSnp>, usize) {
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
    let mut polymorphic = 0;
    for j in 0..width {
        let column = &masks[j * 4 * words..(j + 1) * 4 * words];
        let bases: Vec<&[u64]> = column
            .chunks_exact(words)
            .filter(|m| m.iter().any(|&w| w != 0))
            .collect();
        if bases.len() < 2 {
            continue;
        }
        polymorphic += 1;
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
    (snps, polymorphic)
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

/// Calls `f` with each canonical completion of `pattern` (each missing isolate set to
/// either base, then oriented so that isolate 0 is excluded).
fn for_each_completion(pattern: &Pattern, full: &[u64], mut f: impl FnMut(&[u64])) {
    let mut key: Vec<u64> = pattern.profile.to_vec();
    let mut canonical = vec![0u64; key.len()];
    let u = pattern.missing.len();
    for gray in 0u64..1 << u {
        if gray > 0 {
            // Gray code: consecutive completions differ by one isolate.
            let i = pattern.missing[gray.trailing_zeros() as usize] as usize;
            key[i / 64] ^= 1 << (i % 64);
        }
        if key[0] & 1 == 1 {
            for ((c, &k), &m) in canonical.iter_mut().zip(&key).zip(full) {
                *c = m & !k;
            }
            f(&canonical);
        } else {
            f(&key);
        }
    }
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
        (!pa.missing.is_empty())
            .cmp(&!pb.missing.is_empty())
            // ma / na > mb / nb first, compared exactly.
            .then((mb * na).cmp(&(ma * nb)))
            .then(pa.missing.len().cmp(&pb.missing.len()))
            .then(pa.first_column.cmp(&pb.first_column))
    });
    order
}

/// Greedy clustering of patterns, taken in `order`. Returns the number of clusters and
/// the clusters of each pattern.
fn cluster_patterns(
    patterns: &[Pattern],
    n: usize,
    order: &[usize],
    use_completions: impl Fn(usize, usize) -> bool,
) -> (usize, Vec<Vec<u32>>) {
    let full = all_isolates(n);
    // Isolates with data and canonical profile of each cluster's representative.
    let mut founders: Vec<(Bits, Bits)> = Vec::new();
    // Completions of indexed representatives. Representatives never agree with each
    // other, so each completion belongs to a single cluster.
    let mut index: HashMap<Bits, u32> = HashMap::new();
    // Representatives whose completions are not indexed.
    let mut unindexed: Vec<u32> = Vec::new();
    let mut clusters_of: Vec<Vec<u32>> = vec![Vec::new(); patterns.len()];
    let mut hits: Vec<u32> = Vec::new();
    for &p in order {
        let pattern = &patterns[p];
        let u = pattern.missing.len();
        let agrees = |c: &u32| {
            let (present, profile) = &founders[*c as usize];
            agree(&pattern.present, &pattern.profile, present, profile)
        };
        hits.clear();
        if use_completions(u, founders.len()) {
            for_each_completion(pattern, &full, |key| {
                if let Some(&c) = index.get(key) {
                    hits.push(c);
                }
            });
            hits.extend(unindexed.iter().filter(|c| agrees(c)));
            hits.sort_unstable();
            hits.dedup();
        } else {
            hits.extend((0..founders.len() as u32).filter(|c| agrees(c)));
        }
        if hits.is_empty() {
            let c = founders.len() as u32;
            founders.push((pattern.present.clone(), pattern.profile.clone()));
            if u <= MAX_INDEXED_MISSING {
                for_each_completion(pattern, &full, |key| {
                    index.insert(key.into(), c);
                });
            } else {
                unindexed.push(c);
            }
            hits.push(c);
        }
        clusters_of[p] = hits.clone();
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
        let scanned: Vec<(Vec<ColumnSnp>, usize)> = pool.install(|| {
            (batch..(batch + CHUNKS_PER_BATCH).min(n_chunks))
                .into_par_iter()
                .map(|c| scan_columns(rows, c * CHUNK, ((c + 1) * CHUNK).min(len)))
                .collect()
        });
        for (chunk, polymorphic) in scanned {
            stats.polymorphic += polymorphic;
            for snp in chunk {
                let present = snp.present.unwrap_or_else(|| full.clone());
                let key: Bits = present.iter().chain(snp.profile.iter()).copied().collect();
                let id = *pattern_of_key.entry(key).or_insert_with(|| {
                    let missing = (0..n as u32)
                        .filter(|&i| !is_set(&present, i as usize))
                        .collect();
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
    stats.complete_patterns = patterns.iter().filter(|p| p.missing.is_empty()).count();

    let order = processing_order(&patterns);
    // Enumerate completions when that takes fewer comparisons than scanning clusters.
    let cheaper = |u: usize, clusters: usize| u < 63 && (1usize << u) <= clusters;
    let (clusters, clusters_of) = cluster_patterns(&patterns, n, &order, cheaper);
    stats.clusters = clusters;
    // One entry per (SNP, cluster), in column order.
    let memberships = snps
        .iter()
        .map(|&(_, p)| clusters_of[p as usize].len())
        .sum();
    let mut entries = Vec::with_capacity(memberships);
    let mut weights = Vec::with_capacity(memberships);
    for &(column, p) in &snps {
        let clusters = &clusters_of[p as usize];
        for &c in clusters {
            entries.push((column, c));
            weights.push(1.0 / clusters.len() as f64);
        }
    }
    Ok(SnpClusters {
        snps: entries,
        weights,
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

    /// Number of SNPs in each cluster (split SNPs count fractionally).
    fn sizes(result: &SnpClusters) -> Vec<f64> {
        let mut sizes = vec![0.0; result.stats.clusters];
        for (&(_, cluster), &w) in result.snps.iter().zip(&result.weights) {
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
        let (snps, polymorphic) = scan_columns(&rows, 0, 6);
        assert_eq!(polymorphic, 4);
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
        let (snps, polymorphic) = scan_columns(&rows, 0, 3);
        assert_eq!(polymorphic, 3);
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
        let cluster_of = |column: u32| -> Vec<u32> {
            result
                .snps
                .iter()
                .filter(|s| s.0 == column)
                .map(|s| s.1)
                .collect()
        };
        assert_eq!(result.stats.clusters, 2);
        assert_eq!(cluster_of(4), vec![cluster_of(0)[0], cluster_of(2)[0]]);
        // Each piece is one whole SNP of its cluster.
        assert_eq!(result.weights, vec![1.0; 6]);
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
                polymorphic: 8,
                patterns: 6,
                complete_patterns: 2,
                clusters: 3,
            }
        );
        let entries: Vec<(u32, u32, f64)> = result
            .snps
            .iter()
            .zip(&result.weights)
            .map(|(&(column, cluster), &w)| (column, cluster, w))
            .collect();
        assert_eq!(
            entries,
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
        let (snps, _) = scan_columns(&rows, 0, 3);
        let n = 6;
        let full = bits(n, &(0..n).collect::<Vec<_>>());
        let patterns: Vec<Pattern> = snps
            .into_iter()
            .map(|snp| {
                let present = snp.present.unwrap_or_else(|| full.clone());
                Pattern {
                    missing: (0..n as u32)
                        .filter(|&i| !is_set(&present, i as usize))
                        .collect(),
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

    fn random_alignment(n: usize, columns: usize, seed: u64) -> Vec<Vec<u8>> {
        let mut state = seed;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            state >> 33
        };
        // A few founder splits with random missing data, so that patterns overlap.
        let splits: Vec<Vec<bool>> = (0..12)
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

    #[test]
    fn completion_lookups_and_direct_comparison_give_identical_clusters() {
        for (n, seed) in [(9, 1), (40, 2), (70, 3)] {
            let rows = random_alignment(n, 400, seed);
            let words = n.div_ceil(64);
            let full: Bits = bits(n, &(0..n).collect::<Vec<_>>());
            let mut keyed: HashMap<Bits, usize> = HashMap::new();
            let mut patterns = Vec::new();
            let (snps, _) = scan_columns(&rows, 0, rows[0].len());
            for snp in snps {
                let present = snp.present.unwrap_or_else(|| full.clone());
                let key: Bits = present.iter().chain(snp.profile.iter()).copied().collect();
                let id = *keyed.entry(key).or_insert_with(|| {
                    patterns.push(Pattern {
                        missing: (0..n as u32)
                            .filter(|&i| !is_set(&present, i as usize))
                            .collect(),
                        present,
                        profile: snp.profile,
                        weight: 0,
                        first_column: snp.column,
                    });
                    patterns.len() - 1
                });
                patterns[id].weight += 1;
            }
            assert_eq!(words, patterns[0].present.len());
            let order = processing_order(&patterns);
            let run =
                |method: fn(usize, usize) -> bool| cluster_patterns(&patterns, n, &order, method);
            // Completions only where enumerable; direct comparison everywhere.
            let lookups = run(|u, _| u <= 16);
            let direct = run(|_, _| false);
            assert_eq!(lookups, direct, "n={n}");
        }
    }

    #[test]
    fn results_do_not_depend_on_thread_count() {
        let rows = random_alignment(130, 900, 7);
        let one = cluster_snps(&rows, 1).unwrap();
        let four = cluster_snps(&rows, 4).unwrap();
        assert_eq!(one.stats, four.stats);
        assert_eq!(one.snps, four.snps);
        assert_eq!(one.weights, four.weights);
        // Each SNP's weights sum to 1 (all columns are biallelic here).
        let total: f64 = one.weights.iter().sum();
        let mut columns: Vec<u32> = one.snps.iter().map(|s| s.0).collect();
        columns.dedup();
        assert!((total - columns.len() as f64).abs() < 1e-6);
    }
}

//! Step 2 of recombination detection: mask SNPs concentrated along each genome.
//!
//! Each isolate is processed independently (one isolate per task). For every SNP
//! cluster, the SNPs of that cluster for which the isolate has a base are placed at
//! their variant group's position in the isolate's genome. A SNP split between `n`
//! clusters (see [`crate::snp_clusters`]) counts `1/n` in each of them. If the
//! cluster's SNPs, of total weight `K`, were scattered uniformly over a genome of
//! length `G`, a window of `w` bases would hold `Binomial(K, w / G)` of them; the
//! binomial upper tail is extended to fractional counts with the regularized
//! incomplete beta function, `P(X >= x) = I_p(x, K - x + 1)`, which is the usual tail
//! for whole counts. Windows start at each SNP; a window holding at least two SNPs
//! whose weights reach the critical weight (upper tail below `alpha`, `-a`, 0.05 by
//! default, Bonferroni-corrected for the `G / w` windows of a genome) is significant, and its SNPs are replaced by
//! `N` in that isolate. Flagged SNPs are then removed from `K`, the critical weight is
//! updated, and the remaining SNPs are scanned again until nothing new is flagged. A
//! SNP split between clusters is masked if it is flagged in any of them.
//!
//! Positions are coordinates in the genome's contigs placed end to end. A window never
//! extends past the end of its contig, since the sequences on either side of a contig
//! boundary are not adjacent in the genome: windows near a contig end are shortened to
//! it and tested against the same critical weight. A shortened window can only hold
//! fewer SNPs than a full one, so this keeps the test conservative.
use hashbrown::HashMap;
use rayon::prelude::*;

use crate::group_filtering::NO_POSITION;

/// Default family-wise significance level of the scan, per isolate and cluster (`-a`).
pub(crate) const DEFAULT_ALPHA: f64 = 0.05;

/// A window is significant when its weight is at least the critical weight minus this.
const TOLERANCE: f64 = 1e-9;

/// Counts printed after each run.
#[derive(Debug, Default, PartialEq, Eq)]
pub(crate) struct RecombinationStats {
    /// SNP cells (isolate × SNP column) with a base, before masking.
    pub(crate) cells: usize,
    /// SNP cells replaced by `N`.
    pub(crate) masked_cells: usize,
    /// SNP columns masked in at least one isolate, in ascending order.
    pub(crate) masked_columns: Vec<u32>,
    /// Isolates with at least one masked SNP.
    pub(crate) isolates: usize,
    /// Only when requested (see [`mask_recombinant_snps`]): per isolate, the masked
    /// columns in ascending order with their base before masking.
    pub(crate) original_bases: Vec<Vec<(u32, u8)>>,
}

/// Inputs describing where each SNP lies.
pub(crate) struct SnpLayout<'a> {
    /// Alignment column and cluster of each SNP, in column order (a SNP split
    /// between several clusters appears once per cluster).
    pub(crate) snps: &'a [(u32, u32)],
    /// Weight of each entry of `snps` (1, or a fraction for split SNPs).
    pub(crate) weights: &'a [f64],
    pub(crate) n_clusters: usize,
    /// Alignment column just after the last column of each variant group.
    pub(crate) group_ends: &'a [u32],
    /// Genome position of each variant group in each isolate (`group * n + isolate`).
    pub(crate) positions: &'a [u32],
    /// Per isolate: end of each contig in the genome (contigs placed end to end).
    pub(crate) contig_ends: &'a [Vec<u32>],
}

/// Smallest count `c` such that P(Binomial(k, p) >= c) < `threshold`, or `None` when
/// no count up to `k` is significant (or `k < 2`).
#[cfg(test)]
fn critical_count(k: u32, p: f64, threshold: f64) -> Option<u32> {
    if k < 2 || p >= 1.0 {
        return None;
    }
    let (kf, q) = (k as f64, 1.0 - p);
    // Start just above the mean (and at 2 at least), where the tail is still large.
    let mut c = ((kf * p).floor() as u32 + 1).max(2);
    if c > k {
        return None;
    }
    // ln P(X = c) = ln C(k, c) + c ln p + (k - c) ln q.
    let ln_choose: f64 = (1..=c).map(|i| ((k - c + i) as f64 / i as f64).ln()).sum();
    let mut pmf_c = (ln_choose + c as f64 * p.ln() + (k - c) as f64 * q.ln()).exp();
    // Upper tail from c: terms decrease geometrically above the mean.
    let mut tail = 0.0;
    let mut pmf = pmf_c;
    for j in c..=k {
        tail += pmf;
        if pmf < tail * 1e-16 {
            break;
        }
        pmf *= (k - j) as f64 / (j + 1) as f64 * p / q;
    }
    while tail >= threshold {
        tail -= pmf_c;
        pmf_c *= (k - c) as f64 / (c + 1) as f64 * p / q;
        c += 1;
        if c > k {
            return None;
        }
    }
    Some(c)
}

/// ln Γ(x) for x > 0 (Lanczos approximation, g = 7).
fn ln_gamma(x: f64) -> f64 {
    const G: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    let pi = std::f64::consts::PI;
    if x < 0.5 {
        return (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let t = x + 7.5;
    let a = G[0] + (1..9).map(|i| G[i] / (x + i as f64)).sum::<f64>();
    0.5 * (2.0 * pi).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// Continued fraction for the incomplete beta function (modified Lentz).
fn beta_cf(a: f64, b: f64, x: f64) -> f64 {
    const TINY: f64 = 1e-300;
    let (qab, qap, qam) = (a + b, a + 1.0, a - 1.0);
    let clamp = |v: f64| if v.abs() < TINY { TINY } else { v };
    let mut c = 1.0;
    let mut d = 1.0 / clamp(1.0 - qab * x / qap);
    let mut h = d;
    for m in 1..10_000 {
        let (m, m2) = (m as f64, 2.0 * m as f64);
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 / clamp(1.0 + aa * d);
        c = clamp(1.0 + aa / c);
        h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 / clamp(1.0 + aa * d);
        c = clamp(1.0 + aa / c);
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-15 {
            break;
        }
    }
    h
}

/// Regularized incomplete beta function I_x(a, b).
fn beta_inc(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front =
        (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * beta_cf(a, b, x) / a
    } else {
        1.0 - front * beta_cf(b, a, 1.0 - x) / b
    }
}

/// Smallest weight `x` such that the binomial
/// upper tail, extended to real counts as P(X >= x) = I_p(x, k - x + 1), is below
/// `threshold`; `None` when even `x = k` is not significant. For whole `k` and whole
/// counts this gives the same decisions as `critical_count`.
fn critical_weight(k: f64, p: f64, threshold: f64) -> Option<f64> {
    if k <= 0.0 || p >= 1.0 {
        return None;
    }
    let tail = |x: f64| beta_inc(x, k - x + 1.0, p);
    if tail(k) >= threshold {
        return None;
    }
    let (mut lo, mut hi) = (0.0f64, k);
    while hi - lo > 1e-10 {
        let mid = 0.5 * (lo + hi);
        if tail(mid) < threshold {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some(hi)
}

/// Critical weight for SNPs of total weight `k`, or `None` when no window can be
/// significant: either no weight is, or `heaviest`, the weight of the heaviest window,
/// is too low.
///
/// The bisection of [`critical_weight`] takes about forty evaluations of the tail, and
/// most clusters of an isolate have no significant window. As the tail decreases when
/// the weight increases, these cases are settled with one evaluation, slightly above
/// `heaviest`: if even that weight is not significant, the critical weight is higher
/// still and no window reaches it. `QUICK_EXIT_MARGIN` is far larger than the precision
/// of the bisection and than [`TOLERANCE`], so the windows found significant are the
/// same as when the critical weight is always computed. Critical weights are kept in
/// `known`.
fn critical_weight_if_reached(
    k: f64,
    heaviest: f64,
    p: f64,
    threshold: f64,
    known: &mut HashMap<u64, Option<f64>>,
) -> Option<f64> {
    const QUICK_EXIT_MARGIN: f64 = 1e-6;
    if let Some(&critical) = known.get(&k.to_bits()) {
        return critical;
    }
    let x = heaviest + QUICK_EXIT_MARGIN;
    if beta_inc(x, k - x + 1.0, p) >= threshold {
        return None;
    }
    let critical = critical_weight(k, p, threshold);
    known.insert(k.to_bits(), critical);
    critical
}

/// Flag SNPs in significant windows, iterating until no new SNP is flagged.
/// `snps` holds (position, column, weight) sorted by position; flagged columns are
/// appended to `masked`. Windows are shortened to the end of their contig
/// (`contig_ends`). A window is significant when it holds at least two SNPs whose
/// summed weight reaches the critical weight for the cluster's total weight `K`.
fn scan_cluster(
    snps: &mut Vec<(u32, u32, f64)>,
    contig_ends: &[u32],
    window: u32,
    p: f64,
    threshold: f64,
    critical: &mut HashMap<u64, Option<f64>>,
    masked: &mut Vec<u32>,
) {
    let mut flagged = Vec::new();
    let mut prefix = Vec::new();
    // Windows holding at least two SNPs: index of their first SNP, and of their last
    // SNP plus one.
    let mut spans: Vec<(usize, usize)> = Vec::new();
    loop {
        prefix.clear();
        prefix.push(0.0f64);
        for &(_, _, w) in snps.iter() {
            prefix.push(prefix.last().unwrap() + w);
        }
        let k = *prefix.last().unwrap();
        spans.clear();
        let mut end = 0;
        let mut contig = 0;
        for start in 0..snps.len() {
            // Positions are sorted, so the contig index only moves forward.
            while contig < contig_ends.len() && contig_ends[contig] <= snps[start].0 {
                contig += 1;
            }
            let Some(&contig_end) = contig_ends.get(contig) else {
                break;
            };
            // The window ends at p + w or at the contig end, whichever comes first;
            // this limit never decreases, so `end` only moves forward.
            let limit = (snps[start].0 as u64 + window as u64).min(contig_end as u64);
            end = end.max(start);
            while end < snps.len() && (snps[end].0 as u64) < limit {
                end += 1;
            }
            if end - start >= 2 {
                spans.push((start, end));
            }
        }
        if spans.is_empty() {
            return;
        }
        let weight = |&(start, end): &(usize, usize)| prefix[end] - prefix[start];
        let heaviest = spans.iter().map(weight).fold(0.0, f64::max);
        let Some(c) = critical_weight_if_reached(k, heaviest, p, threshold, critical) else {
            return;
        };
        flagged.clear();
        flagged.resize(snps.len(), false);
        let mut new = 0;
        for span in &spans {
            if weight(span) >= c - TOLERANCE {
                for f in &mut flagged[span.0..span.1] {
                    new += !*f as usize;
                    *f = true;
                }
            }
        }
        if new == 0 {
            return;
        }
        let mut i = 0;
        snps.retain(|&(_, column, _)| {
            let keep = !flagged[i];
            if !keep {
                masked.push(column);
            }
            i += 1;
            keep
        });
    }
}

/// Mask, in each isolate's alignment row, the SNPs found concentrated along its genome.
/// `alpha` is the significance level per isolate and cluster, before Bonferroni
/// correction for the windows of a genome. With `keep_original_bases`, the bases
/// replaced by `N` are kept in [`RecombinationStats::original_bases`].
pub(crate) fn mask_recombinant_snps(
    rows: &mut [Vec<u8>],
    layout: &SnpLayout,
    window: u32,
    alpha: f64,
    num_threads: usize,
    keep_original_bases: bool,
) -> anyhow::Result<RecombinationStats> {
    let n = rows.len();
    // SNPs ordered by cluster (then column), with each SNP's variant group.
    let mut group = 0usize;
    let mut snp_group = Vec::with_capacity(layout.snps.len());
    for &(column, _) in layout.snps {
        while layout.group_ends[group] <= column {
            group += 1;
        }
        snp_group.push(group as u32);
    }
    let mut by_cluster: Vec<u32> = (0..layout.snps.len() as u32).collect();
    by_cluster.sort_by_key(|&s| layout.snps[s as usize].1);
    // Column, variant group and weight of the SNPs of each cluster, stored in the
    // order in which every isolate reads them.
    let by_cluster: Vec<(u32, u32, f64)> = by_cluster
        .iter()
        .map(|&s| {
            let s = s as usize;
            (layout.snps[s].0, snp_group[s], layout.weights[s])
        })
        .collect();
    let mut cluster_bounds = vec![0usize; layout.n_clusters + 1];
    for &(_, cluster) in layout.snps {
        cluster_bounds[cluster as usize + 1] += 1;
    }
    for c in 0..layout.n_clusters {
        cluster_bounds[c + 1] += cluster_bounds[c];
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads.max(1))
        .build()?;
    type IsolateResult = (usize, Vec<u32>, Vec<(u32, u8)>);
    let per_isolate: Vec<IsolateResult> = pool.install(|| {
        rows.par_iter_mut()
            .enumerate()
            .map(|(isolate, row)| {
                let contig_ends = &layout.contig_ends[isolate];
                let genome = contig_ends.last().copied().unwrap_or(0);
                // Genome position of each variant group in this isolate.
                let positions: Vec<u32> = layout
                    .positions
                    .iter()
                    .skip(isolate)
                    .step_by(n)
                    .copied()
                    .collect();
                let mut cells = 0;
                let mut masked = Vec::new();
                let p = window as f64 / genome.max(1) as f64;
                let threshold = alpha / (genome as f64 / window as f64).max(1.0);
                // Critical weight for each total weight K (keyed by its bits).
                let mut critical: HashMap<u64, Option<f64>> = HashMap::new();
                // SNP cells with a base and a position (each column counted once).
                let mut last = u32::MAX;
                for (s, &(column, _)) in layout.snps.iter().enumerate() {
                    if column != last {
                        last = column;
                        let g = snp_group[s] as usize;
                        cells += (matches!(row[column as usize], b'A' | b'C' | b'G' | b'T')
                            && positions[g] != NO_POSITION)
                            as usize;
                    }
                }
                let mut snps: Vec<(u32, u32, f64)> = Vec::new();
                for c in 0..layout.n_clusters {
                    snps.clear();
                    for &(column, g, weight) in
                        &by_cluster[cluster_bounds[c]..cluster_bounds[c + 1]]
                    {
                        if matches!(row[column as usize], b'A' | b'C' | b'G' | b'T') {
                            let pos = positions[g as usize];
                            if pos != NO_POSITION {
                                snps.push((pos, column, weight));
                            }
                        }
                    }
                    if genome == 0 || snps.len() < 2 {
                        continue;
                    }
                    snps.sort_unstable_by_key(|s| (s.0, s.1));
                    scan_cluster(
                        &mut snps,
                        contig_ends,
                        window,
                        p,
                        threshold,
                        &mut critical,
                        &mut masked,
                    );
                }
                // A SNP split between clusters is masked if flagged in any of them.
                masked.sort_unstable();
                masked.dedup();
                let original = if keep_original_bases {
                    masked.iter().map(|&c| (c, row[c as usize])).collect()
                } else {
                    Vec::new()
                };
                for &column in &masked {
                    row[column as usize] = b'N';
                }
                (cells, masked, original)
            })
            .collect()
    });

    let mut stats = RecombinationStats::default();
    let mut column_masked = vec![false; rows.first().map_or(0, Vec::len)];
    for (cells, masked, original) in per_isolate {
        if keep_original_bases {
            stats.original_bases.push(original);
        }
        stats.cells += cells;
        stats.masked_cells += masked.len();
        stats.isolates += !masked.is_empty() as usize;
        for column in masked {
            column_masked[column as usize] = true;
        }
    }
    stats.masked_columns = (0..column_masked.len() as u32)
        .filter(|&c| column_masked[c as usize])
        .collect();
    Ok(stats)
}

/// Re-check the given columns after masking: those that no longer have a base in at
/// least `min_needed` isolates or no longer have two different bases fail, and are
/// returned in the order they were given. They are removed with [`remove_columns`].
pub(crate) fn failing_columns(rows: &[Vec<u8>], columns: &[u32], min_needed: usize) -> Vec<u32> {
    let mut failing = Vec::new();
    for &column in columns {
        let mut seen = [false; 4];
        let mut present = 0;
        for row in rows.iter() {
            let code = crate::dna::encode(row[column as usize]);
            if code < 4 {
                seen[code as usize] = true;
                present += 1;
            }
        }
        let polymorphic = seen.iter().filter(|&&s| s).count() >= 2;
        if present < min_needed || !polymorphic {
            failing.push(column);
        }
    }
    failing
}

/// Remove the given columns from every row.
pub(crate) fn remove_columns(rows: &mut [Vec<u8>], columns: &[u32]) {
    if columns.is_empty() {
        return;
    }
    let mut drop = vec![false; rows.first().map_or(0, Vec::len)];
    for &column in columns {
        drop[column as usize] = true;
    }
    rows.par_iter_mut().for_each(|row| {
        let mut i = 0;
        row.retain(|_| {
            i += 1;
            !drop[i - 1]
        });
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Exact binomial upper tail, for comparison.
    fn tail(k: u32, p: f64, c: u32) -> f64 {
        (c..=k)
            .map(|j| {
                let ln_choose: f64 = (1..=j).map(|i| ((k - j + i) as f64 / i as f64).ln()).sum();
                (ln_choose + j as f64 * p.ln() + (k - j) as f64 * (1.0 - p).ln()).exp()
            })
            .sum()
    }

    #[test]
    fn critical_count_matches_exact_binomial_tail() {
        for &k in &[2u32, 3, 5, 20, 100, 1000, 5000] {
            for &p in &[0.001, 0.009, 0.05, 0.3] {
                for &threshold in &[0.05, 1e-3, 1e-6] {
                    match critical_count(k, p, threshold) {
                        Some(c) => {
                            assert!(tail(k, p, c) < threshold * (1.0 + 1e-9), "k={k} p={p}");
                            if c > ((k as f64 * p).floor() as u32 + 1).max(2) {
                                assert!(
                                    tail(k, p, c - 1) >= threshold * (1.0 - 1e-9),
                                    "k={k} p={p}"
                                );
                            }
                        }
                        None => assert!(k < 2 || tail(k, p, k) >= threshold * (1.0 - 1e-9)),
                    }
                }
            }
        }
        assert_eq!(critical_count(1, 0.01, 0.05), None);
        assert_eq!(critical_count(10, 1.0, 0.05), None);
    }

    #[test]
    fn concentrated_snps_are_flagged_and_scattered_ones_kept() {
        // 1 Mb genome, 20 kb windows: 40 SNPs evenly spread plus 8 within 5 kb.
        let (window, genome) = (20_000u32, 1_000_000u32);
        let mut snps: Vec<(u32, u32, f64)> =
            (0..40).map(|i| (i * 25_000 + 1_000, i, 1.0)).collect();
        snps.extend((0..8).map(|i| (500_000 + i * 600, 100 + i, 1.0)));
        snps.sort_unstable_by_key(|s| (s.0, s.1));
        let p = window as f64 / genome as f64;
        let threshold = DEFAULT_ALPHA / (genome as f64 / window as f64);
        let mut masked = Vec::new();
        scan_cluster(
            &mut snps,
            &[genome],
            window,
            p,
            threshold,
            &mut HashMap::new(),
            &mut masked,
        );
        masked.sort_unstable();
        // The 8 concentrated SNPs (and possibly the evenly spaced one at 501 kb).
        assert!((100..108).all(|c| masked.contains(&c)), "{masked:?}");
        assert!(masked.len() <= 9);
        assert_eq!(snps.len() + masked.len(), 48);
    }

    #[test]
    fn rows_are_masked_per_isolate() {
        // 2 isolates, 1 Mb genomes. Cluster 0: 6 SNPs in one group region (positions
        // close together in isolate 0 only), plus 20 scattered SNPs in cluster 1.
        let n = 2;
        let mut snps = Vec::new();
        let mut positions = Vec::new();
        let mut group_ends = Vec::new();
        for g in 0..26u32 {
            snps.push((g, if g < 6 { 0 } else { 1 }));
            group_ends.push(g + 1);
            let scattered = 10_000 + g * 37_000;
            // Isolate 0: cluster 0 concentrated at ~200 kb; isolate 1: scattered.
            positions.push(if g < 6 {
                200_000 + g * 1_000
            } else {
                scattered
            });
            positions.push(scattered);
        }
        let mut rows = vec![vec![b'A'; 26], vec![b'C'; 26]];
        let weights = vec![1.0; snps.len()];
        let layout = SnpLayout {
            snps: &snps,
            weights: &weights,
            n_clusters: 2,
            group_ends: &group_ends,
            positions: &positions,
            contig_ends: &[vec![1_000_000], vec![1_000_000]],
        };
        let stats =
            mask_recombinant_snps(&mut rows, &layout, 20_000, DEFAULT_ALPHA, 2, false).unwrap();
        assert_eq!(&rows[0][..6], b"NNNNNN");
        assert!(rows[0][6..].iter().all(|&b| b == b'A'));
        assert!(rows[1].iter().all(|&b| b == b'C'));
        assert_eq!(stats.cells, 52);
        assert_eq!(
            (
                stats.masked_cells,
                stats.masked_columns.len(),
                stats.isolates
            ),
            (6, 6, 1)
        );
        assert_eq!(n, rows.len());
        assert!(stats.original_bases.is_empty());

        // Same masking when the original bases are kept, which are then reported.
        let mut kept = vec![vec![b'A'; 26], vec![b'C'; 26]];
        kept[0][3] = b'G';
        let stats =
            mask_recombinant_snps(&mut kept, &layout, 20_000, DEFAULT_ALPHA, 1, true).unwrap();
        assert_eq!(&kept[0][..6], b"NNNNNN");
        let expected: Vec<(u32, u8)> = (0..6)
            .map(|c| (c, if c == 3 { b'G' } else { b'A' }))
            .collect();
        assert_eq!(stats.original_bases, vec![expected, Vec::new()]);

        // A stricter level: the six SNPs within 5 kb (tail 0.02^6, about 6e-11) are
        // no longer significant once the level is divided by the 50 windows.
        let mut strict = vec![vec![b'A'; 26], vec![b'C'; 26]];
        let stats = mask_recombinant_snps(&mut strict, &layout, 20_000, 1e-12, 1, false).unwrap();
        assert_eq!(stats.masked_cells, 0);
        assert!(strict[0].iter().all(|&b| b == b'A'));
    }

    #[test]
    fn windows_are_shortened_at_contig_boundaries() {
        let (window, genome) = (20_000u32, 1_000_000u32);
        let p = window as f64 / genome as f64;
        let threshold = DEFAULT_ALPHA / (genome as f64 / window as f64);
        let run = |cluster: &[(u32, u32)], ends: &[u32]| {
            let mut snps: Vec<(u32, u32, f64)> =
                cluster.iter().map(|&(p, c)| (p, c, 1.0)).collect();
            let mut masked = Vec::new();
            scan_cluster(
                &mut snps,
                ends,
                window,
                p,
                threshold,
                &mut HashMap::new(),
                &mut masked,
            );
            masked.len()
        };
        // K = 4 SNPs: 3 are needed in a window. Two at 498-499 kb, two at 501-502 kb.
        let split: Vec<(u32, u32)> = [498_000, 499_000, 501_000, 502_000]
            .iter()
            .enumerate()
            .map(|(i, &pos)| (pos, i as u32))
            .collect();
        assert_eq!(critical_count(4, p, threshold), Some(3));
        // One contig: all four fall in one window.
        assert_eq!(run(&split, &[genome]), 4);
        // A boundary at 500 kb separates them: the window starting at 498 kb is
        // shortened to 500 kb (2 SNPs), so nothing is significant.
        assert_eq!(run(&split, &[500_000, genome]), 0);
        // Four SNPs 1 kb before the end of a 3 kb contig: the shortened window still
        // holds them, so short contigs and contig ends are tested.
        let near_end: Vec<(u32, u32)> = (0..4).map(|i| (2_000 + i * 200, i)).collect();
        assert_eq!(run(&near_end, &[3_000, genome]), 4);
    }

    #[test]
    fn columns_failing_after_masking_are_removed() {
        // 4 isolates (rows), 5 columns; min_needed = 3.
        let mut rows = vec![
            b"AANAC".to_vec(),
            b"CANAC".to_vec(),
            b"ANCNG".to_vec(),
            b"CNNAG".to_vec(),
        ];
        // Column 0 (A C A C): still polymorphic, kept.
        // Column 1 (A A N N): two isolates with data and one base: fails both checks.
        // Column 2 (N N C N): one isolate with data and one base: fails both checks.
        // Column 3 (A A N A): three isolates with data, one base: not polymorphic only.
        // Column 4: not checked, kept.
        let removed = failing_columns(&rows, &[0, 1, 2, 3], 3);
        assert_eq!(removed, vec![1, 2, 3]);
        // Finding failing columns leaves the rows unchanged.
        assert_eq!(rows[0], b"AANAC");
        remove_columns(&mut rows, &removed);
        assert_eq!(
            rows,
            vec![
                b"AC".to_vec(),
                b"CC".to_vec(),
                b"AG".to_vec(),
                b"CG".to_vec()
            ]
        );
        // Nothing to remove: rows are unchanged.
        let removed = failing_columns(&rows, &[0, 1], 3);
        assert!(removed.is_empty());
        remove_columns(&mut rows, &removed);
        assert_eq!(rows[0], b"AC");
    }

    #[test]
    fn missing_bases_are_ignored() {
        let snps = vec![(0u32, 0u32), (1, 0), (2, 0)];
        let positions = vec![100u32, 200, 300];
        let mut rows = vec![b"N-A".to_vec()];
        let weights = vec![1.0; snps.len()];
        let layout = SnpLayout {
            snps: &snps,
            weights: &weights,
            n_clusters: 1,
            group_ends: &[1, 2, 3],
            positions: &positions,
            contig_ends: &[vec![10_000_000]],
        };
        let stats =
            mask_recombinant_snps(&mut rows, &layout, 20_000, DEFAULT_ALPHA, 1, false).unwrap();
        assert_eq!(stats.cells, 1);
        assert_eq!(stats.masked_cells, 0);
        assert_eq!(rows[0], b"N-A");
    }

    #[test]
    fn critical_weight_matches_critical_count_for_whole_counts() {
        for &k in &[2u32, 3, 5, 20, 100, 1000, 5000] {
            for &p in &[0.0001, 0.001, 0.009, 0.05, 0.3] {
                for &threshold in &[0.05, 1e-3, 1e-6] {
                    let whole = critical_weight(k as f64, p, threshold)
                        .map(|x| (x - 1e-9).ceil().max(2.0) as u32)
                        .filter(|&c| c <= k);
                    assert_eq!(
                        whole,
                        critical_count(k, p, threshold),
                        "k={k} p={p} t={threshold}"
                    );
                }
            }
        }
    }

    #[test]
    fn quick_exit_agrees_with_the_bisection() {
        // Whatever the weight of the heaviest window, including weights right at the
        // critical weight, windows are significant in the same cases as when the
        // critical weight is always computed.
        for &k in &[2.0, 2.5, 3.0, 7.25, 20.0, 100.0, 1000.5] {
            for &p in &[0.0001, 0.001, 0.009, 0.05, 0.3, 1.0] {
                for &threshold in &[0.05, 1e-3, 1e-6] {
                    let full = critical_weight(k, p, threshold);
                    let mut weights: Vec<f64> = (0..=40).map(|i| k * i as f64 / 40.0).collect();
                    if let Some(c) = full {
                        for d in [
                            0.0, 1e-12, 1e-10, 9e-10, 1e-9, 1.1e-9, 3e-9, 1e-7, 1e-6, 2e-6, 1e-4,
                        ] {
                            weights.extend([c - d, c + d]);
                        }
                    }
                    for w in weights {
                        let quick =
                            critical_weight_if_reached(k, w, p, threshold, &mut HashMap::new());
                        let significant = |c: Option<f64>| c.is_some_and(|c| w >= c - TOLERANCE);
                        assert_eq!(significant(quick), significant(full), "k={k} p={p} w={w}");
                        assert!(quick.is_none() || quick == full);
                    }
                }
            }
        }
        // A known critical weight is returned as is.
        let mut known = HashMap::from([(3.0f64.to_bits(), Some(1.5))]);
        assert_eq!(
            critical_weight_if_reached(3.0, 0.0, 0.01, 0.05, &mut known),
            Some(1.5)
        );
    }

    #[test]
    fn split_snps_count_as_fractions() {
        // 1 Mb genome, 10 kb windows. Two whole SNPs plus two halves close together:
        // weight 3 in one window.
        let (window, genome) = (10_000u32, 1_000_000u32);
        let p = window as f64 / genome as f64;
        let threshold = DEFAULT_ALPHA / (genome as f64 / window as f64);
        let x = critical_weight(3.0, p, threshold).unwrap();
        assert!(x > 1.0 && x <= 2.0, "{x}");
        let mut snps = vec![
            (1_000, 0, 1.0),
            (2_000, 1, 0.5),
            (3_000, 2, 0.5),
            (900_000, 3, 1.0),
        ];
        let mut masked = Vec::new();
        scan_cluster(
            &mut snps,
            &[genome],
            window,
            p,
            threshold,
            &mut HashMap::new(),
            &mut masked,
        );
        masked.sort_unstable();
        assert_eq!(masked, vec![0, 1, 2]);
        // A single SNP is never significant on its own, whatever its weight.
        let mut lone = vec![(1_000, 0, 1.0), (900_000, 1, 1.0)];
        let mut masked = Vec::new();
        scan_cluster(
            &mut lone,
            &[genome],
            window,
            p,
            threshold,
            &mut HashMap::new(),
            &mut masked,
        );
        assert!(masked.is_empty());
    }
}

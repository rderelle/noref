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
//!
//! The SNPs of an isolate are sorted by position once; the SNPs of a cluster are then
//! sorted by their rank in that order. Three shortcuts give the same results as
//! scanning every window of every cluster with the exact critical weight: most clusters
//! have no significant window, which a bound of their window weights shows before their
//! SNPs are sorted (see [`WindowBound`]); the critical weight is only computed as
//! precisely as the comparisons need (see [`CriticalWeight`]); and the last round,
//! which flags nothing, is mostly skipped (see [`scan_cluster`]).
use hashbrown::HashMap;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

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
    /// Alignment column and pattern of each SNP, in column order.
    pub(crate) snps: &'a [(u32, u32)],
    /// Clusters of each pattern: a SNP split between `n` clusters counts `1/n` in each.
    pub(crate) pattern_clusters: &'a [Vec<u32>],
    pub(crate) n_clusters: usize,
    /// Alignment column just after the last column of each variant group.
    pub(crate) group_ends: &'a [u32],
    /// Genome position of each variant group in each isolate (`group * n + isolate`).
    pub(crate) positions: &'a [u32],
    /// Per isolate: end of each contig in the genome (contigs placed end to end).
    pub(crate) contig_ends: &'a [Vec<u32>],
}

/// A SNP of a cluster, placed on the genome of an isolate.
#[derive(Clone, Copy, Debug)]
struct Snp {
    pos: u32,
    /// End (excluded) of the window starting at this SNP: `pos + w`, or the end of its
    /// contig if that comes first.
    limit: u32,
    column: u32,
    weight: f64,
}

/// SNPs given as (position, column, weight), sorted by position, placed on a genome
/// whose contigs end at `contig_ends`, with windows of `window` bases.
fn place(
    snps: impl Iterator<Item = (u32, u32, f64)>,
    contig_ends: &[u32],
    window: u32,
) -> Vec<Snp> {
    let mut contig = 0;
    snps.map(|(pos, column, weight)| {
        // Positions are sorted, so the contig index only moves forward.
        while contig < contig_ends.len() && contig_ends[contig] <= pos {
            contig += 1;
        }
        // A position past the last contig starts no window.
        let limit = contig_ends.get(contig).map_or(pos, |&end| {
            (pos as u64 + window as u64).min(end as u64) as u32
        });
        Snp {
            pos,
            limit,
            column,
            weight,
        }
    })
    .collect()
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

/// Critical weight for SNPs of total weight `k`: the smallest weight `x` such that the
/// binomial upper tail, extended to real counts as P(X >= x) = I_p(x, k - x + 1), is
/// below `threshold`. It is found by bisection, to within 1e-10, but the bisection is
/// only carried as far as needed to compare window weights with it: whatever the
/// remaining steps, a weight of at least `hi` (minus [`TOLERANCE`]) reaches the
/// critical weight and a weight below `lo` (minus [`TOLERANCE`]) does not. Windows are
/// therefore found significant exactly as with the critical weight of the complete
/// bisection, in a fraction of its steps.
struct CriticalWeight {
    k: f64,
    p: f64,
    threshold: f64,
    lo: f64,
    hi: f64,
}

impl CriticalWeight {
    /// `None` when even the weight `k` is not significant.
    fn new(k: f64, p: f64, threshold: f64) -> Option<Self> {
        let critical = Self {
            k,
            p,
            threshold,
            lo: 0.0,
            hi: k,
        };
        (k > 0.0 && p < 1.0 && critical.tail(k) < threshold).then_some(critical)
    }

    fn tail(&self, x: f64) -> f64 {
        beta_inc(x, self.k - x + 1.0, self.p)
    }

    /// One step of the bisection; false once it has converged (to `hi`).
    fn narrow(&mut self) -> bool {
        if self.hi - self.lo <= 1e-10 {
            return false;
        }
        let mid = 0.5 * (self.lo + self.hi);
        if self.tail(mid) < self.threshold {
            self.hi = mid;
        } else {
            self.lo = mid;
        }
        true
    }

    /// Whether a window of this weight is significant.
    fn reached_by(&mut self, weight: f64) -> bool {
        loop {
            if weight >= self.hi - TOLERANCE {
                return true;
            }
            if weight < self.lo - TOLERANCE || !self.narrow() {
                return false;
            }
        }
    }
}

/// The critical weight of the complete bisection; `None` when even `x = k` is not
/// significant. For whole `k` and whole counts this gives the same decisions as
/// `critical_count`.
#[cfg(test)]
fn critical_weight(k: f64, p: f64, threshold: f64) -> Option<f64> {
    let mut critical = CriticalWeight::new(k, p, threshold)?;
    while critical.narrow() {}
    Some(critical.hi)
}

/// An upper bound of window weights obtained otherwise than these weights (with sums
/// in another order, or `f32` weights), raised far above its rounding errors, and far
/// above [`TOLERANCE`]. If it is not significant, no window is.
fn raised(bound: f64) -> f64 {
    bound * (1.0 + 1e-6) + 1e-6
}

/// Bound of the weight of the windows of a cluster, computed without sorting its SNPs:
/// a window of `w` bases lies within three consecutive bins of `w / 2` bases (rounded
/// up), the first of which holds the SNP it starts at.
struct WindowBound {
    width: u32,
    bins: Vec<f64>,
    /// Bins holding SNPs: the first `n_filled`, and room for one more.
    filled: Vec<u32>,
    n_filled: usize,
}

impl WindowBound {
    fn new(genome: u32, window: u32) -> Self {
        let width = window.div_ceil(2);
        let n_bins = (genome / width) as usize + 1;
        Self {
            width,
            bins: vec![0.0; n_bins],
            filled: vec![0; n_bins + 1],
            n_filled: 0,
        }
    }

    /// Bin of a position (a SNP past the end of the genome is in no window: any bin
    /// will do).
    fn bin(&self, pos: u32) -> u32 {
        (pos / self.width).min(self.bins.len() as u32 - 1)
    }

    fn add(&mut self, bin: u32, weight: f32) {
        // Without a branch: the bin is kept in `filled` only if it was empty.
        self.filled[self.n_filled] = bin;
        self.n_filled += (self.bins[bin as usize] == 0.0) as usize;
        self.bins[bin as usize] += weight as f64;
    }

    /// Whether a window of the SNPs added since the last call may be significant (if
    /// not, none is). The bins are emptied.
    fn may_be_significant(&mut self, p: f64, threshold: f64) -> bool {
        let filled = &self.filled[..self.n_filled];
        let bins = &self.bins;
        let heaviest = filled
            .iter()
            .map(|&b| {
                bins[b as usize..bins.len().min(b as usize + 3)]
                    .iter()
                    .sum()
            })
            .fold(0.0, f64::max);
        let mut k = 0.0;
        for &b in filled {
            k += self.bins[b as usize];
            self.bins[b as usize] = 0.0;
        }
        self.n_filled = 0;
        // The tail increases with `k`: lowering `k` covers the `f32` weights.
        let (x, k) = (raised(heaviest), k * (1.0 - 1e-6));
        x >= k || beta_inc(x, k - x + 1.0, p) < threshold
    }
}

/// Calls `f(start, end)` for each window holding at least two SNPs: the window
/// starting at SNP `start` holds the SNPs up to `end` (excluded). `snps` are sorted by
/// position.
fn for_each_window(snps: &[Snp], mut f: impl FnMut(usize, usize)) {
    let mut end = 0;
    for (start, snp) in snps.iter().enumerate() {
        // Window limits never decrease, so `end` only moves forward.
        end = end.max(start);
        while end < snps.len() && snps[end].pos < snp.limit {
            end += 1;
        }
        if end - start >= 2 {
            f(start, end);
        }
    }
}

/// Buffers reused by [`scan_cluster`].
#[derive(Default)]
struct Scratch {
    prefix: Vec<f64>,
    flagged: Vec<bool>,
}

/// Flag SNPs in significant windows, round after round until no new SNP is flagged.
/// `snps` are sorted by position; flagged columns are appended to `masked`. A window
/// is significant when it holds at least two SNPs whose summed weight reaches the
/// critical weight for the cluster's total weight `K`, kept for each `K` in `known`.
fn scan_cluster(
    snps: &mut Vec<Snp>,
    p: f64,
    threshold: f64,
    known: &mut HashMap<u64, Option<CriticalWeight>>,
    scratch: &mut Scratch,
    masked: &mut Vec<u32>,
) {
    let Scratch { prefix, flagged } = scratch;
    // Heaviest window left unflagged by the previous round. Removing flagged SNPs only
    // lowers the weight of the remaining windows: a round flags nothing new unless
    // this weight reaches the new critical weight.
    let mut unflagged = f64::INFINITY;
    loop {
        prefix.clear();
        prefix.push(0.0f64);
        for snp in snps.iter() {
            prefix.push(prefix[prefix.len() - 1] + snp.weight);
        }
        let k = prefix[snps.len()];
        let weight = |start: usize, end: usize| prefix[end] - prefix[start];
        let critical = known
            .entry(k.to_bits())
            .or_insert_with(|| CriticalWeight::new(k, p, threshold));
        let Some(critical) = critical else {
            return;
        };
        if !critical.reached_by(raised(unflagged)) {
            return;
        }
        flagged.clear();
        flagged.resize(snps.len(), false);
        // Window ends never decrease: SNPs before `flagged_end` are already flagged.
        let mut flagged_end = 0;
        unflagged = 0.0;
        for_each_window(snps, |start, end| {
            let w = weight(start, end);
            if critical.reached_by(w) {
                flagged[start.max(flagged_end)..end].fill(true);
                flagged_end = end;
            } else {
                unflagged = unflagged.max(w);
            }
        });
        if flagged_end == 0 {
            return;
        }
        let mut i = 0;
        snps.retain(|snp| {
            let keep = !flagged[i];
            if !keep {
                masked.push(snp.column);
            }
            i += 1;
            keep
        });
    }
}

/// Flag the SNPs concentrated along the genome of an isolate, whose alignment row is
/// `row`. `snp_groups` holds the variant group and weight of each SNP, and
/// `cluster_snps[c]` the SNPs of cluster `c`. Returns the SNP cells with a base and a
/// position (each column counted once) and the flagged columns, in ascending order.
fn flag_isolate(
    isolate: usize,
    row: &[u8],
    layout: &SnpLayout,
    snp_groups: &[(u32, f64)],
    cluster_snps: &[&[u32]],
    window: u32,
    alpha: f64,
) -> (usize, Vec<u32>) {
    // Genome position of each variant group in this isolate.
    let positions: Vec<u32> = layout
        .positions
        .iter()
        .skip(isolate)
        .step_by(layout.contig_ends.len())
        .copied()
        .collect();
    let contig_ends = &layout.contig_ends[isolate];
    let genome = contig_ends.last().copied().unwrap_or(0);
    let p = window as f64 / genome.max(1) as f64;
    let threshold = alpha / (genome as f64 / window as f64).max(1.0);
    // SNPs with a base and a position, sorted by position then column (as
    // `(position << 32) | SNP`).
    let mut sorted: Vec<u64> = Vec::new();
    let (mut cells, mut last) = (0, u32::MAX);
    for (s, &(column, _)) in layout.snps.iter().enumerate() {
        let pos = positions[snp_groups[s].0 as usize];
        let with_data =
            matches!(row[column as usize], b'A' | b'C' | b'G' | b'T') && pos != NO_POSITION;
        if with_data {
            sorted.push((pos as u64) << 32 | s as u64);
        }
        if column != last {
            last = column;
            cells += with_data as usize;
        }
    }
    if genome == 0 {
        return (cells, Vec::new());
    }
    sorted.sort_unstable();
    // Rank in genome order of each SNP (the SNPs of a cluster are sorted by sorting
    // their ranks), and its bin and weight for `bound`; `u32::MAX` without data.
    let mut bound = WindowBound::new(genome, window);
    let mut rank = vec![u32::MAX; layout.snps.len()];
    let mut binned = vec![(u32::MAX, 0f32); layout.snps.len()];
    for (r, &key) in sorted.iter().enumerate() {
        let s = key as u32 as usize;
        rank[s] = r as u32;
        binned[s] = (bound.bin((key >> 32) as u32), snp_groups[s].1 as f32);
    }
    let in_genome_order = place(
        sorted.iter().map(|&key| {
            let s = key as u32 as usize;
            ((key >> 32) as u32, layout.snps[s].0, snp_groups[s].1)
        }),
        contig_ends,
        window,
    );
    drop(sorted);

    let mut masked = Vec::new();
    // Critical weight for each total weight K (keyed by its bits).
    let mut known: HashMap<u64, Option<CriticalWeight>> = HashMap::new();
    let mut ranks: Vec<u32> = Vec::new();
    let mut snps: Vec<Snp> = Vec::new();
    let mut scratch = Scratch::default();
    for cluster in cluster_snps {
        for &s in cluster.iter() {
            let (b, w) = binned[s as usize];
            if b != u32::MAX {
                bound.add(b, w);
            }
        }
        // Most clusters have no significant window, which the bound shows without
        // sorting their SNPs.
        if !bound.may_be_significant(p, threshold) {
            continue;
        }
        ranks.clear();
        ranks.extend(
            cluster
                .iter()
                .map(|&s| rank[s as usize])
                .filter(|&r| r != u32::MAX),
        );
        if ranks.len() < 2 {
            continue;
        }
        ranks.sort_unstable();
        snps.clear();
        snps.extend(ranks.iter().map(|&r| in_genome_order[r as usize]));
        scan_cluster(
            &mut snps,
            p,
            threshold,
            &mut known,
            &mut scratch,
            &mut masked,
        );
    }
    // A SNP split between clusters is masked if flagged in any of them.
    masked.sort_unstable();
    masked.dedup();
    (cells, masked)
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
    // Variant group and weight of each SNP.
    let mut group = 0usize;
    let snp_groups: Vec<(u32, f64)> = layout
        .snps
        .iter()
        .map(|&(column, pattern)| {
            while layout.group_ends[group] <= column {
                group += 1;
            }
            let clusters = layout.pattern_clusters[pattern as usize].len();
            (group as u32, 1.0 / clusters as f64)
        })
        .collect();
    // SNPs of each cluster, in column order.
    let clusters_of = |&(_, pattern): &(u32, u32)| &layout.pattern_clusters[pattern as usize];
    let mut bounds = vec![0usize; layout.n_clusters + 1];
    for &c in layout.snps.iter().flat_map(clusters_of) {
        bounds[c as usize + 1] += 1;
    }
    for c in 0..layout.n_clusters {
        bounds[c + 1] += bounds[c];
    }
    let mut all_cluster_snps = vec![0u32; bounds[layout.n_clusters]];
    let mut next = bounds.clone();
    for (s, snp) in layout.snps.iter().enumerate() {
        for &c in clusters_of(snp) {
            all_cluster_snps[next[c as usize]] = s as u32;
            next[c as usize] += 1;
        }
    }
    let cluster_snps: Vec<&[u32]> = bounds
        .windows(2)
        .map(|b| &all_cluster_snps[b[0]..b[1]])
        .collect();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads.max(1))
        .build()?;
    // Columns masked in at least one isolate, one bit per column. Each isolate's
    // masked columns are only kept with `keep_original_bases`.
    let n_columns = rows.first().map_or(0, Vec::len);
    let column_masked: Vec<AtomicU64> = (0..n_columns.div_ceil(64))
        .map(|_| AtomicU64::new(0))
        .collect();
    type IsolateResult = (usize, usize, Vec<(u32, u8)>);
    let per_isolate: Vec<IsolateResult> = pool.install(|| {
        rows.par_iter_mut()
            .enumerate()
            .map(|(isolate, row)| {
                let (cells, masked) = flag_isolate(
                    isolate,
                    row,
                    layout,
                    &snp_groups,
                    &cluster_snps,
                    window,
                    alpha,
                );
                let original = if keep_original_bases {
                    masked.iter().map(|&c| (c, row[c as usize])).collect()
                } else {
                    Vec::new()
                };
                for &column in &masked {
                    row[column as usize] = b'N';
                    column_masked[column as usize / 64].fetch_or(1 << (column % 64), Relaxed);
                }
                (cells, masked.len(), original)
            })
            .collect()
    });

    let mut stats = RecombinationStats::default();
    for (cells, masked, original) in per_isolate {
        if keep_original_bases {
            stats.original_bases.push(original);
        }
        stats.cells += cells;
        stats.masked_cells += masked;
        stats.isolates += (masked > 0) as usize;
    }
    let column_masked: Vec<u64> = column_masked
        .into_iter()
        .map(AtomicU64::into_inner)
        .collect();
    stats.masked_columns = (0..n_columns as u32)
        .filter(|&c| column_masked[c as usize / 64] >> (c % 64) & 1 == 1)
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

    /// Scan one cluster given as (position, column, weight) sorted by position; returns
    /// the flagged columns, sorted, and the number of SNPs left.
    fn scan(
        snps: &[(u32, u32, f64)],
        contig_ends: &[u32],
        window: u32,
        p: f64,
        threshold: f64,
    ) -> (Vec<u32>, usize) {
        let mut snps = place(snps.iter().copied(), contig_ends, window);
        let mut masked = Vec::new();
        let mut known = HashMap::new();
        let mut scratch = Scratch::default();
        scan_cluster(
            &mut snps,
            p,
            threshold,
            &mut known,
            &mut scratch,
            &mut masked,
        );
        masked.sort_unstable();
        (masked, snps.len())
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
        let (masked, left) = scan(&snps, &[genome], window, p, threshold);
        // The 8 concentrated SNPs (and possibly the evenly spaced one at 501 kb).
        assert!((100..108).all(|c| masked.contains(&c)), "{masked:?}");
        assert!(masked.len() <= 9);
        assert_eq!(left + masked.len(), 48);
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
        // Pattern 0 is cluster 0, pattern 1 cluster 1.
        let layout = SnpLayout {
            snps: &snps,
            pattern_clusters: &[vec![0], vec![1]],
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
            let snps: Vec<(u32, u32, f64)> = cluster.iter().map(|&(p, c)| (p, c, 1.0)).collect();
            scan(&snps, ends, window, p, threshold).0.len()
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
        let layout = SnpLayout {
            snps: &snps,
            pattern_clusters: &[vec![0]],
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
    fn lazy_bisection_agrees_with_the_complete_bisection() {
        // Windows are significant in the same cases as with the critical weight of the
        // complete bisection, including weights right at it, in any order.
        for &k in &[2.0, 2.5, 3.0, 7.25, 20.0, 100.0, 1000.5] {
            for &p in &[0.0001, 0.001, 0.009, 0.05, 0.3, 1.0] {
                for &threshold in &[0.05, 1e-3, 1e-6] {
                    let full = critical_weight(k, p, threshold);
                    let Some(mut lazy) = CriticalWeight::new(k, p, threshold) else {
                        assert_eq!(full, None);
                        continue;
                    };
                    let c = full.unwrap();
                    let mut weights: Vec<f64> = (0..=40).map(|i| k * i as f64 / 40.0).collect();
                    for d in [
                        0.0, 1e-12, 1e-10, 9e-10, 1e-9, 1.1e-9, 3e-9, 1e-7, 1e-6, 2e-6, 1e-4,
                    ] {
                        weights.extend([c - d, c + d]);
                    }
                    for i in 0..weights.len() {
                        // Alternately from either end.
                        let w = weights[if i % 2 == 0 {
                            i / 2
                        } else {
                            weights.len() - 1 - i / 2
                        }];
                        assert_eq!(lazy.reached_by(w), w >= c - TOLERANCE, "k={k} p={p} w={w}");
                    }
                }
            }
        }
    }

    /// The scan of one cluster without shortcuts: every window of every round is
    /// compared with the critical weight of the complete bisection.
    fn scan_directly(snps: &[Snp], p: f64, threshold: f64) -> Vec<u32> {
        let mut snps = snps.to_vec();
        let mut masked = Vec::new();
        loop {
            let mut prefix = vec![0.0];
            for snp in &snps {
                prefix.push(prefix[prefix.len() - 1] + snp.weight);
            }
            let Some(c) = critical_weight(prefix[snps.len()], p, threshold) else {
                break;
            };
            let mut flagged = vec![false; snps.len()];
            for_each_window(&snps, |start, end| {
                if prefix[end] - prefix[start] >= c - TOLERANCE {
                    flagged[start..end].fill(true);
                }
            });
            if !flagged.contains(&true) {
                break;
            }
            masked.extend(
                (0..snps.len())
                    .filter(|&i| flagged[i])
                    .map(|i| snps[i].column),
            );
            let mut i = 0;
            snps.retain(|_| {
                i += 1;
                !flagged[i - 1]
            });
        }
        masked.sort_unstable();
        masked
    }

    #[test]
    fn shortcuts_give_the_direct_results() {
        let mut state = 11u64;
        let mut next = move |m: u64| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) % m
        };
        // 1 Mb genome in three contigs, 10 kb windows.
        let (genome, window, ends) = (1_000_000, 10_000, [300_000, 310_000, 1_000_000]);
        let p = window as f64 / genome as f64;
        let mut flagged_clusters = 0;
        for trial in 0..3000 {
            let threshold = [0.05, 1e-3, 1e-6][trial % 3] / (genome as f64 / window as f64);
            // Scattered SNPs and a few denser stretches, whole or split.
            let mut snps: Vec<(u32, u32, f64)> = Vec::new();
            for i in 0..2 + next(60) as u32 {
                let pos = if next(3) == 0 {
                    290_000 + next(30_000) as u32
                } else {
                    next(genome as u64) as u32
                };
                snps.push((pos, i, 1.0 / (1 + next(4) * next(3)) as f64));
            }
            snps.sort_unstable_by_key(|s| (s.0, s.1));
            let placed = place(snps.iter().copied(), &ends, window);
            let direct = scan_directly(&placed, p, threshold);
            assert_eq!(
                scan(&snps, &ends, window, p, threshold).0,
                direct,
                "{trial}"
            );
            flagged_clusters += !direct.is_empty() as usize;
            // The bound rules out only clusters without significant windows.
            let mut bound = WindowBound::new(genome, window);
            for &(pos, _, w) in &snps {
                bound.add(bound.bin(pos), w as f32);
            }
            assert!(bound.may_be_significant(p, threshold) || direct.is_empty());
        }
        assert!(flagged_clusters > 300, "{flagged_clusters}");
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
        let snps = [
            (1_000, 0, 1.0),
            (2_000, 1, 0.5),
            (3_000, 2, 0.5),
            (900_000, 3, 1.0),
        ];
        assert_eq!(
            scan(&snps, &[genome], window, p, threshold).0,
            vec![0, 1, 2]
        );
        // A single SNP is never significant on its own, whatever its weight.
        let lone = [(1_000, 0, 1.0), (900_000, 1, 1.0)];
        assert!(scan(&lone, &[genome], window, p, threshold).0.is_empty());
    }
}

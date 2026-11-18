//! Minimal neighbor-joining implementation used for the optional `--nj` tree.
//!
//! The unrooted NJ tree is written rooted at its midpoint: halfway along the longest
//! path between two leaves.
//!
//! Distances are JC69-corrected pairwise nucleotide differences, ignoring gaps
//! and `N`. The alignment is bit-sliced so each pair comparison handles
//! 64 columns at a time, and invariant gap-free/all-gap columns are folded away.
//! The NJ reduction uses a compact matrix and incrementally maintained row sums.
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::fmt::Write as _;

/// Code for gaps, `N` and anything that is not A, C, G or T.
const GAP: u8 = 4;
/// Two bit planes for the base code plus one validity plane.
const PLANES: usize = 3;
const VALID: usize = 2;

const NT_CODE: [u8; 256] = {
    let mut table = [GAP; 256];
    let order = *b"ACGT";
    let mut i = 0;
    while i < order.len() {
        table[order[i] as usize] = i as u8;
        i += 1;
    }
    table
};

/// Jukes-Cantor (JC69) corrected distance.
fn corrected_distance(valid_sites: u64, mismatches: u64) -> f64 {
    const B: f64 = 0.75;
    if valid_sites == 0 {
        return 0.0;
    }
    let p = mismatches as f64 / valid_sites as f64;
    let pc = (p / B).min(0.999_999_999);
    -B * (1.0 - pc).ln()
}

struct BitAlignment {
    n_taxa: usize,
    words: usize,
    planes: Vec<u64>,
    const_cols: u64,
}

impl BitAlignment {
    #[inline]
    fn taxon(&self, i: usize) -> &[u64] {
        let stride = self.words * PLANES;
        &self.planes[i * stride..(i + 1) * stride]
    }

    #[inline]
    fn counts(&self, i: usize, j: usize) -> (u64, u64) {
        let (a, b) = (self.taxon(i), self.taxon(j));
        let (mut valid_sites, mut mismatches) = (self.const_cols, 0u64);
        let (a_blocks, a_remainder) = a.as_chunks::<PLANES>();
        let (b_blocks, b_remainder) = b.as_chunks::<PLANES>();
        debug_assert!(a_remainder.is_empty() && b_remainder.is_empty());
        for (x, y) in a_blocks.iter().zip(b_blocks) {
            let both = x[VALID] & y[VALID];
            let diff = (x[0] ^ y[0]) | (x[1] ^ y[1]);
            valid_sites += both.count_ones() as u64;
            mismatches += (diff & both).count_ones() as u64;
        }
        (valid_sites, mismatches)
    }
}

fn build_bit_alignment(seqs: &[Vec<u8>], len: usize, pool: &ThreadPool) -> BitAlignment {
    // 0 = all gap, 1 = constant and gap-free, 2 = potentially pair-dependent.
    const BLOCK: usize = 1 << 16;
    let mut class = vec![0u8; len];
    pool.install(|| {
        class
            .par_chunks_mut(BLOCK)
            .enumerate()
            .for_each(|(block, out)| {
                let range = block * BLOCK..block * BLOCK + out.len();
                let first: Vec<u8> = seqs[0][range.clone()]
                    .iter()
                    .map(|&byte| NT_CODE[byte as usize])
                    .collect();
                // Flags: 1 = gap, 2 = difference, 4 = residue.
                let mut flags: Vec<u8> = first
                    .iter()
                    .map(|&code| if code == GAP { 1 } else { 4 })
                    .collect();
                for seq in &seqs[1..] {
                    for ((&byte, &reference), flag) in
                        seq[range.clone()].iter().zip(&first).zip(&mut flags)
                    {
                        let code = NT_CODE[byte as usize];
                        *flag |= if code == GAP {
                            1
                        } else {
                            4 | (((code != reference) as u8) << 1)
                        };
                    }
                }
                for (kind, &flag) in out.iter_mut().zip(&flags) {
                    *kind = if flag & 4 == 0 {
                        0
                    } else if flag & 3 == 0 {
                        1
                    } else {
                        2
                    };
                }
            });
    });

    let kept_cols: Vec<usize> = class
        .iter()
        .enumerate()
        .filter_map(|(column, &kind)| (kind == 2).then_some(column))
        .collect();
    let const_cols = class.iter().filter(|&&kind| kind == 1).count() as u64;
    let words = kept_cols.len().div_ceil(64);
    let stride = words * PLANES;
    let mut planes = vec![0u64; seqs.len() * stride];
    if stride != 0 {
        pool.install(|| {
            planes
                .par_chunks_mut(stride)
                .zip(seqs)
                .for_each(|(out, seq)| {
                    for (block, columns) in kept_cols.chunks(64).enumerate() {
                        let mut acc = [0u64; PLANES];
                        for (bit, &column) in columns.iter().enumerate() {
                            let code = NT_CODE[seq[column] as usize];
                            let valid = (code != GAP) as u64;
                            let value = code as u64 * valid;
                            for (plane, slot) in acc[..VALID].iter_mut().enumerate() {
                                *slot |= ((value >> plane) & 1) << bit;
                            }
                            acc[VALID] |= valid << bit;
                        }
                        out[block * PLANES..(block + 1) * PLANES].copy_from_slice(&acc);
                    }
                });
        });
    }
    BitAlignment {
        n_taxa: seqs.len(),
        words,
        planes,
        const_cols,
    }
}

fn compute_distance_matrix(align: &BitAlignment, pool: &ThreadPool) -> Vec<Vec<f64>> {
    pool.install(|| {
        (0..align.n_taxa)
            .into_par_iter()
            .map(|i| {
                ((i + 1)..align.n_taxa)
                    .map(|j| {
                        let (valid, mismatches) = align.counts(i, j);
                        corrected_distance(valid, mismatches)
                    })
                    .collect()
            })
            .collect()
    })
}

const Q_BLOCK: usize = 8;
const ROW_SUM_REFRESH: usize = 128;

#[inline]
fn row_sum(row: &[f64]) -> f64 {
    let mut acc = [0.0; 4];
    let (blocks, remainder) = row.as_chunks::<4>();
    for block in blocks {
        for (sum, &value) in acc.iter_mut().zip(block) {
            *sum += value;
        }
    }
    (acc[0] + acc[1]) + (acc[2] + acc[3]) + remainder.iter().sum::<f64>()
}

struct Nj {
    dist: Vec<f64>,
    row_sums: Vec<f64>,
    new_row: Vec<f64>,
    slot: Vec<u32>,
    children: Vec<[u32; 2]>,
    lengths: Vec<[f64; 2]>,
}

impl Nj {
    fn new(n: usize) -> Self {
        Self {
            dist: vec![0.0; n * n],
            row_sums: vec![0.0; n],
            new_row: vec![0.0; n],
            slot: Vec::with_capacity(n),
            children: Vec::with_capacity(n),
            lengths: Vec::with_capacity(n),
        }
    }

    fn load(&mut self, n: usize, rows: &[Vec<f64>]) {
        for (i, row) in rows.iter().enumerate().take(n) {
            for (offset, j) in ((i + 1)..n).enumerate() {
                let distance = row[offset];
                self.dist[i * n + j] = distance;
                self.dist[j * n + i] = distance;
            }
        }
    }

    fn reduce(&mut self, n: usize) {
        let (dist, sums, new_row, slot) = (
            &mut self.dist[..],
            &mut self.row_sums[..],
            &mut self.new_row[..],
            &mut self.slot,
        );
        slot.extend(0..n as u32);
        let mut m = n;
        for i in 0..m {
            sums[i] = row_sum(&dist[i * n..i * n + m]);
        }
        let mut since_refresh = 0;

        while m > 2 {
            let factor = (m - 2) as f64;
            let mut best_q = f64::INFINITY;
            let mut best = (0, 1);
            let mut best_key = (u32::MAX, u32::MAX);
            for a in 0..m - 1 {
                let row = &dist[a * n..a * n + m];
                let mut b = a + 1;
                while b + Q_BLOCK <= m {
                    let mut q = [0.0; Q_BLOCK];
                    let mut low = f64::INFINITY;
                    for (offset, value) in q.iter_mut().enumerate() {
                        *value = factor * row[b + offset] - sums[a] - sums[b + offset];
                        low = low.min(*value);
                    }
                    if low <= best_q {
                        for (offset, &value) in q.iter().enumerate() {
                            let key = pair_key(slot[a], slot[b + offset]);
                            if value < best_q || (value == best_q && key < best_key) {
                                (best_q, best, best_key) = (value, (a, b + offset), key);
                            }
                        }
                    }
                    b += Q_BLOCK;
                }
                for b in b..m {
                    let value = factor * row[b] - sums[a] - sums[b];
                    let key = pair_key(slot[a], slot[b]);
                    if value < best_q || (value == best_q && key < best_key) {
                        (best_q, best, best_key) = (value, (a, b), key);
                    }
                }
            }

            let (a, b) = best;
            let dab = dist[a * n + b];
            let la = (0.5 * dab + (sums[a] - sums[b]) / (2.0 * factor)).max(0.0);
            let lb = (dab - la).max(0.0);
            self.children.push([slot[a], slot[b]]);
            self.lengths.push([la, lb]);
            slot[a] = (n + self.children.len() - 1) as u32;

            let mut new_sum = 0.0;
            for t in 0..m {
                if t != a && t != b {
                    let (dat, dbt) = (dist[a * n + t], dist[b * n + t]);
                    let dut = 0.5 * (dat + dbt - dab);
                    new_row[t] = dut;
                    sums[t] = sums[t] - dat - dbt + dut;
                    new_sum += dut;
                }
            }
            for t in 0..m {
                if t != a && t != b {
                    dist[a * n + t] = new_row[t];
                    dist[t * n + a] = new_row[t];
                }
            }
            dist[a * n + a] = 0.0;
            sums[a] = new_sum;

            let last = m - 1;
            if b != last {
                dist.copy_within(last * n..last * n + m, b * n);
                dist[b * n + b] = 0.0;
                for t in 0..m {
                    dist[t * n + b] = dist[b * n + t];
                }
                sums[b] = sums[last];
                slot[b] = slot[last];
            }
            m -= 1;
            since_refresh += 1;
            if since_refresh == ROW_SUM_REFRESH {
                since_refresh = 0;
                for i in 0..m {
                    sums[i] = row_sum(&dist[i * n..i * n + m]);
                }
            }
        }

        let len = (dist[1] * 0.5).max(0.0);
        self.children.push([slot[0], slot[1]]);
        self.lengths.push([len, len]);
    }
}

#[inline]
fn pair_key(a: u32, b: u32) -> (u32, u32) {
    (a.min(b), a.max(b))
}

/// Unrooted tree made by the joins, as adjacency lists of `(neighbour, branch length)`.
///
/// Leaves are `0..n` and the internal node made by join `j` is `n + j`. The last join
/// only links the two remaining clusters, so it adds a branch rather than a node.
fn unrooted_tree(n: usize, nj: &Nj) -> Vec<Vec<(usize, f64)>> {
    let joins = nj.children.len();
    let mut adjacency = vec![Vec::new(); n + joins - 1];
    for (j, (children, lengths)) in nj.children.iter().zip(&nj.lengths).enumerate() {
        let [a, b] = children.map(|child| child as usize);
        if j + 1 == joins {
            let len = lengths[0] + lengths[1];
            adjacency[a].push((b, len));
            adjacency[b].push((a, len));
        } else {
            for (child, len) in [(a, lengths[0]), (b, lengths[1])] {
                adjacency[n + j].push((child, len));
                adjacency[child].push((n + j, len));
            }
        }
    }
    adjacency
}

/// Distance from `start` to every node, and the next node on the way back to `start`.
fn distances_from(adjacency: &[Vec<(usize, f64)>], start: usize) -> (Vec<f64>, Vec<usize>) {
    let mut distance = vec![0.0; adjacency.len()];
    let mut towards_start = vec![start; adjacency.len()];
    let mut stack = vec![start];
    while let Some(node) = stack.pop() {
        for &(next, len) in &adjacency[node] {
            if next != towards_start[node] {
                distance[next] = distance[node] + len;
                towards_start[next] = node;
                stack.push(next);
            }
        }
    }
    (distance, towards_start)
}

/// Leaf other than `start` that is farthest from it (the first one in case of ties).
fn farthest_leaf(distance: &[f64], n: usize, start: usize) -> usize {
    (0..n)
        .filter(|&leaf| leaf != start)
        .reduce(|best, leaf| {
            if distance[leaf] > distance[best] {
                leaf
            } else {
                best
            }
        })
        .expect("a tree has at least 2 leaves")
}

/// Where the tree is rooted: at a node, or on the branch between two nodes, at the
/// given distances from each.
enum Root {
    Node(usize),
    Branch([(usize, f64); 2]),
}

/// Midpoint root: halfway along the longest path between two leaves.
fn midpoint_root(adjacency: &[Vec<(usize, f64)>], n: usize) -> Root {
    let first = farthest_leaf(&distances_from(adjacency, 0).0, n, 0);
    let (distance, towards_first) = distances_from(adjacency, first);
    let last = farthest_leaf(&distance, n, first);
    let half = 0.5 * distance[last];
    // Walk back from the far end of the path to the branch that holds the midpoint.
    let mut node = last;
    while distance[towards_first[node]] > half {
        node = towards_first[node];
    }
    let next = towards_first[node];
    // A midpoint on an internal node (up to rounding) roots the tree at that node.
    let tolerance = 1e-12 * half.max(1.0);
    for candidate in [next, node] {
        if candidate >= n && (distance[candidate] - half).abs() <= tolerance {
            return Root::Node(candidate);
        }
    }
    Root::Branch([
        (next, (half - distance[next]).max(0.0)),
        (node, (distance[node] - half).max(0.0)),
    ])
}

/// Newick tree rooted at the midpoint (see [`midpoint_root`]).
fn midpoint_newick(names: &[String], adjacency: &[Vec<(usize, f64)>]) -> String {
    let mut out = String::new();
    match midpoint_root(adjacency, names.len()) {
        Root::Node(root) => write_subtree(names, adjacency, root, usize::MAX, &mut out),
        Root::Branch([(a, to_a), (b, to_b)]) => {
            out.push('(');
            write_subtree(names, adjacency, a, b, &mut out);
            let _ = write!(out, ":{to_a:.6},");
            write_subtree(names, adjacency, b, a, &mut out);
            let _ = write!(out, ":{to_b:.6})");
        }
    }
    out.push(';');
    out
}

fn escape_name(name: &str) -> String {
    if name
        .chars()
        .any(|character| matches!(character, ' ' | ':' | '(' | ')' | ',' | ';'))
    {
        format!("'{}'", name.replace('\'', "''"))
    } else {
        name.to_string()
    }
}

/// Subtree of `node` seen from its neighbour `parent`.
fn write_subtree(
    names: &[String],
    adjacency: &[Vec<(usize, f64)>],
    node: usize,
    parent: usize,
    out: &mut String,
) {
    if node < names.len() {
        out.push_str(&escape_name(&names[node]));
        return;
    }
    out.push('(');
    let mut first = true;
    for &(child, len) in &adjacency[node] {
        if child != parent {
            if !first {
                out.push(',');
            }
            first = false;
            write_subtree(names, adjacency, child, node, out);
            let _ = write!(out, ":{len:.6}");
        }
    }
    out.push(')');
}

fn validate_alignment(names: &[String], seqs: &[Vec<u8>]) -> anyhow::Result<usize> {
    anyhow::ensure!(
        names.len() == seqs.len(),
        "NJ tree: {} names for {} sequences",
        names.len(),
        seqs.len()
    );
    anyhow::ensure!(names.len() >= 2, "NJ tree: at least 2 sequences are needed");
    let len = seqs[0].len();
    anyhow::ensure!(
        seqs.iter().all(|seq| seq.len() == len),
        "NJ tree: all sequences must have the same length"
    );
    Ok(len)
}

/// Midpoint-rooted neighbor-joining tree of an alignment, in Newick format.
pub fn nj_tree_newick(
    names: &[String],
    seqs: &[Vec<u8>],
    num_threads: usize,
) -> anyhow::Result<String> {
    let len = validate_alignment(names, seqs)?;
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads.max(1))
        .build()?;
    let alignment = build_bit_alignment(seqs, len, &pool);
    let mut nj = Nj::new(names.len());
    nj.load(names.len(), &compute_distance_matrix(&alignment, &pool));
    nj.reduce(names.len());
    Ok(midpoint_newick(names, &unrooted_tree(names.len(), &nj)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_distance(a: &[u8], b: &[u8]) -> f64 {
        let (mut valid, mut mismatches) = (0, 0);
        for (&a, &b) in a.iter().zip(b) {
            let (a, b) = (NT_CODE[a as usize], NT_CODE[b as usize]);
            if a < GAP && b < GAP {
                valid += 1;
                mismatches += (a != b) as u64;
            }
        }
        corrected_distance(valid, mismatches)
    }

    #[test]
    fn bit_slicing_preserves_pairwise_distances() {
        let seqs = vec![
            b"AACG-NTA-GGCAAT".to_vec(),
            b"ACCGANTA-GGTAAT".to_vec(),
            b"GGCGA-TC-GGCAAT".to_vec(),
            b"GGCTA-TC-GGCAAC".to_vec(),
        ];
        let pool = ThreadPoolBuilder::new().num_threads(2).build().unwrap();
        let alignment = build_bit_alignment(&seqs, seqs[0].len(), &pool);
        let rows = compute_distance_matrix(&alignment, &pool);
        for i in 0..seqs.len() {
            for j in i + 1..seqs.len() {
                assert!((rows[i][j - i - 1] - naive_distance(&seqs[i], &seqs[j])).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn jc69_distance_values() {
        assert_eq!(corrected_distance(0, 0), 0.0);
        assert_eq!(corrected_distance(10, 0), 0.0);
        let d = corrected_distance(100, 10);
        assert!((d - (-0.75 * (1.0f64 - 0.4 / 3.0).ln())).abs() < 1e-12);
    }

    fn names(list: &[&str]) -> Vec<String> {
        list.iter().map(|name| name.to_string()).collect()
    }

    /// Unrooted tree from `(node, node, length)` branches.
    fn tree(nodes: usize, branches: &[(usize, usize, f64)]) -> Vec<Vec<(usize, f64)>> {
        let mut adjacency = vec![Vec::new(); nodes];
        for &(a, b, len) in branches {
            adjacency[a].push((b, len));
            adjacency[b].push((a, len));
        }
        adjacency
    }

    #[test]
    fn midpoint_splits_the_branch_halfway_along_the_longest_path() {
        // A and B join at node 4, C and D at node 5; the longest path is A-D (7).
        let adjacency = tree(
            6,
            &[
                (4, 0, 1.0),
                (4, 1, 1.0),
                (5, 2, 1.0),
                (5, 3, 5.0),
                (4, 5, 1.0),
            ],
        );
        assert_eq!(
            midpoint_newick(&names(&["A", "B", "C", "D"]), &adjacency),
            "(D:3.500000,(C:1.000000,(A:1.000000,B:1.000000):1.000000):1.500000);"
        );
    }

    #[test]
    fn a_midpoint_on_an_internal_node_roots_the_tree_there() {
        let adjacency = tree(4, &[(3, 0, 1.0), (3, 1, 1.0), (3, 2, 0.5)]);
        assert_eq!(
            midpoint_newick(&names(&["A", "B", "C"]), &adjacency),
            "(A:1.000000,B:1.000000,C:0.500000);"
        );
    }

    #[test]
    fn identical_sequences_give_a_valid_tree() {
        let seqs = vec![b"ACGT".to_vec(); 3];
        assert_eq!(
            nj_tree_newick(&names(&["A", "B"]), &seqs[..2], 1).unwrap(),
            "(B:0.000000,A:0.000000);"
        );
        let tree = nj_tree_newick(&names(&["A", "B", "C"]), &seqs, 1).unwrap();
        for leaf in ["A", "B", "C"] {
            assert_eq!(tree.matches(leaf).count(), 1, "{tree}");
        }
    }

    #[test]
    fn output_is_independent_of_thread_count() {
        let names = ["A", "B", "C", "D"].map(str::to_string);
        let seqs = vec![
            b"AAAA----ACGT".to_vec(),
            b"AAAC----ACGT".to_vec(),
            b"CCCC----ACGA".to_vec(),
            b"CCCG----ACGA".to_vec(),
        ];
        assert_eq!(
            nj_tree_newick(&names, &seqs, 1).unwrap(),
            nj_tree_newick(&names, &seqs, 4).unwrap()
        );
    }
}

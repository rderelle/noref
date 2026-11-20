//! Recombination plot for [Phandango](https://github.com/jameshadfield/phandango)
//! (`--plot`).
//!
//! The SNP columns removed from the alignment after recombination masking (those
//! failing the second `-f` check, see [`crate::recombination::failing_columns`]) are
//! placed on the genome of a reference isolate, at the position of their variant
//! group in that genome; a column cannot be placed when the reference lacks that
//! group. In each column, every isolate whose base (before masking) differs from the
//! ancestral base is reported as carrying a recombinant SNP (masking covers both
//! sides of a SNP, so it does not tell which isolates carry the change). The ancestral
//! base is the most parsimonious base (Fitch-Hartigan, missing bases matching any) at
//! the root of the midpoint-rooted NJ tree of the filtered alignment, or the most
//! frequent of several equally parsimonious ones.
//!
//! For each isolate, SNPs less than one window (`-w`) apart on the same contig of the
//! reference are merged into one block, from the first to the last SNP. Blocks are
//! written in the GFF format of Gubbins recombination predictions, which Phandango
//! displays next to a tree (for example an `--nj` tree) whose leaves carry the same
//! isolate names.
//!
//! Positions are coordinates in the reference's contigs placed end to end (1-based
//! in the GFF file). They are those of the variant groups, i.e. the middle of the
//! region between their two anchors, so a SNP lies at most half the length of that
//! region (see `-l`) from its reported position.
use anyhow::Result;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::dna::{encode, INVALID};
use crate::group_filtering::NO_POSITION;
use crate::phylo::Postorder;

/// SNPs of the removed columns, placed on the reference genome.
#[derive(Debug, Default, PartialEq, Eq)]
pub(crate) struct MappedSnps {
    /// Per isolate: reference positions of its non-ancestral bases.
    pub(crate) positions: Vec<Vec<u32>>,
    /// Removed columns.
    pub(crate) removed: usize,
    /// Removed columns placed on the reference (it has their variant group).
    pub(crate) mapped: usize,
}

/// Where the variant groups of the reference isolate lie in its genome.
pub(crate) struct Reference<'a> {
    /// Alignment column just after the last column of each variant group.
    pub(crate) group_ends: &'a [u32],
    /// Genome position of each variant group in the reference ([`NO_POSITION`] if
    /// the reference lacks the group).
    pub(crate) positions: &'a [u32],
}

/// The removed columns that can be placed on the reference: their position and their
/// base codes before masking (one per isolate), kept until the tree is built.
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct RemovedSnps {
    removed: usize,
    placed: Vec<(u32, Vec<u8>)>,
}

impl RemovedSnps {
    /// `rows` are the alignment rows after masking, before the removed `columns`
    /// (ascending) are taken out; `original_bases` holds, per isolate, the masked
    /// columns (ascending) with their base before masking.
    pub(crate) fn new(
        rows: &[Vec<u8>],
        columns: &[u32],
        original_bases: &[Vec<(u32, u8)>],
        reference: &Reference,
    ) -> Self {
        let mut placed: Vec<(u32, u32, Vec<u8>)> = columns
            .iter()
            .filter_map(|&column| {
                let group = reference.group_ends.partition_point(|&end| end <= column);
                let position = *reference.positions.get(group)?;
                let bases = rows
                    .iter()
                    .map(|row| encode(row[column as usize]))
                    .collect();
                (position != NO_POSITION).then_some((column, position, bases))
            })
            .collect();
        for (isolate, masked) in original_bases.iter().enumerate() {
            let mut k = 0;
            for &(column, base) in masked {
                while k < placed.len() && placed[k].0 < column {
                    k += 1;
                }
                if k < placed.len() && placed[k].0 == column {
                    placed[k].2[isolate] = encode(base);
                }
            }
        }
        let placed = placed.into_iter().map(|(_, p, codes)| (p, codes)).collect();
        Self {
            removed: columns.len(),
            placed,
        }
    }

    /// Report the non-ancestral bases on `tree`.
    pub(crate) fn place(&self, tree: &Postorder) -> MappedSnps {
        let n = self.placed.first().map_or(0, |(_, codes)| codes.len());
        let mut positions = vec![Vec::new(); n];
        let mut counts = vec![[0; 4]; tree.len()];
        for (position, codes) in &self.placed {
            let ancestral = ancestral_base(codes, tree, &mut counts);
            for (isolate, &code) in codes.iter().enumerate() {
                if code != INVALID && code != ancestral {
                    positions[isolate].push(*position);
                }
            }
        }
        MappedSnps {
            positions,
            removed: self.removed,
            mapped: self.placed.len(),
        }
    }
}

/// Ancestral base of a column (base codes of the leaves of `tree`); `counts` (one per
/// node) is scratch space.
fn ancestral_base(codes: &[u8], tree: &Postorder, counts: &mut [[u8; 4]]) -> u8 {
    counts.fill([0; 4]);
    for &(node, parent) in tree {
        // The bases of a node (bits) are those of the most children.
        let set: u8 = match codes.get(node) {
            Some(&code) if code < 4 => 1 << code,
            Some(_) => 15,
            None => {
                let max = counts[node].iter().max().copied();
                (0..4)
                    .filter(|&b| Some(counts[node][b]) == max)
                    .fold(0, |set, b| set | 1 << b)
            }
        };
        let Some(parent) = parent else {
            // The most frequent of the root's bases (the first one among equals).
            let frequency = |b: u8| codes.iter().filter(|&&c| c == b).count();
            return (0..4)
                .filter(|&b| set >> b & 1 == 1)
                .rev()
                .max_by_key(|&b| frequency(b))
                .unwrap();
        };
        for (b, count) in counts[parent].iter_mut().enumerate() {
            *count += set >> b & 1;
        }
    }
    unreachable!("a tree has a root")
}

/// SNPs of one isolate merged into one block (0-based positions, both included).
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct Block {
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) snps: usize,
}

/// Merge positions into blocks: after sorting, a position starts a new block when it
/// is at least `window` bases after the previous one or on another contig.
fn merge_into_blocks(positions: &mut [u32], contig_ends: &[u32], window: u32) -> Vec<Block> {
    positions.sort_unstable();
    let contig = |p: u32| contig_ends.partition_point(|&end| end <= p);
    let mut blocks: Vec<Block> = Vec::new();
    for &p in positions.iter() {
        match blocks.last_mut() {
            Some(b) if p - b.end < window && contig(p) == contig(b.end) => {
                b.end = p;
                b.snps += 1;
            }
            _ => blocks.push(Block {
                start: p,
                end: p,
                snps: 1,
            }),
        }
    }
    blocks
}

/// Counts printed once the plot is written.
#[derive(Debug, Default, PartialEq, Eq)]
pub(crate) struct PlotStats {
    pub(crate) blocks: usize,
    /// Isolates with at least one block.
    pub(crate) isolates: usize,
    pub(crate) snps: usize,
    /// Isolates with blocks whose name holds a space or a double quote, which
    /// Phandango cannot match with the tree.
    pub(crate) unreadable_names: usize,
}

/// Write the blocks of every isolate in the GFF format of Gubbins.
fn write_gff(
    w: &mut impl Write,
    names: &[String],
    reference: usize,
    contig_ends: &[u32],
    mut mapped: MappedSnps,
    window: u32,
) -> io::Result<PlotStats> {
    let genome = contig_ends.last().copied().unwrap_or(0);
    // Phandango recognises Gubbins output by the source field of the 4th line, so
    // the header has exactly three lines and the first block is on the 4th one.
    writeln!(w, "##gff-version 3")?;
    writeln!(w, "##sequence-region SEQUENCE 1 {genome}")?;
    writeln!(
        w,
        "# noref {}: non-ancestral bases of the alignment columns removed after \
         recombination masking, on the genome of {} (contigs end to end), merged into \
         blocks when less than {} bp apart",
        env!("CARGO_PKG_VERSION"),
        names[reference],
        window
    )?;
    let mut stats = PlotStats::default();
    for (isolate, positions) in mapped.positions.iter_mut().enumerate() {
        let blocks = merge_into_blocks(positions, contig_ends, window);
        if blocks.is_empty() {
            continue;
        }
        let name = &names[isolate];
        stats.isolates += 1;
        stats.unreadable_names += name.contains(|c: char| c.is_whitespace() || c == '"') as usize;
        for b in blocks {
            stats.blocks += 1;
            stats.snps += b.snps;
            writeln!(
                w,
                "SEQUENCE\tGUBBINS\tCDS\t{}\t{}\t0.000\t.\t0\tnode=\"{name}\";\
                 neg_log_likelihood=\"0\";taxa=\"{name}\";snp_count=\"{}\";program=\"noref\";",
                b.start + 1,
                b.end + 1,
                b.snps
            )?;
        }
    }
    Ok(stats)
}

/// Write the recombination plot of `mapped` (see [`RemovedSnps::place`]) to `path`.
/// `contig_ends` are those of the reference genome and `window` is in bases.
pub(crate) fn write_plot(
    path: &Path,
    names: &[String],
    reference: usize,
    contig_ends: &[u32],
    mapped: MappedSnps,
    window: u32,
) -> Result<PlotStats> {
    let mut w = BufWriter::new(File::create(path)?);
    let stats = write_gff(&mut w, names, reference, contig_ends, mapped, window)?;
    w.flush()?;
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(rows: &[&[u8]]) -> Vec<Vec<u8>> {
        rows.iter().map(|r| r.to_vec()).collect()
    }

    #[test]
    fn removed_columns_keep_their_bases_before_masking() {
        // Columns 0-1, 2 and 3 are 3 variant groups; the reference lacks the second.
        let rows = rows(&[b"NTCN", b"ANCA", b"G-CN"]);
        let original = vec![vec![(0, b'C'), (3, b'G')], vec![(1, b'G')], vec![(3, b'T')]];
        let reference = Reference {
            group_ends: &[2, 3, 4],
            positions: &[1_000, NO_POSITION, 50],
        };
        let removed = RemovedSnps::new(&rows, &[0, 1, 2, 3], &original, &reference);
        let codes = |bases: &[u8]| bases.iter().map(|&b| encode(b)).collect::<Vec<_>>();
        assert_eq!(removed.removed, 4);
        assert_eq!(
            removed.placed,
            [
                (1_000, codes(b"CAG")),
                (1_000, codes(b"TG-")),
                (50, codes(b"GAT"))
            ]
        );
    }

    #[test]
    fn bases_that_differ_from_the_root_are_reported() {
        // ((((3, 4), 2), 1), 0): internal nodes 5 to 7, root 8.
        let tree = vec![
            (0, Some(8)),
            (1, Some(7)),
            (2, Some(6)),
            (3, Some(5)),
            (4, Some(5)),
            (5, Some(6)),
            (6, Some(7)),
            (7, Some(8)),
            (8, None),
        ];
        // 10: the C of a clade, although more frequent than the ancestral A (also
        // with a missing base, 20). 30: C or A at the root, A is more frequent. 40:
        // any base at the root, T is the most frequent.
        let columns = [
            (10, b"AACCC"),
            (20, b"AACC-"),
            (30, b"CAAAA"),
            (40, b"ACGTT"),
        ];
        let removed = RemovedSnps {
            removed: 6,
            placed: columns
                .iter()
                .map(|&(p, bases)| (p, bases.iter().map(|&b| encode(b)).collect()))
                .collect(),
        };
        let expected = [
            vec![30, 40],
            vec![40],
            vec![10, 20, 40],
            vec![10, 20],
            vec![10],
        ];
        assert_eq!(removed.place(&tree).positions, expected);
    }

    #[test]
    fn close_snps_on_one_contig_are_merged() {
        // Two contigs: [0, 10_000) and [10_000, 30_000); 5 kb window.
        let mut positions = vec![9_000, 1_000, 4_000, 5_500, 10_500, 10_500, 20_000];
        let blocks = merge_into_blocks(&mut positions, &[10_000, 30_000], 5_000);
        let b = |start, end, snps| Block { start, end, snps };
        assert_eq!(
            blocks,
            vec![
                // 1_000, 4_000, 5_500 and 9_000: each less than 5 kb from the previous.
                b(1_000, 9_000, 4),
                // Same distance but on the next contig.
                b(10_500, 10_500, 2),
                b(20_000, 20_000, 1),
            ]
        );
        assert!(merge_into_blocks(&mut [], &[10_000], 5_000).is_empty());
        // Exactly one window apart: not merged.
        assert_eq!(
            merge_into_blocks(&mut [0, 5_000], &[10_000], 5_000).len(),
            2
        );
    }

    #[test]
    fn gff_is_read_as_gubbins_output_by_phandango() {
        let names: Vec<String> = ["ref", "a", "b c"].map(String::from).to_vec();
        let mapped = MappedSnps {
            positions: vec![vec![], vec![7, 0, 2], vec![900]],
            removed: 5,
            mapped: 4,
        };
        let mut out = Vec::new();
        let stats = write_gff(&mut out, &names, 0, &[500, 1_000], mapped, 5).unwrap();
        assert_eq!(
            stats,
            PlotStats {
                blocks: 3,
                isolates: 2,
                snps: 4,
                unreadable_names: 1
            }
        );
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[1], "##sequence-region SEQUENCE 1 1000");
        assert!(lines[2].starts_with("# noref ") && !lines[2].contains('\t'));
        // Phandango checks the source field of the 4th line.
        assert_eq!(lines[3].split('\t').nth(1), Some("GUBBINS"));
        assert_eq!(
            &lines[3..],
            [
                "SEQUENCE\tGUBBINS\tCDS\t1\t3\t0.000\t.\t0\tnode=\"a\";neg_log_likelihood=\"0\";taxa=\"a\";snp_count=\"2\";program=\"noref\";",
                "SEQUENCE\tGUBBINS\tCDS\t8\t8\t0.000\t.\t0\tnode=\"a\";neg_log_likelihood=\"0\";taxa=\"a\";snp_count=\"1\";program=\"noref\";",
                "SEQUENCE\tGUBBINS\tCDS\t901\t901\t0.000\t.\t0\tnode=\"b c\";neg_log_likelihood=\"0\";taxa=\"b c\";snp_count=\"1\";program=\"noref\";",
            ]
        );
    }
}

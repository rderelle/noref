//! Recombination plot for [Phandango](https://github.com/jameshadfield/phandango)
//! (`--plot`).
//!
//! The SNP columns removed from the alignment after recombination masking (those
//! failing the second `-f` check, see [`crate::recombination::failing_columns`]) are
//! placed on the genome of a reference isolate, at the position of their variant
//! group in that genome. In each of these columns, every other isolate in which the
//! SNP was masked and whose base (before masking) differs from the base of the
//! reference is reported as carrying a recombinant SNP. A column cannot be placed
//! when the reference has no base in it.
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

use crate::group_filtering::NO_POSITION;

/// SNPs of the removed columns, placed on the reference genome.
#[derive(Debug, Default, PartialEq, Eq)]
pub(crate) struct MappedSnps {
    /// Per isolate: reference positions of its masked SNPs that differ from the
    /// reference.
    pub(crate) positions: Vec<Vec<u32>>,
    /// Removed columns.
    pub(crate) removed: usize,
    /// Removed columns in which the reference has a base and a position.
    pub(crate) mapped: usize,
}

/// The reference isolate and where its variant groups lie in its genome.
pub(crate) struct Reference<'a> {
    /// Row of the reference in the alignment.
    pub(crate) index: usize,
    /// Alignment column just after the last column of each variant group.
    pub(crate) group_ends: &'a [u32],
    /// Genome position of each variant group in the reference ([`NO_POSITION`] if
    /// the reference lacks the group).
    pub(crate) positions: &'a [u32],
}

#[inline]
fn is_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Place the SNPs of the removed `columns` (ascending) on the reference genome.
/// `rows` are the alignment rows after masking but before the columns are removed,
/// and `original_bases` holds, per isolate, the masked columns (ascending) with their
/// base before masking. An isolate is reported in a column when it was masked there
/// and its base differs from the base of the reference (before masking).
pub(crate) fn map_removed_snps(
    rows: &[Vec<u8>],
    columns: &[u32],
    original_bases: &[Vec<(u32, u8)>],
    reference: &Reference,
) -> MappedSnps {
    debug_assert!(columns.is_sorted());
    // Position and base (before masking) of the reference in each removed column,
    // when it has both.
    let reference_masked = original_bases
        .get(reference.index)
        .map_or(&[][..], Vec::as_slice);
    let placed: Vec<Option<(u32, u8)>> = columns
        .iter()
        .map(|&column| {
            let group = reference.group_ends.partition_point(|&end| end <= column);
            let position = reference
                .positions
                .get(group)
                .copied()
                .unwrap_or(NO_POSITION);
            let base = match reference_masked.binary_search_by_key(&column, |&(c, _)| c) {
                Ok(i) => reference_masked[i].1,
                Err(_) => rows[reference.index][column as usize],
            };
            (position != NO_POSITION && is_base(base)).then_some((position, base))
        })
        .collect();
    // Masked SNPs of the other isolates in the removed columns.
    let mut positions = vec![Vec::new(); rows.len()];
    for (isolate, masked) in original_bases.iter().enumerate() {
        if isolate == reference.index {
            continue;
        }
        for &(column, base) in masked {
            let Ok(k) = columns.binary_search(&column) else {
                continue;
            };
            if let Some((position, reference_base)) = placed[k] {
                if is_base(base) && base != reference_base {
                    positions[isolate].push(position);
                }
            }
        }
    }
    MappedSnps {
        positions,
        removed: columns.len(),
        mapped: placed.iter().flatten().count(),
    }
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
        "# noref {}: SNPs of the alignment columns removed after recombination \
         masking, on the genome of {} (contigs end to end), merged into blocks when \
         less than {} bp apart",
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

/// Write the recombination plot of `mapped` (see [`map_removed_snps`]) to `path`.
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
    fn masked_isolates_with_another_base_than_the_reference_are_reported() {
        // Isolate 0 is the reference; 3 variant groups: columns 0-1, 2, 3-4.
        let rows = rows(&[b"ANNAC", b"GNNNT", b"NGCAN", b"-NN-T"]);
        // Masked cells with their base before masking; other `N`s are missing data.
        let original = vec![
            vec![(1, b'C')],
            vec![(3, b'G')],
            vec![(0, b'A'), (4, b'T')],
            vec![(1, b'C'), (2, b'C')],
        ];
        let reference = Reference {
            index: 0,
            group_ends: &[2, 3, 5],
            positions: &[1_000, NO_POSITION, 50],
        };
        let mapped = map_removed_snps(&rows, &[0, 1, 2, 3, 4], &original, &reference);
        // Column 2 cannot be placed: the reference lacks its group (and has no base).
        assert_eq!((mapped.removed, mapped.mapped), (5, 4));
        // Isolate 1: masked G in column 3 (A in the reference); its G in column 0 and
        // T in column 4 differ from the reference but were not masked.
        // Isolate 2: masked with the reference's base in column 0, masked T in column 4
        // (C in the reference).
        // Isolate 3: masked with the reference's base in column 1 (where the reference
        // was masked too), and in column 2, which cannot be placed.
        assert_eq!(mapped.positions, vec![vec![], vec![50], vec![50], vec![]]);
        // Masked SNPs of columns that are not removed are not reported.
        let mapped = map_removed_snps(&rows, &[0, 1, 2], &original, &reference);
        assert_eq!(mapped.mapped, 2);
        assert!(mapped.positions.iter().all(Vec::is_empty));
        // A column in which the reference has no base is not placed.
        let mapped = map_removed_snps(&rows, &[1], &[], &reference);
        assert_eq!((mapped.removed, mapped.mapped), (1, 0));
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

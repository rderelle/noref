//! Output writers for the alignments, missing-data report, optional neighbor-joining
//! trees and optional recombination plot.
use anyhow::Result;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

// Derive all output filenames from one prefix so `-o results/foo` writes
// `results/foo_filtered_alignment.fas`, `results/foo_filtered_missing.tsv`, and so on.
fn build_path(base: &Path, suffix: &str, ext: &str) -> PathBuf {
    let base_name = base.file_name().unwrap_or_else(|| OsStr::new("output"));
    let mut name = base_name.to_os_string();
    name.push(suffix);
    name.push(".");
    name.push(ext);
    let mut path = base.parent().map(Path::to_path_buf).unwrap_or_default();
    path.push(name);
    path
}

/// Paths of the output files, derived from the `-o` prefix.
pub struct OutputPaths {
    /// The prefix itself.
    pub prefix: PathBuf,
    /// Alignment before recombination masking.
    pub raw_alignment: PathBuf,
    /// Final alignment (recombinant SNPs masked, failing columns removed).
    pub filtered_alignment: PathBuf,
    /// Missing data per isolate in the raw alignment.
    pub raw_missing: PathBuf,
    /// Missing data per isolate in the filtered alignment.
    pub filtered_missing: PathBuf,
    /// NJ tree of the raw alignment (`--nj`).
    pub raw_tree: PathBuf,
    /// NJ tree of the filtered alignment (`--nj`).
    pub filtered_tree: PathBuf,
    /// Recombination plot (`--plot`).
    pub plot: PathBuf,
    /// Run log.
    pub log: PathBuf,
}

pub fn output_paths(out_base: &Path) -> OutputPaths {
    OutputPaths {
        prefix: out_base.to_path_buf(),
        raw_alignment: build_path(out_base, "_raw_alignment", "fas"),
        filtered_alignment: build_path(out_base, "_filtered_alignment", "fas"),
        raw_missing: build_path(out_base, "_raw_missing", "tsv"),
        filtered_missing: build_path(out_base, "_filtered_missing", "tsv"),
        raw_tree: build_path(out_base, "_raw_NJ", "tree"),
        filtered_tree: build_path(out_base, "_filtered_NJ", "tree"),
        plot: build_path(out_base, "_recombination", "gff"),
        log: build_path(out_base, "", "log"),
    }
}

/// Number of columns of an alignment, and of those without two different bases.
pub struct AlignmentSummary {
    pub length: usize,
    pub constant: usize,
}

/// Write an alignment in FASTA format (one row per isolate) and summarise it.
pub fn write_alignment(
    path: &Path,
    names: &[String],
    rows: &[Vec<u8>],
) -> Result<AlignmentSummary> {
    let total_len = rows.first().map(|v| v.len()).unwrap_or(0);
    let constant = (0..total_len)
        .filter(|&column| {
            let mut base = None;
            rows.iter().all(|row| match row[column] {
                b'-' | b'N' => true,
                b => {
                    base.is_none_or(|previous| previous == b) && {
                        base = Some(b);
                        true
                    }
                }
            })
        })
        .count();
    // FASTA output is wrapped at 60 characters for broad tool compatibility.
    let mut w = BufWriter::new(File::create(path)?);
    for (isolate, name) in names.iter().enumerate() {
        writeln!(w, ">{}", name)?;
        for ch in rows[isolate].chunks(60) {
            w.write_all(ch)?;
            w.write_all(b"\n")?;
        }
    }
    w.flush()?;
    Ok(AlignmentSummary {
        length: total_len,
        constant,
    })
}

/// Write the percentage of missing cells of each isolate in an alignment.
pub fn write_missing(path: &Path, names: &[String], rows: &[Vec<u8>]) -> Result<()> {
    // Per-isolate missingness helps users identify problematic samples.
    let mut tw = BufWriter::new(File::create(path)?);
    writeln!(tw, "isolate\tpct_missing")?;
    for (isolate, name) in names.iter().enumerate() {
        let len = rows[isolate].len().max(1);
        let m = rows[isolate]
            .iter()
            .filter(|&&b| b == b'-' || b == b'N')
            .count();
        writeln!(tw, "{}\t{:.1}", name, 100.0 * m as f64 / len as f64)?;
    }
    tw.flush()?;
    Ok(())
}

/// Build a neighbor-joining tree of an alignment and write it in Newick format. Trees
/// are optional because they add an O(n²·L) distance pass.
pub fn write_nj_tree(
    path: &Path,
    names: &[String],
    rows: &[Vec<u8>],
    threads: usize,
) -> Result<()> {
    let tree = crate::phylo::nj_tree_newick(names, rows, threads)?;
    let mut tw = BufWriter::new(File::create(path)?);
    writeln!(tw, "{}", tree)?;
    tw.flush()?;
    Ok(())
}

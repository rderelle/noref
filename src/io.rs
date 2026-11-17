//! Input discovery and FASTA opening helpers.
//!
//! Inputs can come from a directory of FASTA files or from a two-column TSV table
//! mapping isolate names to FASTA paths. Files ending in `.gz` (in any letter case) are
//! decompressed on the fly by [`open_fasta`].
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

/// Accepted genome FASTA extensions (optionally followed by `.gz`), in any letter case.
const FASTA_EXTENSIONS: &[&str] = &["fa", "fas", "fasta", "fna"];

/// One named genome input consumed by the analysis pipeline.
#[derive(Debug, Clone)]
pub struct IsolateInput {
    /// Isolate name, used in output FASTA headers and trees.
    pub name: String,
    /// FASTA path for this isolate; may be gzip-compressed.
    pub path: PathBuf,
}

fn is_gzipped(path: &Path) -> bool {
    path.extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

fn is_fasta_extension(ext: &str) -> bool {
    FASTA_EXTENSIONS.iter().any(|e| ext.eq_ignore_ascii_case(e))
}

/// Whether a file name ends with a FASTA extension, optionally followed by `.gz`.
fn has_fasta_extension(path: &Path) -> bool {
    let path = if is_gzipped(path) {
        match path.file_stem() {
            Some(stem) => Path::new(stem),
            None => return false,
        }
    } else {
        path
    };
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(is_fasta_extension)
}

/// FASTA files of a directory, sorted by file name. Symbolic links to files are
/// followed.
pub fn collect_fasta_files(input_dir: &Path) -> Result<Vec<PathBuf>> {
    let entries = std::fs::read_dir(input_dir)
        .with_context(|| format!("cannot read the directory {}", input_dir.display()))?;
    let mut files: Vec<PathBuf> = entries
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        // `Path::is_file` follows symbolic links, unlike `DirEntry::file_type`.
        .filter(|p| p.is_file() && has_fasta_extension(p))
        .collect();
    files.sort_unstable_by(|a, b| a.file_name().cmp(&b.file_name()));
    Ok(files)
}

/// Build isolate inputs from all supported FASTA files in a directory.
pub fn collect_isolates_from_dir(input_dir: &Path) -> Result<Vec<IsolateInput>> {
    let files = collect_fasta_files(input_dir)?;
    anyhow::ensure!(
        !files.is_empty(),
        "no FASTA file (.fa, .fas, .fasta or .fna, optionally gzipped) found in {}",
        input_dir.display()
    );
    let mut inputs: Vec<IsolateInput> = files
        .into_iter()
        .map(|path| IsolateInput {
            name: isolate_name_from_path(&path),
            path,
        })
        .collect();
    inputs.sort_unstable_by(|a, b| a.name.cmp(&b.name).then_with(|| a.path.cmp(&b.path)));
    ensure_unique_names(&inputs)?;
    Ok(inputs)
}

/// Build isolate inputs from a two-column `name<TAB>path` table. Relative paths are
/// resolved from the table's directory, which makes tables portable.
pub fn collect_isolates_from_table(table_path: &Path) -> Result<Vec<IsolateInput>> {
    let table = table_path.display();
    let file = File::open(table_path).with_context(|| format!("cannot open {table}"))?;
    let base_dir = table_path.parent().unwrap_or_else(|| Path::new("."));

    let mut inputs: Vec<IsolateInput> = Vec::new();
    for (line_idx, line) in BufReader::new(file).lines().enumerate() {
        let line = line.with_context(|| format!("cannot read {table}"))?;
        let line_no = line_idx + 1;
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }
        let columns: Vec<&str> = trimmed.split('\t').map(str::trim).collect();
        anyhow::ensure!(
            columns.len() == 2 && columns.iter().all(|c| !c.is_empty()),
            "line {line_no} of {table}: expected two tab-separated columns (isolate name and genome path)"
        );
        let mut path = PathBuf::from(columns[1]);
        if path.is_relative() {
            path = base_dir.join(path);
        }
        anyhow::ensure!(
            path.is_file(),
            "line {line_no} of {table}: {} is not a file",
            path.display()
        );
        inputs.push(IsolateInput {
            name: columns[0].to_string(),
            path,
        });
    }
    anyhow::ensure!(!inputs.is_empty(), "no genome listed in {table}");
    inputs.sort_unstable_by(|a, b| a.name.cmp(&b.name).then_with(|| a.path.cmp(&b.path)));
    ensure_unique_names(&inputs)?;
    Ok(inputs)
}

/// Fail if two inputs share a name (inputs are sorted by name, or come from both a
/// directory and a table).
pub fn ensure_unique_names(inputs: &[IsolateInput]) -> Result<()> {
    let mut names: Vec<&str> = inputs.iter().map(|i| i.name.as_str()).collect();
    names.sort_unstable();
    if let Some(pair) = names.windows(2).find(|w| w[0] == w[1]) {
        anyhow::bail!("isolate name '{}' is given more than once", pair[0]);
    }
    Ok(())
}

/// The file name without its compression and FASTA extensions (for example,
/// `isolate.1.fna.gz` becomes `isolate.1`).
pub(crate) fn isolate_name_from_path(p: &Path) -> String {
    let mut name = p
        .file_name()
        .and_then(|x| x.to_str())
        .unwrap_or("unnamed")
        .to_string();
    if is_gzipped(Path::new(&name)) {
        name.truncate(name.len() - 3);
    }
    if let Some((stem, ext)) = name.rsplit_once('.') {
        if is_fasta_extension(ext) {
            name.truncate(stem.len());
        }
    }
    name
}

/// Open a FASTA file, decompressing it if its name ends with `.gz`.
pub(crate) fn open_fasta(path: &Path) -> Result<Box<dyn Read>> {
    let f = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let buffered = BufReader::with_capacity(512 * 1024, f);
    Ok(if is_gzipped(path) {
        Box::new(MultiGzDecoder::new(buffered))
    } else {
        Box::new(buffered)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{write::GzEncoder, Compression};
    use std::io::Write;

    #[test]
    fn names_drop_fasta_and_compression_extensions_in_any_case() {
        for (file, name) in [
            ("isolate.1.fna.gz", "isolate.1"),
            ("G1.FNA.GZ", "G1"),
            ("G2.Fasta", "G2"),
            ("G3.txt.gz", "G3.txt"),
            ("G4", "G4"),
        ] {
            assert_eq!(isolate_name_from_path(Path::new(file)), name);
        }
        assert!(has_fasta_extension(Path::new("a.FA.GZ")));
        assert!(!has_fasta_extension(Path::new("a.txt.gz")));
        assert!(!has_fasta_extension(Path::new("gz")));
    }

    #[test]
    fn upper_case_gz_files_are_decompressed() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("G1.fna.GZ");
        let mut gz = GzEncoder::new(File::create(&path).unwrap(), Compression::default());
        gz.write_all(b">c1\nACGT\n").unwrap();
        gz.finish().unwrap();
        let mut text = String::new();
        open_fasta(&path)
            .unwrap()
            .read_to_string(&mut text)
            .unwrap();
        assert_eq!(text, ">c1\nACGT\n");
    }

    #[cfg(unix)]
    #[test]
    fn symbolic_links_to_genomes_are_followed() {
        let dir = tempfile::tempdir().unwrap();
        let (data, links) = (dir.path().join("data"), dir.path().join("links"));
        std::fs::create_dir_all(&data).unwrap();
        std::fs::create_dir_all(&links).unwrap();
        for name in ["G1.fna", "G2.fa"] {
            std::fs::write(data.join(name), ">c\nACGT\n").unwrap();
            std::os::unix::fs::symlink(data.join(name), links.join(name)).unwrap();
        }
        // A link to a directory is not a genome.
        std::os::unix::fs::symlink(&data, links.join("sub.fa")).unwrap();
        let names: Vec<String> = collect_isolates_from_dir(&links)
            .unwrap()
            .into_iter()
            .map(|i| i.name)
            .collect();
        assert_eq!(names, ["G1", "G2"]);
    }

    #[test]
    fn a_directory_without_fasta_files_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("notes.txt"), "x").unwrap();
        let error = collect_isolates_from_dir(dir.path())
            .unwrap_err()
            .to_string();
        assert!(error.starts_with("no FASTA file"), "{error}");
    }

    #[test]
    fn tables_are_checked() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("g.fna"), ">c\nACGT\n").unwrap();
        let table = dir.path().join("inputs.tsv");
        let read = |text: &str| {
            std::fs::write(&table, text).unwrap();
            collect_isolates_from_table(&table).map_err(|e| e.to_string())
        };
        assert_eq!(read("b\tg.fna\na\tg.fna\n").unwrap()[0].name, "a");
        assert!(read("a g.fna\n")
            .unwrap_err()
            .contains("two tab-separated columns"));
        assert!(read("a\tmissing.fna\n")
            .unwrap_err()
            .contains("is not a file"));
        assert!(read("a\tg.fna\na\tg.fna\n")
            .unwrap_err()
            .contains("more than once"));
        assert!(read("\n").unwrap_err().starts_with("no genome listed"));
    }

    #[test]
    fn names_must_be_unique_across_sources() {
        let input = |name: &str| IsolateInput {
            name: name.to_string(),
            path: PathBuf::from("x"),
        };
        assert!(ensure_unique_names(&[input("a"), input("b")]).is_ok());
        let error = ensure_unique_names(&[input("b"), input("a"), input("b")]).unwrap_err();
        assert_eq!(
            error.to_string(),
            "isolate name 'b' is given more than once"
        );
    }
}

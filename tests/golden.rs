//! Golden integration tests on a small simulated dataset (8 genomes of 20 kb on a
//! known tree, fragmented into contigs of random orientation, with lower-case and
//! ambiguous bases; two files are gzipped).

use std::fs;
use std::path::{Path, PathBuf};

use noref::{run_with_args, Args};

mod common;

const INPUT_DIR: &str = "tests/data/input";
const INPUT_TSV: &str = "tests/data/input.tsv";
const EXPECTED: &str = "tests/data/expected";

fn args(
    input: Option<PathBuf>,
    input_file: Option<PathBuf>,
    output: PathBuf,
    threads: usize,
) -> Args {
    Args {
        input,
        input_file,
        output: Some(output),
        k: 15,
        min_freq: 0.85,
        length_middle: 70,
        mask: 4,
        window: 5,
        alpha: 0.05,
        threads,
        nj: true,
        plot: None,
        skip_recombination: false,
    }
}

fn canonical(path: &str) -> PathBuf {
    Path::new(path).canonicalize().unwrap()
}

#[test]
fn golden_directory_input() {
    let out = tempfile::tempdir().unwrap();
    run_with_args(args(
        Some(canonical(INPUT_DIR)),
        None,
        out.path().join("noref"),
        1,
    ))
    .unwrap();
    common::compare_test_outputs(Path::new(EXPECTED), out.path());
}

#[test]
fn golden_table_input() {
    let out = tempfile::tempdir().unwrap();
    run_with_args(args(
        None,
        Some(canonical(INPUT_TSV)),
        out.path().join("noref"),
        2,
    ))
    .unwrap();
    common::compare_test_outputs(Path::new(EXPECTED), out.path());
}

/// Lines of a recombination plot, without the comment line (it names the version).
fn plot_lines(path: &Path) -> Vec<String> {
    let text = fs::read_to_string(path).unwrap();
    let mut lines: Vec<String> = text.lines().map(String::from).collect();
    assert!(lines[2].starts_with("# noref "), "{}", lines[2]);
    lines.remove(2);
    lines
}

#[test]
fn golden_recombination_plot() {
    const PLOT: &str = "noref_recombination.gff";
    let expected = plot_lines(&Path::new(EXPECTED).join(PLOT));
    for threads in [1, 4] {
        let out = tempfile::tempdir().unwrap();
        let mut a = args(
            Some(canonical(INPUT_DIR)),
            None,
            out.path().join("noref"),
            threads,
        );
        a.plot = Some("G7".to_string());
        run_with_args(a).unwrap();
        // The plot changes no other output.
        common::compare_test_outputs(Path::new(EXPECTED), out.path());
        assert_eq!(
            expected,
            plot_lines(&out.path().join(PLOT)),
            "{threads} threads"
        );
    }
}

#[test]
fn skipping_recombination_stops_at_the_raw_alignment() {
    let out = tempfile::tempdir().unwrap();
    let mut a = args(
        Some(canonical(INPUT_DIR)),
        None,
        out.path().join("noref"),
        2,
    );
    a.skip_recombination = true;
    run_with_args(a).unwrap();
    // Only the raw outputs are written, identical to those of a full run.
    let mut written: Vec<String> = fs::read_dir(out.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().into_string().unwrap())
        .collect();
    written.sort();
    assert_eq!(
        written,
        [
            "noref.log",
            "noref_raw_NJ.tree",
            "noref_raw_alignment.fas",
            "noref_raw_missing.tsv"
        ]
    );
    for name in [
        "noref_raw_alignment.fas",
        "noref_raw_missing.tsv",
        "noref_raw_NJ.tree",
    ] {
        assert_eq!(
            fs::read_to_string(Path::new(EXPECTED).join(name)).unwrap(),
            fs::read_to_string(out.path().join(name)).unwrap(),
            "{name} differs"
        );
    }
    // Missing data recomputed from the raw alignment.
    let fasta = fs::read_to_string(out.path().join("noref_raw_alignment.fas")).unwrap();
    let mut expected = String::from("isolate\tpct_missing\n");
    for record in fasta.split('>').skip(1) {
        let (name, seq) = record.split_once('\n').unwrap();
        let seq: String = seq.lines().collect();
        let missing = seq.bytes().filter(|&b| b == b'-' || b == b'N').count();
        let pct = 100.0 * missing as f64 / seq.len().max(1) as f64;
        expected += &format!("{name}\t{pct:.1}\n");
    }
    assert_eq!(
        fs::read_to_string(out.path().join("noref_raw_missing.tsv")).unwrap(),
        expected
    );
}

#[test]
fn plot_reference_must_be_an_input_isolate() {
    let out = tempfile::tempdir().unwrap();
    let mut a = args(
        Some(canonical(INPUT_DIR)),
        None,
        out.path().join("noref"),
        1,
    );
    a.plot = Some("G9".to_string());
    let error = run_with_args(a).unwrap_err().to_string();
    assert!(error.contains("'G9'"), "{error}");
    // Nothing was written but the log, which ends with the error.
    let written: Vec<_> = fs::read_dir(out.path()).unwrap().collect();
    assert_eq!(written.len(), 1);
    let log = fs::read_to_string(out.path().join("noref.log")).unwrap();
    assert!(
        log.trim_end().ends_with(&format!("Error: {error}")),
        "{log}"
    );
}

#[test]
fn the_run_is_logged() {
    let out = tempfile::tempdir().unwrap();
    let mut a = args(
        Some(canonical(INPUT_DIR)),
        None,
        out.path().join("noref"),
        2,
    );
    a.plot = Some("G7".to_string());
    run_with_args(a).unwrap();
    let log = fs::read_to_string(out.path().join("noref.log")).unwrap();
    let lines: Vec<&str> = log.lines().collect();
    assert!(lines[0].starts_with("command: "), "{log}");
    assert!(lines[1].starts_with("noref "), "{log}");
    assert!(lines[2].starts_with("parameters: ") && lines[2].contains("plot=G7"));
    for section in [
        " . pass 1: estimate k-mer occupancies",
        " . pass 2: extract sequences",
        "# analyse variant groups",
        "# build NJ tree of the raw alignment",
        "# cluster SNPs by isolate profile",
        "# mask concentrated SNPs",
        "# generate recombination file",
        "# build NJ tree of the filtered alignment",
    ] {
        assert!(lines.contains(&section), "{section} missing:\n{log}");
    }
    assert!(lines.last().unwrap().starts_with("total time: "), "{log}");
    // No progress bar.
    assert!(!log.contains('\r') && !log.contains('\x1b') && !log.contains("genome files"));
}

#[test]
fn output_is_identical_across_thread_counts_and_repeated_runs() {
    let temp = tempfile::tempdir().unwrap();
    let mut outputs = Vec::new();
    for (run, threads) in [1, 4, 4].into_iter().enumerate() {
        let dir = temp.path().join(format!("run_{run}"));
        fs::create_dir(&dir).unwrap();
        run_with_args(args(
            Some(canonical(INPUT_DIR)),
            None,
            dir.join("noref"),
            threads,
        ))
        .unwrap();
        outputs.push(common::read_outputs(&dir));
    }
    assert_eq!(
        outputs[0], outputs[1],
        "1-thread and 4-thread output differs"
    );
    assert_eq!(outputs[1], outputs[2], "repeated 4-thread output differs");
}

fn reverse_complement_fasta(text: &str) -> String {
    let mut out = String::new();
    let mut seq = String::new();
    let flush = |seq: &mut String, out: &mut String| {
        let rc: String = seq
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                'a' => 't',
                'c' => 'g',
                'g' => 'c',
                't' => 'a',
                other => other,
            })
            .collect();
        if !rc.is_empty() {
            out.push_str(&rc);
            out.push('\n');
        }
        seq.clear();
    };
    for line in text.lines() {
        if line.starts_with('>') {
            flush(&mut seq, &mut out);
            out.push_str(line);
            out.push('\n');
        } else {
            seq.push_str(line.trim());
        }
    }
    flush(&mut seq, &mut out);
    out
}

#[test]
fn alignment_is_independent_of_contig_orientation() {
    // Reverse-complement every contig of every genome: the alignment, missing-data
    // report and NJ tree must not change.
    let temp = tempfile::tempdir().unwrap();
    let flipped = temp.path().join("flipped");
    fs::create_dir(&flipped).unwrap();
    for entry in fs::read_dir(INPUT_DIR).unwrap() {
        let path = entry.unwrap().path();
        let text = common::read_fasta_text(&path);
        let name = path
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .trim_end_matches(".gz");
        fs::write(flipped.join(name), reverse_complement_fasta(&text)).unwrap();
    }
    let out = temp.path().join("out");
    fs::create_dir(&out).unwrap();
    run_with_args(args(Some(flipped), None, out.join("noref"), 1)).unwrap();
    common::compare_test_outputs(Path::new(EXPECTED), &out);
}

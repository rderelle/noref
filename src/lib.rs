//! Build phylogenomic nucleotide alignments directly from genome assemblies.
//!
//! `noref` finds exact canonical k-mer anchors shared across genomes, extracts the
//! variable regions between adjacent anchors, and concatenates the retained regions
//! into a nucleotide alignment. This avoids choosing a reference genome or defining
//! a fixed set of marker genes.
//!
//! This crate contains the command-line application's reusable entry point. Most users
//! should install and run the `noref` executable; callers embedding it can construct
//! [`Args`] and pass them to [`run_with_args`].
//!
//! # Quick start
//!
//! Put one genome assembly FASTA file per sample in a directory, then run:
//!
//! ```text
//! noref --input-directory genomes --output results/run1 --threads 4
//! ```
//!
//! The output directory (`results` above) must already exist. Run `noref --help` for
//! the complete and authoritative list of options and their defaults.
//!
//! # Input
//!
//! One or both of the following input sources may be supplied:
//!
//! - **Directory input** (`-i`, `--input-directory`) reads every `.fa`, `.fas`,
//!   `.fasta`, or `.fna` file in a directory (optionally gzipped). A filename without
//!   its FASTA and compression extensions becomes the sample name (for example,
//!   `isolate.1.fna.gz` becomes `isolate.1`).
//! - **Manifest input** (`-I`, `--input-file`) reads a headerless, two-column TSV
//!   containing `sample_name<TAB>genome_path` on each line. Relative paths are resolved
//!   from the manifest's directory.
//!
//! Contigs may be in any orientation: k-mers are canonical and variant groups are
//! stored in a strand-independent orientation. Every base other than A, C, G or T
//! (case-insensitive) is treated as `N`, i.e. missing data. A variant group found at
//! several positions in a genome (possible paralogs) is missing data for that genome.
//!
//! # Output
//!
//! `--output` is a filename prefix; when it names an existing directory or ends with a
//! path separator, the files are written in that directory with the default prefix
//! `noref`. With the default prefix, the application writes:
//!
//! - `noref_raw_alignment.fas` — the concatenated nucleotide alignment before
//!   recombination masking;
//! - `noref_filtered_alignment.fas` — the final alignment, with SNPs of recombinant
//!   origin replaced by `N` and the columns failing `-f` afterwards removed (see
//!   [Recombination](#recombination));
//! - `noref_raw_missing.tsv` and `noref_filtered_missing.tsv` — the percentage of
//!   missing (`-`) or ambiguous (`N`) sites per sample in each alignment; and
//! - `noref.log` — what is printed during the run (without the progress bars), after
//!   the command line and followed by the run time or by the error that stopped it.
//!
//! `--nj` additionally writes `noref_raw_NJ.tree` and `noref_filtered_NJ.tree`,
//! neighbor-joining trees in Newick format built from JC69-corrected distances on each
//! alignment. They are intended as a quick overview, not as a replacement for detailed
//! phylogenetic inference.
//!
//! `--skip-recombination` stops after the raw alignment (and its tree with `--nj`):
//! recombination is not detected, and no filtered alignment, tree or missing-data
//! report is written. It cannot be combined with `--plot`.
//!
//! `-p`/`--plot <isolate>` additionally writes `noref_recombination.gff`, a
//! recombination plot for [Phandango](https://github.com/jameshadfield/phandango) on
//! the genome of that isolate (see [Recombination](#recombination)).
//!
//! # Choosing parameters
//!
//! - `-f`/`--min-freq` (above 0.5 and up to 1, default 0.6): lower values retain
//!   sites present in fewer samples, producing a longer alignment with more missing
//!   data.
//! - `-k`/`--k`: anchor k-mer size (1–32, default 17).
//! - `-l`/`--length-middle`: maximum number of bases between two adjacent anchors.
//! - `-m`/`--mask`: minimum run of consecutive differences to the group consensus
//!   for a sequence to be masked (`N`; default 4); `0` disables this filter.
//! - `-w`/`--window`: window size (kb, default 5) used to detect SNPs of recombinant
//!   origin.
//! - `-a`/`--alpha`: p-value threshold for recombinant windows, per isolate and SNP
//!   cluster, before Bonferroni correction for the number of windows (default 0.05);
//!   lower values mask fewer SNPs.
//!
//! # Recombination
//!
//! SNPs are clustered by isolate profile (the set of isolates carrying each base). A
//! position with two bases is one SNP; a position with three or four bases gives one
//! SNP per base other than the most frequent one (isolates with that base vs the other
//! isolates with data), and is left out of the alignment if its two most frequent
//! bases are equally frequent. Every distinct profile without missing data is a
//! cluster; profiles with missing data, taken from the most to the least informative,
//! are split equally between the clusters they agree with, or start a new cluster if
//! they agree with none. Then, for each isolate, the SNPs of each cluster are placed
//! along its genome: windows of `-w` kb holding significantly more of a cluster's SNPs
//! than expected under uniform scattering (binomial test with fractional weights for
//! split SNPs, at `-a` after correction for the number of windows per genome; windows are
//! shortened at contig boundaries) are considered recombinant, and their SNPs are
//! replaced by `N` in that isolate (a position with several SNPs, or a SNP split
//! between clusters, is masked if any of them is). The scan is repeated with the
//! remaining SNPs until nothing new is found.
//! Masked columns are then checked again: those with fewer isolates with a base than
//! required by `-f`, or with a single base left, are removed from the alignment.
//!
//! With `--plot <isolate>`, these removed columns are placed on the genome of that
//! isolate (the reference), at the position of their variant group; columns in which
//! the reference has no base are left out. In each column, every other isolate in
//! which the SNP was masked and whose base before masking differs from the
//! reference's base is reported as carrying a recombinant SNP, and the SNPs of an
//! isolate less than `-w` kb apart on the same reference contig are merged into one
//! block. Blocks are written in the GFF format of Gubbins recombination predictions,
//! which Phandango displays next to a tree.
use clap::Parser;

mod dna;
mod group_extraction;
mod group_filtering;
mod group_sorting;
mod io;
mod output;
mod phylo;
mod plot;
mod proba_filter;
mod recombination;
mod run_log;
mod snp_clusters;

use run_log::{say, RunLog};

use std::path::{Path, PathBuf};

use group_filtering::AlignmentResult;
use io::IsolateInput;
use output::OutputPaths;

/// Prefix of the output files when `-o` is not given or names a directory.
const DEFAULT_PREFIX: &str = "noref";

#[derive(Parser, Debug)]
#[command(name = "noref", author, version, about)]
#[command(group = clap::ArgGroup::new("input_source").required(true).multiple(true).args(["input", "input_file"]))]
/// Parsed command-line options shared by the binary and integration tests.
pub struct Args {
    /// Directory containing genome FASTA files.
    #[arg(short, long = "input-directory", help_heading = "Main parameters")]
    pub input: Option<PathBuf>,

    /// TSV table with `isolate_name<TAB>path_to_fasta` rows.
    #[arg(short = 'I', long = "input-file", help_heading = "Main parameters")]
    pub input_file: Option<PathBuf>,

    /// Prefix for output files, or directory to write them in [o=noref].
    #[arg(short, long, help_heading = "Main parameters")]
    pub output: Option<PathBuf>,

    /// k-mer size [k=17].
    #[arg(
        short,
        long,
        default_value_t = 17,
        hide_default_value = true,
        help_heading = "Main parameters"
    )]
    pub k: usize,

    /// Minimum fraction of isolates present per alignment position [f=0.6].
    #[arg(
        short = 'f',
        long = "min-freq",
        default_value_t = 0.6,
        hide_default_value = true,
        help_heading = "Main parameters"
    )]
    pub min_freq: f64,

    /// Maximum length in bases between two adjacent shared anchors [l=70].
    #[arg(
        short = 'l',
        long = "length-middle",
        default_value_t = 70,
        hide_default_value = true,
        help_heading = "Main parameters"
    )]
    pub length_middle: usize,

    /// Consecutive differences required for masking [m=4]; 0 disables.
    #[arg(
        short = 'm',
        long = "mask",
        default_value_t = 4,
        hide_default_value = true,
        help_heading = "Main parameters"
    )]
    pub mask: usize,

    /// Window size (kb) for detecting concentrated SNPs of recombinant origin [w=5].
    #[arg(
        short = 'w',
        long = "window",
        default_value_t = 5,
        hide_default_value = true,
        help_heading = "Main parameters"
    )]
    pub window: u32,

    /// P-value threshold for recombinant windows, per isolate and SNP cluster, before Bonferroni correction [a=0.05].
    #[arg(
        short = 'a',
        long = "alpha",
        default_value_t = recombination::DEFAULT_ALPHA,
        hide_default_value = true,
        help_heading = "Main parameters"
    )]
    pub alpha: f64,

    /// Number of threads [t=1].
    #[arg(
        short = 't',
        long,
        default_value_t = 1,
        hide_default_value = true,
        help_heading = "Main parameters"
    )]
    pub threads: usize,

    /// Build neighbor-joining trees of the raw and filtered alignments.
    #[arg(long = "nj", help_heading = "Optional output")]
    pub nj: bool,

    /// Write a gff recombination file for Phandango using this isolate genome as reference.
    #[arg(
        short = 'p',
        long = "plot",
        value_name = "ISOLATE",
        help_heading = "Optional output"
    )]
    pub plot: Option<String>,

    /// Stop at the raw alignment: no recombination detection.
    #[arg(
        long = "skip-recombination",
        conflicts_with = "plot",
        help_heading = "Optional output"
    )]
    pub skip_recombination: bool,
}

fn print_startup_banner(args: &Args, paths: &OutputPaths, log: &RunLog) {
    let mut parameters = Vec::new();
    if let Some(input) = args.input.as_ref() {
        parameters.push(format!("input-directory={}", input.display()));
    }
    if let Some(input_file) = args.input_file.as_ref() {
        parameters.push(format!("input-file={}", input_file.display()));
    }
    parameters.push(format!("k={}", args.k));
    parameters.push(format!("min-freq={}", args.min_freq));
    parameters.push(format!("output={}", paths.prefix.display()));
    parameters.push(format!("length-middle={}", args.length_middle));
    parameters.push(format!("mask={}", args.mask));
    parameters.push(format!("window={}kb", args.window));
    parameters.push(format!("alpha={}", args.alpha));
    parameters.push(format!("threads={}", args.threads));
    if args.nj {
        parameters.push("nj=true".to_string());
    }
    if let Some(plot) = args.plot.as_ref() {
        parameters.push(format!("plot={plot}"));
    }
    if args.skip_recombination {
        parameters.push("skip-recombination=true".to_string());
    }
    say!(log, "noref {}", env!("CARGO_PKG_VERSION"));
    say!(log, "parameters: {}", parameters.join(" "));
}

fn validate_args(args: &Args) -> anyhow::Result<()> {
    anyhow::ensure!(
        (1..=dna::MAX_K).contains(&args.k),
        "-k/--k must be between 1 and {}",
        dna::MAX_K
    );
    anyhow::ensure!(
        (1..=u16::MAX as usize).contains(&args.length_middle),
        "-l/--length-middle must be between 1 and {}",
        u16::MAX
    );
    anyhow::ensure!(
        args.min_freq > 0.5 && args.min_freq <= 1.0,
        "-f/--min-freq must be greater than 0.5 and at most 1"
    );
    anyhow::ensure!(args.threads > 0, "-t/--threads must be at least 1");
    anyhow::ensure!(
        (1..=4_000_000).contains(&args.window),
        "-w/--window must be between 1 and 4,000,000 kb"
    );
    anyhow::ensure!(
        args.alpha > 0.0 && args.alpha <= 1.0,
        "-a/--alpha must be greater than 0 and at most 1"
    );
    anyhow::ensure!(
        !(args.skip_recombination && args.plot.is_some()),
        "--plot cannot be used with --skip-recombination"
    );
    Ok(())
}

/// Prefix of the output files: `-o` itself, or `<directory>/noref` when `-o` is an
/// existing directory or ends with a path separator.
fn output_prefix(output: Option<&Path>) -> PathBuf {
    match output {
        None => PathBuf::from(DEFAULT_PREFIX),
        Some(o) if o.is_dir() || o.to_string_lossy().ends_with(std::path::is_separator) => {
            o.join(DEFAULT_PREFIX)
        }
        Some(o) => o.to_path_buf(),
    }
}

/// Minimum number of isolates with a base for a position to be kept: `-f` times the
/// number of isolates, rounded up. The tolerance absorbs the binary representation of
/// `-f`, so that, for example, 0.6 × 25 gives 15 and not 16.
fn min_isolates(min_freq: f64, n: usize) -> usize {
    ((min_freq * n as f64) - 1e-9).ceil().max(1.0) as usize
}

/// Validate arguments, collect inputs, run the analysis, and write output files. What
/// is printed is also written to `<prefix>.log`, with the run time or the error that
/// stopped the run.
pub fn run_with_args(args: Args) -> anyhow::Result<()> {
    let paths = output::output_paths(&output_prefix(args.output.as_deref()));
    let log = RunLog::create(&paths.log)?;
    let result = run(&args, &paths, &log);
    let closed = log.finish(&result);
    result.and(closed)
}

fn run(args: &Args, paths: &OutputPaths, log: &RunLog) -> anyhow::Result<()> {
    print_startup_banner(args, paths, log);
    validate_args(args)?;
    let inputs = collect_inputs(args)?;
    let plot_reference = find_plot_reference(args, &inputs)?;
    let min_needed = min_isolates(args.min_freq, inputs.len());

    let mut alignment = build_alignment(args, &inputs, min_needed, log)?;
    write_raw_outputs(args, paths, &alignment, log)?;
    if args.skip_recombination {
        return Ok(());
    }
    let plot = mask_recombination(args, &mut alignment, min_needed, plot_reference, log)?;
    write_filtered_outputs(args, paths, &alignment, plot, log)
}

/// Isolates of the table (`-I`) and of the directory (`-i`), in that order. Names must
/// be unique, and at least two isolates are needed.
fn collect_inputs(args: &Args) -> anyhow::Result<Vec<IsolateInput>> {
    let mut inputs = Vec::new();
    if let Some(table) = args.input_file.as_ref() {
        inputs.extend(io::collect_isolates_from_table(table)?);
    }
    if let Some(dir) = args.input.as_ref() {
        inputs.extend(io::collect_isolates_from_dir(dir)?);
    }
    io::ensure_unique_names(&inputs)?;
    anyhow::ensure!(
        inputs.len() >= 2,
        "at least 2 genomes are needed, {} found",
        inputs.len()
    );
    Ok(inputs)
}

/// Row of the `--plot` reference (rows follow the input order).
fn find_plot_reference(args: &Args, inputs: &[IsolateInput]) -> anyhow::Result<Option<usize>> {
    let Some(name) = args.plot.as_deref() else {
        return Ok(None);
    };
    let index = inputs
        .iter()
        .position(|input| input.name == name)
        .ok_or_else(|| {
            anyhow::anyhow!("--plot: isolate '{name}' is not among the input genomes")
        })?;
    Ok(Some(index))
}

/// Extract, sort and filter the variant groups into the raw alignment.
fn build_alignment(
    args: &Args,
    inputs: &[IsolateInput],
    min_needed: usize,
    log: &RunLog,
) -> anyhow::Result<AlignmentResult> {
    say!(log, "# process input files");
    let raw_groups = group_extraction::extract_groups(
        inputs,
        args.k,
        min_needed,
        args.length_middle,
        args.threads,
        log,
    )?;
    say!(log, "# analyse variant groups");
    say!(log, " . raw variant groups: {}", raw_groups.groups.len());
    let sorted_groups = group_sorting::sort_and_deduplicate_groups(raw_groups)?;
    let alignment =
        group_filtering::filter_groups(sorted_groups, min_needed, args.mask, args.threads)?;
    say!(log, " . filtered variant groups: {}", alignment.n_groups);
    anyhow::ensure!(
        alignment.n_groups > 0,
        "no SNP found in at least {min_needed} of the {} isolates (see -f, -k and -l)",
        inputs.len()
    );
    Ok(alignment)
}

/// Write the raw alignment, its missing-data report and, with `--nj`, its tree. This is
/// done before recombination masking changes the rows, so that no copy of the
/// alignment is kept.
fn write_raw_outputs(
    args: &Args,
    paths: &OutputPaths,
    alignment: &AlignmentResult,
    log: &RunLog,
) -> anyhow::Result<()> {
    let names = &alignment.isolate_names;
    let raw = output::write_alignment(&paths.raw_alignment, names, &alignment.concat)?;
    output::write_missing(&paths.raw_missing, names, &alignment.concat)?;
    say!(
        log,
        " . raw alignment: {} SNPs and {} constant positions",
        raw.length - raw.constant,
        raw.constant
    );
    if args.nj {
        say!(log, "# build NJ tree of the raw alignment");
        output::write_nj_tree(&paths.raw_tree, names, &alignment.concat, args.threads)?;
    }
    Ok(())
}

/// Cluster the SNPs, mask those concentrated along each genome, and remove the columns
/// that then fail `-f` or are no longer polymorphic. With `--plot`, returns the
/// reference and the SNPs of the removed columns placed on its genome.
fn mask_recombination(
    args: &Args,
    alignment: &mut AlignmentResult,
    min_needed: usize,
    plot_reference: Option<usize>,
    log: &RunLog,
) -> anyhow::Result<Option<(usize, plot::MappedSnps)>> {
    say!(log, "# cluster SNPs by isolate profile");
    let clusters = snp_clusters::cluster_snps(&alignment.concat, args.threads)?;
    let s = &clusters.stats;
    say!(
        log,
        " . distinct SNP patterns: {} ({} without missing data)",
        s.patterns,
        s.complete_patterns
    );
    say!(log, " . clusters: {}", s.clusters);

    say!(log, "# mask concentrated SNPs");
    let layout = recombination::SnpLayout {
        snps: &clusters.snps,
        weights: &clusters.weights,
        n_clusters: clusters.stats.clusters,
        group_ends: &alignment.group_ends,
        positions: &alignment.positions,
        contig_ends: &alignment.contig_ends,
    };
    let r = recombination::mask_recombinant_snps(
        &mut alignment.concat,
        &layout,
        args.window * 1000,
        args.alpha,
        args.threads,
        plot_reference.is_some(),
    )?;
    say!(
        log,
        " . masked SNPs: {} of {} isolate-SNP cells with data ({:.2}%)",
        r.masked_cells,
        r.cells,
        100.0 * r.masked_cells as f64 / r.cells.max(1) as f64
    );
    say!(
        log,
        " . SNP columns masked in at least one isolate: {} of {}",
        r.masked_columns.len(),
        clusters.stats.polymorphic
    );
    say!(log, " . isolates with masked SNPs: {}", r.isolates);
    drop(clusters);

    // Genome position of each variant group in the --plot reference.
    let n = alignment.isolate_names.len();
    let reference_positions: Vec<u32> = plot_reference.map_or_else(Vec::new, |i| {
        alignment
            .positions
            .iter()
            .skip(i)
            .step_by(n)
            .copied()
            .collect()
    });
    alignment.positions = Vec::new();
    // Masked columns must still pass -f and still be polymorphic.
    let removed = recombination::failing_columns(&alignment.concat, &r.masked_columns, min_needed);
    // The SNPs of the removed columns are placed on the reference before removal.
    let plot = plot_reference.map(|index| {
        let reference = plot::Reference {
            index,
            group_ends: &alignment.group_ends,
            positions: &reference_positions,
        };
        let mapped =
            plot::map_removed_snps(&alignment.concat, &removed, &r.original_bases, &reference);
        (index, mapped)
    });
    drop(r);
    recombination::remove_columns(&mut alignment.concat, &removed);
    say!(
        log,
        " . SNPs removed: {} (below -f or no longer polymorphic)",
        removed.len()
    );
    Ok(plot)
}

/// Write the filtered alignment and its missing-data report, then the recombination
/// plot (`--plot`) and the tree of the filtered alignment (`--nj`).
fn write_filtered_outputs(
    args: &Args,
    paths: &OutputPaths,
    alignment: &AlignmentResult,
    plot: Option<(usize, plot::MappedSnps)>,
    log: &RunLog,
) -> anyhow::Result<()> {
    let names = &alignment.isolate_names;
    output::write_alignment(&paths.filtered_alignment, names, &alignment.concat)?;
    output::write_missing(&paths.filtered_missing, names, &alignment.concat)?;
    if let Some((reference, mapped)) = plot {
        write_recombination_plot(args, paths, alignment, reference, mapped, log)?;
    }
    if args.nj {
        say!(log, "# build NJ tree of the filtered alignment");
        output::write_nj_tree(&paths.filtered_tree, names, &alignment.concat, args.threads)?;
    }
    Ok(())
}

fn write_recombination_plot(
    args: &Args,
    paths: &OutputPaths,
    alignment: &AlignmentResult,
    reference: usize,
    mapped: plot::MappedSnps,
    log: &RunLog,
) -> anyhow::Result<()> {
    let names = &alignment.isolate_names;
    say!(log, "# generate recombination file");
    say!(
        log,
        " . removed SNPs placed on {} (--plot): {} of {}",
        names[reference],
        mapped.mapped,
        mapped.removed
    );
    let p = plot::write_plot(
        &paths.plot,
        names,
        reference,
        &alignment.contig_ends[reference],
        mapped,
        args.window * 1000,
    )?;
    say!(
        log,
        " . recombination plot: {} blocks in {} isolates ({} isolate-SNP cells)",
        p.blocks,
        p.isolates,
        p.snps
    );
    if p.blocks == 0 {
        say!(
            log,
            "   nothing to plot: Phandango cannot read a file without blocks"
        );
    }
    if p.unreadable_names > 0 {
        say!(
            log,
            "   warning: {} isolate names hold spaces or double quotes, which Phandango cannot match with tree leaves",
            p.unreadable_names
        );
    }
    Ok(())
}

#[cfg(test)]
mod argument_tests {
    use super::*;

    fn args(k: usize) -> Args {
        Args::parse_from(["noref", "-i", ".", "-k", &k.to_string()])
    }

    #[test]
    fn defaults() {
        let a = Args::parse_from(["noref", "-i", "."]);
        assert_eq!((a.k, a.length_middle, a.mask), (17, 70, 4));
        assert_eq!(a.min_freq, 0.6);
        assert!(validate_args(&a).is_ok());
    }

    #[test]
    fn skip_recombination_excludes_plot() {
        let a = Args::parse_from(["noref", "-i", ".", "--skip-recombination"]);
        assert!(a.skip_recombination && validate_args(&a).is_ok());
        let argv = ["noref", "-i", ".", "--skip-recombination", "-p", "G1"];
        assert!(Args::try_parse_from(argv).is_err());
        // Also checked for arguments built without the parser.
        let mut a = Args::parse_from(["noref", "-i", ".", "-p", "G1"]);
        a.skip_recombination = true;
        assert!(validate_args(&a).is_err());
    }

    #[test]
    fn alpha_bounds() {
        let alpha = |a: &str| Args::parse_from(["noref", "-i", ".", "-a", a]);
        assert_eq!(Args::parse_from(["noref", "-i", "."]).alpha, 0.05);
        assert!(validate_args(&alpha("0.001")).is_ok());
        assert!(validate_args(&alpha("1")).is_ok());
        assert!(validate_args(&alpha("0")).is_err());
        assert!(validate_args(&alpha("1.5")).is_err());
        let mut negative = alpha("0.05");
        negative.alpha = -0.1;
        assert!(validate_args(&negative).is_err());
    }

    #[test]
    fn min_freq_bounds() {
        let freq = |f: &str| Args::parse_from(["noref", "-i", ".", "-f", f]);
        assert!(validate_args(&freq("0.51")).is_ok());
        assert!(validate_args(&freq("1.0")).is_ok());
        assert!(validate_args(&freq("0.5")).is_err());
        assert!(validate_args(&freq("0.49")).is_err());
        assert!(validate_args(&freq("1.01")).is_err());
    }

    #[test]
    fn length_middle_bounds() {
        let l = |l: &str| Args::parse_from(["noref", "-i", ".", "-l", l]);
        assert!(validate_args(&l("1")).is_ok());
        assert!(validate_args(&l("65535")).is_ok());
        assert!(validate_args(&l("0")).is_err());
        assert!(validate_args(&l("65536")).is_err());
    }

    #[test]
    fn min_isolates_is_exact_for_decimal_fractions() {
        // 32-bit arithmetic gave 16 for 0.6 × 25, 61 for 0.6 × 100, 55 for 0.72 × 75.
        for (f, n, needed) in [
            (0.6, 25, 15),
            (0.6, 100, 60),
            (0.72, 75, 54),
            (0.54, 50, 27),
            (0.85, 20, 17),
            (0.85, 21, 18),
            (0.5, 3, 2),
            (0.5, 2, 1),
            (1.0, 7, 7),
        ] {
            assert_eq!(min_isolates(f, n), needed, "f={f} n={n}");
        }
        // Exact for every two-decimal value of -f and up to 3000 isolates.
        for hundredths in 50..=100u64 {
            // Division is correctly rounded, like parsing "0.60": same value as from the CLI.
            let f = hundredths as f64 / 100.0;
            for n in 1..=3000u64 {
                let exact = (hundredths * n).div_ceil(100) as usize;
                assert_eq!(min_isolates(f, n as usize), exact.max(1), "f={f} n={n}");
            }
        }
    }

    #[test]
    fn output_prefix_accepts_a_directory() {
        let dir = tempfile::tempdir().unwrap();
        assert_eq!(output_prefix(None), PathBuf::from("noref"));
        assert_eq!(output_prefix(Some(dir.path())), dir.path().join("noref"));
        let prefix = dir.path().join("run1");
        assert_eq!(output_prefix(Some(&prefix)), prefix);
        // A trailing separator names a directory, even one that does not exist yet.
        let new_dir = format!("{}{}", prefix.display(), std::path::MAIN_SEPARATOR);
        assert_eq!(
            output_prefix(Some(Path::new(&new_dir))),
            prefix.join("noref")
        );
    }

    #[test]
    fn k_bounds() {
        assert!(validate_args(&args(32)).is_ok());
        assert!(validate_args(&args(33)).is_err());
        assert!(validate_args(&args(0)).is_err());
    }

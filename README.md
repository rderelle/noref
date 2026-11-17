# noref

noref (v0.1) is a nucleotide version of [kamino](https://github.com/rderelle/kamino): it builds a nucleotide alignment in a reference-free, alignment-free manner from a set of genome assemblies. Not ‘better’ than traditional marker-based pipelines, but simpler and much faster to run.

## Installation
```bash
cd noref
cargo build --release           # binary: target/release/noref
```

## Usage
```bash
# standard analysis of genome assemblies using 4 CPUs
noref -i <input_dir> -t 4

# also get NJ trees (of the raw and filtered alignments)
noref -i <input_dir> -t 4 --nj

# also get a recombination plot for Phandango, on the genome of isolate 'ncbi_species_1'
noref -i <input_dir> -t 4 --nj -p ncbi_species_1

# stop at the raw alignment (no recombination detection), with its NJ tree
noref -i <input_dir> -t 4 --nj --skip-recombination

# get a shorter alignment with less missing data by raising the minimum sample occupancy per position (default 0.5)
noref -i <input_dir> -t 4 -f 0.85
```

### input
Input consists of genome assemblies in FASTA format (gzipped or not), with one file per sample. Files can be placed in a single directory, specified with the `-i` argument. Files are then recognised by their extension (.fa, .fas, .fasta, .fna; gzipped or not) and filenames minus extensions become isolate names in the final alignment.

Alternatively, the path of input files can be provided in a headerless tab-delimited file using the `-I` argument. The first column should contain the isolate names in the final alignment and the second column the path of the genome files.
```bash
species_soil	new_assemblies/contigs.fas.gz
ncbi_species_1	ncbi/GCA_00076868/genome.fna.gz
ncbi_species_2	ncbi/GCA_01092864/genome.fna.gz
```
Contigs can be in any orientation: k-mers (default k=17) are canonical and variant groups are stored in a strand-independent orientation. Any base other than A, C, G or T (IUPAC ambiguity codes, N...) is encoded as `N` and treated as missing data. A variant group found at several positions in a genome (possible paralogs) is treated as missing data for that genome.

### output
noref outputs two FASTA nucleotide alignments, each with a file containing the percentage of missing data per isolate (`noref_raw_missing.tsv` and `noref_filtered_missing.tsv`). The raw alignment (`noref_raw_alignment.fas`) contains all SNPs, before recombination masking. In the filtered alignment (`noref_filtered_alignment.fas`), SNPs found to be concentrated along a genome (likely recombination) are replaced by `N` in the isolates where they are concentrated (window size set with `-w`, in kb; 5 by default; p-value threshold set with `-a`, 0.05 by default, before Bonferroni correction for the number of windows per genome). Masked SNPs that then fail the `-f` threshold, or are no longer polymorphic, are removed from the filtered alignment. Positions with three or four bases are tested too, as one SNP per base other than the most frequent one; positions whose two most frequent bases are equally frequent are left out of both alignments. A log file (`noref.log`) records the command line, everything printed during the run (without the progress bars) and the total run time. By default, these output files use the prefix 'noref', but this can be changed with the `-o` argument; if `-o` is an existing directory (or ends with `/`), the files are written in it with the prefix 'noref'.

Additionally, Neighbour-joining trees (JC69 distances) of both alignments (`noref_raw_NJ.tree` and `noref_filtered_NJ.tree`) can be generated together with other output files by selecting the `--nj` optional argument.

With `--skip-recombination`, noref stops after the raw alignment (and its NJ tree with `--nj`): recombination is not detected, and no filtered alignment, tree or missing-data file is written. This option cannot be combined with `--plot`.

### recombination plot
With `-p`/`--plot <isolate>`, noref also writes `noref_recombination.gff`, a file to load into [Phandango](https://jameshadfield.github.io/phandango/) together with a tree (e.g. one of the `--nj` trees) to visualise recombination along the genome of the chosen isolate (the reference). It covers the SNP columns removed from the filtered alignment after masking (those failing the `-f` threshold or no longer polymorphic once recombinant SNPs are replaced by `N`):
- each column is placed at the position of its variant group in the reference genome (contigs placed end to end, in file order); columns in which the reference has no base cannot be placed;
- in each column, every other isolate in which the SNP was masked and whose base (before masking) differs from the base of the reference is reported as carrying a recombinant SNP;
- for each isolate, SNPs less than one window (`-w`) apart on the same reference contig are merged into one block (from the first to the last SNP; the number of SNPs is given as `snp_count`).

The file follows the GFF format of Gubbins recombination predictions (its second column must read `GUBBINS` for Phandango to display it; each block is marked `program="noref"`). Positions are those of the variant groups (middle of the region between the two anchors), so a SNP lies at most `-l`/2 bases from its reported position. As SNPs are compared with the reference, the reference's own row is always empty.

## Documentation
The algorithm is described in the manuscript mentioned below (for kamino, the amino-acid version from which noref derives); run `cargo doc --open` for the API documentation of this nucleotide version.

## Citation
> Romain Derelle, John A. Lees and Leonid Chindelevitch. 2026.
> kamino: fast proteome-wide variant calling for amino acid phylogenomics.
> [Molecular Biology and Evolution](https://academic.oup.com/mbe/article/43/10/msag250/8840831)

---
This codebase is provided under the MIT License. Some parts of the code were drafted using AI assistance.

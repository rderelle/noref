<h1 align="center">noref</h1>

noref ('no reference') is designed for within-species phylogenetic analyses (population genomics) from genome assemblies, in a fast and easy way.
It is only recommended for bacterial species and only accepts genome assemblies as input (not sequencing reads).

Reference-free variant callers, such as [SKA](https://github.com/bacpop/ska.rust), have to fall back on a reference genome to deal with recombination: SNPs are mapped onto a reference genome, and recombination tracts are then detected and removed by another tool (e.g. [ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML), [Gubbins](https://github.com/nickjcroucher/gubbins) or [Verticall](https://github.com/rrwick/Verticall)).
noref was developed to perform both steps, SNP calling and detection of recombination tracts, without a reference genome.
This reduces the biases introduced by a reference genome and greatly simplifies phylogenetic workflows: there is no need to choose a closely related reference genome, and one tool does everything.

## under the hood
- SNP calling: this is essentially the algorithm developed in [kamino](https://github.com/rderelle/kamino), but using nucleotides instead of amino acids. Also, the genomic positions in the original input assemblies of each extracted sequence from each variant group are saved in memory instead of protein names, allowing the next step to happen.
- detection of recombination: SNPs are first clustered by isolate profile. Then each input assembly is scanned for local enrichment of each SNP cluster independently, using a window scan and a binomial test (p-value threshold with Bonferroni correction): SNPs in significant windows are masked in the corresponding isolate.

## installation
```bash
git clone https://github.com/rderelle/noref.git
cd noref
cargo build --release           # binary: target/release/noref
```

## Usage
```bash
# standard analysis of genome assemblies using 4 CPUs
noref -i <input_dir> -t 4

# also get NJ trees (of the raw and filtered SNP alignments)
noref -i <input_dir> -t 4 --nj

# also get a recombination gff file for Phandango, on the genome of isolate 'ncbi_species_1'
noref -i <input_dir> -t 4 --plot ncbi_species_1
```

### some examples
[TBC]

### input
[TBC] bacterial genome assemblies in a directory with `-i`, or paths of input files provided in a headerless tab-delimited file using the `-I` argument [kamino](https://github.com/rderelle/kamino#input).

### output
[TBC] SNP alignments, NJ trees with `--nj`, and a recombination gff file for [Phandango](https://github.com/jameshadfield/phandango) with `-p`/`--plot`.
For the latter, the name of the isolate whose genome the recombinant SNPs will be mapped onto should be indicated.


## Documentation
[TBC]

---
This codebase is provided under the MIT License. The code was drafted using AI assistance.

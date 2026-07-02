[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/count_logan_tig_coverage/README.html)

# count_logan_tig_coverage

![Version](https://img.shields.io/badge/version-0.1.0-blue)

Compute **length-weighted mean and median coverage** from a Logan FASTA/FASTQ file.

Coverage values are read from the `ka:f:` tag present in Logan contig and unitig headers.

## Install

```bash
cargo install --git https://github.com/pierrepeterlongo/count_logan_tig_coverage
```

### From source

```bash
git clone https://github.com/pierrepeterlongo/count_logan_tig_coverage
cd count_logan_tig_coverage
cargo build --release
# binary: target/release/count_logan_tig_coverage
```

## Usage

```
count_logan_tig_coverage --in <FASTA/FASTQ file>
```

Supports plain, gzip, and zstd-compressed files (`.fa`, `.fa.gz`, `.fa.zst`, etc.).

### Options

| Option | Description |
|--------|-------------|
| `--in` | Input FASTA/FASTQ file (required) |
| `-t`, `--threads` | Number of threads (default: all available) |

### Example

```bash
count_logan_tig_coverage --in SRR1234567.contigs.fa.zst
```

```
Number of sequences:         42318
Total length:                89 421 304 bp
Mean coverage:               28.741
Median coverage (sequences): 21.500
Median coverage (bases):     24.000
```

## Output

| Field | Description |
|-------|-------------|
| Number of sequences | Total number of sequences with a `ka:f:` tag |
| Total length | Sum of all sequence lengths (bp) |
| Mean coverage | Length-weighted mean of `ka:f:` values |
| Median coverage (sequences) | Median over sequences (each sequence counts once) |
| Median coverage (bases) | Length-weighted median (coverage at which cumulative base count crosses 50%) |

## Context

This tool is designed to work alongside [logan_blaster](https://github.com/pierrepeterlongo/logan_blaster), which uses it to report coverage statistics of downloaded and recruited contigs/unitigs for each processed SRA accession.

## Author

- [Pierre Peterlongo](https://people.rennes.inria.fr/Pierre.Peterlongo/)

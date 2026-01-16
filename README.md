# CodonOpt — Codon Optimization Platform with RNAfold Integration

CodonOpt is a reproducible codon optimization platform for DNA sequences that supports custom codon usage tables, ORF-preserving optimization, GC-content constraints, motif and restriction site avoidance, Kleinbub or unbiased codon sampling, generation of multiple unique optimized sequences, RNA secondary structure evaluation via RNAfold, and a comprehensive TSV report of all computed metrics. The platform is fully containerized using Docker and can be run identically on a local machine, shared workstation, or HPC environment.

CodonOpt is designed for users who want fine-grained control over codon choice while ensuring biological validity, reproducibility, and downstream compatibility with cloning, expression, and RNA structure considerations.

---

## Features

- ORF-preserving codon optimization
- Custom codon usage tables (Excel-based)
- GC-content bounds
- Explicit codon exclusion
- Motif and restriction site avoidance
- Unbiased or Kleinbub-style codon sampling
- Multiple unique optimized sequences per input
- RNA secondary structure metrics via RNAfold
- Detailed TSV metrics report
- Fully containerized execution with Docker

---

## Repository Structure

```text
.
|-- codonopt/
|   |-- main.py
|   |-- optimizer/
|   |-- metrics/
|   |   |-- rnafold.py
|   |-- utils/
|-- Dockerfile
|-- README.md
`-- SelectedCodonTables.xlsx

```
The codonopt/ directory contains all application logic. The Dockerfile defines a fully reproducible runtime environment, including RNAfold. SelectedCodonTables.xlsx is an example codon usage table file.

---

## Codon Usage Tables (Required)

CodonOpt requires a codon usage table file to be provided at runtime. The file can live anywhere on your system, but when running via Docker it must be mounted into the container.

A recommended structure is:

```text
codon_tables/
└── SelectedCodonTables.xlsx
```
The path to this file is supplied using the --codon-xlsx argument.

---

## Codon Table Format

Codon usage tables are provided as an Excel file. Each sheet corresponds to a single codon usage table (for example an organism or expression system). The sheet name is selected using the --sheet argument.

Each sheet must contain the following columns:

AA: Amino acid (1-letter or 3-letter code)  
Codon: DNA codon (uppercase, e.g. GCT)  
Frequency: Relative codon usage (values do not need to sum to 1)

Example:

AA | Codon | Frequency  
A  | GCT   | 0.27  
A  | GCC   | 0.40  
A  | GCA   | 0.23  
A  | GCG   | 0.10  

CodonOpt normalizes frequencies internally as needed.

---

## Installation Using Docker (Recommended)

Docker is the recommended way to run CodonOpt, as it guarantees that RNAfold and all dependencies are installed correctly.

From the repository root, build the Docker image:

docker build -t codonopt:1.0 .

---

## Running CodonOpt Using Docker

When running via Docker, you must mount any input files (DNA sequences and codon tables) into the container.

Minimal example:

docker run --rm \
  -v $(pwd):/data \
  codonopt:1.0 \
  --dna ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG \
  --codon-xlsx /data/SelectedCodonTables.xlsx \
  --sheet ecoli \
  --out /data/output

Typical full pipeline example:

docker run --rm \
  -v $(pwd):/data \
  codonopt:1.0 \
  --dna input.fasta \
  --codon-xlsx /data/SelectedCodonTables.xlsx \
  --sheet ecoli \
  --gc-min 0.45 \
  --gc-max 0.65 \
  --avoid-codons CTA TAG \
  --avoid-motifs GAATTC AAGCTT \
  --optimization kleinbub \
  --n 10 \
  --out /data/output

All input and output paths must be inside the mounted directory (for example /data).

---

## Running CodonOpt Without Docker

CodonOpt can be run directly using Python, but this requires RNAfold to be installed and available on the system PATH.

python -m codonopt.main \
  --dna input.fasta \
  --codon-xlsx SelectedCodonTables.xlsx \
  --sheet ecoli \
  --out output

Docker execution is strongly recommended for reproducibility.

---

## Command-Line Arguments

--dna: Input DNA sequence (FASTA file or raw sequence)  
--codon-xlsx: Path to codon usage Excel file  
--sheet: Sheet name within the Excel file  
--gc-min: Minimum allowed GC content (0–1)  
--gc-max: Maximum allowed GC content (0–1)  
--avoid-codons: Codons to exclude  
--avoid-motifs: DNA motifs or restriction sites to avoid  
--optimization: unbiased or kleinbub  
--n: Number of optimized sequences to generate  
--out: Output directory  
--disable-rnafold: Disable RNAfold evaluation  

---

## Outputs

CodonOpt generates FASTA files for each optimized sequence and a TSV file containing metrics for all sequences.
```text
output/
├── optimized_001.fasta
├── optimized_002.fasta
└── metrics.tsv
```
The metrics TSV contains one row per sequence and includes sequence length, GC content, codon adaptation metrics, codon usage entropy, forbidden codon counts, motif violation counts, RNAfold minimum free energy (ΔG), RNAfold structure string, and ORF validation status.

---

## RNAfold Integration

RNAfold (ViennaRNA) is used to compute RNA secondary structure metrics, including minimum free energy and dot-bracket structures. RNAfold is installed automatically in the Docker image and is enabled by default. It can be disabled using the --disable-rnafold flag.

To verify RNAfold inside the Docker container:

docker run --rm --entrypoint micromamba codonopt:1.0 run -n base RNAfold --version

---

## Notes and Best Practices

- Input DNA must be in frame
- Start and stop codons are preserved
- ORFs are strictly conserved
- Codon tables must use DNA codons (not RNA codons)
- Docker execution is recommended for reproducibility

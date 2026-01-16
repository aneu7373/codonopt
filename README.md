**CodonOpt — Codon Optimization Platform with RNAfold Integration**

CodonOpt is a reproducible codon optimization platform for DNA sequences that supports custom codon usage tables, ORF-preserving optimization, GC-content constraints, motif and restriction site avoidance, Kleinbub or unbiased codon sampling, generation of multiple unique optimized sequences, RNA secondary structure evaluation via RNAfold, and a comprehensive TSV report of all computed metrics. The platform is fully containerized using Docker and can be run identically on a local machine, shared workstation, or HPC environment.

CodonOpt is designed for users who want fine-grained control over codon choice while ensuring biological validity, reproducibility, and downstream compatibility with cloning, expression, and RNA structure considerations.

**Features**

CodonOpt provides ORF-preserving codon optimization, supports custom codon usage tables supplied by the user, enforces GC-content bounds, allows explicit exclusion of specific codons, avoids user-defined motifs and restriction sites, supports unbiased or Kleinbub-style probabilistic codon sampling, generates multiple unique optimized sequences per input, computes RNA secondary structure metrics using RNAfold, outputs a detailed TSV report of all metrics, and runs reproducibly via Docker.

**Repository Structure**

The repository is organized as follows:

.
├── codonopt/
│   ├── main.py
│   ├── optimizer/
│   ├── metrics/
│   │   └── rnafold.py
│   └── utils/
├── Dockerfile
├── README.md
└── SelectedCodonTables.xlsx


The codonopt/ directory contains all application logic. The Dockerfile defines a fully reproducible runtime environment, including RNAfold. SelectedCodonTables.xlsx is an example codon usage table file.

**Codon Usage Tables (Required)**

CodonOpt requires a codon usage table file to be provided at runtime. The file can live anywhere on your system, but when running via Docker it must be mounted into the container. A common and recommended approach is to keep codon tables in a dedicated directory such as:

codon_tables/
└── SelectedCodonTables.xlsx


The path to this file is supplied using the --codon-xlsx argument.

**Codon Table Format**

Codon usage tables are provided as an Excel file. Each sheet corresponds to a single codon usage table, for example an organism or expression system. The sheet name is selected using the --sheet argument.

Each sheet must contain the following columns:

AA: Amino acid (1-letter or 3-letter code)

Codon: DNA codon (uppercase, e.g. GCT)

Frequency: Relative codon usage (values do not need to sum to 1)

An example row layout is:

AA | Codon | Frequency
A | GCT | 0.27
A | GCC | 0.40
A | GCA | 0.23
A | GCG | 0.10

CodonOpt will normalize frequencies internally as needed.

**Installation Using Docker (Recommended)**

Docker is the recommended way to run CodonOpt, as it guarantees that RNAfold and all dependencies are installed correctly.

From the repository root, build the Docker image:

docker build -t codonopt:1.0 .

**Running CodonOpt Using Docker**

When running via Docker, you must mount any input files (DNA sequences and codon tables) into the container. The examples below assume you are running from a directory that contains your input files.

Minimal Example
docker run --rm \
  -v $(pwd):/data \
  codonopt:1.0 \
  --dna ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG \
  --codon-xlsx /data/SelectedCodonTables.xlsx \
  --sheet ecoli \
  --out /data/output

Typical Full Pipeline Example
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


In all Docker runs, input and output paths must be inside the mounted directory (for example /data).

**Running CodonOpt Without Docker**

CodonOpt can also be run directly using Python, although this requires RNAfold to be installed and available on the system PATH.

After installing dependencies, run:

python -m codonopt.main \
  --dna input.fasta \
  --codon-xlsx SelectedCodonTables.xlsx \
  --sheet ecoli \
  --out output


*Using Docker is strongly recommended for reproducibility.*

**Command-Line Arguments**

The main command-line arguments are:

--dna: Input DNA sequence, either as a FASTA file or a raw sequence string

--codon-xlsx: Path to the codon usage Excel file

--sheet: Name of the sheet within the Excel file to use

--gc-min: Minimum allowed GC content (0–1)

--gc-max: Maximum allowed GC content (0–1)

--avoid-codons: One or more codons to exclude from optimization

--avoid-motifs: One or more DNA motifs or restriction sites to avoid

--optimization: Codon sampling strategy (unbiased or kleinbub)

--n: Number of unique optimized sequences to generate

--out: Output directory

--disable-rnafold: Disable RNAfold secondary structure evaluation

**Outputs**

CodonOpt produces FASTA files for each optimized sequence and a TSV file containing metrics for all generated sequences.

The output directory structure is:

output/
├── optimized_001.fasta
├── optimized_002.fasta
└── metrics.tsv


The metrics TSV contains one row per optimized sequence and includes sequence length, GC content, codon adaptation metrics, codon usage entropy, forbidden codon counts, motif violation counts, RNAfold minimum free energy (ΔG), RNAfold structure string, and ORF validation status.

**RNAfold Integration**

RNAfold from the ViennaRNA package is used to compute RNA secondary structure metrics, including minimum free energy and dot-bracket structure strings. RNAfold is installed automatically in the Docker image and is enabled by default. It can be disabled using the --disable-rnafold flag.

To verify that RNAfold is available inside the Docker container, run:

docker run --rm --entrypoint micromamba codonopt:1.0 \
  run -n base RNAfold --version

**Notes and Best Practices**

Input DNA sequences must be in frame. Start and stop codons are preserved, and ORFs are strictly conserved during optimization. Codon usage tables must use DNA codons rather than RNA codons. For consistent and reproducible results across systems, Docker execution is recommended.

CodonOpt — Codon Optimization Platform with RNAfold Integration

CodonOpt is a reproducible codon optimization platform for DNA sequences that supports:

Custom codon usage tables

ORF-preserving optimization

GC-content constraints

Motif and restriction site avoidance

Kleinbub or unbiased codon sampling

Multiple unique optimized sequences per input

RNA secondary structure evaluation via RNAfold

Fully containerized execution via Docker

Features

✅ ORF-preserving codon optimization

✅ Custom codon usage tables (Excel / TSV-based)

✅ GC-content bounds

✅ Codon exclusion lists

✅ Motif & restriction site avoidance

✅ Multiple optimized outputs per sequence

✅ RNAfold secondary structure metrics

✅ TSV report of all metrics

✅ Dockerized & HPC-friendly

Repository Structure
.
├── codonopt/
│   ├── main.py                 # CLI entrypoint
│   ├── optimizer/
│   ├── metrics/
│   │   └── rnafold.py
│   └── utils/
├── Dockerfile
├── README.md
└── SelectedCodonTables.xlsx    # Example codon usage tables

Codon Usage Tables (IMPORTANT)
Where to put your codon table

You must provide a codon usage table file when running CodonOpt.

The file can live anywhere on your system

When using Docker, it must be mounted into the container

Recommended:

codon_tables/
└── SelectedCodonTables.xlsx

Codon table format (Excel)

Each sheet corresponds to a codon usage table (e.g. organism or expression system).

Required columns:

Column	Description
AA	Amino acid (1-letter or 3-letter)
Codon	DNA codon (uppercase, e.g. GCT)
Frequency	Relative usage (does not need to sum to 1)

Example:

AA	Codon	Frequency
A	GCT	0.27
A	GCC	0.40
A	GCA	0.23
A	GCG	0.10

You select the sheet using --sheet.

Installation (Docker — Recommended)
Build the image

From the repository root:

docker build -t codonopt:1.0 .

Running CodonOpt (Docker)
Minimal example
docker run --rm \
  -v $(pwd):/data \
  codonopt:1.0 \
  --dna ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG \
  --codon-xlsx /data/SelectedCodonTables.xlsx \
  --sheet ecoli \
  --out /data/output

Typical full pipeline example
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

Running CodonOpt (Native Python)

Requires Python ≥3.10 and RNAfold installed and on $PATH

pip install -r requirements.txt
python -m codonopt.main \
  --dna input.fasta \
  --codon-xlsx SelectedCodonTables.xlsx \
  --sheet ecoli \
  --out output

Command-line Arguments
Argument	Description
--dna	Input DNA (FASTA file or raw sequence)
--codon-xlsx	Path to codon usage Excel file
--sheet	Sheet name to use
--gc-min	Minimum GC content (0–1)
--gc-max	Maximum GC content (0–1)
--avoid-codons	Codons to exclude
--avoid-motifs	DNA motifs / restriction sites to avoid
--optimization	unbiased or kleinbub
--n	Number of unique optimized sequences
--out	Output directory
--disable-rnafold	Skip RNAfold evaluation
Outputs
FASTA files
output/
├── optimized_001.fasta
├── optimized_002.fasta
└── ...

Metrics TSV (one row per sequence)
output/
└── metrics.tsv


Included metrics:

Sequence length

GC content

Codon Adaptation Index (CAI)

Codon usage entropy

Forbidden codon count

Motif violation count

RNAfold MFE (ΔG)

RNAfold structure string

ORF validation status

RNAfold Integration

RNAfold (ViennaRNA) is used to compute:

Minimum free energy (ΔG)

Secondary structure

RNAfold is:

Automatically enabled if installed

Installed in the Docker image

Optional via --disable-rnafold

To confirm RNAfold in Docker:

docker run --rm --entrypoint micromamba codonopt:1.0 \
  run -n base RNAfold --version

Notes & Best Practices

Input DNA must be in-frame

Stop codons are preserved

ORFs are strictly conserved

Codon usage tables should use DNA codons (not RNA)

Use Docker for full reproducibility

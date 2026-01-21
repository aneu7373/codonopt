# Codon the Barbarian

**Codon the Barbarian** is a flexible, batch-capable codon optimization and back-translation pipeline designed for real-world sequence design workflows. It supports DNA or protein FASTA inputs, organism-specific codon tables (via Excel), tunable optimization strategies, and detailed diagnostic output.

The pipeline is designed to **maximize yield under complex constraints** rather than failing prematurely, while still making failures explicit and debuggable.

---

## Features

- Accepts **DNA CDS FASTA or protein FASTA** inputs
- **Batch mode** via CSV or TSV for large-scale design
- Codon tables loaded from **Excel (.xlsx) files**
  - Optional sheet selection per batch row
- Two optimization modes:
  - **kleinbub (default)**: bounded backtracking, high yield
  - strict: fast greedy optimization
- Rich constraint support:
  - GC content bounds
  - Codons to avoid
  - Motifs to avoid
  - Maximum homopolymer length
- Yield tuning:
  - Automatic retries for failed designs
  - Adjustable search effort
- Deterministic reproducibility via seeds
- Outputs:
  - Single FASTA file with all successful designs
  - TSV report with metrics and explicit failure reasons

---

## Installation

### Option 1: Docker (recommended)

Docker is the easiest way to run codonopt without worrying about dependencies.

    git clone https://github.com/aneu7373/codonopt.git
    cd codonopt
    docker build -t codonopt .

---

### Option 2: Local installation (advanced)

Requires Python 3.9 or newer.

    git clone https://github.com/aneu7373/codonopt.git
    cd codonopt
    pip install -r requirements.txt

---

## Quick Start (Single Sequence)

    docker run --rm \
      -v $(pwd)/data:/data \
      codonopt \
      --sequence /data/example.fasta \
      --codon-table /data/Ecoli.xlsx \
      --out /data/results

Outputs:

- results/optimized_sequences.fasta  
- results/metrics.tsv  

---

## Batch Mode (Recommended)

Batch mode allows you to specify **one row per design job** in a CSV or TSV file.

### Minimal batch CSV

    sequence,codon_table
    /data/input.fasta,/data/Ecoli.xlsx

### Full batch CSV (all supported options)

    sequence,codon_table,codon_table_sheet,optimization_mode,avoid_codons,avoid_motifs,gc_min,gc_max,seed,n,max_tries_per_replicate,kleinbub_search_limit,backtrack_window
    /data/input.fasta,/data/Ecoli.xlsx,Ecoli1,kleinbub,CTA|TAG,GAATTC|AAGCTT,0.4,0.6,42,10,50,400000,25

Run batch mode:

    docker run --rm \
      -v $(pwd)/data:/data \
      codonopt \
      --batch-table /data/batch.csv \
      --out /data/results \
      --verbose

---

## Codon Tables

Codon tables are provided as **Excel (.xlsx) files**.

### Supported formats

#### Codon usage table (recommended)

| AA | Codon | Fraction |
|----|-------|----------|
| A  | GCT   | 0.18     |
| A  | GCC   | 0.27     |
| …  | …     | …        |

Column names are flexible and may include:  
AA, Amino Acid, Codon, Frequency, Fraction, Usage, or RSCU.

#### Wide format

| Amino Acid | Codons                  |
|------------|-------------------------|
| A          | GCT,GCC,GCA,GCG         |
| C          | TGT,TGC                 |
| *          | TAA,TAG,TGA             |

---

## Multiple Sheets

If the Excel file has multiple sheets:

- Specify the desired sheet via:
  - codon_table_sheet (batch mode)
  - --codon-table-sheet (single mode)
- If omitted, the **first non-empty sheet** is used automatically.

---

## Optimization Modes

### kleinbub (default, recommended)

- Bounded backtracking search
- Tries harder before declaring failure
- Best for:
  - Tight constraints
  - Long proteins
  - High-yield design

### strict

- Greedy, per-position retries
- Faster, but may falsely report “impossible”
- Useful for debugging

Set per batch row:

    optimization_mode
    strict

---

## Yield Tuning (Getting N/N Outputs)

By default, codonopt attempts to generate n sequences but does not guarantee success if constraints are tight.

Tune yield using:

| Parameter | Meaning | Typical Values |
|----------|--------|----------------|
| max_tries_per_replicate | Independent restarts per replicate | 25–100 |
| kleinbub_search_limit | Total search steps per attempt | 200k–800k |
| backtrack_window | Backtracking depth | 10–50 |

Example:

    n,max_tries_per_replicate,kleinbub_search_limit,backtrack_window
    10,50,500000,30

---

## Outputs

### FASTA

optimized_sequences.fasta

- Contains only successful designs
- IDs encode job and replicate:

    originalID|job0001|rep003

### TSV

metrics.tsv

Includes successful and failed attempts.

Key columns include:

- optimization_mode
- attempts_used
- failure_reason
- gc_content
- max_homopolymer
- All constraint and tuning parameters

---

## Reproducibility

Set a seed for deterministic behavior:

- Same inputs + same seed → same outputs
- Replicates vary deterministically

---

## Limits and Safety

- Maximum 1000 input FASTA records per run
- Bounded search prevents infinite loops
- No silent failures

---

## When Constraints Are Truly Impossible

codonopt distinguishes between:

- Search exhaustion (fixable by tuning)
- True incompatibility (intrinsic to the protein)

Inspect failure_reason in metrics.tsv.

---

## Development

    git clone https://github.com/aneu7373/codonopt.git
    cd codonopt

Key files:

- codonopt/main.py
- codonopt/core/optimizer.py
- codonopt/io/codon_tables.py

---

## License

MIT License (see LICENSE file).

---

## Citation

If you use this tool academically, please cite:

    codonopt — https://github.com/aneu7373/codonopt

---

Happy designing 🧬

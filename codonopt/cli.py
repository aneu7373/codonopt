
#!/usr/bin/env python3
"""
Codon optimization CLI with FASTA support.

Features:
- Accepts inline DNA (--dna) or a FASTA file (--fasta)
- Batch processing for multi-record FASTA
- GC content bounds: --gc-min / --gc-max
- Avoid lists: --avoid-codons / --avoid-motifs
- Restriction site scanning: --restriction-sites
- Disable RNAfold scoring: --disable-RNAfold
- Reproducibility: --seed
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Tuple

from Bio import SeqIO

from codonopt.io.codon_tables import load_codon_table_xlsx
from codonopt.core.orf import translate_dna
from codonopt.core.optimizer import optimize
from codonopt.io.output import write_outputs


def _parse_csv_list(s: str) -> List[str]:
    return [tok.strip().upper() for tok in s.split(",") if tok.strip()] if s else []


def _validate_gc_bounds(gc_min: float, gc_max: float) -> None:
    if not (0.0 <= gc_min <= 1.0):
        raise ValueError(f"--gc-min must be in [0,1], got {gc_min}")
    if not (0.0 <= gc_max <= 1.0):
        raise ValueError(f"--gc-max must be in [0,1], got {gc_max}")
    if gc_min > gc_max:
        raise ValueError(f"--gc-min ({gc_min}) cannot exceed --gc-max ({gc_max})")


def _validate_codons(codons: List[str]) -> None:
    valid = {"A", "C", "G", "T"}
    for c in codons:
        uc = c.upper()
        if len(uc) != 3 or any(ch not in valid for ch in uc):
            raise ValueError(f"Invalid codon '{c}'. Must be 3 letters from A/C/G/T.")


def _validate_motifs(motifs: List[str]) -> None:
    valid = {"A", "C", "G", "T"}
    for m in motifs:
        um = m.upper()
        if not um or any(ch not in valid for ch in um):
            raise ValueError(f"Invalid motif '{m}'. Motifs must contain only A/C/G/T.")


def _load_fasta_records(path: str) -> List[Tuple[str, str]]:
    """Return a list of (record_id, dna_seq) from a FASTA file."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"FASTA file not found: {path}")
    records = []
    for rec in SeqIO.parse(path, "fasta"):
        seq = str(rec.seq).upper().replace("U", "T")
        records.append((rec.id or "record", seq))
    if not records:
        raise ValueError(f"No sequences found in FASTA: {path}")
    return records


def main() -> None:
    p = argparse.ArgumentParser(
        description="Codon optimization with FASTA support, GC bounds, avoid codons/motifs, restriction sites, and optional RNAfold disabling."
    )

    # Mutually exclusive input sources
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--dna", help="Inline DNA sequence (A/C/G/T).")
    g.add_argument("--fasta", help="Path to a FASTA file with one or more records.")

    # Optional name for single-sequence runs
    p.add_argument("--name", default="sequence_1",
                   help="Name/ID for the single-sequence (--dna) case. Ignored with --fasta.")

    # Codon usage input
    p.add_argument("--codon-xlsx", required=True, help="Path to codon usage Excel file.")
    p.add_argument("--sheet", required=True, help="Sheet name in the Excel file.")

    # Optimization knobs
    p.add_argument("--n", type=int, default=10, help="Number of optimized sequences to output (per input record).")
    p.add_argument("--oversample", type=int, default=10, help="Oversampling factor (pool size = n * oversample).")

    # Constraints & toggles
    p.add_argument("--gc-min", type=float, default=None, help="Minimum GC fraction [0..1].")
    p.add_argument("--gc-max", type=float, default=None, help="Maximum GC fraction [0..1].")
    p.add_argument("--avoid-codons", type=str, default="", help="Comma-separated codons to avoid (e.g., TAG,TGA,ATA).")
    p.add_argument("--avoid-motifs", type=str, default="", help="Comma-separated DNA motifs to avoid (e.g., GGCACC,TTATAA).")
    p.add_argument("--restriction-sites", type=str, default="", help="Comma-separated restriction sites (e.g., GAATTC,GGATCC).")
    p.add_argument("--disable-RNAfold", action="store_true", help="Disable RNA folding metrics.")
    p.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")

    # Outputs
    p.add_argument("--out-fasta", default="output.fasta", help="Output FASTA path for optimized sequences.")
    p.add_argument("--out-report", default="report.tsv", help="Output TSV metrics path.")

    args = p.parse_args()

    # Validate GC bounds (both or none)
    if (args.gc_min is None) ^ (args.gc_max is None):
        print("ERROR: Provide both --gc-min and --gc-max together.", file=sys.stderr)
        sys.exit(2)
    if args.gc_min is not None:
        _validate_gc_bounds(args.gc_min, args.gc_max)  # type: ignore[arg-type]

    # Parse/validate lists
    avoid_codons = _parse_csv_list(args.avoid_codons)
    avoid_motifs = _parse_csv_list(args.avoid_motifs)
    restriction_sites = _parse_csv_list(args.restriction_sites)
    if avoid_codons:
        _validate_codons(avoid_codons)
    if avoid_motifs:
        _validate_motifs(avoid_motifs)

    # Seed (optional)
    if args.seed is not None:
        import random
        random.seed(args.seed)

    # Load codon usage/scoring table
    codon_table = load_codon_table_xlsx(args.codon_xlsx, args.sheet)

    # Scoring weights (RNAfold weights 0.0 if disabled)
    weights = {
        "gc": 1.0,
        "gc_target": 0.5,
        "codon": 1.0,
        "entropy": 0.5,
        "rna_5p": 0.0 if args.disable_RNAfold else 0.2,
        "rna_full": 0.0 if args.disable_RNAfold else 0.1,
        "restriction": 5.0,
        "motif": 2.0,
        "forbidden": 3.0,  # penalize unavoidable forbidden codons
        "gc_min": args.gc_min,
        "gc_max": args.gc_max,
    }

    # Prepare inputs (single sequence vs FASTA batch)
    if args.fasta:
        inputs = _load_fasta_records(args.fasta)  # List[(id, seq)]
    else:
        seq = args.dna.upper().replace("U", "T")  # tolerate RNA input by converting U->T
        inputs = [(args.name, seq)]

    # Run optimization for each input, collect and tag results
    all_results = []
    for rec_id, dna_seq in inputs:
        aa = translate_dna(dna_seq)
        results = optimize(
            dna=dna_seq,
            aa_seq=aa,
            codon_table=codon_table,
            forbidden_codons=set(avoid_codons),
            restriction_sites=restriction_sites,
            motifs=avoid_motifs,
            n_out=args.n,
            oversample=args.oversample,
            weights=weights,
        )
        # Tag record id so downstream reporting includes the source
        for seq, metrics in results:
            metrics = dict(metrics)
            metrics["id"] = rec_id
            all_results.append((seq, metrics))

    # Persist combined outputs (multi-record friendly)
    write_outputs(all_results, args.out_fasta, args.out_report)


if __name__ == "__main__":
    main()

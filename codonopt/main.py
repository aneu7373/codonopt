
import argparse
import sys
from typing import List, Set

from codonopt.io.codon_tables import load_codon_table_xlsx
from codonopt.core.orf import translate_dna
from codonopt.core.optimizer import optimize
from codonopt.io.output import write_outputs


def _parse_csv_list(s: str) -> List[str]:
    """Parse a comma-separated list string into a list of uppercase tokens."""
    if not s:
        return []
    return [tok.strip().upper() for tok in s.split(",") if tok.strip()]


def _validate_gc_bounds(gc_min: float, gc_max: float) -> None:
    """Validate that GC bounds are in [0,1] and min <= max."""
    if not (0.0 <= gc_min <= 1.0):
        raise ValueError(f"--gc-min must be between 0 and 1 (got {gc_min})")
    if not (0.0 <= gc_max <= 1.0):
        raise ValueError(f"--gc-max must be between 0 and 1 (got {gc_max})")
    if gc_min > gc_max:
        raise ValueError(f"--gc-min ({gc_min}) cannot be greater than --gc-max ({gc_max})")


def _validate_codons(codons: List[str]) -> None:
    """Ensure each codon is a length-3 string of A/C/G/T."""
    valid = {"A", "C", "G", "T"}
    for c in codons:
        if len(c) != 3 or any(ch not in valid for ch in c):
            raise ValueError(
                f"Invalid codon '{c}'. Codons must be 3 characters from A/C/G/T."
            )


def _validate_motifs(motifs: List[str]) -> None:
    """Ensure each motif is a non-empty string of A/C/G/T."""
    valid = {"A", "C", "G", "T"}
    for m in motifs:
        if len(m) == 0 or any(ch not in valid for ch in m):
            raise ValueError(
                f"Invalid motif '{m}'. Motifs must be non-empty and contain only A/C/G/T."
            )


def main():
    p = argparse.ArgumentParser(
        description="Codon optimization with GC bounds, avoid-codons/motifs, and optional RNAfold disabling."
    )
    p.add_argument("--dna", required=True, help="Input DNA sequence (A/C/G/T).")
    p.add_argument("--codon-xlsx", required=True, help="Path to codon usage Excel file.")
    p.add_argument("--sheet", required=True, help="Sheet name in the Excel file.")
    p.add_argument("--n", type=int, default=10, help="Number of optimized sequences to output.")

    # New: GC bounds (fractions in [0,1])
    p.add_argument("--gc-min", type=float, default=None,
                   help="Minimum GC fraction (0..1). If provided, must be <= --gc-max.")
    p.add_argument("--gc-max", type=float, default=None,
                   help="Maximum GC fraction (0..1). If provided, must be >= --gc-min.")

    # New: avoid lists
    p.add_argument("--avoid-codons", type=str, default="",
                   help="Comma-separated codons to avoid (e.g., 'TAG,TGA,ATA').")
    p.add_argument("--avoid-motifs", type=str, default="",
                   help="Comma-separated DNA motifs to avoid (e.g., 'GGCACC,TTATAA').")

    # New: disable RNAfold-derived scoring/penalties
    p.add_argument("--disable-RNAfold", action="store_true",
                   help="Disable RNA folding penalties (set RNA structure weights to 0).")

    args = p.parse_args()

    # Translate to amino acids (ensures synonymous-only changes downstream)
    aa = translate_dna(args.dna)

    # Load codon usage table from Excel
    codon_table = load_codon_table_xlsx(args.codon_xlsx, args.sheet)

    # Parse and validate avoid lists
    avoid_codons_list = _parse_csv_list(args.avoid_codons)
    avoid_motifs_list = _parse_csv_list(args.avoid_motifs)
    if avoid_codons_list:
        _validate_codons(avoid_codons_list)
    if avoid_motifs_list:
        _validate_motifs(avoid_motifs_list)

    # Base weights (defaults)
    weights = {
        "gc": 1.0,           # GC scoring weight
        "gc_target": 0.5,    # Default GC target (used if bounds not specified)
        "codon": 1.0,        # Codon usage weight
        "entropy": 0.5,      # Diversity / anti-repeats
        "rna_5p": 0.2,       # 5' structure penalty weight
        "rna_full": 0.1,     # global structure penalty weight
        "restriction": 5.0,  # Restriction-site penalty weight (if sites provided)
        "motif": 2.0,        # Motif penalty weight
    }

    # Apply RNAfold disabling by zeroing RNA weights
    if args.disable_RNAfold:
        weights["rna_5p"] = 0.0
        weights["rna_full"] = 0.0

    # Apply GC bounds if provided
    # - We set 'gc_min' and 'gc_max' for optimizers that support explicit bounds.
    # - We also set 'gc_target' to the midpoint to guide scoring-based approaches.
    gc_min = args.gc_min
    gc_max = args.gc_max
    if gc_min is not None or gc_max is not None:
        if gc_min is None or gc_max is None:
            print("ERROR: --gc-min and --gc-max must be provided together.", file=sys.stderr)
            sys.exit(2)
        _validate_gc_bounds(gc_min, gc_max)
        weights["gc_min"] = gc_min
        weights["gc_max"] = gc_max
        weights["gc_target"] = 0.5 * (gc_min + gc_max)

    # Build inputs for optimize()
    forbidden_codons: Set[str] = set(avoid_codons_list)
    motifs: List[str] = list(avoid_motifs_list)

    # Run the optimization
    results = optimize(
        dna=args.dna,
        aa_seq=aa,
        codon_table=codon_table,
        forbidden_codons=forbidden_codons,
        restriction_sites=[],    # you can add a CLI for these later if needed
        motifs=motifs,
        n_out=args.n,
        oversample=10,           # keep oversampling to select best within constraints
        weights=weights,
    )

    # Write outputs
    write_outputs(results, "output.fasta", "report.tsv")


if __name__ == "__main__":
    main()

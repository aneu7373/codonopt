import argparse
from codonopt.io.codon_tables import load_codon_table_xlsx
from codonopt.core.orf import translate_dna
from codonopt.core.optimizer import optimize
from codonopt.io.output import write_outputs

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dna", required=True)
    p.add_argument("--codon-xlsx", required=True)
    p.add_argument("--sheet", required=True)
    p.add_argument("--n", type=int, default=10)
    args = p.parse_args()

    aa = translate_dna(args.dna)
    codon_table = load_codon_table_xlsx(args.codon_xlsx, args.sheet)

    weights = {
        "gc": 1.0,
        "gc_target": 0.5,
        "codon": 1.0,
        "entropy": 0.5,
        "rna_5p": 0.2,
        "rna_full": 0.1,
        "restriction": 5.0,
        "motif": 2.0,
    }

    results = optimize(
        dna=args.dna,
        aa_seq=aa,
        codon_table=codon_table,
        forbidden_codons=set(),
        restriction_sites=[],
        motifs=[],
        n_out=args.n,
        oversample=10,
        weights=weights,
    )

    write_outputs(results, "output.fasta", "report.tsv")

if __name__ == "__main__":
    main()

from codonopt.core.orf import validate_orf
from codonopt.metrics.basic import gc_content, shannon_entropy
from codonopt.metrics.motifs import count_hits
from codonopt.metrics.rnafold import rnafold_mfe
from codonopt.core.scoring import score_candidate

def optimize(
    dna,
    aa_seq,
    codon_table,
    forbidden_codons,
    restriction_sites,
    motifs,
    n_out,
    oversample,
    weights,
):
    candidates = []

    target = n_out * oversample
    while len(candidates) < target:
        seq, forbidden_hits, codon_freqs = generate_candidate(
            aa_seq, codon_table, forbidden_codons
        )

        if forbidden_hits > 0:
            continue
        if not validate_orf(seq, aa_seq):
            continue

        metrics = {}
        metrics["gc"] = gc_content(seq)
        metrics["mean_codon_freq"] = sum(codon_freqs) / len(codon_freqs)
        metrics["entropy"] = shannon_entropy(codon_freqs)
        metrics["restriction_hits"] = count_hits(seq, restriction_sites)
        metrics["motif_hits"] = count_hits(seq, motifs)
        metrics["rna_5p"] = rnafold_mfe(seq[:150])
        metrics["rna_full"] = rnafold_mfe(seq)

        metrics["score"] = score_candidate(metrics, weights)
        candidates.append((seq, metrics))

    candidates.sort(key=lambda x: x[1]["score"], reverse=True)
    return candidates[:n_out]

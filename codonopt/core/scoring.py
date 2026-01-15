def score_candidate(metrics, weights):
    score = 0.0

    score -= weights["gc"] * abs(metrics["gc"] - weights["gc_target"])
    score += weights["codon"] * metrics["mean_codon_freq"]
    score += weights["entropy"] * metrics["entropy"]

    score += weights["rna_5p"] * metrics["rna_5p"]
    score += weights["rna_full"] * metrics["rna_full"]

    score -= weights["restriction"] * metrics["restriction_hits"]
    score -= weights["motif"] * metrics["motif_hits"]

    return score

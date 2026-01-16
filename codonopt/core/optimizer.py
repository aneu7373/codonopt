
from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

from codonopt.core.orf import validate_orf
from codonopt.metrics.basic import gc_content, shannon_entropy
from codonopt.metrics.motifs import count_hits
from codonopt.metrics.rnafold import rnafold_mfe, rnafold_mfe_5p
from codonopt.core.scoring import score_candidate
from codonopt.core.generation import generate_candidate


def _get_weight(weights: Mapping[str, float], key: str, default: float = 0.0) -> float:
    return float(weights.get(key, default)) if weights is not None else default


def _violates_gc_bounds(gc: float, gc_min: float | None, gc_max: float | None) -> bool:
    if gc_min is not None and gc < gc_min:
        return True
    if gc_max is not None and gc > gc_max:
        return True
    return False


def _mean(values: Sequence[float]) -> float:
    return sum(values) / max(1, len(values))


def optimize(
    dna: str,
    aa_seq: Iterable[str],
    codon_table: Dict[str, Sequence[Tuple[str, float]]],
    forbidden_codons: set[str],
    restriction_sites: List[str],
    motifs: List[str],
    n_out: int,
    oversample: int,
    weights: Mapping[str, float],
    *,
    max_trials: int | None = None,
    rna_5p_window: int = 60,
) -> List[Tuple[str, Dict[str, Any]]]:
    """
    Generate many synonymous candidates, score them, and return the top N.

    - Hard pre-filtering on GC bounds to avoid unnecessary RNAfold calls.
    - Forbidden codons are best-effort avoided; unavoidable hits are penalized via weights["forbidden"].
    - RNAfold calls are skipped entirely when both 'rna_5p' and 'rna_full' weights are 0.0.
    """
    candidates: List[Tuple[str, Dict[str, Any]]] = []

    target = max(1, n_out * max(1, oversample))
    if max_trials is None:
        max_trials = target * 50  # generous cap to prevent infinite loops

    # Weights/knobs
    gc_min = weights.get("gc_min")
    gc_max = weights.get("gc_max")
    use_rnafold = not (
        _get_weight(weights, "rna_5p", 0.0) == 0.0
        and _get_weight(weights, "rna_full", 0.0) == 0.0
    )

    # Normalize AA input
    aa_str = "".join(aa_seq) if not isinstance(aa_seq, str) else aa_seq

    trials = 0
    while len(candidates) < target and trials < max_trials:
        trials += 1

        # 1) Generate candidate (frequency-weighted; avoid forbidden where possible)
        seq, forbidden_hits, codon_freqs = generate_candidate(
            aa_seq=aa_str, codon_table=codon_table, forbidden_codons=forbidden_codons
        )

        # 2) Ensure synonymous integrity
        if not validate_orf(seq, aa_str):
            continue

        # 3) Cheap metrics & early GC pruning
        gc = gc_content(seq)
        if (gc_min is not None or gc_max is not None) and _violates_gc_bounds(gc, gc_min, gc_max):
            continue

        mean_codon_freq = _mean(codon_freqs)
        entropy = shannon_entropy(codon_freqs)

        # 4) Medium-cost metrics
        restriction_hits = count_hits(seq, restriction_sites)
        motif_hits = count_hits(seq, motifs)

        # 5) Expensive metrics (RNAfold) — only if enabled
        if use_rnafold:
            rna_5p = rnafold_mfe_5p(seq, window=rna_5p_window, disabled=False)
            rna_full = rnafold_mfe(seq, disabled=False)
        else:
            rna_5p = 0.0
            rna_full = 0.0

        # 6) Score and keep
        metrics: Dict[str, Any] = {
            "gc": gc,
            "mean_codon_freq": mean_codon_freq,
            "entropy": entropy,
            "restriction_hits": restriction_hits,
            "motif_hits": motif_hits,
            "rna_5p": rna_5p,
            "rna_full": rna_full,
            "forbidden_hits": forbidden_hits,
        }
        metrics["score"] = score_candidate(metrics, weights)
        candidates.append((seq, metrics))

    # Sort by descending score and return top N
    candidates.sort(key=lambda x: x[1]["score"], reverse=True)
    return candidates[:n_out]

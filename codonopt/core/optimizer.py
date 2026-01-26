# codonopt/core/optimizer.py

import random
from typing import List, Tuple, Dict

from codonopt.core.generation import filter_codons_for_constraints
from codonopt.exceptions import ConstraintError
from codonopt.utils import calculate_gc


def _weighted_choice(codons: List[str], weights: List[float]) -> str:
    if not codons:
        raise ConstraintError("No candidate codons available after filtering.")
    if len(codons) != len(weights):
        raise ConstraintError("Internal error: codons/weights length mismatch.")
    # random.choices allows non-normalized weights
    return random.choices(codons, weights=weights, k=1)[0]


def _weighted_permutation(codons: List[str], weights: List[float]) -> List[str]:
    """
    Return a weighted random ordering (without replacement), biased to try higher-weight codons earlier.

    Uses Efraimidis-Spirakis style keys: key = U^(1/w)
    Higher weights tend to yield larger keys, thus appear earlier when sorting descending.
    """
    if not codons:
        return []
    if len(codons) != len(weights):
        raise ConstraintError("Internal error: codons/weights length mismatch.")

    items = []
    for c, w in zip(codons, weights):
        w = float(w)
        if w <= 0:
            # extremely rare/invalid, shove to end deterministically
            key = -1.0
        else:
            u = random.random()
            # u ** (1/w)
            key = u ** (1.0 / w)
        items.append((key, c))

    items.sort(reverse=True, key=lambda x: x[0])
    return [c for _, c in items]


def optimize_sequence(
    dna=None,
    protein=None,
    codon_table=None,
    avoid_codons=None,
    avoid_motifs=None,
    max_homopolymer=5,
    gc_min=None,
    gc_max=None,
    logger=None,
    optimization_mode="kleinbub",
    max_attempts=1000,
    backtrack_window=10,
    min_codon_fraction=0.05,
):
    """
    optimization_mode:
      - 'strict': greedy, retries per position, fast
      - 'kleinbub': bounded backtracking search, higher yield, slower

    codon_table must be:
      dict: AA -> list of (codon, fraction)

    min_codon_fraction:
      excludes codons for an AA if fraction < cutoff.
      default 0.05.
    """
    avoid_codons = avoid_codons or []
    avoid_motifs = avoid_motifs or []

    if protein is None:
        raise ConstraintError("optimize_sequence requires protein input (DNA is translated upstream).")
    if codon_table is None:
        raise ConstraintError("codon_table is required.")

    mode = (optimization_mode or "kleinbub").strip().lower()
    if mode not in ("strict", "kleinbub"):
        raise ConstraintError("optimization_mode must be 'strict' or 'kleinbub'")

    try:
        min_cf = float(min_codon_fraction)
    except Exception:
        raise ConstraintError(f"min_codon_fraction must be a float; got {min_codon_fraction!r}")
    if min_cf < 0 or min_cf > 1:
        raise ConstraintError(f"min_codon_fraction must be between 0 and 1; got {min_cf}")

    # Build AA -> (codons, weights) after rarity cutoff
    aa_codons: Dict[str, Tuple[List[str], List[float]]] = {}
    for aa in set(protein):
        rows = codon_table.get(aa, [])
        if not rows:
            raise ConstraintError(f"No codons available in codon table for amino acid '{aa}'")
        filtered = [(c, float(fr)) for c, fr in rows if float(fr) >= min_cf]
        if not filtered:
            raise ConstraintError(
                f"All codons for amino acid '{aa}' are below min_codon_fraction={min_cf}. "
                f"Lower the cutoff or change the codon table."
            )
        codons = [c for c, _ in filtered]
        weights = [fr for _, fr in filtered]
        aa_codons[aa] = (codons, weights)

    # ---------------------------
    # STRICT MODE
    # ---------------------------
    if mode == "strict":
        seq = ""
        for aa in protein:
            codons, weights = aa_codons[aa]
            attempts = 0
            while attempts < max_attempts:
                # apply motif/homopolymer/codon constraints
                valid_codons = filter_codons_for_constraints(
                    codons, seq, avoid_codons, avoid_motifs, max_homopolymer
                )
                # map weights for valid codons
                wmap = {c: w for c, w in zip(codons, weights)}
                valid_weights = [wmap[c] for c in valid_codons]

                codon = _weighted_choice(valid_codons, valid_weights)
                candidate = seq + codon

                gc = calculate_gc(candidate)
                if ((gc_min is None or gc >= gc_min) and (gc_max is None or gc <= gc_max)):
                    seq = candidate
                    if logger:
                        logger.debug(f"{aa} -> {codon} (strict)")
                    break

                attempts += 1

            else:
                raise ConstraintError(
                    f"Strict optimization failed at amino acid '{aa}'. "
                    f"Try optimization_mode=kleinbub or relax constraints."
                )

        return seq

    # ---------------------------
    # KLEINBUB MODE (bounded backtracking)
    # ---------------------------
    if backtrack_window < 1:
        backtrack_window = 1

    chosen: List[Tuple[str, str]] = []  # (aa, codon)
    seq = ""
    i = 0
    total_steps = 0

    # Treat max_attempts as total search budget
    search_limit = int(max_attempts)

    candidate_lists = [None] * len(protein)
    candidate_index = [0] * len(protein)

    while i < len(protein):
        total_steps += 1
        if total_steps > search_limit:
            raise ConstraintError(
                f"Kleinbub optimization exceeded search limit ({search_limit}). "
                "Increase kleinbub_search_limit or relax constraints."
            )

        aa = protein[i]
        codons, weights = aa_codons[aa]

        # Build candidate list for this position given current prefix
        if candidate_lists[i] is None:
            valid_codons = filter_codons_for_constraints(
                codons, seq, avoid_codons, avoid_motifs, max_homopolymer
            )
            if not valid_codons:
                # force backtrack logic below
                candidate_lists[i] = []
                candidate_index[i] = 0
            else:
                wmap = {c: w for c, w in zip(codons, weights)}
                valid_weights = [wmap[c] for c in valid_codons]
                # weighted order: try higher-weight codons earlier (stochastic)
                ordered = _weighted_permutation(valid_codons, valid_weights)
                candidate_lists[i] = ordered
                candidate_index[i] = 0

        # If exhausted candidates: backtrack
        if candidate_index[i] >= len(candidate_lists[i]):
            candidate_lists[i] = None
            candidate_index[i] = 0

            if i == 0:
                raise ConstraintError(
                    "Kleinbub optimization could not satisfy constraints from the start. "
                    "Constraints likely incompatible."
                )

            back_to = max(0, i - backtrack_window)

            # rewind chosen list
            while len(chosen) > back_to:
                chosen.pop()

            seq = "".join(c for (_, c) in chosen)

            # reset positions from back_to onward
            for j in range(back_to, len(protein)):
                candidate_lists[j] = None
                candidate_index[j] = 0

            i = back_to
            continue

        # Try next candidate codon
        codon = candidate_lists[i][candidate_index[i]]
        candidate_index[i] += 1

        candidate_seq = seq + codon
        gc = calculate_gc(candidate_seq)
        if (gc_min is not None and gc < gc_min) or (gc_max is not None and gc > gc_max):
            continue

        # Accept and advance
        chosen.append((aa, codon))
        seq = candidate_seq
        if logger:
            logger.debug(f"{aa} -> {codon} (kleinbub)")
        i += 1

    return seq

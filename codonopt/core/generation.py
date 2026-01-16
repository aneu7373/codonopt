
from __future__ import annotations

import random
from typing import Dict, Iterable, List, Sequence, Set, Tuple

CodonChoice = Tuple[str, float]               # (codon, frequency)
CodonTable = Dict[str, Sequence[CodonChoice]] # AA -> list of choices


def generate_candidate(
    aa_seq: Iterable[str],
    codon_table: CodonTable,
    forbidden_codons: Set[str] | None = None,
    *,
    max_avoid_attempts: int = 20,
) -> tuple[str, int, List[float]]:
    """
    Generate one DNA candidate by frequency-weighted sampling of synonymous codons.
    Best-effort avoidance of forbidden codons; unavoidable uses are counted.

    Returns
    -------
    dna : str
    forbidden_hits : int
        Positions where a forbidden codon still ended up being chosen.
    codon_freqs : list[float]
        Usage frequency per chosen codon (same length as AA sequence).
    """
    forbidden_codons = {c.upper() for c in (forbidden_codons or set())}

    dna_codons: List[str] = []
    codon_freqs: List[float] = []
    forbidden_hits = 0

    for aa in aa_seq:
        choices = codon_table[aa]
        codons, freqs = zip(*choices)

        selected = None
        attempts = 0
        while attempts < max_avoid_attempts:
            candidate = random.choices(codons, weights=freqs, k=1)[0]
            if candidate.upper() not in forbidden_codons:
                selected = candidate
                break
            attempts += 1

        if selected is None:
            # Could not avoid — accept weighted pick and record a hit
            selected = random.choices(codons, weights=freqs, k=1)[0]
            forbidden_hits += 1

        dna_codons.append(selected)
        freq_map = dict(choices)
        codon_freqs.append(float(freq_map[selected]))

    return "".join(dna_codons), forbidden_hits, codon_freqs


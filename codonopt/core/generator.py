import random

def generate_candidate(aa_seq, codon_table, forbidden_codons):
    dna = []
    forbidden_hits = 0
    codon_freqs = []

    for aa in aa_seq:
        choices = codon_table[aa]
        codons, freqs = zip(*choices)

        codon = random.choices(codons, weights=freqs, k=1)[0]
        if codon in forbidden_codons:
            forbidden_hits += 1

        dna.append(codon)
        codon_freqs.append(dict(choices)[codon])

    return "".join(dna), forbidden_hits, codon_freqs

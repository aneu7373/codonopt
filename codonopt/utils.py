def is_protein_sequence(seq):
    dna_chars = set("ATGCN")
    return any(c not in dna_chars for c in seq.upper())

def calculate_gc(seq):
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return gc_count / len(seq) if len(seq) > 0 else 0

def max_homopolymer_length(seq):
    import re
    runs = re.findall(r"(A+|T+|G+|C+)", seq)
    return max([len(r) for r in runs], default=0)

def is_protein_sequence(seq):
    # simple check: no more than 50% ATGC means likely protein
    seq = seq.upper()
    nt_count = sum(seq.count(n) for n in "ATGC")
    return nt_count / len(seq) < 0.5
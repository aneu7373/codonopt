import re
from codonopt.exceptions import ConstraintError

def violates_homopolymer(seq, max_run):
    runs = re.findall(r"(A+|T+|G+|C+)", seq)
    return any(len(r) > max_run for r in runs)

def filter_codons_for_constraints(codons, current_seq, avoid_codons, avoid_motifs, max_homopolymer):
    valid = []
    for codon in codons:
        test_seq = current_seq[-(max_homopolymer-1):] + codon
        if violates_homopolymer(test_seq, max_homopolymer):
            continue
        if codon in avoid_codons:
            continue
        if any(m in current_seq + codon for m in avoid_motifs):
            continue
        valid.append(codon)
    if not valid:
        raise ConstraintError("No valid codons remaining due to constraints")
    return valid

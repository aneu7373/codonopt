from Bio.Seq import Seq

def translate_dna(dna):
    return str(Seq(dna).translate(to_stop=False))

def validate_orf(candidate_dna, reference_aa):
    aa = translate_dna(candidate_dna)
    return aa == reference_aa

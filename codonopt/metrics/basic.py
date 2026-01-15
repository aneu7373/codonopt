import math

def gc_content(seq):
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)

def shannon_entropy(values):
    entropy = 0.0
    for v in values:
        if v > 0:
            entropy -= v * math.log2(v)
    return entropy

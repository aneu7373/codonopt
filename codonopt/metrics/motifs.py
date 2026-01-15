def count_hits(seq, motifs):
    hits = 0
    for m in motifs:
        hits += seq.count(m)
    return hits

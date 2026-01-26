import random

def choose_random_codon(candidates, seed=None):
    if seed is not None:
        random.seed(seed)
    return random.choice(candidates)


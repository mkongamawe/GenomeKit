import random

from genomekit import GenomeKit


def happy_sequence():
    random.seed(26)
    dna = ''.join(random.choices(["a", "c", "t", "g"], k = 50))
    dna = GenomeKit(dna)
    assert not dna.find_primers().forward.verdict
    assert not dna.find_primers().reverse.verdict

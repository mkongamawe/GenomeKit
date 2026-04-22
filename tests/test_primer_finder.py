import pytest

from genomekit.modules.primer_finder import PrimerFinder


def test_primer_finder_happy_path():
    """A normal 50-bp sequence with moderate GC."""
    finder = PrimerFinder("ccatcttcttcatagattttattactgcgtacggacggattcacggggat")
    result = finder.find()

    assert len(result) == 1
    assert result[0]["forward_primer"] == "CCATCTTCTTCATAGATTTT"
    assert result[0]["forward_verdict"] is False
    assert result[0]["reverse_primer"] == "ACGGACGGATTCACGGGGAT"
    assert result[0]["reverse_verdict"] is True


def test_primer_finder_extreme_gc():
    """Edge case: 100% GC at both ends — should fail the GC window."""
    finder = PrimerFinder("G" * 20 + "A" * 10 + "C" * 20)
    result = finder.find()

    assert result[0]["forward_verdict"] is False
    assert result[0]["reverse_verdict"] is False


def test_primer_finder_short_sequence():
    """Error case: sequence shorter than 20 bp cannot yield a primer."""
    with pytest.raises(ValueError, match="at least 20 bases"):
        finder = PrimerFinder("ATCG")
        finder.find()

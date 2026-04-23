import pytest

from genomekit.modules.primer_finder import Primer


def test_primer_finder_happy_path():
    """A normal 50-bp sequence with moderate GC."""
    finder = Primer("ccatcttcttcatagattttattactgcgtacggacggattcacggggat")
    result = finder.find_primer()

    assert isinstance(result, Primer.PrimerFinderResults)
    assert result.forward == "CCATCTTCTTCATAGATTTT"
    assert result.for_verdict is False
    assert result.reverse == "ACGGACGGATTCACGGGGAT"
    assert result.rev_verdict is True


def test_multiple_primer():
    """Three normal DNA sequences"""
    finder = Primer.from_multiple(
        [
            "tcttgaacattgacaattactaatacctcgtataccataaaggtgtcacc",
            "tttagatagagccctgtcggcggtgagtctcccaccattccggctacgcc",
            "cggtagaaattcaatcatgggaaagggcgctaccgtttcccagagcttca",
        ]
    )
    result = [s.find_primer() for s in finder]

    assert isinstance(result[0], Primer.PrimerFinderResults)
    assert result[0].for_verdict is False
    assert result[1].for_verdict is True
    assert result[2].rev_verdict is True


def test_primer_finder_extreme_gc():
    """Edge case: 100% GC at both ends — should fail the GC window."""
    finder = Primer("G" * 20 + "A" * 10 + "C" * 20)
    result = finder.find_primer()

    assert result.for_verdict is False
    assert result.rev_verdict is False


def test_primer_finder_short_sequence():
    """Error case: sequence shorter than 20 bp cannot yield a primer."""
    with pytest.raises(ValueError, match="at least 20 bases"):
        finder = Primer("ATCG")
        finder.find_primer()

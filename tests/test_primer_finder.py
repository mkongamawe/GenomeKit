import random
import time

import pytest

from genomekit import GenomeKit
from genomekit.modules.primer_finder import BatchSummary, PrimerAnalysis, SequenceResult


def test_happy_sequence():
    random.seed(26)
    dna = "".join(random.choices(["a", "c", "t", "g"], k=50))
    dna = GenomeKit(dna)
    assert not dna.find_primers().forward.verdict
    assert not dna.find_primers().reverse.verdict


def test_happy_sequence_arguments():
    random.seed(26)
    dna = "".join(random.choices(["a", "c", "t", "g"], k=100))
    dna = GenomeKit(dna)
    assert not dna.find_primers(
        primer_length=25, gc_range=(50, 70), tm_range=(50, 55)
    ).forward.verdict
    assert not dna.find_primers(
        primer_length=25, gc_range=(50, 70), tm_range=(50, 55)
    ).reverse.verdict


# -------------------------------------------------------------------------------------------
# Helper lines
# -------------------------------------------------------------------------------------------
def make_sequences(length: int = 50) -> str:
    return "".join(random.choices(["a", "c", "t", "g"], k=length))


sequences = [make_sequences() for _ in range(10)]
ids = [f"sample_{i:03d}" for i in range(10)]

# A high GC sequence
# A sequence guaranteed to have high GC (primers will likely pass gc check)
HIGH_GC = "GCGCGCGCGCGCGCGCGCGC" + "AAAAAAAAAA" + "GCGCGCGCGCGCGCGCGCGC"

# A sequence guaranteed to have low GC (primers will likely fail gc check)
LOW_GC = "ATATATATATATATATATATAT" + "CCCCCCCCCC" + "ATATATATATATATATATATAT"

# -------------------------------------------------------------------------------------------
# construction tests
# -------------------------------------------------------------------------------------------


def test_batch_returns_genomekit_instace():
    kit = GenomeKit.batch(sequences)
    assert isinstance(kit, GenomeKit)


def test_batch_stores_sequences():
    kit = GenomeKit.batch(sequences)
    assert kit._sequences == [s.upper() for s in sequences]


def test_batch_stores_ids_when_provided():
    kit = GenomeKit.batch(sequences, ids)
    assert kit._ids == ids


def test_batch_ids_default_to_none():
    kit = GenomeKit.batch(sequences)
    assert kit._ids is None


def test_batch_mode_is_set():
    kit = GenomeKit.batch(sequences)
    assert kit._mode == "batch"


def test_batch_rejects_non_list():
    with pytest.raises(TypeError):
        GenomeKit.batch(tuple(sequences))


def test_batch_rejects_empty_list():
    with pytest.raises(ValueError):
        GenomeKit.batch([])


def test_batch_rejects_mismatched_ids():
    with pytest.raises(ValueError):
        GenomeKit.batch(sequences, ids=["only_one_id"])


def test_batch_rejects_invalid_bases():
    bad = sequences[:] + ["ATCGHXTZ"]
    with pytest.raises(ValueError, match="invalid characters"):
        GenomeKit.batch(bad)


def test_batch_rejects_empty_sequences_in_list():
    # [ ] Ensure the error message for empty sequences is updated
    with pytest.raises(ValueError):
        GenomeKit.batch(sequences + [""])


def test_batch_uppercase_sequences():
    lower = ["atcgatcgatcgatcgatcgatcgatcg"] * 3
    kit = GenomeKit.batch(lower)
    assert all(s == s.upper() for s in kit._sequences)


# ---------------------------------------------------------------------------
# Mode guard tests
# ---------------------------------------------------------------------------


def test_analyse_raises_in_single_mode():
    kit = GenomeKit(sequences[0])
    with pytest.raises(RuntimeError, match="analyse()"):
        kit.analyse()


def test_find_primers_raises_in_batch_mode():
    kit = GenomeKit.batch(sequences)
    with pytest.raises(RuntimeError, match="find_primers()"):
        kit.find_primers()


# ----------------------------------------------------------------------------
# Return type tests
# ----------------------------------------------------------------------------


def test_analyse_returns_batch_summary():
    kit = GenomeKit.batch(sequences)
    summary = kit.analyse()
    assert isinstance(summary, BatchSummary)


def test_batch_summary_total_matches_input():
    kit = GenomeKit.batch(sequences)
    summary = kit.analyse()
    assert summary.total == len(sequences)


def test_batch_summary_results_are_sequence_results():
    kit = GenomeKit.batch(sequences)
    summary = kit.analyse()
    assert all(isinstance(r, SequenceResult) for r in summary.results)


def test_each_result_has_primer_analysis():
    kit = GenomeKit.batch(sequences)
    summary = kit.analyse()
    for r in summary.results:
        assert isinstance(r.forward, PrimerAnalysis)
        assert isinstance(r.reverse, PrimerAnalysis)


# ----------------------------------------------------------------------------
# ID traceability tests
# ----------------------------------------------------------------------------


def test_results_carry_provided_ids():
    kit = GenomeKit.batch(sequences, ids)
    summary = kit.analyse()
    result_ids = [r.source_id for r in summary.results]
    assert result_ids == ids


def test_results_carry_auto_ids_when_none_provided():
    kit = GenomeKit.batch(sequences)
    summary = kit.analyse()
    result_ids = [r.source_id for r in summary.results]
    expected = [f"seq_{i:04d}" for i in range(len(sequences))]
    assert result_ids == expected


def test_results_carry_source_sequence():
    kit = GenomeKit.batch(sequences, ids=ids)
    summary = kit.analyse()
    for r, original in zip(summary.results, sequences, strict=True):
        assert r.source_sequence == original.upper()


# -----------------------------------------------------------------------
# Filter tests
# -----------------------------------------------------------------------


def test_full_pass_is_subset_of_total():
    kit = GenomeKit.batch(sequences)
    summary = kit.analyse()
    assert len(summary.full_pass) <= summary.total


def test_forward_only_excludes_full_pass():
    kit = GenomeKit.batch(sequences)
    summary = kit.analyse()
    for r in summary.forward_only:
        assert r.forward.verdict is True
        assert r.reverse.verdict is False


def test_reverse_only_excludes_full_pass():
    kit = GenomeKit.batch(sequences)
    summary = kit.analyse()
    for r in summary.reverse_only:
        assert r.reverse.verdict is True
        assert r.forward.verdict is False


def test_no_pass_has_no_verdicts():
    kit = GenomeKit.batch(sequences)
    summary = kit.analyse()
    for r in summary.no_pass:
        assert r.forward.verdict is False
        assert r.reverse.verdict is False


def test_all_categories_are_mutually_exclusive_and_exhaustive():
    kit = GenomeKit.batch(sequences)
    summary = kit.analyse()
    total_accounted = (
        len(summary.full_pass)
        + len(summary.forward_only)
        + len(summary.reverse_only)
        + len(summary.no_pass)
    )
    assert total_accounted == summary.total


# -----------------------------------------------------------------------
# get_sequences and get_ids helpers
# -----------------------------------------------------------------------


def test_get_ids_returns_ids_of_filtered_results():
    kit = GenomeKit.batch(sequences, ids=ids)
    summary = kit.analyse()
    valid = summary.full_pass
    ids_out = summary.get_ids(valid)
    assert all(i in ids for i in ids_out)


def test_get_sequences_returns_sequences_of_filtered_results():
    kit = GenomeKit.batch(sequences, ids=ids)
    summary = kit.analyse()
    valid = summary.full_pass
    seqs_out = summary.get_sequences(valid)
    assert all(s in [seq.upper() for seq in sequences] for s in seqs_out)


def test_get_sequences_and_get_ids_same_length():
    kit = GenomeKit.batch(sequences, ids=ids)
    summary = kit.analyse()
    valid = summary.full_pass
    assert len(summary.get_sequences(valid)) == len(summary.get_ids(valid))


# -----------------------------------------------------------------------
# Performance test
# -----------------------------------------------------------------------


def test_batch_performance_100k_sequences():
    sequences = [make_sequences(100) for _ in range(100_000)]

    start = time.perf_counter()
    summary = GenomeKit.batch(sequences).analyse()
    _ = summary.total  # force materialisation
    elapsed = time.perf_counter() - start

    assert summary.total == 100_000
    assert elapsed < 12.0, f"Batch took too long: {elapsed:.2f}s"

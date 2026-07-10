import random
import time

from genomekit import GenomeKit
from genomekit.modules.gc_calculator import GCBatchSummary

# --------------------------------------------------------
# Simple single sequence tests
# --------------------------------------------------------


def test_happy_sequence():
    random.seed(26)
    dna = "".join(random.choices(["a", "c", "t", "g"], k=50))
    dna = GenomeKit(dna)
    assert not dna.gc_content().gc_verdict


def test_only_gc_sequence():
    random.seed(26)
    dna = "".join(random.choices(["c", "g"], k=50))
    dna = GenomeKit(dna)
    dna.gc_content(50)
    assert dna.gc_content().gc_content == 100
    assert dna.gc_content().gc_ratio == float("inf")


def test_only_ac_sequence():
    random.seed(26)
    dna = "".join(random.choices(["a", "t"], k=50))
    dna = GenomeKit(dna)
    assert not dna.gc_content().gc_verdict
    assert dna.gc_content().gc_content == 0


# ----------------------------------------------------------
# Test helper code
# ----------------------------------------------------------
def make_sequences(length: int = 50) -> str:
    return "".join(random.choices(["a", "c", "t", "g"], k=length))


sequences = [make_sequences() for _ in range(10)]
ids = [f"sample_{i:03d}" for i in range(10)]

# ----------------------------------------------------------
# batch GC analysis
# ----------------------------------------------------------


def test_batch_return_genomekit_instance():
    kit = GenomeKit.batch(sequences)
    summary = kit.gc_content_batch()
    assert isinstance(summary, GCBatchSummary)


def test_batch_summary_total_matches_input():
    kit = GenomeKit.batch(sequences)
    summary = kit.gc_content_batch()
    assert summary.total == len(sequences)


# -----------------------------------------------------------------------
# Performance test
# -----------------------------------------------------------------------


def test_batch_performance_100k_sequences():
    sequences = [make_sequences(100) for _ in range(100_000)]

    start = time.perf_counter()
    summary = GenomeKit.batch(sequences).gc_content_batch()
    _ = summary.total  # force materialisation
    elapsed = time.perf_counter() - start

    assert summary.total == 100_000
    assert elapsed < 12.0, f"Batch took too long: {elapsed:.2f}s"

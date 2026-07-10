# This script has been written by AI, inspiration has been taken from primer_finder
# The goal of this class is to calculate the GC content of the sequence
# But any function can have extended ability i.e GC ratio, region of the
# sequence with high GC and more
# All these have to be done for both a single sequence and a batch of sequences.

from __future__ import annotations

import inspect
import warnings
from collections.abc import Iterator
from dataclasses import dataclass


# First we create a dataclass. A dataclass' work is to store data only
@dataclass(frozen=True)  # Since it is immutable, we add the frozen argument
class GCSequence:
    """
    A value object that represents a sequence instance for gc calculation to occur
    It is Immutable
    """

    id: str
    sequence: str
    source_id: str

    # A dataclass does not need to be instantialised. that occues automatically

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return self.sequence


# Next, we need to store how the results will be stored.
@dataclass  # This is mutable, hence we do not have a frozen argument
class GCResult:
    """
    This contains all the results contained in the GC calculator analysis
    """

    source_id: str
    sequence: str
    gc_content: float
    gc_verdict: bool
    at_content: float
    gc_ratio: float
    # TODO: Implement adding highest GC range - Region in sequence with high GC
    # highest_gc_range: str
    total_length: float


# Now that we have our input class sorted and our output class sorted. We can build the class
# that will analyse the sequence. This is a stateless analyser


class GCCalculator:
    # This is not a dataclass, therefore it needs to be instantiated
    def __init__(self, high_gc_threshold: int = 60):
        self.high_gc_threshold = high_gc_threshold

    # The main function that analyses the GC content of a sequence
    def analyse(self, sequence: str, source_id: str = "unknown") -> GCResult:
        """Analyse a sequence and return the results"""
        sequence = sequence.upper()

        # TODO: Sort out how the warning message appears in terminal.
        if 0 < self.high_gc_threshold > 100:
            warnings.warn(
                f"Warning: {self.high_gc_threshold} is outside the expected"
                f"range [0-100], defaulting to 60",
                UserWarning,
                stacklevel=2,
            )
            self.high_gc_threshold = (
                inspect.signature(self.analyse).parameters["high_gc_content"].default
            )

        seq_len = len(sequence)
        if seq_len == 0:
            raise ValueError("Sequence cannot be empty")

        gc = sequence.count("G") + sequence.count("C")
        at = sequence.count("A") + sequence.count("T")

        gc_content = (gc / len(sequence)) * 100
        at_content = (at / len(sequence)) * 100
        gc_ratio = gc / at if at > 0 else float("inf")

        if gc_content >= self.high_gc_threshold:
            gc_verdict = True
        else:
            gc_verdict = False

        return GCResult(
            source_id=source_id,
            sequence=sequence,
            gc_content=round(gc_content, 2),
            gc_verdict=gc_verdict,
            at_content=round(at_content, 2),
            gc_ratio=round(gc_ratio, 2),
            total_length=seq_len,
        )

    # Now a function that analyses GC for a batch of sequences
    def analyse_batch(
        self, sequences: list[str], ids: list[str] | None = None
    ) -> Iterator[GCResult]:
        """
        Lazily analyses a batch of sequence.
        IDs default to seq_0001, seq_0002... if not provided
        """
        if ids is None:
            ids = [f"seq_{i:04d}" for i in range(len(sequences))]

        if len(ids) != len(sequences):
            raise ValueError("Ids and sequences must be of same length")

        for seq_id, seq in zip(ids, sequences, strict=True):
            yield self.analyse(seq, source_id=seq_id)


# My goal for GenomeKit is that it stops being memory intensive
# So that whole batch of sequences is not loaded to memory, we only return a summary
# Therefore we create a summary class


class GCBatchSummary:
    def __init__(self, results: Iterator[GCResult]):
        self._stream = results
        self._cache: list[GCResult] | None = None

    def _materialise(self) -> list[GCResult]:
        """Consume the stream into cache"""
        if self._cache is None:
            self._cache = list(self._stream)
        return self._cache

    @property
    def total(self) -> int:
        return len(self._materialise())

    @property
    def high_gc(self) -> list[GCResult]:
        return [r for r in self._materialise() if r.gc_verdict]

    @property
    def results(self) -> list[GCResult]:
        return self._materialise()

    def get_sequences(self, results: list[GCResult]) -> list[str]:
        return [r.sequence for r in results]

    def get_ids(self, results: list[GCResult]) -> list[str]:
        return [r.source_id for r in results]

    # Now we set up the __repr_ function for how a batch summary should display
    def __repr__(self):
        return (
            f"\n{'Batch Summary ':=^40}\n"
            f"Total processed   : {self.total}\n"
            f"High GC sequences : {len(self.high_gc)}\n"
            f"{'=' * 40}"
        )

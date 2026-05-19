from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from enum import Enum


class MeltTempMethod(Enum):
    WALLACE = "wallace"  # 4(G+C) + 2(A+T) — fast, short oligos
    NEAREST_NEIGHBOUR = "nn"  # More accurate, longer sequences


@dataclass(frozen=True)
class PrimerSequence:
    """
    A value object representing a single primer candidate.
    Immutable — a primer sequence is a fact, not something that changes.
    """

    id: str  # e.g. "seq_001_forward"
    sequence: str
    source_id: str  # ID of the sequence it was derived from
    orientation: str  # "forward" or "reverse"

    def __post_init__(self):
        invalid = set(self.sequence.upper()) - {"A", "T", "C", "G"}
        if invalid:
            raise ValueError(f"Invalid bases in primer: {sorted(invalid)}")

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return self.sequence


@dataclass
class PrimerAnalysis:
    """
    Holds the full analysis for a single primer candidate.
    Keeps the result coupled to its source — no orphaned results.
    """

    primer: PrimerSequence
    gc_content: float
    melting_temp: float
    has_hairpin: bool
    tm_method: MeltTempMethod

    @property
    def gc_pass(self) -> bool:
        return 40.0 <= self.gc_content <= 60.0

    @property
    def tm_pass(self) -> bool:
        return 50.0 <= self.melting_temp <= 65.0

    @property
    def verdict(self) -> bool:
        return self.gc_pass and self.tm_pass and not self.has_hairpin

    def __repr__(self):
        return (
            f"PrimerAnalysis(id={self.primer.id!r}, "
            f"orientation={self.primer.orientation!r}, "
            f"gc={self.gc_content:.1f}%, "
            f"tm={self.melting_temp:.1f}°C, "
            f"hairpin={self.has_hairpin}, "
            f"verdict={self.verdict})"
        )


@dataclass
class SequenceResult:
    """
    Pairs the forward and reverse analysis for one input sequence.
    This is what you hand off to downstream tools.
    """

    source_id: str
    source_sequence: str
    forward: PrimerAnalysis
    reverse: PrimerAnalysis

    @property
    def full_pass(self) -> bool:
        return self.forward.verdict and self.reverse.verdict

    @property
    def partial_pass(self) -> bool:
        return self.forward.verdict or self.reverse.verdict


class PrimerAnalyser:
    """
    Stateless analyser — takes sequences in, produces results out.
    No state is held between calls.
    """

    def __init__(
        self,
        primer_length: int = 20,
        gc_range: tuple[float, float] = (40.0, 60.0),
        tm_range: tuple[float, float] = (50.0, 65.0),
        tm_method: MeltTempMethod = MeltTempMethod.WALLACE,
    ):
        self.primer_length = primer_length
        self.gc_range = gc_range
        self.tm_range = tm_range
        self.tm_method = tm_method

    def analyse(self, sequence: str, source_id: str = "unknown") -> SequenceResult:
        """Analyse a single sequence and return a fully-labelled result."""
        sequence = sequence.upper()

        if len(sequence) < self.primer_length * 2:
            raise ValueError(
                f"Sequence '{source_id}' is too short for primer length {self.primer_length}. "
                f"Minimum length is {self.primer_length * 2}bp."
            )

        fwd_seq = sequence[: self.primer_length]
        rev_seq = sequence[-self.primer_length :]

        forward = self._analyse_candidate(fwd_seq, source_id, "forward")
        reverse = self._analyse_candidate(rev_seq, source_id, "reverse")

        return SequenceResult(
            source_id=source_id,
            source_sequence=sequence,
            forward=forward,
            reverse=reverse,
        )

    def analyse_batch(
        self,
        sequences: list[str],
        ids: list[str] | None = None,
    ) -> Iterator[SequenceResult]:
        """
        Lazily analyse a batch of sequences.
        IDs default to seq_0001, seq_0002 ... if not provided.
        """
        if ids is None:
            ids = [f"seq_{i:04d}" for i in range(len(sequences))]

        if len(ids) != len(sequences):
            raise ValueError("ids and sequences must be the same length")

        for seq_id, seq in zip(ids, sequences, strict=True):
            yield self.analyse(seq, source_id=seq_id)

    def _analyse_candidate(self, seq: str, source_id: str, orientation: str) -> PrimerAnalysis:
        primer = PrimerSequence(
            id=f"{source_id}_{orientation}",
            sequence=seq,
            source_id=source_id,
            orientation=orientation,
        )
        return PrimerAnalysis(
            primer=primer,
            gc_content=self._gc_content(seq),
            melting_temp=self._melt_temp(seq),
            has_hairpin=self._find_hairpin(seq),
            tm_method=self.tm_method,
        )

    # --- Internal calculations (stateless, could even be module-level) ---

    def _gc_content(self, seq: str) -> float:
        if not seq:
            return 0.0
        return ((seq.count("G") + seq.count("C")) / len(seq)) * 100

    def _melt_temp(self, seq: str) -> float:
        if not seq:
            return 0.0
        if self.tm_method == MeltTempMethod.WALLACE:
            return 4 * (seq.count("G") + seq.count("C")) + 2 * (seq.count("A") + seq.count("T"))
        raise NotImplementedError("Nearest-neighbour Tm not yet implemented")

    def _find_hairpin(self, seq: str) -> bool:
        comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
        for i in range(len(seq)):
            for j in range(i + 6, len(seq) - 2):
                stem_a = seq[i : i + 3]
                stem_b = seq[j : j + 3]
                target = "".join(comp.get(b, "N") for b in stem_a)[::-1]
                if target == stem_b:
                    return True
        return False


class BatchSummary:
    """
    Consumes a stream of SequenceResults and builds a report.
    Retains results for downstream filtering — nothing is thrown away.
    """

    def __init__(self, results: Iterator[SequenceResult]):
        self._results: list[SequenceResult] = list(results)

    @property
    def total(self) -> int:
        return len(self._results)

    @property
    def full_pass(self) -> list[SequenceResult]:
        return [r for r in self._results if r.full_pass]

    @property
    def forward_only(self) -> list[SequenceResult]:
        return [r for r in self._results if r.forward.verdict and not r.reverse.verdict]

    @property
    def reverse_only(self) -> list[SequenceResult]:
        return [r for r in self._results if r.reverse.verdict and not r.forward.verdict]

    @property
    def no_pass(self) -> list[SequenceResult]:
        return [r for r in self._results if not r.partial_pass]

    def filter(
        self,
        require_forward: bool = False,
        require_reverse: bool = False,
        require_both: bool = False,
    ) -> list[SequenceResult]:
        """Flexible filter for downstream pipeline steps."""
        return [
            r
            for r in self._results
            if (not require_forward or r.forward.verdict)
            and (not require_reverse or r.reverse.verdict)
            and (not require_both or r.full_pass)
        ]

    def get_sequences(self, results: list[SequenceResult]) -> list[str]:
        """Extract source sequences from a filtered result list."""
        return [r.source_sequence for r in results]

    def get_ids(self, results: list[SequenceResult]) -> list[str]:
        """Extract source IDs from a filtered result list."""
        return [r.source_id for r in results]

    def __repr__(self):
        return (
            f"\n{'Batch Summary ':=^40}\n"
            f"Total processed : {self.total}\n"
            f"Full pass       : {len(self.full_pass)}\n"
            f"Forward only    : {len(self.forward_only)}\n"
            f"Reverse only    : {len(self.reverse_only)}\n"
            f"No pass         : {len(self.no_pass)}\n"
            f"{'=' * 40}"
        )

from __future__ import annotations

from genomekit.modules.gc_calculator import GCBatchSummary, GCCalculator, GCResult
from genomekit.modules.primer_finder import (
    BatchSummary,
    MeltTempMethod,
    PrimerAnalyser,
    SequenceResult,
)


class GenomeKit:
    """
    Coordinator class that validates DNA sequences and delegates to analysis tools.

    Two modes of operation:
        Single:  GenomeKit(sequence)
        Batch:   GenomeKit.batch(sequences, ids=my_ids)

    Configuration is passed per call, not at construction time.
    """

    _VALID_BASES = {"A", "T", "C", "G"}

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    def __init__(self, sequence: str) -> None:
        """
        Initialise GenomeKit in single sequence mode.

        Validates that the input is a non-empty string containing only
        A, T, C, G characters. Length validation is deferred to the
        analysis methods.

        Args:
            sequence: A DNA sequence string.

        Raises:
            TypeError:  If sequence is not a string.
            ValueError: If sequence is empty or contains invalid characters.
        """
        if not isinstance(sequence, str):
            raise TypeError(f"sequence must be a string, got {type(sequence).__name__!r}")
        if len(sequence) == 0:
            raise ValueError("sequence must not be empty")

        invalid = set(sequence.upper()) - self._VALID_BASES
        if invalid:
            raise ValueError(
                f"sequence contains invalid characters: {sorted(invalid)}. "
                f"Only A, T, C, G are permitted."
            )

        self.sequence = sequence.upper()
        self._mode = "single"

        # IDs and sequence list unused in single mode
        self._sequences: list[str] = []
        self._ids: list[str] | None = None

    @classmethod
    def batch(
        cls,
        sequences: list[str],
        ids: list[str] | None = None,
    ) -> GenomeKit:
        """
        Initialise GenomeKit in batch mode.

        Each sequence is validated individually. IDs are optional — if
        omitted they are auto-generated at analysis time as seq_0001,
        seq_0002, etc.

        Args:
            sequences: A list of DNA sequence strings.
            ids:       An optional list of identifiers, one per sequence.

        Returns:
            A GenomeKit instance in batch mode.

        Raises:
            TypeError:  If sequences is not a list.
            ValueError: If sequences is empty, or if ids length does not
                        match sequences length, or if any sequence is
                        invalid.
        """
        if not isinstance(sequences, list):
            raise TypeError(f"sequences must be a list, got {type(sequences).__name__!r}")
        if len(sequences) == 0:
            raise ValueError("sequences list must not be empty")
        if ids is not None and len(ids) != len(sequences):
            raise ValueError(
                f"ids length ({len(ids)}) must match sequences length ({len(sequences)})"
            )

        # Validate every sequence up front so the user gets a clear error
        # before any processing begins
        for i, seq in enumerate(sequences):
            label = ids[i] if ids else f"index {i}"
            if not isinstance(seq, str) or len(seq) == 0:
                raise ValueError(f"Sequence at {label!r} must be a non-empty string")
            invalid = set(seq.upper()) - cls._VALID_BASES
            if invalid:
                raise ValueError(
                    f"Sequence at {label!r} contains invalid characters: {sorted(invalid)}"
                )

        # Bypass __init__ validation by constructing a blank instance
        # and setting attributes directly
        instance = object.__new__(cls)
        instance.sequence = ""
        instance._mode = "batch"
        instance._sequences = [seq.upper() for seq in sequences]
        instance._ids = ids

        return instance

    # ------------------------------------------------------------------
    # Public analysis methods
    # ------------------------------------------------------------------

    def find_primers(
        self,
        primer_length: int = 20,
        gc_range: tuple[float, float] = (40.0, 60.0),
        tm_range: tuple[float, float] = (50.0, 65.0),
        tm_method: MeltTempMethod = MeltTempMethod.WALLACE,
    ) -> SequenceResult:
        """
        Find forward and reverse primer candidates in the sequence.

        Available in single mode only.

        Args:
            primer_length: Number of bases to evaluate at each end (default 20).
            gc_range:      Acceptable GC content range as (min, max) percent
                           (default (40.0, 60.0)).
            tm_range:      Acceptable melting temperature range as (min, max)
                           in degrees Celsius (default (50.0, 65.0)).
            tm_method:     Melting temperature calculation method
                           (default MeltTempMethod.WALLACE).

        Returns:
            A SequenceResult containing forward and reverse PrimerAnalysis
            objects, each with verdict, gc_content, melting_temp, and
            has_hairpin attributes.

        Raises:
            RuntimeError: If called on a batch-mode instance.
            ValueError:   If the sequence is too short for the given
                          primer_length.
        """
        self._require_mode("single", "find_primers")

        analyser = PrimerAnalyser(
            primer_length=primer_length,
            gc_range=gc_range,
            tm_range=tm_range,
            tm_method=tm_method,
        )
        return analyser.analyse(self.sequence, source_id="single")

    def analyse(
        self,
        primer_length: int = 20,
        gc_range: tuple[float, float] = (40.0, 60.0),
        tm_range: tuple[float, float] = (50.0, 65.0),
        tm_method: MeltTempMethod = MeltTempMethod.WALLACE,
    ) -> BatchSummary:
        """
        Run primer analysis across all sequences in the batch.

        Available in batch mode only. Returns a BatchSummary immediately
        without triggering processing — results are materialised lazily
        when the user accesses filters such as .full_pass or .forward_only.

        Args:
            primer_length: Number of bases to evaluate at each end (default 20).
            gc_range:      Acceptable GC content range as (min, max) percent
                           (default (40.0, 60.0)).
            tm_range:      Acceptable melting temperature range as (min, max)
                           in degrees Celsius (default (50.0, 65.0)).
            tm_method:     Melting temperature calculation method
                           (default MeltTempMethod.WALLACE).

        Returns:
            A BatchSummary wrapping a lazy stream of SequenceResult objects.

        Raises:
            RuntimeError: If called on a single-mode instance.
        """
        self._require_mode("batch", "analyse")

        analyser = PrimerAnalyser(
            primer_length=primer_length,
            gc_range=gc_range,
            tm_range=tm_range,
            tm_method=tm_method,
        )
        stream = analyser.analyse_batch(self._sequences, ids=self._ids)
        return BatchSummary(stream)

    def gc_content(self, high_gc_threshold: int = 60) -> GCResult:
        """
        Return GC statistics for the sequence.

        Available in single mode only.

        Args:
            high_gc_threshold: The user's threshold for high GC content

        Returns:
            A GCResult that contains: gc_content, whether that passes a desired threshold,
            at_content and gc_ratio

        Raises:
            RuntimeError: If called on a batch-mode instance.
        """
        self._require_mode("single", "gc_content")

        analyser = GCCalculator(high_gc_threshold=high_gc_threshold)

        return analyser.analyse(self.sequence, source_id="single")

    def gc_content_batch(self, high_gc_content: int = 60) -> GCBatchSummary:
        """
        Run a GC content analysis across a batch of sequences

        A GCBatchSummary report is created that lazy loads the results when filters
        are accessed e.g. sequences, ids

        Args:
            high_gc_content: The user's threshold for high GC content

        Returns:
            A GCBatchSummary report that contains: Total sequences and those with `high GC`

        Raises:
            RuntimeError: If called on single_mode instance
        """
        self._require_mode("batch", "gc_content")

        analyser = GCCalculator(high_gc_threshold=high_gc_content)

        stream = analyser.analyse_batch(self._sequences, ids=self._ids)
        return GCBatchSummary(stream)

    # ------------------------------------------------------------------
    # Dunder methods
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        if self._mode == "single":
            return (
                f"GenomeKit(mode='single', "
                f"sequence={self.sequence[:20]!r}{'...' if len(self.sequence) > 20 else ''})"
            )
        return (
            f"GenomeKit(mode='batch', "
            f"sequences={len(self._sequences)}, "
            f"ids={'provided' if self._ids else 'auto'})"
        )

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _require_mode(self, required: str, method: str) -> None:
        """
        Guard method calls against the wrong mode.

        Args:
            required: The mode this method requires ('single' or 'batch').
            method:   The name of the calling method, used in the error message.

        Raises:
            RuntimeError: If the current mode does not match required.
        """
        if self._mode != required:
            other = "batch" if required == "single" else "single"
            suggestion = "analyse()" if required == "batch" else "find_primers()"
            raise RuntimeError(
                f"{method}() is not available in {other} mode. "
                f"Use {suggestion} instead, or construct GenomeKit "
                f"with the appropriate constructor."
            )

"""The sequence model every other module reads and writes.

Coordinates are 0-based and half-open: ``Segment(start, end)`` holds bases ``start`` to
``end - 1``. A segment across the origin of a circular record ends past the record's length.
File formats convert at their own boundary; see ``docs/adr/0001-coordinates.md``.
"""

from collections.abc import Mapping
from dataclasses import KW_ONLY, dataclass, field
from enum import IntEnum
from typing import Literal

#: The IUPAC nucleotide codes a DNA sequence may hold.
IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")

_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)

type Topology = Literal["linear", "circular"]


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of an IUPAC DNA sequence, keeping its case.

    Examples
    --------
    >>> reverse_complement("GGTCTC")
    'GAGACC'
    """
    return sequence.translate(_COMPLEMENT)[::-1]


class Strand(IntEnum):
    """The strand a feature or primer binding site lies on."""

    FORWARD = 1
    REVERSE = -1
    NONE = 0
    #: SnapGene's bidirectional.
    BOTH = 2


@dataclass(frozen=True, slots=True)
class Segment:
    """One contiguous span of a feature.

    Parameters
    ----------
    start, end
        0-based, half-open. `end` passes the length of a circular record when the segment
        crosses the origin.
    name
        A label for this segment alone, such as a promoter's ``"-10"``.
    color
        ``"#rrggbb"``, or ``None`` to take the feature's colour.

    Raises
    ------
    ValueError
        Unless ``0 <= start < end``.
    """

    start: int
    end: int
    _: KW_ONLY
    name: str = ""
    color: str | None = None

    def __post_init__(self) -> None:
        """Refuse an empty or negative span."""
        if not 0 <= self.start < self.end:
            raise ValueError(f"need 0 <= start < end, got start={self.start}, end={self.end}")


@dataclass(frozen=True, slots=True)
class Feature:
    """A named, typed annotation over one or more segments.

    Parameters
    ----------
    name
        The label shown on a map.
    type
        A GenBank feature key, such as ``"CDS"`` or ``"promoter"``.
    segments
        In top-strand order, whatever the strand.
    strand
        The strand the feature reads along.
    qualifiers
        GenBank-style qualifiers, each holding one or more values.
    color
        ``"#rrggbb"``, or ``None`` for the writer's default.

    Raises
    ------
    ValueError
        If `segments` is empty.
    """

    name: str
    type: str
    segments: tuple[Segment, ...]
    _: KW_ONLY
    strand: Strand = Strand.NONE
    qualifiers: Mapping[str, tuple[str | int, ...]] = field(default_factory=dict, hash=False)
    color: str | None = None

    def __post_init__(self) -> None:
        """Refuse a feature with no segment."""
        if not self.segments:
            raise ValueError(f"feature {self.name!r} has no segment")


@dataclass(frozen=True, slots=True)
class SequenceRecord:
    """A DNA sequence with its topology and features.

    Parameters
    ----------
    sequence
        Stored upper-case.
    topology
        ``"linear"`` or ``"circular"``.
    features
        Every segment must fit the sequence under its topology.

    Raises
    ------
    ValueError
        If `sequence` holds a letter outside `IUPAC_DNA`, `topology` is unknown, or a segment
        does not fit.
    """

    sequence: str
    _: KW_ONLY
    topology: Topology = "linear"
    features: tuple[Feature, ...] = ()

    def __post_init__(self) -> None:
        """Upper-case the sequence and check every segment fits."""
        if self.topology not in ("linear", "circular"):
            raise ValueError(f"topology must be 'linear' or 'circular', got {self.topology!r}")
        sequence = self.sequence.upper()
        if bad := set(sequence) - IUPAC_DNA:
            raise ValueError(f"not IUPAC DNA: {''.join(sorted(bad))}")
        object.__setattr__(self, "sequence", sequence)
        for feature in self.features:
            for segment in feature.segments:
                self._check_fits(segment, f"feature {feature.name!r}")

    def __len__(self) -> int:
        """Return the number of bases."""
        return len(self.sequence)

    def extract(self, span: Feature | Segment) -> str:
        """Return the bases under a feature or segment, reading across the origin.

        A reverse-strand feature reads as the reverse complement of its joined segments. A bare
        segment reads off the top strand.

        Raises
        ------
        ValueError
            If a segment does not fit this record.

        Examples
        --------
        >>> SequenceRecord("AACCGGTTAC", topology="circular").extract(Segment(8, 12))
        'ACAA'
        """
        segments = span.segments if isinstance(span, Feature) else (span,)
        bases = "".join(self._bases(segment) for segment in segments)
        if isinstance(span, Feature) and span.strand == Strand.REVERSE:
            return reverse_complement(bases)
        return bases

    def _bases(self, segment: Segment) -> str:
        self._check_fits(segment, "extracted")
        n = len(self.sequence)
        if segment.end <= n:
            return self.sequence[segment.start : segment.end]
        return self.sequence[segment.start :] + self.sequence[: segment.end - n]

    def _check_fits(self, segment: Segment, owner: str) -> None:
        n = len(self.sequence)
        span = f"{owner} segment {segment.start}-{segment.end}"
        if self.topology == "linear":
            if segment.end > n:
                raise ValueError(f"{span} runs past the end of a linear record of {n} bases")
        elif segment.start >= n or segment.end - segment.start > n:
            raise ValueError(f"{span} does not fit a circular record of {n} bases")

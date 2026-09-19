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


def _dna(sequence: str) -> str:
    upper = sequence.upper()
    if bad := set(upper) - IUPAC_DNA:
        raise ValueError(f"not IUPAC DNA: {''.join(sorted(bad))}")
    return upper


def _check_span(start: int, end: int) -> None:
    if not 0 <= start < end:
        raise ValueError(f"need 0 <= start < end, got start={start}, end={end}")


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
        _check_span(self.start, self.end)


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
        In the order the feature reads them along the top strand, so a reverse-strand feature
        reads them last to first. One across the origin, or cut apart at the two ends of a linear
        record, may list a higher start first.
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
class BindingSite:
    """Where a primer's 3' part anneals to a record.

    Parameters
    ----------
    start, end
        0-based, half-open, as for `Segment`. The primer's 5' tail lies outside.
    strand
        `Strand.FORWARD` when the primer reads along the top strand, towards higher
        coordinates; `Strand.REVERSE` when it reads along the bottom strand.

    Raises
    ------
    ValueError
        Unless ``0 <= start < end`` and `strand` is forward or reverse.
    """

    start: int
    end: int
    strand: Strand

    def __post_init__(self) -> None:
        """Refuse an empty span or a strand that is neither forward nor reverse."""
        _check_span(self.start, self.end)
        if self.strand not in (Strand.FORWARD, Strand.REVERSE):
            raise ValueError(f"a binding site needs a forward or reverse strand, got {self.strand}")


@dataclass(frozen=True, slots=True)
class Primer:
    """A named oligonucleotide.

    Parameters
    ----------
    name
        The name it is ordered under.
    sequence
        5' to 3', tail included. Stored upper-case.
    binding_sites
        Where it anneals to the record that carries it.
    description
        Free text.

    Raises
    ------
    ValueError
        If `sequence` holds a letter outside `IUPAC_DNA`.
    """

    name: str
    sequence: str
    _: KW_ONLY
    binding_sites: tuple[BindingSite, ...] = ()
    description: str = ""

    def __post_init__(self) -> None:
        """Upper-case and check the sequence."""
        object.__setattr__(self, "sequence", _dna(self.sequence))


@dataclass(frozen=True, slots=True)
class SequenceRecord:
    """A DNA sequence with its topology, features, primers and notes.

    Parameters
    ----------
    sequence
        Stored upper-case.
    topology
        ``"linear"`` or ``"circular"``.
    name
        The name shown on a map.
    features, primers
        Every feature segment and primer binding site must fit the sequence under its topology.
    notes
        Descriptive fields, such as ``"Description"``, keyed by name.
    extras
        Format-specific data a reader keeps for its writer. Ignored by ``==``.

    Raises
    ------
    ValueError
        If `sequence` holds a letter outside `IUPAC_DNA`, `topology` is unknown, or a span does
        not fit.
    """

    sequence: str
    _: KW_ONLY
    topology: Topology = "linear"
    name: str = ""
    features: tuple[Feature, ...] = ()
    primers: tuple[Primer, ...] = ()
    notes: Mapping[str, str] = field(default_factory=dict, hash=False)
    extras: Mapping[str, object] = field(default_factory=dict, hash=False, compare=False)

    def __post_init__(self) -> None:
        """Upper-case the sequence and check every span fits."""
        if self.topology not in ("linear", "circular"):
            raise ValueError(f"topology must be 'linear' or 'circular', got {self.topology!r}")
        object.__setattr__(self, "sequence", _dna(self.sequence))
        for feature in self.features:
            for segment in feature.segments:
                self._check_fits(segment.start, segment.end, f"feature {feature.name!r}")
        for primer in self.primers:
            for site in primer.binding_sites:
                self._check_fits(site.start, site.end, f"primer {primer.name!r}")

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
        self._check_fits(segment.start, segment.end, "extracted")
        n = len(self.sequence)
        if segment.end <= n:
            return self.sequence[segment.start : segment.end]
        return self.sequence[segment.start :] + self.sequence[: segment.end - n]

    def _check_fits(self, start: int, end: int, owner: str) -> None:
        n = len(self.sequence)
        span = f"{owner} span {start}-{end}"
        if self.topology == "linear":
            if end > n:
                raise ValueError(f"{span} runs past the end of a linear record of {n} bases")
        elif start >= n or end - start > n:
            raise ValueError(f"{span} does not fit a circular record of {n} bases")

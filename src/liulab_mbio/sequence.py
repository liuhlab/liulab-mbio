"""The sequence model every other module reads and writes.

Coordinates are 0-based and half-open: ``Segment(start, end)`` holds bases ``start`` to
``end - 1``. A segment across the origin of a circular record ends past the record's length,
and a feature's segments after the origin start past it. File formats convert at their own
boundary; see ``docs/adr/0001-coordinates.md``.
"""

import dataclasses
from collections.abc import Iterable, Mapping
from dataclasses import KW_ONLY, dataclass, field
from enum import IntEnum
from itertools import pairwise
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
        reads them last to first. On a circular record they start in ascending order: a segment
        after the origin starts past the record's length, as one across it ends past it. A linear
        record has no origin to cross, so a feature cut apart at its two ends lists the higher
        start first.
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

    def counted_round(self, length: int) -> "Feature":
        """Return this feature as a circular record of `length` bases holds it.

        Each segment keeps its width and its place in the reading order, starts where it lies in
        the record, and is then counted round past the length if it comes after the origin.

        Examples
        --------
        >>> Feature("f", "CDS", (Segment(91, 98), Segment(2, 8))).counted_round(100).segments[1]
        Segment(start=102, end=108, name='', color=None)
        """
        starts = counted_round((one.start % length for one in self.segments), length)
        return dataclasses.replace(
            self,
            segments=tuple(
                dataclasses.replace(one, start=start, end=start + one.end - one.start)
                for start, one in zip(starts, self.segments, strict=True)
            ),
        )


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
        Every feature and primer binding site must fit the sequence under its topology. On a
        circular record a feature's first segment fits, the rest follow it in ascending order,
        and the whole spans at most one turn.
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
            self._check_feature(feature)
        for primer in self.primers:
            for site in primer.binding_sites:
                self._check_fits(site.start, site.end, f"primer {primer.name!r}")

    def __len__(self) -> int:
        """Return the number of bases."""
        return len(self.sequence)

    def extract(self, span: Feature | Segment) -> str:
        """Return the bases under a feature or segment, reading across the origin.

        A reverse-strand feature reads as the reverse complement of its joined segments. A bare
        segment reads off the top strand. Each segment reads as `bases` reads it.

        Raises
        ------
        ValueError
            If the feature or segment does not fit this record, as an annotation must.

        Examples
        --------
        >>> SequenceRecord("AACCGGTTAC", topology="circular").extract(Segment(8, 12))
        'ACAA'
        """
        if isinstance(span, Feature):
            self._check_feature(span)
        else:
            self._check_fits(span.start, span.end, "extracted")
        segments = span.segments if isinstance(span, Feature) else (span,)
        bases = "".join(self.bases(segment.start, segment.end) for segment in segments)
        if isinstance(span, Feature) and span.strand == Strand.REVERSE:
            return reverse_complement(bases)
        return bases

    def bases(self, start: int, end: int) -> str:
        """Return the top-strand bases in ``[start, end)``, reading across the origin.

        On a circular record a position counts round the circle, so `start` may lie outside the
        record; the span is at most one turn long.

        Raises
        ------
        ValueError
            If `end` comes before `start`, or the span runs off a linear record or is longer than
            a circular one.

        Examples
        --------
        >>> SequenceRecord("AACCGGTTAC", topology="circular").bases(-2, 2)
        'ACAA'
        """
        n = len(self.sequence)
        span = f"span {start}-{end}"
        if end < start:
            raise ValueError(f"{span} ends before it starts")
        if self.topology == "linear":
            if start < 0 or end > n:
                raise ValueError(f"{span} runs past the end of a linear record of {n} bases")
            return self.sequence[start:end]
        if end - start > n:
            raise ValueError(f"{span} is longer than a circular record of {n} bases")
        if start == end:
            return ""
        first = start % n
        last = first + end - start
        if last <= n:
            return self.sequence[first:last]
        return self.sequence[first:] + self.sequence[: last - n]

    def covers(
        self, span: Feature | Segment | BindingSite, inner: int | Segment | BindingSite
    ) -> bool:
        """Whether `span` holds a position, or the whole of a shorter span, across the origin too.

        A feature holds what one of its segments holds. On a circular record a position counts
        round the circle, so one past the length names the base it reaches there.

        Examples
        --------
        >>> SequenceRecord("AACCGGTTAC", topology="circular").covers(Segment(8, 12), 1)
        True
        """
        segments = span.segments if isinstance(span, Feature) else (span,)
        first, last = (inner, inner + 1) if isinstance(inner, int) else (inner.start, inner.end)
        n = len(self.sequence)
        if self.topology == "linear":
            return any(one.start <= first and last <= one.end for one in segments)
        return any(
            (first - one.start) % n + last - first <= one.end - one.start for one in segments
        )

    def fits(self, start: int, end: int) -> bool:
        """Whether ``[start, end)`` lies on this record, as every span it annotates must.

        On a linear record it lies inside the bases. On a circular one it starts inside them and
        is at most one turn long, ending past the length across the origin. An empty span fits
        wherever a base could be inserted.

        Examples
        --------
        >>> SequenceRecord("AACCGGTTAC", topology="circular").fits(8, 12)
        True
        """
        n = len(self.sequence)
        if not 0 <= start <= end:
            return False
        if self.topology == "linear":
            return end <= n
        return start < n and end - start <= n

    def _check_feature(self, feature: Feature) -> None:
        owner, segments = f"feature {feature.name!r}", feature.segments
        if self.topology == "linear":
            for segment in segments:
                self._check_fits(segment.start, segment.end, owner)
            return
        for before, after in pairwise(segments):
            if after.start < before.start:
                raise ValueError(
                    f"{owner} lists span {after.start}-{after.end} after {before.start}-"
                    f"{before.end}: past the origin of a circular record, a span starts past "
                    "its length"
                )
        self._check_fits(segments[0].start, max(one.end for one in segments), owner)

    def _check_fits(self, start: int, end: int, owner: str) -> None:
        if self.fits(start, end):
            return
        n = len(self.sequence)
        span = f"{owner} span {start}-{end}"
        if self.topology == "linear":
            raise ValueError(f"{span} runs past the end of a linear record of {n} bases")
        raise ValueError(f"{span} does not fit a circular record of {n} bases")


def across_the_origin(span: Feature | Segment | BindingSite, length: int) -> bool:
    """Whether `span` runs on from the last of `length` bases into the first.

    A segment or binding site does when it ends past `length`. A feature does when its segments,
    counted round in reading order, do; on a linear record, such a feature lies at both ends.

    Examples
    --------
    >>> across_the_origin(Segment(8, 12), 10)
    True
    """
    if not isinstance(span, Feature):
        return span.end > length
    segments = span.segments
    starts = counted_round((one.start for one in segments), length)
    return any(
        start + one.end - one.start > length for start, one in zip(starts, segments, strict=True)
    )


def counted_round(positions: Iterable[int], length: int) -> tuple[int, ...]:
    """Return positions read in order round a circle of `length` bases, at most one turn.

    A position before the one read ahead of it lies past the origin, and moves on by `length`, as
    a span across the origin ends past the length.

    Examples
    --------
    >>> counted_round((91, 2, 5), 100)
    (91, 102, 105)
    """
    counted: list[int] = []
    for position in positions:
        counted.append(position + length if counted and position < counted[-1] else position)
    return tuple(counted)

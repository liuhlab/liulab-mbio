"""Holding one clone's sequencing results against the record it should be, region by region.

`verify` lays each result on the record, counts a base only where its result trusts it, and
gives every region the record must read true one verdict:

| Verdict | When |
| --- | --- |
| pass | every base is read by a counted base of some result, and nothing disagrees |
| fail | a disagreement every result reading it shows; a mixed base; a mixed result reading it |
| none | part of it went unread, and nothing failed; the detail says which bases |
| warn | results contradict each other at a base; a base read by more than 10 reads and fewer than 20 |

Every disagreement is placed on the record as an ADR 0001 span and named with the regions it
falls in. One outside every region is named with the record's features instead and judges
nothing, because whether it matters depends on what that feature does. The thresholds are read
in ``docs/research/sequencing-read-evidence.md``, cited by section beside each.
"""

from collections.abc import Iterable, Sequence
from dataclasses import KW_ONLY, dataclass, field
from itertools import pairwise
from typing import Literal

from mbio.checks import Check, Status, counted, worst_of
from mbio.sequence import (
    IUPAC_BASES,
    Feature,
    Segment,
    SequenceRecord,
    Strand,
    position_text,
    span_text,
)
from mbio.verification.align import Column, place
from mbio.verification.result import SequencingResult

#: A base counts from this Phred value: Applied Biosystems' "acceptable" band and GENEWIZ's
#: read-length unit, the one threshold the owners agree on (section 2.4).
QUALITY_FLOOR = 20

#: A base counts only above this many reads: GENEWIZ fails a sample at "10x coverage or less"
#: (section 3.4).
DEPTH_FLOOR = 10

#: A base read by fewer reads than this warns: "approximately 20x" supports an accurate
#: consensus, in Plasmidsaurus's and Eurofins' words (section 3.4).
DEPTH_WANTED = 20

#: The share of its trusted bases a result may misread and still read as this product, chosen
#: as ten times the error a Q20 base allows (section 2.4).
MISREAD_SHARE = 0.1

#: The qualifier each cloning plan writes on the junction features of its product. Its value
#: names what begins after the junction, reading along the top strand: the part, or the record
#: a restriction piece was cut from or a Gateway segment moved from.
JUNCTION_TAG = "mbio_junction"

#: What `JUNCTION_TAG` holds on the one junction the vector's own bases follow. The plan knows
#: which that is, so nothing has to work it out from where the product's origin was turned to.
BACKBONE = "backbone"


#: The kinds of disagreement.
type Kind = Literal["substitution", "insertion", "deletion", "mixed"]


@dataclass(frozen=True, slots=True)
class Disagreement:
    """One place the results read otherwise than the record.

    Parameters
    ----------
    kind
        A substitution, an insertion, a deletion or a mixed base.
    start, end
        Where it lies on the record, as an ADR 0001 span. An insertion is empty, at the point
        its bases go in.
    bases
        What the results read there, on the record's top strand; nothing for a deletion.
    results
        The names of the results that show it.
    regions
        The regions it falls in, in record order; empty when it falls outside every one.
    features
        The record's features it falls in, named only when it falls outside every region.
    """

    kind: Kind
    start: int
    end: int
    bases: str
    _: KW_ONLY
    results: tuple[str, ...]
    regions: tuple[str, ...] = ()
    features: tuple[str, ...] = ()

    def said(self, length: int) -> str:
        """Return what it is and where, as a person reads it on a record of `length` bases."""
        if self.kind == "insertion":
            return f"insertion after {position_text(self.start - 1, length)}"
        return f"{self.kind} at {_runs(range(self.start, self.end), length)}"


@dataclass(frozen=True, slots=True)
class Placement:
    """Where one result landed on the record, if anywhere.

    Parameters
    ----------
    result
        The result's name.
    strand
        The record's strand it reads along; ``None`` where nothing was placed, or the placement
        is withheld.
    span
        The record's bases it lies over, as an ADR 0001 span; ``None`` with `strand`.
    trusted
        Its trusted span, half-open in the result's own coordinates.
    """

    result: str
    strand: Strand | None
    span: Segment | None
    trusted: tuple[int, int]


@dataclass(frozen=True, slots=True)
class Verification:
    """What a clone's sequencing results came to against its record.

    Parameters
    ----------
    checks
        One per region, in record order, named for it. Its value is how many of the region's
        bases were read.
    result_checks
        One per result, in the order given, named for it. Its value is the share of the result's
        counted bases that disagree or do not line up.
    disagreements
        Every place the results read otherwise than the record, in record order.
    placements
        Where each result landed, in the order given.
    """

    checks: tuple[Check, ...]
    result_checks: tuple[Check, ...]
    disagreements: tuple[Disagreement, ...]
    placements: tuple[Placement, ...]

    @property
    def status(self) -> Status:
        """The worst verdict of every check, by `mbio.checks.worst_of`."""
        return worst_of((*self.checks, *self.result_checks))

    @property
    def verified(self) -> bool:
        """Whether some region was judged, every region passes and no result check fails."""
        return (
            bool(self.checks)
            and all(one.status == "pass" for one in self.checks)
            and all(one.status != "fail" for one in self.result_checks)
        )


@dataclass(slots=True, eq=False)
class _Reading:
    """One result as the judge counts it."""

    result: SequencingResult
    placement: Placement
    #: Each record base the result reads with a counted base, and whether it reads it thinly.
    read: dict[int, bool] = field(default_factory=dict)
    shown: list[tuple[Kind, int, int, str]] = field(default_factory=list)
    lined_up: int = 0
    misread: int = 0


def _counts(result: SequencingResult, index: int) -> bool:
    if result.bases[index] == "N":
        return False
    if result.quality is not None and result.quality[index] < QUALITY_FLOOR:
        return False
    return result.depth is None or result.depth[index] > DEPTH_FLOOR


def _thin(result: SequencingResult, indices: Sequence[int]) -> bool:
    return result.depth is not None and any(result.depth[i] < DEPTH_WANTED for i in indices)


def _base_kind(read: str, expected: str) -> Kind | None:
    if read not in "ACGT":
        return "mixed"
    return None if read in IUPAC_BASES[expected] else "substitution"


def _read(expected: SequenceRecord, result: SequencingResult) -> _Reading:
    """Lay one result on the record and count what it reads, before any gate."""
    reading = _Reading(result, Placement(result.name, None, None, result.trusted_span))
    first, last = result.trusted_span
    trusted = sum(_counts(result, i) for i in range(first, last))
    if not trusted:
        return reading
    laid = place(expected, result)
    if laid is None:
        reading.misread = reading.lined_up = trusted
        return reading
    reading.placement = Placement(result.name, laid.strand, laid.span, result.trusted_span)
    for column in laid.placed:
        _count(expected, reading, column)
    return reading


def _count(expected: SequenceRecord, reading: _Reading, column: Column) -> None:
    result, n = reading.result, len(expected)
    counting = [i for i in column.behind if _counts(result, i)]
    if column.start == column.end:
        reading.lined_up += len(counting)
        reading.misread += len(counting)
        if counting:
            reading.shown.append(("insertion", column.start, column.end, column.bases))
        return
    if not column.bases:
        if len(counting) == len(column.behind):
            thin = _thin(result, column.behind)
            for position in range(column.start, column.end):
                reading.read[position % n] = reading.read.get(position % n, True) and thin
            reading.shown.append(("deletion", column.start, column.end, ""))
        return
    if not counting:
        return
    reading.lined_up += 1
    position = column.start % n
    reading.read[position] = reading.read.get(position, True) and _thin(result, counting)
    kind = _base_kind(column.bases, expected.sequence[position])
    if kind is not None:
        reading.misread += 1
        reading.shown.append((kind, column.start, column.end, column.bases))


def _gated(reading: _Reading) -> bool:
    return reading.misread > MISREAD_SHARE * reading.lined_up


def _touched(kind: Kind, start: int, end: int, n: int) -> frozenset[int]:
    """Return the record bases a disagreement sits on, or an insertion between."""
    if kind == "insertion":
        return frozenset({(start - 1) % n, start % n})
    return frozenset(position % n for position in range(start, end))


def _positions(feature: Feature) -> list[int]:
    return [p for one in feature.segments for p in range(one.start, one.end)]


def _runs(positions: Sequence[int], n: int) -> str:
    """Return positions as a person reads them, consecutive ones joined into spans."""
    spans: list[list[int]] = []
    for position in positions:
        if spans and position == spans[-1][1]:
            spans[-1][1] += 1
        else:
            spans.append([position, position + 1])
    return ", ".join(
        position_text(start, n) if end - start == 1 else span_text(start, end, n)
        for start, end in spans
    )


def _nothing_trusted(result: SequencingResult) -> str:
    """Return why a result trusts no base, read off its own fields."""
    first, last = result.trusted_span
    if first == last:
        if result.quality and not any(result.quality):
            return "the trace carries no quality values"
        return "the trimmed span is empty"
    quality, depth = result.quality, result.depth
    if quality is not None and all(one < QUALITY_FLOOR for one in quality[first:last]):
        return f"no base in the trimmed span reaches quality {QUALITY_FLOOR}"
    if depth is not None and all(one <= DEPTH_FLOOR for one in depth[first:last]):
        return f"no base in the trimmed span is read by more than {DEPTH_FLOOR} reads"
    return "no base in the trimmed span is called"


def _result_check(reading: _Reading, n: int) -> Check:
    result, placement = reading.result, reading.placement
    if not reading.lined_up:
        return Check(result.name, None, 0.0, _nothing_trusted(result))
    share = reading.misread / reading.lined_up
    if _gated(reading):
        return Check(
            result.name,
            "fail",
            share,
            f"does not read as this product: {reading.misread} of {reading.lined_up} trusted "
            "bases disagree or do not line up",
        )
    if result.others:
        sequences = counted(len(result.others) + 1, "consensus sequence")
        return Check(result.name, "fail", share, f"{sequences} came back: the sample is mixed")
    span = placement.span
    assert span is not None
    strand = "reverse" if placement.strand is Strand.REVERSE else "forward"
    return Check(
        result.name,
        "pass",
        share,
        f"reads {span_text(span.start, span.end, n)} on the {strand} strand",
    )


def verify(
    expected: SequenceRecord,
    results: Sequence[SequencingResult],
    regions: Sequence[Feature],
) -> Verification:
    """Hold one clone's sequencing results against the record it should be.

    A base of a result counts only inside its trusted span, at quality 20 or better where it
    carries a quality, and read by more than 10 reads where it carries a depth: both, where it
    carries both. An ``N`` reads nothing. Bases running past either end of a linear record are
    not listed and count for nothing below.

    A result does not read as this product when more than a tenth of its counted bases are
    substitutions, mixed bases or bases inserted or left unplaced: ten times the error a Q20
    base allows. Such a result's check fails, and its disagreements and placement are withheld,
    so it reads no region; the other results are judged as before. A result carrying a second
    consensus fails as mixed, and one trusting no base carries no verdict.

    Raises
    ------
    ValueError
        If `expected` is empty, or a region does not fit it.

    Examples
    --------
    >>> record = SequenceRecord("GATTACAGGCATTCGACCTA")
    >>> region = Feature("insert", "misc_feature", (Segment(4, 16),))
    >>> verify(record, (SequencingResult("r", "GATTACAGGCATTCGACCTA"),), (region,)).verified
    True
    """
    n = len(expected)
    if not n:
        raise ValueError("a result is held against a record, and this one is empty")
    for region in regions:
        expected.extract(region)
    readings = [_read(expected, one) for one in results]
    judged = [one for one in readings if one.lined_up and not _gated(one)]
    shown: dict[tuple[Kind, int, int, str], list[_Reading]] = {}
    for reading in judged:
        for key in dict.fromkeys(reading.shown):
            shown.setdefault(key, []).append(reading)
    ordered = sorted(regions, key=lambda one: min(part.start for part in one.segments))
    held = [frozenset(p % n for p in _positions(one)) for one in ordered]
    found: list[tuple[Disagreement, frozenset[int], bool]] = []
    for (kind, start, end, bases), showing in sorted(shown.items(), key=lambda item: item[0][1:]):
        touched = _touched(kind, start, end, n)
        names = tuple(
            region.name
            for region, bases_held in zip(ordered, held, strict=True)
            if touched & bases_held
        )
        features = (
            ()
            if names
            else tuple(
                feature.name
                for feature in expected.features
                if touched & {p % n for p in _positions(feature)}
            )
        )
        disagreement = Disagreement(
            kind,
            start,
            end,
            bases,
            results=tuple(reading.result.name for reading in showing),
            regions=names,
            features=features,
        )
        covering = [reading for reading in judged if touched <= reading.read.keys()]
        shared = kind == "mixed" or all(reading in showing for reading in covering)
        found.append((disagreement, touched, shared))
    return Verification(
        tuple(_region_check(*pair, n, judged, found) for pair in zip(ordered, held, strict=True)),
        tuple(_result_check(one, n) for one in readings),
        tuple(one for one, _, _ in found),
        tuple(
            one.placement
            if one in judged
            else Placement(one.result.name, None, None, one.placement.trusted)
            for one in readings
        ),
    )


def _region_check(
    region: Feature,
    held: frozenset[int],
    n: int,
    judged: Sequence[_Reading],
    found: Sequence[tuple[Disagreement, frozenset[int], bool]],
) -> Check:
    positions = _positions(region)
    inside = [(one, shared) for one, touched, shared in found if touched & held]
    failing = [one.said(n) for one, shared in inside if shared]
    failing += [
        f"{one.result.name} is mixed"
        for one in judged
        if one.result.others and held & one.read.keys()
    ]
    contradicted = [one.said(n) for one, shared in inside if not shared]
    unread = [p for p in positions if not any(p % n in one.read for one in judged)]
    thin = [
        p
        for p in positions
        if all(one.read[p % n] for one in judged if p % n in one.read)
        and any(p % n in one.read for one in judged)
    ]
    warning = []
    if contradicted:
        warning.append(f"results disagree with each other: {', '.join(contradicted)}")
    if thin:
        warning.append(f"read by fewer than {DEPTH_WANTED} reads: {_runs(thin, n)}")
    read = float(len(positions) - len(unread))
    if failing:
        return Check(region.name, "fail", read, "; ".join(failing))
    if unread:
        unread_text = f"{counted(len(unread), 'base')} unread: {_runs(unread, n)}"
        return Check(region.name, None, read, "; ".join((unread_text, *warning)))
    if warning:
        return Check(region.name, "warn", read, "; ".join(warning))
    return Check(region.name, "pass", read, "every base read, and none disagrees")


def regions(record: SequenceRecord, names: Iterable[str] = ()) -> tuple[Feature, ...]:
    """Return the features a verification of `record` judges, in record order.

    These are the junctions tagged with `JUNCTION_TAG`, and the insert between each two
    consecutive ones, named for the part the first one's tag names. The stretch the junction
    tagged `BACKBONE` opens is the vector, and is left out. A name that repeats takes an
    ordinal. Given `names`, the features of `record` so named are returned instead.

    Raises
    ------
    ValueError
        If one of `names` names no feature of `record`.

    Examples
    --------
    >>> tags = ((10, "GFP"), (30, BACKBONE))
    >>> joins = [Feature("j", "misc_feature", (Segment(at, at + 2),),
    ...          qualifiers={JUNCTION_TAG: (part,)}) for at, part in tags]
    >>> plasmid = SequenceRecord("A" * 50, topology="circular", features=tuple(joins))
    >>> [(one.name, one.segments[0].start) for one in regions(plasmid)]
    [('j', 10), ('GFP insert', 12), ('j', 30)]
    """
    wanted = tuple(names)
    if wanted:
        missing = sorted(set(wanted) - {one.name for one in record.features})
        if missing:
            raise ValueError(f"{record.name or 'the record'} has no feature named {missing[0]!r}")
        found = (one for one in record.features if one.name in wanted)
        return tuple(sorted(found, key=_place))
    junctions = sorted(
        (one for one in record.features if JUNCTION_TAG in one.qualifiers), key=_place
    )
    return tuple(sorted((*junctions, *_inserts(record, junctions)), key=_place))


def _inserts(record: SequenceRecord, junctions: list[Feature]) -> list[Feature]:
    """Return the stretch between each two consecutive junctions, but the backbone.

    A circular record's last junction is followed by its first, one turn on. Which stretch is
    the vector is the junction's own tag and not a question of geometry: a product keeps its
    vector's origin, and a junction may straddle base 0, so where that base falls says nothing
    about which side of that junction the vector lies on.
    """
    length = len(record)
    gaps: list[tuple[int, int, Feature]] = []
    reach = 0
    for before, after in pairwise(junctions):
        reach = max(reach, _end(before))
        gaps.append((reach, after.segments[0].start, before))
    if junctions and record.topology == "circular":
        last = max(junctions, key=_end)
        gaps.append((_end(last), junctions[0].segments[0].start + length, last))
    # A stretch starting a turn on is brought back into the record.
    kept = sorted(
        (
            (start % length, end - start // length * length, opener)
            for start, end, opener in gaps
            if start < end and opener.qualifiers[JUNCTION_TAG][0] != BACKBONE
        ),
        key=lambda gap: gap[:2],
    )
    taken: set[str] = set()
    inserts = []
    for start, end, opener in kept:
        stem = f"{opener.qualifiers[JUNCTION_TAG][0]} insert"
        chosen, n = stem, 1
        while chosen in taken:
            n += 1
            chosen = f"{stem} {n}"
        taken.add(chosen)
        inserts.append(Feature(chosen, "misc_feature", (Segment(start, end),)))
    return inserts


def _end(feature: Feature) -> int:
    """Where a feature's last base ends, past the length where it runs across the origin."""
    return max(one.end for one in feature.segments)


def _place(feature: Feature) -> tuple[int, int]:
    """Order features by where each begins, as a product lists them."""
    return feature.segments[0].start, feature.segments[0].end

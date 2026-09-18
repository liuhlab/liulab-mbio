"""Edit a `SequenceRecord`, and carry what it annotates into another one.

An edit shifts the features and binding sites it does not remove. `carried` and `annealed`
put them on a different record instead, which is how a simulated product keeps the annotations
of the templates it was built from, and `ordered` sorts what a product ends up with.

The functions here know no file format. Coordinates are the model's: 0-based, half-open, and a
span across the origin of a circular record ends past the record's length.
"""

import dataclasses
from dataclasses import dataclass
from enum import IntEnum

from liulab_mbio.sequence import (
    BindingSite,
    Feature,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)


class _Fate(IntEnum):
    """What an edit did to one span, in order of how much it says."""

    KEPT = 0
    CHANGED = 1
    TRIMMED = 2
    REMOVED = 3


@dataclass(frozen=True, slots=True)
class EditReport:
    """What an edit did to the features and binding sites it did not simply shift.

    Every entry holds the feature or primer as it was **before** the edit.

    Attributes
    ----------
    trimmed
        Features that lost bases, including any that lost a whole segment.
    dropped
        Features left with no bases at all.
    changed
        Features the edit fell inside: kept, and now spanning the new bases.
    dropped_sites
        Binding sites the edit overlapped, with the primer they belong to. A primer no longer
        anneals where its bases changed, so such a site is dropped rather than trimmed.
    """

    trimmed: tuple[Feature, ...] = ()
    dropped: tuple[Feature, ...] = ()
    changed: tuple[Feature, ...] = ()
    dropped_sites: tuple[tuple[Primer, BindingSite], ...] = ()


def insert(record: SequenceRecord, position: int, bases: str) -> tuple[SequenceRecord, EditReport]:
    """Insert `bases` before `position`, shifting everything after it.

    Examples
    --------
    >>> record, report = insert(SequenceRecord("AACCGG"), 2, "TT")
    >>> record.sequence
    'AATTCCGG'
    """
    return replace(record, position, position, bases)


def delete(record: SequenceRecord, start: int, end: int) -> tuple[SequenceRecord, EditReport]:
    """Delete the bases in ``[start, end)``."""
    return replace(record, start, end, "")


def replace(
    record: SequenceRecord, start: int, end: int, bases: str
) -> tuple[SequenceRecord, EditReport]:
    """Replace the bases in ``[start, end)`` with `bases`.

    On a circular record the span may run across the origin, ending past the record's length;
    the bases that remain keep their own indices where they can, so the origin moves only when
    the edit removes it.

    Raises
    ------
    ValueError
        If the span does not fit the record under its topology.
    """
    length = len(record)
    _check_span(record, start, end)
    if record.topology == "linear":
        return _splice(record, record, start, end - start, bases, wrapped=False)
    turned = rotate(record, start)
    edited, report = _splice(turned, record, 0, end - start, bases, wrapped=True)
    origin = len(bases) + (length - end if end <= length else 0)
    return (rotate(edited, origin) if len(edited) else edited), report


def rotate(record: SequenceRecord, origin: int) -> SequenceRecord:
    """Return the circular `record` read from `origin`, which becomes its first base.

    Raises
    ------
    ValueError
        If the record is linear.

    Examples
    --------
    >>> rotate(SequenceRecord("AACCGG", topology="circular"), 2).sequence
    'CCGGAA'
    """
    if record.topology != "linear":
        length = len(record)
        start = origin % length
        return dataclasses.replace(
            record,
            sequence=record.sequence[start:] + record.sequence[:start],
            features=tuple(_turned_feature(feature, start, length) for feature in record.features),
            primers=tuple(_turned_primer(primer, start, length) for primer in record.primers),
        )
    raise ValueError("only a circular record has an origin to rotate")


def flipped(record: SequenceRecord) -> SequenceRecord:
    """Return `record` read from the other strand, features and binding sites turned with it.

    Raises
    ------
    ValueError
        If a span runs across the origin, which has no place on the other strand of a record
        this turns end for end.

    Examples
    --------
    >>> flipped(SequenceRecord("AAAACCCG")).sequence
    'CGGGTTTT'
    """
    length = len(record)
    other = {Strand.FORWARD: Strand.REVERSE, Strand.REVERSE: Strand.FORWARD}
    spans = [
        (segment.start, segment.end) for feature in record.features for segment in feature.segments
    ]
    spans += [(site.start, site.end) for primer in record.primers for site in primer.binding_sites]
    if any(end > length for _, end in spans):
        raise ValueError("a record with a span across its origin cannot be turned end for end")
    features = tuple(
        dataclasses.replace(
            feature,
            segments=tuple(
                sorted(
                    (
                        Segment(
                            length - segment.end,
                            length - segment.start,
                            name=segment.name,
                            color=segment.color,
                        )
                        for segment in feature.segments
                    ),
                    key=lambda segment: (segment.start, segment.end),
                )
            ),
            strand=other.get(feature.strand, feature.strand),
        )
        for feature in record.features
    )
    primers = tuple(
        dataclasses.replace(
            primer,
            binding_sites=tuple(
                BindingSite(length - site.end, length - site.start, other[site.strand])
                for site in primer.binding_sites
            ),
        )
        for primer in record.primers
    )
    return dataclasses.replace(
        record,
        sequence=reverse_complement(record.sequence),
        features=features,
        primers=primers,
        extras={},
    )


def carried(
    record: SequenceRecord, start: int, end: int, *, offset: int
) -> tuple[tuple[Feature, ...], tuple[Primer, ...]]:
    """Carry what `record` annotates over ``[start, end)`` into a record `offset` bases along.

    A base at index `index` of `record` stands at ``index + offset`` in that record, so a span
    lifted to the front of one takes a negative `offset`.

    A feature reaching outside the span is cut down to it, and one meeting it in two places
    keeps a segment for each. A primer is kept only where a whole binding site survives, having
    nowhere to anneal otherwise.

    `end` passes `record`'s length where the span runs across the origin of a circular one.

    Examples
    --------
    >>> feature = Feature("p", "promoter", (Segment(2, 8),))
    >>> features, primers = carried(SequenceRecord("AACCGGTT", features=(feature,)), 4, 8, offset=6)
    >>> features[0].segments
    (Segment(start=10, end=14, name='', color=None),)
    """
    features = []
    for feature in record.features:
        kept = [
            dataclasses.replace(segment, start=first + offset, end=last + offset)
            for segment in feature.segments
            for first, last in _pieces(segment.start, segment.end, start, end, record)
        ]
        if kept:
            kept.sort(key=lambda segment: (segment.start, segment.end))
            features.append(dataclasses.replace(feature, segments=tuple(kept)))
    primers = []
    for primer in record.primers:
        sites = [
            BindingSite(first + offset, last + offset, site.strand)
            for site in primer.binding_sites
            for first, last in _pieces(site.start, site.end, start, end, record)
            if last - first == site.end - site.start
        ]
        if sites:
            primers.append(dataclasses.replace(primer, binding_sites=tuple(sites)))
    return tuple(features), tuple(primers)


def ordered(record: SequenceRecord) -> SequenceRecord:
    """Return `record` with its features in position order, each one's segments in base order.

    A product is built part by part and then turned to its vector's origin, and neither step
    reorders what it moves, so the record it leaves is sorted here before anyone reads it.

    Examples
    --------
    >>> late = Feature("b", "misc_feature", (Segment(4, 6),))
    >>> early = Feature("a", "misc_feature", (Segment(0, 2),))
    >>> [one.name for one in ordered(SequenceRecord("AACCGG", features=(late, early))).features]
    ['a', 'b']
    """
    features = [
        dataclasses.replace(
            feature, segments=tuple(sorted(feature.segments, key=lambda one: (one.start, one.end)))
        )
        for feature in record.features
    ]
    features.sort(key=lambda feature: (feature.segments[0].start, feature.segments[0].end))
    return dataclasses.replace(record, features=tuple(features))


def annealed(primer: Primer, edge: int, strand: Strand) -> Primer:
    """Return `primer` annotated where it anneals on a record it primed.

    A designed primer carries the length of its annealing region but not the coordinates of a
    record it has not been used on yet. `edge` is where that region meets the rest: on the
    forward strand it begins there, on the reverse strand it ends there.

    Examples
    --------
    >>> primer = Primer("p", "GGGGAACCGG", binding_sites=(BindingSite(0, 6, Strand.FORWARD),))
    >>> site = annealed(primer, 4, Strand.FORWARD).binding_sites[0]
    >>> site.start, site.end
    (4, 10)
    """
    site = primer.binding_sites[0]
    length = site.end - site.start
    start = edge if strand is Strand.FORWARD else edge - length
    return dataclasses.replace(primer, binding_sites=(BindingSite(start, start + length, strand),))


def _pieces(
    start: int, end: int, low: int, high: int, record: SequenceRecord
) -> list[tuple[int, int]]:
    """Return where ``[start, end)`` of `record` falls inside the window ``[low, high)``.

    A span of a circular record is tried a turn either way, so one meeting the window twice —
    as a feature either side of the span an outward PCR drops does — gives a piece for each.
    """
    turns = (0, len(record), -len(record)) if record.topology == "circular" else (0,)
    found = []
    for turn in turns:
        first, last = max(start + turn, low), min(end + turn, high)
        if first < last:
            found.append((first, last))
    return sorted(found)


def _check_span(record: SequenceRecord, start: int, end: int) -> None:
    length = len(record)
    span = f"span {start}-{end}"
    if not 0 <= start <= end:
        raise ValueError(f"need 0 <= start <= end, got {span}")
    if record.topology == "linear":
        if end > length:
            raise ValueError(f"{span} runs past the end of a linear record of {length} bases")
    elif start >= length or end - start > length:
        raise ValueError(f"{span} does not fit a circular record of {length} bases")


def _turned_feature(feature: Feature, start: int, length: int) -> Feature:
    return dataclasses.replace(
        feature,
        segments=tuple(_turned_segment(segment, start, length) for segment in feature.segments),
    )


def _turned_segment(segment: Segment, start: int, length: int) -> Segment:
    turned = (segment.start - start) % length
    return dataclasses.replace(segment, start=turned, end=turned + segment.end - segment.start)


def _turned_primer(primer: Primer, start: int, length: int) -> Primer:
    sites = []
    for site in primer.binding_sites:
        turned = (site.start - start) % length
        sites.append(BindingSite(turned, turned + site.end - site.start, site.strand))
    return dataclasses.replace(primer, binding_sites=tuple(sites))


def _splice(
    record: SequenceRecord,
    reported: SequenceRecord,
    start: int,
    span: int,
    bases: str,
    *,
    wrapped: bool,
) -> tuple[SequenceRecord, EditReport]:
    """Replace ``[start, start + span)`` of `record`, reporting `reported`'s own features.

    On a circular record the edit is at the origin, so a segment across the origin meets it
    twice: once at `start` and once a length further on.
    """
    length = len(record)
    edits = [(start, start + span)]
    if wrapped:
        edits.insert(0, (start + length, start + span + length))
    delta = len(bases) - span
    features, trimmed, dropped, changed = [], [], [], []
    for feature, before in zip(record.features, reported.features, strict=True):
        segments, fate = _edited_segments(
            feature.segments, edits, delta, len(bases), length + delta
        )
        if segments:
            features.append(dataclasses.replace(feature, segments=segments))
        if not segments:
            dropped.append(before)
        elif fate is _Fate.TRIMMED or fate is _Fate.REMOVED:
            trimmed.append(before)
        elif fate is _Fate.CHANGED:
            changed.append(before)
    primers, dropped_sites = [], []
    for primer, before in zip(record.primers, reported.primers, strict=True):
        kept = []
        for site, site_before in zip(primer.binding_sites, before.binding_sites, strict=True):
            spans, fate = _edited_segments(
                (Segment(site.start, site.end),), edits, delta, len(bases), length + delta
            )
            if fate is _Fate.KEPT and spans:
                kept.append(BindingSite(spans[0].start, spans[0].end, site.strand))
            else:
                dropped_sites.append((before, site_before))
        primers.append(dataclasses.replace(primer, binding_sites=tuple(kept)))
    edited = dataclasses.replace(
        record,
        sequence=record.sequence[:start] + bases + record.sequence[start + span :],
        features=tuple(features),
        primers=tuple(primers),
    )
    return edited, EditReport(tuple(trimmed), tuple(dropped), tuple(changed), tuple(dropped_sites))


def _edited_segments(
    segments: tuple[Segment, ...],
    edits: list[tuple[int, int]],
    delta: int,
    inserted: int,
    length: int,
) -> tuple[tuple[Segment, ...], _Fate]:
    kept, worst = [], _Fate.KEPT
    for segment in segments:
        span: tuple[int, int] | None = (segment.start, segment.end)
        for edit_start, edit_end in edits:
            if span is None:
                break
            span, fate = _edited_span(span, edit_start, edit_end, delta, inserted)
            worst = max(worst, fate)
        if span is None:
            worst = max(worst, _Fate.REMOVED)
            continue
        first = span[0] % length if length else 0
        kept.append(dataclasses.replace(segment, start=first, end=first + span[1] - span[0]))
    return tuple(kept), worst


def _edited_span(
    span: tuple[int, int], start: int, end: int, delta: int, inserted: int
) -> tuple[tuple[int, int] | None, _Fate]:
    """Apply one edit to one span, in the unrolled coordinates the caller passes."""
    first, last = span
    if last <= start:
        return (first, last), _Fate.KEPT
    if first >= end:
        return (first + delta, last + delta), _Fate.KEPT
    if start <= first and last <= end:
        return None, _Fate.REMOVED
    if first <= start and end <= last:
        if last + delta <= first:
            return None, _Fate.REMOVED
        inside = inserted > 0 or (first < start and end < last)
        return (first, last + delta), _Fate.CHANGED if inside else _Fate.TRIMMED
    if first < start:
        return (first, start), _Fate.TRIMMED
    return (end + delta, last + delta), _Fate.TRIMMED

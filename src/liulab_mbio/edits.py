"""Insert, delete and replace spans of a `SequenceRecord`, shifting what it annotates.

The functions here know no file format. Coordinates are the model's: 0-based, half-open, and a
span across the origin of a circular record ends past the record's length.
"""

import dataclasses
from dataclasses import dataclass
from enum import IntEnum

from liulab_mbio.sequence import BindingSite, Feature, Primer, Segment, SequenceRecord


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

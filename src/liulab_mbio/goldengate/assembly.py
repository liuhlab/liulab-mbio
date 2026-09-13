"""Simulate a Golden Gate assembly: amplify each part, cut it, and ligate the pieces.

Every part reaches the reaction as a PCR product whose primers carry a Type IIS tail, so that
cutting the amplicon leaves the overhang the design asked for. The vector is one of those parts:
`open_vector` points its primers outward from the span the assembly replaces, so the whole
backbone amplifies and the template that templated it is removed with DpnI.

Two rules run through the module:

- **An overhang is written on the top strand, so two ends join when the strings are equal.**
  That is `liulab_mbio.sites`' rule, and it makes a junction simply the bases standing at that
  point of the product.
- **A part's `left_overhang` replaces the first bases of its own span.** A junction that keeps
  what the template already spells passes those bases back; one that changes them changes the
  product there, and nothing is added or lost either way.

Coordinates are the model's, 0-based and half-open, and a span across the origin of a circular
record ends past the record's length.
"""

import dataclasses
from collections.abc import Iterable
from dataclasses import dataclass

from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.primers import (
    Q5,
    THRESHOLDS,
    PairReport,
    Polymerase,
    Thresholds,
    design_pair,
    evaluate_pair,
)
from liulab_mbio.sequence import (
    BindingSite,
    Feature,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)
from liulab_mbio.sites import SPACER_LENGTH, EnzymeLike, primer_tail

#: Dam methylates the adenine of this site, and DpnI cuts only where it has.
DAM_SITE = "GATC"


def dam_sites(record: SequenceRecord) -> int:
    """Count the sites in `record` that DpnI can cut once Dam has methylated them.

    A plasmid grown in a Dam-positive host carries them methylated and a PCR product does not,
    which is what lets DpnI take the template away and leave the amplicon.

    Examples
    --------
    >>> dam_sites(SequenceRecord("AAGATCAA"))
    1
    """
    haystack = record.sequence
    if record.topology == "circular":
        haystack += record.sequence[: len(DAM_SITE) - 1]
    return haystack.count(DAM_SITE)


@dataclass(frozen=True, slots=True)
class Part:
    """One piece of an assembly: the PCR that makes it, and the ends the enzyme leaves on it.

    Parameters
    ----------
    name
        What its tube is labelled.
    template
        The record it is amplified from.
    span
        The template bases it contributes to the product, 0-based and half-open. `end` passes
        the template's length when the span runs across the origin.
    forward, reverse
        The primers, tails included.
    amplicon
        What the PCR makes: both tails, the span, and the template's features and primers
        carried to their new coordinates.
    left_overhang
        The bases the enzyme leaves at the amplicon's left end, which stand in for the span's
        own first bases.
    right_overhang
        The bases at its right end. They belong to the next part round the circle, so they
        reach this part only through the reverse primer's tail.
    dpni
        Whether DpnI should be used to take the template away afterwards.
    report
        What the pair scored on the template, carrying the annealing temperature and the
        extension time the PCR needs.
    """

    name: str
    template: SequenceRecord
    span: tuple[int, int]
    forward: Primer
    reverse: Primer
    amplicon: SequenceRecord
    left_overhang: str
    right_overhang: str
    dpni: bool
    report: PairReport

    @property
    def length(self) -> int:
        """Bases of amplicon, which is what a gel measures and a PCR program times."""
        return len(self.amplicon)

    @property
    def fragment_length(self) -> int:
        """Bases this part puts into the product, its own overhang counted."""
        return self.span[1] - self.span[0]


def amplify(
    template: SequenceRecord,
    enzyme: EnzymeLike,
    start: int,
    end: int,
    *,
    left_overhang: str,
    right_overhang: str,
    name: str = "",
    dpni: bool | None = None,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
    spacer: str | None = None,
    spacer_length: int = SPACER_LENGTH,
    avoid: Iterable[EnzymeLike] = (),
) -> Part:
    """Design the PCR that turns ``template[start:end]`` into one part of an assembly.

    Each primer carries a tail holding the recognition site pointing back into the part, so
    cutting the amplicon leaves `left_overhang` at its left end and `right_overhang` at its
    right. `left_overhang` stands in for the span's own first bases, so the annealing region
    begins after them.

    `end` passes the template's length when the span runs across the origin; `open_vector` is
    that case with the arithmetic done for you.

    Parameters
    ----------
    template
        The record to amplify.
    enzyme
        The Type IIS enzyme the assembly is cut with.
    start, end
        The span of `template` this part puts into the product.
    left_overhang, right_overhang
        The junctions at this part's two ends, written on the top strand.
    name
        Names the part, its amplicon and its two primers.
    dpni
        Whether DpnI removes the template. Decided from the template when not given: a circular
        one is a plasmid from a Dam-positive host, and `dam_sites` says DpnI can cut it.
    polymerase, thresholds, spacer, spacer_length, avoid
        Passed to `liulab_mbio.primers` and `liulab_mbio.sites.primer_tail`.

    Returns
    -------
    Part
        The primers, the amplicon they make, and the ends a digest will leave.

    Raises
    ------
    ValueError
        If an overhang is not one this enzyme leaves, if the span is no longer than its own
        left overhang, or if no annealing region fits either end.
    """
    one = _enzyme(enzyme)
    left, right = left_overhang.upper(), right_overhang.upper()
    anneal = start + len(left)
    if anneal >= end:
        raise ValueError(
            f"the span {start}-{end} is shorter than the overhang {left!r} standing in for "
            "its first bases"
        )
    tails = tuple(
        primer_tail(one, bases, spacer=spacer, spacer_length=spacer_length, avoid=avoid)
        for bases in (left, reverse_complement(right))
    )
    forward, reverse = design_pair(
        template,
        anneal,
        end,
        forward_tail=tails[0],
        reverse_tail=tails[1],
        forward_name=f"{name} forward".strip(),
        reverse_name=f"{name} reverse".strip(),
        polymerase=polymerase,
        thresholds=thresholds,
    )
    bases = tails[0] + template.extract(Segment(anneal, end)) + reverse_complement(tails[1])
    features, carried = _carried(template, start, end, len(tails[0]) - len(left) - start)
    placed = (
        _placed(forward, len(tails[0]), Strand.FORWARD),
        _placed(reverse, len(bases) - len(tails[1]), Strand.REVERSE),
    )
    return Part(
        name,
        template,
        (start, end),
        forward,
        reverse,
        SequenceRecord(
            bases,
            name=name or template.name,
            features=features,
            primers=carried + placed,
        ),
        left,
        right,
        (template.topology == "circular" and dam_sites(template) > 0) if dpni is None else dpni,
        evaluate_pair(forward, reverse, template, polymerase=polymerase, thresholds=thresholds),
    )


def open_vector(
    vector: SequenceRecord,
    enzyme: EnzymeLike,
    start: int,
    end: int,
    *,
    overhangs: tuple[str, str],
    name: str = "",
    dpni: bool | None = None,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
    spacer: str | None = None,
    spacer_length: int = SPACER_LENGTH,
    avoid: Iterable[EnzymeLike] = (),
) -> Part:
    """Linearise a circular vector by PCR, replacing ``vector[start:end]`` with the assembly.

    The primers face outward from that span, so the whole backbone amplifies and the plasmid
    that templated it is taken away with DpnI.

    Parameters
    ----------
    vector
        The circular record to open.
    enzyme
        The Type IIS enzyme the assembly is cut with.
    start, end
        The span the assembly replaces.
    overhangs
        The two junctions, in the vector's own coordinates: the one standing at `start`, which
        the next part round the circle brings, and the one at `end`, which the backbone itself
        begins with.
    name, dpni, polymerase, thresholds, spacer, spacer_length, avoid
        As `amplify`.

    Returns
    -------
    Part
        The backbone, ready to assemble.

    Raises
    ------
    ValueError
        If `vector` is not circular, or for any reason `amplify` refuses.
    """
    if vector.topology != "circular":
        raise ValueError("a vector is opened outward across its origin, so it must be circular")
    at_start, at_end = overhangs
    return amplify(
        vector,
        enzyme,
        end,
        start + len(vector),
        left_overhang=at_end,
        right_overhang=at_start,
        name=name,
        dpni=dpni,
        polymerase=polymerase,
        thresholds=thresholds,
        spacer=spacer,
        spacer_length=spacer_length,
        avoid=avoid,
    )


def _enzyme(enzyme: EnzymeLike) -> Enzyme:
    """Read one enzyme by name or by record."""
    return get_enzyme(enzyme) if isinstance(enzyme, str) else enzyme


def _placed(primer: Primer, edge: int, strand: Strand) -> Primer:
    """Put `primer` where it anneals on a record whose copy of it ends, or begins, at `edge`."""
    annealed = primer.binding_sites[0]
    length = annealed.end - annealed.start
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


def _carried(
    record: SequenceRecord, start: int, end: int, offset: int
) -> tuple[tuple[Feature, ...], tuple[Primer, ...]]:
    """Carry what `record` annotates over ``[start, end)`` into a record `offset` bases along.

    A feature reaching outside the span is cut down to it, and one meeting it in two places
    keeps a segment for each. A primer is kept only where a whole binding site survives, having
    nowhere to anneal otherwise.
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

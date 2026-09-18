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
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from liulab_mbio.checks import Check, Status, worst_of
from liulab_mbio.edits import rotate
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.primers.design import design_pair
from liulab_mbio.primers.evaluation import PairReport, evaluate_pair
from liulab_mbio.primers.polymerase import Q5, Polymerase
from liulab_mbio.primers.thresholds import THRESHOLDS, Thresholds
from liulab_mbio.sequence import (
    BindingSite,
    Feature,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)
from liulab_mbio.sites import (
    SPACER_LENGTH,
    EnzymeLike,
    Fragment,
    digest,
    find_sites,
    primer_tail,
)

#: Dam methylates the adenine of this site, and DpnI cuts only where it has.
DAM_SITE = "GATC"

#: What a junction is drawn in. A feature built in code has no colour of its own, and
#: `liulab_mbio.snapgene` writes SnapGene's default grey for one that has none.
JUNCTION_COLOR = "#ff9900"


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

    @property
    def bases(self) -> str:
        """What this part puts into the product, its own overhang standing first."""
        start, end = self.span
        return self.left_overhang + self.template.extract(
            Segment(start + len(self.left_overhang), end)
        )


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


@dataclass(frozen=True, slots=True)
class Junction:
    """A junction as it came out: the bases two parts' overhangs paired on in the product.

    Parameters
    ----------
    start
        0-based index of the first of those bases in the product.
    overhang
        What they spell, written on the top strand.
    before, after
        The parts either side, named as they were given.
    """

    start: int
    overhang: str
    before: str
    after: str

    @property
    def end(self) -> int:
        """Where the junction's bases end."""
        return self.start + len(self.overhang)

    @property
    def span(self) -> Segment:
        """The junction's bases, as a segment of the product."""
        return Segment(self.start, self.end)


@dataclass(frozen=True, slots=True)
class Assembly:
    """What one Golden Gate reaction makes.

    Parameters
    ----------
    product
        The circular plasmid: every part's features carried to their new coordinates, the
        primers annotated where they anneal, and each junction drawn in `JUNCTION_COLOR`.
    parts
        The parts that went in, in the order they were given.
    junctions
        Where they meet, in the product's own order.
    enzyme
        The enzyme the reaction was cut with.
    """

    product: SequenceRecord
    parts: tuple[Part, ...]
    junctions: tuple[Junction, ...]
    enzyme: Enzyme

    @property
    def junction_positions(self) -> tuple[int, ...]:
        """Where each junction begins, which is what a validation design reads across."""
        return tuple(one.start for one in self.junctions)

    @property
    def checks(self) -> tuple[Check, ...]:
        """Judge the product, as data a protocol can print.

        The sites left for the enzyme come first — one that remains would cut the product open
        again — then a check for each part counting the copies of it the product holds, then how
        many junctions spell the overhang they claim. A part is counted from its template and
        not from the digest, so the count answers for the ligation rather than repeating it.
        """
        left = len(find_sites(self.product, self.enzyme))
        checks = [
            Check(
                "sites",
                "pass" if left == 0 else "fail",
                left,
                f"{self.enzyme.name} no longer cuts the product"
                if left == 0
                else f"{self.enzyme.name} still cuts the product open",
            )
        ]
        for index, part in enumerate(self.parts):
            copies = _copies(self.product, part.bases)
            checks.append(
                Check(
                    part.name or f"part {index}",
                    "pass" if copies == 1 else "fail",
                    copies,
                    f"{part.fragment_length} bases, whole and once"
                    if copies == 1
                    else f"{copies} whole copies of its {part.fragment_length} bases",
                )
            )
        matched = sum(self.product.extract(one.span) == one.overhang for one in self.junctions)
        checks.append(
            Check(
                "junctions",
                "pass" if matched == len(self.junctions) else "fail",
                matched,
                ", ".join(f"{one.overhang} at {one.start}" for one in self.junctions),
            )
        )
        return tuple(checks)

    @property
    def status(self) -> Status:
        """The worst status of any check."""
        return worst_of(self.checks)

    def __getitem__(self, name: str) -> Check:
        """Return the check of that name.

        Raises
        ------
        KeyError
            If no check has it.
        """
        for check in self.checks:
            if check.name == name:
                return check
        raise KeyError(name)


def assemble(parts: Sequence[Part], enzyme: EnzymeLike, *, name: str = "") -> Assembly:
    """Cut every part with `enzyme` and ligate them into one circular product.

    Two ends join where their overhangs are equal, both being written on the top strand, so the
    parts are chained from the first one round until the circle closes. That first part sets the
    origin: the product is turned so its template's own first base keeps the place it had, which
    leaves the vector's coordinates readable and keeps a junction off base zero.

    Parameters
    ----------
    parts
        Two or more, the linearised vector first.
    enzyme
        The Type IIS enzyme, which must be the one the parts carry tails for.
    name
        What to call the product.

    Returns
    -------
    Assembly
        The product, the parts that made it, and the junctions between them.

    Raises
    ------
    ValueError
        If fewer than two parts are given, if cutting one does not leave the fragment its tails
        were designed for, or if the overhangs do not chain into one circle.
    """
    one = _enzyme(enzyme)
    if len(parts) < 2:
        raise ValueError(f"an assembly joins at least two parts, got {len(parts)}")
    order = _chain(parts)
    bases = ""
    features: list[Feature] = []
    primers: list[Primer] = []
    joins: list[tuple[int, Part, Part]] = []
    for place, index in enumerate(order):
        part, piece = parts[index], _cut(parts[index], one)
        at = len(bases)
        bases += part.amplicon.extract(Segment(piece.start, piece.end))
        carried, kept = _carried(part.template, *part.span, at - part.span[0])
        features.extend(carried)
        primers.extend(kept)
        primers.append(_placed(part.forward, at + len(part.left_overhang), Strand.FORWARD))
        primers.append(_placed(part.reverse, at + part.fragment_length, Strand.REVERSE))
        joins.append((at, parts[order[place - 1]], part))
    features.extend(_junction_feature(at, before, after, one) for at, before, after in joins)
    origin = _origin(parts[order[0]])
    product = SequenceRecord(
        bases, topology="circular", name=name, features=tuple(features), primers=tuple(primers)
    )
    return Assembly(
        _ordered(rotate(product, origin) if origin else product),
        tuple(parts),
        tuple(
            sorted(
                (
                    Junction(
                        (at - origin) % len(bases), after.left_overhang, before.name, after.name
                    )
                    for at, before, after in joins
                ),
                key=lambda junction: junction.start,
            )
        ),
        one,
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


def _copies(record: SequenceRecord, bases: str) -> int:
    """Count where `bases` reads whole in `record`, on either strand and across the origin."""
    if not bases or len(bases) > len(record):
        return 0
    haystack = record.sequence
    if record.topology == "circular":
        haystack += record.sequence[: len(bases) - 1]
    return _occurrences(haystack, bases) + _occurrences(haystack, reverse_complement(bases))


def _occurrences(haystack: str, needle: str) -> int:
    """Count where `needle` reads in `haystack`, overlapping copies counted separately."""
    found, at = 0, haystack.find(needle)
    while at >= 0:
        found, at = found + 1, haystack.find(needle, at + 1)
    return found


def _ordered(record: SequenceRecord) -> SequenceRecord:
    """Put `record`'s features in position order, each one's segments in top-strand order.

    Turning a record moves its segments without reordering them, and parts are joined in the
    order they ligate rather than the order they end up in.
    """
    features = [
        dataclasses.replace(
            feature, segments=tuple(sorted(feature.segments, key=lambda one: (one.start, one.end)))
        )
        for feature in record.features
    ]
    features.sort(key=lambda feature: (feature.segments[0].start, feature.segments[0].end))
    return dataclasses.replace(record, features=tuple(features))


def _junction_feature(at: int, before: Part, after: Part, enzyme: Enzyme) -> Feature:
    """Draw the junction that begins at `at`, where `before` gives way to `after`."""
    overhang = after.left_overhang
    return Feature(
        f"{overhang} junction",
        "misc_feature",
        (Segment(at, at + len(overhang)),),
        color=JUNCTION_COLOR,
        qualifiers={"note": (f"{before.name} to {after.name}, {enzyme.name} overhang",)},
    )


def _cut(part: Part, enzyme: Enzyme) -> Fragment:
    """Return the one piece a digest of `part` releases with the ends it was designed for."""
    wanted = (part.left_overhang, part.right_overhang)
    pieces = [
        piece
        for piece in digest(part.amplicon, enzyme)
        if (piece.left_overhang, piece.right_overhang) == wanted
    ]
    if len(pieces) != 1:
        raise ValueError(
            f"cutting {part.name or 'a part'} with {enzyme.name} leaves {len(pieces)} pieces "
            f"ending in {wanted[0]} and {wanted[1]}, the overhangs its tails were designed for"
        )
    return pieces[0]


def _chain(parts: Sequence[Part]) -> list[int]:
    """Return the order the parts ligate in, starting from the first one given."""
    order, used = [0], {0}
    while len(order) < len(parts):
        before = parts[order[-1]]
        found = [
            index
            for index, part in enumerate(parts)
            if index not in used and part.left_overhang == before.right_overhang
        ]
        if len(found) != 1:
            raise ValueError(
                f"{len(found)} parts begin with the overhang {before.right_overhang}, which "
                f"{before.name or 'the part before'} ends with"
            )
        order.append(found[0])
        used.add(found[0])
    last, first = parts[order[-1]], parts[0]
    if last.right_overhang != first.left_overhang:
        raise ValueError(
            f"the circle does not close: {last.name or 'the last part'} ends with "
            f"{last.right_overhang} and {first.name or 'the first'} begins with "
            f"{first.left_overhang}"
        )
    return order


def _origin(part: Part) -> int:
    """Where the first part's template origin falls in the product, or 0 when it is not there."""
    start, end = part.span
    turns = (0, len(part.template)) if part.template.topology == "circular" else (0,)
    for turn in turns:
        if start <= turn < end:
            return turn - start
    return 0

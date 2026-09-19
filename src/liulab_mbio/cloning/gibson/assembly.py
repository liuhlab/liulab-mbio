"""Simulate a Gibson assembly: amplify each part with its overlaps, and join them on those.

Most parts reach the reaction as a PCR product. What a part contributes to the product is its
own span; the overlaps at its two ends are the neighbours' bases, carried as 5' tails so that
the two fragments spell the same thing where they meet and the reaction anneals them there.
`open_vector` points a circular vector's primers outward from the span the inserts replace, so
the whole backbone amplifies and DpnI takes the plasmid that templated it away. A backbone that
is already linear needs none of that and goes in as it is: `given_vector` makes the part.

Two more parts reach it without a PCR. A part short enough is **stitched** -- ordered as
overlapping oligos tiling both strands and assembled in the same reaction -- which `stitch`
makes from the oligos `liulab_mbio.cloning.gibson.design` lays out. And where two neighbours
share no homology at all, one **bridging oligo** carrying homology to both ends joins them, so
neither carries a tail: that is `Bridge`, and `Part.bridge` is where it sits.

Two rules run through the module:

- **A tail is the neighbour's sequence, never the part's own.** So a part's amplicon is longer
  than what it puts into the product, and the junction's bases are counted once.
- **Exactly one molecule carries a junction.** The part that does not already spell the bases
  carries them as a primer tail, or a bridging oligo carries them and neither part does.
  `taken_from` records which part already spelt them and `bridge` which oligo carries them.

Coordinates are the model's, 0-based and half-open, and a span across the origin of a circular
record ends past the record's length.
"""

import dataclasses
from collections.abc import Sequence
from dataclasses import dataclass

from liulab_mbio.bench.steps import dam_sites
from liulab_mbio.bench.validation import insert_order
from liulab_mbio.checks import Check, Status, worst_of
from liulab_mbio.edits import annealed, carried, ordered, rotate
from liulab_mbio.primers.design import design_pair
from liulab_mbio.primers.evaluation import PairReport, evaluate_pair
from liulab_mbio.primers.polymerase import Q5, Polymerase
from liulab_mbio.primers.thresholds import THRESHOLDS, Thresholds
from liulab_mbio.sequence import (
    Feature,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)

#: What an overlap is drawn in. A feature built in code has no colour of its own, and
#: `liulab_mbio.snapgene` writes SnapGene's default grey for one that has none.
OVERLAP_COLOR = "#ff9900"


@dataclass(frozen=True, slots=True)
class StitchOligo:
    """One oligo of the set that tiles a stitched part, on the strand it is written for.

    Parameters
    ----------
    name
        What to order it under.
    sequence
        Its bases 5' to 3', so a reverse-strand oligo reads the complement of the part.
    strand
        Which strand of the part it lies on. Neighbours alternate, so the set tiles both.
    start
        Where it begins on the part's own molecule, tails included, 0-based.
    """

    name: str
    sequence: str
    strand: Strand
    start: int

    @property
    def end(self) -> int:
        """Where it ends on that molecule."""
        return self.start + len(self.sequence)


@dataclass(frozen=True, slots=True)
class Bridge:
    """One oligo joining two parts that share no homology, so neither needs a tail.

    Parameters
    ----------
    name
        What to order it under.
    sequence
        Its bases 5' to 3': the last `homology_bp` of the part before, then the first
        `homology_bp` of the part after.
    homology_bp
        Bases it shares with each of the two.
    """

    name: str
    sequence: str
    homology_bp: int


@dataclass(frozen=True, slots=True)
class Part:
    """One fragment of an assembly: the PCR that makes it, and the overlaps it is tailed with.

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
        The primers, tails included, or ``None`` for a part handed in ready to assemble.
    amplicon
        What the PCR makes: both tails, the span, and the template's features and primers
        carried to their new coordinates. The record itself, for a part that is not amplified.
    left_tail, right_tail
        The overlaps this part carries, which belong to the part before it and the part after
        it round the product. Either is empty where that neighbour carries the overlap instead.
    dpni
        Whether DpnI should be used to take the template away afterwards.
    report
        What the pair scored on the template, carrying the annealing temperature and the
        extension time the PCR needs, or ``None`` for a part that is not amplified.
    oligos
        The oligos that tile this part where it is stitched rather than amplified, in the order
        they run along it. Empty for every other part.
    bridge
        The oligo joining this part to the one before it round the product, where the two share
        no homology. ``None`` where a primer tail carries that junction instead.
    """

    name: str
    template: SequenceRecord
    span: tuple[int, int]
    forward: Primer | None
    reverse: Primer | None
    amplicon: SequenceRecord
    left_tail: str
    right_tail: str
    dpni: bool
    report: PairReport | None
    oligos: tuple[StitchOligo, ...] = ()
    bridge: Bridge | None = None

    @property
    def length(self) -> int:
        """Bases of amplicon, which is what a gel measures and a PCR program times."""
        return len(self.amplicon)

    @property
    def fragment_length(self) -> int:
        """Bases this part puts into the product, which is fewer than its amplicon holds."""
        return self.span[1] - self.span[0]

    @property
    def bases(self) -> str:
        """What this part puts into the product, its neighbours' overlaps left off."""
        return self.template.extract(Segment(*self.span))

    @property
    def amplified(self) -> bool:
        """Whether a PCR makes this part, or it was handed in ready to assemble."""
        return self.report is not None

    @property
    def stitched(self) -> bool:
        """Whether oligos tiling both strands make this part rather than a PCR."""
        return bool(self.oligos)


def amplify(
    template: SequenceRecord,
    start: int,
    end: int,
    *,
    left_tail: str = "",
    right_tail: str = "",
    name: str = "",
    dpni: bool | None = None,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
) -> Part:
    """Design the PCR that turns ``template[start:end]`` into one part of an assembly.

    `left_tail` joins the forward primer's 5' end and `right_tail`, reverse-complemented, the
    reverse primer's, so the amplicon reads ``left_tail + template[start:end] + right_tail``.
    Both are the neighbours' bases, so neither moves where the annealing regions sit and the Tm
    of each primer is read from its annealing region alone.

    `end` passes the template's length when the span runs across the origin; `open_vector` is
    that case with the arithmetic done for you.

    Parameters
    ----------
    template
        The record to amplify.
    start, end
        The span of `template` this part puts into the product.
    left_tail, right_tail
        The overlaps this part carries, written on the top strand of the product.
    name
        Names the part, its amplicon and its two primers.
    dpni
        Whether DpnI removes the template. Decided from the template when not given: a circular
        one is a plasmid from a Dam-positive host, and `dam_sites` says DpnI can cut it.
    polymerase, thresholds
        Passed to `liulab_mbio.primers`.

    Returns
    -------
    Part
        The primers, the amplicon they make, and the overlaps at its two ends.

    Raises
    ------
    ValueError
        If the span is empty, or if no annealing region fits either end.
    """
    if start >= end:
        raise ValueError(f"a part needs bases to contribute, got the span {start}-{end}")
    left, right = left_tail.upper(), right_tail.upper()
    reverse_tail = reverse_complement(right)
    forward, reverse = design_pair(
        template,
        start,
        end,
        forward_tail=left,
        reverse_tail=reverse_tail,
        forward_name=f"{name} forward".strip(),
        reverse_name=f"{name} reverse".strip(),
        polymerase=polymerase,
        thresholds=thresholds,
    )
    bases = left + template.extract(Segment(start, end)) + right
    features, kept = carried(template, start, end, offset=len(left) - start)
    placed = (
        annealed(forward, len(left), Strand.FORWARD),
        annealed(reverse, len(bases) - len(right), Strand.REVERSE),
    )
    return Part(
        name,
        template,
        (start, end),
        forward,
        reverse,
        SequenceRecord(bases, name=name or template.name, features=features, primers=kept + placed),
        left,
        right,
        (template.topology == "circular" and dam_sites(template) > 0) if dpni is None else dpni,
        evaluate_pair(forward, reverse, template, polymerase=polymerase, thresholds=thresholds),
    )


def stitch(
    record: SequenceRecord,
    oligos: Sequence[StitchOligo],
    *,
    left_tail: str = "",
    right_tail: str = "",
    name: str = "",
) -> Part:
    """Make one part out of `oligos` tiling it, rather than out of a PCR.

    The oligos assemble into ``left_tail + record + right_tail`` in the same reaction as
    everything else, so the part they make is the part an amplicon would have made: it puts the
    same bases into the product and carries the same two overlaps. There is no template behind
    it for DpnI to take away, and no primer to judge.
    `liulab_mbio.cloning.gibson.design.stitch_oligos` lays the oligos out.

    Raises
    ------
    ValueError
        If `record` is circular, which is opened rather than stitched, or if no oligo is given.
    """
    if record.topology != "linear":
        raise ValueError(f"{record.name or 'this part'} is circular, so it cannot be stitched")
    if not oligos:
        raise ValueError(f"{record.name or 'a stitched part'} needs the oligos that tile it")
    left, right = left_tail.upper(), right_tail.upper()
    features, kept = carried(record, 0, len(record), offset=len(left))
    return Part(
        name or record.name,
        record,
        (0, len(record)),
        None,
        None,
        SequenceRecord(
            left + record.sequence + right,
            name=name or record.name,
            features=features,
            primers=kept,
        ),
        left,
        right,
        False,
        None,
        tuple(oligos),
    )


def open_vector(
    vector: SequenceRecord,
    start: int,
    end: int,
    *,
    left_tail: str = "",
    right_tail: str = "",
    name: str = "",
    dpni: bool | None = None,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
) -> Part:
    """Linearise a circular vector by PCR, replacing ``vector[start:end]`` with the inserts.

    The primers face outward from that span, so the whole backbone amplifies and the plasmid
    that templated it is taken away with DpnI. A vector whose junctions take their overlaps
    from the vector itself needs no tail at all, which is what NEB asks for where the same
    opened backbone is reused.

    Parameters
    ----------
    vector
        The circular record to open.
    start, end
        The span the inserts replace.
    left_tail, right_tail, name, dpni, polymerase, thresholds
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
    return amplify(
        vector,
        end,
        start + len(vector),
        left_tail=left_tail,
        right_tail=right_tail,
        name=name,
        dpni=dpni,
        polymerase=polymerase,
        thresholds=thresholds,
    )


def given_vector(vector: SequenceRecord, *, name: str = "") -> Part:
    """Take a vector that is already linear as the opened part, with no PCR at all.

    A backbone cut or amplified last week goes in as it is: it spells what it spells, so it
    carries no tail and its neighbours' primers carry its bases instead, and there is no
    template behind it for DpnI to take away.

    Raises
    ------
    ValueError
        If `vector` is circular, which is opened by `open_vector` instead.
    """
    if vector.topology != "linear":
        raise ValueError(
            f"{vector.name or 'this vector'} is circular, so it is opened by PCR rather than "
            "taken as the opened part"
        )
    return Part(
        name or vector.name,
        vector,
        (0, len(vector)),
        None,
        None,
        vector,
        "",
        "",
        False,
        None,
    )


@dataclass(frozen=True, slots=True)
class Junction:
    """A junction as it came out: the bases the part before it spells where the two meet.

    They are the bases something else has to anneal to. Where a primer tail carries them the
    part after spells them too, which is the overlap the two share; where a bridging oligo
    carries them it spells them and as many again of the part after, and the two parts share
    nothing.

    Parameters
    ----------
    start
        0-based index of the first of those bases in the product.
    overlap
        What they spell, written on the top strand.
    before, after
        The parts either side, named as they were given.
    taken_from
        Which of the two already spelt them.
    bridge
        The oligo carrying them, where neither part does: it spells these bases and as many
        again of the part after, so the two parts share nothing. Empty where the part that did
        not spell them carries them as a primer tail instead.
    """

    start: int
    overlap: str
    before: str
    after: str
    taken_from: str
    bridge: str = ""

    @property
    def end(self) -> int:
        """Where the junction's bases end."""
        return self.start + len(self.overlap)

    @property
    def boundary(self) -> int:
        """Where the part before gives way to the part after, which is not where the bases begin.

        The bases belong to one of the two parts, so they lie on that part's side of the
        boundary: `end` where the part before them spells them, `start` where the part after
        does. A bridged junction is the first case, the oligo reaching as far again past it.
        """
        return self.end if self.taken_from == self.before else self.start

    @property
    def length(self) -> int:
        """Bases carried here, which is what a bridging oligo takes from each side."""
        return len(self.overlap)

    @property
    def span(self) -> Segment:
        """The junction's bases, as a segment of the product."""
        return Segment(self.start, self.end)


@dataclass(frozen=True, slots=True)
class Assembly:
    """What one Gibson reaction makes.

    Parameters
    ----------
    product
        The circular plasmid: every part's features carried to their new coordinates, the
        primers annotated where they anneal, and each junction drawn in `OVERLAP_COLOR`.
    parts
        The parts that went in, in the order they go round the product.
    junctions
        Where they meet, in the product's own order.
    """

    product: SequenceRecord
    parts: tuple[Part, ...]
    junctions: tuple[Junction, ...]

    @property
    def junction_positions(self) -> tuple[int, ...]:
        """Where each junction's shared bases begin, which is not where a part gives way."""
        return tuple(one.start for one in self.junctions)

    @property
    def boundaries(self) -> tuple[int, ...]:
        """Where each part gives way to the next, which is what a validation design reads across.

        `Junction.boundary` is the rule and `insert_order` the order: from the boundary the
        first part gives way at, so the last passes the product's length where the inserts run
        across the origin.
        """
        opened = self.parts[0].name
        at = next(index for index, one in enumerate(self.junctions) if one.before == opened)
        return insert_order([one.boundary for one in self.junctions], at, len(self.product))

    @property
    def insert_span(self) -> tuple[int, int]:
        """The product bases the inserts own: everything the first part does not.

        From the boundary the first part gives way at to the one it takes over at, which are the
        two ends of `boundaries`.
        """
        boundaries = self.boundaries
        return boundaries[0], boundaries[-1]

    @property
    def checks(self) -> tuple[Check, ...]:
        """Judge the product, as data a protocol can print.

        One check for each part, counting the copies of it the product holds, then one for how
        many junctions spell the overlap they claim. A part is counted from its template, so
        the count answers for the assembly rather than repeating it.
        """
        checks = [
            Check(
                part.name or f"part {index}",
                "pass" if copies == 1 else "fail",
                copies,
                f"{part.fragment_length} bases, whole and once"
                if copies == 1
                else f"{copies} whole copies of its {part.fragment_length} bases",
            )
            for index, part in enumerate(self.parts)
            for copies in (_copies(self.product, part.bases),)
        ]
        matched = sum(self.product.extract(one.span) == one.overlap for one in self.junctions)
        checks.append(
            Check(
                "junctions",
                "pass" if matched == len(self.junctions) else "fail",
                matched,
                ", ".join(
                    f"{one.length} bp at {one.start}, "
                    + (f"bridged by {one.bridge}" if one.bridge else f"from {one.taken_from}")
                    for one in self.junctions
                ),
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


def assemble(parts: Sequence[Part], *, name: str = "") -> Assembly:
    """Join the parts into one circular product, in the order they are given.

    Each part contributes its own span, and the overlap at each junction is the bases one of
    the two parts already spells there. The first part sets the origin: the product is turned
    so its template's own first base keeps the place it had, which leaves the vector's
    coordinates readable and keeps a junction off base zero.

    Parameters
    ----------
    parts
        Two or more, the linearised vector first.
    name
        What to call the product.

    Returns
    -------
    Assembly
        The product, the parts that made it, and the junctions between them.

    Raises
    ------
    ValueError
        If fewer than two parts are given, or two neighbours share no overlap at all.
    """
    if len(parts) < 2:
        raise ValueError(f"an assembly joins at least two parts, got {len(parts)}")
    bases = ""
    features: list[Feature] = []
    primers: list[Primer] = []
    starts: list[int] = []
    for part in parts:
        at = len(bases)
        starts.append(at)
        bases += part.bases
        over, kept = carried(part.template, *part.span, offset=at - part.span[0])
        features.extend(over)
        primers.extend(kept)
        if part.forward is not None and part.reverse is not None:
            primers.append(annealed(part.forward, at, Strand.FORWARD))
            primers.append(annealed(part.reverse, at + part.fragment_length, Strand.REVERSE))
    length = len(bases)
    joins = [
        _junction(parts[index - 1], part, starts[index], length) for index, part in enumerate(parts)
    ]
    features.extend(_overlap_feature(one) for one in joins)
    product = SequenceRecord(
        bases, topology="circular", name=name, features=tuple(features), primers=tuple(primers)
    )
    origin = _origin(parts[0])
    return Assembly(
        ordered(rotate(product, origin) if origin else product),
        tuple(parts),
        tuple(
            sorted(
                (dataclasses.replace(one, start=(one.start - origin) % length) for one in joins),
                key=lambda one: one.start,
            )
        ),
    )


def _junction(before: Part, after: Part, at: int, length: int) -> Junction:
    """Return the junction where `before` gives way to `after` at `at` in the product.

    Raises
    ------
    ValueError
        If neither part carries an overlap there and no oligo bridges it, so nothing anneals
        the two.
    """
    if after.bridge is not None:
        homology = after.bridge.homology_bp
        return Junction(
            (at - homology) % length,
            before.bases[-homology:],
            before.name,
            after.name,
            before.name,
            after.bridge.name,
        )
    if after.left_tail:
        return Junction(
            (at - len(after.left_tail)) % length,
            after.left_tail,
            before.name,
            after.name,
            before.name,
        )
    if before.right_tail:
        return Junction(at, before.right_tail, before.name, after.name, after.name)
    raise ValueError(
        f"{before.name or 'a part'} and {after.name or 'the next'} share no overlap, so nothing "
        "joins them: one of the two has to carry the other's bases as a primer tail, or a "
        "bridging oligo has to span them"
    )


def _overlap_feature(junction: Junction) -> Feature:
    """Draw what carries the junction, so a map shows where it is and what holds it together."""
    if junction.bridge:
        return Feature(
            f"{junction.before}-{junction.after} bridge",
            "misc_feature",
            (Segment(junction.start, junction.end + junction.length),),
            color=OVERLAP_COLOR,
            qualifiers={
                "note": (
                    f"{junction.bridge}, {junction.length} bp of {junction.before} and as many "
                    f"of {junction.after}; neither carries a tail",
                )
            },
        )
    return Feature(
        f"{junction.before}-{junction.after} overlap",
        "misc_feature",
        (Segment(junction.start, junction.end),),
        color=OVERLAP_COLOR,
        qualifiers={
            "note": (
                f"{junction.length} bp shared with {junction.before}, taken from "
                f"{junction.taken_from}",
            )
        },
    )


def _origin(part: Part) -> int:
    """Where the first part's template origin falls in the product, or 0 when it is not there."""
    start, end = part.span
    turns = (0, len(part.template)) if part.template.topology == "circular" else (0,)
    for turn in turns:
        if start <= turn < end:
            return turn - start
    return 0


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

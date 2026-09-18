"""Simulate one Gateway recombination: the segment between two att sites moves into a backbone.

Both reactions have one shape. One record carries the segment that moves, between its two att
sites; the other is a circular vector whose backbone outside its own two sites takes it. LR
moves an entry clone's insert into a destination vector; BP moves an attB substrate into a
donor vector. What each junction spells afterwards is `liulab_mbio.cloning.gateway.att`, and a
`PlannedReaction` is one of these with the bench numbers a plan works out beside it.

The product is built rather than asserted: every feature of both records is carried to its new
coordinates, each primer is kept where it still anneals, and each att junction is marked. The
junction is not scarless -- a whole 25 bp att site stands at each end of the moved segment --
so the product is what says what it spells.

Coordinates are the model's, 0-based and half-open, and a span across the origin of a circular
record ends past the record's length.
"""

from dataclasses import dataclass

from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.phenotype import Phenotype
from liulab_mbio.checks import Check, Status, worst_of
from liulab_mbio.cloning.gateway.att import (
    CROSSOVER,
    REGION_BP,
    SUBSTRATES,
    AttSite,
    Reaction,
    att_pair,
    find_att_sites,
    identify,
    joined,
    writes,
)
from liulab_mbio.edits import carried, flipped, ordered, rotate
from liulab_mbio.sequence import Feature, Primer, Segment, SequenceRecord, Strand

#: What an att junction is drawn in. A feature built in code has no colour of its own, and
#: `liulab_mbio.snapgene` writes SnapGene's default grey for one that has none.
JUNCTION_COLOR = "#ff9900"


@dataclass(frozen=True, slots=True)
class Junction:
    """One att site a recombination wrote, where it stands in the product.

    Not `liulab_mbio.overhangs.Junction`, which is a pair of ligated ends, nor
    `liulab_mbio.cloning.goldengate.assembly.Junction`, which is the bases two overhangs paired
    on. This one is an att site, and the product gains all 25 of its bases.

    Parameters
    ----------
    name
        Which att site the arithmetic says it is, such as ``"attB1"``.
    start
        0-based index of the region's first base on the product's top strand.
    bases
        The 25 bases in the site's own orientation.
    strand
        `Strand.FORWARD` where the site reads along the product's top strand.
    before, after
        The records that gave the first `CROSSOVER` bases and the rest, named as they were
        given.
    """

    name: str
    start: int
    bases: str
    strand: Strand
    before: str
    after: str

    @property
    def end(self) -> int:
        """Where the region ends, past the product's length when it runs across the origin."""
        return self.start + REGION_BP

    @property
    def span(self) -> Segment:
        """The region, as a segment of the product."""
        return Segment(self.start, self.end)


@dataclass(frozen=True, slots=True)
class Piece:
    """What one substrate put into the product.

    Parameters
    ----------
    name
        The record's, as it was given.
    record
        The substrate, on the strand that went in.
    span
        The bases it contributed, in its own coordinates. `end` passes the record's length
        where the span runs across the origin.
    """

    name: str
    record: SequenceRecord
    span: tuple[int, int]

    @property
    def bases(self) -> str:
        """What it contributed, read off its own record."""
        return self.record.extract(Segment(*self.span))

    @property
    def length(self) -> int:
        """How many bases that is."""
        return self.span[1] - self.span[0]


@dataclass(frozen=True, slots=True)
class Recombination:
    """What one Gateway reaction makes.

    Parameters
    ----------
    reaction
        ``"BP"`` or ``"LR"``.
    product
        The circular plasmid: both records' features carried to their new coordinates, their
        primers kept where they still anneal, and each att junction drawn in `JUNCTION_COLOR`.
    moved
        The segment that moved, and the record it came from.
    backbone
        The vector bases outside its att sites, and the record they came from.
    junctions
        The two att sites the reaction wrote, site 1 first.
    """

    reaction: Reaction
    product: SequenceRecord
    moved: Piece
    backbone: Piece
    junctions: tuple[Junction, Junction]

    @property
    def boundaries(self) -> tuple[int, int]:
        """Where the backbone gives way to the moved segment and back, which validation reads across.

        Not the junctions' own edges: the first `CROSSOVER` bases of each att region are the
        acceptor's, so the segment starts inside the first region and ends inside the second.
        """
        first, second = self.junctions
        return first.start + CROSSOVER, second.end - CROSSOVER

    @property
    def cassette(self) -> tuple[int, int]:
        """The acceptor bases the reaction throws away, in the acceptor's own coordinates.

        The ccdB and chloramphenicol genes that stood between its att sites, which leave with
        the by-product. `end` passes the record's length where the span runs across the origin.
        """
        length = len(self.backbone.record)
        start = self.backbone.span[1] % length
        return start, start + (length - self.backbone.length)

    @property
    def checks(self) -> tuple[Check, ...]:
        """Judge the product, as data a protocol can print.

        Both are read off the product's own bases rather than off the arithmetic that built it:
        that each junction reads as the att site it is named for, and that the product carries
        those two att sites and no other. A third site would recombine where nobody meant it to.
        Each carries the reaction's name, so a plan running both keeps its verdicts apart.
        """
        spelled = sum(identify(one.bases) == one.name for one in self.junctions)
        found = find_att_sites(self.product)
        places = {(one.start, one.name) for one in self.junctions}
        extra = [one for one in found if (one.start, one.name) not in places]
        return (
            Check(
                f"{self.reaction} junctions",
                "pass" if spelled == len(self.junctions) else "fail",
                spelled,
                ", ".join(f"{one.name} at {one.start}" for one in self.junctions),
            ),
            Check(
                f"{self.reaction} att sites",
                "pass" if len(found) == len(self.junctions) and not extra else "fail",
                len(found),
                "the two the reaction wrote, and no other"
                if not extra
                else "also " + ", ".join(f"{one.name} at {one.start}" for one in extra),
            ),
        )

    @property
    def status(self) -> Status:
        """The worst status of any check."""
        return worst_of(self.checks)


@dataclass(frozen=True, slots=True)
class PlannedReaction:
    """One recombination as a plan holds it: what it makes, what it takes, and what grows.

    Parameters
    ----------
    recombination
        The simulated reaction.
    amounts
        What to put in it, the record whose segment moves first.
    phenotype
        What the product says about itself.
    """

    recombination: Recombination
    amounts: tuple[Amount, ...]
    phenotype: Phenotype

    @property
    def product(self) -> SequenceRecord:
        """The plasmid the reaction makes."""
        return self.recombination.product

    @property
    def junctions(self) -> tuple[Junction, Junction]:
        """The two att sites it wrote, site 1 first."""
        return self.recombination.junctions


def recombine(
    carrier: SequenceRecord, acceptor: SequenceRecord, *, reaction: Reaction, name: str = ""
) -> Recombination:
    """Recombine the segment between `carrier`'s att sites into `acceptor`'s backbone.

    Both records are read for the pair of att sites the reaction takes, on either strand, and
    each is turned so that its site 1 reads along the top strand. The product is then the
    acceptor's bases outside its own sites, an att junction, the carrier's segment, and the
    other junction.

    Parameters
    ----------
    carrier
        The record whose segment moves: an entry clone for LR, an attB substrate for BP.
    acceptor
        The circular vector whose backbone takes it: a destination vector for LR, a donor
        vector for BP.
    reaction
        Which reaction to run, which is what fixes the kind of site looked for in each record.
    name
        What to call the product.

    Returns
    -------
    Recombination
        The product, the two pieces that made it, and the junctions between them.

    Raises
    ------
    KeyError
        If `reaction` is neither ``"BP"`` nor ``"LR"``.
    ValueError
        If either record does not carry one pair of the att sites the reaction takes pointing
        at each other, if the acceptor is not circular, or if either record's two sites leave
        no bases between them.
    """
    moving, taking = SUBSTRATES[reaction]
    if acceptor.topology != "circular":
        raise ValueError(
            f"a {reaction} reaction closes its product round the {taking} vector's backbone, so "
            f"{acceptor.name or 'it'} must be circular"
        )
    carrier, (c1, c2) = _oriented(carrier, moving)
    acceptor, (a1, a2) = _oriented(acceptor, taking)
    moved = Piece(carrier.name, carrier, _between(carrier, c1, c2))
    backbone = Piece(acceptor.name, acceptor, _outside(acceptor, a1, a2))
    bases = backbone.bases + moved.bases
    junctions = (
        Junction(
            writes(a1.name, c1.name),
            backbone.length - CROSSOVER,
            joined(a1.bases, c1.bases),
            Strand.FORWARD,
            acceptor.name,
            carrier.name,
        ),
        Junction(
            writes(a2.name, c2.name),
            backbone.length + moved.length - (REGION_BP - CROSSOVER),
            joined(a2.bases, c2.bases),
            Strand.REVERSE,
            acceptor.name,
            carrier.name,
        ),
    )
    features: list[Feature] = []
    primers: list[Primer] = []
    for piece, at in ((backbone, 0), (moved, backbone.length)):
        over, kept = carried(piece.record, *piece.span, offset=at - piece.span[0])
        features.extend(over)
        primers.extend(kept)
    features.extend(_junction_feature(one) for one in junctions)
    product = SequenceRecord(
        bases, topology="circular", name=name, features=tuple(features), primers=tuple(primers)
    )
    origin = _origin(acceptor, *backbone.span)
    first, second = (_turned(one, origin, len(bases)) for one in junctions)
    return Recombination(
        reaction,
        ordered(rotate(product, origin) if origin else product),
        moved,
        backbone,
        (first, second),
    )


def _oriented(record: SequenceRecord, kind: str) -> tuple[SequenceRecord, tuple[AttSite, AttSite]]:
    """Return the record read so that its site 1 lies on the top strand, and that pair.

    A record whose site 1 reads the other way is the same molecule written on the other strand,
    so it is turned over rather than refused.
    """
    pair = att_pair(record, kind)
    if pair[0].strand is Strand.REVERSE:
        record = flipped(record)
        pair = att_pair(record, kind)
    return record, pair


def _between(record: SequenceRecord, first: AttSite, second: AttSite) -> tuple[int, int]:
    """Return the span that moves: from inside site 1's overlap to inside site 2's.

    Raises
    ------
    ValueError
        If the two sites leave nothing between them.
    """
    start = (first.start + CROSSOVER) % len(record)
    length = _length(record, first.start + CROSSOVER, second.end - CROSSOVER, "moves")
    return start, start + length


def _outside(record: SequenceRecord, first: AttSite, second: AttSite) -> tuple[int, int]:
    """Return the span that stays: from inside site 2's overlap round to inside site 1's.

    Raises
    ------
    ValueError
        If the two sites leave no backbone outside them.
    """
    start = (second.end - CROSSOVER) % len(record)
    length = _length(record, second.end - CROSSOVER, first.start + CROSSOVER + len(record), "stays")
    return start, start + length


def _length(record: SequenceRecord, start: int, end: int, what: str) -> int:
    """Return how many bases lie from `start` to `end`, going round a circular record once."""
    length = end - start
    if record.topology == "circular":
        length %= len(record)
    if length <= 0:
        raise ValueError(
            f"{record.name or 'this record'} carries its two att sites with no segment between "
            f"them, so nothing {what}"
        )
    return length


def _origin(record: SequenceRecord, start: int, end: int) -> int:
    """Where `record`'s own first base falls in the piece taken from ``[start, end)``, or 0.

    Turning the product to put it back at zero is what keeps the vector's coordinates readable,
    and it moves the junction that would otherwise straddle the origin.
    """
    turns = (0, len(record)) if record.topology == "circular" else (0,)
    for turn in turns:
        if start <= turn < end:
            return turn - start
    return 0


def _turned(junction: Junction, origin: int, length: int) -> Junction:
    """Return the junction as the turned product carries it."""
    return Junction(
        junction.name,
        (junction.start - origin) % length,
        junction.bases,
        junction.strand,
        junction.before,
        junction.after,
    )


def _junction_feature(junction: Junction) -> Feature:
    """Draw one att junction, saying what it spells and which record gave which side."""
    return Feature(
        junction.name,
        "misc_recomb",
        (junction.span,),
        strand=junction.strand,
        color=JUNCTION_COLOR,
        qualifiers={
            "note": (
                f"{junction.bases}; the first {CROSSOVER} bases from {junction.before} and the "
                f"rest from {junction.after}",
            )
        },
    )

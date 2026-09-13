"""What confirms a clone: a colony PCR reading its junctions, and Sanger reads across them.

The distances are the ones ``docs/research/primer-design-and-pcr.md`` sets out for colony PCR
bands and for sequencing primers.
"""

import dataclasses
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise

from liulab_mbio import edits
from liulab_mbio.bench.gels import agarose_percent, choose_ladder
from liulab_mbio.primers import (
    ONETAQ,
    Q5,
    THRESHOLDS_FOR,
    Polymerase,
    PrimerReport,
    Thresholds,
    amplicon_sizes,
    design_pair,
    design_primer,
    evaluate_primer,
)
from liulab_mbio.protocol import Gel, Ladder, Lane
from liulab_mbio.sequence import Primer, SequenceRecord, Strand, reverse_complement

#: Vector kept either side of the junctions by a designed colony PCR pair, bases. Twice this is
#: the empty-vector band, and the note asks for every band to stay at 100 bp or more.
COLONY_FLANK = 60

#: How far into the insert a junction primer anneals, bases, for the same reason.
JUNCTION_OFFSET = 100

#: Genewiz asks for a sequencing primer 100 bases from what it reads, 50 to 60 at the closest.
SANGER_FLANK = 100

#: What the candidate plasmids are called, correct first. A reversed lane is numbered after
#: `REVERSED_CLONE` wherever an assembly holds more than one insert to turn round.
CORRECT_CLONE = "Correct clone"
EMPTY_CLONE = "Empty vector"
REVERSED_CLONE = "Reversed insert"

#: What a primer reading out of one insert is called, numbered the same way.
_JUNCTION_PRIMER = "Junction reverse"


@dataclass(frozen=True, slots=True)
class Clone:
    """One plasmid a colony may carry, and the bands a colony PCR gives from it."""

    name: str
    bands_bp: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class ColonyCheck:
    """A colony PCR reading across an assembly's junctions.

    Parameters
    ----------
    primers
        The flanking pair first, then one junction primer per insert where they were asked for.
    reports
        What each primer scored on the assembled plasmid.
    clones
        The candidates a colony can hold: the correct one, the empty vector, and one carrying
        each insert the other way round.
    annealing_temperature
        By the polymerase's rule over the two lowest Tms of the set, °C.
    extension_seconds
        For the largest band any candidate gives.
    ladder, agarose_percent
        Chosen for that band range.
    """

    primers: tuple[Primer, ...]
    reports: tuple[PrimerReport, ...]
    clones: tuple[Clone, ...]
    annealing_temperature: float
    extension_seconds: int
    ladder: Ladder
    agarose_percent: float

    @property
    def gel(self) -> Gel:
        """The gel these clones should give, one lane each."""
        return Gel(
            self.ladder,
            tuple(Lane(clone.name, clone.bands_bp) for clone in self.clones),
            title="Colony PCR",
        )

    @property
    def tells_orientation(self) -> bool:
        """Whether the gel separates every reversed insert from the correct clone.

        Two vector primers flanking the inserts never do: they amplify them whichever way round
        they sit. A junction primer does, unless the two vector primers happen to lie the same
        distance from their own junctions.
        """
        correct = next(
            (clone.bands_bp for clone in self.clones if clone.name == CORRECT_CLONE), None
        )
        turned = [clone.bands_bp for clone in self.clones if clone.name.startswith(REVERSED_CLONE)]
        return bool(turned) and all(lane != correct for lane in turned)


def colony_pcr_check(
    product: SequenceRecord,
    junctions: Sequence[int],
    *,
    vector: SequenceRecord,
    primers: tuple[Primer, ...] | None = None,
    insert_primer: bool = False,
    flank: int = COLONY_FLANK,
    junction_offset: int = JUNCTION_OFFSET,
    polymerase: Polymerase = ONETAQ,
    thresholds: Thresholds = THRESHOLDS_FOR["colony PCR"],
) -> ColonyCheck:
    """Return what a colony PCR across these junctions should show.

    An assembly of n inserts has n + 1 junctions and the inserts are the spans between them, so
    a product whose inserts cross the origin is rotated first. Without `primers`, a pair is
    designed in the vector `flank` bases outside the first and the last junction;
    `insert_primer` adds one primer per insert, annealing `junction_offset` bases into it. Those
    are what tell a reversed insert apart and what put a band of their own on each junction.
    Each candidate plasmid is amplified on its own, so the bands are simulated rather than
    derived.

    Raises
    ------
    ValueError
        If the junctions are not two or more separate positions inside the product, if fewer
        than two primers are given, if an insert is too short for a junction primer, or if the
        primers amplify nothing at all.
    """
    places = _junction_span(junctions, len(product))
    start, end = places[0], places[-1]
    inserts = tuple(pairwise(places))
    chosen = (
        list(primers)
        if primers is not None
        else list(_flanking_pair(product, start, end, flank, polymerase, thresholds))
    )
    if insert_primer:
        chosen.extend(
            _junction_primer(product, first, last, junction_offset, polymerase, thresholds, name)
            for name, (first, last) in zip(
                _numbered(_JUNCTION_PRIMER, inserts), inserts, strict=True
            )
        )
    if len(chosen) < 2:
        raise ValueError("a colony PCR needs at least two primers")
    placed = tuple(chosen)
    candidates = [(CORRECT_CLONE, product), (EMPTY_CLONE, vector)]
    candidates += [
        (name, _reversed_insert(product, first, last))
        for name, (first, last) in zip(_numbered(REVERSED_CLONE, inserts), inserts, strict=True)
    ]
    clones = tuple(Clone(name, _bands(placed, record, thresholds)) for name, record in candidates)
    sizes = tuple(sorted({bp for clone in clones for bp in clone.bands_bp}))
    if not sizes:
        raise ValueError("these primers amplify nothing on any of the candidate plasmids")
    reports = tuple(
        evaluate_primer(primer, product, polymerase=polymerase, thresholds=thresholds)
        for primer in placed
    )
    tms = sorted(report["tm"].value for report in reports)
    return ColonyCheck(
        placed,
        reports,
        clones,
        polymerase.annealing_temperature(tms[0], tms[1]),
        polymerase.extension_seconds(max(sizes)),
        choose_ladder(sizes),
        agarose_percent(sizes),
    )


@dataclass(frozen=True, slots=True)
class SangerRead:
    """A sequencing primer and the read it has to give.

    Parameters
    ----------
    primer
        Reading towards the junction it sits outside.
    distance_bp
        From its 3' end to that junction.
    read_bp
        From its 3' end to the far junction, which is what the read must cover for the insert to
        be confirmed at both ends.
    """

    primer: Primer
    distance_bp: int
    read_bp: int


def sanger_primers(
    product: SequenceRecord,
    junctions: Sequence[int],
    *,
    flank: int = SANGER_FLANK,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS_FOR["sequencing"],
) -> tuple[SangerRead, SangerRead]:
    """Return a sequencing primer reading into the inserts from outside the first and last junction.

    Every 3' end lands at least `flank` bases from its own junction, near enough for the read to
    be clean there, and `SangerRead.read_bp` is how far it must carry to reach the far junction.
    A provider whose reads are shorter needs a primer inside the inserts as well.

    Raises
    ------
    ValueError
        If the junctions are not two or more separate positions inside the product, or no primer
        fits outside the first or the last.
    """
    places = _junction_span(junctions, len(product))
    start, end = places[0], places[-1]
    longest = _longest_annealing(thresholds)
    forward = design_primer(
        product,
        start - flank - longest,
        Strand.FORWARD,
        name="Sequencing forward",
        polymerase=polymerase,
        thresholds=thresholds,
    )
    reverse = design_primer(
        product,
        end + flank + longest,
        Strand.REVERSE,
        name="Sequencing reverse",
        polymerase=polymerase,
        thresholds=thresholds,
    )
    length = len(product)
    ahead = forward.binding_sites[0].end
    behind = reverse.binding_sites[0].start
    return (
        SangerRead(forward, (start - ahead) % length, (end - ahead) % length),
        SangerRead(reverse, (behind - end) % length, (behind - start) % length),
    )


def _junction_span(junctions: Sequence[int], length: int) -> tuple[int, ...]:
    """Return the junctions in rising order, refusing a set no assembly could leave."""
    if len(junctions) < 2:
        raise ValueError(
            f"an assembly has a junction at each end of every insert, got {len(junctions)}"
        )
    places = tuple(sorted(junctions))
    if len(set(places)) != len(places):
        raise ValueError(f"two junctions share a position: {places}")
    if places[0] < 0 or places[-1] > length:
        raise ValueError(f"junctions {places} do not lie inside {length} bases")
    return places


def _numbered(label: str, inserts: Sequence[tuple[int, int]]) -> tuple[str, ...]:
    """Label one thing per insert, numbered only where there is more than one to tell apart."""
    if len(inserts) == 1:
        return (label,)
    return tuple(f"{label} {number}" for number in range(1, len(inserts) + 1))


def _flanking_pair(
    product: SequenceRecord,
    start: int,
    end: int,
    flank: int,
    polymerase: Polymerase,
    thresholds: Thresholds,
) -> tuple[Primer, Primer]:
    return design_pair(
        product,
        start - flank,
        end + flank,
        forward_name="Colony PCR forward",
        reverse_name="Colony PCR reverse",
        polymerase=polymerase,
        thresholds=thresholds,
    )


def _junction_primer(
    product: SequenceRecord,
    start: int,
    end: int,
    offset: int,
    polymerase: Polymerase,
    thresholds: Thresholds,
    name: str = _JUNCTION_PRIMER,
) -> Primer:
    if end - start <= offset:
        raise ValueError(f"an insert is shorter than the {offset} bases a junction primer needs")
    return design_primer(
        product,
        start + offset,
        Strand.REVERSE,
        name=name,
        polymerase=polymerase,
        thresholds=thresholds,
    )


def _reversed_insert(product: SequenceRecord, start: int, end: int) -> SequenceRecord:
    flipped, _ = edits.replace(product, start, end, reverse_complement(product.sequence[start:end]))
    return flipped


def _bands(
    primers: tuple[Primer, ...], record: SequenceRecord, thresholds: Thresholds
) -> tuple[int, ...]:
    """Return every band this primer set gives on one candidate plasmid.

    Binding sites are cleared first: they were found on the product, and each candidate has to
    be searched on its own.
    """
    free = [dataclasses.replace(primer, binding_sites=()) for primer in primers]
    sizes: set[int] = set()
    for index, one in enumerate(free):
        for other in free[index + 1 :]:
            sizes.update(amplicon_sizes(one, other, record, thresholds=thresholds))
    return tuple(sorted(sizes))


def _longest_annealing(thresholds: Thresholds) -> int:
    """Return the longest annealing region `design_primer` will choose."""
    band = thresholds.length
    return int(band.high if band.warn_high == float("inf") else band.warn_high)

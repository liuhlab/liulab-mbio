"""The DMX stage of DAD-GGA-DMX: how a well is marked, how deeply it is read, and what passes.

This is the **validating** experiment, not the pooled one. A design sits one per well, the well
is marked, sequenced and called on its own, and identity stays with well position throughout.
`liulab_synbio.igga` is iGGA, whose product is a pool nothing re-identifies per member, so
nothing here is a gate on a library.

Two routes mark a well and one judgement reads them. Route A ligates four DMX barcodes into the
construct in lysate; Route B amplifies each well with one barcoded primer pair. The picking,
the pass rule and the reformat are shared; the marking step, the plate and the depth floor are
the route's own, because each floor was measured on its own library prep.

**A well's marks are arithmetic, not a recorded draw.** The address factorises across the
barcode axes and one axis is the plate, so two plates on one flow cell are told apart by
construction and a demultiplexer can check an address instead of trusting a file.

The 96 barcode sequences are not shipped. They are read from a copy the user holds, the way
`liulab_mbio.bench.prices` and `liulab_mbio.ligase` read theirs; `read_kit` finds one.
"""

import csv
import itertools
import os
from collections import Counter
from collections.abc import Sequence
from dataclasses import KW_ONLY, dataclass
from pathlib import Path

from liulab_mbio.bench import plates
from liulab_mbio.checks import Check, Status, worst
from liulab_mbio.protocol.model import (
    Citation,
    Component,
    Hole,
    Incubation,
    Material,
    Plate,
    ReactionTable,
    Reference,
    Source,
    Stage,
    Step,
    ThermocyclerProgram,
    Transfer,
    Troubleshooting,
    Vessel,
    Well,
)

#: The kit: 96 plasmids in four groups of 24, used as supplied.
#: ``docs/research/synthesis-and-assembly-barcode-kit.md``.
GROUPS = 4
GROUP_SIZE = 24

#: The overhangs the four groups chain on, read on the strand the cargo reads on, so the four
#: barcodes assemble in one order and none can be left out.
CHAIN: tuple[str, ...] = ("AGGA", "GTTC", "CCTT", "TCAG", "TTCC")

#: The two universal primers flanking the design, which exist only after barcoding because they
#: come from the barcodes. Group 4's 3' constant and the reverse complement of group 1's 5'.
DMX0 = "TTTGATACGCAGGAAGATGGCCCAC"
DMX7 = "ATCGGTGACGGCGATTCTCACATTT"

#: Colonies picked per design, and the fragment counts Lund measured a clean colony at. Four is
#: the only measured anchor and is a project input, not a constant of the method.
COLONIES_PER_DESIGN = 4
CLEAN_COLONY_CURVE: tuple[tuple[int, float], ...] = (
    (2, 1.0),
    (3, 0.938),
    (5, 0.846),
    (8, 0.667),
    (12, 0.40),
    (16, 0.0),
)

#: The plate each stage uses. Colonies are picked into 384-well plates and four of those are
#: compressed into one 1536-well plate, which is why the format parameter has to reach 1536.
#: Qian SI Day 3 and Day 4.1. Route B amplifies in the half-skirted 96-well PCR plate its
#: barcodes address, LevSeq SI step 2.
PICKED_WELLS = 384
COMPRESSED_WELLS = 1536
INDEX_WELLS = 96
PICKED_CATALOG = "Beckman Coulter #c74290"
COMPRESSED_CATALOG = "Greiner #782270"
INDEX_CATALOG = "USA Scientific #1402-9700"

#: How many picked plates one compressed plate holds.
PLATES_COMPRESSED = COMPRESSED_WELLS // PICKED_WELLS

#: Medium in a picked well, µL, and what grows in it. Qian SI Day 3.
CULTURE_UL = 60.0

#: One well of the barcoding reaction, µL: lysate, barcodes, water and mastermix to 3.5 µL.
#: Qian SI Day 4.1.
LYSATE_UL = 2.0
BARCODE_UL = 0.5
WATER_UL = 0.5
MASTERMIX_UL = 0.5
WELL_UL = LYSATE_UL + BARCODE_UL + WATER_UL + MASTERMIX_UL

#: What moves the liquid, and what the colonies are picked with, both Qian's. The handler is
#: cited on the move it makes; the picker is equipment, which carries no citation, so the page
#: reaches it through the Qian reference alone.
ACOUSTIC = "ECHO 525 acoustic liquid handler"
PICKER = "QPix XE Microbial Colony Picker"

#: What a quarter of a picked plate is sampled into a Route B index plate with, and how much of
#: each well goes. LevSeq SI step 4. An acoustic handler does the same move where one is booked.
SAMPLE_UL = 1.0
MULTICHANNEL = "multichannel pipette"

#: Route B's index PCR: a colony PCR straight from overnight culture, cycled with a touchdown.
#: The master mix is one well's share of the mix LevSeq states per plate, and the cycling is its
#: thermal-cycler table with the two loop lines spelled out: ten touchdown cycles and 25 more,
#: 35 in all. ``docs/research/route-b-index-pcr.md`` section 2.
INDEX_PCR_UL = 10.0
INDEX_MIX_UL = 7.0
INDEX_PRIMER_UL = 2.0
INDEX_DENATURE_C = 95.0
INDEX_TOUCHDOWN_C = (68.0, 63.5)
INDEX_TOUCHDOWN_STEP_C = 0.5
INDEX_PLATEAU_CYCLES = 25

#: The cycles the touchdown takes, which is every step from one end of it to the other.
INDEX_TOUCHDOWN_CYCLES = (
    round((INDEX_TOUCHDOWN_C[0] - INDEX_TOUCHDOWN_C[1]) / INDEX_TOUCHDOWN_STEP_C) + 1
)

#: The mix is made 1.5 times what a full plate takes, which is what LevSeq's per-plate table
#: spells out. Nothing published says what the surplus is for, or how to scale a part plate.
INDEX_OVERAGE = 0.5

#: Extension at the top of the touchdown, seconds: what the authors ran for a gene below 1 kb.
#: The rule is minimally one minute a kilobase, so a longer design wants it raised. SI note 5b.
INDEX_EXTENSION_SECONDS = 120

#: What the index PCR is pipetted from, LevSeq SI step 1. A catalogue number is ordered as
#: written, so each is the SI's own.
THERMOPOL_CATALOG = "NEB #B9004"
DNTP_CATALOG = "NEB #N0447"
TAQ_CATALOG = "NEB #M0267"
DMSO_CATALOG = "MP Biomedicals #194819"

#: Colonies a 25 cm BioAssay plate carries before picking gets hard. Qian SI Day 2.
BIOASSAY_COLONIES = 2500
BIOASSAY_CATALOG = "Corning #431111"

#: Columns a pooled 1536-well plate takes, because one saturates. Qian SI Day 4.1.
POOL_COLUMNS = 2

#: Where the 96 barcode sequences are read from when a caller names no file. The package ships
#: none; this points at a copy the user holds.
KIT_ENV = "LIULAB_SYNBIO_DMX_BARCODES"

#: The columns one of those files holds, in any order.
KIT_COLUMNS = ("name", "group", "index", "overhang5", "umi", "overhang3", "final_seq")

#: What a file has to hold to be one of these, said in the refusal so a caller need not guess.
KIT_EXPECTED = (
    f"a tab-separated file with the columns {', '.join(KIT_COLUMNS)}, holding "
    f"{GROUPS * GROUP_SIZE} rows: {GROUPS} groups of {GROUP_SIZE}, indexed 1 to {GROUP_SIZE} "
    "within each group, every UMI distinct"
)

#: The documents these numbers were read from.
SOURCES: dict[str, Source] = {
    "Qian SI": Source(
        "Qian et al. 2026, Supplementary Information",
        edition="Nat. Commun. 10.1038/s41467-026-76740-5",
        read_as="held under reference_docs/",
        date="2026-10-06",
    ),
    "LevSeq": Source(
        "Long, Y. et al. 2025, LevSeq, Supporting Information",
        read_as="held under reference_docs/",
        date="2026-10-06",
    ),
}


@dataclass(frozen=True, slots=True)
class Barcode:
    """One member of the kit: which group it belongs to, what it spells, and how it chains.

    Parameters
    ----------
    name
        What the kit calls it, such as ``"DMX_1_1"``.
    group, index
        Which of the `GROUPS` groups it is in, counting from one, and which of `GROUP_SIZE`.
    overhang5, overhang3
        The two overhangs it chains on, which are its group's and not its own.
    umi
        The bases that tell it from the rest of its group.
    sequence
        The whole barcode as it is ordered.
    """

    name: str
    _: KW_ONLY
    group: int
    index: int
    overhang5: str
    overhang3: str
    umi: str
    sequence: str


@dataclass(frozen=True, slots=True)
class Kit:
    """The 96 barcodes a user holds, grouped as the kit groups them.

    Parameters
    ----------
    path
        The file they were read from, so a report can say which copy was used.
    barcodes
        Every barcode, in group and then index order.
    """

    path: Path
    barcodes: tuple[Barcode, ...]

    def group(self, group: int) -> tuple[Barcode, ...]:
        """Return one group's barcodes, in index order.

        Raises
        ------
        ValueError
            If the kit has no such group.
        """
        if not 1 <= group <= GROUPS:
            raise ValueError(f"the kit has groups 1 to {GROUPS}, and {group} is not one of them")
        return tuple(one for one in self.barcodes if one.group == group)

    def at(self, group: int, index: int) -> Barcode:
        """Return the barcode at one group and index, both counting from one.

        Raises
        ------
        ValueError
            If the kit holds no barcode there.
        """
        found = [one for one in self.group(group) if one.index == index]
        if not found:
            raise ValueError(f"the kit holds no barcode {index} of group {group}")
        return found[0]


@dataclass(frozen=True, slots=True)
class Route:
    """One way of marking a well, and how deep its read has to be before anyone calls it.

    The two routes differ in their marking step, their plate and their floor; the picking, the
    pass rule and the reformat are shared. Which one a project runs is the project's choice.

    Parameters
    ----------
    name
        ``"A"`` or ``"B"``.
    marking
        What the step does, in a few words.
    well_axes
        How many marks each axis addressing the well carries, in chain order. The well's address
        factorises across them.
    plate_axis
        How many marks the one axis carrying the plate carries.
    wanted_reads
        The mark a well's depth has to meet, by `meets`, before it is called. A project may
        raise it.
    tolerable_reads
        The mark that same comparison holds a well to for any verdict at all; under it the well
        carries none. Equal to `wanted_reads` where the route publishes one number rather than
        two.
    inclusive
        Whether a depth equal to either of a route's marks meets it. Qian publishes a consensus
        depth *above* 150, and LevSeq's SI checklist an alignment count *above* 20, so neither
        route calls a well that lands exactly on a mark. LevSeq's article reads the same 20 as a
        minimum instead; the SI governs, because its checklist decides whether a well's data may
        be used, which is what this decides, while the article's number is where its software
        warns. Where two sources still contest a boundary the stricter reading stands: a well
        wrongly failed is re-sequenced, a well wrongly passed contaminates a result.
    source
        The key of the `SOURCES` entry the two depths were read from.
    """

    name: str
    _: KW_ONLY
    marking: str
    well_axes: tuple[int, ...]
    plate_axis: int
    wanted_reads: int
    tolerable_reads: int
    source: str
    inclusive: bool = True

    def meets(self, reads: int, mark: int) -> bool:
        """Whether this many reads meet `mark`, by the comparison this route publishes."""
        return reads >= mark if self.inclusive else reads > mark

    @property
    def wells_per_plate(self) -> int:
        """How many wells the axes addressing a well can tell apart."""
        total = 1
        for axis in self.well_axes:
            total *= axis
        return total

    @property
    def capacity(self) -> int:
        """How many wells the route addresses in all, over every plate it tells apart."""
        return self.wells_per_plate * self.plate_axis


#: Four DMX barcodes a well, three groups addressing the well and the fourth the plate, one
#: barcode across a whole plate from a reservoir. 24 x 24 x 24 is 13,824 against the 1,536 a
#: compressed plate needs. Depth from Qian SI, "Generation of consensus sequences".
ROUTE_A = Route(
    "A",
    marking="ligate one barcode from each of the four kit groups into the construct, in lysate",
    well_axes=(GROUP_SIZE, GROUP_SIZE, GROUP_SIZE),
    plate_axis=GROUP_SIZE,
    wanted_reads=150,
    tolerable_reads=150,
    source="Qian SI",
    inclusive=False,
)

#: One barcoded primer pair a well, 96 forward marks addressing the well and 96 reverse the
#: plate: 9,216 wells on the 192 primers already held. Twenty reads wanted, ten tolerable, both
#: depths to exceed, so a well at exactly twenty is called and warns, where the SI's suboptimal
#: branch puts it. LevSeq pairs the floor with a second criterion, a mean error below 10%, which
#: nothing here judges — a named hole: it is a mean over per-position base counts, and a well
#: reaches this package as a read count and a consensus call.
ROUTE_B = Route(
    "B",
    marking="amplify each well with one barcoded primer pair",
    well_axes=(INDEX_WELLS,),
    plate_axis=INDEX_WELLS,
    wanted_reads=20,
    tolerable_reads=10,
    source="LevSeq",
    inclusive=False,
)

#: Both routes, by the name a project names one with.
ROUTES: dict[str, Route] = {ROUTE_A.name: ROUTE_A, ROUTE_B.name: ROUTE_B}


@dataclass(frozen=True, slots=True)
class Address:
    """Which marks name one well, worked out from where the well is rather than looked up.

    Parameters
    ----------
    route
        The route whose axes this was worked out on.
    plate, well
        Which plate and which well of it, both counting from zero.
    well_marks
        One mark an axis addressing the well, each counting from one, in chain order.
    plate_mark
        The mark the whole plate carries, counting from one.
    """

    route: Route
    _: KW_ONLY
    plate: int
    well: int
    well_marks: tuple[int, ...]
    plate_mark: int

    @property
    def marks(self) -> tuple[int, ...]:
        """Every mark this well carries, the well's axes first and the plate's last."""
        return (*self.well_marks, self.plate_mark)


def address(route: Route, *, plate: int, well: int) -> Address:
    """Return the marks naming one well, factorising its address across the route's axes.

    The combination is arithmetic and no file records it, which is what keeps two plates apart
    on one flow cell and lets a demultiplexer check an address rather than trust one.

    Raises
    ------
    ValueError
        If either index is negative, or runs past what the route addresses.

    Examples
    --------
    >>> address(ROUTE_A, plate=0, well=25).marks
    (2, 2, 1, 1)
    >>> address(ROUTE_B, plate=2, well=0).marks
    (1, 3)
    """
    if plate < 0 or well < 0:
        raise ValueError(f"a well's address counts from zero, got plate {plate} and well {well}")
    if well >= route.wells_per_plate:
        raise ValueError(
            f"route {route.name} tells {route.wells_per_plate} wells of a plate apart, and well "
            f"{well} is past the last of them"
        )
    if plate >= route.plate_axis:
        raise ValueError(
            f"route {route.name} tells {route.plate_axis} plates apart, and plate {plate} is "
            "past the last of them"
        )
    marks: list[int] = []
    left = well
    for axis in route.well_axes:
        marks.append(left % axis + 1)
        left //= axis
    return Address(route, plate=plate, well=well, well_marks=tuple(marks), plate_mark=plate + 1)


def barcodes_for(kit: Kit, one: Address) -> tuple[Barcode, ...]:
    """Return the four kit barcodes an address names, in chain order.

    Raises
    ------
    ValueError
        If the address was not worked out on Route A, which is the route the kit marks.
    """
    if one.route is not ROUTE_A:
        raise ValueError(
            f"the DMX barcode kit marks a well on route {ROUTE_A.name}, and this address was "
            f"worked out on route {one.route.name}"
        )
    return tuple(kit.at(group, mark) for group, mark in enumerate(one.marks, start=1))


def read_kit(path: str | os.PathLike[str] | None = None) -> Kit:
    """Read the 96 barcode sequences from a copy the user holds.

    The package ships none. `path` names the file, or `KIT_ENV` does.

    Raises
    ------
    ValueError
        If no file is named, or the file is not one of these, saying what one is.

    Examples
    --------
    >>> read_kit("dmx-barcodes.tsv").barcodes[0].group  # doctest: +SKIP
    1
    """
    named = path if path is not None else os.environ.get(KIT_ENV)
    if not named:
        raise ValueError(
            f"the DMX barcode sequences are not shipped: name the file, or set {KIT_ENV}. "
            f"It is {KIT_EXPECTED}"
        )
    one = Path(named)
    try:
        barcodes = _kit_rows(one.read_text(encoding="utf-8-sig").splitlines())
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            f"{one.name} is not the DMX barcode kit ({error}); expected {KIT_EXPECTED}"
        ) from error
    return Kit(one, barcodes)


def depth_check(route: Route, reads: int, *, wanted: int | None = None) -> Check:
    """Return the verdict on whether a well was read deeply enough to call.

    A well below the route's tolerable depth carries **no verdict**, not a failure: it is read
    again or picked again, and reformatting does not compact it out. Between tolerable and
    wanted it is called and warns. `wanted` raises the route's own mark, which a project may do
    and may not lower.

    Raises
    ------
    ValueError
        If `reads` is negative, or `wanted` is below the route's own mark.

    Examples
    --------
    >>> depth_check(ROUTE_B, 12).status, depth_check(ROUTE_B, 4).status
    ('warn', None)
    """
    if reads < 0:
        raise ValueError(f"a well cannot carry {reads} reads")
    mark = route.wanted_reads if wanted is None else wanted
    if mark < route.wanted_reads:
        raise ValueError(
            f"route {route.name} wants {route.wanted_reads} reads a well and a project may only "
            f"raise that, not lower it to {mark}"
        )
    status: Status | None = None
    if route.meets(reads, mark):
        status = "pass"
    elif route.meets(reads, route.tolerable_reads):
        status = "warn"
    return Check(
        "reads_per_well",
        status,
        float(reads),
        f"route {route.name} wants {mark} reads a well and tolerates "
        f"{route.tolerable_reads}; below that no read is deep enough to call",
    )


def identity_check(called: Sequence[str], designed: str) -> Check:
    """Return the verdict on whether a well's call is the design, base for base.

    A pass is an exact match across the whole designed region — both entry overhangs, the
    fragment, the stuffer and the barcode. A silent mismatch is not a pass: it leaves the
    protein right and the barcode wrong, and a barcode that no longer names its member cannot be
    put right by the linkage read. More than one consensus is mixed, and mixed fails. No
    consensus at all is nothing to judge, so it carries no verdict.

    Raises
    ------
    ValueError
        If `designed` is empty.

    Examples
    --------
    >>> identity_check(("ACGT",), "ACGT").status
    'pass'
    >>> identity_check(("ACGT", "ACGA"), "ACGT").status
    'fail'
    """
    if not designed:
        raise ValueError("a well is judged against a designed region, and this one is empty")
    consensus = [one.upper() for one in called]
    if not consensus:
        return Check("well_identity", None, 0.0, "no consensus was called from this well")
    if len(consensus) > 1:
        return Check(
            "well_identity",
            "fail",
            float(len(consensus)),
            f"{len(consensus)} consensus sequences: the well is mixed",
        )
    if consensus[0] != designed.upper():
        return Check("well_identity", "fail", 1.0, "the call is not the design, base for base")
    return Check("well_identity", "pass", 1.0, "the call is the design across its whole length")


@dataclass(frozen=True, slots=True)
class WellVerdict:
    """What one well's read came to, and whether reformatting carries it forward.

    Parameters
    ----------
    well
        Where the well is, as a plate name and a well name.
    checks
        The depth and the identity, in that order.
    """

    well: Well
    checks: tuple[Check, ...]

    @property
    def status(self) -> Status:
        """The worst verdict the well carries, by `liulab_mbio.checks.worst`."""
        return worst(one.status for one in self.checks)

    @property
    def called(self) -> bool:
        """Whether anything judged this well at all."""
        return any(one.status is not None for one in self.checks)


def judge_well(
    well: Well,
    *,
    route: Route,
    reads: int,
    called: Sequence[str],
    designed: str,
    wanted: int | None = None,
) -> WellVerdict:
    """Judge one well: deep enough to call, and then an exact match or not.

    The two questions run in that order, which is why a shallow well is never a failure.

    Raises
    ------
    ValueError
        For any reason `depth_check` or `identity_check` refuses.
    """
    depth = depth_check(route, reads, wanted=wanted)
    if depth.status is None:
        return WellVerdict(well, (depth,))
    return WellVerdict(well, (depth, identity_check(called, designed)))


def reformat(verdicts: Sequence[WellVerdict]) -> tuple[WellVerdict, ...]:
    """Return the wells reformatting carries forward: everything that did not fail.

    A well nobody could call is kept. Compacting it out would throw away a design that may be
    clean and has only been read too thinly.

    Examples
    --------
    >>> reformat(())
    ()
    """
    return tuple(one for one in verdicts if one.status != "fail")


def clean_colony_chance(fragments: int) -> float:
    """Return the chance one picked colony carries a clean copy of a design of this many pieces.

    Linear interpolation between Lund's six measured anchors, `CLEAN_COLONY_CURVE`: a design in
    two pieces came out clean every time, one in sixteen never. Fewer pieces than the first
    anchor takes the first anchor's value and more than the last takes the last's, because
    nothing was measured outside them.

    Raises
    ------
    ValueError
        If `fragments` is not positive.

    Examples
    --------
    >>> clean_colony_chance(8), clean_colony_chance(16)
    (0.667, 0.0)
    """
    if fragments < 1:
        raise ValueError(f"a design is built from at least one fragment, got {fragments}")
    anchors = CLEAN_COLONY_CURVE
    if fragments <= anchors[0][0]:
        return anchors[0][1]
    for (low, left), (high, right) in itertools.pairwise(anchors):
        if fragments <= high:
            return left + (right - left) * (fragments - low) / (high - low)
    return anchors[-1][1]


@dataclass(frozen=True, slots=True)
class Design:
    """One design a well may be read back for, and how many pieces it was assembled from.

    Parameters
    ----------
    name
        What the design is called.
    fragments
        How many pieces were joined to build it, which is what its chance of a clean colony
        falls with.
    """

    name: str
    fragments: int

    def __post_init__(self) -> None:
        """Refuse a design built from no piece at all."""
        if self.fragments < 1:
            raise ValueError(
                f"{self.name} is built from at least one fragment, got {self.fragments}"
            )

    @property
    def chance(self) -> float:
        """The chance one picked colony carries a clean copy, by `clean_colony_chance`."""
        return clean_colony_chance(self.fragments)


def validated(designs: Sequence[Design], floor: int | None) -> tuple[Design, ...]:
    """Return the designs a fragment-count floor reads back, in the order they were given.

    No floor reads nothing, because validation is optional and a library headed for a pooled
    screen takes its identity from that screen. A floor of zero reads every design. Nothing here
    supplies a default: `clean_colony_chance` gives a design's chance of a clean colony, not the
    chance worth paying to check, and that is the project's own call.

    Raises
    ------
    ValueError
        If the floor is negative.

    Examples
    --------
    >>> some = (Design("short", 2), Design("long", 8))
    >>> [one.name for one in validated(some, 5)], validated(some, None)
    (['long'], ())
    """
    if floor is None:
        return ()
    if floor < 0:
        raise ValueError(f"a fragment-count floor counts fragments, and {floor} is below none")
    return tuple(one for one in designs if one.fragments >= floor)


def picked_plate(name: str, colonies: int) -> Plate:
    """Return one plate of picked colonies: `colonies` wells of selective medium.

    Picking fills one quarter of the plate at a time, in the order a head built for the index
    format covers it. A part-filled plate then gives full index plates rather than part-filled
    ones, no reverse mark is spent on a plate that is mostly empty, and a design's colonies stay
    together.

    Raises
    ------
    ValueError
        If more colonies are picked than the format holds.

    Examples
    --------
    >>> picked = picked_plate("picked 1", 288)
    >>> picked.wells, len(picked.seating), picked.seating["B1"]
    (384, 288, 'quarter 3')
    """
    if colonies > PICKED_WELLS:
        raise ValueError(f"{colonies} colonies do not fit a {PICKED_WELLS}-well plate")
    seating: dict[str, str] = {}
    left = colonies
    for number, quarter in enumerate(plates.interleave(PICKED_WELLS, INDEX_WELLS), start=1):
        seating |= dict.fromkeys(quarter[: max(left, 0)], f"quarter {number}")
        left -= len(quarter)
    return plates.plate(
        name,
        PICKED_WELLS,
        catalog=PICKED_CATALOG,
        holds=f"one picked colony each in {CULTURE_UL:g} µL low-salt LB with carbenicillin",
        seating=seating,
        note=f"{colonies} of {PICKED_WELLS} wells picked, a quarter at a time",
    )


def compressed_plate(name: str) -> Plate:
    """Return the plate four picked plates are compressed into, where the barcoding runs."""
    return plates.plate(
        name,
        COMPRESSED_WELLS,
        catalog=COMPRESSED_CATALOG,
        holds=f"{WELL_UL:g} µL: {LYSATE_UL:g} µL lysate, barcodes, water and mastermix",
        note=f"{PLATES_COMPRESSED} picked plates to one",
    )


def index_plate(name: str, samples: int, *, plate: int = 0) -> Plate:
    """Return one Route B plate: each well amplified with the pair its address names.

    Raises
    ------
    ValueError
        If more samples are given than the plate holds, or `plate` is past the last the route
        tells apart.

    Examples
    --------
    >>> index_plate("index 1", 3, plate=1).seating["A2"]
    'forward 2, reverse 2'
    """
    if samples > INDEX_WELLS:
        raise ValueError(f"{samples} samples do not fit a {INDEX_WELLS}-well plate")
    seating = {}
    names = plates.plate(name, INDEX_WELLS).well_names
    for well in range(samples):
        one = address(ROUTE_B, plate=plate, well=well)
        seating[names[well]] = f"forward {one.well_marks[0]}, reverse {one.plate_mark}"
    return plates.plate(
        name,
        INDEX_WELLS,
        catalog=INDEX_CATALOG,
        holds="one barcoded PCR each, the pair its own address names",
        seating=seating,
        note=f"plate mark {plate + 1}; {samples} of {INDEX_WELLS} wells used",
    )


def bioassay_plate(name: str) -> Vessel:
    """Return the 25 cm plate the colonies are picked from.

    It is a vessel and not a plate: its colonies land where they land, so they have no
    positions to seat.
    """
    return Vessel(
        name,
        kind="25 cm BioAssay plate",
        catalog=BIOASSAY_CATALOG,
        holds=f"about {BIOASSAY_COLONIES:,} colonies, which is the density picking wants",
        note="100 µg/mL carbenicillin, overnight at 37 °C",
    )


def compression(picked: Sequence[Plate], compressed: Plate) -> Transfer:
    """Return the acoustic move of every picked well into the compressed plate.

    A well that failed is simply left out of `picked`'s seating, so the destination is dense:
    compaction and a plain acoustic transfer are one shape.

    Raises
    ------
    ValueError
        If the picked wells do not fit the compressed plate.

    Examples
    --------
    >>> one = compression([picked_plate("picked 1", 2)], compressed_plate("lysate"))
    >>> one.moves[0].volume_ul, one.instrument
    (2.0, 'ECHO 525 acoustic liquid handler')
    """
    wells = [
        Well(plate.name, name)
        for plate in picked
        for name in (plate.seating or dict.fromkeys(plate.well_names))
    ]
    return plates.compact(
        wells,
        compressed,
        LYSATE_UL,
        title=f"Compress {len(picked)} picked plates into {compressed.name}",
        instrument=ACOUSTIC,
        note="Invert the picked plates for 30 minutes first, so the cells gather at the meniscus",
        citation=Citation("Qian SI", "Day 4.1"),
    )


def sampling(picked: Plate, index: Sequence[Plate]) -> tuple[Transfer, ...]:
    """Return Route B's move out of one picked plate: one quarter of it into each index plate.

    A 384-well plate's wells sit at half a 96-well plate's spacing, so one well in four lines up
    under a standard multichannel head and the plate is covered in four passes. A quarter nobody
    picked into is no index plate and no move, which is what keeps a part-filled picked plate
    giving full index plates.

    Raises
    ------
    ValueError
        If the picked plate holds no colony, or `index` is not one plate per filled quarter.

    Examples
    --------
    >>> moves = sampling(picked_plate("picked 1", 288), [index_plate("index 1", 96)] * 3)
    >>> len(moves), len(moves[0].moves), moves[0].instrument
    (3, 96, 'multichannel pipette')
    """
    filled = picked.seating or dict.fromkeys(picked.well_names)
    quarters = [
        tuple(name for name in quarter if name in filled)
        for quarter in plates.interleave(picked.wells, INDEX_WELLS)
    ]
    taken = [one for one in quarters if one]
    if not taken:
        raise ValueError(f"{picked.name} holds no picked colony to sample")
    if len(taken) != len(index):
        raise ValueError(
            f"{picked.name} was picked into {len(taken)} quarter(s) and {len(index)} index "
            "plate(s) were given; one index plate covers one quarter"
        )
    return tuple(
        plates.compact(
            [Well(picked.name, name) for name in quarter],
            plate,
            SAMPLE_UL,
            title=f"Sample a quarter of {picked.name} into {plate.name}",
            instrument=MULTICHANNEL,
            note="one well in four, which is the spacing a head built for the index format reads",
            citation=Citation("LevSeq", "step 4"),
        )
        for quarter, plate in zip(taken, index, strict=True)
    )


@dataclass(frozen=True, slots=True)
class Validation:
    """Which designs are read back, on which route, and the plates that takes.

    Picking is shared and sized once, because the picking format no longer differs by route: one
    well a colony a design, filled a quarter of a plate at a time. Each route's own plates follow
    from those wells.

    Parameters
    ----------
    route
        The route the project named.
    designs
        The designs read back, which `validated` chose from the project's floor.
    floor
        That floor, carried so the protocol can print it beside each design's chance.
    colonies
        Colonies picked per design.
    """

    route: Route
    designs: tuple[Design, ...]
    _: KW_ONLY
    floor: int
    colonies: int = COLONIES_PER_DESIGN

    def __post_init__(self) -> None:
        """Refuse a read of no design, or of no colony per design."""
        if not self.designs:
            raise ValueError("a validation reads at least one design; a floor reading none is None")
        if self.colonies < 1:
            raise ValueError(f"{self.colonies} colonies a design reads nothing back")

    @property
    def wells(self) -> int:
        """How many wells the read takes: one a colony a design."""
        return len(self.designs) * self.colonies

    @property
    def picked(self) -> tuple[Plate, ...]:
        """The picked plates, each filled a quarter at a time and the last one part-filled."""
        full, rest = divmod(self.wells, PICKED_WELLS)
        sizes = [PICKED_WELLS] * full + ([rest] if rest else [])
        return tuple(
            picked_plate(f"picked {number}", size) for number, size in enumerate(sizes, start=1)
        )

    @property
    def compressed(self) -> tuple[Plate, ...]:
        """Route A's 1536-well plates, `PLATES_COMPRESSED` picked plates compressed into one."""
        self._only(ROUTE_A)
        many = -(-len(self.picked) // PLATES_COMPRESSED)
        return tuple(compressed_plate(f"lysate {number}") for number in range(1, many + 1))

    @property
    def index(self) -> tuple[Plate, ...]:
        """Route B's 96-well plates, one a filled quarter, marked in one run across the run."""
        self._only(ROUTE_B)
        made: list[Plate] = []
        for one in self.picked:
            left = len(one.seating)
            while left > 0:
                made.append(
                    index_plate(f"index {len(made) + 1}", min(left, INDEX_WELLS), plate=len(made))
                )
                left -= INDEX_WELLS
        return tuple(made)

    def _only(self, route: Route) -> None:
        """Refuse a plate the other route never pours.

        Raises
        ------
        ValueError
            If this validation runs the other route.
        """
        if self.route is not route:
            raise ValueError(
                f"only route {route.name} pours this plate, and this read runs route "
                f"{self.route.name}"
            )


def validation(
    route: Route,
    designs: Sequence[Design],
    floor: int | None,
    *,
    colonies: int = COLONIES_PER_DESIGN,
) -> Validation | None:
    """Return what reading `designs` back on `route` takes, or `None` where the floor reads none.

    `None` is the answer for a project that states no floor and for one whose floor is above
    every design: either way nothing is read, and a protocol then carries no validation at all.

    Raises
    ------
    ValueError
        If the floor is negative, or fewer than one colony a design is picked.

    Examples
    --------
    >>> validation(ROUTE_B, (Design("one", 4),), 2).wells
    4
    >>> validation(ROUTE_B, (Design("one", 4),), None) is None
    True
    """
    if floor is None:
        return None
    read = validated(designs, floor)
    if not read:
        return None
    return Validation(route, read, floor=floor, colonies=colonies)


def pooling(compressed: Plate, reservoir: Vessel) -> Transfer:
    """Return every barcoded well run into one reservoir, which is what a pool is.

    Examples
    --------
    >>> pooled = pooling(compressed_plate("lysate"), Vessel("reservoir"))
    >>> len(pooled.moves)
    1536
    """
    return plates.pool(
        plates.wells_of(compressed),
        Well(reservoir.name, "1"),
        WELL_UL,
        title=f"Pool {compressed.name}",
        note=f"Invert and spin at 200 x g; {POOL_COLUMNS} miniprep columns, because one saturates",
        citation=Citation("Qian SI", "Day 4.1"),
    )


#: Route B's marks, which the package holds none of. The two annealing regions are published and
#: bind the DMX vector verbatim; which 192 index sequences sit on their 5' ends is a plate the lab
#: buys and holds, as Route A's kit is. ``docs/research/route-b-index-primers.md`` section 7.
INDEX_MARKS = Hole(
    "B1",
    "no index mark set is named for the barcoded primer pairs",
    "lab",
    where="route B, the pair marking one well",
    filled_by="the prepared primer plate the lab holds",
    issue="liuhlab/liulab-mbio#225",
)

#: What 0.05 µL of Taq is worth. LevSeq states the volume and never the enzyme's concentration,
#: and NEB's specification for M0267 could not be read, so the table prints the volume and the
#: unit count stands empty. ``docs/research/route-b-index-pcr.md`` section 8.
INDEX_TAQ_UNITS = Hole(
    "B2",
    "the units of Taq one index PCR takes",
    "unread",
    where="route B, the index PCR's polymerase",
    filled_by="NEB's own specification for M0267, which gives the stock concentration",
    issue="liuhlab/liulab-mbio#225",
)


def index_pcr_reaction(reactions: int = 1) -> ReactionTable:
    """Return LevSeq's index PCR, one `INDEX_PCR_UL` µL reaction a well.

    The first five lines are the master mix, each one well's share of the mix LevSeq states per
    plate; at `INDEX_OVERAGE` they scale back to that table. The pair and the culture go in a
    well at a time, because each well takes its own pair. How many units 0.05 µL of Taq is
    nobody published, so `INDEX_TAQ_UNITS` stands where the unit count would.

    Examples
    --------
    >>> round(sum(one.volume_ul for one in index_pcr_reaction().components), 2)
    10.0
    >>> index_pcr_reaction().mix_volumes(96)[:4]
    (144.0, 28.8, 7.2, 57.6)
    """
    mix = Citation("LevSeq", "step 1")
    return ReactionTable(
        (
            Component("ThermoPol Reaction Buffer", 1.0, stock="10X", final="1X", citation=mix),
            Component("dNTP mix", 0.2, stock="10 mM each", final="0.2 mM each", citation=mix),
            Component("Taq DNA Polymerase", 0.05, citation=mix),
            Component("DMSO", 0.4, stock="100%", final="4% (v/v)", citation=mix),
            Component("Nuclease-free water", 5.35, citation=mix),
            Component(
                "Barcoded primer mix",
                INDEX_PRIMER_UL,
                stock="1 µM each",
                final="0.2 µM each",
                master_mix=False,
                citation=Citation("LevSeq", "step 3"),
            ),
            Component(
                "Overnight culture",
                SAMPLE_UL,
                master_mix=False,
                citation=Citation("LevSeq", "step 4"),
            ),
        ),
        title="Index PCR, one well a sample",
        reactions=reactions,
        overage=INDEX_OVERAGE,
    )


def index_pcr_program() -> ThermocyclerProgram:
    """Return LevSeq's touchdown program, the SI's two loop lines spelled out.

    Ten cycles drop the annealing temperature `INDEX_TOUCHDOWN_STEP_C` each, over
    `INDEX_TOUCHDOWN_C`; `INDEX_PLATEAU_CYCLES` more anneal and extend together at the top of it.
    The touchdown is one stage a cycle, because a stage holds one temperature.

    Examples
    --------
    >>> stages = index_pcr_program().stages
    >>> stages[1].incubations[1].temperature_c, stages[10].incubations[1].temperature_c
    (68.0, 63.5)
    """
    cite = Citation("LevSeq", "thermal cycler table")
    top = INDEX_TOUCHDOWN_C[0]
    denature = Incubation("Denature", INDEX_DENATURE_C, 20, cite)
    extend = Incubation("Extend", top, INDEX_EXTENSION_SECONDS, cite)
    return ThermocyclerProgram(
        (
            Stage((Incubation("Initial denaturation", INDEX_DENATURE_C, 300, cite),)),
            *(
                Stage(
                    (
                        denature,
                        Incubation("Anneal", top - at * INDEX_TOUCHDOWN_STEP_C, 20, cite),
                        extend,
                    )
                )
                for at in range(INDEX_TOUCHDOWN_CYCLES)
            ),
            Stage(
                (denature, Incubation("Anneal and extend", top, INDEX_EXTENSION_SECONDS, cite)),
                cycles=INDEX_PLATEAU_CYCLES,
            ),
            Stage((Incubation("Final extension", top, 300, cite),)),
            Stage((Incubation("Hold", 4.0, None, cite),)),
        ),
        title="Index PCR",
    )


def chances(designs: Sequence[Design]) -> tuple[str, ...]:
    """Return one line a fragment count: how many designs it covers, and each one's chance.

    Grouped rather than listed, because every design of one fragment count carries the same
    chance and a library has more designs than a page has room for.

    Examples
    --------
    >>> chances((Design("a", 2), Design("b", 2), Design("c", 8)))
    ('2 fragment(s): 2 design(s), 100.0% of picks clean', '8 fragment(s): 1 design(s), 66.7% of picks clean')
    """
    counted = Counter(one.fragments for one in designs)
    return tuple(
        f"{pieces} fragment(s): {number} design(s), "
        f"{clean_colony_chance(pieces):.1%} of picks clean"
        for pieces, number in sorted(counted.items())
    )


def validation_materials(one: Validation) -> tuple[Material, ...]:
    """Return what reading these designs back consumes, beyond the designs themselves."""
    made = [
        Material(
            "25 cm BioAssay plate",
            supplier="Corning",
            catalog=BIOASSAY_CATALOG.split("#")[-1],
            amount="one spot a design",
            note=f"carbenicillin at 100 µg/mL; about {BIOASSAY_COLONIES:,} colonies a plate",
            citation=Citation("Qian SI", "Day 2"),
        ),
        Material(
            f"{PICKED_WELLS}-well culture plate",
            supplier="Beckman Coulter",
            catalog=PICKED_CATALOG.split("#")[-1],
            amount=f"{len(one.picked)}, one a {PICKED_WELLS} wells",
            note=f"{CULTURE_UL:g} µL low-salt LB with carbenicillin a well",
            citation=Citation("Qian SI", "Day 3"),
        ),
    ]
    if one.route is ROUTE_A:
        made += [
            Material(
                f"{COMPRESSED_WELLS}-well barcoding plate",
                supplier="Greiner",
                catalog=COMPRESSED_CATALOG.split("#")[-1],
                amount=f"{len(one.compressed)}, one a {PLATES_COMPRESSED} picked plates",
                citation=Citation("Qian SI", "Day 4.1"),
            ),
            Material(
                "DMX barcode kit",
                storage="-20 °C",
                amount=f"{GROUPS * GROUP_SIZE} plasmids, {GROUPS} groups of {GROUP_SIZE}",
                note="the lab's own stock; its sequences are read from the copy you hold",
            ),
            Material(
                "Barcoding master mix",
                storage="-20 °C",
                amount=f"{MASTERMIX_UL:g} µL a well",
                note=f"a well is {WELL_UL:g} µL: {LYSATE_UL:g} µL lysate, {BARCODE_UL:g} µL "
                f"barcodes, {WATER_UL:g} µL water and {MASTERMIX_UL:g} µL master mix",
                citation=Citation("Qian SI", "Day 4.1"),
            ),
        ]
    else:
        made += [
            Material(
                f"{INDEX_WELLS}-well half-skirted PCR plate",
                supplier="USA Scientific",
                catalog=INDEX_CATALOG.split("#")[-1],
                amount=f"{len(one.index)}, one a quarter of a picked plate",
                citation=Citation("LevSeq", "step 2"),
            ),
            Material(
                "Barcoded index primer plate",
                storage="-20 °C",
                amount=f"{INDEX_PRIMER_UL:g} µL a well at 1 µM each, {INDEX_WELLS} forward and "
                f"{INDEX_WELLS} reverse",
                note="prepared once as lab stock, by its own protocol; a run calls for it. The "
                "SI's plate preparation ends at 0.1 µM, where its own protocol and the article "
                "both stamp 1 µM, which is what this takes",
                citation=Citation("LevSeq", "step 3"),
            ),
            Material(
                "ThermoPol Reaction Buffer",
                supplier="NEB",
                catalog=THERMOPOL_CATALOG.split("#")[-1],
                storage="-20 °C",
                amount="1 µL a reaction, 10X",
                citation=Citation("LevSeq", "step 1"),
            ),
            Material(
                "dNTP mix",
                supplier="NEB",
                catalog=DNTP_CATALOG.split("#")[-1],
                storage="-20 °C",
                amount="0.2 µL a reaction, 10 mM each",
                citation=Citation("LevSeq", "step 1"),
            ),
            Material(
                "Taq DNA Polymerase",
                supplier="NEB",
                catalog=TAQ_CATALOG.split("#")[-1],
                storage="-20 °C",
                amount="0.05 µL a reaction",
                note="LevSeq gives the volume and no source gives the stock, so the reaction "
                "carries no unit count",
                citation=Citation("LevSeq", "step 1"),
            ),
            Material(
                "DMSO, molecular biology grade",
                supplier="MP Biomedicals",
                catalog=DMSO_CATALOG.split("#")[-1],
                amount="0.4 µL a reaction, 4% (v/v) in the well",
                citation=Citation("LevSeq", "step 1"),
            ),
        ]
    return tuple(made)


def validation_equipment(one: Validation) -> tuple[str, ...]:
    """Return the hardware reading these designs back needs and no reagent table covers."""
    route = (
        (ACOUSTIC,)
        if one.route is ROUTE_A
        else (
            f"{MULTICHANNEL} on {INDEX_WELLS}-well spacing",
            f"Thermocycler taking a {INDEX_WELLS}-well plate",
        )
    )
    return (
        PICKER,
        *route,
        "Incubator at 37 °C",
        "A sequencer, and a demultiplexer that can check an address",
    )


def validation_steps(one: Validation) -> tuple[Step, ...]:
    """Return the steps that read these designs back, the route's own in the middle.

    The picking is shared and the calling is shared; between them sits the route's own marking.

    Examples
    --------
    >>> one = validation(ROUTE_B, (Design("a", 4),), 0)
    >>> for step in validation_steps(one):
    ...     print(step.title)
    Array 1 design(s) and grow
    Pick 4 colonies of each design
    Sample the picked plates into index plates
    Amplify each well with its own pair
    Pool and sequence
    Call every well
    """
    marking = _route_a_steps(one) if one.route is ROUTE_A else _route_b_steps(one)
    return (_array_step(one), _pick_step(one), *marking, _call_step(one))


def _array_step(one: Validation) -> Step:
    """Spot every design the floor reads back, which is where the read-back set is chosen."""
    floor = (
        "every design, which a floor of zero does"
        if one.floor == 0
        else f"every design in {one.floor} fragment(s) or more"
    )
    return Step(
        f"Array {len(one.designs)} design(s) and grow",
        instructions=(
            "Spot each design from its archive plate as its own spot on a 25 cm BioAssay plate.",
            "Grow overnight at 37 °C on 100 µg/mL carbenicillin.",
        ),
        expected=(
            f"One spot a design, at about {BIOASSAY_COLONIES:,} colonies a plate, which is the "
            "density picking wants.",
        ),
        notes=(
            f"This project reads back {floor}: {len(one.designs)} design(s). The rest stay "
            "polyclonal and are never read one design at a time.",
            "The archive is untouched: this reads a copy of it.",
        ),
        troubleshooting=(
            Troubleshooting(
                "A spot grows nothing",
                "That design has no clone to read. Re-transform it from the archive before "
                "anyone re-synthesises it.",
            ),
        ),
    )


def _pick_step(one: Validation) -> Step:
    """Pick the colonies, a quarter of a plate at a time, and print what each pick is worth."""
    sizes = ", ".join(f"{len(plate.seating)}" for plate in one.picked)
    return Step(
        f"Pick {one.colonies} colonies of each design",
        instructions=(
            f"Pick {one.colonies} colonies a design into {CULTURE_UL:g} µL low-salt LB with "
            f"carbenicillin, with the {PICKER}.",
            f"Fill one quarter of each {PICKED_WELLS}-well plate before starting the next.",
            "Grow overnight at 37 °C.",
        ),
        expected=(
            f"{one.wells} wells over {len(one.picked)} plate(s): {sizes} picked.",
            "Every colony of one design sits on one plate.",
        ),
        notes=(
            f"{one.colonies} colonies a design is Lund's anchor and the only measured one; four "
            "gave a clean copy of 343 of 458 genes.",
            *chances(one.designs),
            "A quarter at a time is what makes a part-filled plate give full plates downstream, "
            "and what keeps a mark off a plate that is mostly empty.",
        ),
        troubleshooting=(
            Troubleshooting(
                "A design in many fragments gives no clean colony",
                "The chances above say which designs that is likeliest for. Pick more colonies "
                "of those designs before the plates are poured, not after.",
            ),
        ),
    )


def _route_a_steps(one: Validation) -> tuple[Step, ...]:
    """Compress into 1536, barcode in lysate, then pool and sequence."""
    moves = tuple(
        compression(one.picked[at : at + PLATES_COMPRESSED], plate)
        for at, plate in zip(
            range(0, len(one.picked), PLATES_COMPRESSED), one.compressed, strict=True
        )
    )
    return (
        Step(
            f"Compress the picked plates into {len(one.compressed)} barcoding plate(s)",
            instructions=(
                "Invert the picked plates for 30 minutes so the cells gather at the meniscus.",
                f"Move {LYSATE_UL:g} µL of each well into the {COMPRESSED_WELLS}-well plate.",
            ),
            transfers=moves,
            expected=(f"{one.wells} wells of lysate, {PLATES_COMPRESSED} picked plates to one.",),
        ),
        Step(
            "Barcode each well in lysate",
            instructions=(
                f"Add one barcode from each of the {GROUPS} kit groups to every well, by the "
                "address that well's position gives.",
                f"Make each well up to {WELL_UL:g} µL with {BARCODE_UL:g} µL barcodes, "
                f"{WATER_UL:g} µL water and {MASTERMIX_UL:g} µL master mix.",
                "Run the ligation, then pool.",
            ),
            expected=(
                "One barcoded construct a well, carrying four barcodes chained head to tail.",
            ),
            notes=(
                "The address is worked out from where the well is, not looked up: three groups "
                "name the well and the fourth goes across the whole plate from a reservoir.",
                "The kit's sequences are not shipped. Read them from the copy you hold.",
            ),
        ),
        _sequencing_step(
            one,
            f"Pool every well and amplify the three primer pairs separately. {POOL_COLUMNS} "
            "miniprep columns, because one saturates.",
        ),
    )


def _route_b_steps(one: Validation) -> tuple[Step, ...]:
    """Sample a quarter at a time into index plates, amplify on the pair each well's address names."""
    at = 0
    moves: list[Transfer] = []
    for plate in one.picked:
        many = -(-len(plate.seating) // INDEX_WELLS)
        moves += sampling(plate, one.index[at : at + many])
        at += many
    return (
        Step(
            "Sample the picked plates into index plates",
            instructions=(
                f"Move {SAMPLE_UL:g} µL of each well into its {INDEX_WELLS}-well plate, one "
                "quarter of the picked plate a pass.",
            ),
            transfers=tuple(moves),
            expected=(f"{len(one.index)} index plate(s), {one.wells} reactions in all.",),
            notes=(
                "One well in four lines up under a head built for the smaller format, so a "
                "quarter moves in one pass.",
            ),
        ),
        Step(
            "Amplify each well with its own pair",
            instructions=(
                f"Add {INDEX_MIX_UL:g} µL of the master mix below to each well, which already "
                f"holds its {SAMPLE_UL:g} µL of culture.",
                f"Add {INDEX_PRIMER_UL:g} µL of the pair its address names from the prepared "
                f"primer plate, for {INDEX_PCR_UL:g} µL a well.",
                "Seal the plate, spin it down, and run the program below.",
            ),
            cautions=("Keep the polymerase on ice.",),
            tables=(index_pcr_reaction(one.wells),),
            programs=(index_pcr_program(),),
            expected=(
                "One barcoded amplicon a well. Well *n* of a plate takes forward mark *n*, and "
                "every well of one plate takes that plate's own reverse mark.",
                "One band a well, the design plus about 100 bp: the two marks, and the stretch "
                "between the primer sites and the reading frame.",
            ),
            notes=(
                f"{INDEX_WELLS} forward marks and {INDEX_WELLS} reverse reach "
                f"{ROUTE_B.capacity:,} wells, so the pairs already held cover far more than this "
                "run needs.",
                "The primer plate is built once as lab stock and a run calls for it; this "
                "protocol does not build one.",
                f"The first {INDEX_TOUCHDOWN_CYCLES} cycles touch down from "
                f"{INDEX_TOUCHDOWN_C[0]:g} to {INDEX_TOUCHDOWN_C[1]:g} °C, "
                f"{INDEX_TOUCHDOWN_STEP_C:g} °C a cycle, and {INDEX_PLATEAU_CYCLES} more run at "
                f"{INDEX_TOUCHDOWN_C[0]:g} °C.",
                f"Extension is {INDEX_EXTENSION_SECONDS // 60} minutes, what the authors ran for "
                "a gene below 1 kb. The published rule is at least one minute a kilobase, so a "
                "longer design wants it raised.",
                f"The mix is made {1 + INDEX_OVERAGE:g} times what the wells take, which is what "
                "LevSeq's own per-plate table spells out. Nothing says what the surplus is for.",
            ),
            troubleshooting=(
                Troubleshooting(
                    "A plate reads back far fewer wells than the others",
                    "LevSeq traces that to the amplification, not the sequencing: on a run of "
                    "ten plates three returned under 60% of their variants. Check each pool on "
                    "a gel before the library prep.",
                ),
            ),
            holes=(INDEX_MARKS, INDEX_TAQ_UNITS),
        ),
        _sequencing_step(one, "Pool each index plate on its own and clean the pool up."),
    )


def _sequencing_step(one: Validation, pooling_instruction: str) -> Step:
    """Pool the marked wells and sequence them, which both routes end their own stretch on."""
    return Step(
        "Pool and sequence",
        instructions=(pooling_instruction, "Sequence the pool."),
        expected=(
            f"Reads for {one.wells} wells, every well told from the rest by the marks it carries.",
        ),
        notes=(
            "Two plates on one flow cell are told apart by construction, because one axis of "
            "the address is the plate.",
        ),
    )


def _call_step(one: Validation) -> Step:
    """Demultiplex, judge each well on depth and then identity, and compact out what failed."""
    mark = "at least" if one.route.inclusive else "more than"
    return Step(
        "Call every well",
        instructions=(
            "Demultiplex the reads by address, checking each address rather than trusting a file.",
            "Call a consensus a well, then compare it with that well's design base for base.",
            "Reformat, compacting out the wells that failed.",
        ),
        expected=(
            f"A well read {mark} {one.route.wanted_reads} times is called; route "
            f"{one.route.name} tolerates {one.route.tolerable_reads} and warns between the two.",
            "A pass matches across the whole designed region: both entry overhangs, the "
            "fragment, the stuffer and the barcode.",
        ),
        notes=(
            "A well read too thinly carries no verdict and is not a failure, so reformatting "
            "does not compact it out.",
            "More than one consensus is mixed, and mixed fails.",
            "A design with no passing well is picked again from the same archive spot before "
            "anyone re-synthesises it.",
        ),
        troubleshooting=(
            Troubleshooting(
                "A well's call is the right protein and the wrong barcode",
                "It fails. A barcode that no longer names its member cannot be put right by the "
                "linkage read later.",
            ),
        ),
    )


def _kit_rows(lines: Sequence[str]) -> tuple[Barcode, ...]:
    """Parse the barcode rows, refusing a file that is not the kit.

    Raises
    ------
    ValueError
        Naming what is wrong with it.
    """
    reader = csv.DictReader(lines, delimiter="\t")
    if missing := sorted(set(KIT_COLUMNS) - set(reader.fieldnames or ())):
        raise ValueError(f"missing column(s) {', '.join(missing)}")
    found = [
        Barcode(
            row["name"],
            group=int(row["group"]),
            index=int(row["index"]),
            overhang5=row["overhang5"].upper(),
            overhang3=row["overhang3"].upper(),
            umi=row["umi"].upper(),
            sequence=row["final_seq"].upper(),
        )
        for row in reader
    ]
    _check_kit(found)
    return tuple(sorted(found, key=lambda one: (one.group, one.index)))


def _check_kit(found: Sequence[Barcode]) -> None:
    """Refuse a kit that is not four groups of 24 with distinct UMIs chaining in one order.

    Raises
    ------
    ValueError
        Naming what is wrong with it.
    """
    if len(found) != GROUPS * GROUP_SIZE:
        raise ValueError(f"{len(found)} rows, not {GROUPS * GROUP_SIZE}")
    if len({one.umi for one in found}) != len(found):
        raise ValueError("two rows share a UMI")
    for group in range(1, GROUPS + 1):
        members = [one for one in found if one.group == group]
        if sorted(one.index for one in members) != list(range(1, GROUP_SIZE + 1)):
            raise ValueError(f"group {group} is not indexed 1 to {GROUP_SIZE}")
        ends = {(one.overhang5, one.overhang3) for one in members}
        if len(ends) != 1:
            raise ValueError(f"group {group} chains on {len(ends)} different overhang pairs")
    # Head to tail, so a barcode cannot land in the wrong position and none can be left out.
    # Checked strand-agnostically: a copy held on the map strand spells the chain reversed and
    # complemented, and the adjacency is the same either way.
    for group in range(1, GROUPS):
        left = next(one for one in found if one.group == group)
        right = next(one for one in found if one.group == group + 1)
        if left.overhang3 != right.overhang5:
            raise ValueError(
                f"group {group} ends on {left.overhang3} and group {group + 1} begins on "
                f"{right.overhang5}, so the four do not chain head to tail"
            )


#: Where the numbers above come from, ready for a protocol's reference list.
REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Qian, Z. et al. (2026) Accelerating protein design by scaling experimental "
        "characterization. Nat. Commun., for the DMX barcode kit, the 1536-well barcoding in "
        "lysate, and a consensus called above 150 reads",
        url="https://doi.org/10.1038/s41467-026-76740-5",
    ),
    Reference(
        "Long, Y. et al. (2025) LevSeq: rapid generation of sequence-function data for "
        "directed evolution and machine learning, for index PCR marking a well on one "
        "barcoded primer pair, its reaction and touchdown cycling, twenty reads wanted and "
        "ten tolerable",
        url="https://doi.org/10.1021/acssynbio.4c00625",
    ),
    Reference(
        "Lund, S. et al. (2024) Highly parallelized construction of DNA from low-cost "
        "oligonucleotide mixtures using Data-optimized Assembly Design and Golden Gate. ACS "
        "Synth. Biol. 13, 745-751, for four colonies giving a clean copy of 343 of 458 genes, "
        "and the clean rate against the number of fragments a design is built from",
        url="https://doi.org/10.1021/acssynbio.3c00694",
    ),
)

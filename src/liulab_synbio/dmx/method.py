"""How a well is marked, how deeply it is read, and what reading a plate back takes.

Two routes mark a well and one judgement reads them, and `liulab_synbio.dmx` says what each
one is. The picking, the pass rule and the reformat are shared; the marking step, the plate and
the depth floor are the route's own. The two floors are not a strict and a lenient setting of
one scale: one is where consensus calling starts, the other where a reader stops trusting a
well, over different amplification and different read filters.
`docs/research/route-choice.md` section 4 compares them.

Every number a read-back prints is here, with the document it was read from beside it.
`liulab_synbio.dmx.steps` writes them up as steps.
"""

from collections.abc import Sequence
from dataclasses import KW_ONLY, dataclass

from liulab_mbio.bench import plates
from liulab_mbio.bench.materials import material
from liulab_mbio.bench.readback import WellVerdict, clean_colony_chance, identity_check
from liulab_mbio.checks import Check, Status
from liulab_mbio.protocol.model import (
    Citation,
    Component,
    Hole,
    Incubation,
    Item,
    Material,
    Plate,
    ReactionTable,
    Reference,
    Source,
    Stage,
    ThermocyclerProgram,
    Topic,
    Transfer,
    Well,
)
from liulab_synbio.dmx.kit import GROUP_SIZE, GROUPS, Kit, KitBarcode

#: Colonies picked per design. Four is Lund's only measured anchor, and it is a build's input
#: rather than a constant of the method.
COLONIES_PER_DESIGN = 4

#: The plate each stage uses. Colonies are picked into 384-well plates and four of those are
#: compressed into one 1536-well plate, which is why the format parameter has to reach 1536.
#: Qian SI Day 3 and Day 4.1. Index PCR amplifies in the half-skirted 96-well PCR plate its
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

#: What a quarter of a picked plate is sampled into an index plate with, and how much of
#: each well goes. LevSeq SI step 4. An acoustic handler does the same move where one is booked.
SAMPLE_UL = 1.0
MULTICHANNEL = "multichannel pipette"

#: The index PCR itself: a colony PCR straight from overnight culture, cycled with a touchdown.
#: The master mix is one well's share of the mix LevSeq states per plate, and the cycling is its
#: thermal-cycler table with the two loop lines spelled out: ten touchdown cycles and 25 more,
#: 35 in all. ``docs/research/route-b-index-pcr.md`` section 2.
INDEX_PCR_UL = 10.0
INDEX_MIX_UL = 7.0
INDEX_PRIMER_UL = 2.0
INDEX_TAQ_UL = 0.05
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

#: What the Taq stock is worth, which LevSeq never states and NEB's own specification for M0267
#: does. ``docs/research/route-b-index-pcr.md`` section 9.
TAQ_STOCK_UNITS_UL = 5.0
TAQ_SOURCE_KEY = "PS-M0267"

#: Colonies a 25 cm BioAssay plate carries before picking gets hard. Qian SI Day 2.
BIOASSAY_COLONIES = 2500
BIOASSAY_CATALOG = "Corning #431111"

#: Columns a pooled 1536-well plate takes, because one saturates. Qian SI Day 4.1.
POOL_COLUMNS = 2

#: The documents these numbers were read from.
SOURCES: dict[str, Source] = {
    "Qian SI": Source(
        "Qian et al. 2026, Supplementary Information",
        edition="Nat. Commun. 10.1038/s41467-026-76740-9",
        read_as="held under reference_docs/",
        date="2026-10-06",
    ),
    "LevSeq": Source(
        "Long, Y. et al. 2025, LevSeq, Supporting Information",
        read_as="held under reference_docs/",
        date="2026-10-06",
    ),
    TAQ_SOURCE_KEY: Source(
        "New England Biolabs, Product Specification: Taq DNA Polymerase with ThermoPol Buffer",
        edition="PS-M0267S/L/X/E v2.0, effective 12 Feb 2020",
        url="https://www.neb.com/en/-/media/catalog/specifications/m/0/m0267s_l_x_e_v2.pdf",
        read_as="plain curl",
        date="2026-10-08",
        note="docs/research/route-b-index-pcr.md",
    ),
}


@dataclass(frozen=True, slots=True)
class Route:
    """One way of marking a well, and how deep its read has to be before anyone calls it.

    The two routes differ in their marking step, their plate and their floor; the picking, the
    pass rule and the reformat are shared. Which one a build runs is the build's choice.

    Parameters
    ----------
    name
        What the marking step does: ``"barcode ligation"`` or ``"index PCR"``. It is both the
        name the page prints and the key a build names a route by.
    marking
        What the step does, in a few words.
    well_axes
        How many marks each axis addressing the well carries, in chain order. The well's address
        factorises across them.
    plate_axis
        How many marks the one axis carrying the plate carries.
    wanted_reads
        A well is called **above** this mark. Qian publishes a consensus depth above 150 and
        LevSeq's SI checklist an alignment count above 20, so neither route calls a well landing
        exactly on it. A build may raise it.
    tolerable_reads
        The mark at or above which a well still carries a verdict; `None` where the route
        publishes one number and has no warn band. LevSeq's SI routes a well at ``<=20`` into a
        proceed-with-validation branch and its Figure S1E puts the floor of detection at 10
        reads, so a well at exactly 10 warns rather than going unjudged.
    source
        The key of the `SOURCES` entry the depths were read from.
    """

    name: str
    _: KW_ONLY
    marking: str
    well_axes: tuple[int, ...]
    plate_axis: int
    wanted_reads: int
    tolerable_reads: int | None
    source: str

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
ROUTE_LIGATION = Route(
    "barcode ligation",
    marking="ligate one barcode from each of the four kit groups into the construct, in lysate",
    well_axes=(GROUP_SIZE, GROUP_SIZE, GROUP_SIZE),
    plate_axis=GROUP_SIZE,
    wanted_reads=150,
    tolerable_reads=None,
    source="Qian SI",
)

#: One barcoded primer pair a well, 96 forward marks addressing the well and 96 reverse the
#: plate: 9,216 wells on the 192 primers already held. Twenty reads wanted, ten tolerable,
#: twenty to exceed and ten to reach, so a well at exactly twenty is called and warns, where the
#: SI's suboptimal branch puts it. LevSeq pairs the floor with a second criterion, a mean error
#: below 10%, which nothing here judges — a named hole: it is a mean over per-position base
#: counts, and a well reaches this package as a read count and a consensus call.
ROUTE_INDEX_PCR = Route(
    "index PCR",
    marking="amplify each well with one barcoded primer pair",
    well_axes=(INDEX_WELLS,),
    plate_axis=INDEX_WELLS,
    wanted_reads=20,
    tolerable_reads=10,
    source="LevSeq",
)

#: Both routes, by the name a build names one with.
ROUTES: dict[str, Route] = {
    ROUTE_LIGATION.name: ROUTE_LIGATION,
    ROUTE_INDEX_PCR.name: ROUTE_INDEX_PCR,
}

#: The job both routes are a way of doing, in the words a bench reader uses for it. A run
#: offering both titles its guidance with this, which is what puts each way's page one link
#: from what to weigh.
READ_BACK = "read every well back"


def route_choice() -> Topic:
    """Return what to weigh before doing one of the two routes, for a run offering both.

    Every figure is `docs/research/route-choice.md`'s, which is also what says which comparisons
    no source carries. Those are written as cautions and never as figures: there is no published
    cost for index PCR, no cost for either route per well, and no library size at which one
    overtakes the other.
    """
    return Topic(
        READ_BACK,
        (
            "Both routes mark every well so that sequencing says which well a read came from, "
            "and both end with a pass or a fail for each well. Do one of them, never both.",
            "How many wells each reaches is not what picks between them, and it does not rank "
            "them the way the bench's rule of thumb does. Barcode ligation addresses more than "
            "330,000 wells from one plate of barcodes; index PCR addresses 9,216. Most runs sit "
            "well below either.",
            "Barcode ligation pays a fixed price once, whatever the library's size: $695 of "
            "sequencing, and five days. Shared over 100 designs that is $19.52 a design, and "
            "over 2,000 designs $3.38. The more designs share it, the less each one carries. "
            "Nobody has published the same table for index PCR, so no price for it, and no "
            "price per well for either route, is stated anywhere here.",
            "Index PCR runs one thermocycled reaction in every well, so thousands of wells need "
            "many thermocyclers. Barcode ligation holds one temperature, and a whole 1536-well "
            "plate marks in a single incubator. That is the only direct comparison of the two "
            "in print, and Qian and colleagues wrote it about the route their own replaced; no "
            "independent measurement of one against the other was found.",
            "Doing neither is also a choice. Every well is then enriched for the design it was "
            "built from, and nothing firmer than that can be said: the one published purity "
            "figure, 89.3% of 929 wells at 90% clonal purity or better, was measured on short "
            "arrayed fragments under 314 bases rather than on cargo from an oligo pool, and that "
            "paper says polyclonality grows with length. Nothing measured says the material is "
            "even. To know what a well holds, read it back.",
        ),
    )


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
    >>> address(ROUTE_LIGATION, plate=0, well=25).marks
    (2, 2, 1, 1)
    >>> address(ROUTE_INDEX_PCR, plate=2, well=0).marks
    (1, 3)
    """
    if plate < 0 or well < 0:
        raise ValueError(f"a well's address counts from zero, got plate {plate} and well {well}")
    if well >= route.wells_per_plate:
        raise ValueError(
            f"the {route.name} route tells {route.wells_per_plate} wells of a plate apart, and "
            f"well "
            f"{well} is past the last of them"
        )
    if plate >= route.plate_axis:
        raise ValueError(
            f"the {route.name} route tells {route.plate_axis} plates apart, and plate {plate} is "
            "past the last of them"
        )
    marks: list[int] = []
    left = well
    for axis in route.well_axes:
        marks.append(left % axis + 1)
        left //= axis
    return Address(route, plate=plate, well=well, well_marks=tuple(marks), plate_mark=plate + 1)


def barcodes_for(kit: Kit, one: Address) -> tuple[KitBarcode, ...]:
    """Return the four kit barcodes an address names, in chain order.

    Raises
    ------
    ValueError
        If the address was not worked out on the barcode ligation route, which the kit marks.
    """
    if one.route is not ROUTE_LIGATION:
        raise ValueError(
            f"the DMX barcode kit marks a well on the {ROUTE_LIGATION.name} route, and this "
            f"address was worked out on the {one.route.name} route"
        )
    return tuple(kit.at(group, mark) for group, mark in enumerate(one.marks, start=1))


def depth_check(route: Route, reads: int, *, wanted: int | None = None) -> Check:
    """Return the verdict on whether a well was read deeply enough to call.

    A well below the route's tolerable depth carries **no verdict**, not a failure: it is read
    again or picked again, and reformatting does not compact it out. From the tolerable mark up
    to the wanted one it warns; above the wanted one it passes. A route with no tolerable mark
    has no warn band. `wanted` raises the route's own mark, which a build may do and may not
    lower.

    Raises
    ------
    ValueError
        If `reads` is negative, or `wanted` is below the route's own mark.

    Examples
    --------
    >>> depth_check(ROUTE_INDEX_PCR, 12).status, depth_check(ROUTE_INDEX_PCR, 4).status
    ('warn', None)
    >>> depth_check(ROUTE_INDEX_PCR, 10).status
    'warn'
    """
    if reads < 0:
        raise ValueError(f"a well cannot carry {reads} reads")
    mark = route.wanted_reads if wanted is None else wanted
    if mark < route.wanted_reads:
        raise ValueError(
            f"the {route.name} route wants {route.wanted_reads} reads a well, and that mark "
            f"may be raised, not lowered to {mark}"
        )
    status: Status | None = None
    if reads > mark:
        status = "pass"
    elif route.tolerable_reads is not None and reads >= route.tolerable_reads:
        status = "warn"
    tolerated = "" if route.tolerable_reads is None else f" and tolerates {route.tolerable_reads}"
    return Check(
        "reads_per_well",
        status,
        float(reads),
        f"the {route.name} route wants more than {mark} reads a well{tolerated}; below that no "
        f"is deep enough to call",
    )


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
    chance worth paying to check, and that is the build's own call.

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


def refuse_unclonal(cargo: str, *, clonal: bool) -> None:
    """Refuse to read back cargo nobody can pick a colony from.

    DMX takes any cargo and a design's name does not say which kind it is, so the caller names
    its cargo and says whether it arrives as colonies.

    Raises
    ------
    ValueError
        If the cargo is not clonal.
    """
    if clonal:
        return
    raise ValueError(
        "DMX reads clonal material: it spots each design from its archive plate, picks colonies "
        "off it and grows every pick under the destination's own selection. There is no colony "
        f"to pick and no marker to select on in {cargo}"
    )


def selected_on(selection: str) -> str:
    """Return what a plate of this read's transformants carries, named or left to the record.

    The drug is the vector's, read off its marker by the caller: this method's own DMX vector is
    not the one the published protocol was written for, and naming that one sends a reader to an
    empty plate.
    """
    return selection or "the vector's own antibiotic"


def picked_plate(name: str, colonies: int, selection: str = "") -> Plate:
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
    >>> picked.wells, len(picked.labels), picked.labels["B1"]
    (384, 288, 'quarter 3')
    """
    if colonies > PICKED_WELLS:
        raise ValueError(f"{colonies} colonies do not fit a {PICKED_WELLS}-well plate")
    labels: dict[str, str] = {}
    left = colonies
    for number, quarter in enumerate(plates.interleave(PICKED_WELLS, INDEX_WELLS), start=1):
        labels |= dict.fromkeys(quarter[: max(left, 0)], f"quarter {number}")
        left -= len(quarter)
    return plates.plate(
        name,
        PICKED_WELLS,
        catalog=PICKED_CATALOG,
        holds=f"one picked colony each in {CULTURE_UL:g} µL low-salt LB with "
        f"{selected_on(selection)}",
        labels=labels,
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
    """Return one index plate: each well amplified with the pair its address names.

    Raises
    ------
    ValueError
        If more samples are given than the plate holds, or `plate` is past the last the route
        tells apart.

    Examples
    --------
    >>> index_plate("index 1", 3, plate=1).labels["A2"]
    'forward 2, reverse 2'
    """
    if samples > INDEX_WELLS:
        raise ValueError(f"{samples} samples do not fit a {INDEX_WELLS}-well plate")
    labels = {}
    names = plates.plate(name, INDEX_WELLS).well_names
    for well in range(samples):
        one = address(ROUTE_INDEX_PCR, plate=plate, well=well)
        labels[names[well]] = f"forward {one.well_marks[0]}, reverse {one.plate_mark}"
    return plates.plate(
        name,
        INDEX_WELLS,
        catalog=INDEX_CATALOG,
        holds="one barcoded PCR each, the pair its own address names",
        labels=labels,
        note=f"plate mark {plate + 1}; {samples} of {INDEX_WELLS} wells used",
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
        for name in (plate.labels or dict.fromkeys(plate.well_names))
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
    """Return the index PCR move out of one picked plate: a quarter of it into each index plate.

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
    filled = picked.labels or dict.fromkeys(picked.well_names)
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
        The route the build named.
    designs
        The designs read back, which `validated` chose from the build's floor.
    floor
        That floor, carried so the protocol can print it beside each design's chance.
    colonies
        Colonies picked per design.
    selection
        What to plate on, read off the vector's own marker by the caller. Empty where the record
        annotates none, and every plate then says so rather than naming a drug.
    index_plate
        What the lab calls the prepared plate of barcoded primer pairs the index PCR route takes
        one pair of per well. Empty, the step says no more than that a prepared plate is wanted,
        and `INDEX_MARKS` stands.
    """

    route: Route
    designs: tuple[Design, ...]
    _: KW_ONLY
    floor: int
    colonies: int = COLONIES_PER_DESIGN
    selection: str = ""
    index_plate: str = ""

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
            picked_plate(f"picked {number}", size, self.selection)
            for number, size in enumerate(sizes, start=1)
        )

    @property
    def compressed(self) -> tuple[Plate, ...]:
        """The barcode ligation plates, `PLATES_COMPRESSED` picked plates compressed into one."""
        self._only(ROUTE_LIGATION)
        many = -(-len(self.picked) // PLATES_COMPRESSED)
        return tuple(compressed_plate(f"lysate {number}") for number in range(1, many + 1))

    @property
    def index(self) -> tuple[Plate, ...]:
        """The index plates, one a filled quarter, marked in one run across the run."""
        self._only(ROUTE_INDEX_PCR)
        made: list[Plate] = []
        for one in self.picked:
            left = len(one.labels)
            while left > 0:
                made.append(
                    index_plate(f"index {len(made) + 1}", min(left, INDEX_WELLS), plate=len(made))
                )
                left -= INDEX_WELLS
        return tuple(made)

    @property
    def plates(self) -> tuple[Plate, ...]:
        """Every plate this read pours, so each well a transfer names has one drawn for it."""
        rest = self.index if self.route is ROUTE_INDEX_PCR else self.compressed
        return (*self.picked, *rest)

    def _only(self, route: Route) -> None:
        """Refuse a plate the other route never pours.

        Raises
        ------
        ValueError
            If this validation runs the other route.
        """
        if self.route is not route:
            raise ValueError(
                f"only the {route.name} route pours this plate, and this read runs the "
                f"{self.route.name} route"
            )


def validation(
    route: Route,
    designs: Sequence[Design],
    floor: int | None,
    *,
    colonies: int = COLONIES_PER_DESIGN,
    selection: str = "",
    index_plate: str = "",
) -> Validation | None:
    """Return what reading `designs` back on `route` takes, or `None` where the floor reads none.

    `None` is the answer for a build that states no floor and for one whose floor is above
    every design: either way nothing is read, and a protocol then carries no validation at all.

    `selection` is what the caller read off the vector's marker. Left empty, every plate says
    the vector's own antibiotic rather than naming one this read cannot know.

    `index_plate` is what the lab calls its prepared plate of barcoded primer pairs, which only
    the index PCR route takes.

    Raises
    ------
    ValueError
        If the floor is negative, or fewer than one colony a design is picked.

    Examples
    --------
    >>> validation(ROUTE_INDEX_PCR, (Design("one", 4),), 2).wells
    4
    >>> validation(ROUTE_INDEX_PCR, (Design("one", 4),), None) is None
    True
    """
    if floor is None:
        return None
    read = validated(designs, floor)
    if not read:
        return None
    return Validation(
        route,
        read,
        floor=floor,
        colonies=colonies,
        selection=selection,
        index_plate=index_plate,
    )


#: What a read-back hands on, under the names a chain matches by. They are this method's and not
#: one caller's: every run that reads a plate back leaves one clone a well and a verdict a well.
PICKED_PLATE = Item("clonal picked plate", "one well a picked colony, each well one clone")
WELL_CALLS = Item("well calls", "a pass or a fail a well, and which design each well holds")


def marking_stock(one: Validation | None) -> Item | None:
    """Return the lab stock the chosen route marks with, which no protocol of a run makes.

    `None` where nothing is read back and no well is marked at all. Both routes take stock the
    lab prepares once and a run calls for, so it is a run's input wherever a read-back stands.
    The index plate a build named is said beside it; the item keeps its own name, which is what
    a chain hands it over by.
    """
    if one is None:
        return None
    if one.route is ROUTE_LIGATION:
        return Item(
            "DMX barcode kit",
            "the lab's own barcoding plasmids, one group a picked plate",
            storage="-20 °C",
        )
    named = f"{one.index_plate}, " if one.index_plate else "the lab's own index primers, "
    return Item(
        "Barcoded index primer plate",
        f"{named}prepared once and called for by a run",
        spec=("1 µM each",),
        storage="-20 °C",
    )


#: The index PCR marks, which the package holds none of. The published annealing regions bind the
#: DMX vector verbatim; which 192 index sequences sit on their 5' ends is a plate the lab holds,
#: as the barcode ligation kit is. ``docs/research/route-b-index-primers.md`` §7.
INDEX_MARKS = Hole(
    "IDX1",
    "no index mark set is named for the barcoded primer pairs",
    "lab",
    where="index PCR, the pair marking one well",
    filled_by="the prepared primer plate the lab holds",
)


def index_pcr_reaction(reactions: int = 1) -> ReactionTable:
    """Return LevSeq's index PCR, one `INDEX_PCR_UL` µL reaction a well.

    The first five lines are the master mix, each one well's share of the mix LevSeq states per
    plate; at `INDEX_OVERAGE` they scale back to that table. The pair and the culture go in a
    well at a time, because each well takes its own pair. LevSeq never gives the Taq stock, so
    the unit count comes from the supplier's own specification.

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
            Component(
                "Taq DNA Polymerase",
                INDEX_TAQ_UL,
                stock=f"{TAQ_STOCK_UNITS_UL:g} U/µL",
                final=f"{INDEX_TAQ_UL * TAQ_STOCK_UNITS_UL:g} units",
                citation=mix,
            ),
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
    """Return LevSeq's touchdown program, the SI's two loop lines as two stages.

    The touchdown is one stage of `INDEX_TOUCHDOWN_CYCLES` cycles whose annealing step drops
    `INDEX_TOUCHDOWN_STEP_C` a cycle, across `INDEX_TOUCHDOWN_C`; `INDEX_PLATEAU_CYCLES` more
    anneal and extend together at the top of it.

    Examples
    --------
    >>> touchdown = index_pcr_program().stages[1]
    >>> anneal = touchdown.incubations[1]
    >>> anneal.temperature_c, anneal.delta_c, anneal.last_c(touchdown.cycles)
    (68.0, -0.5, 63.5)
    """
    cite = Citation("LevSeq", "thermal cycler table")
    top = INDEX_TOUCHDOWN_C[0]
    denature = Incubation("Denature", INDEX_DENATURE_C, 20, citation=cite)
    return ThermocyclerProgram(
        (
            Stage((Incubation("Initial denaturation", INDEX_DENATURE_C, 300, citation=cite),)),
            Stage(
                (
                    denature,
                    Incubation("Anneal", top, 20, delta_c=-INDEX_TOUCHDOWN_STEP_C, citation=cite),
                    Incubation("Extend", top, INDEX_EXTENSION_SECONDS, citation=cite),
                ),
                cycles=INDEX_TOUCHDOWN_CYCLES,
            ),
            Stage(
                (
                    denature,
                    Incubation("Anneal and extend", top, INDEX_EXTENSION_SECONDS, citation=cite),
                ),
                cycles=INDEX_PLATEAU_CYCLES,
            ),
            Stage((Incubation("Final extension", top, 300, citation=cite),)),
            Stage((Incubation("Hold", 4.0, None, citation=cite),)),
        ),
        title="Index PCR",
    )


def validation_materials(one: Validation) -> tuple[Material, ...]:
    """Return what reading these designs back consumes, beyond the designs themselves."""
    made = [
        Material(
            "25 cm BioAssay plate",
            supplier="Corning",
            catalog=BIOASSAY_CATALOG.split("#")[-1],
            amount="one spot a design",
            note=f"{selected_on(one.selection)}; about {BIOASSAY_COLONIES:,} colonies a plate",
            citation=Citation("Qian SI", "Day 2"),
        ),
        Material(
            f"{PICKED_WELLS}-well culture plate",
            supplier="Beckman Coulter",
            catalog=PICKED_CATALOG.split("#")[-1],
            amount=f"{len(one.picked)}, one a {PICKED_WELLS} wells",
            note=f"{CULTURE_UL:g} µL low-salt LB with {selected_on(one.selection)} a well",
            citation=Citation("Qian SI", "Day 3"),
        ),
    ]
    if one.route is ROUTE_LIGATION:
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
            material(
                "Taq DNA Polymerase",
                supplier="NEB",
                catalog=TAQ_CATALOG.split("#")[-1],
                storage="-20 °C",
                amount=f"{INDEX_TAQ_UL:g} µL a reaction",
                note=f"{TAQ_STOCK_UNITS_UL:g} U/µL, so one reaction takes "
                f"{INDEX_TAQ_UL * TAQ_STOCK_UNITS_UL:g} units",
                citation=Citation(TAQ_SOURCE_KEY, "Concentration"),
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
        if one.route is ROUTE_LIGATION
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


#: Where the numbers above come from, ready for a protocol's reference list.
REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Qian, Z. et al. (2026) Accelerating protein design by scaling experimental "
        "characterization. Nat. Commun., for the DMX barcode kit, the 1536-well barcoding in "
        "lysate, and a consensus called above 150 reads",
        url="https://doi.org/10.1038/s41467-026-76740-9",
    ),
    Reference(
        "Long, Y. et al. (2025) LevSeq: rapid generation of sequence-function data for "
        "directed evolution and machine learning, for index PCR marking a well on one "
        "barcoded primer pair, its reaction and touchdown cycling, twenty reads wanted and "
        "ten tolerable, twenty to exceed and ten to reach",
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

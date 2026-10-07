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
from collections.abc import Sequence
from dataclasses import KW_ONLY, dataclass
from pathlib import Path

from liulab_mbio.bench import plates
from liulab_mbio.checks import Check, Status, worst
from liulab_mbio.protocol.model import Citation, Plate, Reference, Source, Transfer, Vessel, Well

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
#: Qian SI Day 3 and Day 4.1. Route B amplifies in the 96-well plate its barcodes address.
PICKED_WELLS = 384
COMPRESSED_WELLS = 1536
INDEX_WELLS = 96
PICKED_CATALOG = "Beckman Coulter #c74290"
COMPRESSED_CATALOG = "Greiner #782270"

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

#: What moves the liquid, and what the colonies are picked with. Qian SI Days 3 and 4.1.
ACOUSTIC = "ECHO 525 acoustic liquid handler"
PICKER = "QPix XE Microbial Colony Picker"

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
        The depth at and above which a well is called. A project may raise it.
    tolerable_reads
        The depth below which no read is deep enough to call, so the well carries no verdict.
        Equal to `wanted_reads` where the route publishes one number rather than two.
    inclusive
        Whether a depth equal to a mark meets it. Qian publishes a consensus depth *above* 150
        and LevSeq twenty reads wanted, so the two routes read their own marks differently and
        neither is rewritten to suit the other.
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
#: plate: 9,216 wells on the 192 primers already held. Twenty reads wanted, ten tolerable.
ROUTE_B = Route(
    "B",
    marking="amplify each well with one barcoded primer pair",
    well_axes=(INDEX_WELLS,),
    plate_axis=INDEX_WELLS,
    wanted_reads=20,
    tolerable_reads=10,
    source="LevSeq",
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


def picked_plate(name: str, colonies: int) -> Plate:
    """Return one plate of picked colonies: `colonies` wells of selective medium.

    Raises
    ------
    ValueError
        If more colonies are picked than the format holds.

    Examples
    --------
    >>> picked_plate("picked 1", 300).wells
    384
    """
    if colonies > PICKED_WELLS:
        raise ValueError(f"{colonies} colonies do not fit a {PICKED_WELLS}-well plate")
    return plates.plate(
        name,
        PICKED_WELLS,
        catalog=PICKED_CATALOG,
        holds=f"one picked colony each in {CULTURE_UL:g} µL low-salt LB with carbenicillin",
        note=f"{colonies} of {PICKED_WELLS} wells picked",
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
        "barcoded primer pair, twenty reads wanted and ten tolerable",
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

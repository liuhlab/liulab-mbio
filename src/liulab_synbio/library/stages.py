"""What this method asks of the bench model: its plates, its moves, its rules and its holes.

Every number here is quoted from ``docs/research/bench-numbers.md``, which names the document
each came from. Nothing is invented: where the method needs a number nobody published, a `Hole`
stands in its place and says which ticket it is routed to.

The plate model, the moves and the rules are all `liulab_mbio.bench`'s. This module only says
which format each stage uses, what it moves, and which materials carry which rule — a different
method would call the same functions with its own numbers.
"""

from collections.abc import Sequence

from liulab_mbio.bench import materials, plates
from liulab_mbio.protocol.model import (
    Citation,
    Hole,
    Material,
    Plate,
    Source,
    Transfer,
    Vessel,
    Well,
)

#: The plate each stage uses. Colonies are picked into 384-well plates and four of those are
#: compressed into one 1536-well plate, which is why the format parameter has to reach 1536.
#: Qian SI Day 3 and Day 4.1.
PICKED_WELLS = 384
COMPRESSED_WELLS = 1536
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

#: The documents these numbers were read from, beside the ones a material brings with it.
SOURCES: dict[str, Source] = {
    "Qian SI": Source(
        "Qian et al. 2026, Supplementary Information",
        edition="Nat. Commun. 10.1038/s41467-026-76740-5",
        read_as="held under reference_docs/",
        date="2026-10-06",
    ),
    **materials.SOURCES,
}

#: The ligase this method's round runs on, and the buffer it runs in. Both carry their own
#: rules, so neither can be used in a protocol that does not show them: never add PEG to a T7
#: ligase reaction, and do not heat-inactivate it in PEG.
LIGASE = materials.material(
    "T7 DNA Ligase",
    supplier="NEB",
    catalog="#M0318L",
    storage="-20 °C",
    note="3,000,000 U/mL. Refuses blunt ends, which is what removes escapees.",
    citation=Citation("M0318", "reaction conditions"),
)
LIGASE_BUFFER = materials.material(
    "StickTogether DNA Ligase Buffer",
    supplier="NEB",
    catalog="#B0535S",
    storage="-20 °C",
    note="Supplied 2x. 1x is 7.5% PEG 6000, well under the 20% that would force blunt ligation.",
    citation=Citation("M0318", "reaction conditions"),
)

#: What this method cannot write completely. Each is a number nobody published, routed to the
#: ticket that would decide it; the ids are the research note's own, so a reference still
#: resolves.
HOLES: tuple[Hole, ...] = (
    Hole(
        "H14",
        "no reads-per-well pass mark anywhere",
        "unpublished",
        where="read-out, what counts as a pass",
        filled_by="the lab, from a first run",
        issue="liuhlab/liulab-mbio#264",
    ),
    Hole(
        "H22",
        "no colony count and no selection antibiotic for a round",
        "undecided",
        where="each round, plating",
        filled_by="the vector's own marker, and the coverage the package computes",
        issue="liuhlab/liulab-mbio#264",
    ),
    Hole(
        "H23",
        "no published document describes the split digest; the one source is the paper itself",
        "unpublished",
        where="the split digest, units per reaction",
        filled_by="nothing; it is the method's own",
        issue="liuhlab/liulab-mbio#264",
    ),
    Hole(
        "H24",
        "every mass in the final assembly is unsourced, the SPRI ratio with them",
        "undecided",
        where="final assembly, DNA in",
        filled_by="a pilot, or carrying the round's own numbers across and saying so",
        issue="liuhlab/liulab-mbio#264",
    ),
)


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


def method_materials() -> tuple[Material, ...]:
    """Return the materials whose own parameters and rules this method depends on."""
    return (LIGASE, LIGASE_BUFFER)

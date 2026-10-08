"""The steps that read a plate of wells back, the route's own marking in the middle.

The picking and the calling are shared, and between them sits the step that marks a well:
barcode ligation in lysate on one route, index PCR on the other. The numbers these steps print
are `liulab_synbio.dmx.method`'s, so a page and the plates it pours cannot disagree.

A caller supplies the figure the marking step draws, because the records it names belong to the
run and this module holds none.
"""

from collections import Counter
from collections.abc import Sequence

from liulab_mbio.bench.readback import clean_colony_chance
from liulab_mbio.protocol.model import Figure, Step, Transfer, Troubleshooting
from liulab_synbio.dmx.kit import GROUPS
from liulab_synbio.dmx.method import (
    BARCODE_UL,
    BIOASSAY_COLONIES,
    COMPRESSED_WELLS,
    CULTURE_UL,
    INDEX_EXTENSION_SECONDS,
    INDEX_MARKS,
    INDEX_MIX_UL,
    INDEX_OVERAGE,
    INDEX_PCR_UL,
    INDEX_PLATEAU_CYCLES,
    INDEX_PRIMER_UL,
    INDEX_TOUCHDOWN_C,
    INDEX_TOUCHDOWN_CYCLES,
    INDEX_TOUCHDOWN_STEP_C,
    INDEX_WELLS,
    LYSATE_UL,
    MASTERMIX_UL,
    PICKED_WELLS,
    PICKER,
    PLATES_COMPRESSED,
    POOL_COLUMNS,
    ROUTE_INDEX_PCR,
    ROUTE_LIGATION,
    SAMPLE_UL,
    WATER_UL,
    WELL_UL,
    Design,
    Validation,
    compression,
    index_pcr_program,
    index_pcr_reaction,
    sampling,
    selected_on,
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


def validation_steps(one: Validation, *, marking: Figure | None = None) -> tuple[Step, ...]:
    """Return the steps that read these designs back, the route's own in the middle.

    The picking is shared and the calling is shared; between them sits the route's own marking.

    `marking` is drawn on the one step of the route that changes a molecule: the lysate ligation
    on one route, the index PCR on the other. The records it names belong to the run, which this
    module does not hold, so the caller chooses it.

    Examples
    --------
    >>> from liulab_synbio.dmx.method import validation
    >>> one = validation(ROUTE_INDEX_PCR, (Design("a", 4),), 0)
    >>> for step in validation_steps(one):
    ...     print(step.title)
    Array 1 design(s) and grow
    Pick 4 colonies of each design
    Sample the picked plates into index plates
    Amplify each well with its own pair
    Pool and sequence
    Call every well
    """
    route = (
        _ligation_steps(one, marking)
        if one.route is ROUTE_LIGATION
        else _index_pcr_steps(one, marking)
    )
    return (_array_step(one), _pick_step(one), *route, _call_step(one))


def _array_step(one: Validation) -> Step:
    """Spot every design the floor reads back, which is where the read-back set is chosen."""
    floor = (
        "every design, which a floor of zero does"
        if one.floor == 0
        else f"every design in {one.floor} fragment(s) or more"
    )
    return Step(
        f"Array {len(one.designs)} design(s) and grow",
        key="array-designs",
        instructions=(
            "Spot each design from its archive plate as its own spot on a 25 cm BioAssay plate.",
            f"Grow overnight at 37 °C on {selected_on(one.selection)}.",
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
    sizes = ", ".join(f"{len(plate.labels)}" for plate in one.picked)
    return Step(
        f"Pick {one.colonies} colonies of each design",
        key="pick-colonies",
        instructions=(
            f"Pick {one.colonies} colonies a design into {CULTURE_UL:g} µL low-salt LB with "
            f"{selected_on(one.selection)}, with the {PICKER}.",
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


def _ligation_steps(one: Validation, marking: Figure | None = None) -> tuple[Step, ...]:
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
            key="compress-plates",
            instructions=(
                "Invert the picked plates for 30 minutes so the cells gather at the meniscus.",
                f"Move {LYSATE_UL:g} µL of each well into the {COMPRESSED_WELLS}-well plate.",
            ),
            transfers=moves,
            expected=(f"{one.wells} wells of lysate, {PLATES_COMPRESSED} picked plates to one.",),
        ),
        Step(
            "Barcode each well in lysate",
            key="barcode-wells",
            figures=() if marking is None else (marking,),
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


def _index_pcr_steps(one: Validation, marking: Figure | None = None) -> tuple[Step, ...]:
    """Sample a quarter at a time into index plates, amplify on the pair each well's address names."""
    at = 0
    moves: list[Transfer] = []
    for plate in one.picked:
        many = -(-len(plate.labels) // INDEX_WELLS)
        moves += sampling(plate, one.index[at : at + many])
        at += many
    return (
        Step(
            "Sample the picked plates into index plates",
            key="sample-into-index-plates",
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
            key="index-pcr",
            figures=() if marking is None else (marking,),
            instructions=(
                f"Add {INDEX_MIX_UL:g} µL of the master mix below to each well, which already "
                f"holds its {SAMPLE_UL:g} µL of culture.",
                f"Add {INDEX_PRIMER_UL:g} µL of the pair its address names from the prepared "
                f"primer plate, for {INDEX_PCR_UL:g} µL a well.",
                "Seal the plate, spin it down, and run the program below.",
            ),
            tables=(index_pcr_reaction(one.wells),),
            programs=(index_pcr_program(),),
            expected=(
                "One barcoded amplicon a well. Each well of a plate takes the forward mark its "
                "own place in the plate names, and every well of one plate takes that plate's "
                "own reverse mark.",
                "One band a well, the design plus about 100 bp: the two marks, and the stretch "
                "between the primer sites and the reading frame.",
            ),
            notes=(
                f"{INDEX_WELLS} forward marks and {INDEX_WELLS} reverse reach "
                f"{ROUTE_INDEX_PCR.capacity:,} wells, so the pairs already held cover far more than this "
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
            holes=(INDEX_MARKS,),
        ),
        _sequencing_step(one, "Pool each index plate on its own and clean the pool up."),
    )


def _sequencing_step(one: Validation, pooling_instruction: str) -> Step:
    """Pool the marked wells and sequence them, which both routes end their own stretch on."""
    return Step(
        "Pool and sequence",
        key="pool-and-sequence",
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
    tolerated = (
        "."
        if one.route.tolerable_reads is None
        else f"; the {one.route.name} route tolerates {one.route.tolerable_reads} reads and "
        "warns between the two."
    )
    return Step(
        "Call every well",
        key="call-wells",
        instructions=(
            "Demultiplex the reads by address, checking each address rather than trusting a file.",
            "Call a consensus a well, then compare it with that well's design base for base.",
            "Reformat, compacting out the wells that failed.",
        ),
        expected=(
            f"A well read more than {one.route.wanted_reads} times is called{tolerated}",
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

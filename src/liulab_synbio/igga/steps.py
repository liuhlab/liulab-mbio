"""The library-assembly bench protocol: its reactions, its programs, and the step order.

The protocol covers every round as one experiment. A round is two digests, a clean-up, a
ligation, a second clean-up, an electroporation, a growth and a prep, and the rounds run in the
order the scheme fills its positions.

Every number comes from `liulab_synbio.igga.bench` and `liulab_synbio.igga.coverage`, each
sourced in ``docs/research/protein-library-assembly.md``, or from the design the plan computed.
What the method leaves unpublished -- the ligase's units, the buffer's strength -- is one line
saying so rather than a number invented here. The pulse is the cells' own, read from their
catalogue number and printed with the manual it came from.

The shared step builders in `liulab_mbio.bench.steps` are shaped for a heat-shock
transformation, which this method does not run. What it does share is the reaction table, the
PCR that pulls a block out of the oligo pool, and the smaller pieces of a protocol.
"""

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import KW_ONLY, dataclass

from liulab_mbio import checks as judged
from liulab_mbio.barcodes import deletion_ambiguity
from liulab_mbio.bench.amounts import REFERENCES as AMOUNT_REFERENCES
from liulab_mbio.bench.amounts import Amount, to_nanograms
from liulab_mbio.bench.gels import REFERENCES as GEL_REFERENCES
from liulab_mbio.bench.gels import choose_ladder
from liulab_mbio.bench.materials import Electroporation, electroporation
from liulab_mbio.bench.pcr import REFERENCES as PCR_REFERENCES
from liulab_mbio.bench.pcr import pcr_program, pcr_reaction
from liulab_mbio.bench.phenotype import selection_marker
from liulab_mbio.bench.plates import plate, seat
from liulab_mbio.bench.prices import Item, PriceRecord
from liulab_mbio.bench.prices import bill as priced
from liulab_mbio.bench.reactions import reaction_table
from liulab_mbio.bench.steps import badges, card, catalogued, enzyme_material, listed
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.primers.polymerase import Q5, Polymerase, melting_temperature
from liulab_mbio.protocol.model import (
    Bill,
    Citation,
    Component,
    Gel,
    Incubation,
    Lane,
    Material,
    Plate,
    Protocol,
    ReactionTable,
    Reference,
    Source,
    Stage,
    Step,
    ThermocyclerProgram,
    Timer,
    Troubleshooting,
    citing,
    number,
)
from liulab_mbio.sequence import Segment, SequenceRecord
from liulab_mbio.sites import find_sites
from liulab_synbio import dmx
from liulab_synbio.igga import stages
from liulab_synbio.igga.bench import (
    DIGEST_CELSIUS,
    DIGEST_SECONDS,
    DIGEST_VOLUME_UL,
    ENZYME_UL,
    GROWTH_CELSIUS,
    LIGATION_SECONDS,
    LIGATION_VOLUME_UL,
    MOLAR_RATIO,
    OUTGROWTH_SECONDS,
    RECOVERY_SECONDS,
    SPRI_AFTER_DIGEST,
    SPRI_AFTER_LIGATION,
    TRANSFORMATION_NG,
    digest_amount,
    pool_floor_ng_ul,
)
from liulab_synbio.igga.bench import REFERENCES as BENCH_REFERENCES
from liulab_synbio.igga.cargo import Batch, PoolPlan
from liulab_synbio.igga.coverage import REFERENCES as COVERAGE_REFERENCES
from liulab_synbio.igga.coverage import (
    REPRESENTATION_MARKS,
    RepresentationMarks,
    RoundCoverage,
    absent_probability,
    colonies_for_completeness,
    reads_for_representation,
)
from liulab_synbio.igga.method import SYNTHESIS_ENZYME, Scheme
from liulab_synbio.igga.parts import Part
from liulab_synbio.igga.reads import ReadPair, ReadPairs
from liulab_synbio.igga.rounds import Round
from liulab_synbio.igga.standard import PartList, Standard
from liulab_synbio.igga.vector import Destination, Working, released_cargo

#: The buffer both digests run in, and the ligase and buffer the ligation runs in. The method
#: names all three and publishes neither the ligase's units nor either buffer's strength.
CUTSMART = "CutSmart Buffer"
LIGASE = "T7 DNA Ligase"
LIGASE_BUFFER = "StickTogether DNA Ligase Buffer"

#: What each clean-up is done with. The method gives the two volume ratios and not the product.
SPRI_BEADS = "SPRI paramagnetic beads"

#: The electrocompetent strain the method names, and who sells it.
STRAIN = "Endura ElectroCompetent Cells"
STRAIN_SUPPLIER = "Lucigen"
STRAIN_CATALOG = "60242-2"


def _pulse() -> Electroporation:
    """Return the strain's own program, which belongs to the cells and not to the step.

    Raises
    ------
    LookupError
        If the package stops shipping a program for these cells, rather than printing none.
    """
    found = electroporation(STRAIN_CATALOG)
    if found is None:
        raise LookupError(f"no electroporation program ships for {STRAIN_CATALOG}")
    return found


#: What amplifies the pool. The package's high-fidelity default, named here so the two PCRs
#: and the material agree about which buffer the annealing temperatures were computed in.
POOL_POLYMERASE: Polymerase = Q5
POOL_POLYMERASE_PRODUCT = "Q5 High-Fidelity DNA Polymerase (M0491)"

#: Twist's own cycle counts for amplifying an oligo pool: each band's longest oligo in nt, then
#: the fewest and the most cycles it allows. FRM-001034 REV 8 p. 2 and DOC-4060 REV 1.0 give the
#: same three. ``docs/research/oligo-pool-pcr-cycles.md`` section 2.
POOL_CYCLE_BANDS: tuple[tuple[int, int, int], ...] = ((100, 6, 10), (150, 10, 12), (350, 12, 14))

#: The document PCR1's cycle count is cited to, and where in it the count stands.
POOL_CYCLE_SOURCE_KEY = "FRM-001034"
POOL_CYCLE_CITATION = Citation(POOL_CYCLE_SOURCE_KEY, "p. 2, cycle chart")
POOL_SOURCES: dict[str, Source] = {
    POOL_CYCLE_SOURCE_KEY: Source(
        "Twist Bioscience, Amplifying Twist Oligo Pools",
        edition="REV 8",
        url="https://www.twistbioscience.com/content/dam/twistbioscience/resources/2026-01/"
        "FRM-001034-AmplifyingOligoPools-REV8%20singles.pdf",
        read_as="plain curl",
        date="2026-10-07",
    )
}

#: How many wells the plate PCR2 runs in holds. One batch is one plate of PCR2, which is what
#: fixes `liulab_synbio.igga.method.ORTHOGONAL_SPLIT`'s 96 inner primers.
PCR2_WELLS = 96
PCR2_PLATE = "PCR2 plate"

#: The hardware a round needs, which no reagent table covers.
EQUIPMENT: tuple[str, ...] = (
    f"Incubator or heat block at {DIGEST_CELSIUS:g} °C",
    "Magnetic rack for the bead clean-ups",
    "Electroporator and cuvettes",
    f"Shaking incubator at {GROWTH_CELSIUS:g} °C",
    "Spectrophotometer or fluorometer",
    "Long-read sequencer, for the two reads of the finished library",
)


@dataclass(frozen=True, slots=True)
class RoundBench:
    """What one round takes at the bench, computed from that round's own lengths.

    Parameters
    ----------
    number
        Which round it is, counting from one.
    position
        The position it fills, which names its part list.
    destination_digest, donor_digest
        What goes into each of the two digests: the library built so far, and the part list.
    ligation
        What the ligation takes, the opened destination first and the released part list second.
    transformation
        The most of the purified ligation one electroporation takes.
    coverage
        What this round has to cover, and what the colonies asked for leave out.
    """

    number: int
    position: str
    _: KW_ONLY
    destination_digest: Amount
    donor_digest: Amount
    ligation: tuple[Amount, Amount]
    transformation: Amount
    coverage: RoundCoverage


def choppers(scheme: Scheme) -> tuple[tuple[Enzyme, ...], tuple[Enzyme, ...]]:
    """Return the blunt enzymes each digest carries, the internal digest's first.

    A chopper belongs to the digest whose discarded piece it cuts: the internal stuffer core the
    internal digest excises, or the external stuffers the external digest leaves behind. One that
    reads a site in both is named by both, and the method has already refused one that reads a
    site in neither.
    """
    core = SequenceRecord(scheme.internal_stuffer_core)
    flanks = [
        SequenceRecord(bases) for bases in (scheme.external_stuffer_5, scheme.external_stuffer_3)
    ]
    inside = tuple(one for one in scheme.blunt if find_sites(core, one))
    outside = tuple(one for one in scheme.blunt if any(find_sites(flank, one) for flank in flanks))
    return inside, outside


def digest_reaction(
    dna: Amount,
    enzymes: Sequence[Enzyme],
    *,
    volume_ul: float = DIGEST_VOLUME_UL,
    reactions: int = 1,
) -> ReactionTable:
    """Return one digest: the DNA, each enzyme in turn, and the buffer and water that fill it.

    The buffer is one line with the water because the method names `CUTSMART` and not the
    strength it is supplied at, so its own volume is the supplier's to set.

    Raises
    ------
    ValueError
        If the DNA and the enzymes do not fit `volume_ul`.
    """
    return reaction_table(
        (dna,),
        tuple(Component(one.supplier_label, ENZYME_UL) for one in enzymes),
        volume_ul=volume_ul,
        title=f"Digest with {listed([one.name for one in enzymes])}",
        filler=f"{CUTSMART} and nuclease-free water",
        reactions=reactions,
    )


def digest_program(enzymes: Sequence[Enzyme]) -> ThermocyclerProgram:
    """Return the two hours a digest runs: the first enzyme alone, then the rest added to it."""
    first, rest = enzymes[0], tuple(enzymes[1:])
    stages = [Stage((Incubation(first.name, DIGEST_CELSIUS, DIGEST_SECONDS),))]
    if rest:
        added = listed([one.name for one in rest])
        stages.append(
            Stage((Incubation(f"{first.name} and {added}", DIGEST_CELSIUS, DIGEST_SECONDS),))
        )
    return ThermocyclerProgram(
        tuple(stages), title=f"Digest with {listed([one.name for one in enzymes])}"
    )


def ligation_reaction(
    amounts: tuple[Amount, Amount],
    *,
    volume_ul: float = LIGATION_VOLUME_UL,
    reactions: int = 1,
) -> ReactionTable:
    """Return one round's ligation: the two digest fragments, then the ligase, buffer and water.

    The last line carries all three together: the method publishes neither the ligase's units nor
    its volume, so the split between them is the supplier's to set and is not invented here.

    Raises
    ------
    ValueError
        If the two fragments do not fit `volume_ul`.
    """
    return reaction_table(
        amounts,
        volume_ul=volume_ul,
        title="Ligation",
        filler=f"{LIGASE} in {LIGASE_BUFFER}, and nuclease-free water",
        reactions=reactions,
    )


def growth_program() -> ThermocyclerProgram:
    """Return the recovery and the outgrowth, both at `GROWTH_CELSIUS` and not at 37 °C.

    The outgrowth carries the shorter of the two times the method gives; the step says the range.
    """
    return ThermocyclerProgram(
        (
            Stage((Incubation("Recovery, shaking", GROWTH_CELSIUS, RECOVERY_SECONDS),)),
            Stage((Incubation("Outgrowth", GROWTH_CELSIUS, OUTGROWTH_SECONDS[0]),)),
        ),
        title=f"Recovery and outgrowth at {GROWTH_CELSIUS:g} °C",
    )


def protocol(
    *,
    scheme: Scheme,
    positions: Sequence[str],
    barcode_length: int,
    vector: SequenceRecord,
    destination: Destination,
    part_lists: Sequence[PartList],
    standard: Standard,
    parts: Sequence[Part],
    rounds: Sequence[Round],
    bench: Sequence[RoundBench],
    constructs: int,
    checks: Sequence[judged.Check],
    host: str,
    sheet: str,
    barcodes: str,
    validation: dmx.Validation | None = None,
    prices: PriceRecord | None = None,
    pool: PoolPlan | None = None,
    pool_sheet: str = "",
    primer_sheet: str = "",
    working: Working | None = None,
    reads: ReadPairs | None = None,
    marks: RepresentationMarks = REPRESENTATION_MARKS,
    linkage_fidelity: float | None = None,
    block_vectors: Sequence[tuple[str, str]] = (),
) -> Protocol:
    """Return the bench protocol for one planned library, ready to render.

    Each argument is the `liulab_synbio.igga.plan.LibraryPlan` field or property of that name;
    `sheet` and `barcodes` are what the plan calls the two files a step points at. The steps run
    in the order someone does them: order the blocks, read back the designs the project asks for,
    pool each part list, then every round in turn, and finally read the finished library twice,
    for linkage and for representation.

    `validation` is `None` for a project that states no fragment-count floor, and the protocol
    then carries no validation at all: the library stays polyclonal, which is the default.

    The bill is always there, because its quantities come from the design. Money comes only from
    `prices`, and every row it does not price carries a hole.

    `working` is the vector the finished library is moved into, and `None` for a project naming
    none. Without one the final assembly's steps still run, because it is a stage of the method,
    and what the vector would have fixed is a hole instead.

    `pool` is the oligo pool the blocks are built from, and `None` for a project that writes
    none. With one, the blocks are not ordered: they are amplified out of the pool and
    assembled, so the first steps and the bill say that instead. `pool_sheet` and `primer_sheet`
    are what the plan calls the two files those steps point at.

    `block_vectors` is what a block closes into, one a position: what each is called and the file
    the plan wrote it to. A part enters on its own position's entry overhang, so one destination
    serves one position. Given none, the steps name the destination the first round opens and
    say nothing about the rest.
    """
    inside, outside = choppers(scheme)
    selection = stages.selection_for(vector)
    one = Protocol(
        f"Library assembly: {len(part_lists)} part lists into {vector.name or 'the vector'}",
        summary=(
            f"Join {len(parts)} synthesised parts into {constructs} distinct constructs in "
            f"{len(rounds)} rounds. Each round opens the library with {scheme.internal.name}, "
            f"releases one part list with {scheme.external.name}, ligates the two, and "
            "transforms, grows and preps the result for the round after it."
        ),
        overview=_overview(
            scheme,
            positions,
            barcode_length,
            standard,
            rounds,
            bench,
            constructs,
            part_lists,
            host,
            validation,
            pool,
        ),
        highlights=_highlights(
            scheme,
            positions,
            barcode_length,
            destination,
            standard,
            rounds,
            bench,
            constructs,
            barcodes,
        ),
        checks=badges(checks),
        materials=(
            *_materials(
                scheme,
                positions,
                part_lists,
                vector,
                sheet,
                pool,
                pool_sheet,
                primer_sheet,
                working,
                block_vectors,
            ),
            *(dmx.validation_materials(validation) if validation else ()),
        ),
        # Deduplicated: the pool route and the read-back route each ask for a thermocycler.
        equipment=tuple(
            dict.fromkeys(
                (
                    *EQUIPMENT,
                    *(POOL_EQUIPMENT if pool else ()),
                    *(dmx.validation_equipment(validation) if validation else ()),
                )
            )
        ),
        plates=(
            *((_pcr2_plate(pool),) if pool else ()),
            *(_validation_plates(validation) if validation else ()),
        ),
        steps=_steps(
            scheme,
            positions,
            barcode_length,
            parts,
            part_lists,
            rounds,
            bench,
            inside,
            outside,
            constructs,
            sheet,
            barcodes,
            validation,
            pool,
            pool_sheet,
            primer_sheet,
            working,
            selection,
            reads,
            marks,
            linkage_fidelity,
            block_vectors,
        ),
        references=(
            *_references(scheme),
            *(POOL_REFERENCES if pool else ()),
            *(dmx.REFERENCES if validation else ()),
        ),
        sources=_sources(prices, validation, pool),
        holes=stages.holes_for(vector),
        bill=_consumed(scheme, parts, rounds, inside, outside, prices, pool),
    )
    return citing(one)


def _sources(
    prices: PriceRecord | None, validation: dmx.Validation | None, pool: PoolPlan | None
) -> dict[str, Source]:
    """Return every document this run could cite; `citing` drops the ones it did not."""
    found = dict(stages.SOURCES)
    if validation:
        found |= dmx.SOURCES
    if pool:
        found |= POOL_SOURCES
    if prices:
        found[PRICES_SOURCE] = prices.source
    return found


def _overview(
    scheme: Scheme,
    positions: Sequence[str],
    barcode_length: int,
    standard: Standard,
    rounds: Sequence[Round],
    bench: Sequence[RoundBench],
    constructs: int,
    part_lists: Sequence[PartList],
    host: str,
    validation: dmx.Validation | None,
    pool: PoolPlan | None,
) -> dict[str, str]:
    """Return the facts to check before starting, each short enough to be a card."""
    last = bench[-1]
    product = rounds[-1].product
    blocks = sum(len(one) for one in part_lists)
    sizes = ", ".join(
        f"{position} {len(parts)}" for position, parts in zip(positions, part_lists, strict=True)
    )
    return {
        "Method": card(scheme.name, f"{len(positions)} positions"),
        "Part lists": card(sizes, f"{len(part_lists)} lists"),
        "Parts": (
            f"{blocks} blocks, assembled from {pool.pool.count} oligos"
            if pool
            else f"{blocks} synthesised blocks"
        ),
        "Constructs": f"{constructs:,} distinct",
        "Rounds": f"{len(rounds)}, one a part list",
        "Opened with": scheme.internal.supplier_label,
        "Released with": scheme.external.supplier_label,
        "Blunt enzymes": listed([one.name for one in scheme.blunt]),
        "Entry overhangs": card(listed(standard.entry_overhangs), "one a position"),
        "Cloning scar": standard.scar_overhang,
        "Barcode": (
            f"{barcode_length} bp, block "
            f"{scheme.barcode_block_length(barcode_length, len(positions))} bp"
        ),
        "Codon usage": host,
        "Product": card(f"{product.name}, {len(product)} bp", f"{len(product)} bp"),
        "Completeness": card(
            f"P {last.coverage.completeness:g}, {last.coverage.colonies:,} colonies at the end",
            f"{last.coverage.coverage:.0f}x",
        ),
        "Amino acids changed": f"{standard.cost} over {len(standard.changes)} part end(s)",
        "Designs read back": (
            card(
                f"{len(validation.designs)}, every one"
                if validation.floor == 0
                else f"{len(validation.designs)}, from {validation.floor} fragment(s)",
                f"route {validation.route.name}",
            )
            if validation
            else "none; the library stays polyclonal"
        ),
    }


def _highlights(
    scheme: Scheme,
    positions: Sequence[str],
    barcode_length: int,
    destination: Destination,
    standard: Standard,
    rounds: Sequence[Round],
    bench: Sequence[RoundBench],
    constructs: int,
    barcodes: str,
) -> tuple[str, ...]:
    """Return what the facts mean, a sentence each."""
    last = bench[-1]
    said = [
        f"One round appends one part list to every member of the library at once, so "
        f"{len(rounds)} rounds make {constructs} constructs out of "
        f"{sum(one.coverage.part_list_size for one in bench)} synthesised parts.",
        f"The product keeps {scheme.internal.name}'s sites, which is what lets the next round "
        f"open it, and loses {scheme.external.name}'s with the external stuffers.",
    ]
    if destination.edit is None:
        said.append(
            f"Your vector already carried an internal stuffer, so every part enters position "
            f"{positions[0]} on {standard.entry_overhangs[0]}, which is the overhang "
            "that stuffer spells: DNA that exists cannot be re-chosen."
        )
    else:
        said.append(
            f"Your vector carried no internal stuffer, so one was put in at "
            f"{destination.stuffer.start} carrying {standard.entry_overhangs[0]}, the overhang "
            "this design chose. Read the edit before you order the vector."
        )
    if standard.cost:
        said.append(
            f"The overhang standard changes {standard.cost} amino acid(s) across "
            f"{len(standard.changes)} part end(s); changes.tsv has wild type beside synthesised, "
            "and that is what you are about to pay for."
        )
    else:
        said.append(
            "The overhang standard changes no amino acid: every part list already spells its "
            "junctions."
        )
    said.append(
        f"The last round needs {last.coverage.colonies:,} colonies for a "
        f"{last.coverage.completeness:g} chance that none of its {last.coverage.products:,} "
        f"products is missing, which is {last.coverage.coverage:.0f}x its products. A round "
        "short of that floor loses members no later round can put back."
    )
    said.append(
        f"The finished barcode block is "
        f"{scheme.barcode_block_length(barcode_length, len(positions))} bp and reads in the "
        f"reverse of the order the rounds ran; {barcodes} says which barcode names which part."
    )
    return tuple(said)


def _destination_note(block_vectors: Sequence[tuple[str, str]]) -> str:
    """Say what the destination is for, and that a later position has one of its own."""
    opened = "opened by the internal digest in round 1"
    if len(block_vectors) < 2:
        return opened
    return (
        f"{opened}, and to take each block's cargo. A part enters on its own position's entry "
        f"overhang, so there is one a position: {listed(_vector_names(block_vectors))}"
    )


def _enzyme_material(enzyme: Enzyme, note: str) -> Material:
    """Return one enzyme as a material, carrying what every digest takes of it."""
    return enzyme_material(enzyme, amount=f"{ENZYME_UL:g} µL per digest", note=note)


def _selection_note(record: SequenceRecord, which: str) -> str:
    """Return what this vector's transformants are selected on, read off its own marker.

    A record annotating no marker this method can name a drug for leaves the plate to the reader
    rather than naming one, and `stages.ROUND_SELECTION` then stands as the protocol's hole.
    """
    plate = stages.selection_for(record)
    marker = selection_marker(record)
    if not plate or marker is None:
        return f"the {which} vector's own antibiotic, which this plan does not name"
    return f"LB with {plate}, the {which} vector's own marker ({marker.name})"


def _materials(
    scheme: Scheme,
    positions: Sequence[str],
    part_lists: Sequence[PartList],
    vector: SequenceRecord,
    sheet: str,
    pool: PoolPlan | None = None,
    pool_sheet: str = "",
    primer_sheet: str = "",
    working: Working | None = None,
    block_vectors: Sequence[tuple[str, str]] = (),
) -> tuple[Material, ...]:
    """Every reagent and consumable the protocol asks for.

    Where there is a pool the blocks are not bought, so the part lists say what they are built
    from and the pool, its primers, its polymerase and its assembly enzyme are the materials
    bought instead.
    """
    inside, outside = choppers(scheme)
    came_from = (
        f"assembled from the oligo pool; {pool_sheet} says which oligos"
        if pool
        else f"synthesised blocks, ordered from {sheet}"
    )
    made = [
        Material(
            f"{position} part list",
            storage="-20 °C",
            amount=f"{len(one)} members, pooled",
            note=came_from,
        )
        for position, one in zip(positions, part_lists, strict=True)
    ]
    if pool:
        made += list(_pool_materials(pool, pool_sheet, primer_sheet))
    made.append(
        Material(
            f"{vector.name or 'destination'} vector",
            storage="-20 °C",
            amount="" if len(block_vectors) < 2 else f"{len(block_vectors)}, one a position",
            note=_destination_note(block_vectors),
        )
    )
    if working is not None:
        made.append(
            Material(
                f"{working.record.name or 'Working'} vector",
                storage="-20 °C",
                note=f"the library's final home; {working.enzyme.name} releases its ccdB "
                "cassette to admit the cargo",
            )
        )
        made.append(
            _enzyme_material(working.enzyme, "admits the finished cargo to the working vector")
        )
    made.append(_enzyme_material(scheme.internal, "opens the library, excising its stuffer"))
    made.append(_enzyme_material(scheme.external, "releases a part from its own block"))
    for one in scheme.blunt:
        jobs = [
            what
            for what, named in (("the excised stuffer", inside), ("the donor backbone", outside))
            if one in named
        ]
        made.append(_enzyme_material(one, f"cuts {listed(jobs)}, so it cannot ligate back"))
    made.append(Material(CUTSMART, storage="-20 °C", note="both digests run in it"))
    # Both carry their own rules, so neither can appear in a protocol that does not show them.
    made += list(stages.method_materials())
    made.append(
        Material(
            SPRI_BEADS,
            amount=f"{SPRI_AFTER_DIGEST:g} volumes after a digest, "
            f"{SPRI_AFTER_LIGATION:g} after a ligation",
            note="the two ratios differ; both elute in water",
        )
    )
    pulse = _pulse()
    made.append(
        Material(
            STRAIN,
            supplier=STRAIN_SUPPLIER,
            catalog=STRAIN_CATALOG,
            storage="-80 °C",
            amount=f"one aliquot per round, {pulse.cells_ul:g} µL a pulse",
            note=f"pulsed at {pulse.volts:g} V, {pulse.ohms:g} Ω and {pulse.microfarads:g} µF "
            f"in a {pulse.cuvette_mm:g} mm cuvette",
            citation=pulse.citation,
        )
    )
    made.append(Material("Recovery medium", amount="one outgrowth per round"))
    made.append(Material("Selective broth and plates", note=_selection_note(vector, "destination")))
    if working is not None:
        made.append(
            Material(
                "Selective plates for the final transfer",
                note=_selection_note(working.record, "working"),
            )
        )
    made.append(Material("Plasmid prep kit", amount="one prep per round"))
    made.append(Material("Electroporation cuvettes", amount=f"{len(part_lists)}, one per round"))
    return tuple(made)


#: What a price record prices the synthesis order and the bill's own source by. Neither has a
#: catalogue number, so each is keyed by the name the protocol's own row already carries. The
#: pool's own key is the project's, through `liulab_synbio.igga.cargo.design_pool`.
BLOCKS_KEY = "synthesised blocks"
POOL_PRIMER_KEY = "pool-primers"
PRICES_SOURCE = "prices"


def _consumed(
    scheme: Scheme,
    parts: Sequence[Part],
    rounds: Sequence[Round],
    inside: Sequence[Enzyme],
    outside: Sequence[Enzyme],
    record: PriceRecord | None,
    pool: PoolPlan | None,
) -> Bill:
    """Return what this build buys, in the quantities the design computes.

    Only what the design fixes is billed. The buffer, the beads and the medium scale with volumes
    the method leaves to the supplier, so a row for one would be a quantity nobody computed.

    **The blocks are bought once or not at all.** A project with a pool buys oligos and the
    primers that amplify them, and assembles its blocks from those; one without buys the blocks
    themselves. Billing both would charge the same DNA twice.
    """
    digests = [(scheme.internal, len(rounds)), (scheme.external, len(rounds))]
    digests += [(one, len(rounds) * ((one in inside) + (one in outside))) for one in scheme.blunt]
    items = [
        *(
            (
                Item(
                    "Synthesised blocks",
                    len(parts),
                    unit="blocks",
                    key=BLOCKS_KEY,
                    quantities={
                        "count": len(parts),
                        "length_nt": max(one.length for one in parts),
                    },
                ),
            )
            if pool is None
            else (pool.item, _primer_item(pool))
        ),
        *(
            Item(
                one.commercial_name or one.name,
                ENZYME_UL * uses,
                unit="µL",
                key=one.catalog_number or "",
                quantities={"volume_ul": ENZYME_UL * uses},
            )
            for one, uses in digests
        ),
        Item(
            STRAIN,
            len(rounds),
            unit="aliquots",
            key=STRAIN_CATALOG,
            quantities={"count": len(rounds)},
        ),
        Item(
            "Electroporation cuvettes",
            len(rounds),
            unit="cuvettes",
            key="cuvettes",
            quantities={"count": len(rounds)},
        ),
        Item(
            "Plasmid preps",
            len(rounds),
            unit="preps",
            key="plasmid prep",
            quantities={"count": len(rounds)},
        ),
    ]
    return priced(items, record, source_key=PRICES_SOURCE)


def _references(scheme: Scheme) -> tuple[Reference, ...]:
    """Where the numbers come from, and where the scheme itself came from."""
    from liulab_mbio.cloning.goldengate.bench import REFERENCES as GOLDEN_GATE_REFERENCES

    items = [
        *BENCH_REFERENCES,
        *COVERAGE_REFERENCES,
        *AMOUNT_REFERENCES,
        *READOUT_REFERENCES,
        *POOL_RESUSPENSION_REFERENCES,
        *GOLDEN_GATE_REFERENCES,
    ]
    if scheme.source:
        items.append(Reference(f"The method this build was planned by: {scheme.source}"))
    return tuple(items)


#: What a synthesised block is resuspended in, and the floor both vendors publish for it. Each
#: gives 10 ng/µL as a least, not a value: at it a 1,000 ng pool fills 100 µL, more than twice
#: what the round's digest leaves, so the pool's own floor is computed and these set the buffer.
POOL_RESUSPENSION_REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Twist Bioscience, How should Multiplexed Gene Fragments be resuspended? — nuclease-free "
        "TE pH 8.0 or 10 mM Tris-HCl pH 8.0, at least 10 ng/µL for the stock dilution",
        url="https://www.twistbioscience.com/faq/multiplexed-gene-fragments/how-should-multiplexed-gene-fragments-be-resuspended",
    ),
    Reference(
        "Integrated DNA Technologies, gBlocks Gene Fragments resuspension — spin down, add IDTE "
        "or molecular-grade water to 10 ng/µL, vortex, 50 °C for 15-20 min, then verify",
        url="https://www.idtdna.com/page/?p=890",
    ),
)


#: What the two readout cautions of the confirming step are measured by.
READOUT_REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Van Nieuwerburgh, F. et al. (2011) Quantitative bias in Illumina TruSeq and a novel "
        "post amplification barcoding strategy for multiplexed DNA and small RNA deep "
        "sequencing. PLoS ONE 6, e26969, for a barcode 3 bp from the insert giving up to "
        "100-fold differences in read counts, where one introduced during the PCR 34 bp away "
        "gave R2 = 0.9977",
        url="https://doi.org/10.1371/journal.pone.0026969",
    ),
    Reference(
        "Alon, S. et al. (2011) Barcoding bias in high-throughput multiplex sequencing of "
        "miRNA. Genome Res. 21, 1506-1511, for ligated barcodes spreading read counts about "
        "twofold where the same barcodes introduced during the PCR did not",
        url="https://doi.org/10.1101/gr.121715.111",
    ),
)


def _steps(
    scheme: Scheme,
    positions: Sequence[str],
    barcode_length: int,
    parts: Sequence[Part],
    part_lists: Sequence[PartList],
    rounds: Sequence[Round],
    bench: Sequence[RoundBench],
    inside: Sequence[Enzyme],
    outside: Sequence[Enzyme],
    constructs: int,
    sheet: str,
    barcodes: str,
    validation: dmx.Validation | None,
    pool: PoolPlan | None = None,
    pool_sheet: str = "",
    primer_sheet: str = "",
    working: Working | None = None,
    selection: str = "",
    reads: ReadPairs | None = None,
    marks: RepresentationMarks = REPRESENTATION_MARKS,
    linkage_fidelity: float | None = None,
    block_vectors: Sequence[tuple[str, str]] = (),
) -> tuple[Step, ...]:
    """Return every step in the order it happens, the rounds one after another.

    A pool replaces the order step with the four it takes to get the same blocks: order the
    pool, pull each batch out of it, pull each block out of its batch, and clone that block's
    cargo into its own position's destination.
    """
    made = (
        list(
            _pool_steps(
                pool,
                parts,
                sheet,
                pool_sheet,
                primer_sheet,
                scheme,
                block_vectors or ((rounds[0].destination.name or "the destination", ""),),
                bench[0].ligation[0],
                selection,
            )
        )
        if pool
        else [_order_step(parts, sheet)]
    )
    if validation:
        made += dmx.validation_steps(validation)
    made.append(_pool_step(bench, parts, pool))
    for one, row in zip(rounds, bench, strict=True):
        made.extend(_round_steps(scheme, one, row, inside, outside, len(rounds), selection))
    made.append(
        _linkage_step(
            scheme,
            positions,
            barcode_length,
            rounds,
            parts,
            barcodes,
            None if reads is None else reads.linkage,
            linkage_fidelity,
        )
    )
    made.append(
        _representation_step(
            scheme,
            positions,
            barcode_length,
            rounds,
            constructs,
            barcodes,
            None if reads is None else reads.representation,
            marks,
        )
    )
    made += _final_steps(
        scheme,
        rounds,
        constructs,
        bench[-1].coverage.completeness,
        barcodes,
        working,
        None if reads is None else reads.final_representation,
        marks,
    )
    return tuple(made)


def _order_step(parts: Sequence[Part], sheet: str) -> Step:
    """Order every block, which is what the whole design comes down to."""
    lengths = [part.length for part in parts]
    return Step(
        "Order the synthesised parts",
        instructions=(
            f"Order every row of {sheet} as a double-stranded synthesised block.",
            "Ask for sequence-verified material: a block carries the scheme's sites at fixed "
            "offsets, and a base out of place stops it being cut where the design says.",
        ),
        expected=(
            f"{len(parts)} blocks, {min(lengths)} to {max(lengths)} bp.",
            "Each block reads: 5' external stuffer, coding bases, internal stuffer, barcode, "
            "3' external stuffer.",
        ),
        notes=(
            f"{sheet} carries each block's barcode on its own row, so the sheet you order from "
            "is also what decodes the sequencing afterwards.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The vendor cannot synthesise a block",
                "It is usually a repeat or a GC run in the coding bases. Plan again with another "
                "codon usage table; the stuffers and the overhangs are fixed by the scheme.",
            ),
        ),
    )


#: The hardware the pool route needs on top of a round's, which no reagent table covers.
POOL_EQUIPMENT: tuple[str, ...] = (
    "Thermocycler taking a 96-well plate",
    "Gel tank and a transilluminator",
)

#: Where PCR1's cycle count is read. Two independently revised Twist documents give the same
#: three length bands, and the second's appendix answers what more cycles cost.
POOL_CYCLE_REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Twist Bioscience, Amplifying Twist Oligo Pools, FRM-001034 REV 8, p. 2, for the cycle "
        "count banded by the pool's length"
    ),
    Reference(
        "Twist Bioscience, Twist Oligo Pools Amplification Protocol, DOC-4060 REV 1.0, for the "
        "same three bands, and for the FAQ answering that more cycles give worse uniformity"
    ),
)

#: What the two PCRs and their gel cite, beside the round's own references.
POOL_REFERENCES: tuple[Reference, ...] = (
    *PCR_REFERENCES,
    *GEL_REFERENCES,
    *POOL_CYCLE_REFERENCES,
)


def _validation_plates(one: dmx.Validation) -> tuple[Plate, ...]:
    """Return the plates the read-back fills, so every well a transfer names has one."""
    return (*one.picked, *(one.index if one.route is dmx.ROUTE_B else ()))


def _pcr2_plate(pool: PoolPlan) -> Plate:
    """Return the plate PCR2 runs in, one block a well, drawn as the first batch fills it.

    The wells are labelled rather than seated: a block is no material, oligo, vessel or plate,
    so its name describes the well instead of naming an occupant the protocol declares.
    """
    batches = pool.batches
    return plate(
        PCR2_PLATE,
        PCR2_WELLS,
        holds="one block a well, in the order its batch allotted the inner primers",
        labels=seat(batches[0].blocks, PCR2_WELLS) if batches else {},
        note=(
            f"{_counted(len(batches), 'plate')}, one a batch; "
            f"{', '.join(f'batch {one.number} holds {len(one.blocks)}' for one in batches)}"
        ),
    )


def _pool_materials(pool: PoolPlan, pool_sheet: str, primer_sheet: str) -> tuple[Material, ...]:
    """Return what the pool route buys that a project ordering its blocks does not."""
    layout = pool.pool.layout
    roles = ", ".join(f"{count} {role}" for role, count in _primer_roles(pool).items())
    return (
        Material(
            "Oligo pool",
            storage="-20 °C",
            amount=f"{pool.pool.count} oligos, every one {layout.length} nt",
            note=f"ordered from {pool_sheet}, which says which block each oligo is a piece of",
        ),
        Material(
            "Pool amplification primers",
            storage="-20 °C",
            amount=f"{len(pool.pool.primers)} primers: {roles}",
            note=f"ordered from {primer_sheet}; a pair a batch, an inner primer a block",
        ),
        catalogued(
            POOL_POLYMERASE_PRODUCT,
            supplier="NEB",
            storage="-20 °C",
            amount="one reaction per batch, then one per block",
            note="the buffer the annealing temperatures below were computed in",
        ),
        catalogued(
            "NEBridge Golden Gate Assembly Kit, BsmBI-v2 (E1602)",
            supplier="NEB",
            storage="-20 °C",
            amount="one assembly per block",
            note=f"its mix carries the {SYNTHESIS_ENZYME} that cuts every oligo back to its own "
            "fragment, which this method reserves, and the T4 DNA Ligase that joins them",
        ),
    )


def _primer_roles(pool: PoolPlan) -> dict[str, int]:
    """How many primers of the pool take each role, in the order the set allots them."""
    counted: dict[str, int] = {}
    for one in pool.pool.primers:
        counted[one.role] = counted.get(one.role, 0) + 1
    return counted


def _primer_item(pool: PoolPlan) -> Item:
    """Return the amplification primers as one line of the bill, beside the pool's own."""
    primers = pool.pool.primers
    return Item(
        "Pool amplification primers",
        len(primers),
        unit="primers",
        key=POOL_PRIMER_KEY,
        quantities={"count": len(primers), "length": max(len(one) for one in primers)},
    )


def _annealing(
    pairs: Iterable[tuple[str, str]], sequences: Mapping[str, str]
) -> tuple[float, float]:
    """Return the lowest and the highest annealing temperature these primer pairs ask for, °C."""
    found = [
        POOL_POLYMERASE.annealing_temperature(
            melting_temperature(sequences[one], POOL_POLYMERASE),
            melting_temperature(sequences[other], POOL_POLYMERASE),
        )
        for one, other in pairs
    ]
    return min(found), max(found)


def _counted(number: int, word: str) -> str:
    """Say a count and its noun, the noun plural only where the count is not one."""
    return f"{number} {word}" if number == 1 else f"{number} {word}s"


def _one_band(name: str, length_bp: int, *, title: str) -> Gel:
    """Return the gel one lane of a PCR makes, which is the same lane in every tube of it."""
    return Gel(choose_ladder((length_bp,)), (Lane(name, (length_bp,)),), title=title)


def _band_note(low: float, high: float) -> str:
    """Say what annealing temperature a block of tubes runs at, and the spread it covers."""
    if low == high:
        return f"Every pair anneals at {low:g} °C, so one block of tubes takes them all."
    return (
        f"The pairs anneal between {low:g} and {high:g} °C. The program runs at {low:g}, the "
        "lowest of them, so one block of tubes takes them all."
    )


def _pool_steps(
    pool: PoolPlan,
    parts: Sequence[Part],
    sheet: str,
    pool_sheet: str,
    primer_sheet: str,
    scheme: Scheme,
    destinations: Sequence[tuple[str, str]],
    opened: Amount,
    selection: str,
) -> tuple[Step, ...]:
    """Return the steps that turn an oligo pool into the blocks a round's part list holds."""
    sequences = {one.name: one.sequence for one in pool.pool.primers}
    layout = pool.pool.layout
    batches = pool.batches
    first = _annealing(((one.forward, one.outer) for one in batches), sequences)
    second = _annealing(pool.inner_pairs, sequences)
    inner_length = layout.length - layout.primer_length
    return (
        _pool_order_step(pool, pool_sheet, primer_sheet),
        _pcr1_step(pool, batches, first, layout.length),
        _pcr2_step(pool, batches, parts, second, inner_length),
        _assembly_step(pool, parts, sheet, scheme, destinations, opened, selection),
    )


def _pool_order_step(pool: PoolPlan, pool_sheet: str, primer_sheet: str) -> Step:
    """Order the pool and its primers, which is what the whole design comes down to."""
    layout = pool.pool.layout
    roles = ", ".join(f"{count} {role}" for role, count in _primer_roles(pool).items())
    over = pool.over_floor
    spent = (
        f"{pool.floor} oligos is what the length budget allows for; this design spends "
        f"{pool.pool.count}, one more on {listed(list(over))}, whose forced cuts spell no legal "
        "overhang."
        if over
        else f"{pool.pool.count} oligos is the fewest the length budget allows for."
    )
    return Step(
        "Order the oligo pool and the primers that amplify it",
        instructions=(
            f"Order every row of {pool_sheet} as one synthesised oligo pool, {pool.pool.count} "
            f"members at {layout.length} nt.",
            f"Order every row of {primer_sheet} as an ordinary oligo: {roles}.",
        ),
        expected=(
            f"One pool, {pool.pool.count} oligos, every one {layout.length} nt and no length "
            "spread, which is what the padding is for.",
            f"{len(pool.pool.primers)} primers, each "
            f"{min(len(one) for one in pool.pool.primers)} nt.",
        ),
        notes=(
            f"{pool_sheet} names each oligo's block, which piece of it that is, and the primers "
            "that pull it out, so an oligo is traceable to its protein.",
            spent,
        ),
        troubleshooting=(
            Troubleshooting(
                "The vendor bins the pool by length",
                "Every oligo is padded to one length, so the spread is zero. A vendor asking "
                "for a different length wants the project's oligo length changed and the "
                "design run again.",
            ),
        ),
    )


def pool_cycles(length_nt: int) -> tuple[int, int]:
    """Return the fewest and the most cycles Twist allows a pool of this length.

    The count is banded by length, so it follows the project's own oligo length rather than
    sitting fixed in the source. A pool longer than the last band takes that band's count: the
    table stops where the product does.

    Examples
    --------
    >>> pool_cycles(350)
    (12, 14)
    >>> pool_cycles(60)
    (6, 10)
    >>> pool_cycles(400)
    (12, 14)
    """
    for longest, fewest, most in POOL_CYCLE_BANDS:
        if length_nt <= longest:
            return (fewest, most)
    _, fewest, most = POOL_CYCLE_BANDS[-1]
    return (fewest, most)


def _pcr1_step(
    pool: PoolPlan,
    batches: Sequence[Batch],
    annealing: tuple[float, float],
    length_bp: int,
) -> Step:
    """Pull one batch of blocks out of the whole pool, which is what PCR1 is for."""
    low, high = annealing
    fewest, most = pool_cycles(length_bp)
    pairs = "; ".join(f"batch {one.number}: {one.forward} with {one.outer}" for one in batches)
    return Step(
        f"PCR1: pull {_counted(len(batches), 'batch')} out of the pool",
        instructions=(
            f"Set up {_counted(len(batches), 'reaction')}, one a batch, with the pool as template.",
            f"Give each its own pair: {pairs}.",
            "Run the program below. On a real-time instrument, add an intercalating dye and stop "
            "before the curve plateaus; the printed count is what to run without one.",
        ),
        cautions=("Keep the polymerase on ice.",),
        tables=(
            pcr_reaction(POOL_POLYMERASE, reactions=len(batches), title="PCR1, one tube a batch"),
        ),
        programs=(
            pcr_program(
                POOL_POLYMERASE,
                annealing_temperature=low,
                amplicon_length=length_bp,
                cycles=fewest,
                cycles_citation=POOL_CYCLE_CITATION,
                title="PCR1",
            ),
        ),
        gels=(_one_band("PCR1 product", length_bp, title="PCR1, any batch"),),
        expected=(
            f"One band at {length_bp} bp in every batch, which is the whole oligo.",
            f"Batch sizes: "
            f"{', '.join(f'{one.number} holds {len(one.blocks)}' for one in batches)}.",
        ),
        notes=(
            _band_note(low, high),
            f"Twist's band for a {length_bp} nt pool is {fewest} to {most} cycles; Twist Oligo "
            "Pools Amplification Protocol DOC-4060 REV 1.0 gives the same three bands. Its FAQ "
            f"answers that more cycles give worse uniformity, so {fewest} is what prints.",
            "The outer primer is what makes a batch a batch: it is dropped at PCR2, so a block "
            "cannot be pulled out of a batch it does not sit in.",
        ),
        holes=(stages.PCR1_POLYMERASE,),
        troubleshooting=(
            Troubleshooting(
                "No band",
                "Drop the annealing temperature by 3 °C and check the pool went in.",
            ),
            Troubleshooting(
                "A hump after the peak on a Bioanalyzer or TapeStation trace",
                "Heteroduplexes, which is what over-amplification leaves. Run it again with "
                "fewer cycles. What too many cycles cost is dropout, chimeras and polymerase "
                "error; pooling the blocks equimolar later restores the evenness but none of "
                "those.",
            ),
        ),
    )


def _pcr2_step(
    pool: PoolPlan,
    batches: Sequence[Batch],
    parts: Sequence[Part],
    annealing: tuple[float, float],
    length_bp: int,
) -> Step:
    """Pull one block out of its batch, after which a block is named by the well it sits in.

    The program's cycle count is blank. Twist's band is for amplifying the pool as it arrives,
    and this reaction's template is PCR1's product, so the reader is given the rule that stops
    the reaction and `liulab_synbio.igga.stages.PCR2_CYCLES` stands where the number would be.
    """
    low, high = annealing
    return Step(
        f"PCR2: pull each of the {len(parts)} blocks out of its batch",
        instructions=(
            f"Set up one reaction a block, {len(parts)} in all, in {PCR2_PLATE}.",
            "Give each its batch's PCR1 product as template, that batch's forward primer, and "
            "the block's own inner primer.",
            "Run the program below. Nobody published a cycle count for this reaction, so run it "
            "on a real-time instrument with an intercalating dye and stop before the curve "
            "plateaus.",
        ),
        cautions=("Keep the polymerase on ice.",),
        tables=(
            pcr_reaction(POOL_POLYMERASE, reactions=len(parts), title="PCR2, one well a block"),
        ),
        programs=(
            pcr_program(
                POOL_POLYMERASE,
                annealing_temperature=low,
                amplicon_length=length_bp,
                # A blank count has nothing to cite; `stages.PCR2_CYCLES` stands where it would.
                cycles=None,
                cycles_citation=None,
                title="PCR2",
            ),
        ),
        gels=(_one_band("PCR2 product", length_bp, title="PCR2, any well"),),
        expected=(
            f"One band at {length_bp} bp in every well, the outer primer's "
            f"{pool.pool.layout.primer_length} nt shorter than PCR1's.",
            f"{len(parts)} wells filled across {_counted(len(batches), 'plate')}.",
        ),
        notes=(
            _band_note(low, high),
            "A block is named by its well from here on, not by anything in the tube.",
        ),
        holes=(stages.PCR2_CYCLES,),
        troubleshooting=(
            Troubleshooting(
                "A band at PCR1's length",
                "The outer primer region was not dropped. Check the inner primer went in, and "
                "that the template is PCR1's product and not the pool.",
            ),
            Troubleshooting(
                "A hump after the peak on a Bioanalyzer or TapeStation trace",
                "Heteroduplexes, which is what over-amplification leaves. Run it again with "
                "fewer cycles; the forward primer is the whole batch's, so over-cycling here "
                "also pulls a neighbour's product into the well.",
            ),
        ),
    )


def _vector_names(destinations: Sequence[tuple[str, str]]) -> list[str]:
    """Say each destination as the bench reads it: what it is called, and the file holding it."""
    return [f"{name} ({file})" if file else name for name, file in destinations]


def _assembly_step(
    pool: PoolPlan,
    parts: Sequence[Part],
    sheet: str,
    scheme: Scheme,
    destinations: Sequence[tuple[str, str]],
    opened: Amount,
    selection: str,
) -> Step:
    """Clone each block's cargo into its position's opened destination, one well a block.

    The destination supplies the stuffers the cargo is no longer synthesised with, and its own
    marker is what selects a well that closed. A part enters on its own position's entry
    overhang, so there is one destination a position and a block goes into its own. The reaction
    and its cycling are NEB's kit table for the most pieces any block takes, so one master mix
    covers the plate.
    """
    from liulab_mbio.cloning.goldengate.bench import (
        KIT,
        assembly_amounts,
        assembly_program,
        assembly_reaction,
    )

    enzyme = get_enzyme(SYNTHESIS_ENZYME)
    pieces = {split.pieces for split in pool.splits}
    most = max(pieces)
    piece_bp = pool.pool.layout.length - pool.pool.layout.primer_length
    cargo = [len(split.cargo.sequence) for split in pool.splits]
    amounts = assembly_amounts(
        (opened.name, opened.length_bp),
        tuple((f"PCR2 piece {number}", piece_bp) for number in range(1, most + 1)),
    )
    named = listed(_vector_names(destinations))
    plated = selection or "the destinations' own marker, which this plan does not name"
    return Step(
        f"Assemble each cargo into its position's destination, from its {min(pieces)} to "
        f"{most} pieces",
        instructions=(
            f"Open {named} with {scheme.internal.name} and "
            f"{listed([one.name for one in choppers(scheme)[0]])}, the digest a round opens the "
            "library with, and clean each up.",
            "Set up one assembly a well, holding that block's own PCR2 pieces, the opened "
            "destination of the position that block fills, and nothing from another well.",
            f"Transform, plate on {plated}, and pick one colony a block.",
        ),
        tables=(assembly_reaction(enzyme, amounts, system=KIT, reactions=len(parts)),),
        programs=(assembly_program(enzyme, fragments=most + 1, system=KIT),),
        expected=(
            f"{len(parts)} plasmids, one a block: its position's destination carrying that "
            f"block's cargo, {min(cargo):,} to {max(cargo):,} bp of it.",
            f"Each cargo reads: the overhang its part enters on, its coding bases, the internal "
            f"stuffer, its barcode, and the {scheme.cloning_scar} cloning scar. The external "
            f"stuffers {sheet} spells either side of it are the destination's own bases.",
        ),
        notes=(
            f"{SYNTHESIS_ENZYME} cuts each oligo back to its fragment, so the primer sites, the "
            "recognition sites and the padding all stay outside the cargo.",
            "The table is sized at the most pieces any block takes, so one master mix covers "
            "the plate; a well with fewer pieces fills fewer of its DNA rows.",
            f"A part enters on its own position's entry overhang, so there is one destination a "
            f"position. They are one vector differing in the {len(scheme.entry_overhang)} bases "
            "of that overhang; a block put in the wrong one cannot close.",
        ),
        troubleshooting=(
            Troubleshooting(
                "A block comes out short",
                "A piece was missing from the well. Check that well's PCR2 lane before "
                "assembling it again; the pieces of one block are not interchangeable.",
            ),
            Troubleshooting(
                "Colonies carrying the destination with no cargo",
                "The destination was not opened to completion. Run a little of its digest on a "
                "gel before setting the plate up again.",
            ),
        ),
    )


def _pool_step(bench: Sequence[RoundBench], parts: Sequence[Part], pool: PoolPlan | None) -> Step:
    """Pool each part list, which is what a round joins in one tube.

    Every number here is the round's own donor digest read backwards. That digest takes
    `DIGEST_NG` of the pool in the volume its two enzymes leave it, so the pool's total, its
    concentration floor and each member's share are fixed downstream rather than chosen.
    """
    floor = _stated_floor()
    left = DIGEST_VOLUME_UL - 2 * ENZYME_UL
    members = {position: _at(parts, position) for position in {row.position for row in bench}}
    first = (
        "Clean up each assembly and measure each concentration."
        if pool
        else "Spin each tube down, resuspend in TE pH 8.0 or 10 mM Tris-HCl pH 8.0, 50 °C for "
        "15-20 min, and measure each concentration."
    )
    return Step(
        "Pool each part list",
        instructions=(
            first,
            "Pool the members of each part list in equal picomoles, one tube a position: each "
            "member gives the position's total divided by its member count, which is unequal "
            "masses because the members differ in length.",
            f"Bring each pool to at least {floor:g} ng/µL. The round's {DIGEST_VOLUME_UL:g} µL "
            f"digest leaves {left:g} µL for the DNA after its two {ENZYME_UL:g} µL enzymes, and "
            "the whole total has to arrive in it.",
        ),
        expected=tuple(
            f"{row.position}: one tube, {len(members[row.position])} member(s), at least "
            f"{row.donor_digest.nanograms:,.0f} ng at {floor:g} ng/µL or above, "
            f"{number(row.donor_digest.pmol / len(members[row.position]))} pmol of each member."
            for row in bench
            if members[row.position]
        ),
        notes=(
            "Library coverage is counted on equally represented members, so an uneven pool loses "
            "members that no later round can put back.",
            "Equal picomoles are unequal masses: a short member weighs less than a long one for "
            "the same number of molecules, and weighing them equally would not pool them equally.",
            *_pool_masses(bench, members),
        ),
        troubleshooting=(
            Troubleshooting(
                "One member is much more dilute than the rest",
                "Pool to the lowest member rather than to the mean; a member short here is short "
                "in every round after it.",
            ),
            Troubleshooting(
                f"The pool is below {floor:g} ng/µL",
                f"Concentrate it, by SPRI at the round's own {SPRI_AFTER_DIGEST:g}x ratio, or "
                "scale the digest up so the same mass arrives in a larger volume.",
            ),
        ),
    )


def _at(parts: Sequence[Part], position: str) -> tuple[Part, ...]:
    """Return the parts filling one position, in the order they were designed."""
    return tuple(one for one in parts if one.position == position)


def _stated_floor() -> float:
    """Return the pool's floor as a page states it, rounded up so it is never below the real one."""
    return math.ceil(pool_floor_ng_ul() * 10) / 10


def _pool_masses(
    bench: Sequence[RoundBench], members: Mapping[str, Sequence[Part]]
) -> tuple[str, ...]:
    """Return what a position's equal picomoles weigh, shortest member to longest."""
    said = []
    for row in bench:
        each = members[row.position]
        if not each:
            continue
        pmol = row.donor_digest.pmol / len(each)
        shortest, longest = min(each, key=_len_of), max(each, key=_len_of)
        said.append(
            f"{row.position}: {number(pmol)} pmol a member is "
            f"{number(to_nanograms(pmol, shortest.length))} ng of its shortest at "
            f"{shortest.length:,} bp and {number(to_nanograms(pmol, longest.length))} ng of its "
            f"longest at {longest.length:,} bp."
        )
    return tuple(said)


def _len_of(one: Part) -> int:
    """How many bases the part is ordered as."""
    return one.length


def _round_steps(
    scheme: Scheme,
    one: Round,
    row: RoundBench,
    inside: Sequence[Enzyme],
    outside: Sequence[Enzyme],
    total: int,
    selection: str = "",
) -> list[Step]:
    """Return the eight steps of one round, in the order they happen."""
    number = row.number
    opened, released = row.ligation
    internal = (scheme.internal, *inside)
    external = (scheme.external, *outside)
    return [
        _open_step(scheme, one, row, internal, opened),
        _release_step(scheme, row, external, released),
        _digest_cleanup_step(row, opened, released),
        _ligation_step(row, opened, released),
        _ligation_cleanup_step(row),
        _electroporation_step(row),
        _growth_step(row, selection),
        _prep_step(scheme, one, row, number == total),
    ]


def _digest_instructions(enzymes: Sequence[Enzyme]) -> tuple[str, ...]:
    """Return the two hours a digest runs, the second enzyme going in after the first."""
    first, rest = enzymes[0], tuple(enzymes[1:])
    said = [f"Set the digest below up with {first.name} alone, and start the program."]
    if rest:
        said.append(
            f"After the first hour add {listed([one.name for one in rest])} and run the second."
        )
    return tuple(said)


def _open_step(
    scheme: Scheme, one: Round, row: RoundBench, enzymes: Sequence[Enzyme], opened: Amount
) -> Step:
    """Open the library built so far, which is what the round's part enters."""
    chopped = listed([enzyme.name for enzyme in enzymes[1:]]) or "nothing else"
    return Step(
        f"Round {row.number}: open {row.destination_digest.name} with {scheme.internal.name}",
        instructions=_digest_instructions(enzymes),
        tables=(digest_reaction(row.destination_digest, enzymes),),
        programs=(digest_program(enzymes),),
        expected=(
            f"{row.destination_digest.name} opens on {one.entry_overhang} and "
            f"{one.scar_overhang}, leaving {opened.length_bp} bp of backbone.",
            f"The {one.excised.length} bp internal stuffer comes out, and {chopped} cuts it so "
            "it cannot go back in.",
        ),
        notes=(
            f"{scheme.internal.name} is the enzyme the product keeps sites for. That is the "
            "design working, not a site left behind.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Mostly empty library later",
                "The destination was not opened to completion. Run a little of the digest on a "
                "gel before ligating the next round.",
            ),
        ),
    )


def _release_step(
    scheme: Scheme, row: RoundBench, enzymes: Sequence[Enzyme], released: Amount
) -> Step:
    """Cut the round's whole part list out of its blocks, in one tube."""
    chopped = listed([enzyme.name for enzyme in enzymes[1:]]) or "nothing else"
    return Step(
        f"Round {row.number}: release the {row.position} part list with {scheme.external.name}",
        instructions=_digest_instructions(enzymes),
        tables=(digest_reaction(row.donor_digest, enzymes),),
        programs=(digest_program(enzymes),),
        expected=(
            f"Every member of the {row.position} part list is released, about "
            f"{released.length_bp} bp, on the same two overhangs the destination now offers.",
            f"{chopped} cuts the two external stuffers left behind, so neither can ligate back.",
        ),
        notes=(
            "One tube takes the whole part list: every member carries the same stuffers and "
            "differs only in its coding bases and its barcode.",
            "The amounts are weighed at the pool's mean block length, "
            f"{row.donor_digest.length_bp} bp.",
        ),
        troubleshooting=(
            Troubleshooting(
                "A member is missing from the library later",
                "Its block was under-represented in the pool, or its own digest was incomplete. "
                "Check the pool before repeating the round.",
            ),
        ),
    )


def _digest_cleanup_step(row: RoundBench, opened: Amount, released: Amount) -> Step:
    """Take the enzymes and the shredded pieces away, and measure what is left."""
    return Step(
        f"Round {row.number}: clean both digests up",
        instructions=(
            f"Add {SPRI_AFTER_DIGEST:g} volumes of {SPRI_BEADS} to each digest and elute in water.",
            "Measure both concentrations; the ligation table asks for picomoles, not nanograms.",
        ),
        expected=(
            f"{opened.name}: {opened.pmol:g} pmol is {opened.nanograms:g} ng at "
            f"{opened.length_bp} bp.",
            f"{released.name}: {released.pmol:g} pmol is {released.nanograms:g} ng at "
            f"{released.length_bp} bp.",
        ),
        notes=(
            f"{SPRI_AFTER_DIGEST:g} volumes here and {SPRI_AFTER_LIGATION:g} after the ligation. "
            "The two ratios differ; they are not one number used twice.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Too dilute for the ligation",
                "Elute in a smaller volume, or scale the ligation up.",
            ),
        ),
    )


def _ligation_step(row: RoundBench, opened: Amount, released: Amount) -> Step:
    """Join the opened library and the released part list, at a molar ratio and not a mass one."""
    return Step(
        f"Round {row.number}: ligate the {row.position} part list into the library",
        instructions=(
            "Pipette the two DNAs into the tube first, then the ligase, its buffer and water.",
            f"Hold at room temperature for {LIGATION_SECONDS // 60} minutes.",
        ),
        tables=(ligation_reaction(row.ligation),),
        timers=(Timer("Ligation", LIGATION_SECONDS),),
        expected=(
            f"Nothing visible. {released.pmol:g} pmol of the released part list against "
            f"{opened.pmol:g} pmol of the opened library is the {MOLAR_RATIO:g}:1 molar ratio.",
        ),
        notes=(
            "Picomoles, not nanograms: the shorter fragment weighs less at the same ratio.",
            f"The method names {LIGASE} and {LIGASE_BUFFER} and publishes neither the units nor "
            "the volume, so the last line of the table is the supplier's own.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Few colonies after this round",
                "Check that both digests went to completion; an uncut end cannot ligate.",
            ),
        ),
    )


def _ligation_cleanup_step(row: RoundBench) -> Step:
    """Desalt the ligation, which is what stops the cuvette arcing."""
    return Step(
        f"Round {row.number}: clean the ligation up",
        instructions=(
            f"Add {SPRI_AFTER_LIGATION:g} volume of {SPRI_BEADS} and elute in water.",
            "Measure the concentration.",
        ),
        expected=(
            f"Enough for one electroporation: {row.transformation.nanograms:g} ng is "
            f"{row.transformation.pmol:g} pmol at {row.transformation.length_bp} bp.",
        ),
        notes=(
            f"One volume here, not the {SPRI_AFTER_DIGEST:g} the digests took. Eluting in water "
            "rather than a buffer is what keeps the salt out of the cuvette.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The eluate is still salty",
                "Repeat the clean-up; salt carried into the cuvette arcs the pulse.",
            ),
        ),
    )


def _electroporation_step(row: RoundBench) -> Step:
    """Get the whole ligation into cells, which is where the library's size is won or lost."""
    pulse = _pulse()
    low, high = pulse.time_constant_ms
    return Step(
        f"Round {row.number}: electroporate into {STRAIN}",
        instructions=(
            f"Thaw one aliquot of {STRAIN} on ice, {pulse.cells_ul:g} µL a pulse.",
            f"Add at most {TRANSFORMATION_NG:g} ng of the purified ligation and mix without "
            "making bubbles.",
            f"Pulse at {pulse.volts:g} V, {pulse.ohms:g} Ω and {pulse.microfarads:g} µF in a "
            f"{pulse.cuvette_mm:g} mm cuvette.",
        ),
        cautions=("Keep the cells and the cuvette on ice; a warm cuvette arcs.",),
        expected=(f"A pulse with no arc, and a time constant of {low:g} to {high:g} ms.",),
        notes=(
            "The method names the instrument and none of the settings, so the program is the "
            f"cells' own: it is keyed by {STRAIN_CATALOG} and changes when the cells do. The "
            "materials table says where it was read.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The cuvette arcs",
                "The DNA carries salt. Clean it up again, elute in water, and use a fresh "
                "aliquot of cells.",
            ),
        ),
    )


def _growth_step(row: RoundBench, selection: str = "") -> Step:
    """Recover and grow, both at 30 °C, and bound the round on net colonies.

    Two plates, both growing during the outgrowth: a measured dilution of the recovery, and the
    same cut destination carried through the ligation with no donor added. What grows on the
    control is parental destination that survived, so the round is judged on the difference.
    """
    coverage = row.coverage
    dilution, control = stages.titre_plates(row.number)
    return Step(
        f"Round {row.number}: recover and grow at {GROWTH_CELSIUS:g} °C",
        instructions=(
            "Add recovery medium straight away and shake for the first hour.",
            f"Plate a measured dilution of the recovery on {selection or 'selection'} as "
            f"{dilution.name}, and grow the rest in selective broth.",
            f"Plate the no-donor ligation from the same digest as {control.name}, at the same "
            "dilution.",
        ),
        programs=(growth_program(),),
        expected=(
            f"At least {coverage.colonies:,} net colonies, scaled up from the dilution: the "
            f"floor for the {coverage.completeness:g} chance the project asked for that none of "
            f"its {coverage.products:,} distinct products is missing, equally represented.",
            f"At that count the chance a named product is missing is "
            f"{number(coverage.absent_probability)}.",
            f"{control.name} should be near empty beside it; its colonies come off the count.",
        ),
        notes=(
            f"Both steps run at {GROWTH_CELSIUS:g} °C and not at 37 °C. That is a library "
            "precaution rather than an oversight.",
            f"The program carries the shorter outgrowth; {OUTGROWTH_SECONDS[0] // 3600} to "
            f"{OUTGROWTH_SECONDS[1] // 3600} hours is the range the method gives.",
            f"The {coverage.completeness:g} is the project's own, and follows from the "
            "representation the screen downstream asks for. No source sets it, and nothing here "
            f"defaults it; it works out at {coverage.coverage:.0f}x this round's products.",
            "The control measures the chain the design rests on — two cuts, a blunt chopper, a "
            "ligase that refuses blunt ends and the 2x clean-up — rather than assuming it.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Fewer net colonies than the count above",
                "The round has lost library members and no later round can put them back. "
                "Electroporate more of the ligation, or run the round again.",
            ),
            Troubleshooting(
                "The no-donor control is not empty",
                "The destination was not opened to completion, or the blunt chopper missed. "
                "What grows is the round before's library, one position short, and it reaches "
                "the linkage read as truncated members.",
            ),
        ),
    )


def _prep_step(scheme: Scheme, one: Round, row: RoundBench, last: bool) -> Step:
    """Prep the round's plasmid, which is either the next destination or the finished library."""
    where = (
        "This is the finished library."
        if last
        else f"This is the destination round {row.number + 1} opens."
    )
    return Step(
        f"Round {row.number}: prep the library",
        instructions=(
            "Harvest the whole culture rather than a single colony.",
            "Prep the plasmid and measure the concentration.",
        ),
        expected=(
            f"One pool of plasmid, {len(one.product)} bp per molecule. {where}",
            f"{scheme.internal.name} still cuts it and {scheme.external.name} no longer does.",
        ),
        notes=(
            "Harvesting the whole culture is what keeps the library a library; picking colonies "
            "here throws away everything not picked.",
        ),
        troubleshooting=(
            Troubleshooting(
                "A low yield",
                "The culture was grown too briefly, or the antibiotic was wrong for this "
                "backbone. Check the vector's own marker.",
            ),
        ),
    )


def _with(pair: ReadPair | None) -> str:
    """Name the designed pair and its amplicon, or say nothing where none was designed."""
    if pair is None:
        return ""
    return (
        f", on {pair.forward.name} {pair.forward.sequence} and {pair.reverse.name} "
        f"{pair.reverse.sequence}, which give a {pair.amplicon_length} bp amplicon"
    )


def _as_platform(pair: ReadPair | None) -> str:
    """Say what the amplicon asks of the sequencing, which follows from what it has to carry."""
    return "" if pair is None else f" as a {pair.platform}"


def _marks_sentence(constructs: int, marks: RepresentationMarks, what: str) -> str:
    """State the three marks the counts are judged on, and the depth they are judged at."""
    depth = reads_for_representation(constructs, marks.reads_per_member) if constructs else 0
    at = f", which is {depth:,} reads over {constructs:,} combinations" if depth else ""
    return (
        f"At least {marks.seen:.1%} {what} seen; a 90th/10th percentile skew ratio below "
        f"{marks.skew:g}, judged at {marks.reads_per_member} or more reads a member{at}."
    )


def _linkage_step(
    scheme: Scheme,
    positions: Sequence[str],
    barcode_length: int,
    rounds: Sequence[Round],
    parts: Sequence[Part],
    barcodes: str,
    pair: ReadPair | None = None,
    fidelity: float | None = None,
) -> Step:
    """Read the whole cargo back, which is what says the barcode block still names its parts.

    Read once. A barcode that names the wrong member is the one fault no later round repairs,
    so this is where the library is carried forward or a round is sent back.
    """
    final = rounds[-1]
    ambiguous = max(
        deletion_ambiguity([one.barcode for one in parts if one.index == index])
        for index in range(len(positions))
    )
    order = listed([one.position for one in reversed(rounds)])
    return Step(
        "Read linkage",
        instructions=(
            f"Amplify the whole cargo out of the finished library, from the vector before the "
            f"first {rounds[0].entry_overhang} to the vector past the final "
            f"{rounds[0].scar_overhang}, so one read carries a member's parts and its barcode "
            f"block together{_with(pair)}.",
            f"Sequence the amplicon{_as_platform(pair)}.",
            f"Decode each read's block against {barcodes}, then read the coding bases beside it "
            "against the member that block names.",
        ),
        expected=(
            "A table from barcode combination to cargo: each read decodes to one member of each "
            "part list, and the coding bases it carries are that member's.",
            f"The {scheme.barcode_block_length(barcode_length, len(positions))} bp block reads "
            f"{order}, each barcode separated from the last by the cloning scar "
            f"{scheme.cloning_scar}, at {final.block.start}-{final.block.end} of the "
            "representative construct.",
        ),
        notes=(
            "The block reads in the reverse of the order the rounds ran: each round inserted its "
            "barcode ahead of the ones already there.",
            f"{ambiguous:.1%} of the single-base deletions a barcode can carry leave a read "
            "another barcode of the same part list could leave, which no read can be assigned "
            "through.",
            (
                f"This build passes the linkage read at {fidelity:.1%} of reads carrying a "
                "barcode that still names its part."
                if fidelity is not None
                else "Takacsi-Nagy's Figures 1D and 1E read about 95% of their reads carrying a "
                "valid barcode at every position, and nearly 90% of the library correctly "
                "linked. That is what one source reached, not a mark this library is held to."
            ),
        ),
        troubleshooting=(
            Troubleshooting(
                "A read carries fewer barcodes than there are rounds",
                "A round did not go in. Check that round's prep against its own record before "
                "blaming the sequencing.",
            ),
            Troubleshooting(
                "A combination decodes to a member the coding bases are not",
                f"That round joined a block to the wrong barcode, so {barcodes} no longer names "
                "what the library holds. Rebuild that round rather than carrying the table "
                "forward; no later round repairs it.",
            ),
        ),
        holes=(stages.READ_PASS_MARK,),
    )


def _final_steps(
    scheme: Scheme,
    rounds: Sequence[Round],
    constructs: int,
    completeness: float,
    barcodes: str,
    working: Working | None,
    pair: ReadPair | None = None,
    marks: RepresentationMarks = REPRESENTATION_MARKS,
) -> list[Step]:
    """Return the five steps that move the finished library into the working vector.

    One tube, staged: the cargo is freed and the enzymes that freed it are killed before the
    working vector, the cargo enzyme and the ligase go in. The assembly alone is `H24`, and a
    build naming no working vector carries `H31` in place of the enzyme and its cycling.
    """
    product = rounds[-1].product
    span = released_cargo(product, scheme)
    freeing = (scheme.external, *_shredders(scheme, product, span))
    return [
        _pick_working_step(scheme, working),
        _release_step_final(scheme, product, span, freeing),
        _assemble_step(scheme, product, span, working),
        _final_growth_step(constructs, completeness, working),
        _final_representation_step(barcodes, working, pair, constructs, marks),
    ]


def _shredders(scheme: Scheme, product: SequenceRecord, span: Segment | None) -> tuple[Enzyme, ...]:
    """Return the blunt enzymes cutting the backbone the cargo leaves, read off the record."""
    if span is None:
        return ()
    return tuple(
        one
        for one in scheme.blunt
        if any(not product.covers(span, site.span) for site in find_sites(product, one))
    )


def _pick_working_step(scheme: Scheme, working: Working | None) -> Step:
    """Pick the vector the library moves into, which is what fixes the cargo enzyme."""
    entry, scar = scheme.entry_overhang, scheme.scar_overhang
    if working is None:
        return Step(
            "Pick the working vector",
            instructions=(
                "Take one tube of the working vector stock for the application this library is "
                "built for.",
                f"Confirm its own cargo enzyme opens it on {entry} and {scar}, giving up the "
                "ccdB cassette the library displaces.",
            ),
            expected=(
                "One opened backbone and the ccdB cassette beside it, and nothing else cut.",
            ),
            notes=(
                "The cargo enzyme is chosen against every molecule in this pot, so it has to be "
                "settled before the blocks are designed: the pipeline reads it off the working "
                "vector before it designs any block, so no block spells it.",
            ),
            holes=(stages.WORKING_VECTOR,),
        )
    record, cargo = working.record, working.enzyme
    stuffer = working.destination.stuffer
    return Step(
        f"Pick the working vector and confirm {cargo.name} opens it",
        instructions=(
            f"Take one tube of {record.name or 'the working vector'} stock.",
            f"Digest a little of it with {cargo.supplier_label} and run it on a gel.",
        ),
        expected=(
            f"{len(record)} bp opens on {entry} and {scar}, giving up its "
            f"{stuffer.end - stuffer.start} bp ccdB cassette.",
            f"Two bands and no more: {cargo.name} reads this vector nowhere else.",
        ),
        notes=(
            f"{working.cargo.check.detail}.",
            "ccdB is what makes the assembly self-selecting: a vector that took no cargo keeps "
            "the cassette and kills its host.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The stock will not open",
                "The cassette is not where the record says, or the enzyme has lost activity. "
                "Sequence the stock before assembling into it.",
            ),
        ),
    )


def _release_step_final(
    scheme: Scheme,
    product: SequenceRecord,
    span: Segment | None,
    enzymes: Sequence[Enzyme],
) -> Step:
    """Free the cargo from the backbone the rounds ran in, then kill what freed it."""
    if span is None:
        return Step(
            "Release the cargo from the library backbone",
            instructions=(
                f"Digest the finished library with {scheme.external.name} and the blunt enzyme "
                "that shreds the backbone it leaves.",
                "Heat-kill both. Nothing is purified: the working vector goes into this tube.",
            ),
            expected=(
                f"The whole cargo free, on {scheme.entry_overhang} and {scheme.scar_overhang}.",
            ),
            notes=(
                "The sites that free the cargo belong to the vector the rounds ran in, not to "
                "the cargo, so a backbone without them cannot release it.",
            ),
            holes=(stages.CARGO_RELEASE,),
        )
    named = listed([one.name for one in enzymes])
    return Step(
        f"Release the cargo with {named}",
        instructions=(
            f"Digest the finished library with {named} at {DIGEST_CELSIUS:g} °C.",
            "Heat-kill, then leave the tube alone: nothing is purified between the two stages.",
        ),
        tables=(
            digest_reaction(
                digest_amount((product.name or "the finished library", len(product))), enzymes
            ),
        ),
        programs=(
            ThermocyclerProgram(
                (Stage((Incubation(named, DIGEST_CELSIUS, DIGEST_SECONDS),)),),
                title=f"Digest with {named}",
            ),
            *_kill(enzymes),
        ),
        expected=(
            f"The whole {span.end - span.start} bp cargo free, on {scheme.entry_overhang} and "
            f"{scheme.scar_overhang}, out of {len(product)} bp of library.",
            f"{listed([one.name for one in enzymes[1:]]) or 'Nothing else'} cuts the backbone it "
            "came out of, so that backbone cannot close again.",
        ),
        notes=(
            "One tube, two stages. Killing the releasing enzymes before the working vector goes "
            "in is what keeps them off its backbone.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The library backbone survives the assembly",
                "The blunt enzyme missed, or the heat kill was short. Both show up as colonies "
                "carrying the round's own vector rather than the working one.",
            ),
        ),
    )


def _kill(enzymes: Sequence[Enzyme]) -> tuple[ThermocyclerProgram, ...]:
    """Return the suppliers' heat inactivations, enzymes agreeing on one sharing a program."""
    shared: dict[tuple[int, int], list[Enzyme]] = {}
    for one in enzymes:
        celsius, minutes = one.heat_inactivation_celsius, one.heat_inactivation_minutes
        if celsius is not None and minutes is not None:
            shared.setdefault((celsius, minutes), []).append(one)
    return tuple(
        ThermocyclerProgram(
            (Stage((Incubation("Heat inactivation", float(celsius), minutes * 60),)),),
            title=f"{listed([one.name for one in named])} heat inactivation",
        )
        for (celsius, minutes), named in shared.items()
    )


def _assemble_step(
    scheme: Scheme, product: SequenceRecord, span: Segment | None, working: Working | None
) -> Step:
    """Join the freed cargo to the opened working vector, in the tube the release left."""
    from liulab_mbio.cloning.goldengate.bench import assembly_program

    if working is None:
        return Step(
            "Assemble the cargo into the working vector",
            instructions=(
                "Add the working vector, its cargo enzyme and the ligase to the release tube, "
                "and run the enzyme's own Golden Gate cycling.",
            ),
            expected=(
                "One circular final vector a member, the ccdB cassette displaced by the cargo.",
            ),
            notes=(
                "This is the one reaction where the working vector meets material the rounds "
                "made; the rounds all finish first.",
            ),
            holes=(stages.WORKING_VECTOR, stages.FINAL_MASSES),
        )
    cargo = working.enzyme
    joined = (
        f"about {len(working.record) - _cassette_length(working) + span.end - span.start} bp, "
        if span is not None
        else ""
    )
    return Step(
        f"Assemble the cargo into {working.record.name or 'the working vector'} with {cargo.name}",
        instructions=(
            f"Add the working vector, {cargo.supplier_label} and {LIGASE} in "
            f"{LIGASE_BUFFER} to the release tube.",
            "Run the cycling below without purifying anything first.",
        ),
        programs=(assembly_program(cargo, fragments=2, library=True),),
        expected=(
            f"One circular final vector a member, {joined}joined on "
            f"{scheme.entry_overhang} and {scheme.scar_overhang}.",
            "The ccdB cassette is displaced, so a vector that took no cargo kills its host.",
        ),
        notes=(
            f"{len(product)} bp of library goes in and the cargo alone comes out: the backbone "
            "the rounds ran in is shredded and stays behind.",
            "The cycling is NEB's longer single-insert program, which it gives for library "
            "preparation rather than for cloning one gene.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Colonies that are still ccdB-positive",
                "The working vector was not opened to completion. Check the pick step's gel "
                "before repeating.",
            ),
        ),
        holes=(stages.FINAL_MASSES,),
    )


def _cassette_length(working: Working) -> int:
    """How many bases the ccdB cassette takes out of the working vector."""
    stuffer = working.destination.stuffer
    return stuffer.end - stuffer.start


def _final_growth_step(constructs: int, completeness: float, working: Working | None) -> Step:
    """Clean the assembly up and get all of it into cells, which is the library's last bottleneck."""
    pulse = _pulse()
    colonies = colonies_for_completeness(constructs, completeness)
    return Step(
        f"Clean the assembly up and electroporate into {STRAIN}",
        instructions=(
            f"Add {SPRI_BEADS} at {SPRI_AFTER_LIGATION:g}x the volume and elute in water, as "
            "every round did after its ligation.",
            f"Pulse at {pulse.volts:g} V, {pulse.ohms:g} Ω and {pulse.microfarads:g} µF in a "
            f"{pulse.cuvette_mm:g} mm cuvette.",
            f"Recover and grow at {GROWTH_CELSIUS:g} °C, as every round did.",
        ),
        programs=(growth_program(),),
        cautions=("Keep the cells and the cuvette on ice; a warm cuvette arcs.",),
        expected=(
            f"At least {colonies:,} net colonies: the floor for the {completeness:g} chance the "
            f"project asked for that none of its {constructs:,} distinct members is missing, "
            "equally represented.",
            f"At that count the chance a named member is missing is "
            f"{number(absent_probability(constructs, colonies))}.",
            "Near-empty plates from a no-cargo control beside it; what grows there is working "
            "vector that kept its ccdB cassette.",
        ),
        notes=(
            "This is a bottleneck like a round's, and the library can only lose members here. "
            "Electroporate all of the assembly rather than a measured part of it.",
            "The settings are the cells' own and are keyed by their catalogue number, not set "
            "by this method.",
        )
        + (
            ()
            if working is None
            else (
                f"ccdB does the selecting here, so the plates carry "
                f"{working.record.name or 'the working vector'}'s own marker.",
            )
        ),
        troubleshooting=(
            Troubleshooting(
                "Fewer net colonies than the count above",
                "The library has lost members in the transfer. Nothing downstream puts them "
                "back; repeat the assembly from more of the released cargo.",
            ),
        ),
    )


def _final_representation_step(
    barcodes: str,
    working: Working | None,
    pair: ReadPair | None = None,
    constructs: int = 0,
    marks: RepresentationMarks = REPRESENTATION_MARKS,
) -> Step:
    """Read the library again on the other side of the move, which is the only way to size the loss."""
    where = working.record.name if working is not None else "the working vector"
    return Step(
        "Read representation in the final vector",
        instructions=(
            f"Amplify across the barcode block again{_with(pair)}.",
            f"Sequence, decode each read against {barcodes}, and compare the counts with the "
            "read taken in the library backbone.",
        ),
        expected=(
            "The same combinations, at a similar evenness. A combination seen before the move "
            f"and not after it was lost in the transfer into {where or 'the working vector'}.",
            _marks_sentence(constructs, marks, "of what survived the move"),
        ),
        notes=(
            "Linkage is read once, in the library backbone; a barcode still names the same part "
            "after the move, because the move carries the whole cargo in one piece.",
            "The forward anchor is the same retained internal stuffer, which travels with the "
            "cargo; the reverse anchor moved with the vector, so this pair is designed against "
            "the final record rather than reused from the read before the move.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Markedly fewer combinations than before the move",
                "The transfer was the bottleneck, not the rounds. Compare the colony count with "
                "the step above before rebuilding anything.",
            ),
        ),
    )


def _representation_step(
    scheme: Scheme,
    positions: Sequence[str],
    barcode_length: int,
    rounds: Sequence[Round],
    constructs: int,
    barcodes: str,
    pair: ReadPair | None = None,
    marks: RepresentationMarks = REPRESENTATION_MARKS,
) -> Step:
    """Count which combinations the library holds and how evenly, over the barcode block alone.

    This is the read a bottleneck repeats, so it spans the block and nothing else: the forward
    anchor is the internal stuffer every member keeps, which is a method constant.
    """
    block = scheme.barcode_block_length(barcode_length, len(positions))
    return Step(
        "Read representation",
        instructions=(
            f"Amplify across the {block} bp barcode block alone, forward from the "
            f"{len(scheme.internal_stuffer)} bp internal stuffer every member keeps and back "
            f"from the vector past the final {rounds[0].scar_overhang}{_with(pair)}.",
            f"Sequence the amplicon{_as_platform(pair)}, decode each read against {barcodes}, "
            "and count the reads each barcode combination gets.",
        ),
        expected=(
            f"A count for each of up to {constructs:,} distinct combinations: how many of them "
            "are seen at all, and how evenly they are read.",
            "A combination with no reads is a member the library has lost.",
            _marks_sentence(constructs, marks, "of the combinations"),
        ),
        notes=(
            "Read representation again after every later bottleneck — the final assembly, and "
            "anything downstream that resamples the library. Linkage is read once; this one is "
            "read at each.",
            "The short amplicon is what makes repeating it cheap: linkage spans the whole cargo, "
            "this spans the block.",
            "Where you amplify the block to read it, carry any sample index on a primer "
            "rather than ligating it on, and keep the barcodes away from where a primer "
            "anneals: both cost more read counts than what a barcode spells does.",
        ),
        troubleshooting=(
            Troubleshooting(
                "A combination is missing",
                "It was lost at a round or at a bottleneck since, and no later step puts it "
                "back. Read the rounds' own titre plates before re-reading this one.",
            ),
            Troubleshooting(
                "The counts are heavily skewed",
                "Members differ in length and a bottleneck can favour the short ones. Imkeller's "
                "Table 2 prices the skew in screen coverage: a 90th/10th ratio of 2.5 wants "
                "200-fold, 5 wants 300-fold and 10 wants 400-fold, so a skewed library costs "
                "cells downstream rather than failing here.",
            ),
        ),
    )

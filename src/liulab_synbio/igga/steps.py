"""The library-assembly bench protocol: its reactions, its programs, and the step order.

The protocol covers every round as one experiment. A round is two digests, a clean-up, a
ligation, a second clean-up, an electroporation, a growth and a prep, and the rounds run in the
order the scheme fills its positions.

Every number comes from `liulab_synbio.igga.bench` and `liulab_synbio.igga.coverage`, each
sourced in ``docs/research/protein-library-assembly.md``, or from the design the plan computed.
What the method leaves unpublished -- the ligase's units, the buffer's strength, the
electroporation settings -- is one line saying so rather than a number invented here.

The shared step builders in `liulab_mbio.bench.steps` are shaped for a heat-shock
transformation, which this method does not run. What it does share is the reaction table, the
PCR that pulls a block out of the oligo pool, and the smaller pieces of a protocol.
"""

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import KW_ONLY, dataclass

from liulab_mbio import checks as judged
from liulab_mbio.barcodes import deletion_ambiguity
from liulab_mbio.bench.amounts import REFERENCES as AMOUNT_REFERENCES
from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.gels import REFERENCES as GEL_REFERENCES
from liulab_mbio.bench.gels import choose_ladder
from liulab_mbio.bench.pcr import REFERENCES as PCR_REFERENCES
from liulab_mbio.bench.pcr import pcr_program, pcr_reaction
from liulab_mbio.bench.plates import plate
from liulab_mbio.bench.prices import Item, PriceRecord
from liulab_mbio.bench.prices import bill as priced
from liulab_mbio.bench.reactions import reaction_table
from liulab_mbio.bench.steps import badges, card, catalogued, enzyme_material, listed
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.primers.polymerase import Q5, Polymerase, melting_temperature
from liulab_mbio.protocol.model import (
    Bill,
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
)
from liulab_mbio.sequence import SequenceRecord
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
)
from liulab_synbio.igga.bench import REFERENCES as BENCH_REFERENCES
from liulab_synbio.igga.cargo import Batch, PoolPlan
from liulab_synbio.igga.coverage import REFERENCES as COVERAGE_REFERENCES
from liulab_synbio.igga.coverage import RoundCoverage
from liulab_synbio.igga.method import SYNTHESIS_ENZYME, Scheme
from liulab_synbio.igga.parts import Part
from liulab_synbio.igga.rounds import Round
from liulab_synbio.igga.standard import PartList, Standard
from liulab_synbio.igga.vector import Destination

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

#: What amplifies the pool. The package's high-fidelity default, named here so the two PCRs
#: and the material agree about which buffer the annealing temperatures were computed in.
POOL_POLYMERASE: Polymerase = Q5
POOL_POLYMERASE_PRODUCT = "Q5 High-Fidelity DNA Polymerase (M0491)"

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

    `pool` is the oligo pool the blocks are built from, and `None` for a project that writes
    none. With one, the blocks are not ordered: they are amplified out of the pool and
    assembled, so the first steps and the bill say that instead. `pool_sheet` and `primer_sheet`
    are what the plan calls the two files those steps point at.
    """
    inside, outside = choppers(scheme)
    return Protocol(
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
                scheme, positions, part_lists, vector, sheet, pool, pool_sheet, primer_sheet
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
        plates=(_pcr2_plate(pool),) if pool else (),
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
        ),
        references=(
            *_references(scheme),
            *(POOL_REFERENCES if pool else ()),
            *(dmx.REFERENCES if validation else ()),
        ),
        sources=_sources(prices, validation),
        holes=stages.HOLES,
        bill=_consumed(scheme, parts, rounds, inside, outside, prices, pool),
    )


def _sources(prices: PriceRecord | None, validation: dmx.Validation | None) -> dict[str, Source]:
    """Return every document this run's citations resolve against, and no document it never cites."""
    found = dict(stages.SOURCES)
    if validation:
        found |= dmx.SOURCES
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
        "Coverage": f"{last.coverage.coverage:g}x, {last.coverage.colonies:,} colonies at the end",
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
        f"The last round needs {last.coverage.colonies:,} colonies for "
        f"{last.coverage.coverage:g}x coverage of its {last.coverage.products:,} products. A "
        "round short of that loses members no later round can put back."
    )
    said.append(
        f"The finished barcode block is "
        f"{scheme.barcode_block_length(barcode_length, len(positions))} bp and reads in the "
        f"reverse of the order the rounds ran; {barcodes} says which barcode names which part."
    )
    return tuple(said)


def _enzyme_material(enzyme: Enzyme, note: str) -> Material:
    """Return one enzyme as a material, carrying what every digest takes of it."""
    return enzyme_material(enzyme, amount=f"{ENZYME_UL:g} µL per digest", note=note)


def _materials(
    scheme: Scheme,
    positions: Sequence[str],
    part_lists: Sequence[PartList],
    vector: SequenceRecord,
    sheet: str,
    pool: PoolPlan | None = None,
    pool_sheet: str = "",
    primer_sheet: str = "",
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
            note="opened by the internal digest in round 1",
        )
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
    made.append(
        Material(
            STRAIN,
            supplier=STRAIN_SUPPLIER,
            catalog=STRAIN_CATALOG,
            storage="-80 °C",
            amount="one aliquot per round",
        )
    )
    made.append(Material("Recovery medium", amount="one outgrowth per round"))
    made.append(
        Material(
            "Selective broth and plates",
            note="the destination vector's own antibiotic, which this plan does not name",
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
    items = [*BENCH_REFERENCES, *COVERAGE_REFERENCES, *AMOUNT_REFERENCES, *READOUT_REFERENCES]
    if scheme.source:
        items.append(Reference(f"The method this build was planned by: {scheme.source}"))
    return tuple(items)


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
) -> tuple[Step, ...]:
    """Return every step in the order it happens, the rounds one after another.

    A pool replaces the order step with the four it takes to get the same blocks: order the
    pool, pull each batch out of it, pull each block out of its batch, and assemble the block.
    """
    made = (
        list(_pool_steps(pool, parts, sheet, pool_sheet, primer_sheet))
        if pool
        else [_order_step(parts, sheet)]
    )
    if validation:
        made += dmx.validation_steps(validation)
    made.append(_pool_step(positions, part_lists))
    for one, row in zip(rounds, bench, strict=True):
        made.extend(_round_steps(scheme, one, row, inside, outside, len(rounds)))
    made.append(_linkage_step(scheme, positions, barcode_length, rounds, parts, barcodes))
    made.append(
        _representation_step(scheme, positions, barcode_length, rounds, constructs, barcodes)
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

#: What the two PCRs and their gel cite, beside the round's own references.
POOL_REFERENCES: tuple[Reference, ...] = (*PCR_REFERENCES, *GEL_REFERENCES)


def _pcr2_plate(pool: PoolPlan) -> Plate:
    """Return the plate PCR2 runs in, one block a well.

    The seating is not declared. Which block sits in which well follows from the order its
    batch allotted the inner primers, and `liulab_mbio.bench.plates.seat` would write it, but
    `liulab_mbio.protocol.model.Protocol.audit` resolves every seated well against a declared
    material, oligo, vessel or plate, and a block is none of those. Issue 330 is open on that
    rule, so the plate says what it holds and the wells stay unnamed rather than being named
    wrongly.
    """
    batches = pool.batches
    return plate(
        PCR2_PLATE,
        PCR2_WELLS,
        holds="one block a well, in the order its batch allotted the inner primers",
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
        enzyme_material(
            get_enzyme(SYNTHESIS_ENZYME),
            amount="one assembly per block",
            note="cuts every oligo back to its own fragment, and is reserved for that",
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
        _assembly_step(pool, parts, sheet),
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


def _pcr1_step(
    pool: PoolPlan,
    batches: Sequence[Batch],
    annealing: tuple[float, float],
    length_bp: int,
) -> Step:
    """Pull one batch of blocks out of the whole pool, which is what PCR1 is for."""
    low, high = annealing
    pairs = "; ".join(f"batch {one.number}: {one.forward} with {one.outer}" for one in batches)
    return Step(
        f"PCR1: pull {_counted(len(batches), 'batch')} out of the pool",
        instructions=(
            f"Set up {_counted(len(batches), 'reaction')}, one a batch, with the pool as template.",
            f"Give each its own pair: {pairs}.",
            "Run the program below.",
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
            "The outer primer is what makes a batch a batch: it is dropped at PCR2, so a block "
            "cannot be pulled out of a batch it does not sit in.",
        ),
        troubleshooting=(
            Troubleshooting(
                "No band",
                "Drop the annealing temperature by 3 °C and check the pool went in.",
            ),
            Troubleshooting(
                "A smear rather than a band",
                "Too many cycles over a pool loses the evenness the coverage is counted on. "
                "Take the fewest cycles that give a visible band.",
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
    """Pull one block out of its batch, after which a block is named by the well it sits in."""
    low, high = annealing
    return Step(
        f"PCR2: pull each of the {len(parts)} blocks out of its batch",
        instructions=(
            f"Set up one reaction a block, {len(parts)} in all, in {PCR2_PLATE}.",
            "Give each its batch's PCR1 product as template, that batch's forward primer, and "
            "the block's own inner primer.",
            "Run the program below.",
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
        troubleshooting=(
            Troubleshooting(
                "A band at PCR1's length",
                "The outer primer region was not dropped. Check the inner primer went in, and "
                "that the template is PCR1's product and not the pool.",
            ),
        ),
    )


def _assembly_step(pool: PoolPlan, parts: Sequence[Part], sheet: str) -> Step:
    """Assemble each block out of its own pieces, which is the step the model cannot finish."""
    pieces = {split.pieces for split in pool.splits}
    lengths = [part.length for part in parts]
    return Step(
        f"Assemble each block from its {min(pieces)} to {max(pieces)} pieces",
        instructions=(
            f"Set up one {SYNTHESIS_ENZYME} assembly a well, holding that block's own PCR2 "
            "pieces and nothing from another well.",
        ),
        expected=(
            f"{len(parts)} blocks, {min(lengths)} to {max(lengths)} bp, as {sheet} spells them.",
            "Each block reads: 5' external stuffer, coding bases, internal stuffer, barcode, "
            "3' external stuffer.",
        ),
        notes=(
            f"{SYNTHESIS_ENZYME} cuts each oligo back to its fragment, so the primer sites, the "
            "recognition sites and the padding all stay outside the block.",
        ),
        holes=stages.POOL_HOLES,
        troubleshooting=(
            Troubleshooting(
                "A block comes out short",
                "A piece was missing from the well. Check that well's PCR2 lane before "
                "assembling it again; the pieces of one block are not interchangeable.",
            ),
        ),
    )


def _pool_step(positions: Sequence[str], part_lists: Sequence[PartList]) -> Step:
    """Pool each part list, which is what a round joins in one tube."""
    return Step(
        "Pool each part list",
        instructions=(
            "Resuspend every block and measure each concentration.",
            "Pool the members of each part list in equal molar amounts, one tube per position.",
        ),
        expected=tuple(
            f"{position}: one tube holding {len(one)} member(s)."
            for position, one in zip(positions, part_lists, strict=True)
        ),
        notes=(
            "Library coverage is counted on equally represented members, so an uneven pool loses "
            "members that no later round can put back.",
        ),
        troubleshooting=(
            Troubleshooting(
                "One member is much more dilute than the rest",
                "Pool to the lowest member rather than to the mean; a member short here is short "
                "in every round after it.",
            ),
        ),
    )


def _round_steps(
    scheme: Scheme,
    one: Round,
    row: RoundBench,
    inside: Sequence[Enzyme],
    outside: Sequence[Enzyme],
    total: int,
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
        _growth_step(row),
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
    return Step(
        f"Round {row.number}: electroporate into {STRAIN}",
        instructions=(
            f"Thaw one aliquot of {STRAIN} on ice.",
            f"Add at most {TRANSFORMATION_NG:g} ng of the purified ligation and mix without "
            "making bubbles.",
            "Pulse with the settings the cell supplier gives for your cuvette.",
        ),
        cautions=("Keep the cells and the cuvette on ice; a warm cuvette arcs.",),
        expected=(
            "A pulse with no arc, and a time constant in the range the cell supplier's manual "
            "gives.",
        ),
        notes=(
            "The method names the instrument and not the voltage, the capacitance, the "
            "resistance or the cuvette gap, so the settings are the cell supplier's.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The cuvette arcs",
                "The DNA carries salt. Clean it up again, elute in water, and use a fresh "
                "aliquot of cells.",
            ),
        ),
    )


def _growth_step(row: RoundBench) -> Step:
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
            f"Plate a measured dilution of the recovery on selection as {dilution.name}, and "
            "grow the rest in selective broth.",
            f"Plate the no-donor ligation from the same digest as {control.name}, at the same "
            "dilution.",
        ),
        programs=(growth_program(),),
        expected=(
            f"At least {coverage.colonies:,} net colonies, scaled up from the dilution: "
            f"{coverage.products:,} distinct products this round can make, times the "
            f"{coverage.coverage:g}x the project asked for, rounded up.",
            f"At that count the chance a named product is missing is "
            f"{coverage.absent_probability:.3g}.",
            f"{control.name} should be near empty beside it; its colonies come off the count.",
        ),
        notes=(
            f"Both steps run at {GROWTH_CELSIUS:g} °C and not at 37 °C. That is a library "
            "precaution rather than an oversight.",
            f"The program carries the shorter outgrowth; {OUTGROWTH_SECONDS[0] // 3600} to "
            f"{OUTGROWTH_SECONDS[1] // 3600} hours is the range the method gives.",
            f"{coverage.coverage:g}x is the project's own, and follows from the representation "
            "the screen downstream asks for. No source sets it, and nothing here defaults it.",
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


def _linkage_step(
    scheme: Scheme,
    positions: Sequence[str],
    barcode_length: int,
    rounds: Sequence[Round],
    parts: Sequence[Part],
    barcodes: str,
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
            "block together.",
            "Sequence the amplicon as a long read.",
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
            "This plan designs no sequencing primers.",
            f"{ambiguous:.1%} of the single-base deletions a barcode can carry leave a read "
            "another barcode of the same part list could leave, which no read can be assigned "
            "through.",
            "Takacsi-Nagy's Figures 1D and 1E read about 95% of their reads carrying a valid "
            "barcode at every position, and nearly 90% of the library correctly linked. That is "
            "what one source reached, not a mark this library is held to.",
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
        holes=(stages.READ_PRIMERS,),
    )


def _representation_step(
    scheme: Scheme,
    positions: Sequence[str],
    barcode_length: int,
    rounds: Sequence[Round],
    constructs: int,
    barcodes: str,
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
            f"from the vector past the final {rounds[0].scar_overhang}.",
            f"Sequence the amplicon, decode each read against {barcodes}, and count the reads "
            "each barcode combination gets.",
        ),
        expected=(
            f"A count for each of up to {constructs:,} distinct combinations: how many of them "
            "are seen at all, and how evenly they are read.",
            "A combination with no reads is a member the library has lost.",
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
                "Members differ in length and a bottleneck can favour the short ones. This read "
                "measures the skew; nothing here says how much of it is tolerable.",
            ),
        ),
        holes=(stages.READ_PASS_MARK,),
    )

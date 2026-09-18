"""The Gateway bench protocol: its own reaction steps, its own notes, and the step order.

The steps any bench shares are `liulab_mbio.bench.steps`. This module runs them around the LR
reaction and adds what only Gateway has to say. Every number is computed by
`liulab_mbio.cloning.gateway.plan` or is one `liulab_mbio.cloning.gateway.bench` cites.
"""

from collections.abc import Sequence

from liulab_mbio import checks as judged
from liulab_mbio.bench import REFERENCES as BENCH_REFERENCES
from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.phenotype import Phenotype
from liulab_mbio.bench.steps import badges, card, listed, phenotype_sentences, transform_step
from liulab_mbio.cloning.gateway.bench import (
    CELL_EFFICIENCY_CFU_UG,
    DESTINATION_NG,
    ENTRY_MIN_NG,
    ENTRY_NG,
    GATEWAY_TRANSFORMATION,
    LR_CELSIUS,
    LR_CLONASE,
    LR_CLONASE_CATALOG,
    LR_CLONASE_UL,
    LR_COLONIES,
    LR_LONG_BP,
    LR_LONG_SECONDS,
    LR_SECONDS,
    LR_VOLUME_UL,
    PROTEINASE_K_UG_UL,
    PROTEINASE_K_UL,
    REFERENCES,
    STOP_CELSIUS,
    STOP_SECONDS,
    SUPPLIER,
    TE_BUFFER,
    lr_reaction,
)
from liulab_mbio.cloning.gateway.recombination import Junction, Recombination
from liulab_mbio.protocol.model import Material, Protocol, Step, Timer, Troubleshooting
from liulab_mbio.sequence import SequenceRecord

#: The hardware a run needs, which no reagent table covers.
EQUIPMENT: tuple[str, ...] = (
    f"Water bath or heat block at {LR_CELSIUS:g} °C",
    f"Water bath or heat block at {STOP_CELSIUS:g} °C",
    f"Heat block or water bath at {GATEWAY_TRANSFORMATION.heat_shock_celsius:g} °C",
    f"Shaking incubator and a plate incubator at {GATEWAY_TRANSFORMATION.outgrowth_celsius:g} °C",
    "Microcentrifuge",
)


def protocol(
    *,
    entry: SequenceRecord,
    destination: SequenceRecord,
    recombination: Recombination,
    amounts: Sequence[Amount],
    phenotype: Phenotype,
    checks: Sequence[judged.Check],
    host: str,
) -> Protocol:
    """Return the bench protocol for one planned LR reaction, ready to render.

    Each argument is the `liulab_mbio.cloning.gateway.plan.Plan` field or property of that name.
    The steps run in the order someone does them: set the reaction up, run it, stop it with
    proteinase K, then transform and plate.
    """
    insert = recombination.moved
    return Protocol(
        f"Gateway LR: {entry.name} into {destination.name}",
        summary=(
            f"Recombine the {insert.length} bp segment between {entry.name}'s attL sites into "
            f"{destination.name}, in one LR reaction, and select the expression clone on the "
            "destination vector's own marker."
        ),
        overview=_overview(entry, destination, recombination, phenotype),
        highlights=_highlights(recombination, phenotype),
        checks=badges(checks),
        materials=_materials(entry=entry, destination=destination, host=host, phenotype=phenotype),
        equipment=EQUIPMENT,
        steps=_steps(
            recombination=recombination,
            amounts=amounts,
            phenotype=phenotype,
            entry=entry,
            host=host,
        ),
        references=(*REFERENCES, *BENCH_REFERENCES),
    )


def _overview(
    entry: SequenceRecord,
    destination: SequenceRecord,
    recombination: Recombination,
    phenotype: Phenotype,
) -> dict[str, str]:
    """Return the facts to check before starting, each short enough to be a card."""
    product = recombination.product
    facts = {
        "Entry clone": f"{entry.name}, {len(entry)} bp",
        "Destination vector": f"{destination.name}, {len(destination)} bp",
        "Reaction": f"{LR_CLONASE}, {LR_VOLUME_UL:g} µL at {LR_CELSIUS:g} °C",
        "Insert": f"{recombination.moved.length} bp between the att sites",
        "Junctions": card(listed([one.name for one in recombination.junctions]), "two att sites"),
        "Product": f"{product.name}, {len(product)} bp",
    }
    if phenotype.antibiotic:
        facts["Selection"] = phenotype.antibiotic
    elif phenotype.marker is not None:
        facts["Selection"] = f"{phenotype.marker.name} marker"
    return facts


def _highlights(recombination: Recombination, phenotype: Phenotype) -> tuple[str, ...]:
    """Return what the facts mean, a sentence each: what moves, what the scar is, what grows."""
    junctions = recombination.junctions
    return (
        f"One reaction moves the {recombination.moved.length} bp between the entry clone's att "
        f"sites into the {recombination.backbone.length} bp of destination vector outside its "
        "own, and the ccdB cassette leaves with the by-product.",
        f"The junction is not scarless: the clone gains a whole att site at each end of the "
        f"insert, {junctions[0].name} spelling {junctions[0].bases} and {junctions[1].name} "
        f"spelling {junctions[1].bases}. A fusion reads through both.",
        *phenotype_sentences(phenotype, [recombination.moved.name or "the insert"]),
    )


def _materials(
    *,
    entry: SequenceRecord,
    destination: SequenceRecord,
    host: str,
    phenotype: Phenotype,
) -> tuple[Material, ...]:
    """Every reagent and consumable the protocol asks for."""
    return (
        Material(
            f"{entry.name} entry clone",
            storage="-20 °C",
            amount=f"{ENTRY_NG:g} ng per reaction",
            note="supercoiled, which is the substrate the manual calls most efficient",
        ),
        Material(
            f"{destination.name} destination vector",
            storage="-20 °C",
            amount=f"{DESTINATION_NG:g} ng per reaction",
            note="grown in a ccdB-resistant strain, which is the only kind it grows in",
        ),
        Material(
            LR_CLONASE,
            supplier=SUPPLIER,
            catalog=LR_CLONASE_CATALOG,
            storage="-20 °C",
            amount=f"{LR_CLONASE_UL:g} µL per reaction",
            note="thaw on ice and return it to the freezer at once",
        ),
        Material(TE_BUFFER, amount=f"to {LR_VOLUME_UL - LR_CLONASE_UL:g} µL per reaction"),
        Material(
            "Proteinase K solution",
            supplier=SUPPLIER,
            storage="-20 °C",
            amount=f"{PROTEINASE_K_UL:g} µL per reaction",
            note=f"{PROTEINASE_K_UG_UL:g} µg/µL, supplied with the enzyme mix",
        ),
        Material(
            host,
            storage="-80 °C",
            amount=f"{GATEWAY_TRANSFORMATION.cells_ul:g} µL per transformation",
            note=f"{CELL_EFFICIENCY_CFU_UG:,} cfu/µg or better, and no F' episome",
        ),
        Material(
            "S.O.C. medium",
            amount=f"{GATEWAY_TRANSFORMATION.outgrowth_ul:g} µL per transformation",
        ),
        Material(_plate(phenotype), amount="one plate per transformation"),
    )


def _plate(phenotype: Phenotype) -> str:
    """Return what to pour the selection plates with."""
    antibiotic = phenotype.antibiotic or "the destination vector's own antibiotic"
    return f"LB agar plates with {antibiotic}"


def _steps(
    *,
    recombination: Recombination,
    amounts: Sequence[Amount],
    phenotype: Phenotype,
    entry: SequenceRecord,
    host: str,
) -> tuple[Step, ...]:
    """Return the steps in the order they happen."""
    return (
        _reaction_step(amounts),
        _incubation_step(recombination),
        _stop_step(),
        transform_step(
            host,
            phenotype,
            inserts=[recombination.moved.name or "the insert"],
            colonies=(
                f"More than {LR_COLONIES:,} colonies where the whole reaction is transformed "
                f"and plated, with cells at {CELL_EFFICIENCY_CFU_UG:,} cfu/µg or better. What "
                "fraction of them is correct is not published; Hartley 2000 counted 96 of 102."
            ),
            protocol=GATEWAY_TRANSFORMATION,
            notes=(
                "Unreacted destination vector and the by-product both keep the ccdB gene, "
                f"which kills {host}, so they do not grow. A strain carrying F' would supply "
                "ccdA and cancel that.",
                f"Small colonies beside large ones are usually unreacted {entry.name} "
                "co-transforming; restreak them on the entry clone's own antibiotic to tell.",
            ),
        ),
    )


def _reaction_step(amounts: Sequence[Amount]) -> Step:
    """Set the one-tube recombination up."""
    return Step(
        "Set up the LR reaction",
        instructions=(
            "Thaw the enzyme mix on ice and vortex it briefly twice.",
            "Pipette the two plasmids and the TE buffer into a tube at room temperature.",
            f"Add {LR_CLONASE_UL:g} µL of {LR_CLONASE}, mix well and spin down.",
        ),
        tables=(lr_reaction(amounts),),
        expected=(f"A {LR_VOLUME_UL:g} µL reaction holding both plasmids.",),
        notes=(
            f"The volumes above take each plasmid at the concentration the manual's own table "
            f"assumes; pipette what your prep needs for the nanograms and make the difference "
            f"up with {TE_BUFFER}.",
            f"Do not go over {ENTRY_NG:g} ng of entry clone: the manual reports colonies "
            f"carrying several molecules above it, and fewer colonies below {ENTRY_MIN_NG:g} ng.",
            "Supercoiled plasmids are the substrates the manual calls most efficient for LR.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The DNA does not fit the reaction volume",
                "Concentrate either plasmid, or scale the whole reaction up keeping the enzyme "
                "mix at its stated fraction of the volume.",
            ),
        ),
    )


def _incubation_step(recombination: Recombination) -> Step:
    """Run it, and say what the product carries when it is done."""
    return Step(
        "Run the LR reaction",
        instructions=(f"Incubate at {LR_CELSIUS:g} °C for {LR_SECONDS // 3600} hour.",),
        timers=(Timer("LR incubation", LR_SECONDS),),
        expected=(
            "Nothing visible: the reaction is a recombination, not a digest.",
            *(_junction_sentence(one) for one in recombination.junctions),
        ),
        notes=(
            f"A plasmid of {LR_LONG_BP:,} bp or more runs up to "
            f"{LR_LONG_SECONDS // 3600} hours instead; efficiency falls as the DNA gets longer.",
            "Junction positions are 0-based, on the product.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Few or no colonies later, though the transformation control worked",
                "Use an entry clone (attL) with a destination vector (attR): the LR reaction "
                "takes those and no others. Do not freeze and thaw the enzyme mix more than "
                "ten times.",
            ),
        ),
    )


def _junction_sentence(junction: Junction) -> str:
    """Say what one junction spells and which record gave which side of it."""
    return (
        f"The expression clone carries {junction.name} at {junction.start}, spelling "
        f"{junction.bases}: its first bases from {junction.before} and the rest from "
        f"{junction.after}."
    )


def _stop_step() -> Step:
    """End the reaction, which the manual requires before transforming."""
    return Step(
        "Stop the reaction with proteinase K",
        instructions=(
            f"Add {PROTEINASE_K_UL:g} µL of proteinase K at {PROTEINASE_K_UG_UL:g} µg/µL.",
            f"Incubate at {STOP_CELSIUS:g} °C for {STOP_SECONDS // 60} minutes.",
        ),
        timers=(Timer("Proteinase K", STOP_SECONDS),),
        expected=("A reaction that can go straight into competent cells.",),
        notes=(
            "The manual lists an untreated reaction as a cause of few or no colonies, so this "
            "step is not optional.",
        ),
    )

"""The Gateway bench protocol: its own reaction steps, its own notes, and the step order.

The steps any bench shares are `liulab_mbio.bench.steps`. This module runs them around each
recombination and adds what only Gateway has to say. A plan that runs BP reads as two staged
reactions: BP, its plate, the miniprep that purifies the entry clone, then LR and its own
plate, with the attB PCR and its cleanup in front of them where one was designed. Every number
is computed by `liulab_mbio.cloning.gateway.plan` or is one
`liulab_mbio.cloning.gateway.bench` cites.
"""

from collections.abc import Sequence

from liulab_mbio import checks as judged
from liulab_mbio.bench import REFERENCES as BENCH_REFERENCES
from liulab_mbio.bench.oligos import oligo_row
from liulab_mbio.bench.pcr import PRIMER_STOCK_UM
from liulab_mbio.bench.phenotype import Phenotype
from liulab_mbio.bench.steps import (
    badges,
    card,
    cleanup_step,
    listed,
    pcr_step,
    pcr_title,
    phenotype_sentences,
    transform_step,
)
from liulab_mbio.cloning.gateway.att import REGION_BP
from liulab_mbio.cloning.gateway.bench import (
    BP_CELSIUS,
    BP_CLONASE,
    BP_CLONASE_CATALOG,
    BP_CLONASE_UL,
    BP_COLONIES,
    BP_DONOR_MAX_NG,
    BP_LONG_BP,
    BP_LONG_SECONDS,
    BP_PMOL,
    BP_SECONDS,
    BP_SUBSTRATE_MIN_PMOL,
    BP_TOTAL_MAX_NG,
    BP_TRANSFORMATION,
    BP_VOLUME_UL,
    CELL_EFFICIENCY_CFU_UG,
    CHLORAMPHENICOL_UG_ML,
    DESTINATION_NG,
    ENTRY_MIN_NG,
    ENTRY_NG,
    LR_CELSIUS,
    LR_CLONASE,
    LR_CLONASE_CATALOG,
    LR_CLONASE_UL,
    LR_COLONIES,
    LR_LONG_BP,
    LR_LONG_SECONDS,
    LR_SECONDS,
    LR_TRANSFORMATION,
    LR_VOLUME_UL,
    ONE_TUBE_YIELD,
    PROPAGATION_HOST,
    PROTEINASE_K_UG_UL,
    PROTEINASE_K_UL,
    REFERENCES,
    STOP_CELSIUS,
    STOP_SECONDS,
    SUPPLIER,
    TE_BUFFER,
    bp_reaction,
    lr_reaction,
)
from liulab_mbio.cloning.gateway.design import SPACER, Amplicon
from liulab_mbio.cloning.gateway.recombination import Junction, PlannedReaction
from liulab_mbio.protocol.model import Material, Oligo, Protocol, Step, Timer, Troubleshooting
from liulab_mbio.sequence import SequenceRecord

#: The hardware a run needs, which no reagent table covers.
EQUIPMENT: tuple[str, ...] = (
    f"Water bath or heat block at {LR_CELSIUS:g} °C",
    f"Water bath or heat block at {STOP_CELSIUS:g} °C",
    f"Heat block or water bath at {LR_TRANSFORMATION.heat_shock_celsius:g} °C",
    f"Shaking incubator and a plate incubator at {LR_TRANSFORMATION.outgrowth_celsius:g} °C",
    "Microcentrifuge",
)

#: What a run amplifying its own insert needs on top of `EQUIPMENT`.
PCR_EQUIPMENT: tuple[str, ...] = (
    "Thermocycler with a heated lid",
    "Agarose gel rig and power supply",
    "Spectrophotometer or fluorometer",
)


def protocol(
    *,
    lr: PlannedReaction,
    bp: PlannedReaction | None,
    amplicon: Amplicon | None = None,
    checks: Sequence[judged.Check],
    host: str,
) -> Protocol:
    """Return the bench protocol for one planned Gateway experiment, ready to render.

    Each argument is the `liulab_mbio.cloning.gateway.plan.Plan` field or property of that name.
    The steps run in the order someone does them, one reaction at a time: set it up, run it,
    stop it with proteinase K, transform and plate. A planned BP reaction puts its own four
    steps and the miniprep that follows them in front of LR's, and an attB PCR puts its own
    two in front of those.
    """
    return Protocol(
        _title(lr, bp),
        summary=_summary(lr, bp),
        overview=_overview(lr, bp, amplicon),
        highlights=_highlights(lr, bp, amplicon),
        checks=badges(checks),
        materials=_materials(lr=lr, bp=bp, amplicon=amplicon, host=host),
        oligos=_oligos(amplicon),
        equipment=EQUIPMENT if amplicon is None else (*PCR_EQUIPMENT, *EQUIPMENT),
        steps=(
            *_pcr_steps(amplicon),
            *_bp_steps(bp, host=host, entry=_carrier(lr)),
            *_lr_steps(lr, host=host),
        ),
        references=(*REFERENCES, *BENCH_REFERENCES),
    )


def _carrier(reaction: PlannedReaction) -> SequenceRecord:
    """Return the record whose segment moved."""
    return reaction.recombination.moved.record


def _acceptor(reaction: PlannedReaction) -> SequenceRecord:
    """Return the vector whose backbone took it."""
    return reaction.recombination.backbone.record


def _title(lr: PlannedReaction, bp: PlannedReaction | None) -> str:
    """Name the experiment by what goes in and what it ends in."""
    start = _carrier(lr) if bp is None else _carrier(bp)
    return f"Gateway {'LR' if bp is None else 'BP then LR'}: {start.name} into {_acceptor(lr).name}"


def _summary(lr: PlannedReaction, bp: PlannedReaction | None) -> str:
    """One paragraph: what is recombined into what, and in how many stages."""
    entry, destination = _carrier(lr), _acceptor(lr)
    if bp is None:
        return (
            f"Recombine the {lr.recombination.moved.length} bp segment between {entry.name}'s "
            f"attL sites into {destination.name}, in one LR reaction, and select the expression "
            "clone on the destination vector's own marker."
        )
    return (
        f"Recombine the {bp.recombination.moved.length} bp between {_carrier(bp).name}'s attB "
        f"sites into {_acceptor(bp).name}, giving the entry clone {entry.name}; grow that up "
        f"and recombine it into {destination.name}. Two staged reactions, each with its own "
        "incubation, stop, transformation and plate."
    )


def _overview(
    lr: PlannedReaction, bp: PlannedReaction | None, amplicon: Amplicon | None
) -> dict[str, str]:
    """Return the facts to check before starting, each short enough to be a card."""
    entry, product = _carrier(lr), lr.product
    facts = {}
    if amplicon is not None:
        forward, reverse = amplicon.tails
        facts["attB PCR"] = (
            f"{amplicon.length} bp product, {len(forward)} and {len(reverse)} bp tails"
        )
    if bp is not None:
        facts["attB DNA"] = f"{_carrier(bp).name}, {bp.recombination.moved.length} bp insert"
        facts["Donor vector"] = f"{_acceptor(bp).name}, {len(_acceptor(bp))} bp"
    facts["Entry clone"] = f"{entry.name}, {len(entry)} bp"
    facts["Destination vector"] = f"{_acceptor(lr).name}, {len(_acceptor(lr))} bp"
    facts["Reaction"] = (
        f"{LR_CLONASE}, {LR_VOLUME_UL:g} µL at {LR_CELSIUS:g} °C"
        if bp is None
        else f"BP then LR, {LR_VOLUME_UL:g} µL each at {LR_CELSIUS:g} °C"
    )
    facts["Insert"] = f"{lr.recombination.moved.length} bp between the att sites"
    facts["Junctions"] = card(listed([one.name for one in lr.junctions]), "two att sites")
    facts["Product"] = f"{product.name}, {len(product)} bp"
    if lr.phenotype.antibiotic:
        facts["Selection"] = lr.phenotype.antibiotic
    elif lr.phenotype.marker is not None:
        facts["Selection"] = f"{lr.phenotype.marker.name} marker"
    return facts


def _highlights(
    lr: PlannedReaction, bp: PlannedReaction | None, amplicon: Amplicon | None
) -> tuple[str, ...]:
    """Return what the facts mean, a sentence each: what moves, what the scar is, what grows."""
    junctions = lr.junctions
    said = []
    if amplicon is not None:
        forward, reverse = amplicon.tails
        said.append(
            f"{amplicon.template.name or 'The insert'} carries no att site, so it is amplified "
            f"onto attB ends first: each primer is a whole tail -- {len(SPACER)} G residues, "
            f"the {REGION_BP} bp att site and the frame bases this fusion needs, {len(forward)} "
            f"bases forward and {len(reverse)} reverse -- over an annealing region read off the "
            "insert's own end. The G residues leave with the BP by-product, so the entry clone "
            "is the same record as one made from an insert that arrived attB-flanked."
        )
    if bp is not None:
        said.append(
            f"BP runs first: the {bp.recombination.moved.length} bp between the attB sites move "
            f"into {_acceptor(bp).name}, making {_carrier(lr).name}, which is picked, grown up "
            "and purified before LR takes it."
        )
    said.append(
        f"{'LR then moves' if bp is not None else 'One reaction moves'} the "
        f"{lr.recombination.moved.length} bp between the entry clone's att sites into the "
        f"{lr.recombination.backbone.length} bp of destination vector outside its own, and the "
        "ccdB cassette leaves with the by-product."
    )
    said.append(
        f"The junction is not scarless: the clone gains a whole att site at each end of the "
        f"insert, {junctions[0].name} spelling {junctions[0].bases} and {junctions[1].name} "
        f"spelling {junctions[1].bases}. A fusion reads through both."
    )
    said.extend(phenotype_sentences(lr.phenotype, [lr.recombination.moved.name or "the insert"]))
    return tuple(said)


def _materials(
    *,
    lr: PlannedReaction,
    bp: PlannedReaction | None,
    amplicon: Amplicon | None,
    host: str,
) -> tuple[Material, ...]:
    """Every reagent and consumable the protocol asks for, the first reaction's first."""
    entry, destination = _carrier(lr), _acceptor(lr)
    materials: list[Material] = []
    if amplicon is not None:
        materials.extend(
            (
                Material(
                    f"{amplicon.name} primers",
                    storage="-20 °C",
                    amount=f"{PRIMER_STOCK_UM:g} µM each",
                    note="ordered off the sheet below; the manual asks for HPLC- or "
                    "PAGE-purified oligos where colonies are few",
                ),
                Material(
                    f"{amplicon.polymerase.name} DNA Polymerase",
                    storage="-20 °C",
                    note="the manual's answer to an entry clone made of primer-dimers is a "
                    "hot-start enzyme, so use one",
                ),
                Material("PCR and gel cleanup spin columns"),
                Material("Agarose and 1X TAE or TBE"),
            )
        )
    if bp is not None:
        substrate, donor = bp.amounts
        materials.extend(
            (
                Material(
                    f"{_carrier(bp).name} attB DNA",
                    storage="-20 °C",
                    amount=f"{substrate.pmol:g} pmol ({substrate.nanograms:g} ng) per reaction",
                    note="purified, which is what takes the attB primers and their dimers away",
                ),
                Material(
                    f"{_acceptor(bp).name} donor vector",
                    storage="-20 °C",
                    amount=f"{donor.pmol:g} pmol ({donor.nanograms:g} ng) per reaction",
                    note=f"supercoiled, and grown in {PROPAGATION_HOST}",
                ),
                _clonase(BP_CLONASE, BP_CLONASE_CATALOG, BP_CLONASE_UL),
                Material(_plate(bp.phenotype, "donor"), amount="one plate per transformation"),
                Material("Plasmid miniprep kit", note="for the entry clone the LR reaction takes"),
            )
        )
    materials.extend(
        (
            Material(
                f"{entry.name} entry clone",
                storage="-20 °C",
                amount=f"{ENTRY_NG:g} ng per reaction",
                note="supercoiled, which is the substrate the manual calls most efficient"
                if bp is None
                else "the miniprep from the BP plate, which is what this reaction takes",
            ),
            Material(
                f"{destination.name} destination vector",
                storage="-20 °C",
                amount=f"{DESTINATION_NG:g} ng per reaction",
                note=f"grown in {PROPAGATION_HOST}, which is the only kind of strain it grows in",
            ),
            _clonase(LR_CLONASE, LR_CLONASE_CATALOG, LR_CLONASE_UL),
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
                amount=f"{LR_TRANSFORMATION.cells_ul:g} µL per transformation",
                note=f"{CELL_EFFICIENCY_CFU_UG:,} cfu/µg or better, and no F' episome",
            ),
            Material(
                "S.O.C. medium",
                amount=f"{LR_TRANSFORMATION.outgrowth_ul:g} µL per transformation",
            ),
            Material(_plate(lr.phenotype, "destination"), amount="one plate per transformation"),
        )
    )
    return tuple(materials)


def _clonase(name: str, catalog: str, volume_ul: float) -> Material:
    """Return one enzyme mix as a material."""
    return Material(
        name,
        supplier=SUPPLIER,
        catalog=catalog,
        storage="-20 °C",
        amount=f"{volume_ul:g} µL per reaction",
        note="thaw on ice and return it to the freezer at once",
    )


def _plate(phenotype: Phenotype, vector: str) -> str:
    """Return what to pour the selection plates with, named from the vector's own marker."""
    if phenotype.antibiotic:
        antibiotic = phenotype.antibiotic
    elif phenotype.marker is not None:
        antibiotic = f"the antibiotic {phenotype.marker.name} selects"
    else:
        antibiotic = f"the {vector} vector's own antibiotic"
    return f"{phenotype.medium} agar plates with {antibiotic}"


def _oligos(amplicon: Amplicon | None) -> tuple[Oligo, ...]:
    """Return the order sheet the page prints, which is the sheet the plan writes."""
    if amplicon is None:
        return ()
    return tuple(
        oligo_row(report, purpose=pcr_title(amplicon.name), thresholds=amplicon.thresholds)
        for report in amplicon.reports
    )


def _pcr_steps(amplicon: Amplicon | None) -> tuple[Step, ...]:
    """Return the attB PCR and its cleanup, or nothing where the insert arrived attB-flanked."""
    if amplicon is None:
        return ()
    forward, reverse = amplicon.tails
    return (
        pcr_step(
            amplicon.name,
            amplicon.template.name or "the insert",
            amplicon.length,
            polymerase=amplicon.polymerase,
            annealing_temperature=amplicon.report.annealing_temperature,
            extension_seconds=amplicon.report.extension_seconds,
            notes=(
                f"Each primer carries a whole attB tail: {len(SPACER)} G residues, the "
                f"{REGION_BP} bp att site and the frame bases the fusion needs, {len(forward)} "
                f"bases forward and {len(reverse)} reverse. The manual's own cause of few or no "
                "colonies is a tail short of that.",
                "The annealing temperature above is read from the annealing regions alone; a "
                "tail pairs with nothing on the template in the first cycles.",
                "Over 70 bp of primer the manual switches to a two-step adapter PCR, which "
                "installs the same whole tail in two rounds and is not planned here.",
            ),
        ),
        cleanup_step(
            notes=(
                "The BP reaction takes purified attB DNA: gel-purifying the product is the "
                "manual's fix for few or no colonies, and it takes the attB primers and their "
                "dimers away.",
                "An entry clone running as a 2.2 kb supercoiled plasmid is a BP reaction that "
                "cloned attB primer-dimers instead.",
            )
        ),
    )


def _bp_steps(bp: PlannedReaction | None, *, host: str, entry: SequenceRecord) -> tuple[Step, ...]:
    """Return the first stage, or nothing where an entry clone was given."""
    if bp is None:
        return ()
    return (
        Step(
            "Set up the BP reaction",
            instructions=(
                "Thaw the enzyme mix on ice and vortex it briefly twice.",
                "Pipette the attB DNA, the donor vector and the TE buffer into a tube at room "
                "temperature.",
                f"Add {BP_CLONASE_UL:g} µL of {BP_CLONASE}, mix well and spin down.",
            ),
            tables=(bp_reaction(bp.amounts),),
            expected=(f"A {BP_VOLUME_UL:g} µL reaction holding both DNAs.",),
            notes=(
                _pipetting_note("picomoles"),
                f"The manual fixes the picomoles and not the weight: {BP_PMOL * 1000:g} fmol of "
                f"each, the attB DNA down to {BP_SUBSTRATE_MIN_PMOL * 1000:g} fmol, weighed "
                "here from each record's own length.",
                f"Do not go over {BP_DONOR_MAX_NG:g} ng of donor vector or "
                f"{BP_TOTAL_MAX_NG:g} ng of DNA altogether: excess DNA inhibits the reaction.",
                "A linear attB product and a supercoiled donor vector are the substrates the "
                "manual calls most efficient for BP.",
            ),
            troubleshooting=(_volume_trouble(),),
        ),
        Step(
            "Run the BP reaction",
            instructions=(f"Incubate at {BP_CELSIUS:g} °C for {BP_SECONDS // 3600} hour.",),
            timers=(Timer("BP incubation", BP_SECONDS),),
            expected=(
                "Nothing visible: the reaction is a recombination, not a digest.",
                *(_junction_sentence(one, "entry clone") for one in bp.junctions),
            ),
            notes=(
                f"An attB substrate of {BP_LONG_BP:,} bp or more runs up to "
                f"{BP_LONG_SECONDS // 3600} hours instead; efficiency falls as the DNA gets "
                "longer.",
                "Junction positions are 0-based, on the entry clone.",
            ),
            troubleshooting=(
                Troubleshooting(
                    "Few or no colonies later, though the transformation control worked",
                    "Use an attB substrate with a donor vector (attP): the BP reaction takes "
                    "those and no others. Do not freeze and thaw the enzyme mix more than ten "
                    "times.",
                ),
            ),
        ),
        _stop_step("BP"),
        transform_step(
            host,
            bp.phenotype,
            title="Transform and plate the BP reaction",
            inserts=[bp.recombination.moved.name or "the insert"],
            colonies=(
                f"More than {BP_COLONIES:,} colonies where the whole reaction is transformed "
                f"and plated, with cells at {CELL_EFFICIENCY_CFU_UG:,} cfu/µg or better. What "
                "fraction of them is correct is not published; Hartley 2000 counted 195 of 197."
            ),
            protocol=BP_TRANSFORMATION,
            notes=(
                "Unreacted donor vector and the by-product both keep the ccdB gene, which "
                f"kills {host}, so they do not grow. A strain carrying F' would supply ccdA "
                "and cancel that.",
                "Two sizes of colony here mean the donor vector's ccdB gene has mutated or "
                "been deleted; the negative control then gives a similar count.",
            ),
        ),
        _miniprep_step(entry),
    )


def _miniprep_step(entry: SequenceRecord) -> Step:
    """Grow one colony up and purify it, because that is what the LR reaction takes."""
    return Step(
        "Pick and miniprep the entry clone",
        instructions=(
            "Pick single colonies into overnight cultures on the plate's own antibiotic.",
            "Miniprep each culture and measure what it yielded.",
            f"Take {ENTRY_MIN_NG:g}-{ENTRY_NG:g} ng of one prep into the LR reaction below.",
        ),
        expected=(
            f"Supercoiled {entry.name}, concentrated enough to weigh {ENTRY_NG:g} ng into the "
            "volume the LR table leaves for it.",
        ),
        notes=(
            "The LR reaction takes purified entry clone and not the stopped BP reaction: "
            f"{ENTRY_MIN_NG:g}-{ENTRY_NG:g} ng is a weight of miniprep DNA, and no plan can "
            "weigh a yield nobody has measured yet.",
            "Supercoiled miniprep DNA is the substrate the manual calls most efficient for LR.",
            f"The vendor's one-tube protocol chains the two reactions without this step and "
            f"gives {ONE_TUBE_YIELD}, each of which it says to sequence. It runs on its own "
            "timings and is not the protocol below.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The prep is too dilute for the reaction",
                "Concentrate it, or pipette more of it and less TE buffer.",
            ),
        ),
    )


def _lr_steps(lr: PlannedReaction, *, host: str) -> tuple[Step, ...]:
    """Return the LR stage, which every plan runs."""
    return (
        Step(
            "Set up the LR reaction",
            instructions=(
                "Thaw the enzyme mix on ice and vortex it briefly twice.",
                "Pipette the two plasmids and the TE buffer into a tube at room temperature.",
                f"Add {LR_CLONASE_UL:g} µL of {LR_CLONASE}, mix well and spin down.",
            ),
            tables=(lr_reaction(lr.amounts),),
            expected=(f"A {LR_VOLUME_UL:g} µL reaction holding both plasmids.",),
            notes=(
                _pipetting_note("nanograms"),
                f"Do not go over {ENTRY_NG:g} ng of entry clone: the manual reports colonies "
                f"carrying several molecules above it, and fewer colonies below "
                f"{ENTRY_MIN_NG:g} ng.",
                "Supercoiled plasmids are the substrates the manual calls most efficient for LR.",
            ),
            troubleshooting=(_volume_trouble(),),
        ),
        Step(
            "Run the LR reaction",
            instructions=(f"Incubate at {LR_CELSIUS:g} °C for {LR_SECONDS // 3600} hour.",),
            timers=(Timer("LR incubation", LR_SECONDS),),
            expected=(
                "Nothing visible: the reaction is a recombination, not a digest.",
                *(_junction_sentence(one, "expression clone") for one in lr.junctions),
            ),
            notes=(
                f"A plasmid of {LR_LONG_BP:,} bp or more runs up to "
                f"{LR_LONG_SECONDS // 3600} hours instead; efficiency falls as the DNA gets "
                "longer.",
                "Junction positions are 0-based, on the expression clone.",
            ),
            troubleshooting=(
                Troubleshooting(
                    "Few or no colonies later, though the transformation control worked",
                    "Use an entry clone (attL) with a destination vector (attR): the LR "
                    "reaction takes those and no others. Do not freeze and thaw the enzyme mix "
                    "more than ten times.",
                ),
            ),
        ),
        _stop_step("LR"),
        transform_step(
            host,
            lr.phenotype,
            title="Transform and plate the LR reaction",
            inserts=[lr.recombination.moved.name or "the insert"],
            colonies=(
                f"More than {LR_COLONIES:,} colonies where the whole reaction is transformed "
                f"and plated, with cells at {CELL_EFFICIENCY_CFU_UG:,} cfu/µg or better. What "
                "fraction of them is correct is not published; Hartley 2000 counted 96 of 102."
            ),
            protocol=LR_TRANSFORMATION,
            notes=(
                "Unreacted destination vector and the by-product both keep the ccdB gene, "
                f"which kills {host}, so they do not grow. A strain carrying F' would supply "
                "ccdA and cancel that.",
                f"The expression clone lost the chloramphenicol cassette with the by-product, so "
                f"restreaking a colony on {CHLORAMPHENICOL_UG_ML} µg/mL chloramphenicol "
                "confirms it: a true expression clone does not grow there, and one carrying a "
                "mutated ccdB gene does.",
                f"Small colonies beside large ones are usually unreacted "
                f"{_carrier(lr).name} co-transforming; restreak them on the entry clone's own "
                "antibiotic to tell.",
            ),
        ),
    )


def _stop_step(reaction: str) -> Step:
    """End one reaction, which the manual requires before transforming."""
    return Step(
        f"Stop the {reaction} reaction with proteinase K",
        instructions=(
            f"Add {PROTEINASE_K_UL:g} µL of proteinase K at {PROTEINASE_K_UG_UL:g} µg/µL.",
            f"Incubate at {STOP_CELSIUS:g} °C for {STOP_SECONDS // 60} minutes.",
        ),
        timers=(Timer(f"Proteinase K, {reaction}", STOP_SECONDS),),
        expected=("A reaction that can go straight into competent cells.",),
        notes=(
            "The manual lists an untreated reaction as a cause of few or no colonies, so this "
            "step is not optional.",
        ),
    )


def _pipetting_note(units: str) -> str:
    """Say that the table's volumes assume the manual's own concentrations."""
    return (
        f"The volumes above take each DNA at the concentration the manual's own table assumes; "
        f"pipette what your prep needs for the {units} and make the difference up with "
        f"{TE_BUFFER}."
    )


def _volume_trouble() -> Troubleshooting:
    """Return what to do when the DNA will not fit the reaction."""
    return Troubleshooting(
        "The DNA does not fit the reaction volume",
        "Concentrate either DNA, or scale the whole reaction up keeping the enzyme mix at its "
        "stated fraction of the volume.",
    )


def _junction_sentence(junction: Junction, clone: str) -> str:
    """Say what one junction spells and which record gave which side of it."""
    return (
        f"The {clone} carries {junction.name} at {junction.start}, spelling {junction.bases}: "
        f"its first bases from {junction.before} and the rest from {junction.after}."
    )

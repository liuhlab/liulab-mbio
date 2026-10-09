"""The Gateway bench protocol: its own reaction steps, its own notes, and the step order.

The steps any bench shares are `mbio.bench.steps`. This module runs them around each
recombination and adds what only Gateway has to say. A plan that runs BP reads as two staged
reactions: BP, its plate, the miniprep that purifies the entry clone, then LR and its own
plate, with the attB PCR and its cleanup in front of them where one was designed, and the
colony PCR and the sequencing that confirm the expression clone after them. Every number is
computed by `mbio.cloning.gateway.plan` or is one `mbio.cloning.gateway.bench`
cites.
"""

from collections.abc import Mapping, Sequence
from types import MappingProxyType

from mbio import checks as judged
from mbio.bench import REFERENCES as BENCH_REFERENCES
from mbio.bench.oligos import oligo_row
from mbio.bench.pcr import (
    COLONY_PCR_MASTER_MIX,
    PRIMER_STOCK_UM,
    colony_pcr_master_mix_component,
    cycle_citation,
)
from mbio.bench.pcr import SOURCES as PCR_SOURCES
from mbio.bench.phenotype import Phenotype
from mbio.bench.steps import (
    COLONY_PCR_TITLE,
    SEQUENCING_TITLE,
    badges,
    card,
    catalogued,
    cleanup_step,
    colony_pcr_step,
    listed,
    pcr_step,
    pcr_title,
    phenotype_sentences,
    sequencing_step,
    transform_step,
)
from mbio.bench.steps import SOURCES as BENCH_SOURCES
from mbio.bench.validation import ColonyCheck, SangerRead
from mbio.cloning.gateway.att import REGION_BP
from mbio.cloning.gateway.bench import (
    BP_CELSIUS,
    BP_CLONASE,
    BP_CLONASE_CATALOG,
    BP_CLONASE_UL,
    BP_COLONIES,
    BP_LOAD_CAUTION,
    BP_LONG_BP,
    BP_LONG_SECONDS,
    BP_PMOL,
    BP_SECONDS,
    BP_SUBSTRATE_MIN_PMOL,
    BP_TRANSFORMATION,
    BP_VOLUME_UL,
    CELL_EFFICIENCY_CFU_UG,
    CHLORAMPHENICOL_UG_ML,
    DESTINATION_NG,
    DIMER_ENTRY_BP,
    ENTRY_MIN_NG,
    ENTRY_NG,
    LR_CELSIUS,
    LR_CLONASE,
    LR_CLONASE_CATALOG,
    LR_CLONASE_UL,
    LR_COLONIES,
    LR_LOAD_CAUTION,
    LR_LONG_BP,
    LR_LONG_SECONDS,
    LR_SECONDS,
    LR_TRANSFORMATION,
    LR_VOLUME_UL,
    M13_VECTOR_BP,
    ONE_TUBE_YIELD,
    PROPAGATION_HOST,
    PROTEINASE_K_UG_UL,
    PROTEINASE_K_UL,
    REFERENCES,
    SEQUENCING_MAX_PMOL,
    SEQUENCING_MIN_PMOL,
    SEQUENCING_NG,
    STOP_CELSIUS,
    STOP_SECONDS,
    SUPPLIER,
    TE_BUFFER,
    bp_reaction,
    lr_reaction,
)
from mbio.cloning.gateway.bench import SOURCES as REACTION_SOURCES
from mbio.cloning.gateway.design import SPACER, Amplicon, Fusion
from mbio.cloning.gateway.oligos import DesignedOligo
from mbio.cloning.gateway.recombination import Junction, PlannedReaction
from mbio.primers.thresholds import THRESHOLDS_FOR, PrimerRole, Thresholds
from mbio.protocol.model import (
    Citation,
    Material,
    Note,
    Oligo,
    Protocol,
    Source,
    Step,
    Timer,
    Troubleshooting,
    citing,
    number,
)
from mbio.sequence import SequenceRecord, position_text

#: The documents this method's own rows and notes cite, read into
#: `docs/research/gateway-cloning.md`, which names where each was fetched from. The two a
#: caution cites are in `bench`, beside the numbers whose pages they are. The same documents
#: stand in `REFERENCES`, which lists what a run read rather than what a row names.
SOURCES: Mapping[str, Source] = MappingProxyType(
    {
        "11789": Source(
            "Thermo Fisher Scientific Gateway BP Clonase II enzyme mix product sheet",
            edition="11789.II.pps, revision 31 October 2010",
            date="2026-09-18",
            note="docs/research/gateway-cloning.md",
        ),
        "MAN0000291": Source(
            "Thermo Fisher Scientific #MAN0000291 Gateway pDONR Vectors user guide",
            edition="part 25-0531, revised 29 March 2012",
            date="2026-09-18",
            note="docs/research/gateway-cloning.md",
        ),
        "MAN0000437": Source(
            "Thermo Fisher Scientific #MAN0000437 pCR8/GW/TOPO TA Cloning Kit user guide",
            edition="part 25-0706, revised 23 March 2012",
            date="2026-09-18",
            note="docs/research/gateway-cloning.md",
        ),
    }
)

#: The hardware a run needs, which no reagent table covers. Every run screens its colonies by
#: PCR, so the thermocycler and the gel rig are here rather than beside the attB PCR.
EQUIPMENT: tuple[str, ...] = (
    f"Water bath or heat block at {LR_CELSIUS:g} °C",
    f"Water bath or heat block at {STOP_CELSIUS:g} °C",
    f"Heat block or water bath at {LR_TRANSFORMATION.heat_shock_celsius:g} °C",
    f"Shaking incubator and a plate incubator at {LR_TRANSFORMATION.outgrowth_celsius:g} °C",
    "Thermocycler with a heated lid",
    "Agarose gel rig and power supply",
    "Microcentrifuge",
)

#: What a run amplifying its own insert needs on top of `EQUIPMENT`.
PCR_EQUIPMENT: tuple[str, ...] = ("Spectrophotometer or fluorometer",)


def protocol(
    *,
    lr: PlannedReaction,
    bp: PlannedReaction | None,
    amplicon: Amplicon | None = None,
    colony: ColonyCheck,
    reads: Sequence[SangerRead],
    oligos: Sequence[DesignedOligo],
    checks: Sequence[judged.Check],
    host: str,
    fusion: Fusion = "none",
    thresholds: Mapping[PrimerRole, Thresholds] = THRESHOLDS_FOR,
) -> Protocol:
    """Return the bench protocol for one planned Gateway experiment, ready to render.

    Each argument is the `mbio.cloning.gateway.plan.Plan` field or property of that name.
    The steps run in the order someone does them, one reaction at a time: set it up, run it,
    stop it with proteinase K, transform and plate. A planned BP reaction puts its own four
    steps and the miniprep that follows them in front of LR's, an attB PCR puts its own two in
    front of those, and the colony PCR and the sequencing that confirm the clone come last.
    """
    one = Protocol(
        _title(lr, bp),
        summary=_summary(lr, bp),
        overview=_overview(lr, bp, amplicon),
        highlights=_highlights(lr, bp, amplicon),
        checks=badges(checks),
        materials=_materials(lr=lr, bp=bp, amplicon=amplicon, colony=colony, host=host),
        oligos=_oligos(oligos, amplicon, thresholds),
        equipment=EQUIPMENT if amplicon is None else (*PCR_EQUIPMENT, *EQUIPMENT),
        steps=(
            *_pcr_steps(amplicon),
            *_bp_steps(bp, host=host, entry=_carrier(lr)),
            *_lr_steps(lr, host=host),
            *_validation_steps(lr, colony, reads, host=host, fusion=fusion),
        ),
        references=(*REFERENCES, *BENCH_REFERENCES),
        sources={**BENCH_SOURCES, **PCR_SOURCES, **REACTION_SOURCES, **SOURCES},
    )
    return citing(one)


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
    colony: ColonyCheck,
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
                    note="ordered off the oligo sheet on this page; HPLC- or PAGE-purified "
                    "oligos are the published answer where colonies are few",
                    citation=Citation("MAN0000470", "pp. 43-44"),
                ),
                Material(
                    f"{amplicon.polymerase.name} DNA Polymerase",
                    storage="-20 °C",
                    note="a hot-start enzyme is the published answer to an entry clone made "
                    "of primer-dimers, so use one",
                    citation=Citation("MAN0000470", "pp. 43-44"),
                ),
                Material("PCR and gel cleanup spin columns"),
            )
        )
    if bp is not None:
        substrate, donor = bp.amounts
        materials.extend(
            (
                Material(
                    f"{_carrier(bp).name} attB DNA",
                    storage="-20 °C",
                    amount=f"{number(substrate.pmol)} pmol "
                    f"({substrate.nanograms:g} ng) per reaction",
                    note="purified, which is what takes the attB primers and their dimers away",
                ),
                Material(
                    f"{_acceptor(bp).name} donor vector",
                    storage="-20 °C",
                    amount=f"{number(donor.pmol)} pmol ({donor.nanograms:g} ng) per reaction",
                    note=f"supercoiled, and grown in {PROPAGATION_HOST}",
                ),
                _clonase(BP_CLONASE, BP_CLONASE_CATALOG, BP_CLONASE_UL),
                Material(_plate(bp.phenotype, "donor"), amount="one plate per transformation"),
            )
        )
    materials.extend(
        (
            Material(
                f"{entry.name} entry clone",
                storage="-20 °C",
                amount=f"{ENTRY_NG:g} ng per reaction",
                note="supercoiled, which is the most efficient substrate"
                if bp is None
                else "the miniprep from the BP plate, which is what this reaction takes",
                citation=Citation("MAN0001032", "p. 2") if bp is None else None,
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
            catalogued(
                COLONY_PCR_MASTER_MIX,
                storage="-20 °C",
                amount=f"{colony_pcr_master_mix_component().volume_ul:g} µL per reaction",
            ),
            Material("Agarose, 1X TAE or TBE, and a DNA stain"),
            catalogued(colony.ladder.name),
            Material(
                "Plasmid miniprep kit",
                note="for the clones that go to sequencing"
                if bp is None
                else "for the entry clone the LR reaction takes, and for the clones that go to "
                "sequencing",
            ),
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


def _oligos(
    oligos: Sequence[DesignedOligo],
    amplicon: Amplicon | None,
    thresholds: Mapping[PrimerRole, Thresholds],
) -> tuple[Oligo, ...]:
    """Return the order sheet the page prints, which is the sheet the plan writes."""
    return tuple(
        oligo_row(one.report, purpose=_purpose(one, amplicon), thresholds=thresholds[one.role])
        for one in oligos
    )


def _purpose(oligo: DesignedOligo, amplicon: Amplicon | None) -> str:
    """Return the title of the step that uses this oligo."""
    if oligo.role == "amplification":
        return pcr_title(amplicon.name) if amplicon is not None else ""
    return COLONY_PCR_TITLE if oligo.role == "colony PCR" else SEQUENCING_TITLE


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
            cycles=amplicon.polymerase.pcr.cycles,
            cycles_citation=cycle_citation(amplicon.polymerase),
            notes=(
                Note(
                    f"Each primer carries a whole attB tail: {len(SPACER)} G residues, the "
                    f"{REGION_BP} bp att site and the frame bases the fusion needs, "
                    f"{len(forward)} bases forward and {len(reverse)} reverse. A tail short of "
                    "that is a published cause of few or no colonies.",
                    citation=Citation("MAN0000470", "pp. 43-44"),
                ),
                "The annealing temperature above is read from the annealing regions alone; a "
                "tail pairs with nothing on the template in the first cycles.",
                Note(
                    "Over 70 bp of primer a two-step adapter PCR installs the same whole tail "
                    "in two rounds; it is not planned here.",
                    citation=Citation("MAN0000470", "pp. 47-48"),
                ),
            ),
        ),
        cleanup_step(
            notes=(
                Note(
                    "The BP reaction takes purified attB DNA: gel-purifying the product is the "
                    "published fix for few or no colonies, and it takes the attB primers and "
                    "their dimers away.",
                    citation=Citation("MAN0000470", "pp. 43-44"),
                ),
            ),
            troubleshooting=(
                Troubleshooting(
                    f"An entry clone later runs as a {DIMER_ENTRY_BP / 1000:g} kb supercoiled "
                    "plasmid",
                    "The BP reaction cloned attB primer-dimers instead; gel-purify this "
                    "product before running BP again.",
                    citation=Citation("MAN0000470", "pp. 43-44"),
                ),
            ),
        ),
    )


def _bp_steps(bp: PlannedReaction | None, *, host: str, entry: SequenceRecord) -> tuple[Step, ...]:
    """Return the first stage, or nothing where an entry clone was given."""
    if bp is None:
        return ()
    return (
        Step(
            "Set up the BP reaction",
            key="set-up-bp",
            instructions=(
                "Thaw the enzyme mix on ice and vortex it briefly twice.",
                "Pipette the attB DNA, the donor vector and the TE buffer into a tube at room "
                "temperature.",
                f"Add {BP_CLONASE_UL:g} µL of {BP_CLONASE}, mix well and spin down.",
            ),
            cautions=(BP_LOAD_CAUTION,),
            tables=(bp_reaction(bp.amounts),),
            expected=(f"A {BP_VOLUME_UL:g} µL reaction holding both DNAs.",),
            notes=(
                _pipetting_note("picomoles"),
                Note(
                    f"What is fixed is the picomoles and not the weight: {BP_PMOL * 1000:g} "
                    f"fmol of each, the attB DNA down to {BP_SUBSTRATE_MIN_PMOL * 1000:g} fmol, "
                    "weighed here from each record's own length.",
                    citation=Citation("MAN0000470", "pp. 21-22"),
                ),
                Note(
                    "A linear attB product and a supercoiled donor vector are the most "
                    "efficient substrates for BP.",
                    citation=Citation("11789", "p. 2"),
                ),
            ),
            troubleshooting=(_volume_trouble(),),
        ),
        Step(
            "Run the BP reaction",
            key="run-bp",
            instructions=(f"Incubate at {BP_CELSIUS:g} °C for {BP_SECONDS // 3600} hour.",),
            timers=(Timer("BP incubation", BP_SECONDS),),
            expected=(
                "Nothing visible: the reaction is a recombination, not a digest.",
                *(_junction_sentence(one, "entry clone", len(bp.product)) for one in bp.junctions),
            ),
            notes=(
                Note(
                    f"An attB substrate of {BP_LONG_BP:,} bp or more runs up to "
                    f"{BP_LONG_SECONDS // 3600} hours instead; efficiency falls as the DNA gets "
                    "longer.",
                    citation=Citation("MAN0000470", "p. 23"),
                ),
                "Junction positions are 1-based, on the entry clone.",
            ),
            troubleshooting=(
                Troubleshooting(
                    "Few or no colonies later, though the transformation control worked",
                    "Use an attB substrate with a donor vector (attP): the BP reaction takes "
                    "those and no others. Do not freeze and thaw the enzyme mix more than ten "
                    "times.",
                    citation=Citation("MAN0000470", "p. 40"),
                ),
            ),
        ),
        _stop_step("BP"),
        transform_step(
            host,
            bp.phenotype,
            title="Transform and plate the BP reaction",
            key="transform-bp",
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
            ),
            troubleshooting=(
                Troubleshooting(
                    "Two sizes of colony on the plate",
                    "The donor vector's ccdB gene has mutated or been deleted; the negative "
                    "control then gives a similar count.",
                    citation=Citation("MAN0000470", "p. 41"),
                ),
            ),
        ),
        _miniprep_step(entry),
    )


def _miniprep_step(entry: SequenceRecord) -> Step:
    """Grow one colony up and purify it, because that is what the LR reaction takes."""
    return Step(
        "Pick and miniprep the entry clone",
        key="miniprep-entry-clone",
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
            Note(
                "Supercoiled miniprep DNA is the most efficient substrate for LR.",
                citation=Citation("MAN0001032", "p. 2"),
            ),
            Note(
                f"The one-tube protocol chains the two reactions without this step and gives "
                f"{ONE_TUBE_YIELD}, each of which it says to sequence. It runs on its own "
                "timings and is not the protocol below.",
                citation=Citation("MAN0000470", "pp. 45-46"),
            ),
        ),
        troubleshooting=(
            Troubleshooting(
                "The prep is too dilute for the reaction",
                "Concentrate it, or pipette more of it and less TE buffer.",
            ),
            Troubleshooting(
                f"The prep runs as a {DIMER_ENTRY_BP / 1000:g} kb supercoiled plasmid",
                "That size is a BP reaction that cloned attB primer-dimers. Gel-purify the "
                "attB DNA and amplify it with a hot-start polymerase before running BP again.",
                citation=Citation("MAN0000470", "pp. 43-44"),
            ),
        ),
    )


def _lr_steps(lr: PlannedReaction, *, host: str) -> tuple[Step, ...]:
    """Return the LR stage, which every plan runs."""
    return (
        Step(
            "Set up the LR reaction",
            key="set-up-lr",
            instructions=(
                "Thaw the enzyme mix on ice and vortex it briefly twice.",
                "Pipette the two plasmids and the TE buffer into a tube at room temperature.",
                f"Add {LR_CLONASE_UL:g} µL of {LR_CLONASE}, mix well and spin down.",
            ),
            cautions=(LR_LOAD_CAUTION,),
            tables=(lr_reaction(lr.amounts),),
            expected=(f"A {LR_VOLUME_UL:g} µL reaction holding both plasmids.",),
            notes=(
                _pipetting_note("nanograms"),
                Note(
                    f"Above {ENTRY_NG:g} ng of entry clone the colonies carry several "
                    f"molecules, and below {ENTRY_MIN_NG:g} ng there are fewer of them.",
                    citation=Citation("MAN0001032", "p. 2"),
                ),
                Note(
                    "Supercoiled plasmids are the most efficient substrates for LR.",
                    citation=Citation("MAN0001032", "p. 2"),
                ),
            ),
            troubleshooting=(_volume_trouble(),),
        ),
        Step(
            "Run the LR reaction",
            key="run-lr",
            instructions=(f"Incubate at {LR_CELSIUS:g} °C for {LR_SECONDS // 3600} hour.",),
            timers=(Timer("LR incubation", LR_SECONDS),),
            expected=(
                "Nothing visible: the reaction is a recombination, not a digest.",
                *(
                    _junction_sentence(one, "expression clone", len(lr.product))
                    for one in lr.junctions
                ),
            ),
            notes=(
                Note(
                    f"A plasmid of {LR_LONG_BP:,} bp or more runs up to "
                    f"{LR_LONG_SECONDS // 3600} hours instead; efficiency falls as the DNA gets "
                    "longer.",
                    citation=Citation("MAN0000470", "p. 32"),
                ),
                "Junction positions are 1-based, on the expression clone.",
            ),
            troubleshooting=(
                Troubleshooting(
                    "Few or no colonies later, though the transformation control worked",
                    "Use an entry clone (attL) with a destination vector (attR): the LR "
                    "reaction takes those and no others. Do not freeze and thaw the enzyme mix "
                    "more than ten times.",
                    citation=Citation("MAN0000470", "p. 40"),
                ),
            ),
        ),
        _stop_step("LR"),
        transform_step(
            host,
            lr.phenotype,
            title="Transform and plate the LR reaction",
            key="transform-lr",
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
            ),
            troubleshooting=(
                Troubleshooting(
                    "Small colonies beside large ones",
                    f"Usually unreacted {_carrier(lr).name} co-transforming; restreak them on "
                    "the entry clone's own antibiotic to tell.",
                    citation=Citation("MAN0000470", "p. 41"),
                ),
            ),
        ),
    )


def _validation_steps(
    lr: PlannedReaction,
    colony: ColonyCheck,
    reads: Sequence[SangerRead],
    *,
    host: str,
    fusion: Fusion,
) -> tuple[Step, ...]:
    """Return what tells a correct expression clone from the plate's other colonies.

    Both are read on the LR product: the attB junctions are what a colony PCR crosses and what
    a sequencing read has to cover, and neither exists until LR has run.
    """
    entry, insert = _carrier(lr).name, lr.recombination.moved.name or "the insert"
    return (
        colony_pcr_step(
            colony,
            junctions=len(lr.junctions),
            notes=(
                "Both primers sit in the destination vector's backbone, so one product crosses "
                "both attB junctions and an unrecombined vector gives a band of its own.",
                f"An unrecombined destination vector still carries ccdB, which kills {host}, so "
                "its lane is what a strain supplying ccdA or a damaged ccdB gene would put on "
                "the plate.",
                Note(
                    "A restriction digest beside this screen is recommended the first time it "
                    "is run: mispriming and contaminating template both give artefacts.",
                    citation=Citation("MAN0000291", "p. 10"),
                ),
            ),
            troubleshooting=(
                Troubleshooting(
                    "The plate carries large colonies and small ones",
                    f"The small ones are usually unreacted {entry} co-transforming. It carries "
                    "neither of these primers, so it adds no band here; restreak on the entry "
                    "clone's own antibiotic to tell.",
                    citation=Citation("MAN0000470", "p. 41"),
                ),
            ),
        ),
        sequencing_step(
            reads,
            junctions=[one.bases for one in lr.junctions],
            inserts=[insert],
            instructions=(
                f"Send at least {SEQUENCING_NG:g} ng of plasmid with "
                f"{SEQUENCING_MIN_PMOL:g}-{SEQUENCING_MAX_PMOL:g} pmol of each primer.",
            ),
            notes=(
                Note(
                    "No kit primer reads these junctions: GW1 and GW2 suit pCR8/GW/TOPO alone, "
                    f"and an M13 primer crosses at least {M13_VECTOR_BP} bp of vector first. "
                    "Both primers above are designed against this record.",
                    citation=Citation("MAN0000437", "p. 13"),
                ),
                *(
                    ()
                    if fusion == "none"
                    else (
                        f"This is a {fusion} fusion, so read the att junction it runs through "
                        "and check the frame the plan judged is the frame on the trace.",
                    )
                ),
            ),
        ),
    )


def _stop_step(reaction: str) -> Step:
    """End one reaction, which the manual requires before transforming."""
    return Step(
        f"Stop the {reaction} reaction with proteinase K",
        key=f"stop-{reaction}",
        instructions=(
            f"Add {PROTEINASE_K_UL:g} µL of proteinase K at {PROTEINASE_K_UG_UL:g} µg/µL.",
            f"Incubate at {STOP_CELSIUS:g} °C for {STOP_SECONDS // 60} minutes.",
        ),
        timers=(Timer(f"Proteinase K, {reaction}", STOP_SECONDS),),
        expected=("A reaction that can go straight into competent cells.",),
        notes=(
            Note(
                "An untreated reaction is a published cause of few or no colonies, so this "
                "step is not optional.",
                citation=Citation("MAN0000470", "p. 40"),
            ),
        ),
    )


def _pipetting_note(units: str) -> str:
    """Say that the table's volumes assume the concentrations the table itself prints."""
    return (
        f"The volumes above take each DNA at the concentration the table assumes; pipette what "
        f"your prep needs for the {units} and make the difference up with {TE_BUFFER}."
    )


def _volume_trouble() -> Troubleshooting:
    """Return what to do when the DNA will not fit the reaction."""
    return Troubleshooting(
        "The DNA does not fit the reaction volume",
        "Concentrate either DNA, or scale the whole reaction up keeping the enzyme mix at its "
        "stated fraction of the volume.",
    )


def _junction_sentence(junction: Junction, clone: str, length: int) -> str:
    """Say what one junction spells and which record gave which side of it."""
    return (
        f"The {clone} carries {junction.name} at {position_text(junction.start, length)}, "
        f"spelling {junction.bases}: "
        f"its first bases from {junction.before} and the rest from {junction.after}."
    )

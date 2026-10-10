"""The Golden Gate bench protocol: its own two steps, its own notes, and the step order.

The steps any bench shares are `mbio.bench.steps`. This module runs them around the
assembly and its cycling and adds what only Golden Gate has to say. Every number is computed by
`mbio.cloning.goldengate.plan` or by the modules it calls; the constants below are the choices no
table of NEB's covers, and each says where it comes from.
"""

from collections.abc import Mapping, Sequence

from mbio import checks as judged
from mbio.bench.amounts import Amount
from mbio.bench.gels import choose_ladder
from mbio.bench.goldengate import (
    GOLDEN_GATE_PCR_CYCLES,
    GOLDEN_GATE_PCR_CYCLES_CITATION,
    assembly_program,
    assembly_reaction,
    enzyme_component,
    fidelity_citation,
    fidelity_sources,
    golden_gate_temperature,
    ligase_master_mix_component,
)
from mbio.bench.goldengate import SOURCES as GOLDEN_GATE_SOURCES
from mbio.bench.inactivation import heat_inactivation
from mbio.bench.oligos import oligo_row
from mbio.bench.pcr import (
    COLONY_PCR_MASTER_MIX,
    DNTP_STOCK_MM,
    colony_pcr_master_mix_component,
)
from mbio.bench.pcr import SOURCES as PCR_SOURCES
from mbio.bench.phenotype import Phenotype
from mbio.bench.steps import (
    CELLS_UL,
    COLONY_PCR_TITLE,
    DPNI_UNITS,
    HEAT_SHOCK_CELSIUS,
    MINIPREP_KIT,
    OUTGROWTH_CELSIUS,
    OUTGROWTH_UL,
    QUANTIFY_EQUIPMENT,
    SEQUENCING_TITLE,
    badges,
    card,
    catalogued,
    cleanup_step,
    colony_pcr_step,
    dam_sites,
    dpni_step,
    enzyme_material,
    gel_step,
    listed,
    pcr_step,
    pcr_title,
    phenotype_sentences,
    plate_material,
    quantify_step,
    sequencing_step,
    transform_step,
)
from mbio.bench.steps import SOURCES as BENCH_SOURCES
from mbio.bench.validation import ColonyCheck, SangerRead
from mbio.cloning.goldengate.assembly import (
    Assembly,
    Junction,
    Part,
    overhang_name,
    site_name,
)
from mbio.cloning.goldengate.design import OverhangSet
from mbio.cloning.goldengate.oligos import DesignedOligo
from mbio.cloning.plan import (
    JOIN_SECTION,
    MAKE_SECTION,
    PRODUCT_FILE,
    SCREEN_SECTION,
    amplicon_files,
)
from mbio.enzymes import Enzyme
from mbio.overhangs import FidelityReport
from mbio.primers.polymerase import Polymerase
from mbio.primers.thresholds import PrimerRole, Thresholds
from mbio.protocol.figures import tail_figure
from mbio.protocol.model import (
    Citation,
    Component,
    Figure,
    Material,
    Note,
    Protocol,
    Step,
    Troubleshooting,
    citing,
    figured,
    sectioned,
)
from mbio.sequence import SequenceRecord, position_text

#: Correct colonies NEB counts from a single-insert assembly with 2.5 µL of the outgrowth
#: plated (E1601 technical note, `docs/research/golden-gate-assembly.md` §5). An
#: order-of-magnitude check, not a target.
NEB_COLONIES = 687

#: Who sells the products this protocol names. Every catalogue number it prints comes out of an
#: enzyme record or a product name a supplier wrote; none is written here.
SUPPLIER = "New England Biolabs"

#: The strain a protocol names unless the caller picks one. Blue/white screening needs a host
#: that supplies the rest of the lacZ fragment the vector carries, which this one does.
DEFAULT_HOST = "NEB 5-alpha Competent E. coli (C2987)"

#: The catalogue number NEB sells NEBridge Ligase Master Mix under.
LIGASE_MIX_CATALOG = "M1100"

#: The hardware a run needs, which no reagent table covers.
EQUIPMENT: tuple[str, ...] = (
    "Thermocycler with a heated lid",
    "Agarose gel rig and power supply",
    "Microcentrifuge",
    QUANTIFY_EQUIPMENT,
    f"Heat block or water bath at {HEAT_SHOCK_CELSIUS:g} °C",
    f"Shaking incubator and a plate incubator at {OUTGROWTH_CELSIUS:g} °C",
)


def protocol(
    *,
    vector: SequenceRecord,
    span: tuple[int, int],
    overhangs: OverhangSet,
    linearised_vector: Part,
    insert_parts: Sequence[Part],
    assembly: Assembly,
    colony: ColonyCheck,
    reads: Sequence[SangerRead],
    amounts: tuple[Amount, ...],
    phenotype: Phenotype,
    oligos: Sequence[DesignedOligo],
    checks: Sequence[judged.Check],
    host: str,
    polymerase: Polymerase,
    cleanup_kit: Material,
    thresholds: Mapping[PrimerRole, Thresholds],
) -> Protocol:
    """Return the bench protocol for one planned assembly, ready to render.

    Each argument is the `mbio.cloning.goldengate.plan.Plan` field or property of that name, and
    `oligos` is `Plan.designed_oligos`. The steps run in the order someone does them: one PCR
    per part, the gel that checks them, the DpnI digest and cleanup, quantification, the
    assembly, transformation and plating, colony PCR, and sequencing. Every step is per
    experiment, however many parts there are.
    """
    parts = (linearised_vector, *insert_parts)
    names = tuple(part.name for part in insert_parts)
    enzyme = assembly.enzyme
    inserts = listed(names)
    one = Protocol(
        f"Golden Gate assembly: {inserts} into {vector.name}",
        summary=(
            f"Open {vector.name} by PCR across {span[0]}-{span[1]}, amplify "
            f"{inserts} with {enzyme.name} tails, join the {len(parts)} fragments in "
            "one Golden Gate reaction, and confirm the clone by colony PCR and sequencing."
        ),
        overview=_overview(vector, insert_parts, assembly, overhangs, phenotype),
        highlights=_highlights(parts, phenotype, names),
        checks=badges(checks),
        materials=_materials(
            vector=vector,
            parts=parts,
            inserts=names,
            colony=colony,
            enzyme=enzyme,
            host=host,
            polymerase=polymerase,
            phenotype=phenotype,
            cleanup_kit=cleanup_kit,
        ),
        oligos=tuple(
            oligo_row(oligo.report, purpose=_purpose(oligo), thresholds=thresholds[oligo.role])
            for oligo in oligos
        ),
        equipment=EQUIPMENT,
        steps=_steps(
            parts=parts,
            inserts=names,
            assembly=assembly,
            overhangs=overhangs,
            amounts=amounts,
            colony=colony,
            reads=reads,
            phenotype=phenotype,
            host=host,
            polymerase=polymerase,
            cleanup_kit=cleanup_kit,
        ),
        sources={
            **BENCH_SOURCES,
            **PCR_SOURCES,
            **GOLDEN_GATE_SOURCES,
            **fidelity_sources(overhangs.fidelity),
        },
    )
    return citing(one)


def _overview(
    vector: SequenceRecord,
    insert_parts: Sequence[Part],
    assembly: Assembly,
    overhangs: OverhangSet,
    phenotype: Phenotype,
) -> dict[str, str]:
    """Return the facts to check before starting, each short enough to be a card."""
    fidelity = overhangs.fidelity
    enzyme = assembly.enzyme
    facts = {
        "Vector": f"{vector.name}, {len(vector)} bp",
        "Insert" if len(insert_parts) == 1 else "Inserts": _insert_fact(insert_parts),
        "Enzyme": f"{enzyme.supplier_label} at {golden_gate_temperature(enzyme):g} °C",
        "Fragments": f"{len(insert_parts) + 1} in one reaction",
        "Overhangs": _brief(overhangs.overhangs, "junctions"),
        "Fidelity": f"{fidelity.value:.0%}, {fidelity.label}",
        "Product": f"{assembly.product.name}, {len(assembly.product)} bp",
    }
    selection = _selection(phenotype)
    if selection:
        facts["Selection"] = selection
    return facts


def _insert_fact(insert_parts: Sequence[Part]) -> str:
    """Return the insert card: the one insert with its length, or what the inserts are called."""
    if len(insert_parts) == 1:
        return f"{insert_parts[0].name}, {len(insert_parts[0].template)} bp"
    return _brief([part.name for part in insert_parts], "inserts")


def _selection(phenotype: Phenotype) -> str:
    """Return what to select transformants on, or nothing where the vector annotates no marker."""
    if phenotype.antibiotic:
        return phenotype.antibiotic
    return f"{phenotype.marker.name} marker" if phenotype.marker is not None else ""


def _brief(items: Sequence[str], noun: str) -> str:
    """List these where a card holds the list, and count them where it does not.

    Examples
    --------
    >>> _brief(("ATGA", "TGGC"), "junctions")
    'ATGA and TGGC'
    """
    return card(listed(items), f"{len(items)} {noun}")


def _highlights(
    parts: Sequence[Part], phenotype: Phenotype, inserts: Sequence[str]
) -> tuple[str, ...]:
    """Return what the facts mean, a sentence each: the fragments, then the phenotype."""
    joined = listed([f"{part.name} ({part.fragment_length} bp)" for part in parts])
    return (
        f"One reaction joins {len(parts)} fragments: {joined}.",
        *phenotype_sentences(phenotype, inserts),
    )


def _materials(
    *,
    vector: SequenceRecord,
    parts: Sequence[Part],
    inserts: Sequence[str],
    colony: ColonyCheck,
    enzyme: Enzyme,
    host: str,
    polymerase: Polymerase,
    phenotype: Phenotype,
    cleanup_kit: Material,
) -> tuple[Material, ...]:
    """Every reagent and consumable the protocol asks for. The oligos are the order sheet."""
    fragments = len(parts)
    mix = ligase_master_mix_component(fragments)
    ladders = dict.fromkeys(
        (choose_ladder(tuple(part.length for part in parts)).name, colony.ladder.name)
    )
    return (
        Material(f"{vector.name} plasmid", storage="-20 °C", note="PCR template"),
        *(Material(f"{name} template", storage="-20 °C", note="PCR template") for name in inserts),
        Material(
            f"{polymerase.name} DNA Polymerase and its reaction buffer",
            supplier=SUPPLIER,
            storage="-20 °C",
        ),
        Material("dNTP mix", storage="-20 °C", note=f"{DNTP_STOCK_MM:g} mM of each base"),
        Material(
            "DpnI",
            supplier=SUPPLIER,
            storage="-20 °C",
            amount=f"{DPNI_UNITS} units per PCR",
            note="cuts the methylated plasmid template only",
        ),
        cleanup_kit,
        Material(
            mix.name,
            supplier=SUPPLIER,
            catalog=LIGASE_MIX_CATALOG,
            storage="-20 °C",
            amount=_per_reaction(mix),
        ),
        enzyme_material(enzyme, amount=_per_reaction(enzyme_component(enzyme, fragments))),
        catalogued(
            host,
            supplier=SUPPLIER if host == DEFAULT_HOST else "",
            storage="-80 °C",
            amount=f"{CELLS_UL:g} µL per transformation",
        ),
        Material(
            "SOC or NEB 10-beta/Stable Outgrowth Medium",
            amount=f"{OUTGROWTH_UL:g} µL per transformation",
        ),
        plate_material(phenotype),
        MINIPREP_KIT,
        catalogued(
            COLONY_PCR_MASTER_MIX,
            supplier=SUPPLIER,
            storage="-20 °C",
            amount=_per_reaction(colony_pcr_master_mix_component()),
        ),
        Material("Agarose and 1X TAE or TBE"),
        *(catalogued(name, supplier=SUPPLIER) for name in ladders),
    )


def _per_reaction(component: Component) -> str:
    """Return what one reaction takes of this component."""
    return f"{component.volume_ul:g} µL per reaction"


def _purpose(oligo: DesignedOligo) -> str:
    """Return the title of the step that uses this oligo."""
    if oligo.part is not None:
        return pcr_title(oligo.part.name)
    return COLONY_PCR_TITLE if oligo.role == "colony PCR" else SEQUENCING_TITLE


def _steps(
    *,
    parts: Sequence[Part],
    inserts: Sequence[str],
    assembly: Assembly,
    overhangs: OverhangSet,
    amounts: tuple[Amount, ...],
    colony: ColonyCheck,
    reads: Sequence[SangerRead],
    phenotype: Phenotype,
    host: str,
    polymerase: Polymerase,
    cleanup_kit: Material,
) -> tuple[Step, ...]:
    """Return the steps in the order they happen, the shared ones carrying Golden Gate's notes."""
    enzyme = assembly.enzyme
    cut = [part for part in parts if part.dpni]
    files = amplicon_files(part.name for part in parts)
    made = [
        figured(_pcr_step(part, polymerase), _tail_figure(part, enzyme, path))
        for part, path in zip(parts, files, strict=True)
    ]
    made.append(gel_step([(part.name, part.length) for part in parts]))
    if cut:
        made.append(
            dpni_step(
                [part.name for part in cut],
                [(part.template.name, dam_sites(part.template)) for part in cut],
                notes=(
                    Note(
                        "The kit's own protocol has no step for a methylated template, so the "
                        "incubation is this package's choice.",
                        citation=Citation("E1601"),
                    ),
                ),
            )
        )
    made.append(
        cleanup_step(
            kit=cleanup_kit,
            notes=(
                Note(
                    "The reaction takes purified amplicons: polymerase carried over from the "
                    "PCR fills in the four-base overhangs, which blunts the ends and "
                    "mis-assembles them.",
                    citation=Citation("E1601"),
                ),
            ),
        )
    )
    made.append(quantify_step(amounts))
    return (
        *sectioned(MAKE_SECTION, *made),
        *sectioned(
            JOIN_SECTION,
            _assembly_step(enzyme, amounts),
            figured(
                _cycling_step(
                    enzyme,
                    len(parts),
                    assembly.junctions,
                    len(assembly.product),
                    overhangs.fidelity,
                ),
                _assembly_figure(assembly),
            ),
        ),
        *sectioned(
            SCREEN_SECTION,
            transform_step(
                host,
                phenotype,
                inserts=inserts,
                colonies=f"Hundreds of colonies; NEB counts about {NEB_COLONIES} correct ones "
                "from a single-insert assembly with 2.5 µL of the outgrowth plated.",
            ),
            colony_pcr_step(
                colony,
                junctions=len(assembly.junctions),
                troubleshooting=(
                    Troubleshooting(
                        "Every colony reads as empty vector",
                        "The template survived the DpnI digest, or the vector re-closed; check "
                        "the 60 °C soak ran.",
                        citation=Citation("E1601", "FAQ 10"),
                    ),
                ),
            ),
            sequencing_step(reads, junctions=overhangs.overhangs, inserts=inserts),
        ),
    )


def _pcr_step(part: Part, polymerase: Polymerase) -> Step:
    """Amplify one part with the tails that carry the enzyme site."""
    return pcr_step(
        part.name,
        part.template.name,
        part.length,
        polymerase=polymerase,
        annealing_temperature=part.report.annealing_temperature,
        extension_seconds=part.report.extension_seconds,
        cycles=GOLDEN_GATE_PCR_CYCLES,
        cycles_citation=GOLDEN_GATE_PCR_CYCLES_CITATION,
        notes=(
            "The cycle count is the fewest enough for an amplicon going into an assembly; "
            "fewer cycles means fewer PCR errors.",
        ),
    )


def _tail_figure(part: Part, enzyme: Enzyme, path: str) -> Figure:
    """Return the end of the part's amplicon its forward primer makes, and what the cut leaves.

    That end's overhang is the part's own, the other belonging to the next part round the
    circle, so across the parts' steps every junction is drawn once.
    """
    left, right = part.left_overhang, part.right_overhang
    return tail_figure(
        part.amplicon,
        path=path,
        primer=part.forward.name,
        enzymes=(enzyme.name,),
        caption=(
            f"The end of the {part.name} amplicon that {part.forward.name} makes. Its tail "
            f"carries a {enzyme.name} site pointing into the fragment, and {enzyme.name} cuts the "
            f"two strands {len(left)} bases apart, leaving {left} single-stranded. The other "
            f"end leaves {right} the same way."
        ),
        highlight=(site_name(enzyme), overhang_name(left)),
    )


def _assembly_step(enzyme: Enzyme, amounts: tuple[Amount, ...]) -> Step:
    """Set the one-tube digest and ligation up."""
    table = assembly_reaction(enzyme, amounts, measured=True)
    total = sum(component.volume_ul for component in table.components)
    return Step(
        "Set up the Golden Gate reaction",
        key="set-up-assembly",
        instructions=(
            "Thaw the master mix on ice and mix it well; it is viscous.",
            "Pipette the volume that gives each fragment its picomoles, and make the "
            "difference up with water.",
            "Pipette the DNA into the tube first, then the rest.",
            "Mix gently and spin down.",
        ),
        tables=(table,),
        expected=(f"A {total:g} µL reaction holding every fragment.",),
        troubleshooting=(
            Troubleshooting(
                "The DNA does not fit the reaction volume",
                "Concentrate the fragments, or scale the whole reaction up.",
                citation=Citation("E1601", "assembly reaction, note 6"),
            ),
        ),
    )


def _cycling_step(
    enzyme: Enzyme,
    fragments: int,
    junctions: Sequence[Junction],
    length: int,
    fidelity: FidelityReport,
) -> Step:
    """Run it, and the heat inactivation the supplier gives."""
    scored = fidelity_citation(fidelity)
    programs = [assembly_program(enzyme, fragments=fragments)]
    kill = heat_inactivation(enzyme)
    if kill is not None:
        programs.append(kill)
    expected = [
        "Nothing visible. The 60 °C soak at the end is a digest, not heat inactivation: it "
        "cuts vector that never opened or has closed again, so fewer empty colonies grow.",
        *(
            f"The product carries {one.overhang} at "
            f"{position_text(one.start, length)}, joining {one.before} to "
            f"{one.after}."
            for one in junctions
        ),
    ]
    return Step(
        "Run the Golden Gate program",
        key="run-assembly",
        instructions=("Put the tube in the thermocycler and run the program below.",),
        programs=tuple(programs),
        expected=tuple(expected),
        notes=(
            *(
                (
                    Note(
                        "The fidelity this design was scored on is read from measured "
                        "ligation counts.",
                        citation=scored,
                    ),
                )
                if scored
                else ()
            ),
            "Junction positions are 1-based, on the product.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Mostly empty vector later",
                "Keep the 60 °C soak, and check the template was digested with DpnI.",
                citation=Citation("E1601", "FAQ 10"),
            ),
            Troubleshooting(
                "Few colonies later",
                "Raise the cycle count, or plate more of the outgrowth.",
                citation=Citation("E1601", "FAQ 14"),
            ),
        ),
    )


def _assembly_figure(assembly: Assembly) -> Figure:
    """Return the plasmid the reaction closes, every junction it spells lit."""
    enzyme = assembly.enzyme.name
    spelled = listed([one.overhang for one in assembly.junctions])
    return Figure(
        (PRODUCT_FILE,),
        f"{assembly.product.name} as the reaction closes it: the {spelled} junctions it spells, "
        f"and every {enzyme} site left on it.",
        enzymes=(enzyme,),
        highlight=tuple(dict.fromkeys(one.feature_name for one in assembly.junctions)),
    )

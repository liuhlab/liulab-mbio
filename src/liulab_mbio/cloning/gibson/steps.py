"""The Gibson bench protocol: its own two steps, its own notes, and the step order.

The steps any bench shares are `liulab_mbio.bench.steps`. This module runs them around the
assembly and its incubation and adds what only this method has to say. Every number is computed
by `liulab_mbio.cloning.gibson.plan` or by the modules it calls; the sentences below that carry
a supplier's number name the section of ``docs/research/gibson-assembly.md`` it comes from.
"""

from collections.abc import Mapping, Sequence

from liulab_mbio import checks as judged
from liulab_mbio.bench import REFERENCES as BENCH_REFERENCES
from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.gels import choose_ladder
from liulab_mbio.bench.oligos import oligo_row
from liulab_mbio.bench.pcr import (
    COLONY_PCR_MASTER_MIX,
    DNTP_STOCK_MM,
    colony_pcr_master_mix_component,
)
from liulab_mbio.bench.phenotype import Phenotype
from liulab_mbio.bench.steps import (
    CELLS_UL,
    COLONY_PCR_TITLE,
    DPNI_REFERENCE,
    DPNI_UNITS,
    HEAT_SHOCK_CELSIUS,
    IPTG_UM,
    OUTGROWTH_CELSIUS,
    OUTGROWTH_UL,
    PLATE_REFERENCE,
    SEQUENCING_TITLE,
    XGAL_UG_ML,
    badges,
    card,
    catalogued,
    cleanup_step,
    colony_pcr_step,
    dam_sites,
    dpni_step,
    gel_step,
    listed,
    pcr_step,
    pcr_title,
    phenotype_sentences,
    quantify_step,
    sequencing_step,
    transform_step,
)
from liulab_mbio.bench.validation import ColonyCheck, SangerRead
from liulab_mbio.cloning.gibson.assembly import Assembly, Junction, Part
from liulab_mbio.cloning.gibson.bench import (
    AssemblyProduct,
    assembly_program,
    assembly_reaction,
)
from liulab_mbio.cloning.gibson.oligos import DesignedOligo
from liulab_mbio.primers.polymerase import Polymerase
from liulab_mbio.primers.thresholds import PrimerRole, Thresholds
from liulab_mbio.protocol.model import (
    Incubation,
    Material,
    Protocol,
    Reference,
    Step,
    Timer,
    Troubleshooting,
)
from liulab_mbio.sequence import SequenceRecord

#: The strain a protocol names unless the caller picks one. Blue/white screening needs a host
#: that supplies the rest of the lacZ fragment the vector carries, which this one does.
DEFAULT_HOST = "NEB 5-alpha Competent E. coli (C2987)"

#: The DpnI digest NEB prescribes for this method, unlike Golden Gate's, which is this package's
#: own: 20 units, 30 minutes at 37 °C, then 20 minutes at 80 °C (note §11).
DPNI_DIGEST_SECONDS = 1800
DPNI_INACTIVATION = Incubation("Heat-inactivate DpnI", 80.0, 1200)

#: What NEB's own lot test yields, which is the firmest colony count the note carries: more than
#: 100 colonies from six fragments at 0.05 pmol each with a tenth of the outgrowth plated
#: (note §16). An order-of-magnitude check, not a target.
RELEASE_COLONIES = 100

#: Where NEB stops calling a column optional: this many PCR fragments, or a fragment this many
#: kilobases long, at which it puts the gain at two- to tenfold (note §11).
CLEANUP_FRAGMENTS = 3
CLEANUP_KB = 5

#: Colonies to screen, and what Takara read of that many across its In-Fusion series: all of
#: them correct at two fragments, falling to `CORRECT_AT_FIVE` at five. The only source that
#: measures the fall with fragment count, and In-Fusion's number and not NEB's (note §16).
SCREENED_COLONIES = 10
CORRECT_AT_FIVE = 4

#: What Gibson 2009 found of the junctions it sequenced: about one error per this many molecules
#: joined, which is why a clone is sequenced whatever the screen said (note §16).
MOLECULES_PER_ERROR = 50

#: What a protocol cites for how many screened colonies read correct and for how often a
#: junction is misjoined (note §16). Every plan screens and sequences, so both are always cited.
SCREENING_REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Takara Bio, In-Fusion Cloning FAQs, for the correct clones out of ten screened from "
        "two to five fragments",
        url="https://www.takarabio.com/learning-centers/cloning/in-fusion-cloning-faqs",
    ),
    Reference(
        "Gibson, D.G. et al. (2009) Enzymatic assembly of DNA molecules up to several hundred "
        "kilobases. Nat. Methods 6, 343-345, for the junction error rate",
        url="https://doi.org/10.1038/nmeth.1318",
    ),
)

#: The hardware a run needs, which no reagent table covers.
EQUIPMENT: tuple[str, ...] = (
    "Thermocycler with a heated lid",
    "Agarose gel rig and power supply",
    "Microcentrifuge",
    "Spectrophotometer or fluorometer",
    f"Heat block or water bath at {HEAT_SHOCK_CELSIUS:g} °C",
    f"Shaking incubator and a plate incubator at {OUTGROWTH_CELSIUS:g} °C",
)


def protocol(
    *,
    vector: SequenceRecord,
    span: tuple[int, int],
    product: AssemblyProduct,
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
    thresholds: Mapping[PrimerRole, Thresholds],
) -> Protocol:
    """Return the bench protocol for one planned assembly, ready to render.

    Each argument is the `liulab_mbio.cloning.gibson.plan.Plan` field or property of that name,
    and `oligos` is `Plan.designed_oligos`. The steps run in the order someone does them: one
    PCR per part, the gel that checks them, the DpnI digest, the cleanup and the quantification,
    the assembly and its incubation, the transformation, then the colony PCR and the sequencing
    that settle it. Every step is per experiment, however many parts there are.
    """
    parts = (linearised_vector, *insert_parts)
    names = tuple(part.name for part in insert_parts)
    inserts = listed(names)
    opening = (
        f"Open {vector.name} by PCR across {span[0]}-{span[1]}"
        if linearised_vector.amplified
        else f"Take {vector.name}, which is already linear, as the backbone"
    )
    return Protocol(
        f"Gibson assembly: {inserts} into {vector.name}",
        summary=(
            f"{opening}, amplify {inserts} with the overlap at each junction carried as a "
            f"primer tail, join the {len(parts)} fragments in one {product.name} reaction, and "
            "confirm the clone by colony PCR and sequencing."
        ),
        overview=_overview(vector, insert_parts, assembly, product, phenotype),
        highlights=_highlights(parts, assembly, phenotype, names),
        checks=badges(checks),
        materials=_materials(
            vector=vector,
            parts=parts,
            inserts=names,
            product=product,
            colony=colony,
            host=host,
            polymerase=polymerase,
            phenotype=phenotype,
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
            product=product,
            amounts=amounts,
            colony=colony,
            reads=reads,
            phenotype=phenotype,
            host=host,
            polymerase=polymerase,
        ),
        references=_references(parts, product, phenotype),
    )


def _overview(
    vector: SequenceRecord,
    insert_parts: Sequence[Part],
    assembly: Assembly,
    product: AssemblyProduct,
    phenotype: Phenotype,
) -> dict[str, str]:
    """Return the facts to check before starting, each short enough to be a card."""
    tier = product.tier(len(assembly.parts))
    facts = {
        "Vector": f"{vector.name}, {len(vector)} bp",
        "Insert" if len(insert_parts) == 1 else "Inserts": _insert_fact(insert_parts),
        "Fragments": f"{len(assembly.parts)} in one reaction",
        "Overlaps": _brief([f"{one.length} bp" for one in assembly.junctions], "junctions"),
        "Assembly": f"{product.name}, {product.celsius:g} °C for "
        f"{tier.incubation_seconds // 60} min",
        "Product": f"{assembly.product.name}, {len(assembly.product)} bp",
    }
    selection = _selection(phenotype)
    if selection:
        facts["Selection"] = selection
    return facts


def _insert_fact(insert_parts: Sequence[Part]) -> str:
    """Return the insert card: the one insert with its length, or what the inserts are called."""
    if len(insert_parts) == 1:
        return f"{insert_parts[0].name}, {insert_parts[0].fragment_length} bp"
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
    >>> _brief(("18 bp", "20 bp"), "junctions")
    '18 bp and 20 bp'
    """
    return card(listed(items), f"{len(items)} {noun}")


def _highlights(
    parts: Sequence[Part], assembly: Assembly, phenotype: Phenotype, inserts: Sequence[str]
) -> tuple[str, ...]:
    """Return what the facts mean, a sentence each: the fragments, the overlaps, the phenotype."""
    joined = listed([f"{part.name} ({part.fragment_length} bp)" for part in parts])
    return (
        f"One reaction joins {len(parts)} fragments: {joined}.",
        *(
            f"{one.before} and {one.after} share {one.length} bp at {one.start}, taken from "
            f"{one.taken_from} and carried by the other as a primer tail."
            for one in assembly.junctions
        ),
        *phenotype_sentences(phenotype, inserts),
    )


def _materials(
    *,
    vector: SequenceRecord,
    parts: Sequence[Part],
    inserts: Sequence[str],
    product: AssemblyProduct,
    colony: ColonyCheck,
    host: str,
    polymerase: Polymerase,
    phenotype: Phenotype,
) -> tuple[Material, ...]:
    """Every reagent and consumable the protocol asks for. The oligos are the order sheet."""
    ladders = dict.fromkeys(
        (
            choose_ladder(tuple(part.length for part in parts if part.amplified)).name,
            colony.ladder.name,
        )
    )
    return (
        Material(
            f"{vector.name} plasmid" if parts[0].amplified else vector.name,
            storage="-20 °C",
            note="PCR template" if parts[0].amplified else "the opened backbone, used as given",
        ),
        *(Material(f"{name} template", storage="-20 °C", note="PCR template") for name in inserts),
        Material(
            f"{polymerase.name} DNA Polymerase and its reaction buffer",
            supplier=product.supplier,
            storage="-20 °C",
        ),
        Material("dNTP mix", storage="-20 °C", note=f"{DNTP_STOCK_MM:g} mM of each base"),
        *(
            (
                Material(
                    "DpnI",
                    supplier=product.supplier,
                    storage="-20 °C",
                    amount=f"{DPNI_UNITS} units per PCR",
                    note="cuts the methylated plasmid template only",
                ),
            )
            if any(part.dpni for part in parts)
            else ()
        ),
        Material("PCR and gel cleanup spin columns"),
        Material(
            product.name,
            supplier=product.supplier,
            catalog=product.catalog,
            storage="-20 °C",
            amount=f"{product.master_mix_ul:g} µL per reaction",
        ),
        catalogued(
            host,
            supplier=product.supplier if host == DEFAULT_HOST else "",
            storage="-80 °C",
            amount=f"{CELLS_UL:g} µL per transformation",
        ),
        Material(
            "SOC or NEB 10-beta/Stable Outgrowth Medium",
            amount=f"{OUTGROWTH_UL:g} µL per transformation",
        ),
        Material(_plate(phenotype), amount="one plate per transformation"),
        catalogued(
            COLONY_PCR_MASTER_MIX,
            supplier=product.supplier,
            storage="-20 °C",
            amount=f"{colony_pcr_master_mix_component().volume_ul:g} µL per reaction",
        ),
        Material("Agarose and 1X TAE or TBE"),
        *(catalogued(name, supplier=product.supplier) for name in ladders),
    )


def _purpose(oligo: DesignedOligo) -> str:
    """Return the title of the step that uses this oligo."""
    if oligo.part is not None:
        return pcr_title(oligo.part.name)
    return COLONY_PCR_TITLE if oligo.role == "colony PCR" else SEQUENCING_TITLE


def _plate(phenotype: Phenotype) -> str:
    """Return what to pour the selection plates with."""
    antibiotic = phenotype.antibiotic or "the vector's own antibiotic"
    if phenotype.blue_white:
        return (
            f"{phenotype.medium} agar plates with {antibiotic}, {XGAL_UG_ML} µg/mL X-gal "
            f"and {IPTG_UM} µM IPTG"
        )
    return f"{phenotype.medium} agar plates with {antibiotic}"


def _steps(
    *,
    parts: Sequence[Part],
    inserts: Sequence[str],
    assembly: Assembly,
    product: AssemblyProduct,
    amounts: tuple[Amount, ...],
    colony: ColonyCheck,
    reads: Sequence[SangerRead],
    phenotype: Phenotype,
    host: str,
    polymerase: Polymerase,
) -> tuple[Step, ...]:
    """Return the steps in the order they happen, the shared ones carrying this method's notes.

    A part handed in ready to assemble has no PCR and nothing to run on the gel, so neither
    step names it.
    """
    cut = [part for part in parts if part.dpni]
    made = [part for part in parts if part.amplified]
    steps = [_pcr_step(part, assembly, polymerase) for part in made]
    steps.append(gel_step([(part.name, part.length) for part in made]))
    if cut:
        steps.append(
            dpni_step(
                [part.name for part in cut],
                [(part.template.name, dam_sites(part.template)) for part in cut],
                seconds=DPNI_DIGEST_SECONDS,
                inactivation=DPNI_INACTIVATION,
                notes=(
                    "This digest is NEB's own, not this package's: its dose, its time and its "
                    "heat inactivation are the ones the assembly manual prints.",
                    "NEB's other answer to template background is to keep it low in the first "
                    "place: 0.1-0.5 ng of plasmid per 50 µL PCR.",
                ),
            )
        )
    steps.append(
        cleanup_step(
            notes=(
                f"NEB calls a column optional below {CLEANUP_FRAGMENTS} PCR fragments: a "
                "product more than 90% pure goes into the reaction unpurified, within the "
                "fraction the assembly step gives.",
                f"At {CLEANUP_FRAGMENTS} fragments or more, or a fragment over {CLEANUP_KB} kb, "
                "NEB calls it highly recommended and puts the gain in assembly and "
                "transformation at two- to tenfold. A PCR showing anything but one band is "
                "gel-purified whatever the count.",
            )
        )
    )
    steps.append(quantify_step(amounts))
    steps.append(_assembly_step(product, amounts))
    steps.append(_incubation_step(product, assembly.junctions, len(parts)))
    steps.append(
        transform_step(
            host,
            phenotype,
            inserts=inserts,
            colonies=f"Hundreds of colonies. NEB's own lot test asks for more than "
            f"{RELEASE_COLONIES} from a six-fragment assembly with a tenth of the outgrowth "
            f"plated; this reaction joins {len(parts)}.",
        )
    )
    steps.append(_colony_pcr_step(colony, assembly))
    steps.append(
        sequencing_step(
            reads,
            junctions=[one.overlap for one in assembly.junctions],
            inserts=inserts,
            notes=(
                f"Gibson 2009 sequenced 210 repaired junctions and found about one error per "
                f"{MOLECULES_PER_ERROR} molecules joined, so a gel that reads right is not the "
                "same as a junction that is right.",
            ),
        )
    )
    return tuple(steps)


def _pcr_step(part: Part, assembly: Assembly, polymerase: Polymerase) -> Step:
    """Amplify one part, with whatever overlaps its primers carry.

    Raises
    ------
    ValueError
        If the part was handed in ready to assemble, so no PCR makes it.
    """
    if part.report is None:
        raise ValueError(f"{part.name or 'this part'} is not amplified, so it has no PCR step")
    tails = [
        f"{len(tail)} bp of {_neighbour(assembly, part, before=first)} on the {end} primer"
        for tail, first, end in (
            (part.left_tail, True, "forward"),
            (part.right_tail, False, "reverse"),
        )
        if tail
    ]
    carried = (
        f"Its 5' tails are {listed(tails)}, so the amplicon already spells what its neighbours "
        "spell where they meet."
        if tails
        else "Its primers carry no tail: both its junctions take their overlap from this "
        "fragment, and the neighbours' primers carry it."
    )
    return pcr_step(
        part.name,
        part.template.name,
        part.length,
        polymerase=polymerase,
        annealing_temperature=part.report.annealing_temperature,
        extension_seconds=part.report.extension_seconds,
        notes=(
            carried,
            "The annealing temperature above is read from the annealing regions alone; a tail "
            "pairs with nothing on the template in the first cycles.",
        ),
    )


def _neighbour(assembly: Assembly, part: Part, *, before: bool) -> str:
    """Name the part whose bases this one's tail spells, or the tail's own length otherwise."""
    for one in assembly.junctions:
        if before and one.after == part.name:
            return one.taken_from
        if not before and one.before == part.name:
            return one.taken_from
    return "its neighbour"


def _assembly_step(product: AssemblyProduct, amounts: tuple[Amount, ...]) -> Step:
    """Set the one-tube assembly up, at the molar ratio this product asks for."""
    table = assembly_reaction(product, amounts)
    tier = product.tier(len(amounts))
    total = sum(component.volume_ul for component in table.components)
    ratio = (
        f"Each insert goes in at {tier.insert_ratio:g} times the vector's moles, and any insert "
        f"under {product.short_insert_bp} bp at {product.short_insert_ratio:g} times, which is "
        f"what {product.supplier} asks for at this fragment count. The table gives each in "
        "picomoles and in nanograms, so it can be pipetted at whatever concentration it is."
    )
    unpurified = (
        f"{product.supplier} takes unpurified PCR product straight from the tube for up to "
        f"{product.unpurified_fraction:.0%} of the reaction, which is {product.unpurified_ul:g} "
        "µL here, as long as the product is a single band."
        if product.unpurified_fraction is not None
        else f"{product.supplier} asks for purified PCR product and documents no allowance for "
        "unpurified DNA in the reaction."
    )
    return Step(
        f"Set up the {product.name} reaction",
        instructions=(
            "Thaw the master mix on ice and mix it well.",
            "Pipette the DNA into the tube first, then the master mix.",
            "Mix gently and spin down.",
        ),
        tables=(table,),
        expected=(f"A {total:g} µL reaction holding every fragment.",),
        notes=(
            ratio,
            "The volumes above assume each fragment is concentrated enough to carry its "
            "picomoles in the microlitre the table gives; make the difference up with water.",
            unpurified,
        ),
        troubleshooting=(
            Troubleshooting(
                "The DNA does not fit the reaction volume",
                "Concentrate the fragments, or scale the whole reaction up and add master mix "
                "in proportion.",
            ),
            Troubleshooting(
                "Colonies later carry the empty vector",
                "The backbone PCR carried its plasmid template through; digest it with DpnI "
                "again, or gel-purify the backbone.",
            ),
        ),
    )


def _incubation_step(
    product: AssemblyProduct, junctions: Sequence[Junction], fragments: int
) -> Step:
    """Run the one isothermal incubation, which is the whole reaction."""
    tier = product.tier(fragments)
    return Step(
        "Incubate the assembly",
        instructions=(
            "Put the tube in the thermocycler and run the program below.",
            "Put the reaction on ice afterwards, or keep it at -20 °C.",
        ),
        programs=(assembly_program(product, fragments=fragments),),
        timers=(Timer("Assembly", tier.incubation_seconds),),
        expected=(
            "Nothing visible. One exonuclease, one polymerase and one ligase work together at "
            "this one temperature; there is nothing to cycle.",
            *(
                f"The product spells {one.overlap} at {one.start}, where {one.before} meets "
                f"{one.after}."
                for one in junctions
            ),
        ),
        notes=(
            f"{tier.incubation_seconds // 60} minutes is what {product.supplier} asks for at "
            f"this fragment count. {product.incubation_note}",
            "Junction positions are 0-based, on the product.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Few colonies later",
                "Run the kit's positive control beside the assembly; it tells a bad master mix "
                "from a bad design.",
            ),
            Troubleshooting(
                "No colonies later",
                "Run the reaction on a gel: an efficient assembly shows the fragments gone and "
                "a product of the right size.",
            ),
        ),
    )


def _colony_pcr_step(colony: ColonyCheck, assembly: Assembly) -> Step:
    """Screen the colonies, with what the note measured of how many should read correct."""
    return colony_pcr_step(
        colony,
        junctions=len(assembly.junctions),
        notes=(
            f"Pick {SCREENED_COLONIES}. Takara's In-Fusion series is the only one that measures "
            f"how the correct fraction falls with fragment count: {SCREENED_COLONIES} of "
            f"{SCREENED_COLONIES} correct at two fragments, {CORRECT_AT_FIVE} of "
            f"{SCREENED_COLONIES} at five. That is In-Fusion's number and not this product's.",
            "Where nothing grew at all, NEB asks for the same PCR on the assembly reaction "
            "itself, with primers flanking the assembled product.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Every colony reads as empty vector",
                "NEB's answer is the template: a PCR-generated vector carries uncut plasmid "
                "through, so digest it with DpnI again or gel-purify the backbone.",
            ),
            Troubleshooting(
                "A colony gives a band of the wrong size",
                "The PCR that made a part was not a single band; gel-purify it and assemble "
                "again. NEB suggests NEB Stable Competent E. coli (#C3040) for an insert "
                "carrying repeats.",
            ),
        ),
    )


def _references(
    parts: Sequence[Part], product: AssemblyProduct, phenotype: Phenotype
) -> tuple[Reference, ...]:
    """Where the numbers come from."""
    items = [*product.references, *SCREENING_REFERENCES, *BENCH_REFERENCES]
    if any(part.dpni for part in parts):
        items.append(DPNI_REFERENCE)
    if phenotype.blue_white:
        items.append(PLATE_REFERENCE)
    return tuple(items)

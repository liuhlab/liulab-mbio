"""The bench protocol a plan writes, step by step.

Every number here is computed by `liulab_mbio.goldengate.plan` or by the modules it calls, and
every sentence about the phenotype is read off the product's own features. The constants below
are the bench choices that no table of NEB's covers; each says where it comes from.
"""

import re
from collections.abc import Sequence
from typing import TYPE_CHECKING

from liulab_mbio.enzymes import Enzyme
from liulab_mbio.goldengate.assembly import Part, dam_sites
from liulab_mbio.goldengate.bench import (
    COLONY_PCR_MASTER_MIX,
    DNA_VOLUME_UL,
    DNTP_STOCK_MM,
    GOLDEN_GATE_PCR_CYCLES,
    PRIMER_STOCK_UM,
    REFERENCES,
    agarose_percent,
    assembly_program,
    assembly_reaction,
    choose_ladder,
    colony_pcr_program,
    colony_pcr_reaction,
    golden_gate_temperature,
    heat_inactivation,
    pcr_program,
    pcr_reaction,
)
from liulab_mbio.goldengate.plan import DEFAULT_HOST
from liulab_mbio.protocol import (
    OVERVIEW_CHARS,
    Check,
    Gel,
    Lane,
    Material,
    Oligo,
    Protocol,
    ReactionTable,
    Reference,
    Step,
    Timer,
    Troubleshooting,
)

if TYPE_CHECKING:
    from liulab_mbio.goldengate.plan import Plan

#: The DpnI digest that takes the plasmid template away. NEB's Golden Gate pages prescribe no
#: such step -- `docs/research/golden-gate-assembly.md` §3 justifies it from REBASE's record of
#: DpnI as methyl-directed -- so these are this package's choices, not a supplier's table.
DPNI_UNITS = 20
DPNI_CELSIUS = 37.0
DPNI_SECONDS = 3600

#: Transformation and plating, from the NEB #E1601 and #E1602 manuals (§5 of the same note):
#: microlitres, degrees Celsius and seconds.
CELLS_UL = 50.0
THAW_SECONDS = 600
ASSEMBLY_UL = 2.0
ICE_SECONDS = 1800
HEAT_SHOCK_CELSIUS = 42.0
HEAT_SHOCK_SECONDS = 30
RECOVER_SECONDS = 300
OUTGROWTH_UL = 950.0
OUTGROWTH_CELSIUS = 37.0
OUTGROWTH_SECONDS = 3600
PLATE_UL = 50.0
PLATE_DILUTION = 5

#: The indicator plate, from Potapov et al. 2018's recipe (§5): micrograms per millilitre of
#: X-gal and micromolar IPTG.
XGAL_UG_ML = 80
IPTG_UM = 200

#: Correct colonies NEB counts from a single-insert assembly with 2.5 µL of the outgrowth
#: plated (E1601 technical note, §5). An order-of-magnitude check, not a target.
NEB_COLONIES = 687

#: Who sells the products this protocol names. Every catalogue number it prints comes out of an
#: enzyme record or a product name a supplier wrote; none is written here.
SUPPLIER = "New England Biolabs"

#: NEBridge Ligase Master Mix: as `assembly_reaction` names the line, and as NEB sells it.
LIGASE_MIX = "NEBridge Ligase Master Mix"
LIGASE_MIX_PRODUCT = f"{LIGASE_MIX} (M1100)"

#: The steps an oligo's row points at, written once so a row and its step cannot drift.
COLONY_STEP = "Screen colonies by PCR"
SEQUENCING_STEP = "Confirm the clone by sequencing"

#: The hardware a run needs, which no reagent table covers.
EQUIPMENT: tuple[str, ...] = (
    "Thermocycler with a heated lid",
    "Agarose gel rig and power supply",
    "Microcentrifuge",
    "Spectrophotometer or fluorometer",
    f"Heat block or water bath at {HEAT_SHOCK_CELSIUS:g} °C",
    f"Shaking incubator and a plate incubator at {OUTGROWTH_CELSIUS:g} °C",
)

#: A catalogue number at the end of a product name, such as ``"(M1100)"``. A letter and then a
#: digit, so a bracketed enzyme name is not read as one.
_CATALOG_RE = re.compile(r"^(?P<name>.*?)\s*\((?P<catalog>[A-Z]\d[\w./-]*)\)$")


def protocol(plan: "Plan") -> Protocol:
    """Return the bench protocol for `plan`, ready to render.

    The steps run in the order someone does them: one PCR per part, the gel that checks them,
    the DpnI digest and cleanup, quantification, the assembly, transformation and plating,
    colony PCR, and sequencing. Every step is per experiment, however many parts there are.
    """
    vector = plan.vector
    inserts = _listed(_insert_names(plan))
    return Protocol(
        f"Golden Gate assembly: {inserts} into {vector.name}",
        summary=(
            f"Open {vector.name} by PCR across {plan.span[0]}-{plan.span[1]}, amplify "
            f"{inserts} with {plan.enzyme.name} tails, join the {len(plan.parts)} fragments in "
            "one Golden Gate reaction, and confirm the clone by colony PCR and sequencing."
        ),
        overview=_overview(plan),
        highlights=_highlights(plan),
        checks=_checks(plan),
        materials=_materials(plan),
        oligos=_oligos(plan),
        equipment=EQUIPMENT,
        steps=_steps(plan),
        references=_references(plan),
    )


def _overview(plan: "Plan") -> dict[str, str]:
    """Return the facts to check before starting, each short enough to be a card."""
    fidelity = plan.overhangs.fidelity
    facts = {
        "Vector": f"{plan.vector.name}, {len(plan.vector)} bp",
        "Insert" if len(plan.inserts) == 1 else "Inserts": _insert_fact(plan),
        "Enzyme": f"{_label(plan.enzyme)} at {golden_gate_temperature(plan.enzyme):g} °C",
        "Fragments": f"{len(plan.parts)} in one reaction",
        "Overhangs": _brief(plan.overhangs.overhangs, "junctions"),
        "Fidelity": f"{fidelity.value:.0%}, {fidelity.label}",
        "Product": f"{plan.product.name}, {len(plan.product)} bp",
    }
    selection = _selection(plan)
    if selection:
        facts["Selection"] = selection
    return facts


def _insert_fact(plan: "Plan") -> str:
    """Return the insert card: the one insert with its length, or what the inserts are called."""
    names = _insert_names(plan)
    if len(names) == 1:
        return f"{names[0]}, {len(plan.inserts[0])} bp"
    return _brief(names, "inserts")


def _selection(plan: "Plan") -> str:
    """Return what to select transformants on, or nothing where the vector annotates no marker."""
    phenotype = plan.phenotype
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
    listed = _listed(items)
    return listed if len(listed) <= OVERVIEW_CHARS else f"{len(items)} {noun}"


def _highlights(plan: "Plan") -> tuple[str, ...]:
    """Return what the facts mean, a sentence each: the fragments, then the phenotype."""
    joined = _listed([f"{part.name} ({part.fragment_length} bp)" for part in plan.parts])
    return (
        f"One reaction joins {len(plan.parts)} fragments: {joined}.",
        *_phenotype_sentences(plan),
    )


def _phenotype_sentences(plan: "Plan") -> tuple[str, ...]:
    """Return what the product's own features say about the insert and about the plate."""
    phenotype = plan.phenotype
    coding = phenotype.coding.name if phenotype.coding is not None else _listed(_insert_names(plan))
    lines: list[str] = []
    if phenotype.promoter is not None:
        way = (
            f"reads on the same strand as {phenotype.promoter.name}, "
            f"{phenotype.gap_bp} bp downstream of it"
            if phenotype.driven
            else (
                f"reads on the opposite strand from {phenotype.promoter.name}, "
                f"{phenotype.gap_bp} bp away, so that promoter does not transcribe it"
            )
        )
        lines.append(f"{coding} {way}.")
    site = "is" if phenotype.ribosome_binding_site else "is no"
    lines.append(
        f"There {site} ribosome binding site annotated ahead of {coding}, so the clone "
        f"{'may make' if phenotype.expressed else 'is not expected to make'} its protein."
    )
    if phenotype.blue_white and phenotype.reporter is not None:
        lines.append(
            f"The insertion interrupts {phenotype.reporter.name}, so correct clones are white "
            "and empty vector is blue on X-gal and IPTG."
        )
    return tuple(lines)


def _checks(plan: "Plan") -> tuple[Check, ...]:
    """Return the plan's verdicts, one badge each, so a warning is seen and not read."""
    return tuple(Check(check.name, check.status, detail=check.detail) for check in plan.checks)


def _materials(plan: "Plan") -> tuple[Material, ...]:
    """Every reagent and consumable the protocol asks for. The oligos are `_oligos`."""
    enzyme = plan.enzyme
    assembly = assembly_reaction(enzyme, plan.amounts)
    ladders = dict.fromkeys(
        (choose_ladder(tuple(part.length for part in plan.parts)).name, plan.colony.ladder.name)
    )
    return (
        Material(f"{plan.vector.name} plasmid", storage="-20 °C", note="PCR template"),
        *(
            Material(f"{name} template", storage="-20 °C", note="PCR template")
            for name in _insert_names(plan)
        ),
        Material(
            f"{plan.polymerase.name} DNA Polymerase and its reaction buffer",
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
        Material("PCR and gel cleanup spin columns"),
        _catalogued(
            LIGASE_MIX_PRODUCT,
            supplier=SUPPLIER,
            storage="-20 °C",
            amount=_per_reaction(assembly, LIGASE_MIX),
        ),
        Material(
            enzyme.commercial_name or enzyme.name,
            supplier=enzyme.supplier or "",
            catalog=enzyme.catalog_number or "",
            storage="-20 °C",
            amount=_per_reaction(assembly, _label(enzyme)),
        ),
        _catalogued(
            plan.host,
            supplier=SUPPLIER if plan.host == DEFAULT_HOST else "",
            storage="-80 °C",
            amount=f"{CELLS_UL:g} µL per transformation",
        ),
        Material(
            "SOC or NEB 10-beta/Stable Outgrowth Medium",
            amount=f"{OUTGROWTH_UL:g} µL per transformation",
        ),
        Material(_plate(plan), amount="one plate per transformation"),
        _catalogued(
            COLONY_PCR_MASTER_MIX,
            supplier=SUPPLIER,
            storage="-20 °C",
            amount=_per_reaction(colony_pcr_reaction(), COLONY_PCR_MASTER_MIX),
        ),
        Material("Agarose and 1X TAE or TBE"),
        *(_catalogued(name, supplier=SUPPLIER) for name in ladders),
    )


def _catalogued(name: str, *, supplier: str = "", storage: str = "", amount: str = "") -> Material:
    """Return a material, taking the catalogue number out of a product name that carries one.

    A name carrying none leaves the cell empty; nothing here invents one.

    Examples
    --------
    >>> _catalogued("NEB 100 bp DNA Ladder (N3231)").catalog
    'N3231'
    >>> _catalogued("Agarose and 1X TAE or TBE").catalog
    ''
    """
    found = _CATALOG_RE.match(name)
    if found is None:
        return Material(name, supplier=supplier, storage=storage, amount=amount)
    return Material(
        found["name"],
        supplier=supplier,
        catalog=found["catalog"],
        storage=storage,
        amount=amount,
    )


def _per_reaction(table: ReactionTable, name: str) -> str:
    """Return what one reaction takes of a component, or nothing where no line names it."""
    for component in table.components:
        if component.name.startswith(name):
            return f"{component.volume_ul:g} µL per reaction"
    return ""


def _oligos(plan: "Plan") -> tuple[Oligo, ...]:
    """Every designed oligo, in the order the primer sheet lists them."""
    stock = f"{PRIMER_STOCK_UM:g} µM"
    return tuple(
        Oligo(
            report.primer.name,
            report.primer.sequence,
            purpose=purpose,
            tm_c=round(report["tm"].value, 1),
            stock=stock,
        )
        for report, purpose in zip(plan.reports, _purposes(plan), strict=True)
    )


def _purposes(plan: "Plan") -> tuple[str, ...]:
    """Which step uses each oligo, in the order `Plan.reports` lists them."""
    return (
        *(f"Amplify {part.name}" for part in plan.parts for _ in ("forward", "reverse")),
        *(COLONY_STEP for _ in plan.colony.reports),
        *(SEQUENCING_STEP for _ in plan.reads),
    )


def _plate(plan: "Plan") -> str:
    """Return what to pour the selection plates with."""
    antibiotic = plan.phenotype.antibiotic or "the vector's own antibiotic"
    if plan.phenotype.blue_white:
        return f"LB agar plates with {antibiotic}, {XGAL_UG_ML} µg/mL X-gal and {IPTG_UM} µM IPTG"
    return f"LB agar plates with {antibiotic}"


def _steps(plan: "Plan") -> tuple[Step, ...]:
    """Return the steps, in the order they happen."""
    steps = [_pcr_step(plan, part) for part in plan.parts]
    steps.append(_gel_step(plan))
    if any(part.dpni for part in plan.parts):
        steps.append(_dpni_step(plan))
    steps.append(_cleanup_step(plan))
    steps.append(_quantify_step(plan))
    steps.append(_assembly_step(plan))
    steps.append(_cycling_step(plan))
    steps.append(_transform_step(plan))
    steps.append(_colony_step(plan))
    steps.append(_sequencing_step(plan))
    return tuple(steps)


def _pcr_step(plan: "Plan", part: Part) -> Step:
    """Amplify one part with the tails that carry the enzyme site."""
    report = part.report
    return Step(
        f"Amplify {part.name}",
        instructions=(
            "Thaw the buffer, dNTPs and primers on ice, then vortex and spin them down.",
            f"Mix the master mix and put it in each tube, then add the "
            f"{part.template.name} template.",
            f"Run the program below: {report.annealing_temperature:g} °C annealing and "
            f"{report.extension_seconds} s extension for a {part.length} bp product.",
        ),
        cautions=("Keep the polymerase on ice.",),
        tables=(pcr_reaction(plan.polymerase, title=f"{part.name} PCR"),),
        programs=(
            pcr_program(
                plan.polymerase,
                annealing_temperature=report.annealing_temperature,
                amplicon_length=part.length,
                cycles=GOLDEN_GATE_PCR_CYCLES,
                title=f"{part.name} PCR",
            ),
        ),
        expected=(f"One band at {part.length} bp.",),
        notes=(
            f"{GOLDEN_GATE_PCR_CYCLES} cycles, which is the fewest NEB finds enough for an "
            "amplicon going into an assembly; fewer cycles means fewer PCR errors.",
            f"The primers carry a {plan.enzyme.name} site pointing back into the part, so "
            f"cutting the amplicon leaves {part.left_overhang} and {part.right_overhang}.",
        ),
        troubleshooting=(
            Troubleshooting(
                "No band",
                f"Drop the annealing temperature by 3 °C and check the "
                f"{part.template.name} template is there.",
            ),
            Troubleshooting(
                "Several bands",
                "Raise the annealing temperature, or gel-purify the band of the right size.",
            ),
        ),
    )


def _gel_step(plan: "Plan") -> Step:
    """Check every PCR before anything is spent on them."""
    sizes = tuple(part.length for part in plan.parts)
    percent = agarose_percent(sizes)
    return Step(
        "Check the PCRs on a gel",
        instructions=(
            f"Pour a {percent:g}% agarose gel.",
            "Load 5 µL of each reaction beside the ladder.",
            "Run until the dye front is two thirds down the gel.",
        ),
        gels=(
            Gel(
                choose_ladder(sizes),
                tuple(Lane(part.name, (part.length,)) for part in plan.parts),
                title="PCR products",
            ),
        ),
        expected=tuple(f"{part.name}: one band at {part.length} bp." for part in plan.parts),
        troubleshooting=(
            Troubleshooting(
                "A smear or an extra band",
                "Gel-purify the band of the right size; a wrong template in the assembly "
                "gives wrong clones.",
            ),
        ),
    )


def _dpni_step(plan: "Plan") -> Step:
    """Take the plasmid template away so it cannot transform as itself."""
    cut = [part for part in plan.parts if part.dpni]
    counted = ", ".join(
        f"{part.template.name} ({dam_sites(part.template)} Dam sites)" for part in cut
    )
    return Step(
        "Digest the plasmid template with DpnI",
        instructions=(
            *(f"Add {DPNI_UNITS} units of DpnI to the {part.name} PCR and mix." for part in cut),
            f"Incubate at {DPNI_CELSIUS:g} °C for {DPNI_SECONDS // 60} minutes.",
        ),
        timers=(Timer("DpnI digest", DPNI_SECONDS),),
        expected=(
            "Nothing visible. The digest shows up later as fewer colonies carrying the "
            "template plasmid.",
        ),
        notes=(
            f"DpnI cuts GATC only where Dam has methylated it, so it cuts {counted} and "
            "leaves the PCR product, which carries no methylation.",
            "This step is not in NEB's Golden Gate protocol; the incubation is this "
            "package's choice.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Many colonies on the no-insert control",
                "The template survived: digest longer, or use more DpnI.",
            ),
        ),
    )


def _cleanup_step(plan: "Plan") -> Step:
    """Purify every amplicon, which NEB asks for."""
    return Step(
        "Purify every amplicon",
        instructions=(
            "Run each reaction over a spin column and elute in the smallest volume the kit allows.",
        ),
        expected=("Clean DNA, free of polymerase, primers and dNTPs.",),
        notes=(
            "NEB asks for purified amplicons: polymerase carried over from the PCR fills in "
            "the four-base overhangs, which blunts the ends and mis-assembles them.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Low recovery",
                "Elute twice through the same column, or pool two reactions before purifying.",
            ),
        ),
    )


def _quantify_step(plan: "Plan") -> Step:
    """Measure what the reaction is about to take."""
    wanted = tuple(
        f"{amount.name}: {amount.pmol:g} pmol is {amount.nanograms:g} ng, so "
        f"{amount.nanograms / DNA_VOLUME_UL:.0f} ng/µL or more fits in {DNA_VOLUME_UL:g} µL."
        for amount in plan.amounts
    )
    return Step(
        "Measure every concentration",
        instructions=(
            "Measure each purified amplicon by A260 or with a fluorometer.",
            "Work out the volume that carries the picomoles the next table asks for.",
        ),
        expected=wanted,
        notes=(
            "Picomoles, not nanograms: the shorter fragment weighs less at the same molar "
            "ratio. Mass to moles here is NEBioCalculator's 36.04 + 615.94 per base pair, "
            "which is about 5% off the 650 Da per base pair of NEB's manuals.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Too dilute to fit in the reaction",
                "Concentrate the amplicon, or scale the reaction up.",
            ),
        ),
    )


def _assembly_step(plan: "Plan") -> Step:
    """Set the one-tube digest and ligation up."""
    table = assembly_reaction(plan.enzyme, plan.amounts)
    total = sum(component.volume_ul for component in table.components)
    return Step(
        "Set up the Golden Gate reaction",
        instructions=(
            "Thaw the master mix on ice and mix it well; it is viscous.",
            "Pipette the DNA into the tube first, then the rest.",
            "Mix gently and spin down.",
        ),
        tables=(table,),
        expected=(f"A {total:g} µL reaction holding every fragment.",),
        notes=(
            "The volumes above assume the concentrations measured in the step before; "
            "pipette the volume that gives the picomoles, and make the difference up with "
            "water.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The DNA does not fit the reaction volume",
                "Concentrate the fragments, or scale the whole reaction up.",
            ),
        ),
    )


def _cycling_step(plan: "Plan") -> Step:
    """Run it, and the heat inactivation the supplier gives."""
    programs = [assembly_program(plan.enzyme, fragments=len(plan.parts))]
    kill = heat_inactivation(plan.enzyme)
    if kill is not None:
        programs.append(kill)
    expected = [
        "Nothing visible. The 60 °C soak at the end is a digest, not heat inactivation: it "
        "cuts vector that never opened or has closed again, so fewer empty colonies grow.",
        *(
            f"The product carries {one.overhang} at {one.start}, joining {one.before} to "
            f"{one.after}."
            for one in plan.assembly.junctions
        ),
    ]
    return Step(
        "Run the Golden Gate program",
        instructions=("Put the tube in the thermocycler and run the program below.",),
        programs=tuple(programs),
        expected=tuple(expected),
        notes=("Junction positions are 0-based, on the product.",),
        troubleshooting=(
            Troubleshooting(
                "Mostly empty vector later",
                "Keep the 60 °C soak, and check the template was digested with DpnI.",
            ),
            Troubleshooting(
                "Few colonies later",
                "Raise the cycle count, or plate more of the outgrowth.",
            ),
        ),
    )


def _transform_step(plan: "Plan") -> Step:
    """Transform and plate, with the colour the plate should show."""
    phenotype = plan.phenotype
    expected = [
        f"Hundreds of colonies; NEB counts about {NEB_COLONIES} correct ones from a "
        "single-insert assembly with 2.5 µL of the outgrowth plated.",
    ]
    if phenotype.blue_white and phenotype.reporter is not None:
        expected.append(
            f"Correct clones are white and empty vector is blue: the insertion interrupts "
            f"{phenotype.reporter.name}, which is then not there to complete the host's own."
        )
    notes = [
        f"The plate reads colour only with an alpha-complementing host, such as {plan.host}. "
        "A host that cannot complement gives white colonies whatever the clone carries."
        if phenotype.blue_white
        else "Colour does not report this insertion; screen every colony by PCR.",
    ]
    if not phenotype.expressed:
        notes.append(_expression_note(plan))
    return Step(
        "Transform and plate",
        instructions=(
            f"Thaw {CELLS_UL:g} µL of {plan.host} on ice for {THAW_SECONDS // 60} minutes.",
            f"Add {ASSEMBLY_UL:g} µL of the assembly and flick the tube four or five times.",
            f"Hold on ice for {ICE_SECONDS // 60} minutes.",
            f"Heat shock at {HEAT_SHOCK_CELSIUS:g} °C for {HEAT_SHOCK_SECONDS} seconds.",
            f"Return to ice for {RECOVER_SECONDS // 60} minutes.",
            f"Add {OUTGROWTH_UL:g} µL of outgrowth medium and shake at "
            f"{OUTGROWTH_CELSIUS:g} °C for {OUTGROWTH_SECONDS // 60} minutes at 250 rpm.",
            f"Spread {PLATE_UL:g} µL of a 1:{PLATE_DILUTION} dilution on a warmed plate and "
            "grow overnight at 37 °C.",
        ),
        cautions=("Competent cells die if they warm up; keep them on ice until the shock.",),
        timers=(
            Timer("On ice", ICE_SECONDS),
            Timer("Heat shock", HEAT_SHOCK_SECONDS),
            Timer("Outgrowth", OUTGROWTH_SECONDS),
        ),
        expected=tuple(expected),
        notes=tuple(notes),
        troubleshooting=(
            Troubleshooting(
                "No colonies",
                "Check the antibiotic and the cells' efficiency, and plate the rest of the "
                "outgrowth.",
            ),
            Troubleshooting(
                "A lawn",
                "Plate a smaller volume or a greater dilution next time.",
            ),
        ),
    )


def _expression_note(plan: "Plan") -> str:
    """One sentence on whether the clone should make the insert's protein."""
    phenotype = plan.phenotype
    coding = phenotype.coding.name if phenotype.coding is not None else _listed(_insert_names(plan))
    reasons = []
    if phenotype.promoter is not None and not phenotype.driven:
        reasons.append(f"it reads on the opposite strand from {phenotype.promoter.name}")
    if not phenotype.ribosome_binding_site:
        reasons.append("no ribosome binding site is annotated ahead of it")
    because = " and ".join(reasons) if reasons else "the product annotates nothing to drive it"
    return f"The clone is not expected to make {coding}: {because}."


def _colony_step(plan: "Plan") -> Step:
    """Screen the colonies, and say which band means what."""
    check = plan.colony
    sizes = tuple(bp for clone in check.clones for bp in clone.bands_bp)
    expected = [
        f"Each of the {len(plan.assembly.junctions)} junctions is read: the flanking pair "
        "crosses them all, and each junction primer stops inside its own insert.",
        *(
            f"{clone.name}: {', '.join(f'{bp} bp' for bp in clone.bands_bp) or 'no band'}."
            for clone in check.clones
        ),
    ]
    expected.append(
        "A reversed insert is told from a correct one, because the two vector primers sit at "
        "different distances from their own junctions."
        if check.tells_orientation
        else "These primers cannot tell a reversed insert from a correct one."
    )
    return Step(
        COLONY_STEP,
        instructions=(
            "Touch a well-separated colony with a sterile toothpick and stir it into the "
            "tube until the liquid clouds.",
            "Streak the same toothpick onto a numbered plate, so the clone survives the PCR.",
            f"Run the program below and load 5 µL on a {check.agarose_percent:g}% gel.",
        ),
        tables=(colony_pcr_reaction(),),
        programs=(
            colony_pcr_program(
                annealing_temperature=check.annealing_temperature,
                amplicon_length=max(sizes),
            ),
        ),
        gels=(check.gel,),
        expected=tuple(expected),
        notes=("The long first step at 94 °C lyses the cells; there is no purified template.",),
        troubleshooting=(
            Troubleshooting(
                "No band in any lane",
                "The colony was too much material: touch a smaller one, or dilute it.",
            ),
            Troubleshooting(
                "Every colony reads as empty vector",
                "The template survived the DpnI digest, or the vector re-closed; check the "
                "60 °C soak ran.",
            ),
        ),
    )


def _sequencing_step(plan: "Plan") -> Step:
    """Confirm the junctions, which is the only thing that settles it."""
    reads = tuple(
        f"{read.primer.name} anneals {read.distance_bp} bp from its own junction and has to "
        f"read {read.read_bp} bp to cover the far one."
        for read in plan.reads
    )
    return Step(
        SEQUENCING_STEP,
        instructions=(
            "Miniprep two or three colonies that read as correct.",
            "Send each with both sequencing primers.",
            "Check the read across every junction and the whole of each insert.",
        ),
        expected=(
            *reads,
            f"The junctions read as {_listed(plan.overhangs.overhangs)}, and the parts match "
            f"{_listed(_insert_names(plan))}.",
        ),
        notes=(
            "NEB asks for the assembly to be confirmed by sequencing across the junctions "
            "whatever the screen said.",
            "A provider whose read is shorter than the lengths above needs a further primer "
            "inside the inserts.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The read starts too close to the junction",
                "Move the primer further out; the first bases after a primer are unreadable.",
            ),
        ),
    )


def _references(plan: "Plan") -> tuple[Reference, ...]:
    """Where the numbers come from."""
    items = list(REFERENCES)
    items.append(Reference(plan.overhangs.fidelity.source))
    if any(part.dpni for part in plan.parts):
        items.append(
            Reference(
                "REBASE record for DpnI, which cuts G6mATC and so only methylated template",
                url="http://rebase.neb.com/rebase/rebase.html",
            )
        )
    if plan.phenotype.blue_white:
        items.append(
            Reference(
                "Potapov, V. et al. (2018) Comprehensive profiling of four base overhang "
                "ligation fidelity by T4 DNA Ligase and application to DNA assembly. ACS "
                "Synth. Biol. 7, 2665-2674, for the X-gal and IPTG plate",
                url="https://doi.org/10.1021/acssynbio.8b00333",
            )
        )
    return tuple(items)


def _insert_names(plan: "Plan") -> tuple[str, ...]:
    """Return what each insert is called, which is what its part and its tube are labelled."""
    return tuple(part.name for part in plan.parts[1:])


def _listed(items: Sequence[str]) -> str:
    """Join names the way a sentence does, with `and` before the last."""
    if len(items) < 3:
        return " and ".join(items)
    return f"{', '.join(items[:-1])} and {items[-1]}"


def _label(enzyme: Enzyme) -> str:
    """Return the enzyme as a supplier sells it."""
    name = enzyme.commercial_name or enzyme.name
    return f"{name} ({enzyme.catalog_number})" if enzyme.catalog_number else name

"""The protocol steps any cloning pipeline reuses, each built from plain facts.

A pipeline runs these in its own order around its own steps, and passes its own notes where it
has something of its own to say; they follow the step's. Every sentence about the phenotype is
read off `liulab_mbio.bench.phenotype`. The smaller pieces every pipeline shapes the same way --
a verdict's badge, an overview card and a material row for an enzyme or a catalogued product --
are here too.
"""

import re
from collections.abc import Sequence
from dataclasses import KW_ONLY, dataclass

from liulab_mbio import checks as judged
from liulab_mbio.bench.amounts import DNA_VOLUME_UL, Amount
from liulab_mbio.bench.gels import agarose_percent, choose_ladder
from liulab_mbio.bench.pcr import colony_pcr_program, colony_pcr_reaction, pcr_program, pcr_reaction
from liulab_mbio.bench.phenotype import Phenotype
from liulab_mbio.bench.validation import ColonyCheck, SangerRead
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.primers.polymerase import Polymerase
from liulab_mbio.protocol.model import (
    OVERVIEW_CHARS,
    Check,
    Gel,
    Incubation,
    Lane,
    Material,
    Reference,
    Step,
    Timer,
    Troubleshooting,
)

#: The DpnI digest that takes the plasmid template away. No supplier's table sets these, so they
#: are this package's choices; `docs/research/golden-gate-assembly.md` §3 justifies the digest
#: from REBASE's record of DpnI as methyl-directed.
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

#: The titles of the steps an oligo's row points at, written once so a row and its step cannot
#: drift. `pcr_title` gives the third.
COLONY_PCR_TITLE = "Screen colonies by PCR"
SEQUENCING_TITLE = "Confirm the clone by sequencing"

#: What a protocol cites when it runs the DpnI digest.
DPNI_REFERENCE = Reference(
    "REBASE record for DpnI, which cuts G6mATC and so only methylated template",
    url="http://rebase.neb.com/rebase/rebase.html",
)

#: What a protocol cites when it pours the indicator plate.
PLATE_REFERENCE = Reference(
    "Potapov, V. et al. (2018) Comprehensive profiling of four base overhang ligation fidelity by "
    "T4 DNA Ligase and application to DNA assembly. ACS Synth. Biol. 7, 2665-2674, for the X-gal "
    "and IPTG plate",
    url="https://doi.org/10.1021/acssynbio.8b00333",
)


@dataclass(frozen=True, slots=True)
class Transformation:
    """One supplier's heat-shock protocol, which no two of them state the same way.

    A method brings the numbers its own kit's manual gives; `NEB_TRANSFORMATION` is the one
    `transform_step` uses where a method names none.

    Parameters
    ----------
    cells_ul, reaction_ul
        Competent cells per tube, and how much of the reaction goes into them.
    source
        What the step calls that reaction, such as ``"the assembly"``.
    thaw_seconds, ice_seconds, heat_shock_celsius, heat_shock_seconds, recover_seconds
        Thawing the cells, and the shock itself. `thaw_seconds` is ``None`` where the kit's
        manual states no time, and the step then asks for a thaw without one.
    outgrowth_ul, outgrowth_celsius, outgrowth_seconds
        The medium and the recovery.
    plate_ul, dilution
        What is spread on one plate, and the dilution it is spread from.
    """

    _: KW_ONLY
    cells_ul: float
    reaction_ul: float
    source: str
    thaw_seconds: int | None
    ice_seconds: int
    heat_shock_celsius: float
    heat_shock_seconds: int
    recover_seconds: int
    outgrowth_ul: float
    outgrowth_celsius: float
    outgrowth_seconds: int
    plate_ul: float
    dilution: int


#: What `transform_step` asks for where a method names no protocol of its own: the constants
#: above, which are NEB's.
NEB_TRANSFORMATION = Transformation(
    cells_ul=CELLS_UL,
    reaction_ul=ASSEMBLY_UL,
    source="the assembly",
    thaw_seconds=THAW_SECONDS,
    ice_seconds=ICE_SECONDS,
    heat_shock_celsius=HEAT_SHOCK_CELSIUS,
    heat_shock_seconds=HEAT_SHOCK_SECONDS,
    recover_seconds=RECOVER_SECONDS,
    outgrowth_ul=OUTGROWTH_UL,
    outgrowth_celsius=OUTGROWTH_CELSIUS,
    outgrowth_seconds=OUTGROWTH_SECONDS,
    plate_ul=PLATE_UL,
    dilution=PLATE_DILUTION,
)


def listed(items: Sequence[str]) -> str:
    """Join names the way a sentence does, with `and` before the last.

    Examples
    --------
    >>> listed(("GFP", "Linker", "Tag"))
    'GFP, Linker and Tag'
    """
    if len(items) < 3:
        return " and ".join(items)
    return f"{', '.join(items[:-1])} and {items[-1]}"


def badges(checks: Sequence[judged.Check]) -> tuple[Check, ...]:
    """Return a plan's verdicts, one badge each, so a warning is seen and not read.

    A check no sourced threshold judges keeps its badge and shows no verdict on it, which is
    what stops it being read as a pass.
    """
    return tuple(Check(check.name, check.status, detail=check.detail) for check in checks)


def card(value: str, short: str) -> str:
    """Return the fact where an overview card holds it, and `short` where it does not.

    Examples
    --------
    >>> card("2 positions", "shorter")
    '2 positions'
    """
    return value if len(value) <= OVERVIEW_CHARS else short


def enzyme_material(enzyme: Enzyme, *, amount: str = "", note: str = "") -> Material:
    """Return one enzyme as a material, its own record carrying the supplier and catalogue number.

    `amount` is what one reaction takes of it, and `note` what it is there to cut.
    """
    return Material(
        enzyme.commercial_name or enzyme.name,
        supplier=enzyme.supplier or "",
        catalog=enzyme.catalog_number or "",
        storage="-20 °C",
        amount=amount,
        note=note,
    )


#: A catalogue number at the end of a product name, such as ``"(M1100)"``. A letter and then a
#: digit, so a bracketed enzyme name is not read as one.
_CATALOG_RE = re.compile(r"^(?P<name>.*?)\s*\((?P<catalog>[A-Z]\d[\w./-]*)\)$")


def catalogued(name: str, *, supplier: str = "", storage: str = "", amount: str = "") -> Material:
    """Return one product as a material, taking the catalogue number out of a name carrying one.

    A name carrying none leaves the cell empty; nothing here invents one.

    Examples
    --------
    >>> catalogued("NEB 100 bp DNA Ladder (N3231)").catalog
    'N3231'
    >>> catalogued("Agarose and 1X TAE or TBE").catalog
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


def pcr_title(name: str) -> str:
    """Return the title of the step that amplifies `name`, which its oligos name as their purpose."""
    return f"Amplify {name}"


def pcr_step(
    name: str,
    template: str,
    length_bp: int,
    *,
    polymerase: Polymerase,
    annealing_temperature: float,
    extension_seconds: int | None,
    cycles: int | None = None,
    notes: Sequence[str] = (),
) -> Step:
    """Return the step that makes one amplicon by PCR.

    Parameters
    ----------
    name
        What the amplicon is called, which titles the step, its reaction and its program.
    template
        What its template is called.
    length_bp
        The amplicon's length.
    polymerase, annealing_temperature, extension_seconds
        The PCR's, as the primer pair's report gives them.
    cycles
        Replaces the polymerase profile's own count.
    notes
        The caller's own.
    """
    return Step(
        pcr_title(name),
        instructions=(
            "Thaw the buffer, dNTPs and primers on ice, then vortex and spin them down.",
            f"Mix the master mix and put it in each tube, then add the {template} template.",
            f"Run the program below: {annealing_temperature:g} °C annealing and "
            f"{extension_seconds} s extension for a {length_bp} bp product.",
        ),
        cautions=("Keep the polymerase on ice.",),
        tables=(pcr_reaction(polymerase, title=f"{name} PCR"),),
        programs=(
            pcr_program(
                polymerase,
                annealing_temperature=annealing_temperature,
                amplicon_length=length_bp,
                cycles=cycles,
                title=f"{name} PCR",
            ),
        ),
        expected=(f"One band at {length_bp} bp.",),
        notes=tuple(notes),
        troubleshooting=(
            Troubleshooting(
                "No band",
                f"Drop the annealing temperature by 3 °C and check the {template} template is "
                "there.",
            ),
            Troubleshooting(
                "Several bands",
                "Raise the annealing temperature, or gel-purify the band of the right size.",
            ),
        ),
    )


def gel_step(amplicons: Sequence[tuple[str, int]]) -> Step:
    """Return the gel that checks every PCR before anything is spent on them.

    `amplicons` gives each amplicon's name and its length in base pairs, one lane each.
    """
    sizes = tuple(length_bp for _, length_bp in amplicons)
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
                tuple(Lane(name, (length_bp,)) for name, length_bp in amplicons),
                title="PCR products",
            ),
        ),
        expected=tuple(f"{name}: one band at {length_bp} bp." for name, length_bp in amplicons),
        troubleshooting=(
            Troubleshooting(
                "A smear or an extra band",
                "Gel-purify the band of the right size; a wrong template in the assembly gives "
                "wrong clones.",
            ),
        ),
    )


def dpni_step(
    pcrs: Sequence[str],
    templates: Sequence[tuple[str, int]],
    *,
    seconds: int = DPNI_SECONDS,
    inactivation: Incubation | None = None,
    notes: Sequence[str] = (),
) -> Step:
    """Return the DpnI digest that takes the plasmid template away, so it cannot transform.

    Parameters
    ----------
    pcrs
        The PCRs to digest, by name.
    templates
        The plasmid each was amplified from, and the Dam sites it carries.
    seconds
        How long to digest. `DPNI_SECONDS` is this package's own choice; a method whose
        supplier prescribes the digest passes that supplier's time instead.
    inactivation
        The heat inactivation the supplier asks for after it, where one is prescribed.
    notes
        The caller's own, after the step's.
    """
    counted = ", ".join(f"{name} ({sites} Dam sites)" for name, sites in templates)
    kill = (
        ()
        if inactivation is None
        else (
            f"{inactivation.label} at {inactivation.temperature_c:g} °C for "
            f"{(inactivation.seconds or 0) // 60:g} minutes.",
        )
    )
    return Step(
        "Digest the plasmid template with DpnI",
        instructions=(
            *(f"Add {DPNI_UNITS} units of DpnI to the {name} PCR and mix." for name in pcrs),
            f"Incubate at {DPNI_CELSIUS:g} °C for {seconds // 60} minutes.",
            *kill,
        ),
        timers=(
            Timer("DpnI digest", seconds),
            *(
                ()
                if inactivation is None or inactivation.seconds is None
                else (Timer(inactivation.label, inactivation.seconds),)
            ),
        ),
        expected=(
            "Nothing visible. The digest shows up later as fewer colonies carrying the "
            "template plasmid.",
        ),
        notes=(
            f"DpnI cuts GATC only where Dam has methylated it, so it cuts {counted} and "
            "leaves the PCR product, which carries no methylation.",
            *notes,
        ),
        troubleshooting=(
            Troubleshooting(
                "Many colonies on the no-insert control",
                "The template survived: digest longer, or use more DpnI.",
            ),
        ),
    )


def cleanup_step(*, notes: Sequence[str] = ()) -> Step:
    """Return the spin-column cleanup of every amplicon, carrying the caller's own notes."""
    return Step(
        "Purify every amplicon",
        instructions=(
            "Run each reaction over a spin column and elute in the smallest volume the kit allows.",
        ),
        expected=("Clean DNA, free of polymerase, primers and dNTPs.",),
        notes=tuple(notes),
        troubleshooting=(
            Troubleshooting(
                "Low recovery",
                "Elute twice through the same column, or pool two reactions before purifying.",
            ),
        ),
    )


def quantify_step(amounts: Sequence[Amount]) -> Step:
    """Return the step that measures what the next reaction is about to take."""
    wanted = tuple(
        f"{amount.name}: {amount.pmol:g} pmol is {amount.nanograms:g} ng, so "
        f"{amount.nanograms / DNA_VOLUME_UL:.0f} ng/µL or more fits in {DNA_VOLUME_UL:g} µL."
        for amount in amounts
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


def transform_step(
    host: str,
    phenotype: Phenotype,
    *,
    inserts: Sequence[str],
    colonies: str,
    protocol: Transformation = NEB_TRANSFORMATION,
    title: str = "Transform and plate",
    expected: Sequence[str] = (),
    notes: Sequence[str] = (),
) -> Step:
    """Return the transformation and plating, with the colour the plate should show.

    Parameters
    ----------
    host
        The competent strain.
    phenotype
        What the product says about itself.
    inserts
        What the inserts are called, for a product that annotates no coding sequence among them.
    colonies
        How many colonies to expect, as a sentence: a count belongs to the pipeline's reaction.
    protocol
        The volumes and times to run it by, which are the kit manufacturer's.
    title
        The step's, for a method that transforms more than once and has to tell them apart.
    expected, notes
        The caller's own, after the step's.
    """
    results = [colonies]
    if phenotype.blue_white and phenotype.reporter is not None:
        results.append(
            f"Correct clones are white and empty vector is blue: the insertion interrupts "
            f"{phenotype.reporter.name}, which is then not there to complete the host's own."
        )
    results.extend(expected)
    said = [
        f"The plate reads colour only with an alpha-complementing host, such as {host}. "
        "A host that cannot complement gives white colonies whatever the clone carries."
        if phenotype.blue_white
        else "Colour does not report this insertion; screen every colony by PCR.",
    ]
    if not phenotype.expressed:
        said.append(_expression_note(phenotype, inserts))
    said.extend(notes)
    return Step(
        title,
        instructions=(
            f"Thaw {protocol.cells_ul:g} µL of {host} on ice"
            + (
                "."
                if protocol.thaw_seconds is None
                else f" for {protocol.thaw_seconds // 60} minutes."
            ),
            f"Add {protocol.reaction_ul:g} µL of {protocol.source} and flick the tube four or "
            "five times.",
            f"Hold on ice for {protocol.ice_seconds // 60} minutes.",
            f"Heat shock at {protocol.heat_shock_celsius:g} °C for "
            f"{protocol.heat_shock_seconds} seconds.",
            f"Return to ice for {protocol.recover_seconds // 60} minutes.",
            f"Add {protocol.outgrowth_ul:g} µL of outgrowth medium and shake at "
            f"{protocol.outgrowth_celsius:g} °C for {protocol.outgrowth_seconds // 60} minutes "
            "at 250 rpm.",
            f"Spread {protocol.plate_ul:g} µL of a 1:{protocol.dilution} dilution on a warmed "
            "plate and grow overnight at 37 °C.",
        ),
        cautions=("Competent cells die if they warm up; keep them on ice until the shock.",),
        timers=(
            Timer("On ice", protocol.ice_seconds),
            Timer("Heat shock", protocol.heat_shock_seconds),
            Timer("Outgrowth", protocol.outgrowth_seconds),
        ),
        expected=tuple(results),
        notes=tuple(said),
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


def colony_pcr_step(
    check: ColonyCheck,
    *,
    junctions: int,
    notes: Sequence[str] = (),
    troubleshooting: Sequence[Troubleshooting] = (),
) -> Step:
    """Return the colony PCR screen, saying which band means what.

    Parameters
    ----------
    check
        The colony PCR, its expected clones included.
    junctions
        How many junctions it reads.
    notes, troubleshooting
        The caller's own, after the step's. How many colonies read correct is one of them: it
        is the method's own measurement and not this builder's.
    """
    sizes = tuple(bp for clone in check.clones for bp in clone.bands_bp)
    expected = [
        f"Each of the {junctions} junctions is read: the flanking pair crosses them all, and "
        "each junction primer stops inside its own insert.",
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
        COLONY_PCR_TITLE,
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
        notes=(
            "The long first step at 94 °C lyses the cells; there is no purified template.",
            *notes,
        ),
        troubleshooting=(
            Troubleshooting(
                "No band in any lane",
                "The colony was too much material: touch a smaller one, or dilute it.",
            ),
            *troubleshooting,
        ),
    )


def sequencing_step(
    reads: Sequence[SangerRead],
    *,
    junctions: Sequence[str],
    inserts: Sequence[str],
    notes: Sequence[str] = (),
) -> Step:
    """Return the sequencing that confirms the junctions, which is the only thing that settles it.

    Parameters
    ----------
    reads
        The sequencing primers and the reads they have to give.
    junctions
        The bases each junction spells.
    inserts
        What the inserts are called.
    notes
        The caller's own, after the step's, such as how often its method misjoins a junction.
    """
    lengths = tuple(
        f"{read.primer.name} anneals {read.distance_bp} bp from its own junction and has to "
        f"read {read.read_bp} bp to cover the far one."
        for read in reads
    )
    return Step(
        SEQUENCING_TITLE,
        instructions=(
            "Miniprep two or three colonies that read as correct.",
            "Send each with both sequencing primers.",
            "Check the read across every junction and the whole of each insert.",
        ),
        expected=(
            *lengths,
            f"The junctions read as {listed(junctions)}, and the parts match {listed(inserts)}.",
        ),
        notes=(
            "NEB asks for the assembly to be confirmed by sequencing across the junctions "
            "whatever the screen said.",
            "A provider whose read is shorter than the lengths above needs a further primer "
            "inside the inserts.",
            *notes,
        ),
        troubleshooting=(
            Troubleshooting(
                "The read starts too close to the junction",
                "Move the primer further out; the first bases after a primer are unreadable.",
            ),
        ),
    )


def phenotype_sentences(phenotype: Phenotype, inserts: Sequence[str]) -> tuple[str, ...]:
    """Return what the product's own features say about the insert and about the plate.

    `inserts` names the inserts for a product that annotates no coding sequence among them.
    """
    coding = _coding_name(phenotype, inserts)
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


def _expression_note(phenotype: Phenotype, inserts: Sequence[str]) -> str:
    """One sentence on whether the clone should make the insert's protein."""
    reasons = []
    if phenotype.promoter is not None and not phenotype.driven:
        reasons.append(f"it reads on the opposite strand from {phenotype.promoter.name}")
    if not phenotype.ribosome_binding_site:
        reasons.append("no ribosome binding site is annotated ahead of it")
    because = " and ".join(reasons) if reasons else "the product annotates nothing to drive it"
    return f"The clone is not expected to make {_coding_name(phenotype, inserts)}: {because}."


def _coding_name(phenotype: Phenotype, inserts: Sequence[str]) -> str:
    """Name the coding sequence the insert carries, or the inserts where none is annotated."""
    return phenotype.coding.name if phenotype.coding is not None else listed(inserts)

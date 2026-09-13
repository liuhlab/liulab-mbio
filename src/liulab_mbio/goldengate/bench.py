"""What a Golden Gate experiment takes at the bench: quantities, reactions and programs.

Every number is NEB's, through ``docs/research/golden-gate-assembly.md`` for the assembly and
``docs/research/primer-design-and-pcr.md`` for the conversions, the PCR programs and the
ladders. Functions return `liulab_mbio.protocol` values, so a protocol prints them unchanged.

NEB ships two Golden Gate systems and their tables do not mix. `LIGASE_MASTER_MIX` is NEBridge
Ligase Master Mix (M1100), which takes any NEB Type IIS enzyme; `KIT` is one of the kits, which
carry their own enzyme mix and so exist only for BsaI-HFv2 and BsmBI-v2.
"""

import dataclasses
from collections.abc import Mapping, Sequence
from dataclasses import KW_ONLY, dataclass
from itertools import pairwise
from typing import Literal

from liulab_mbio import edits
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.primers import (
    ONETAQ,
    Q5,
    THRESHOLDS,
    Polymerase,
    PrimerReport,
    Thresholds,
    amplicon_sizes,
    design_pair,
    design_primer,
    evaluate_primer,
)
from liulab_mbio.protocol import (
    Component,
    Gel,
    Incubation,
    Ladder,
    Lane,
    ReactionTable,
    Reference,
    Stage,
    ThermocyclerProgram,
)
from liulab_mbio.sequence import Primer, SequenceRecord, Strand, reverse_complement

#: NEBioCalculator's double-stranded DNA weight, g/mol: `_DUPLEX_ENDS + bp * _BASE_PAIR`. NEB's
#: manuals use 650 Da per base pair instead, which differs by about 5%, so a protocol says which.
_DUPLEX_ENDS = 36.04
_BASE_PAIR = 615.94

#: Picomoles of each fragment NEB's Golden Gate reactions take, whichever system is used.
FRAGMENT_PMOL = 0.05

#: Insert to vector molar ratio, from the E1601 manual revision 5.0_6/26 for amplicon inserts.
#: An earlier revision of the same page asked for 2:1.
INSERT_RATIO = 1.0

#: What to pipette of a fragment whose concentration is not known yet, µL.
DNA_VOLUME_UL = 1.0


def molecular_weight(length_bp: int) -> float:
    """Return the weight of a double-stranded DNA molecule, g/mol."""
    return _DUPLEX_ENDS + length_bp * _BASE_PAIR


def to_pmol(nanograms: float, length_bp: int) -> float:
    """Return how many picomoles a mass of double-stranded DNA of this length is.

    Examples
    --------
    >>> round(to_pmol(1000, 2686), 3)
    0.604
    """
    return 1000.0 * nanograms / molecular_weight(length_bp)


def to_nanograms(pmol: float, length_bp: int) -> float:
    """Return what those picomoles of double-stranded DNA of this length weigh, ng."""
    return pmol * molecular_weight(length_bp) / 1000.0


@dataclass(frozen=True, slots=True)
class Fragment:
    """One piece of DNA going into an assembly.

    Parameters
    ----------
    name
        What its tube is labelled.
    length_bp
        Base pairs.
    concentration_ng_ul
        Of that tube, or ``None`` when it is not measured yet.

    Raises
    ------
    ValueError
        If either number is not positive.
    """

    name: str
    length_bp: int
    _: KW_ONLY
    concentration_ng_ul: float | None = None

    def __post_init__(self) -> None:
        """Refuse a length or a concentration that is not positive."""
        if self.length_bp <= 0:
            raise ValueError(f"fragment {self.name!r}: length_bp must be positive")
        if self.concentration_ng_ul is not None and self.concentration_ng_ul <= 0:
            raise ValueError(f"fragment {self.name!r}: concentration_ng_ul must be positive")


@dataclass(frozen=True, slots=True)
class Amount:
    """How much of one fragment a reaction takes, in the three units a bench needs.

    Parameters
    ----------
    name, length_bp
        The fragment's.
    pmol
        What the protocol asks for.
    nanograms
        The same amount weighed.
    volume_ul
        What to pipette at the fragment's concentration, or `DNA_VOLUME_UL` without one.

    Raises
    ------
    ValueError
        If any of the three is not positive.
    """

    name: str
    length_bp: int
    _: KW_ONLY
    pmol: float
    nanograms: float
    volume_ul: float

    def __post_init__(self) -> None:
        """Refuse an amount no one can pipette."""
        for field, value in (
            ("pmol", self.pmol),
            ("nanograms", self.nanograms),
            ("volume_ul", self.volume_ul),
        ):
            if value <= 0:
                raise ValueError(f"amount {self.name!r}: {field} must be positive")


def assembly_amounts(
    vector: Fragment,
    inserts: tuple[Fragment, ...],
    *,
    vector_pmol: float = FRAGMENT_PMOL,
    insert_ratio: float = INSERT_RATIO,
) -> tuple[Amount, ...]:
    """Return what to put in the assembly, vector first.

    NEB asks for `FRAGMENT_PMOL` of the destination plasmid and the same of each precloned
    insert; `insert_ratio` is the insert to vector molar ratio for amplicon inserts.
    """
    return tuple(
        _amount(fragment, pmol)
        for fragment, pmol in (
            (vector, vector_pmol),
            *((insert, vector_pmol * insert_ratio) for insert in inserts),
        )
    )


def _amount(fragment: Fragment, pmol: float) -> Amount:
    nanograms = to_nanograms(pmol, fragment.length_bp)
    volume = (
        DNA_VOLUME_UL
        if fragment.concentration_ng_ul is None
        else nanograms / fragment.concentration_ng_ul
    )
    return Amount(
        fragment.name,
        fragment.length_bp,
        pmol=pmol,
        nanograms=round(nanograms, 2),
        volume_ul=round(volume, 2),
    )


# --------------------------------------------------------------------------------------
# The assembly reaction and its program
# --------------------------------------------------------------------------------------

#: Which of NEB's two Golden Gate systems a reaction follows.
type System = Literal["ligase master mix", "kit"]

#: NEBridge Ligase Master Mix (M1100), with a Type IIS enzyme of your own.
LIGASE_MASTER_MIX: System = "ligase master mix"

#: A NEBridge Golden Gate Assembly Kit, whose enzyme mix already holds T4 DNA Ligase.
KIT: System = "kit"


@dataclass(frozen=True, slots=True)
class Dose:
    """How much Type IIS enzyme one reaction takes, and the units that is."""

    volume_ul: float
    units: float


#: M1100: the Type IIS enzyme per reaction, for two fragments, three to six, and seven or more.
#: The stock a dose implies is its units over its volume, which is how BbsI-HF reaches 50 units
#: in one microlitre: that tier needs the R3539M concentrate.
_ENZYME_DOSE: Mapping[str, tuple[Dose, Dose, Dose]] = {
    "BbsI": (Dose(1.0, 20), Dose(1.0, 20), Dose(1.0, 50)),
    "BsaI": (Dose(1.0, 20), Dose(1.0, 20), Dose(1.0, 20)),
    "BsmBI": (Dose(3.0, 30), Dose(3.0, 30), Dose(6.0, 60)),
    "BspQI": (Dose(1.0, 10), Dose(1.0, 10), Dose(2.0, 20)),
    "Esp3I": (Dose(2.0, 20), Dose(3.0, 30), Dose(4.0, 40)),
    "PaqCI": (Dose(1.0, 10), Dose(1.0, 10), Dose(2.5, 25)),
    "SapI": (Dose(1.0, 10), Dose(1.0, 10), Dose(2.0, 20)),
}

#: The temperature NEB's cycling tables run each enzyme at, °C. It is not always the supplier's
#: digestion temperature, which is why `Enzyme.incubation_celsius` is not read here.
_GOLDEN_GATE_CELSIUS: Mapping[str, int] = {
    "BbsI": 37,
    "BsaI": 37,
    "BsmBI": 42,
    "BspQI": 42,
    "Esp3I": 37,
    "PaqCI": 37,
    "SapI": 37,
}

#: The kit built around each enzyme, by catalogue number.
_KIT_CATALOG: Mapping[str, str] = {"BsaI": "E1601", "BsmBI": "E1602"}

#: PaqCI is a tetramer that has to engage several copies of its site, so NEB supplies a short
#: activator duplex: µL per reaction for two fragments, three to six, and seven or more.
_ACTIVATOR_UL = (0.5, 0.5, 1.25)
_NEEDS_ACTIVATOR = frozenset({"PaqCI"})

#: PaqCI Activator as NEB supplies it, µM.
ACTIVATOR_UM = 20.0

#: M1100 totals, µL: the reaction and the 3X ligase master mix in it, up to six fragments and
#: from seven.
_MASTER_MIX_VOLUMES = ((15.0, 5.0), (30.0, 10.0))

#: The kit reaction, µL: the total and its 10X T4 DNA Ligase Buffer.
_KIT_VOLUME_UL = 20.0
_KIT_BUFFER_UL = 2.0

#: The background-reducing digest every Golden Gate program ends with. It is not heat
#: inactivation: it cuts destination plasmid that was never cut or has re-ligated.
END_SOAK_CELSIUS = 60.0
END_SOAK_SECONDS = 300

#: NEB asks for the fewest cycles that work when a Golden Gate insert is an amplicon.
GOLDEN_GATE_PCR_CYCLES = 20


def golden_gate_temperature(enzyme: Enzyme) -> float:
    """Return the temperature NEB's Golden Gate cycling runs this enzyme at, °C.

    Raises
    ------
    ValueError
        If NEB publishes no Golden Gate protocol for it.

    Examples
    --------
    >>> from liulab_mbio.enzymes import get_enzyme
    >>> golden_gate_temperature(get_enzyme("BbsI-HF"))
    37.0
    """
    if enzyme.name not in _GOLDEN_GATE_CELSIUS:
        raise ValueError(f"NEB publishes no Golden Gate protocol for {enzyme.name}")
    return float(_GOLDEN_GATE_CELSIUS[enzyme.name])


def assembly_reaction(
    enzyme: Enzyme,
    amounts: tuple[Amount, ...],
    *,
    system: System = LIGASE_MASTER_MIX,
    reactions: int = 1,
) -> ReactionTable:
    """Return the Golden Gate reaction for these fragments, the vector counted among them.

    DNA goes in each tube, everything else into the master mix, which is the order NEB asks for.
    The PaqCI activator line appears only for PaqCI.

    Raises
    ------
    ValueError
        If fewer than two fragments are given, if the enzyme has no NEB table for `system`, or
        if the DNA does not fit the reaction volume NEB's table sets.
    """
    if len(amounts) < 2:
        raise ValueError("a Golden Gate reaction joins at least two fragments")
    golden_gate_temperature(enzyme)
    total, components = (
        _kit_components(enzyme, amounts)
        if system == KIT
        else _master_mix_components(enzyme, amounts)
    )
    used = sum(component.volume_ul for component in components)
    if used >= total:
        raise ValueError(
            f"the DNA and enzyme take {used:g} µL of a {total:g} µL reaction; "
            "concentrate the fragments or scale the reaction up"
        )
    components.append(
        Component("Nuclease-free water", round(total - used, 2), final=f"to {total:g} µL")
    )
    return ReactionTable(
        tuple(components), title=_reaction_title(enzyme, system), reactions=reactions
    )


def _reaction_title(enzyme: Enzyme, system: System) -> str:
    if system == KIT:
        return f"Golden Gate assembly, NEBridge kit {_KIT_CATALOG[enzyme.name]}"
    return "Golden Gate assembly, NEBridge Ligase Master Mix (M1100)"


def _dna_components(amounts: tuple[Amount, ...]) -> list[Component]:
    return [
        Component(
            amount.name,
            amount.volume_ul,
            final=f"{amount.pmol:g} pmol ({amount.nanograms:g} ng)",
            master_mix=False,
        )
        for amount in amounts
    ]


def _master_mix_components(
    enzyme: Enzyme, amounts: tuple[Amount, ...]
) -> tuple[float, list[Component]]:
    fragments = len(amounts)
    total, _ = _master_mix_volumes(fragments)
    tier = _dose_tier(fragments)
    components = _dna_components(amounts)
    components.append(ligase_master_mix_component(fragments))
    components.append(enzyme_component(enzyme, fragments))
    if enzyme.name in _NEEDS_ACTIVATOR:
        components.append(
            Component(
                f"{enzyme.name} Activator",
                _ACTIVATOR_UL[tier],
                stock=f"{ACTIVATOR_UM:g} µM",
                final=f"{_ACTIVATOR_UL[tier] * ACTIVATOR_UM:g} pmol",
            )
        )
    return total, components


def ligase_master_mix_component(fragments: int) -> Component:
    """Return the NEBridge Ligase Master Mix component of a reaction joining `fragments`."""
    _, volume = _master_mix_volumes(fragments)
    return Component("NEBridge Ligase Master Mix", volume, stock="3X", final="1X")


def enzyme_component(enzyme: Enzyme, fragments: int) -> Component:
    """Return the Type IIS enzyme component of a Ligase Master Mix reaction joining `fragments`.

    Raises
    ------
    ValueError
        If NEB's Ligase Master Mix table has no row for the enzyme.
    """
    if enzyme.name not in _ENZYME_DOSE:
        raise ValueError(f"NEB's Ligase Master Mix table has no row for {enzyme.name}")
    dose = _ENZYME_DOSE[enzyme.name][_dose_tier(fragments)]
    return Component(
        enzyme.supplier_label,
        dose.volume_ul,
        stock=f"{dose.units / dose.volume_ul:g} U/µL",
        final=f"{dose.units:g} units",
    )


def _master_mix_volumes(fragments: int) -> tuple[float, float]:
    return _MASTER_MIX_VOLUMES[1 if fragments >= 7 else 0]


def _kit_components(enzyme: Enzyme, amounts: tuple[Amount, ...]) -> tuple[float, list[Component]]:
    if enzyme.name not in _KIT_CATALOG:
        raise ValueError(f"no NEBridge kit carries {enzyme.name}; use the Ligase Master Mix system")
    components = _dna_components(amounts)
    components.append(Component("T4 DNA Ligase Buffer", _KIT_BUFFER_UL, stock="10X", final="1X"))
    components.append(
        Component("NEBridge Golden Gate Enzyme Mix", 1.0 if len(amounts) - 1 <= 10 else 2.0)
    )
    return _KIT_VOLUME_UL, components


def _dose_tier(fragments: int) -> int:
    if fragments == 2:
        return 0
    return 1 if fragments <= 6 else 2


def assembly_program(
    enzyme: Enzyme,
    *,
    fragments: int,
    system: System = LIGASE_MASTER_MIX,
    library: bool = False,
) -> ThermocyclerProgram:
    """Return the cycling for this enzyme and fragment count, ending with the 60 °C end soak.

    `fragments` counts the vector. `library` picks NEB's longer single-insert incubation, which
    it recommends for library construction rather than for cloning one gene.

    Raises
    ------
    ValueError
        If fewer than two fragments are given, or the enzyme has no NEB protocol.
    """
    if fragments < 2:
        raise ValueError("a Golden Gate reaction joins at least two fragments")
    celsius = golden_gate_temperature(enzyme)
    stages = (
        _kit_stages(celsius, fragments, library=library)
        if system == KIT
        else _master_mix_stages(celsius, fragments, library=library)
    )
    end_soak = Stage((Incubation("End soak", END_SOAK_CELSIUS, END_SOAK_SECONDS),))
    return ThermocyclerProgram((*stages, end_soak), title="Golden Gate assembly")


def _cycle(celsius: float, seconds: int, cycles: int) -> Stage:
    return Stage(
        (Incubation("Digest", celsius, seconds), Incubation("Ligate", 16.0, seconds)),
        cycles=cycles,
    )


def _master_mix_stages(celsius: float, fragments: int, *, library: bool) -> tuple[Stage, ...]:
    if fragments == 2:
        if celsius == 37.0:
            return (Stage((Incubation("Assembly", celsius, 3600 if library else 900),)),)
        return (_cycle(celsius, 60, 30 if library else 15),)
    if fragments <= 6:
        return (_cycle(celsius, 60, 30),)
    return (_cycle(celsius, 300, 30 if fragments <= 13 else 60),)


def _kit_stages(celsius: float, fragments: int, *, library: bool) -> tuple[Stage, ...]:
    inserts = fragments - 1
    if inserts == 1:
        return (Stage((Incubation("Assembly", celsius, 3600 if library else 300),)),)
    return (_cycle(celsius, 60 if inserts <= 10 else 300, 30),)


def heat_inactivation(enzyme: Enzyme) -> ThermocyclerProgram | None:
    """Return the supplier's heat inactivation for this enzyme, or ``None`` where it gives none.

    This is not the 60 °C end soak `assembly_program` already ends with, which is a digest.
    """
    celsius, minutes = enzyme.heat_inactivation_celsius, enzyme.heat_inactivation_minutes
    if celsius is None or minutes is None:
        return None
    return ThermocyclerProgram(
        (Stage((Incubation("Heat inactivation", float(celsius), minutes * 60),)),),
        title=f"{enzyme.name} heat inactivation",
    )


# --------------------------------------------------------------------------------------
# PCR and colony PCR
# --------------------------------------------------------------------------------------

#: The tubes a PCR is pipetted from: each primer at 10 µM, and a dNTP mix at 10 mM of each base.
PRIMER_STOCK_UM = 10.0
DNTP_STOCK_MM = 10.0


#: NEB's colony PCR: a 2X master mix, a colony picked with a toothpick, and a lysis step long
#: enough to open the cells.
COLONY_PCR_MASTER_MIX = "OneTaq Quick-Load 2X Master Mix with Standard Buffer (M0486)"
COLONY_PCR_VOLUME_UL = 25.0
COLONY_LYSIS_SECONDS = 300
COLONY_HOLD_CELSIUS = 10.0


def pcr_reaction(
    polymerase: Polymerase = Q5,
    *,
    volume_ul: float = 50.0,
    reactions: int = 1,
    template_volume_ul: float = DNA_VOLUME_UL,
    title: str = "PCR",
) -> ReactionTable:
    """Return the PCR NEB's protocol sets up for this polymerase, in its own order.

    Volumes follow from the concentrations NEB's table gives: the buffer, the dNTP mix at
    `DNTP_STOCK_MM`, each primer at `PRIMER_STOCK_UM`, and the polymerase from its own tube.
    Template is added per tube.

    Raises
    ------
    ValueError
        If the components do not fit `volume_ul`.
    """
    profile = polymerase.pcr
    components = [
        Component(
            profile.buffer_name,
            round(volume_ul / profile.buffer_fold, 2),
            stock=f"{profile.buffer_fold:g}X",
            final="1X",
        ),
        Component(
            "dNTP mix",
            round(volume_ul * profile.dntp_um_each / (DNTP_STOCK_MM * 1000), 2),
            stock=f"{DNTP_STOCK_MM:g} mM each",
            final=f"{profile.dntp_um_each:g} µM each",
        ),
        *(
            Component(
                f"{end} primer",
                round(volume_ul * polymerase.primer_nm / 1000 / PRIMER_STOCK_UM, 2),
                stock=f"{PRIMER_STOCK_UM:g} µM",
                final=f"{polymerase.primer_nm:g} nM",
            )
            for end in ("Forward", "Reverse")
        ),
        Component("Template DNA", template_volume_ul, master_mix=False),
        Component(
            f"{polymerase.name} DNA Polymerase",
            round(volume_ul * profile.units_per_ul / profile.stock_units_ul, 2),
            stock=f"{profile.stock_units_ul:g} U/µL",
            final=f"{volume_ul * profile.units_per_ul:g} units",
        ),
    ]
    return _filled(components, volume_ul, title=title, reactions=reactions)


def colony_pcr_reaction(
    polymerase: Polymerase = ONETAQ,
    *,
    volume_ul: float = COLONY_PCR_VOLUME_UL,
    reactions: int = 1,
) -> ReactionTable:
    """Return NEB's colony PCR reaction, which is a 2X master mix and the two primers.

    The colony itself is no line of the table: it is picked with a toothpick and stirred into
    the tube until the solution clouds.

    Raises
    ------
    ValueError
        If the primers and master mix do not fit `volume_ul`.
    """
    components = [
        colony_pcr_master_mix_component(volume_ul),
        *(
            Component(
                f"{end} primer",
                round(volume_ul * polymerase.primer_nm / 1000 / PRIMER_STOCK_UM, 2),
                stock=f"{PRIMER_STOCK_UM:g} µM",
                final=f"{polymerase.primer_nm:g} nM",
            )
            for end in ("Forward", "Reverse")
        ),
    ]
    return _filled(components, volume_ul, title="Colony PCR", reactions=reactions)


def colony_pcr_master_mix_component(volume_ul: float = COLONY_PCR_VOLUME_UL) -> Component:
    """Return the 2X master mix component of a colony PCR of `volume_ul`."""
    return Component(COLONY_PCR_MASTER_MIX, round(volume_ul / 2, 2), stock="2X", final="1X")


def _filled(
    components: list[Component], volume_ul: float, *, title: str, reactions: int
) -> ReactionTable:
    used = sum(component.volume_ul for component in components)
    if used >= volume_ul:
        raise ValueError(f"the components take {used:g} µL of a {volume_ul:g} µL reaction")
    components.append(
        Component("Nuclease-free water", round(volume_ul - used, 2), final=f"to {volume_ul:g} µL")
    )
    return ReactionTable(tuple(components), title=title, reactions=reactions)


def pcr_program(
    polymerase: Polymerase = Q5,
    *,
    annealing_temperature: float,
    amplicon_length: int,
    cycles: int | None = None,
    title: str = "PCR",
) -> ThermocyclerProgram:
    """Return the program for this polymerase, annealing temperature and amplicon.

    Annealing and extension are combined into one step at the extension temperature once the
    annealing temperature reaches the polymerase's `PcrProfile.two_step_celsius`. An amplicon
    that is a Golden Gate insert wants `GOLDEN_GATE_PCR_CYCLES`, the fewest NEB finds enough.
    """
    profile = polymerase.pcr
    initial = Incubation(
        "Initial denaturation",
        profile.initial_denaturation_c,
        profile.initial_denaturation_seconds,
    )
    return _program(
        polymerase,
        initial,
        annealing_temperature=annealing_temperature,
        amplicon_length=amplicon_length,
        cycles=profile.cycles if cycles is None else cycles,
        hold_c=profile.hold_c,
        title=title,
    )


def colony_pcr_program(
    polymerase: Polymerase = ONETAQ,
    *,
    annealing_temperature: float,
    amplicon_length: int,
    cycles: int | None = None,
) -> ThermocyclerProgram:
    """Return the colony PCR program, which opens the cells before it denatures anything."""
    profile = polymerase.pcr
    lysis = Incubation("Lysis", profile.initial_denaturation_c, COLONY_LYSIS_SECONDS)
    return _program(
        polymerase,
        lysis,
        annealing_temperature=annealing_temperature,
        amplicon_length=amplicon_length,
        cycles=profile.cycles if cycles is None else cycles,
        hold_c=COLONY_HOLD_CELSIUS,
        title="Colony PCR",
    )


def _program(
    polymerase: Polymerase,
    first: Incubation,
    *,
    annealing_temperature: float,
    amplicon_length: int,
    cycles: int,
    hold_c: float,
    title: str,
) -> ThermocyclerProgram:
    profile = polymerase.pcr
    extension = polymerase.extension_seconds(amplicon_length)
    denature = Incubation("Denature", profile.denaturation_c, profile.denaturation_seconds)
    if annealing_temperature >= profile.two_step_celsius:
        inside = (
            denature,
            Incubation("Anneal and extend", polymerase.extension_temperature, extension),
        )
    else:
        inside = (
            denature,
            Incubation("Anneal", annealing_temperature, profile.annealing_seconds),
            Incubation("Extend", polymerase.extension_temperature, extension),
        )
    return ThermocyclerProgram(
        (
            Stage((first,)),
            Stage(inside, cycles=cycles),
            Stage(
                (
                    Incubation(
                        "Final extension",
                        polymerase.extension_temperature,
                        profile.final_extension_seconds,
                    ),
                )
            ),
            Stage((Incubation("Hold", hold_c, None),)),
        ),
        title=title,
    )


# --------------------------------------------------------------------------------------
# Colony PCR validation
# --------------------------------------------------------------------------------------

#: NEB's two ladders, each with 500 and 517 counted as the one band NEB counts them as.
LADDER_100_BP = Ladder(
    "NEB 100 bp DNA Ladder (N3231)",
    (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1517),
)
LADDER_1_KB_PLUS = Ladder(
    "NEB 1 kb Plus DNA Ladder (N3200)",
    (
        100,
        200,
        300,
        400,
        500,
        600,
        700,
        800,
        900,
        1000,
        1200,
        1517,
        2017,
        3001,
        4001,
        5001,
        6001,
        8001,
        10002,
    ),
)

#: Where the 100 bp ladder at 2% agarose gives way to the 1 kb Plus ladder at 1%, bp.
_LADDER_LIMIT = 1000

#: Vector kept either side of the junctions by a designed colony PCR pair, bases. Twice this is
#: the empty-vector band, and the note asks for every band to stay at 100 bp or more.
COLONY_FLANK = 60

#: How far into the insert a junction primer anneals, bases, for the same reason.
JUNCTION_OFFSET = 100

#: Genewiz asks for a sequencing primer 100 bases from what it reads, 50 to 60 at the closest.
SANGER_FLANK = 100

#: What the candidate plasmids are called, correct first. A reversed lane is numbered after
#: `REVERSED_CLONE` wherever an assembly holds more than one insert to turn round.
CORRECT_CLONE = "Correct clone"
EMPTY_CLONE = "Empty vector"
REVERSED_CLONE = "Reversed insert"

#: What a primer reading out of one insert is called, numbered the same way.
_JUNCTION_PRIMER = "Junction reverse"


def choose_ladder(bands_bp: tuple[int, ...]) -> Ladder:
    """Return the ladder covering these bands.

    Raises
    ------
    ValueError
        If no band is given.
    """
    return LADDER_100_BP if _largest(bands_bp) < _LADDER_LIMIT else LADDER_1_KB_PLUS


def agarose_percent(bands_bp: tuple[int, ...]) -> float:
    """Return the agarose percentage these bands resolve on.

    Raises
    ------
    ValueError
        If no band is given.
    """
    return 2.0 if _largest(bands_bp) < _LADDER_LIMIT else 1.0


def _largest(bands_bp: tuple[int, ...]) -> int:
    if not bands_bp:
        raise ValueError("a gel needs at least one expected band")
    return max(bands_bp)


@dataclass(frozen=True, slots=True)
class Clone:
    """One plasmid a colony may carry, and the bands a colony PCR gives from it."""

    name: str
    bands_bp: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class ColonyCheck:
    """A colony PCR reading across an assembly's junctions.

    Parameters
    ----------
    primers
        The flanking pair first, then one junction primer per insert where they were asked for.
    reports
        What each primer scored on the assembled plasmid.
    clones
        The candidates a colony can hold: the correct one, the empty vector, and one carrying
        each insert the other way round.
    annealing_temperature
        By the polymerase's rule over the two lowest Tms of the set, °C.
    extension_seconds
        For the largest band any candidate gives.
    ladder, agarose_percent
        Chosen for that band range.
    """

    primers: tuple[Primer, ...]
    reports: tuple[PrimerReport, ...]
    clones: tuple[Clone, ...]
    annealing_temperature: float
    extension_seconds: int
    ladder: Ladder
    agarose_percent: float

    @property
    def gel(self) -> Gel:
        """The gel these clones should give, one lane each."""
        return Gel(
            self.ladder,
            tuple(Lane(clone.name, clone.bands_bp) for clone in self.clones),
            title="Colony PCR",
        )

    @property
    def tells_orientation(self) -> bool:
        """Whether the gel separates every reversed insert from the correct clone.

        Two vector primers flanking the inserts never do: they amplify them whichever way round
        they sit. A junction primer does, unless the two vector primers happen to lie the same
        distance from their own junctions.
        """
        correct = next(
            (clone.bands_bp for clone in self.clones if clone.name == CORRECT_CLONE), None
        )
        turned = [clone.bands_bp for clone in self.clones if clone.name.startswith(REVERSED_CLONE)]
        return bool(turned) and all(lane != correct for lane in turned)


def colony_pcr_check(
    product: SequenceRecord,
    junctions: Sequence[int],
    *,
    vector: SequenceRecord,
    primers: tuple[Primer, ...] | None = None,
    insert_primer: bool = False,
    flank: int = COLONY_FLANK,
    junction_offset: int = JUNCTION_OFFSET,
    polymerase: Polymerase = ONETAQ,
    thresholds: Thresholds = THRESHOLDS,
) -> ColonyCheck:
    """Return what a colony PCR across these junctions should show.

    An assembly of n inserts has n + 1 junctions and the inserts are the spans between them, so
    a product whose inserts cross the origin is rotated first. Without `primers`, a pair is
    designed in the vector `flank` bases outside the first and the last junction;
    `insert_primer` adds one primer per insert, annealing `junction_offset` bases into it. Those
    are what tell a reversed insert apart and what put a band of their own on each junction.
    Each candidate plasmid is amplified on its own, so the bands are simulated rather than
    derived.

    Raises
    ------
    ValueError
        If the junctions are not two or more separate positions inside the product, if fewer
        than two primers are given, if an insert is too short for a junction primer, or if the
        primers amplify nothing at all.
    """
    places = _junction_span(junctions, len(product))
    start, end = places[0], places[-1]
    inserts = tuple(pairwise(places))
    chosen = (
        list(primers)
        if primers is not None
        else list(_flanking_pair(product, start, end, flank, polymerase, thresholds))
    )
    if insert_primer:
        chosen.extend(
            _junction_primer(product, first, last, junction_offset, polymerase, thresholds, name)
            for name, (first, last) in zip(
                _numbered(_JUNCTION_PRIMER, inserts), inserts, strict=True
            )
        )
    if len(chosen) < 2:
        raise ValueError("a colony PCR needs at least two primers")
    placed = tuple(chosen)
    candidates = [(CORRECT_CLONE, product), (EMPTY_CLONE, vector)]
    candidates += [
        (name, _reversed_insert(product, first, last))
        for name, (first, last) in zip(_numbered(REVERSED_CLONE, inserts), inserts, strict=True)
    ]
    clones = tuple(Clone(name, _bands(placed, record, thresholds)) for name, record in candidates)
    sizes = tuple(sorted({bp for clone in clones for bp in clone.bands_bp}))
    if not sizes:
        raise ValueError("these primers amplify nothing on any of the candidate plasmids")
    reports = tuple(
        evaluate_primer(primer, product, polymerase=polymerase, thresholds=thresholds)
        for primer in placed
    )
    tms = sorted(report["tm"].value for report in reports)
    return ColonyCheck(
        placed,
        reports,
        clones,
        polymerase.annealing_temperature(tms[0], tms[1]),
        polymerase.extension_seconds(max(sizes)),
        choose_ladder(sizes),
        agarose_percent(sizes),
    )


@dataclass(frozen=True, slots=True)
class SangerRead:
    """A sequencing primer and the read it has to give.

    Parameters
    ----------
    primer
        Reading towards the junction it sits outside.
    distance_bp
        From its 3' end to that junction.
    read_bp
        From its 3' end to the far junction, which is what the read must cover for the insert to
        be confirmed at both ends.
    """

    primer: Primer
    distance_bp: int
    read_bp: int


def sanger_primers(
    product: SequenceRecord,
    junctions: Sequence[int],
    *,
    flank: int = SANGER_FLANK,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
) -> tuple[SangerRead, SangerRead]:
    """Return a sequencing primer reading into the inserts from outside the first and last junction.

    Every 3' end lands at least `flank` bases from its own junction, near enough for the read to
    be clean there, and `SangerRead.read_bp` is how far it must carry to reach the far junction.
    A provider whose reads are shorter needs a primer inside the inserts as well.

    Raises
    ------
    ValueError
        If the junctions are not two or more separate positions inside the product, or no primer
        fits outside the first or the last.
    """
    places = _junction_span(junctions, len(product))
    start, end = places[0], places[-1]
    longest = _longest_annealing(thresholds)
    forward = design_primer(
        product,
        start - flank - longest,
        Strand.FORWARD,
        name="Sequencing forward",
        polymerase=polymerase,
        thresholds=thresholds,
    )
    reverse = design_primer(
        product,
        end + flank + longest,
        Strand.REVERSE,
        name="Sequencing reverse",
        polymerase=polymerase,
        thresholds=thresholds,
    )
    length = len(product)
    ahead = forward.binding_sites[0].end
    behind = reverse.binding_sites[0].start
    return (
        SangerRead(forward, (start - ahead) % length, (end - ahead) % length),
        SangerRead(reverse, (behind - end) % length, (behind - start) % length),
    )


def _junction_span(junctions: Sequence[int], length: int) -> tuple[int, ...]:
    """Return the junctions in rising order, refusing a set no assembly could leave."""
    if len(junctions) < 2:
        raise ValueError(
            f"an assembly has a junction at each end of every insert, got {len(junctions)}"
        )
    places = tuple(sorted(junctions))
    if len(set(places)) != len(places):
        raise ValueError(f"two junctions share a position: {places}")
    if places[0] < 0 or places[-1] > length:
        raise ValueError(f"junctions {places} do not lie inside {length} bases")
    return places


def _numbered(label: str, inserts: Sequence[tuple[int, int]]) -> tuple[str, ...]:
    """Label one thing per insert, numbered only where there is more than one to tell apart."""
    if len(inserts) == 1:
        return (label,)
    return tuple(f"{label} {number}" for number in range(1, len(inserts) + 1))


def _flanking_pair(
    product: SequenceRecord,
    start: int,
    end: int,
    flank: int,
    polymerase: Polymerase,
    thresholds: Thresholds,
) -> tuple[Primer, Primer]:
    return design_pair(
        product,
        start - flank,
        end + flank,
        forward_name="Colony PCR forward",
        reverse_name="Colony PCR reverse",
        polymerase=polymerase,
        thresholds=thresholds,
    )


def _junction_primer(
    product: SequenceRecord,
    start: int,
    end: int,
    offset: int,
    polymerase: Polymerase,
    thresholds: Thresholds,
    name: str = _JUNCTION_PRIMER,
) -> Primer:
    if end - start <= offset:
        raise ValueError(f"an insert is shorter than the {offset} bases a junction primer needs")
    return design_primer(
        product,
        start + offset,
        Strand.REVERSE,
        name=name,
        polymerase=polymerase,
        thresholds=thresholds,
    )


def _reversed_insert(product: SequenceRecord, start: int, end: int) -> SequenceRecord:
    flipped, _ = edits.replace(product, start, end, reverse_complement(product.sequence[start:end]))
    return flipped


def _bands(
    primers: tuple[Primer, ...], record: SequenceRecord, thresholds: Thresholds
) -> tuple[int, ...]:
    """Return every band this primer set gives on one candidate plasmid.

    Binding sites are cleared first: they were found on the product, and each candidate has to
    be searched on its own.
    """
    free = [dataclasses.replace(primer, binding_sites=()) for primer in primers]
    sizes: set[int] = set()
    for index, one in enumerate(free):
        for other in free[index + 1 :]:
            sizes.update(amplicon_sizes(one, other, record, thresholds=thresholds))
    return tuple(sorted(sizes))


def _longest_annealing(thresholds: Thresholds) -> int:
    """Return the longest annealing region `design_primer` will choose."""
    band = thresholds.length
    return int(band.high if band.warn_high == float("inf") else band.warn_high)


#: Where the numbers above come from, ready for a protocol's reference list.
REFERENCES: tuple[Reference, ...] = (
    Reference(
        "NEB, NEBridge Golden Gate Assembly Kit (BsaI-HFv2) instruction manual, NEB #E1601S/L, "
        "version 5.0_6/26",
        url="https://www.neb.com/-/media/nebus/files/manuals/manuale1601.pdf",
    ),
    Reference(
        "NEB, Protocol for NEBridge Ligase Master Mix (NEB #M1100)",
        url="https://web.archive.org/web/20230331002719id_/https://www.neb.com/protocols/2021/09/14/protocol-for-nebridge-ligase-master-mix-neb-m1100",
    ),
    Reference(
        "NEB, NEBridge Ligase Master Mix Protocol Guidelines",
        url="https://web.archive.org/web/20250713174145id_/https://www.neb.com/en-us/tools-and-resources/usage-guidelines/nebridge-ligase-master-mix-protocol-guidelines",
    ),
    Reference(
        "NEB, Usage Guidelines for Golden Gate Assembly with PaqCI",
        url="https://web.archive.org/web/20210615031818id_/https://www.neb.com/tools-and-resources/usage-guidelines/usage-guidelines-for-golden-gate-assembly-with-paqci",
    ),
    Reference(
        "NEB, Robust Colony PCR from Multiple E. coli Strains using OneTaq Quick-Load Master "
        "Mixes (Y. Xu, 11/13)"
    ),
    Reference("NEB, Nucleic Acid Data, and NEBioCalculator for the ng to pmol conversion"),
    Reference("NEB product pages: 1 kb Plus DNA Ladder (N3200), 100 bp DNA Ladder (N3231)"),
    Reference("NEB, Agarose Gel Resolution, for the percentage a band range resolves on"),
)

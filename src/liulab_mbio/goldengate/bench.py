"""What a Golden Gate experiment takes at the bench: quantities, reactions and programs.

Every number is NEB's, through ``docs/research/golden-gate-assembly.md`` for the assembly and
``docs/research/primer-design-and-pcr.md`` for the conversions, the PCR programs and the
ladders. Functions return `liulab_mbio.protocol` values, so a protocol prints them unchanged.

NEB ships two Golden Gate systems and their tables do not mix. `LIGASE_MASTER_MIX` is NEBridge
Ligase Master Mix (M1100), which takes any NEB Type IIS enzyme; `KIT` is one of the kits, which
carry their own enzyme mix and so exist only for BsaI-HFv2 and BsmBI-v2.
"""

from collections.abc import Mapping
from dataclasses import KW_ONLY, dataclass
from typing import Literal

from liulab_mbio.enzymes import Enzyme
from liulab_mbio.primers import ONETAQ, Q5, Polymerase
from liulab_mbio.protocol import (
    Component,
    Incubation,
    ReactionTable,
    Stage,
    ThermocyclerProgram,
)

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
    if enzyme.name not in _ENZYME_DOSE:
        raise ValueError(f"NEB's Ligase Master Mix table has no row for {enzyme.name}")
    total, master_mix = _MASTER_MIX_VOLUMES[1 if fragments >= 7 else 0]
    tier = _dose_tier(fragments)
    dose = _ENZYME_DOSE[enzyme.name][tier]
    components = _dna_components(amounts)
    components.append(Component("NEBridge Ligase Master Mix", master_mix, stock="3X", final="1X"))
    components.append(
        Component(
            _enzyme_label(enzyme),
            dose.volume_ul,
            stock=f"{dose.units / dose.volume_ul:g} U/µL",
            final=f"{dose.units:g} units",
        )
    )
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


def _kit_components(enzyme: Enzyme, amounts: tuple[Amount, ...]) -> tuple[float, list[Component]]:
    if enzyme.name not in _KIT_CATALOG:
        raise ValueError(f"no NEBridge kit carries {enzyme.name}; use the Ligase Master Mix system")
    components = _dna_components(amounts)
    components.append(Component("T4 DNA Ligase Buffer", _KIT_BUFFER_UL, stock="10X", final="1X"))
    components.append(
        Component("NEBridge Golden Gate Enzyme Mix", 1.0 if len(amounts) - 1 <= 10 else 2.0)
    )
    return _KIT_VOLUME_UL, components


def _enzyme_label(enzyme: Enzyme) -> str:
    name = enzyme.commercial_name or enzyme.name
    return f"{name} ({enzyme.catalog_number})" if enzyme.catalog_number else name


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


@dataclass(frozen=True, slots=True)
class PcrProfile:
    """What NEB's protocol puts in one polymerase's PCR, and the program around it.

    Parameters
    ----------
    buffer_name, buffer_fold
        The reaction buffer as supplied, so its volume is the reaction over `buffer_fold`.
    units_per_ul
        Polymerase in the reaction.
    stock_units_ul
        Polymerase in the tube it is pipetted from.
    initial_denaturation_c, initial_denaturation_seconds
        The step before the cycles.
    denaturation_c, denaturation_seconds, annealing_seconds
        Inside each cycle; the annealing temperature is the primer pair's.
    final_extension_seconds
        At the polymerase's own extension temperature.
    two_step_celsius
        The lowest annealing temperature that gets a two-step program, annealing and extension
        combined. Taq's rule is "above 65 °C", which at a tenth of a degree is 65.1.
    cycles, dntp_um_each, hold_c
        The rest of NEB's table.
    """

    buffer_name: str
    _: KW_ONLY
    buffer_fold: float
    units_per_ul: float
    stock_units_ul: float
    initial_denaturation_c: float
    denaturation_c: float
    denaturation_seconds: int
    annealing_seconds: int
    final_extension_seconds: int
    two_step_celsius: float
    cycles: int = 30
    initial_denaturation_seconds: int = 30
    dntp_um_each: float = 200.0
    hold_c: float = 4.0


#: One profile per polymerase `liulab_mbio.primers` ships, keyed by its name.
PCR_PROFILES: Mapping[str, PcrProfile] = {
    "Q5": PcrProfile(
        "Q5 Reaction Buffer",
        buffer_fold=5.0,
        units_per_ul=0.02,
        stock_units_ul=2.0,
        initial_denaturation_c=98.0,
        denaturation_c=98.0,
        denaturation_seconds=10,
        annealing_seconds=20,
        final_extension_seconds=120,
        two_step_celsius=72.0,
    ),
    "Phusion": PcrProfile(
        "Phusion HF Buffer",
        buffer_fold=5.0,
        units_per_ul=0.02,
        stock_units_ul=2.0,
        initial_denaturation_c=98.0,
        denaturation_c=98.0,
        denaturation_seconds=10,
        annealing_seconds=20,
        final_extension_seconds=300,
        two_step_celsius=72.0,
    ),
    "OneTaq": PcrProfile(
        "OneTaq Standard Reaction Buffer",
        buffer_fold=5.0,
        units_per_ul=0.025,
        stock_units_ul=5.0,
        initial_denaturation_c=94.0,
        denaturation_c=94.0,
        denaturation_seconds=30,
        annealing_seconds=30,
        final_extension_seconds=300,
        two_step_celsius=68.0,
    ),
    "Taq": PcrProfile(
        "Standard Taq Reaction Buffer",
        buffer_fold=10.0,
        units_per_ul=0.025,
        stock_units_ul=5.0,
        initial_denaturation_c=95.0,
        denaturation_c=95.0,
        denaturation_seconds=30,
        annealing_seconds=30,
        final_extension_seconds=300,
        two_step_celsius=65.1,
    ),
}

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
        If the polymerase has no profile, or the components do not fit `volume_ul`.
    """
    profile = _profile(polymerase)
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
        Component(COLONY_PCR_MASTER_MIX, round(volume_ul / 2, 2), stock="2X", final="1X"),
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

    Raises
    ------
    ValueError
        If the polymerase has no profile.
    """
    profile = _profile(polymerase)
    initial = Incubation(
        "Initial denaturation",
        profile.initial_denaturation_c,
        profile.initial_denaturation_seconds,
    )
    return _program(
        profile,
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
    """Return the colony PCR program, which opens the cells before it denatures anything.

    Raises
    ------
    ValueError
        If the polymerase has no profile.
    """
    profile = _profile(polymerase)
    lysis = Incubation("Lysis", profile.initial_denaturation_c, COLONY_LYSIS_SECONDS)
    return _program(
        profile,
        polymerase,
        lysis,
        annealing_temperature=annealing_temperature,
        amplicon_length=amplicon_length,
        cycles=profile.cycles if cycles is None else cycles,
        hold_c=COLONY_HOLD_CELSIUS,
        title="Colony PCR",
    )


def _program(
    profile: PcrProfile,
    polymerase: Polymerase,
    first: Incubation,
    *,
    annealing_temperature: float,
    amplicon_length: int,
    cycles: int,
    hold_c: float,
    title: str,
) -> ThermocyclerProgram:
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


def _profile(polymerase: Polymerase) -> PcrProfile:
    if polymerase.name not in PCR_PROFILES:
        raise ValueError(f"no NEB PCR protocol is recorded for {polymerase.name}")
    return PCR_PROFILES[polymerase.name]

"""What a Golden Gate assembly takes at the bench: its reaction and its cycling.

Every number is NEB's, through ``docs/research/golden-gate-assembly.md``. Functions return
`liulab_mbio.protocol` values, so a protocol prints them unchanged. What any cloning pipeline
shares -- DNA amounts, PCR, gels and validation -- is `liulab_mbio.bench`.

NEB ships two Golden Gate systems and their tables do not mix. `LIGASE_MASTER_MIX` is NEBridge
Ligase Master Mix (M1100), which takes any NEB Type IIS enzyme; `KIT` is one of the kits, which
carry their own enzyme mix and so exist only for BsaI-HFv2 and BsmBI-v2.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal

from liulab_mbio.bench.amounts import Amount, dna_amount
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.protocol import (
    Component,
    Incubation,
    ReactionTable,
    Reference,
    Stage,
    ThermocyclerProgram,
)

#: Picomoles of each fragment NEB's Golden Gate reactions take, whichever system is used.
FRAGMENT_PMOL = 0.05

#: Insert to vector molar ratio, from the E1601 manual revision 5.0_6/26 for amplicon inserts.
#: An earlier revision of the same page asked for 2:1.
INSERT_RATIO = 1.0


def assembly_amounts(
    vector: tuple[str, int],
    inserts: tuple[tuple[str, int], ...],
    *,
    vector_pmol: float = FRAGMENT_PMOL,
    insert_ratio: float = INSERT_RATIO,
) -> tuple[Amount, ...]:
    """Return what to put in the assembly, vector first.

    The vector and each insert are a name and a length in base pairs. NEB asks for
    `FRAGMENT_PMOL` of the destination plasmid and the same of each precloned insert;
    `insert_ratio` is the insert to vector molar ratio for amplicon inserts.

    Raises
    ------
    ValueError
        If a length is not positive.
    """
    return tuple(
        dna_amount(name, length_bp, pmol=pmol)
        for (name, length_bp), pmol in (
            (vector, vector_pmol),
            *((insert, vector_pmol * insert_ratio) for insert in inserts),
        )
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
)

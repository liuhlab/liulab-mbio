"""What a Golden Gate assembly takes at the bench: its enzymes, its reaction and its cycling.

Every number is NEB's, through ``docs/research/golden-gate-assembly.md``, and none of it is one
method's choice: `mbio.cloning.goldengate` and `synbio.igga` run the same tables.
Functions return `mbio.protocol` values, so a protocol prints them unchanged.

NEB ships two Golden Gate systems and their tables do not mix. `LIGASE_MASTER_MIX` is NEBridge
Ligase Master Mix (M1100), which takes any NEB Type IIS enzyme; `KIT` is one of the kits, which
carry their own enzyme mix and so exist only for BsaI-HFv2 and BsmBI-v2.
"""

from collections.abc import Mapping
from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import Literal

from mbio.bench.amounts import Amount, dna_amount
from mbio.bench.reactions import reaction_table
from mbio.enzymes import Enzyme
from mbio.overhangs import FidelityReport, ligation_source
from mbio.protocol.model import (
    Citation,
    Component,
    Incubation,
    ReactionTable,
    Source,
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

#: The enzymes a Golden Gate design may still reach for although nothing above covers them:
#: BtgZI cuts outside its site but NEB publishes no protocol for it, and its cut ends re-ligate
#: poorly.
LAST_RESORT = frozenset({"BtgZI"})

#: The Type IIS enzymes a Golden Gate design ranks by default, the last resorts last. Derived
#: from the cycling table, so an enzyme this module stops covering leaves the list with it.
GOLDEN_GATE_ENZYMES: tuple[str, ...] = (*_GOLDEN_GATE_CELSIUS, *sorted(LAST_RESORT))

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

#: NEB asks for the fewest cycles that work when a Golden Gate insert is an amplicon, and where
#: it asks. The count is the kit manual's, not the polymerase's, so it travels with its own
#: citation rather than through `PcrProfile.cycles`.
GOLDEN_GATE_PCR_CYCLES = 20
GOLDEN_GATE_PCR_CYCLES_CITATION = Citation("E1601", "FAQ 11")

#: The documents a citation in a Golden Gate protocol resolves against, each cited by the row
#: that rests on it, so a protocol borrowing one part of this module lists that part's own.
SOURCES: Mapping[str, Source] = MappingProxyType(
    {
        "E1601": Source(
            "New England Biolabs #E1601S/L NEBridge Golden Gate Assembly Kit (BsaI-HFv2) "
            "instruction manual",
            edition="version 5.0_6/26",
            url="https://www.neb.com/-/media/nebus/files/manuals/manuale1601.pdf",
            date="2026-09-12",
        ),
        "E1602": Source(
            "New England Biolabs #E1602S/L NEBridge Golden Gate Assembly Kit (BsmBI-v2) "
            "instruction manual",
            edition="version 3.0_6/26",
            url="https://www.neb.com/-/media/nebus/files/manuals/manuale1602.pdf",
            date="2026-09-12",
            note="docs/research/golden-gate-assembly.md",
        ),
        "M1100": Source(
            "New England Biolabs, Protocol for NEBridge Ligase Master Mix (NEB #M1100)",
            edition="capture 2023-03-31",
            url="https://web.archive.org/web/20230331002719id_/https://www.neb.com/protocols/2021/09/14/protocol-for-nebridge-ligase-master-mix-neb-m1100",
            read_as="Wayback Machine",
            note="docs/research/golden-gate-assembly.md",
        ),
        "M1100-guidelines": Source(
            "New England Biolabs, NEBridge Ligase Master Mix Protocol Guidelines",
            edition="capture 2025-07-13",
            url="https://web.archive.org/web/20250713174145id_/https://www.neb.com/en-us/tools-and-resources/usage-guidelines/nebridge-ligase-master-mix-protocol-guidelines",
            read_as="Wayback Machine",
            note="docs/research/golden-gate-assembly.md",
        ),
        "pryor-2020": Source(
            ligation_source().citation,
            url=ligation_source().doi_url,
            note="docs/research/ligation-fidelity.md",
        ),
    }
)

#: Where each system's reaction and cycling are read. A kit's numbers come from that kit's own
#: manual, so the citation is keyed by the kit the enzyme is sold in.
_MIX_REACTION = Citation("M1100", "reaction")
_MIX_DOSE = Citation("M1100", "enzyme amounts")
_MIX_CYCLING = Citation("M1100", "cycling")
_ACTIVATOR_CITATION = Citation("M1100-guidelines", "PaqCI Activator")


#: The key a ligase profile the user holds is cited under, where one scored a set.
PROFILE_KEY = "ligase-profile"


def fidelity_sources(report: FidelityReport) -> dict[str, Source]:
    """Return the documents a fidelity score may be read from: the paper, and the user's file.

    A protocol merges them into its `sources`; `citing` keeps the one `fidelity_citation` names.
    """
    return {
        "pryor-2020": SOURCES["pryor-2020"],
        PROFILE_KEY: Source(f"Ligase profile, {report.source}"),
    }


def fidelity_citation(report: FidelityReport) -> Citation | None:
    """Return where a fidelity score was read, or ``None`` where the rules scored the set.

    Examples
    --------
    >>> from mbio.overhangs import fidelity
    >>> fidelity_citation(fidelity(["AATG", "GCTT"], "BsaI-HFv2")).source
    'pryor-2020'
    """
    if not report.measured:
        return None
    if report.enzyme_specific or report.stand_in:
        return Citation("pryor-2020", report.source)
    return Citation(PROFILE_KEY)


def kit_citation(enzyme: Enzyme) -> Citation:
    """Return where the kit reaction for this enzyme is read: that kit's own manual.

    Examples
    --------
    >>> from mbio.enzymes import get_enzyme
    >>> kit_citation(get_enzyme("BsmBI-v2")).source
    'E1602'
    """
    return Citation(_KIT_CATALOG[enzyme.name], "assembly protocol")


def golden_gate_temperature(enzyme: Enzyme) -> float:
    """Return the temperature NEB's Golden Gate cycling runs this enzyme at, °C.

    Raises
    ------
    ValueError
        If NEB publishes no Golden Gate protocol for it.

    Examples
    --------
    >>> from mbio.enzymes import get_enzyme
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
    measured: bool = False,
) -> ReactionTable:
    """Return the Golden Gate reaction for these fragments, the vector counted among them.

    DNA goes in each tube, everything else into the master mix, which is the order NEB asks for.
    The PaqCI activator line appears only for PaqCI. `measured` is `reaction_table`'s.

    Raises
    ------
    ValueError
        If fewer than two fragments are given, if the enzyme has no NEB table for `system`, or
        if the DNA does not fit the reaction volume NEB's table sets.
    """
    if len(amounts) < 2:
        raise ValueError("a Golden Gate reaction joins at least two fragments")
    golden_gate_temperature(enzyme)
    fragments = len(amounts)
    total, components = (
        _kit_components(enzyme, fragments)
        if system == KIT
        else _master_mix_components(enzyme, fragments)
    )
    return reaction_table(
        amounts,
        components,
        volume_ul=total,
        title=_reaction_title(enzyme, system),
        reactions=reactions,
        measured=measured,
    )


def _reaction_title(enzyme: Enzyme, system: System) -> str:
    if system == KIT:
        return f"Golden Gate assembly, NEBridge kit {_KIT_CATALOG[enzyme.name]}"
    return "Golden Gate assembly, NEBridge Ligase Master Mix (M1100)"


def _master_mix_components(enzyme: Enzyme, fragments: int) -> tuple[float, list[Component]]:
    total, _ = _master_mix_volumes(fragments)
    tier = _dose_tier(fragments)
    components = [ligase_master_mix_component(fragments), enzyme_component(enzyme, fragments)]
    if enzyme.name in _NEEDS_ACTIVATOR:
        components.append(
            Component(
                f"{enzyme.name} Activator",
                _ACTIVATOR_UL[tier],
                stock=f"{ACTIVATOR_UM:g} µM",
                final=f"{_ACTIVATOR_UL[tier] * ACTIVATOR_UM:g} pmol",
                citation=_ACTIVATOR_CITATION,
            )
        )
    return total, components


def ligase_master_mix_component(fragments: int) -> Component:
    """Return the NEBridge Ligase Master Mix component of a reaction joining `fragments`."""
    _, volume = _master_mix_volumes(fragments)
    return Component(
        "NEBridge Ligase Master Mix", volume, stock="3X", final="1X", citation=_MIX_REACTION
    )


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
        citation=_MIX_DOSE,
    )


def _master_mix_volumes(fragments: int) -> tuple[float, float]:
    return _MASTER_MIX_VOLUMES[1 if fragments >= 7 else 0]


def _kit_components(enzyme: Enzyme, fragments: int) -> tuple[float, list[Component]]:
    if enzyme.name not in _KIT_CATALOG:
        raise ValueError(f"no NEBridge kit carries {enzyme.name}; use the Ligase Master Mix system")
    cited = kit_citation(enzyme)
    components = [
        Component("T4 DNA Ligase Buffer", _KIT_BUFFER_UL, stock="10X", final="1X", citation=cited),
        Component(
            "NEBridge Golden Gate Enzyme Mix",
            1.0 if fragments - 1 <= 10 else 2.0,
            citation=cited,
        ),
    ]
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
    cited = kit_citation(enzyme) if system == KIT else _MIX_CYCLING
    return ThermocyclerProgram(
        tuple(replace(stage, citation=cited) for stage in (*stages, end_soak)),
        title="Golden Gate assembly",
    )


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

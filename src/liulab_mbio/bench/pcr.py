"""PCR and colony PCR at the bench: the reaction and the thermocycler program for each.

Every number is NEB's, through ``docs/research/primer-design-and-pcr.md``; each polymerase
carries its own profile. Functions return `liulab_mbio.protocol` values, so a protocol prints
them unchanged.
"""

from collections.abc import Mapping
from types import MappingProxyType

from liulab_mbio.bench.amounts import DNA_VOLUME_UL
from liulab_mbio.primers.polymerase import ONETAQ, Q5, Polymerase
from liulab_mbio.protocol.model import (
    Citation,
    Component,
    Incubation,
    ReactionTable,
    Reference,
    Source,
    Stage,
    ThermocyclerProgram,
)

#: The tubes a PCR is pipetted from: each primer at 10 µM, and a dNTP mix at 10 mM of each base.
PRIMER_STOCK_UM = 10.0
DNTP_STOCK_MM = 10.0

#: NEB's protocol for each polymerase, keyed as its `PcrProfile.cycles_source` names it. Read
#: through ``docs/research/primer-design-and-pcr.md`` §3.1.
SOURCES: Mapping[str, Source] = MappingProxyType(
    {
        "M0491": Source(
            "New England Biolabs #M0491 Q5 High-Fidelity DNA Polymerase protocol page",
            edition="capture 2024-08-04",
            read_as="Wayback Machine",
            date="2026-09-12",
        ),
        "E0553": Source(
            "New England Biolabs #E0553 Phusion High-Fidelity PCR Kit manual",
            date="2026-09-12",
        ),
        "M0273": Source(
            "New England Biolabs #M0273 Taq DNA Polymerase protocol",
            edition="protocols.io version 1, 2015-01-29",
            url="https://dx.doi.org/10.17504/protocols.io.ch7t9m",
            date="2026-09-12",
        ),
        "M0480": Source(
            "New England Biolabs #M0480 OneTaq DNA Polymerase protocol",
            edition="protocols.io version 2, 2022-02-21",
            url="https://dx.doi.org/10.17504/protocols.io.bd24i8gw",
            date="2026-09-12",
        ),
    }
)

#: Where in each of those a cycle count stands.
CYCLES_LOCATOR = "thermocycling conditions"


def cycle_citation(polymerase: Polymerase = Q5) -> Citation:
    """Return where this polymerase's own `PcrProfile.cycles` was read.

    A caller taking the profile's count takes this with it, so the page says whose protocol the
    count is.

    Examples
    --------
    >>> cycle_citation(Q5).source
    'M0491'
    """
    return Citation(polymerase.pcr.cycles_source, CYCLES_LOCATOR)


def polymerase_name(polymerase: Polymerase = Q5) -> str:
    """Return what the tube is called: the row `pcr_reaction` pipettes, and its material's name.

    One name, so a step pipetting it names the material it comes from.

    Examples
    --------
    >>> polymerase_name(Q5)
    'Q5 DNA Polymerase'
    """
    return f"{polymerase.name} DNA Polymerase"


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
            polymerase_name(polymerase),
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
    cycles: int | None,
    cycles_citation: Citation | None,
    title: str = "PCR",
) -> ThermocyclerProgram:
    """Return the program for this polymerase, annealing temperature and amplicon.

    Annealing and extension are combined into one step at the extension temperature once the
    annealing temperature reaches the polymerase's `PcrProfile.two_step_celsius`.

    Parameters
    ----------
    cycles
        How many cycles to run. `PcrProfile.cycles` is NEB's own count for this polymerase, and
        ``None`` leaves the count blank where nothing sources it. It has no default: a count
        nobody chose prints on the page as if someone had.
    cycles_citation
        Where `cycles` was read, `cycle_citation` for the polymerase's own count. It has no
        default either, for the same reason `cycles` has none: a default lets a caller leave
        provenance off without saying so, where ``None`` says the count has no source to give.
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
        cycles=cycles,
        cycles_citation=cycles_citation,
        hold_c=profile.hold_c,
        title=title,
    )


def colony_pcr_program(
    polymerase: Polymerase = ONETAQ,
    *,
    annealing_temperature: float,
    amplicon_length: int,
    cycles: int | None,
    cycles_citation: Citation | None,
) -> ThermocyclerProgram:
    """Return the colony PCR program, which opens the cells before it denatures anything.

    `cycles` and `cycles_citation` read as `pcr_program`'s do, neither with a default.
    """
    profile = polymerase.pcr
    lysis = Incubation("Lysis", profile.initial_denaturation_c, COLONY_LYSIS_SECONDS)
    return _program(
        polymerase,
        lysis,
        annealing_temperature=annealing_temperature,
        amplicon_length=amplicon_length,
        cycles=cycles,
        cycles_citation=cycles_citation,
        hold_c=COLONY_HOLD_CELSIUS,
        title="Colony PCR",
    )


def _program(
    polymerase: Polymerase,
    first: Incubation,
    *,
    annealing_temperature: float,
    amplicon_length: int,
    cycles: int | None,
    cycles_citation: Citation | None,
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
            Stage(inside, cycles=cycles, citation=cycles_citation),
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


#: Where the numbers above come from, ready for a protocol's reference list.
REFERENCES: tuple[Reference, ...] = (
    Reference(
        "NEB, Robust Colony PCR from Multiple E. coli Strains using OneTaq Quick-Load Master "
        "Mixes (Y. Xu, 11/13)"
    ),
)

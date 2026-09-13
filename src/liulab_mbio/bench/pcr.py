"""PCR and colony PCR at the bench: the reaction and the thermocycler program for each.

Every number is NEB's, through ``docs/research/primer-design-and-pcr.md``; each polymerase
carries its own profile. Functions return `liulab_mbio.protocol` values, so a protocol prints
them unchanged.
"""

from liulab_mbio.bench.amounts import DNA_VOLUME_UL
from liulab_mbio.primers import ONETAQ, Q5, Polymerase
from liulab_mbio.protocol import Component, Incubation, ReactionTable, Stage, ThermocyclerProgram

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

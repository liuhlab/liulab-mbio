"""PCR and colony PCR pinned against NEB's own tables.

The source is `docs/research/primer-design-and-pcr.md` for the PCR programs.
"""

import pytest

from liulab_mbio.bench.pcr import (
    COLONY_PCR_MASTER_MIX,
    colony_pcr_master_mix_component,
    colony_pcr_program,
    colony_pcr_reaction,
    pcr_program,
    pcr_reaction,
)
from liulab_mbio.primers import POLYMERASES, Q5, Polymerase

from ..reactions import total, volumes


def test_the_colony_pcr_master_mix_component_is_the_one_its_reaction_prints() -> None:
    assert colony_pcr_master_mix_component() in colony_pcr_reaction().components


def test_a_high_annealing_temperature_combines_annealing_and_extension() -> None:
    program = pcr_program(Q5, annealing_temperature=72.0, amplicon_length=800)
    cycled = program.stages[1]
    assert [i.label for i in cycled.incubations] == ["Denature", "Anneal and extend"]
    assert cycled.incubations[1].temperature_c == 72.0


#: NEB's 50 µL PCR for each shipped polymerase, line by line as (name, µL, stock, final); its
#: program for 1.5 kb annealing at 55 °C as (label, °C, seconds); and the lowest two-step Ta.
SHIPPED_PCR = {
    "Q5": (
        [
            ("Q5 Reaction Buffer", 10.0, "5X", "1X"),
            ("dNTP mix", 1.0, "10 mM each", "200 µM each"),
            ("Forward primer", 2.5, "10 µM", "500 nM"),
            ("Reverse primer", 2.5, "10 µM", "500 nM"),
            ("Template DNA", 1.0, "", ""),
            ("Q5 DNA Polymerase", 0.5, "2 U/µL", "1 units"),
            ("Nuclease-free water", 32.5, "", "to 50 µL"),
        ],
        [
            ("Initial denaturation", 98.0, 30),
            ("Denature", 98.0, 10),
            ("Anneal", 55.0, 20),
            ("Extend", 72.0, 40),
            ("Final extension", 72.0, 120),
            ("Hold", 4.0, None),
        ],
        72.0,
    ),
    "Phusion": (
        [
            ("Phusion HF Buffer", 10.0, "5X", "1X"),
            ("dNTP mix", 1.0, "10 mM each", "200 µM each"),
            ("Forward primer", 2.5, "10 µM", "500 nM"),
            ("Reverse primer", 2.5, "10 µM", "500 nM"),
            ("Template DNA", 1.0, "", ""),
            ("Phusion DNA Polymerase", 0.5, "2 U/µL", "1 units"),
            ("Nuclease-free water", 32.5, "", "to 50 µL"),
        ],
        [
            ("Initial denaturation", 98.0, 30),
            ("Denature", 98.0, 10),
            ("Anneal", 55.0, 20),
            ("Extend", 72.0, 30),
            ("Final extension", 72.0, 300),
            ("Hold", 4.0, None),
        ],
        72.0,
    ),
    "Taq": (
        [
            ("Standard Taq Reaction Buffer", 5.0, "10X", "1X"),
            ("dNTP mix", 1.0, "10 mM each", "200 µM each"),
            ("Forward primer", 1.0, "10 µM", "200 nM"),
            ("Reverse primer", 1.0, "10 µM", "200 nM"),
            ("Template DNA", 1.0, "", ""),
            ("Taq DNA Polymerase", 0.25, "5 U/µL", "1.25 units"),
            ("Nuclease-free water", 40.75, "", "to 50 µL"),
        ],
        [
            ("Initial denaturation", 95.0, 30),
            ("Denature", 95.0, 30),
            ("Anneal", 55.0, 30),
            ("Extend", 68.0, 120),
            ("Final extension", 68.0, 300),
            ("Hold", 4.0, None),
        ],
        65.1,
    ),
    "OneTaq": (
        [
            ("OneTaq Standard Reaction Buffer", 10.0, "5X", "1X"),
            ("dNTP mix", 1.0, "10 mM each", "200 µM each"),
            ("Forward primer", 1.0, "10 µM", "200 nM"),
            ("Reverse primer", 1.0, "10 µM", "200 nM"),
            ("Template DNA", 1.0, "", ""),
            ("OneTaq DNA Polymerase", 0.25, "5 U/µL", "1.25 units"),
            ("Nuclease-free water", 35.75, "", "to 50 µL"),
        ],
        [
            ("Initial denaturation", 94.0, 30),
            ("Denature", 94.0, 30),
            ("Anneal", 55.0, 30),
            ("Extend", 68.0, 120),
            ("Final extension", 68.0, 300),
            ("Hold", 4.0, None),
        ],
        68.0,
    ),
}


@pytest.mark.parametrize("polymerase", POLYMERASES, ids=lambda one: one.name)
def test_each_shipped_polymerase_gets_nebs_reaction_and_program(polymerase: Polymerase) -> None:
    lines, program, two_step = SHIPPED_PCR[polymerase.name]
    table = pcr_reaction(polymerase)
    assert [(one.name, one.volume_ul, one.stock, one.final) for one in table.components] == lines
    cycled = pcr_program(polymerase, annealing_temperature=55.0, amplicon_length=1500)
    assert [stage.cycles for stage in cycled.stages] == [1, 30, 1, 1]
    assert [
        (one.label, one.temperature_c, one.seconds)
        for stage in cycled.stages
        for one in stage.incubations
    ] == program
    for annealing, steps in ((round(two_step - 0.1, 1), 3), (two_step, 2)):
        split = pcr_program(polymerase, annealing_temperature=annealing, amplicon_length=1500)
        assert len(split.stages[1].incubations) == steps


def test_the_colony_pcr_reaction_is_half_master_mix() -> None:
    table = colony_pcr_reaction()
    assert total(table) == pytest.approx(25.0)
    assert volumes(table)[COLONY_PCR_MASTER_MIX] == 12.5
    assert volumes(table)["Forward primer"] == 0.5
    assert volumes(table)["Nuclease-free water"] == 11.5


def test_the_colony_pcr_program_lyses_the_cells_and_holds_at_ten() -> None:
    program = colony_pcr_program(annealing_temperature=51.5, amplicon_length=137)
    lysis = program.stages[0].incubations[0]
    assert (lysis.label, lysis.temperature_c, lysis.seconds) == ("Lysis", 94.0, 300)
    assert program.stages[1].incubations[2].seconds == 60
    assert program.stages[-1].incubations[0].temperature_c == 10.0

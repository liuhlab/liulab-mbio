"""The shared protocol steps, each built from plain facts with no assembly plan."""

from liulab_mbio.bench.steps import dpni_step, pcr_step
from liulab_mbio.primers import Q5


def test_a_pcr_step_is_built_from_a_parts_name_and_its_reaction() -> None:
    step = pcr_step(
        "GFP",
        "GFP plasmid",
        749,
        polymerase=Q5,
        annealing_temperature=61.0,
        extension_seconds=20,
        cycles=20,
    )

    assert step.title == "Amplify GFP"
    assert step.instructions[-1] == (
        "Run the program below: 61 °C annealing and 20 s extension for a 749 bp product."
    )
    assert [table.title for table in step.tables] == ["GFP PCR"]
    (program,) = step.programs
    assert program.title == "GFP PCR"
    cycled = program.stages[1]
    assert cycled.cycles == 20
    assert [(one.label, one.temperature_c) for one in cycled.incubations][1] == ("Anneal", 61.0)
    assert step.expected == ("One band at 749 bp.",)
    assert "GFP plasmid template" in step.troubleshooting[0].solution
    assert step.notes == ()


def test_a_pcr_step_takes_the_profiles_cycles_unless_told_and_carries_the_callers_notes() -> None:
    step = pcr_step(
        "GFP",
        "GFP plasmid",
        749,
        polymerase=Q5,
        annealing_temperature=61.0,
        extension_seconds=20,
        notes=("Why this PCR is special.",),
    )

    assert step.programs[0].stages[1].cycles == 30
    assert step.notes == ("Why this PCR is special.",)


def test_a_callers_note_sits_after_the_steps_own() -> None:
    step = dpni_step(
        ("backbone",), (("pUC19", 19),), notes=("This pipeline's own reason for the digest.",)
    )

    assert step.instructions[0] == "Add 20 units of DpnI to the backbone PCR and mix."
    assert step.notes == (
        "DpnI cuts GATC only where Dam has methylated it, so it cuts pUC19 (19 Dam sites) and "
        "leaves the PCR product, which carries no methylation.",
        "This pipeline's own reason for the digest.",
    )

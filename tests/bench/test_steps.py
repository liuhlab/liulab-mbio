"""The shared protocol steps, each built from plain facts with no assembly plan."""

import pytest

from liulab_mbio.bench.phenotype import Phenotype
from liulab_mbio.bench.steps import dpni_step, pcr_step, phenotype_sentences
from liulab_mbio.primers import Q5
from liulab_mbio.sequence import Feature, Segment, Strand


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


@pytest.mark.parametrize(
    ("annotated", "sentence"),
    [
        (True, "There is a ribosome binding site annotated ahead of GFP, so the clone may make "),
        (False, "There is no ribosome binding site annotated ahead of GFP, so the clone is not "),
    ],
)
def test_the_ribosome_binding_site_reads_as_english_either_way(
    annotated: bool, sentence: str
) -> None:
    coding = Feature("GFP", "CDS", (Segment(100, 800),), strand=Strand.FORWARD)
    promoter = Feature("lac promoter", "promoter", (Segment(0, 50),), strand=Strand.FORWARD)
    phenotype = Phenotype((100, 800), coding, promoter, 50, True, annotated, None, None)

    assert phenotype_sentences(phenotype, ["GFP"])[1].startswith(sentence)

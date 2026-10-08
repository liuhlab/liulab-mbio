"""The shared protocol steps, each built from plain facts with no assembly plan."""

import pytest

from liulab_mbio.bench.pcr import cycle_citation
from liulab_mbio.bench.phenotype import Phenotype
from liulab_mbio.bench.plates import primer_plates
from liulab_mbio.bench.steps import (
    WORKING_PLATE_SINGLE_USE,
    dpni_step,
    pcr_step,
    phenotype_sentences,
    primer_plate_protocol,
    primer_plate_steps,
)
from liulab_mbio.primers import Q5
from liulab_mbio.protocol.model import Citation, Oligo, Protocol
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
        cycles_citation=cycle_citation(Q5),
    )

    assert step.title == "Amplify GFP"
    # Keyed by the amplicon and not by the title, so rewording the title keeps the step's mark.
    assert step.key == "pcr-gfp"
    assert step.instructions[-1] == (
        "Run the program below: 61 °C annealing and 20 s extension for a 749 bp product."
    )
    assert [table.title for table in step.tables] == ["GFP PCR"]
    (program,) = step.programs
    assert program.title == "GFP PCR"
    cycled = program.stages[1]
    assert cycled.cycles == 20
    assert cycled.citation == Citation("M0491", "thermocycling conditions")
    assert [(one.label, one.temperature_c) for one in cycled.incubations][1] == ("Anneal", 61.0)
    assert step.expected == ("One band at 749 bp.",)
    assert "GFP plasmid template" in step.troubleshooting[0].solution
    assert step.notes == ()


def test_a_pcr_step_leaves_an_unsourced_count_blank_and_carries_the_callers_notes() -> None:
    step = pcr_step(
        "GFP",
        "GFP plasmid",
        749,
        polymerase=Q5,
        annealing_temperature=61.0,
        extension_seconds=20,
        cycles=None,
        cycles_citation=None,
        notes=("Why this PCR is special.",),
    )

    assert step.programs[0].stages[1].cycles is None
    assert step.programs[0].stages[1].citation is None
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


def _index_oligos(count: int = 3) -> tuple[Oligo, ...]:
    """Routine primers to order once, as an order sheet's rows."""
    return tuple(Oligo(f"IDX{n}", "ACGTACGTACGTACGTACGT") for n in range(1, count + 1))


def _plate_protocol(copies: int = 2) -> Protocol:
    return primer_plate_protocol(
        _index_oligos(),
        wells=96,
        copies=copies,
        nanomoles=25.0,
        stock_um=100.0,
        working_um=1.0,
        working_ul=100.0,
    )


def test_a_primer_plate_protocol_produces_one_item_per_plate() -> None:
    one = _plate_protocol()

    assert one.consumes == ()
    assert [item.name for item in one.produces] == [
        "primer stock plate",
        "primer working plate 1",
        "primer working plate 2",
    ]
    assert one.produces[0].spec == ("100 µM",)
    assert one.produces[1].spec[0] == "1 µM"
    assert {item.storage for item in one.produces} == {"-20 °C"}
    assert [plate.name for plate in one.plates] == [item.name for item in one.produces]


def test_a_working_plate_carries_the_rule_that_it_is_never_put_back() -> None:
    """It reaches both steps that name a plate: the one splitting them and the one storing them."""
    one = _plate_protocol()
    split, stored = one.steps[-2], one.steps[-1]

    assert [(carrier, rule.subject) for carrier, rule in one.rules_for(stored)] == [
        ("primer working plate 1", "return to the freezer"),
        ("primer working plate 2", "return to the freezer"),
    ]
    assert [carrier for carrier, _ in one.rules_for(split)] == [
        "primer working plate 1",
        "primer working plate 2",
    ]
    assert {rule.kind for _, rule in one.rules_for(stored)} == {"forbids"}
    assert [check.status for check in one.audit()] == ["pass", "pass", "pass", "pass"]

    # The plate carries it, because the plate is what this protocol makes: a materials table
    # listing an output would have a reader order it.
    assert "primer working plate 1" not in [material.name for material in one.materials]
    made = next(plate for plate in one.plates if plate.name == "primer working plate 1")
    assert made.rules == (WORKING_PLATE_SINGLE_USE,)


def test_the_primers_ordered_once_are_the_protocols_own_order_sheet() -> None:
    one = _plate_protocol()

    assert [oligo.name for oligo in one.oligos] == ["IDX1", "IDX2", "IDX3"]
    assert one.steps[0].title == "Order the primers"


def test_the_split_moves_every_stock_well_into_the_same_well_of_each_copy() -> None:
    one = _plate_protocol()
    (split,) = [step for step in one.steps if step.transfers]

    assert [transfer.title for transfer in split.transfers] == [
        "Stock into primer working plate 1",
        "Stock into primer working plate 2",
    ]
    moves = split.transfers[0].moves
    assert len(moves) == 3
    assert all(move.source.well == move.destination.well for move in moves)
    assert {move.volume_ul for move in moves} == {1.0}
    assert "99 µL of nuclease-free water" in split.instructions[0]


def test_the_resuspension_volume_reaches_the_stock_concentration() -> None:
    made = primer_plates(["IDX1"], wells=96, copies=1)
    steps = primer_plate_steps(
        made, nanomoles=25.0, stock_um=100.0, working_um=1.0, working_ul=100.0, diluent="water"
    )

    assert "250 µL of water" in steps[1].instructions[0]


def test_a_primer_plate_protocol_refuses_a_working_plate_no_weaker_than_the_stock() -> None:
    with pytest.raises(ValueError, match="dilutes nothing"):
        primer_plate_protocol(
            _index_oligos(),
            wells=96,
            copies=1,
            nanomoles=25.0,
            stock_um=10.0,
            working_um=10.0,
            working_ul=100.0,
        )

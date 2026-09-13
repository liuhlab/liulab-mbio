from collections.abc import Callable
from pathlib import Path

import pytest

from liulab_mbio.protocol import (
    OVERVIEW_CHARS,
    Check,
    Component,
    Gel,
    Incubation,
    Ladder,
    Lane,
    Protocol,
    ReactionTable,
    Reference,
    Stage,
    Step,
    ThermocyclerProgram,
    Timer,
    read_protocol,
)

EXAMPLE = Path(__file__).parent / "data" / "pcr-protocol.json"


def test_an_unknown_key_is_refused_and_located() -> None:
    data = {"title": "t", "steps": [{"title": "s", "instruction": ["typo"]}]}
    with pytest.raises(ValueError, match=r"steps\[0\].*instruction"):
        Protocol.from_dict(data)


def test_a_missing_required_key_is_refused_and_located() -> None:
    data = {"title": "t", "steps": [{"title": "s", "tables": [{"components": [{"name": "x"}]}]}]}
    with pytest.raises(ValueError, match=r"components\[0\].*volume_ul"):
        Protocol.from_dict(data)


def test_the_master_mix_scales_by_reaction_count_with_overage() -> None:
    table = ReactionTable(
        (Component("Buffer", 2.5), Component("Template", 1, master_mix=False)),
        overage=0.1,
    )
    # 2.5 µL x 8 reactions x 1.1 = 22 µL; template goes into each tube, not the mix.
    assert table.mix_volumes(8) == (22.0, None)
    assert table.mix_volumes(1) == (2.75, None)


def test_a_program_duration_counts_cycles_and_skips_an_indefinite_hold() -> None:
    program = ThermocyclerProgram(
        (
            Stage((Incubation("Denature", 98, 30),)),
            Stage((Incubation("Denature", 98, 10), Incubation("Anneal", 60, 20)), cycles=30),
            Stage((Incubation("Hold", 4, None),)),
        )
    )
    assert program.duration_seconds == 30 + 30 * (10 + 20)


def test_gel_migration_is_linear_in_the_log_of_band_size() -> None:
    gel = Gel(Ladder("marker", (100, 1000, 10000)), (Lane("sample", (500,)),))
    assert gel.migration(10000) == pytest.approx(0.0)
    assert gel.migration(1000) == pytest.approx(0.5)
    assert gel.migration(100) == pytest.approx(1.0)
    assert gel.migration(500) > gel.migration(1000)


def test_gel_migration_spans_sample_bands_beyond_the_ladder() -> None:
    gel = Gel(Ladder("marker", (1000,)), (Lane("sample", (100, 10000)),))
    assert gel.migration(1000) == pytest.approx(0.5)


@pytest.mark.parametrize(
    "build",
    [
        lambda: Component("water", 0),
        lambda: ReactionTable(()),
        lambda: ReactionTable((Component("water", 1),), reactions=0),
        lambda: ReactionTable((Component("water", 1),), overage=-0.1),
        lambda: Stage((Incubation("x", 37, 60),), cycles=0),
        lambda: Incubation("x", 37, 0),
        lambda: Timer("x", 0),
        lambda: Ladder("marker", ()),
        lambda: Lane("sample", (0,)),
        lambda: Reference("x", url="javascript:alert(1)"),
        lambda: Step(""),
        lambda: Protocol(""),
    ],
)
def test_the_model_refuses_values_no_bench_could_follow(build: Callable[[], object]) -> None:
    with pytest.raises(ValueError):  # noqa: PT011
        build()


def test_master_mix_volumes_round_to_a_hundredth_of_a_microlitre() -> None:
    table = ReactionTable((Component("Polymerase", 0.33),), overage=0.1)
    # 0.33 x 2 x 1.1 = 0.726
    assert table.mix_volumes(2) == (0.73,)


def test_a_sentence_is_refused_a_place_in_the_card_grid() -> None:
    sentence = "The insert reads on the opposite strand, so that promoter does not transcribe it."
    assert len(sentence) > OVERVIEW_CHARS
    with pytest.raises(ValueError, match="highlights"):
        Protocol("t", overview={"Orientation": sentence})


def test_a_fact_of_a_few_words_is_a_card() -> None:
    assert Protocol("t", overview={"Vector": "pUC19, 2686 bp"}).overview["Vector"]


def test_a_verdict_outside_the_three_is_refused() -> None:
    with pytest.raises(ValueError, match="status"):
        Protocol.from_dict({"title": "t", "checks": [{"name": "junctions", "status": "ok"}]})


def test_a_protocol_reads_from_its_json_file() -> None:
    protocol = read_protocol(EXAMPLE)
    assert protocol.title == "Colony check by PCR"
    assert protocol.overview["Expected product"] == "500 bp"
    assert protocol.highlights[0].startswith("The reaction is a colony check")
    assert protocol.checks == (
        Check("primers", "pass", detail="2 designed, 0 with a warning"),
        Check("controls", "warn", detail="no positive control is set up"),
    )
    assert [m.name for m in protocol.materials][3:5] == ["M13 fwd", "M13 rev"]
    assert protocol.materials[3].sequence == "GTAAAACGACGGCCAGT"
    assert [s.title for s in protocol.steps] == [
        "Set up the PCR",
        "Run the thermocycler",
        "Check the product on a gel",
    ]
    table = protocol.steps[0].tables[0]
    assert table.reactions == 4
    assert table.components[-1].name == "Template DNA"
    assert table.components[-1].master_mix is False
    stage = protocol.steps[1].programs[0].stages[1]
    assert stage.cycles == 30
    assert stage.incubations[2].temperature_c == 68
    assert protocol.steps[1].programs[0].stages[-1].incubations[0].seconds is None
    gel = protocol.steps[2].gels[0]
    assert gel.ladder.name == "1 kb ladder"
    assert gel.lanes[1].bands_bp == ()
    assert protocol.steps[2].troubleshooting[0].problem == "No band"
    assert protocol.references[0].url == "https://example.org/pcr"

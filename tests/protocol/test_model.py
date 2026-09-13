import re
from collections.abc import Callable, Mapping
from pathlib import Path

import pytest

from liulab_mbio.goldengate import Plan
from liulab_mbio.protocol import (
    OVERVIEW_CHARS,
    Check,
    Component,
    Gel,
    Incubation,
    Ladder,
    Lane,
    Oligo,
    Protocol,
    ReactionTable,
    Reference,
    Stage,
    Step,
    ThermocyclerProgram,
    Timer,
    read_protocol,
    write_protocol,
)


def test_an_unknown_key_is_refused_and_located() -> None:
    data = {"title": "t", "steps": [{"title": "s", "instruction": ["typo"]}]}
    with pytest.raises(ValueError, match=r"steps\[0\].*instruction"):
        Protocol.from_dict(data)


def test_a_missing_required_key_is_refused_and_located() -> None:
    data = {"title": "t", "steps": [{"title": "s", "tables": [{"components": [{"name": "x"}]}]}]}
    with pytest.raises(ValueError, match=r"components\[0\].*volume_ul"):
        Protocol.from_dict(data)


def _one_table(
    *,
    step: Mapping[str, object] | None = None,
    table: Mapping[str, object] | None = None,
    component: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """A protocol of one step holding one reaction table of one component, with values changed."""
    one_component = {"name": "water", "volume_ul": 1, **(component or {})}
    one_table = {"components": [one_component], **(table or {})}
    return {"title": "t", "steps": [{"title": "s", "tables": [one_table], **(step or {})}]}


STEP = "protocol.steps[0]"
TABLE = f"{STEP}.tables[0]"
COMPONENT = f"{TABLE}.components[0]"


@pytest.mark.parametrize(
    ("data", "where", "expected"),
    [
        (_one_table(step={"instructions": "Mix well."}), f"{STEP}.instructions", "a list"),
        (_one_table(component={"volume_ul": "5"}), f"{COMPONENT}.volume_ul", "a number"),
        (_one_table(component={"volume_ul": True}), f"{COMPONENT}.volume_ul", "a number"),
        (_one_table(component={"volume_ul": None}), f"{COMPONENT}.volume_ul", "a number"),
        (_one_table(component={"master_mix": "false"}), f"{COMPONENT}.master_mix", "true or false"),
        (_one_table(table={"reactions": 2.5}), f"{TABLE}.reactions", "a whole number"),
        (_one_table(step={"title": 5}), f"{STEP}.title", "a string"),
        ({"title": "t", "overview": ["Vector"]}, "protocol.overview", "an object"),
    ],
)
def test_a_value_of_the_wrong_type_is_refused_and_located(
    data: Mapping[str, object], where: str, expected: str
) -> None:
    with pytest.raises(ValueError, match=f"^{re.escape(where)}: expected {expected}"):
        Protocol.from_dict(data)


def test_a_whole_number_is_a_measurement() -> None:
    protocol = Protocol.from_dict(_one_table(component={"volume_ul": 5}))
    assert protocol.steps[0].tables[0].components[0].volume_ul == 5


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
        lambda: Oligo("M13 fwd", "  "),
        lambda: Oligo("M13 fwd", "GTAAAACG", checks=(Check("length", "warn", "17"),)),
        lambda: Oligo(
            "M13 fwd", "GTAAAACG", status="pass", checks=(Check("length", "warn", "17"),)
        ),
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


def test_an_oligo_carries_its_own_verdict_and_the_checks_that_fired() -> None:
    row = Oligo(
        "Sequencing forward",
        "AACTGTTGGGAAGGGC",
        status="warn",
        checks=(
            Check("length", "warn", "16 (band 18-30)"),
            Check("GC clamp", "warn", "4 (band 1-3)"),
        ),
    )
    assert row.status == "warn"
    assert [(check.name, check.detail) for check in row.checks] == [
        ("length", "16 (band 18-30)"),
        ("GC clamp", "4 (band 1-3)"),
    ]


def test_an_oligo_nothing_judged_says_so_rather_than_reading_as_a_pass() -> None:
    assert Oligo("M13 fwd", "GTAAAACGACGGCCAGT").status is None
    assert Oligo("M13 fwd", "GTAAAACGACGGCCAGT", status="pass").checks == ()


def test_a_verdict_outside_the_three_is_refused() -> None:
    with pytest.raises(ValueError, match="status"):
        Protocol.from_dict({"title": "t", "checks": [{"name": "junctions", "status": "ok"}]})
    with pytest.raises(ValueError, match="status"):
        Protocol.from_dict(
            {"title": "t", "oligos": [{"name": "M13 fwd", "sequence": "ACGT", "status": "ok"}]}
        )


def test_a_protocol_reads_from_its_json_file(data_dir: Path) -> None:
    protocol = read_protocol(data_dir / "pcr-protocol.json")
    assert protocol.title == "Colony check by PCR"
    assert protocol.overview["Expected product"] == "500 bp"
    assert protocol.highlights[0].startswith("The reaction is a colony check")
    assert protocol.checks == (
        Check(
            "primers",
            "warn",
            detail="2 designed, 2 with a warning, 0 failing: length on 2; "
            "not judged: full-primer Tm, 3' end stability",
        ),
        Check("controls", "warn", detail="no positive control is set up"),
    )
    # Oligos are their own list, and a material the file gives no catalogue number keeps none.
    assert [m.name for m in protocol.materials][:1] == ["Taq DNA Polymerase"]
    assert (protocol.materials[0].supplier, protocol.materials[0].catalog) == ("NEB", "M0273")
    assert (protocol.materials[-1].supplier, protocol.materials[-1].catalog) == ("", "")
    assert [o.name for o in protocol.oligos] == ["M13 fwd", "M13 rev"]
    assert protocol.oligos[0] == Oligo(
        "M13 fwd",
        "GTAAAACGACGGCCAGT",
        purpose="Set up the PCR",
        tm_c=55.4,
        stock="10 µM",
        status="warn",
        checks=(Check("length", "warn", "17 (band 18-30)"),),
    )
    assert protocol.equipment[0] == "Thermocycler with a heated lid"
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


def test_a_protocol_written_as_json_reads_back_equal(
    data_dir: Path, plan: Plan, tmp_path: Path
) -> None:
    for name, protocol in (
        ("example", read_protocol(data_dir / "pcr-protocol.json")),
        ("golden-gate", plan.protocol()),
    ):
        path = write_protocol(protocol, tmp_path / f"{name}.json")
        assert read_protocol(path) == protocol


def test_every_field_is_written_in_its_declared_order_even_when_empty(tmp_path: Path) -> None:
    path = write_protocol(Protocol("Spin at 4 °C"), tmp_path / "protocol.json")
    assert path.read_bytes() == "\n".join(
        [
            "{",
            '  "title": "Spin at 4 °C",',
            '  "summary": "",',
            '  "overview": {},',
            '  "highlights": [],',
            '  "checks": [],',
            '  "materials": [],',
            '  "oligos": [],',
            '  "equipment": [],',
            '  "steps": [],',
            '  "references": []',
            "}",
            "",
        ]
    ).encode("utf-8")

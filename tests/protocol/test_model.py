import ast
import re
from collections.abc import Callable, Mapping
from pathlib import Path

import pytest

from liulab_mbio.protocol import (
    OVERVIEW_CHARS,
    Check,
    Citation,
    Component,
    Figure,
    Gel,
    Incubation,
    Ladder,
    Lane,
    Oligo,
    Project,
    Protocol,
    ReactionTable,
    Reference,
    Stage,
    Step,
    ThermocyclerProgram,
    Timer,
    Wait,
    read_protocol,
    write_protocol,
)

#: The checkout a `Source.note` is relative to.
REPO = Path(__file__).resolve().parents[2]


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


def test_a_blank_cycle_count_leaves_the_run_time_unknown() -> None:
    program = ThermocyclerProgram(
        (
            Stage((Incubation("Denature", 98, 30),)),
            Stage((Incubation("Denature", 98, 10),), cycles=None),
        )
    )
    assert program.duration_seconds is None


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
        lambda: Figure((), "The product"),
        lambda: Figure(("a.dna",), " "),
        lambda: Figure(("a.dna",), "The product", span=(400, 100)),
        lambda: Step(""),
        lambda: Protocol(""),
        lambda: Project(""),
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


def test_a_protocol_written_as_json_reads_back_equal(data_dir: Path, tmp_path: Path) -> None:
    protocol = read_protocol(data_dir / "pcr-protocol.json")
    path = write_protocol(protocol, tmp_path / "example.json")
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
            '  "consumes": [],',
            '  "produces": [],',
            '  "materials": [],',
            '  "oligos": [],',
            '  "order_sheet": "",',
            '  "equipment": [],',
            '  "vessels": [],',
            '  "plates": [],',
            '  "steps": [],',
            '  "references": [],',
            '  "sources": {},',
            '  "holes": [],',
            '  "bill": null',
            "}",
            "",
        ]
    ).encode("utf-8")


def test_a_step_waits_on_a_vendor_nobody_attends() -> None:
    wait = Wait("the oligo pool to arrive", duration="10-15 working days")
    step = Step("Order the oligo pool", waits=(wait,))
    assert step.waits == (wait,)
    # A turnaround nobody will state is admitted, never written as a zero.
    assert Wait("sequencing to come back").duration == ""


def test_a_wait_that_says_nothing_is_waited_on_is_refused() -> None:
    with pytest.raises(ValueError, match="waits on"):
        Wait(" ")


def test_the_hands_on_share_of_a_step_is_unknown_until_it_is_stated() -> None:
    assert Step("Thaw the cells").hands_on_seconds is None
    assert Step("Thaw the cells", hands_on_seconds=0).hands_on_seconds == 0
    with pytest.raises(ValueError, match="hands_on_seconds"):
        Step("Thaw the cells", hands_on_seconds=-1)


def test_a_wait_cites_its_turnaround_like_any_other_row() -> None:
    cited = Citation("vendor")
    one = Protocol("t", steps=(Step("Order it", waits=(Wait("the pool", citation=cited),)),))
    assert one.cited == frozenset({"vendor"})
    (check,) = [c for c in one.audit() if c.name == "sources"]
    assert check.status == "fail"


def test_a_steps_time_round_trips_through_json(tmp_path: Path) -> None:
    one = Protocol(
        "t",
        steps=(
            Step(
                "Order the oligo pool",
                waits=(Wait("the pool to arrive", duration="10-15 working days"),),
                hands_on_seconds=900,
            ),
        ),
    )
    assert read_protocol(write_protocol(one, tmp_path / "protocol.json")) == one


def test_a_figure_names_its_record_by_path_and_round_trips_through_json(tmp_path: Path) -> None:
    one = Protocol(
        "t",
        steps=(
            Step(
                "Assemble the vector and the insert",
                figures=(
                    Figure(
                        ("product.dna",),
                        "The assembled plasmid, opened at the first junction",
                        span=(2683, 2689),
                        linear=True,
                        sequence_view=True,
                        enzymes=("BsaI",),
                        highlight=("GFP",),
                    ),
                ),
            ),
        ),
    )
    assert read_protocol(write_protocol(one, tmp_path / "protocol.json")) == one


def test_a_figure_cites_its_source_as_a_note_does() -> None:
    one = Protocol(
        "t",
        steps=(
            Step("Draw it", figures=(Figure(("v.dna",), "The vector", citation=Citation("k")),)),
        ),
    )
    assert one.cited == frozenset({"k"})
    (check,) = [c for c in one.audit() if c.name == "sources"]
    assert check.status == "fail"


def _notes() -> set[str]:
    """Every note a `Source(...)` call in `src/` spells as a literal `note=`.

    Read off the source text, so a source built inside a function counts; a path a function
    computes does not, and no scan would see one.
    """
    found: set[str] = set()
    for path in (REPO / "src").rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        if "Source(" not in text:
            continue
        for node in ast.walk(ast.parse(text)):
            if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name):
                continue
            if node.func.id != "Source":
                continue
            for keyword in node.keywords:
                if keyword.arg == "note" and isinstance(keyword.value, ast.Constant):
                    found.add(str(keyword.value.value))
    return found


def test_every_note_a_source_names_is_on_disk() -> None:
    """`Source.note` is a repo-relative path, so a note that moves or goes fails here.

    A source naming no note is not a failure, so naming none anywhere is not one either.
    """
    assert sorted(note for note in _notes() if not (REPO / note).is_file()) == []


def test_a_span_that_is_not_a_pair_of_numbers_is_refused_where_it_stands() -> None:
    figure = {"records": ["v.dna"], "caption": "c", "span": [1, 2, 3]}
    data = {"title": "t", "steps": [{"title": "s", "figures": [figure]}]}
    with pytest.raises(ValueError, match=r"figures\[0\]\.span: expected a list of 2, got 3"):
        Protocol.from_dict(data)

import re
from dataclasses import replace
from pathlib import Path

import pytest

from mbio.bench import plates
from mbio.plot import layers
from mbio.protocol import (
    OVERVIEW_CHARS,
    AmountToVolume,
    Caution,
    Check,
    Citation,
    Component,
    Figure,
    Folder,
    Incubation,
    Item,
    Material,
    Move,
    Note,
    Oligo,
    Page,
    Plate,
    Protocol,
    ReactionTable,
    Source,
    Stage,
    Step,
    ThermocyclerProgram,
    Timer,
    Transfer,
    Troubleshooting,
    Well,
    read_protocol,
    render_html,
    write_html,
)
from mbio.protocol.render import HIGHLIGHTS_HEADING, NO_NUMBER, page_key

from ..html import Node, parse


@pytest.fixture(scope="module")
def page_html(data_dir: Path) -> str:
    return render_html(read_protocol(data_dir / "pcr-protocol.json"))


@pytest.fixture(scope="module")
def page(page_html: str) -> Node:
    return parse(page_html)


def test_the_page_loads_nothing_from_outside_itself(page_html: str, page: Node) -> None:
    assert not page.find_all("link")
    assert not [n for n in page.iter() if "src" in n.attrs]
    assert "@import" not in page_html
    links = {a.attrs["href"] for a in page.find_all("a")}
    # The only addresses on the page are links a reader may choose to follow.
    assert set(re.findall(r"https?://[^\s\"'<>]+", page_html)) <= links


def test_text_from_the_protocol_is_escaped() -> None:
    protocol = Protocol(
        '<script>alert("t")</script>',
        summary="A & B",
        oligos=(Oligo("oligo", 'AC"GT'),),
        steps=(Step("<b>mix</b>", instructions=("5 µL < 10 µL",)),),
    )
    html = render_html(protocol)
    page = parse(html)
    assert len(page.find_all("script")) == 1
    assert not page.find_all("b")
    assert page.find_all("h1")[0].text == '<script>alert("t")</script>'
    assert "5 µL &lt; 10 µL" in html
    assert page.find_all("button", cls="copy")[0].attrs["data-copy"] == 'AC"GT'


def test_the_header_reads_in_bands_and_its_sentences_carry_a_heading(page: Node) -> None:
    header = page.find_all("header", cls="intro")[0]
    blocks = [n.attrs.get("class") or n.tag for n in header.children if isinstance(n, Node)]
    assert blocks == ["h1", "summary", "overview", "highlights", "status", "toolbar"]
    assert header.find_all(cls="highlights")[0].find_all("h3")[0].text == HIGHLIGHTS_HEADING


def test_a_card_holds_a_fact_and_a_sentence_is_prose(page: Node) -> None:
    cards = page.find_all("dl", cls="overview")[0].find_all("div")
    assert [card.find_all("dt")[0].text for card in cards] == ["Template", "Expected product"]
    assert all(len(card.find_all("dd")[0].text) <= OVERVIEW_CHARS for card in cards)
    prose = [p.text for p in page.find_all(cls="highlights")[0].find_all("p")]
    assert len(prose) == 2
    assert prose[0].startswith("The reaction is a colony check")
    assert prose[0] not in page.find_all("dl", cls="overview")[0].text


def test_each_check_is_a_badge_and_a_warn_shows_without_being_read(page: Node) -> None:
    badges = page.find_all("li", cls="check")
    assert [badge.attrs["class"] for badge in badges] == ["check is-warn", "check is-warn"]
    assert [badge.find_all(cls="verdict")[0].text for badge in badges] == ["warn", "warn"]
    assert badges[1].find_all(cls="check-name")[0].text == "controls"
    # Only a verdict that is not a pass spells its detail out, and it names the kinds.
    details = [p.text for p in page.find_all(cls="check-detail")]
    assert details[0].startswith("primers warn: 2 designed, 2 with a warning")
    assert "0 failing: length on 2" in details[0]
    assert details[1] == "controls warn: no positive control is set up"


def test_a_check_no_threshold_judges_shows_as_unjudged_and_never_as_a_pass() -> None:
    protocol = Protocol("Plan", checks=(Check("buffer", None, detail="look the pair up"),))
    page = parse(render_html(protocol))
    badge = page.find_all("li", cls="check")[0]
    assert badge.attrs["class"] == "check is-none"
    assert badge.find_all(cls="verdict")[0].text == "not judged"
    assert page.find_all(cls="check-detail")[0].text == "buffer not judged: look the pair up"


def test_every_step_and_instruction_has_its_own_checkbox(page: Node) -> None:
    """`protocol.js` ticks a step from the boxes in its instructions list, and they from it."""
    steps = page.find_all("section", cls="step")
    assert [len(s.find_all("input", type="checkbox")) for s in steps] == [1 + 3, 1 + 1, 1 + 2]
    assert [len(s.find_all("input", cls="done")) for s in steps] == [1, 1, 1]
    lists = [s.find_all("ol", cls="instructions")[0] for s in steps]
    assert [len(one.find_all("input", type="checkbox")) for one in lists] == [3, 1, 2]
    keys = [box.attrs["data-key"] for box in page.find_all("input", type="checkbox")]
    assert len(set(keys)) == len(keys)


def reworded() -> tuple[Protocol, Protocol]:
    """One protocol and the same one with every step's title rewritten, as an agent edits it."""
    before = Protocol(
        "Assemble",
        key="assemble",
        steps=(
            Step("Set the reaction up", key="set-up", instructions=("Thaw the mix.",)),
            Step("Run the thermocycler", key="cycle", timers=(Timer("ligate", 60),)),
        ),
    )
    after = replace(
        before,
        steps=(
            replace(before.steps[0], title="Set up the Golden Gate reaction"),
            replace(before.steps[1], title="Cycle it"),
        ),
    )
    return before, after


def test_a_step_is_addressed_by_its_key_so_a_reworded_title_keeps_its_ticks() -> None:
    """ADR 0002's own case: an agent rewords the steps, renders again, nothing is ticked twice."""
    marks = [
        sorted(
            box.attrs["data-key"]
            for box in parse(render_html(one)).find_all("input", type="checkbox")
        )
        for one in reworded()
    ]
    assert marks[0] == ["step-cycle", "step-set-up", "step-set-up.1"]
    assert marks[1] == marks[0]


def test_a_page_remembers_under_the_key_its_protocol_carries() -> None:
    before, after = reworded()
    stores = [
        parse(render_html(one)).find_all("body")[0].attrs["data-protocol"]
        for one in (before, after)
    ]
    assert stores == ["assemble", "assemble"]
    # A protocol nobody has keyed falls back to a digest of its content, which an edit changes.
    unkeyed = [replace(one, key="") for one in (before, after)]
    assert page_key(unkeyed[0]) != page_key(unkeyed[1])


def test_no_two_marks_of_one_page_are_alike_however_its_steps_are_keyed() -> None:
    """A duplicate is a defect the package's own tests catch; the page renders regardless.

    The keys below are the ones that collide where a repeat is numbered once and left: ``cut``
    twice beside a step already called ``cut-2``.
    """
    protocol = Protocol(
        "Digest",
        steps=(
            Step("Digest", key="cut"),
            Step("Digest again", key="cut"),
            Step("Digest once more", key="cut-2"),
            Step("Set up", key="set-up", instructions=("Thaw the mix.",)),
            Step("Set up again", key="set-up-1"),
        ),
    )
    page = parse(render_html(protocol))
    assert [one.attrs["id"] for one in page.find_all("section", cls="step")] == [
        "step-cut",
        "step-cut-2",
        "step-cut-2-3",
        "step-set-up",
        "step-set-up-1",
    ]
    marks = [box.attrs["data-key"] for box in page.find_all("input", type="checkbox")]
    assert len(set(marks)) == len(marks)


def rounds() -> Protocol:
    """A protocol of two like rounds between a step that opens it and one that closes it."""
    return Protocol(
        "Assemble in rounds",
        steps=(
            Step("Pool each part list", key="pool"),
            *(
                Step(f"Round {n}: {what}", key=f"round-{n}-{what}", section=f"Round {n}")
                for n in (1, 2)
                for what in ("open", "ligate")
            ),
            Step("Read the library back", key="read", section="Read the library back"),
        ),
    )


def test_the_navigation_groups_the_steps_under_the_section_each_belongs_to() -> None:
    [nav] = parse(render_html(rounds())).find_all("nav", cls="within")
    groups = nav.find_all("details")
    assert [g.attrs["data-steps"] for g in groups] == [
        "round-1-open round-1-ligate",
        "round-2-open round-2-ligate",
        "read",
    ]
    counts = [g.find_all("summary")[0].text for g in groups]
    assert counts == [
        "Round 1 0 of 2 done",
        "Round 2 0 of 2 done",
        "Read the library back 0 of 1 done",
    ]
    # The bench arrives at the top of the run, so the first section is the one standing open.
    assert [("open" in g.attrs) for g in groups] == [True, False, False]
    # A step naming no section is not forced into one, and the numbering runs through every list.
    assert [a.text for a in nav.find_all("a")][:2] == ["1 Pool each part list", "2 Round 1: open"]
    assert nav.find_all("a")[-1].text == "6 Read the library back"


def test_a_page_lists_its_steps_in_the_column_beside_it_alone_or_in_a_run() -> None:
    """A page lists its steps in the column beside it; a page alone has no run's left column."""
    one = rounds()
    alone = parse(render_html(one))
    [frame] = alone.find_all("div", cls="alone")
    kinds = [n.attrs.get("class", "") for n in frame.children if isinstance(n, Node)]
    assert kinds == ["page", "column within"]
    [column] = parse(
        render_html(one, folder=Folder((Page.of(1, one),)), here="01-x.html")
    ).find_all("nav", cls="within")
    [beside] = frame.find_all("nav", cls="within")
    assert [g.attrs["data-steps"] for g in column.find_all("details")] == [
        g.attrs["data-steps"] for g in beside.find_all("details")
    ]
    # One list of steps a page, and it is the column's.
    assert len(alone.find_all("nav")) == 1


def test_steps_naming_no_section_stay_the_one_list_they_were(page: Node) -> None:
    [nav] = page.find_all("nav", cls="within")
    assert not nav.find_all("details")
    assert next(a.attrs["href"] for a in nav.find_all("a")).startswith("#step-")


def test_the_header_states_what_the_bench_is_handed_and_what_it_is_left_with() -> None:
    one = Protocol(
        "Transform",
        overview={"Strain": "DH10B"},
        consumes=(Item("ligation", "the joined plasmid", spec=("≥10 ng/µL",)),),
        produces=(Item("colonies", "transformed cells on a plate"),),
    )
    header = parse(render_html(one)).find_all("header", cls="intro")[0]
    blocks = [n.attrs.get("class") or n.tag for n in header.children if isinstance(n, Node)]
    assert blocks[:3] == ["h1", "overview", "handover"]
    [handover] = header.find_all("div", cls="handover")
    assert [h.text for h in handover.find_all("h3")] == ["Have in hand", "Leaves you with"]
    assert handover.find_all("li")[0].text == "ligation the joined plasmid · ≥10 ng/µL"
    assert handover.find_all("li")[1].text == "colonies transformed cells on a plate"


def test_a_protocol_declaring_neither_states_no_handover(page: Node) -> None:
    assert not page.find_all("div", cls="handover")


def test_a_protocol_declaring_one_side_states_that_side_alone() -> None:
    one = Protocol("Transform", produces=(Item("colonies", "transformed cells on a plate"),))
    [handover] = parse(render_html(one)).find_all("div", cls="handover")
    assert [h.text for h in handover.find_all("h3")] == ["Leaves you with"]


def test_the_toolbar_offers_to_reset_everything_the_page_remembers(page: Node) -> None:
    """One button for one store: `protocol.js` empties it, marks, counts and timers alike."""
    [button] = page.find_all("button", cls="reset")
    assert button.text == "Reset page"


def test_an_oligo_is_an_order_sheet_row_and_never_a_material(page: Node) -> None:
    materials = page.find_all("section", cls="materials")[0]
    assert "M13 fwd" not in materials.text
    assert "GTAAAACGACGGCCAGT" not in materials.text
    oligos = page.find_all("section", cls="oligos")[0]
    assert [th.text for th in oligos.find_all("th")] == [
        "Name",
        "Sequence (5'→3')",
        "Length",
        "Tm (°C)",
        "For",
        "Working stock",
        "Checks",
    ]
    row = next(r for r in oligos.find_all("tr") if "M13 fwd" in r.text)
    cells = [cell.text for cell in row.find_all("td")]
    assert cells[0] == "M13 fwd"
    assert cells[2] == str(len("GTAAAACGACGGCCAGT"))
    assert cells[3:] == ["55.4", "Set up the PCR", "10 µM", "warn"]


def test_a_warned_row_says_its_verdict_in_a_word_and_why_behind_a_toggle(page: Node) -> None:
    oligos = page.find_all("section", cls="oligos")[0]
    row = next(r for r in oligos.find_all("tr") if "M13 fwd" in r.text)
    # The word carries the verdict, so a colour fill is never what the reader has to see.
    assert row.find_all("td", cls="verdict-cell")[0].find_all(cls="verdict")[0].text == "warn"
    toggle = oligos.find_all("details")[0]
    assert toggle.find_all("summary")[0].text == "Checks on 2 of 2 oligos"
    assert [li.text for li in toggle.find_all("li")] == [
        "M13 fwd warn: length 17 (band 18-30)",
        "M13 rev warn: length 17 (band 18-30)",
    ]


def test_a_clean_row_says_pass_and_an_unjudged_one_says_so() -> None:
    protocol = Protocol(
        "t",
        oligos=(Oligo("clean", "ACGTACGTACGTACGTACGT", status="pass"), Oligo("bare", "ACGTACGT")),
    )
    page = parse(render_html(protocol))
    assert [td.text for td in page.find_all("td", cls="verdict-cell")] == ["pass", "not judged"]
    # Nothing fired, so the sheet carries no toggle at all.
    assert not page.find_all("details")


def test_a_material_with_no_catalogue_number_gets_an_empty_cell(page: Node) -> None:
    materials = page.find_all("section", cls="materials")[0]
    head = [th.text for th in materials.find_all("th")]
    column = head.index("Catalogue")
    rows = [r for r in materials.find_all("tr") if r.find_all("td")]
    numbers = [r.find_all("td")[column].text for r in rows]
    assert numbers[0] == "M0273"
    # Nothing is invented for the rest.
    assert numbers[1:] == ["", "", ""]


def test_the_equipment_is_one_light_line_under_the_materials(page: Node) -> None:
    line = page.find_all(cls="equipment")[0].text
    assert line.startswith("Equipment:")
    assert "Thermocycler with a heated lid" in line


def test_each_oligo_sequence_has_a_copy_button(page: Node) -> None:
    copies = [b.attrs["data-copy"] for b in page.find_all("button", cls="copy")]
    assert "GTAAAACGACGGCCAGT" in copies
    assert "CAGGAAACAGCTATGAC" in copies
    assert "M13 fwd\tGTAAAACGACGGCCAGT\nM13 rev\tCAGGAAACAGCTATGAC" in copies


def test_a_reaction_table_opens_scaled_to_its_reaction_count(page: Node) -> None:
    table = page.find_all(cls="reaction")[0]
    assert table.find_all("input", type="number")[0].attrs["value"] == "4"
    # Per-reaction volume x 4 reactions x 1.1 for overage, to three figures; the last cell is
    # the mix total, which says so because a component is added per tube and not to the mix.
    assert [c.text for c in table.find_all("td", cls="mix")] == [
        "87.5",
        "11",
        "2.2",
        "2.2",
        "2.2",
        "0.55",
        "106mix only",
    ]
    template = next(r for r in table.find_all("tr") if "Template DNA" in r.text)
    assert "each tube" in template.text


def test_a_measured_row_opens_at_the_concentration_its_volume_assumes() -> None:
    """The plan's numbers stand until the reader types; the warning is another step's own entry."""
    calculator = AmountToVolume(81.98, made_up_by="Water", too_dilute="Too dilute")
    row = Component("pUC19", 2.0, master_mix=False, calculator=calculator)
    measure = Step("Measure", troubleshooting=(Troubleshooting("Too dilute", "Concentrate it."),))
    protocol = Protocol(
        "t", steps=(measure, Step("Mix", tables=(ReactionTable((row, Component("Water", 13.0))),)))
    )
    figure = parse(render_html(protocol)).find_all(cls="reaction")[0]

    [measured] = figure.find_all("tr", cls="measured")
    assert (measured.attrs["data-ng-ul"], measured.attrs["data-fill"]) == ("40.99", "1")
    assert measured.find_all("input", cls="calc-value")[0].attrs["value"] == "41"
    [warning] = figure.find_all(cls="calc-warning")
    assert "hidden" in warning.attrs
    assert warning.text == "Too dilute. Concentrate it."
    # Out of the digest, so a calculator re-keys no page and the bench keeps its ticks.
    plain = ReactionTable((replace(row, calculator=None), Component("Water", 13.0)))
    assert page_key(protocol) == page_key(
        replace(protocol, steps=(measure, Step("Mix", tables=(plain,))))
    )


def volume_cells(volume_ul: float) -> list[str]:
    """Every number a one-component reaction table of `volume_ul` prints, and its dispense line."""
    table = ReactionTable((Component("Water", volume_ul),), reactions=1, overage=0.1)
    figure = parse(render_html(Protocol("t", steps=(Step("Mix", tables=(table,)),))))
    figure = figure.find_all(cls="reaction")[0]
    return [cell.text for cell in figure.find_all("td", cls="num")] + [
        figure.find_all(cls="dispense")[0].text
    ]


@pytest.mark.parametrize(
    ("volume_ul", "printed"),
    [
        (2.5e-9, "2.5 × 10⁻⁹"),
        (1e7, "1 × 10⁷"),
        (0.0009, "9 × 10⁻⁴"),
        (0.001, "0.001"),
        (999999, "999,999"),
        (1e6, "1 × 10⁶"),
        (1161.6, "1,162"),
        (198.0, "198"),
        (0.05, "0.05"),
    ],
)
def test_a_volume_prints_to_three_figures_and_far_from_the_bench_as_a_power(
    volume_ul: float, printed: str
) -> None:
    assert volume_cells(volume_ul)[0] == printed


def test_a_volume_under_half_a_thousandth_of_a_microlitre_never_prints_as_zero() -> None:
    # The mix is the volume times 1.1, the total row repeats both, and the line says it again.
    assert volume_cells(4e-4) == [
        "4 × 10⁻⁴",
        "4.4 × 10⁻⁴",
        "4 × 10⁻⁴",
        "4.4 × 10⁻⁴",
        "Put 4 × 10⁻⁴ µL of mix in each tube.",
    ]


def test_a_thermocycler_program_lists_temperatures_times_and_cycles(page: Node) -> None:
    program = page.find_all(cls="program")[0]
    rows = [[c.text for c in r.find_all(("td", "th"))] for r in program.find_all("tr")][1:]
    assert rows == [
        ["Initial denaturation", "95 °C", "30 s", "1"],
        ["Denaturation", "95 °C", "15 s", "×30"],
        ["Annealing", "55 °C", "15 s"],
        ["Extension", "68 °C", "1 min"],
        ["Final extension", "68 °C", "5 min", "1"],
        ["Hold", "4 °C", "∞", "1"],
    ]
    caption = program.find_all("figcaption")[0].text
    assert "105 °C" in caption
    assert "50 min 30 s" in caption
    # The run the caption times is the program's own timer, so no step has to add one.
    [timer] = program.find_all(cls="timer")
    assert timer.attrs["data-seconds"] == "3030.0"


def _program_rows(program: ThermocyclerProgram) -> tuple[Node, list[list[str]]]:
    """One program rendered on a page of its own, as its figure and its body rows."""
    page = parse(render_html(Protocol("Run", steps=(Step("Cycle", programs=(program,)),))))
    figure = page.find_all(cls="program")[0]
    rows = [[c.text for c in r.find_all("td")] for r in figure.find_all("tr")][1:]
    return figure, rows


def test_a_touchdown_is_one_cycled_stage_printed_from_its_start_to_its_derived_end() -> None:
    """The end follows from the step and the count, so the page cannot contradict the model."""
    cite = Citation("LevSeq", "thermal cycler table")
    anneal = Incubation("Anneal", 68.0, 20, delta_c=-0.5, citation=cite)
    figure, rows = _program_rows(
        ThermocyclerProgram(
            (
                Stage((Incubation("Denature", 95.0, 20, citation=cite), anneal), cycles=10),
                Stage((anneal,), cycles=None),
            ),
            title="Index PCR",
        )
    )
    assert rows == [
        ["Denature", "95 °C", "20 s", "×10"],
        ["Anneal", "68 → 63.5 °C-0.5 °C a cycle", "20 s"],
        ["Anneal", "68 °C-0.5 °C a cycle", "20 s", NO_NUMBER],
    ]
    # One source behind every row is the program's, so the rows carry no citation of their own.
    assert [cited.attrs["title"] for cited in figure.find_all("a", cls="cite")] == [
        "LevSeq · thermal cycler table"
    ]
    # A blank count bounds no run, so no timer stands in for one.
    assert not figure.find_all(cls="timer")


def test_a_program_held_at_one_temperature_is_not_given_ramps_it_does_not_run() -> None:
    """A water bath and an incubator are drawn as programs, and their rows take the whole time."""
    figure, rows = _program_rows(
        ThermocyclerProgram(
            (
                Stage((Incubation("Recovery, shaking", 30.0, 3600),)),
                Stage((Incubation("Outgrowth", 30.0, 43200),)),
            )
        )
    )
    assert rows == [["Recovery, shaking", "30 °C", "1 h", "1"], ["Outgrowth", "30 °C", "12 h", "1"]]
    assert "plus ramps" not in figure.find_all("figcaption")[0].text


def test_a_gel_draws_each_band_at_a_height_set_by_log_size(page: Node) -> None:
    gel = page.find_all(cls="gel")[0]
    lanes = {
        g.attrs["data-lane"]: {
            int(b.attrs["data-bp"]): float(b.attrs["y"]) for b in g.find_all("rect", cls="band")
        }
        for g in gel.find_all("g", cls="lane")
    }
    assert list(lanes) == ["1 kb ladder", "Sample", "No template"]
    ladder = lanes["1 kb ladder"]
    heights = [ladder[bp] for bp in sorted(ladder, reverse=True)]
    assert heights == sorted(set(heights))
    # 1000 bp is one decade from both 10000 and 100 bp.
    assert ladder[1000] == pytest.approx((ladder[10000] + ladder[100]) / 2, abs=0.1)
    assert lanes["Sample"] == {500: pytest.approx(ladder[500])}
    assert lanes["No template"] == {}
    legend = gel.find_all(cls="gel-legend")[0].text
    assert "Sample: 500 bp" in legend
    assert "No template: no band" in legend


def test_a_timer_starts_from_the_plans_time_in_a_field_the_reader_may_type_over(
    page: Node,
) -> None:
    """`data-seconds` keeps the plan's time, which the page offers back once the reader's differs."""
    [timer] = [one for one in page.find_all(cls="timer") if "Gel run" in one.text]
    assert timer.attrs["data-seconds"] == "1800"
    [field] = timer.find_all("input", cls="timer-time")
    assert field.attrs["value"] == "30:00"
    assert "readonly" not in field.attrs


def test_each_timer_is_keyed_so_a_running_one_survives_a_page_turn() -> None:
    """`protocol.js` keeps a deadline under this key, so a key is a timer's own and no other's."""
    ligate = ThermocyclerProgram((Stage((Incubation("ligate", 25.0, 60),)),))
    protocol = Protocol(
        "Incubate",
        steps=(
            Step("Digest", timers=(Timer("digest", 60), Timer("heat", 120))),
            Step("Ligate", programs=(ligate,), timers=(Timer("ligate", 60),)),
        ),
    )
    page = parse(render_html(protocol))
    keys = [timer.attrs["data-key"] for timer in page.find_all(cls="timer")]
    assert keys == [
        "step-digest.timer.1",
        "step-digest.timer.2",
        "step-ligate.program.1",
        "step-ligate.timer.1",
    ]


def test_a_duration_of_an_hour_or_more_is_printed_to_the_minute() -> None:
    """Nothing a bench reads in hours is planned to the second, so the seconds are not printed."""
    program = ThermocyclerProgram(
        (
            Stage(
                (Incubation("hold", 72.0, 3661), Incubation("rest", 4.0, 59)),
            ),
        )
    )
    page = parse(render_html(Protocol("Run", steps=(Step("Cycle", programs=(program,)),))))
    times = [cell.text for cell in page.find_all("td", cls="num")]
    assert "1 h 1 min" in times
    assert "59 s" in times


def test_expected_results_troubleshooting_and_references_are_shown(page: Node) -> None:
    expected = page.find_all(cls="expected")[0].text
    assert "One band at 500 bp" in expected
    assert "Band in the no-template lane" in page.find_all(cls="trouble")[0].text
    assert page.find_all("a", href="https://example.org/pcr")


def test_a_page_inlining_a_drawing_embeds_its_faces_once_and_keeps_it_on_white(page: Node) -> None:
    """The plate is laid out on white in the faces the package measured it in, in either scheme."""
    drawn = parse(render_html(Protocol("Pick", plates=(plates.plate("picked", 96),))))
    style = drawn.find_all("style")[0].text
    [figure] = drawn.find_all("figure", cls="drawing")

    families = {text.attrs["font-family"] for text in figure.find_all("text")}

    faces = re.findall(r"@font-face \{([^}]*)\}", style)

    assert families
    assert all(any(f"font-family: {family};" in face for face in faces) for family in families)
    assert len(faces) == len(set(faces))
    # One rule, under no scheme's media query, so the dark ground never shows through.
    assert style.count(".drawing svg") == 1
    [declarations] = re.findall(r"\.drawing svg \{([^}]*)\}", style)
    assert "background: #ffffff" in declarations
    # A page with no drawing carries no face: they cost more than the rest of the page.
    assert "@font-face" not in page.find_all("style")[0].text


def test_the_page_fits_a_phone_prints_and_follows_dark_mode(page: Node) -> None:
    assert page.find_all("meta", name="viewport")
    style = page.find_all("style")[0].text
    assert "@media print" in style
    assert "prefers-color-scheme: dark" in style


def test_write_html_writes_the_rendered_page(data_dir: Path, tmp_path: Path) -> None:
    protocol = read_protocol(data_dir / "pcr-protocol.json")
    path = write_html(protocol, tmp_path / "protocol.html")
    assert path.read_text(encoding="utf-8") == render_html(protocol)


def test_a_plate_says_its_kind_count_and_drops_the_legend_past_ten() -> None:
    """A page with a 96-entry legend under it is unreadable, so the count stands instead."""
    seating = {f"A{n + 1}": f"mark {n}" for n in range(11)}
    page = parse(render_html(Protocol("x", plates=(Plate("index", 96, seating=seating),))))

    figure = page.find_all("figure", cls="plate")[0]

    assert not figure.find_all("ul", cls="plate-legend")
    assert "11 kinds" in figure.find_all("figcaption")[0].text


def _drawn(caption: str) -> Protocol:
    figure = Figure(("pUC19.dna",), caption, span=(400, 700), highlight=("MCS",))
    return Protocol("Clone", steps=(Step("Cut the vector", figures=(figure,)),))


def test_a_steps_figure_draws_the_record_it_names(data_dir: Path) -> None:
    """A figure is a spec the page draws, so the map follows the record it names."""
    page = parse(render_html(_drawn("pUC19, cut at the MCS"), base=data_dir))

    [figure] = page.find_all("figure", cls="map")

    assert figure.find_all("svg")
    assert "drawing" in figure.attrs["class"]
    assert figure.find_all("figcaption")[0].text == "pUC19, cut at the MCS"


def test_a_figure_reads_its_record_from_beside_the_page(data_dir: Path, tmp_path: Path) -> None:
    """A pipeline writes the records, the data and the page into one directory."""
    (tmp_path / "pUC19.dna").write_bytes((data_dir / "pUC19.dna").read_bytes())
    path = write_html(_drawn("pUC19"), tmp_path / "protocol.html")
    assert "<svg" in path.read_text(encoding="utf-8")


def test_a_figure_whose_record_is_not_there_names_the_step_and_the_path(tmp_path: Path) -> None:
    one = Protocol("Clone", steps=(Step("Cut", figures=(Figure(("gone.dna",), "The vector"),)),))
    with pytest.raises(FileNotFoundError, match=r"step 1 'Cut'.*gone\.dna"):
        render_html(one, base=tmp_path)


def test_several_records_stack_as_rows_each_labelled_by_its_own_name(data_dir: Path) -> None:
    """A figure of more than one record is the rows the iGGA assembly figure stacks."""
    figure = Figure(("pUC19.dna", "GFP.dna"), "The vector, then the insert", linear=True)
    one = Protocol("Clone", steps=(Step("Join them", figures=(figure,)),))

    page = parse(render_html(one, base=data_dir))

    [drawn] = page.find_all("figure", cls="rows")
    rows = drawn.find_all("div", cls="row")
    assert len(rows) == 2
    assert [row.find_all("p", cls="row-name")[0].text for row in rows] == ["pUC19", "GFP"]
    assert all(row.find_all("svg") for row in rows)


def test_one_record_draws_as_it_did_with_no_row_around_it(data_dir: Path) -> None:
    """A one-record figure is the common case and keeps the markup it had."""
    page = parse(render_html(_drawn("pUC19, cut at the MCS"), base=data_dir))

    [drawn] = page.find_all("figure", cls="map")

    assert not drawn.find_all("div", cls="row")
    assert "rows" not in drawn.attrs["class"]


def test_a_highlight_lights_the_rows_answering_to_it_and_dims_the_rest(data_dir: Path) -> None:
    """Which row is lit is the whole reason the assembly figure takes a round."""
    figure = Figure(("pUC19.dna", "GFP.dna"), "The insert lit", linear=True, highlight=("GFP",))
    one = Protocol("Clone", steps=(Step("Join them", figures=(figure,)),))

    rows = parse(render_html(one, base=data_dir)).find_all("div", cls="row")

    unlit, named = (_paints(row) for row in rows)
    assert set(unlit.values()) == {layers.DIM}
    assert named["GFP"] != layers.DIM


def _paints(row: Node) -> dict[str, str]:
    """What each feature of a row is drawn in: `layers.DIM` wherever the highlight left it unlit."""
    return {
        group.attrs["data-name"]: next(
            n.attrs["fill"] for n in group.iter() if n.attrs.get("fill", "none") != "none"
        )
        for group in row.find_all("g", cls="feature")
    }


def test_a_highlight_no_record_answers_to_names_the_step(data_dir: Path) -> None:
    figure = Figure(("pUC19.dna", "GFP.dna"), "Nothing lit", linear=True, highlight=("mCherry",))
    one = Protocol("Clone", steps=(Step("Join them", figures=(figure,)),))
    with pytest.raises(ValueError, match=r"step 1 'Join them'.*'mCherry'"):
        render_html(one, base=data_dir)


def _sampled() -> Transfer:
    return Transfer(
        "Sample a quarter",
        tuple(
            Move(Well("picked", at), Well("index", to), 1.0)
            for at, to in (("A2", "A1"), ("A4", "A2"), ("C2", "B1"))
        ),
        instrument="multichannel pipette",
    )


def test_a_stamp_is_drawn_as_its_two_plates_and_keeps_its_moves_behind_a_toggle() -> None:
    """One pattern draws once: the move list spells the same thing out a row at a time."""
    one = Protocol(
        "Index",
        plates=(Plate("picked", 384), Plate("index", 96)),
        steps=(Step("Sample", transfers=(_sampled(),)),),
    )

    figure = parse(render_html(one)).find_all("figure", cls="transfer")[0]

    assert len(figure.find_all("svg")) == 2
    caption = figure.find_all("figcaption")[0].text
    assert "every other row and column, starting A2" in caption
    assert "3 moves" in caption
    assert "1 µL each" in caption
    toggle = figure.find_all("details")[0]
    assert "open" not in toggle.attrs
    assert "picked C2" in toggle.text


def test_a_transfer_naming_a_plate_the_protocol_does_not_declare_stays_a_table() -> None:
    """Nothing says what those wells look like, so the moves are all the page has."""
    one = Protocol("Index", steps=(Step("Sample", transfers=(_sampled(),)),))

    figure = parse(render_html(one)).find_all("figure", cls="transfer")[0]

    assert not figure.find_all("svg")
    assert not figure.find_all("details")
    assert "picked C2" in figure.find_all("table")[0].text


def _ordered(count: int, *, blocks: bool = False) -> Protocol:
    """A sheet of `count` primers for two purposes, seated in the order it lists them.

    The purposes alternate down the plate, or take half of it each where `blocks`.
    """

    def purpose(n: int) -> str:
        return ("gene" if n * 2 <= count else "index") if blocks else ("index" if n % 2 else "gene")

    oligos = tuple(
        Oligo(f"OP{n}", "ACGT" * (5 + n % 2), purpose=purpose(n), tm_c=60.0 + n)
        for n in range(1, count + 1)
    )
    stock = Plate("stock", 96)
    seating = dict(zip(stock.well_names, (oligo.name for oligo in oligos), strict=False))
    return Protocol(
        "Order",
        oligos=oligos,
        order_sheet="../primers.tsv",
        plates=(Plate("stock", 96, seating=seating),),
    )


def _summary_row(one: Protocol, purpose: str) -> list[str]:
    """Return the summary row for `purpose`, cell by cell."""
    section = parse(render_html(one)).find_all("section", cls="oligos")[0]
    summary = section.find_all("table")[0]
    row = next(r for r in summary.find_all("tr") if purpose in r.text)
    return [cell.text for cell in row.find_all("td")]


def test_a_long_order_sheet_is_summarised_and_the_rows_go_behind_a_toggle() -> None:
    """A page says what 20 rows have in common; the rows themselves are what the file is for."""
    section = parse(render_html(_ordered(20))).find_all("section", cls="oligos")[0]

    assert section.find_all("a", href="../primers.tsv")
    assert "20 oligos · 20 to 24 bases · Tm 61.0 to 80.0 °C" in section.text
    toggle = section.find_all("details", cls="listing")[0]
    assert toggle.find_all("summary")[0].text == "All 20 rows"
    assert "OP20" in toggle.text


def test_a_purpose_seated_among_others_claims_neither_a_run_of_names_nor_of_wells() -> None:
    """Ten of the twenty sit every other well, so counting ten off A1 stops short of B7."""
    assert _summary_row(_ordered(20), "index") == [
        "index",
        "some of OP1 to OP19",
        "10",
        "stock, 10 wells from A1 to B7",
    ]


def test_a_purpose_filling_its_own_block_of_the_plate_reads_as_a_range() -> None:
    """Every name and every well between the ends is the group's, so the ends are a range."""
    assert _summary_row(_ordered(20, blocks=True), "index") == [
        "index",
        "OP11 to OP20",
        "10",
        "stock A11 to B8",
    ]


def test_a_note_that_cites_a_document_anchors_it_where_a_troubleshooting_row_does() -> None:
    """A reader follows the *why* back to the document without leaving the page."""
    one = Protocol(
        "Digest",
        sources={"NEB": Source("NEB Technical Guide")},
        steps=(
            Step(
                "Set up the digest",
                notes=(
                    "This run chose 2 hours.",
                    Note("Glycerol above 5% is what stars.", citation=Citation("NEB", "§2")),
                ),
            ),
        ),
    )
    notes = parse(render_html(one)).find_all(cls="notes")[0]
    assert [item.text for item in notes.find_all("li")] == [
        "This run chose 2 hours.",
        "Glycerol above 5% is what stars.[1]",
    ]
    assert [a.attrs["href"] for a in notes.find_all("a", cls="cite")] == ["#source-neb"]


def test_a_citation_is_one_number_a_document_and_its_entry_names_where_the_page_cites_it() -> None:
    """Numbered by first citation, so the list reads in order and lands each click on its entry."""
    one = Protocol(
        "Digest",
        sources={"MAN": Source("Kit manual"), "NEB": Source("Guide"), "OLD": Source("Unread")},
        steps=(
            Step(
                "Set up the digest",
                notes=(
                    Note("Glycerol above 5% is what stars.", citation=Citation("NEB", "§2")),
                    Note("The buffer comes with the kit.", citation=Citation("MAN")),
                    Note("Heat stops the enzyme.", citation=Citation("NEB", "§5")),
                ),
            ),
        ),
    )
    page = parse(render_html(one))
    assert [(a.text, a.attrs["href"]) for a in page.find_all("a", cls="cite")] == [
        ("[1]", "#source-neb"),
        ("[2]", "#source-man"),
        ("[1]", "#source-neb"),
    ]
    [listed] = page.find_all("section", cls="sources")[0].find_all("ol")
    assert [
        (item.attrs["id"], [at.text for at in item.find_all(cls="cited-at")])
        for item in listed.find_all("li")
    ] == [("source-neb", ["cited at §2 · §5"]), ("source-man", []), ("source-old", [])]


def test_a_caution_that_cites_a_document_anchors_it_where_a_note_does() -> None:
    """A hazard the bench is warned of says whose warning it is, and links the document."""
    one = Protocol(
        "Recombine",
        sources={"MAN": Source("Gateway Technology with Clonase II")},
        steps=(
            Step(
                "Set up the BP reaction",
                cautions=(
                    "Keep the enzyme mix on ice.",
                    Caution("Excess DNA inhibits it.", citation=Citation("MAN", "p. 21")),
                ),
            ),
        ),
    )
    page = parse(render_html(one))
    said = [one.text for one in page.find_all(cls="caution")]
    assert said == [
        "Caution: Keep the enzyme mix on ice.",
        "Caution: Excess DNA inhibits it.[1]",
    ]
    assert [a.attrs["href"] for a in page.find_all("a", cls="cite")] == ["#source-man"]


def _naming_files() -> Protocol:
    """A protocol whose text says the sheets its run writes, as the bench says them."""
    return Protocol(
        "Order",
        summary="Order pool.tsv, then amplify it.",
        overview={"Ordered from": "pool.tsv"},
        highlights=("pool-primers.tsv pairs a primer to a block.",),
        files=("../pool.tsv", "../pool-primers.tsv", "../changes.tsv"),
        materials=(
            Material(
                "Oligo pool",
                note="ordered from pool.tsv",
                cautions=("Thaw the tubes pool.tsv names on ice.",),
            ),
        ),
        steps=(
            Step(
                "Order the pool",
                instructions=("Order every row of pool.tsv.",),
                notes=("pool.tsv names each oligo's block.",),
                expected=("One pool, as pool.tsv has it.",),
                troubleshooting=(Troubleshooting("A short row", "Order pool.tsv again."),),
            ),
        ),
    )


def test_a_file_the_run_writes_is_linked_wherever_the_page_says_its_name() -> None:
    """A filename a gloved reader cannot open is not a filename, so every mention is a link."""
    html = render_html(_naming_files())
    page = parse(html)

    assert {a.attrs["href"] for a in page.find_all("a") if a.text.endswith(".tsv")} == {
        "../pool.tsv",
        "../pool-primers.tsv",
    }
    # No field says a file name as bare text, whichever field the pipeline put it in.
    assert ".tsv" not in re.sub(r"<[^>]+>", " ", re.sub(r"<a [^>]*>[^<]*</a>", "", html))
    # `changes.tsv` is written by the run and said by no page of it, so it links nowhere.
    assert "changes.tsv" not in page.find_all("main", cls="page")[0].text


def test_a_file_name_spelled_inside_a_longer_one_is_not_linked_on_its_own() -> None:
    """`pool.tsv` and `pool-primers.tsv` are two files, and a reader must reach the right one."""
    one = Protocol(
        "Order",
        files=("../pool.tsv", "../pool-primers.tsv"),
        steps=(Step("Amplify", instructions=("Use pool-primers.tsv on the pool.",)),),
    )

    page = parse(render_html(one))

    linked = [(a.attrs["href"], a.text) for a in page.find_all("a") if a.text.endswith(".tsv")]
    assert linked == [("../pool-primers.tsv", "pool-primers.tsv")]


def test_text_beside_a_linked_file_name_is_still_escaped() -> None:
    one = Protocol(
        "Order",
        files=("../pool.tsv",),
        steps=(Step("Order", instructions=("Order <b>pool.tsv</b> & nothing else.",)),),
    )

    html = render_html(one)

    assert 'Order &lt;b&gt;<a href="../pool.tsv">pool.tsv</a>&lt;/b&gt; &amp; nothing' in html


def test_a_short_order_sheet_stays_the_sheet_it_is() -> None:
    section = parse(render_html(_ordered(19))).find_all("section", cls="oligos")[0]

    assert not section.find_all("details", cls="listing")
    assert "19 oligos" not in section.text


def test_a_warned_row_is_seen_without_opening_the_summarised_sheet() -> None:
    """A verdict a reader has to open the sheet to find is a verdict they do not see."""
    one = _ordered(20)
    warned = replace(one.oligos[2], status="warn", checks=(Check("length", "warn", "17 bases"),))
    section = parse(render_html(replace(one, oligos=(*one.oligos[:2], warned, *one.oligos[3:]))))

    toggles = section.find_all("section", cls="oligos")[0].find_all("details")

    assert [t.attrs.get("class") for t in toggles] == ["listing", "oligo-checks"]


def test_a_move_through_a_well_the_declared_plate_has_not_got_stays_a_table() -> None:
    """A drawing that silently left out the wells it cannot place would claim the wrong move."""
    moved = Transfer(
        "Sample",
        tuple(
            Move(Well("picked", at), Well("index", to), 1.0)
            for at, to in (("A1", "A1"), ("E5", "B2"))
        ),
    )
    one = Protocol(
        "Index",
        plates=(Plate("picked", 12), Plate("index", 12)),
        steps=(Step("Sample", transfers=(moved,)),),
    )

    figure = parse(render_html(one)).find_all("figure", cls="transfer")[0]

    assert moved.stamp is not None
    assert not figure.find_all("svg")
    assert "picked E5" in figure.find_all("table")[0].text

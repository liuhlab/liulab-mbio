import re
from pathlib import Path

import pytest

from liulab_mbio.bench import plates
from liulab_mbio.plot import layers
from liulab_mbio.protocol import (
    OVERVIEW_CHARS,
    Check,
    Component,
    Figure,
    Oligo,
    Plate,
    Protocol,
    ReactionTable,
    Step,
    read_protocol,
    render_html,
    write_html,
)

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


def test_the_header_reads_title_summary_facts_sentences_then_badges(page: Node) -> None:
    header = page.find_all("header", cls="intro")[0]
    blocks = [n.attrs.get("class") or n.tag for n in header.children if isinstance(n, Node)]
    assert blocks[:5] == ["h1", "summary", "overview", "highlights", "checks"]


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
    steps = page.find_all("section", cls="step")
    assert [len(s.find_all("input", type="checkbox")) for s in steps] == [1 + 3, 1 + 1, 1 + 2]
    keys = [box.attrs["data-key"] for box in page.find_all("input", type="checkbox")]
    assert len(set(keys)) == len(keys)


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
    # the mix total.
    assert [c.text for c in table.find_all("td", cls="mix")] == [
        "87.5",
        "11",
        "2.2",
        "2.2",
        "2.2",
        "0.55",
        "106",
    ]
    template = next(r for r in table.find_all("tr") if "Template DNA" in r.text)
    assert "each tube" in template.text


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
        ["Denaturation", "95 °C", "15 s", "30"],
        ["Annealing", "55 °C", "15 s"],
        ["Extension", "68 °C", "1 min"],
        ["Final extension", "68 °C", "5 min", "1"],
        ["Hold", "4 °C", "∞", "1"],
    ]
    caption = program.find_all("figcaption")[0].text
    assert "105 °C" in caption
    assert "50 min 30 s" in caption


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


def test_a_timer_starts_from_its_duration(page: Node) -> None:
    timer = page.find_all("button", cls="timer")[0]
    assert timer.attrs["data-seconds"] == "1800"
    assert "30:00" in timer.text


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

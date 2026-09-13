import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from html.parser import HTMLParser
from pathlib import Path

import pytest

from liulab_mbio.protocol import (
    OVERVIEW_CHARS,
    Oligo,
    Protocol,
    Step,
    read_protocol,
    render_html,
    write_html,
)

VOID = {"area", "base", "br", "col", "embed", "hr", "img", "input", "link", "meta", "source", "wbr"}


@dataclass
class Node:
    tag: str
    attrs: dict[str, str]
    children: list["Node | str"] = field(default_factory=list)

    def iter(self) -> Iterator["Node"]:
        for child in self.children:
            if isinstance(child, Node):
                yield child
                yield from child.iter()

    def find_all(
        self, tag: str | tuple[str, ...] = (), cls: str = "", **attrs: str
    ) -> list["Node"]:
        tags = (tag,) if isinstance(tag, str) else tag
        return [
            n
            for n in self.iter()
            if (not tags or n.tag in tags)
            and (not cls or cls in n.attrs.get("class", "").split())
            and all(n.attrs.get(k.replace("_", "-")) == v for k, v in attrs.items())
        ]

    @property
    def text(self) -> str:
        return "".join(c if isinstance(c, str) else c.text for c in self.children).strip()


class _Builder(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.root = Node("#document", {})
        self.stack = [self.root]

    def handle_starttag(self, tag, attrs) -> None:
        node = Node(tag, {k: v or "" for k, v in attrs})
        self.stack[-1].children.append(node)
        if tag not in VOID:
            self.stack.append(node)

    def handle_startendtag(self, tag, attrs) -> None:
        self.stack[-1].children.append(Node(tag, {k: v or "" for k, v in attrs}))

    def handle_endtag(self, tag) -> None:
        for i in range(len(self.stack) - 1, 0, -1):
            if self.stack[i].tag == tag:
                del self.stack[i:]
                return

    def handle_data(self, data) -> None:
        self.stack[-1].children.append(data)


def parse(html: str) -> Node:
    builder = _Builder()
    builder.feed(html)
    return builder.root


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
    # Per-reaction volume x 4 reactions x 1.1 for overage; the last cell is the mix total.
    assert [c.text for c in table.find_all("td", cls="mix")] == [
        "87.45",
        "11",
        "2.2",
        "2.2",
        "2.2",
        "0.55",
        "105.6",
    ]
    template = next(r for r in table.find_all("tr") if "Template DNA" in r.text)
    assert "each tube" in template.text


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


def test_the_page_fits_a_phone_prints_and_follows_dark_mode(page: Node) -> None:
    assert page.find_all("meta", name="viewport")
    style = page.find_all("style")[0].text
    assert "@media print" in style
    assert "prefers-color-scheme: dark" in style


def test_write_html_writes_the_rendered_page(data_dir: Path, tmp_path: Path) -> None:
    protocol = read_protocol(data_dir / "pcr-protocol.json")
    path = write_html(protocol, tmp_path / "protocol.html")
    assert path.read_text(encoding="utf-8") == render_html(protocol)

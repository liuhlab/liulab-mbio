"""What the three pages a run shares hold: the index, the reagents and the references.

The frame around them — the nav bar, the two columns and the links between pages — is
`test_folder.py`'s. This file reads what fills them.
"""

from pathlib import Path

import pytest

from liulab_mbio.protocol import (
    REAGENTS_FILE,
    Bill,
    BillRow,
    Check,
    Citation,
    Folder,
    Hole,
    Incubation,
    Item,
    Material,
    Page,
    Plate,
    Project,
    Protocol,
    Reference,
    Source,
    Stage,
    Step,
    ThermocyclerProgram,
    Timer,
    Topic,
    Vessel,
    Wait,
    page_key,
    page_name,
    render_html,
    render_index,
    render_reagents,
    render_references,
    write_project_files,
)

from ..html import Node, parse

POOL = Item("oligo pool", "the synthesised pool, resuspended", spec=("10 ng/µL",))
BLOCKS = Item("block plasmids", "one plasmid per part, miniprepped")
LIBRARY = Item("library", "the pooled plasmid library")
VECTOR = Item("destination vector", "the retrofitted vector, miniprepped")
WATER = Material("Nuclease-free water", supplier="Thermo", catalog="AM9937", storage="room")
SMITH = Reference("Smith 2020. Pooled assembly.", url="https://example.org/smith")


def chain() -> Project:
    """A run of three protocols with enough on it to fill every page the run shares."""
    return Project(
        "AP-1 library",
        summary="Build the library in three sittings.",
        background=(
            Topic(
                "Why one pool",
                ("One pool orders every part at once.", "A second pool doubles the wait."),
            ),
        ),
        inputs=(VECTOR,),
        checks=(
            Check("coverage", "pass", "every combination is covered"),
            Check("fidelity", "warn", "one junction ligates at 94%"),
        ),
        bill=Bill(
            (
                BillRow("oligo pool", 1, unit="pool", key="S-1", charge="1200.00"),
                BillRow(
                    "sequencing",
                    2,
                    unit="run",
                    key="SEQ",
                    hole=Hole("H9", "what one sequencing run costs", "price"),
                ),
            ),
            currency="USD",
            total="1200.00",
        ),
        protocols=(
            Protocol(
                "Order the pool",
                produces=(POOL,),
                materials=(WATER,),
                equipment=("Thermocycler",),
                references=(SMITH,),
                sources={"NEB": Source("NEB catalogue")},
                steps=(
                    Step(
                        "Order the pool",
                        section="Day 1",
                        waits=(Wait("the pool to arrive", "10-15 working days", Citation("NEB")),),
                    ),
                ),
            ),
            Protocol(
                "Build the blocks",
                consumes=(POOL, VECTOR),
                produces=(BLOCKS,),
                materials=(
                    Material("Nuclease-free water", supplier="Thermo", catalog="AM9937"),
                    Material("BsaI-HFv2", supplier="NEB", catalog="R3733"),
                ),
                equipment=("Thermocycler", "Plate reader"),
                vessels=(Vessel("pool tube", kind="1.5 mL tube"),),
                plates=(Plate("blocks", 96, catalog="AB-1400"),),
                references=(SMITH, Reference("Jones 2019. Golden Gate at scale.")),
                sources={"NEB": Source("NEB catalogue"), "M0491": Source("NEB M0491 manual")},
                steps=(
                    Step("Set the reaction up", timers=(Timer("ligate", 3600),)),
                    Step(
                        "Plate the colonies", holes=(Hole("H3", "how long colonies grow", "lab"),)
                    ),
                ),
            ),
            Protocol(
                "Pool and sequence",
                consumes=(BLOCKS, Item("barcode map", "the table naming each barcode")),
                produces=(LIBRARY,),
                holes=(Hole("H7", "the depth one library needs", "unpublished"),),
                steps=(
                    Step(
                        "Run the thermocycler",
                        hands_on_seconds=60,
                        programs=(
                            ThermocyclerProgram(
                                (
                                    Stage(
                                        (
                                            Incubation("denature", 98, 30),
                                            Incubation("extend", 72, 60),
                                        ),
                                        cycles=2,
                                    ),
                                )
                            ),
                        ),
                    ),
                    Step("Send the plate away", waits=(Wait("the sequencing to come back"),)),
                ),
            ),
        ),
    )


def folder_of(project: Project) -> Folder:
    """The folder `write_project_files` would compute for `project`, without writing it."""
    return Folder(
        tuple(
            Page(one.title, page_name(n, one.title), len(one.steps), page_key(one))
            for n, one in enumerate(project.protocols, 1)
        )
    )


@pytest.fixture(scope="module")
def project() -> Project:
    return chain()


@pytest.fixture(scope="module")
def index(project: Project) -> Node:
    return parse(render_index(project, folder_of(project)))


def main_of(page: Node) -> Node:
    [main] = page.find_all("main", cls="page")
    return main


def test_the_index_explains_the_run_before_any_protocol_does(index: Node, project: Project) -> None:
    [topic] = index.find_all("section", cls="topic")
    assert topic.find_all("h2")[0].text == "Why one pool"
    assert [p.text for p in topic.find_all("p")] == list(project.background[0].body)


def test_the_background_renders_on_the_index_and_on_no_protocol_page(project: Project) -> None:
    page = parse(render_html(project.protocols[0], folder=folder_of(project), here="01.html"))
    assert not page.find_all("section", cls="topic")


def test_every_name_on_the_flow_chart_is_an_item_some_protocol_declared(
    index: Node, project: Project
) -> None:
    declared = {item.name for item in project.inputs}
    for one in project.protocols:
        declared |= {item.name for item in (*one.consumes, *one.produces)}
    drawn = [span.text for span in index.find_all("span", cls="flow-item")]
    assert drawn
    assert set(drawn) <= declared


def test_the_flow_chart_boxes_every_protocol_in_order_and_links_each_to_its_page(
    index: Node, project: Project
) -> None:
    boxes = index.find_all("a", cls="flow-box")
    assert [box.find_all("span", cls="flow-title")[0].text for box in boxes] == [
        one.title for one in project.protocols
    ]
    assert [box.attrs["href"] for box in boxes] == [page.href for page in folder_of(project).pages]
    assert not any("step" in box.text for box in boxes)


def test_the_flow_chart_says_where_each_handed_name_comes_from(index: Node) -> None:
    said = {
        item.find_all("span", cls="flow-item")[0].text: item.find_all("span", cls="flow-from")[
            0
        ].text
        for item in index.find_all("li", cls="flow-hand-item")
    }
    assert said["destination vector"] == "the bench already holds it"
    assert said["oligo pool"] == "from Order the pool"
    assert said["barcode map"] == "nothing in the run hands this over"
    assert said["library"] == "the run ends holding it"


def schedule(index: Node) -> tuple[list[list[str]], list[str]]:
    """The schedule's protocol rows and its total, and the wait row under each protocol."""
    [table] = index.find_all("table", cls="schedule")
    rows = [[cell.text for cell in row.find_all(("th", "td"))] for row in table.find_all("tr")]
    waiting = [row.text for row in table.find_all("tr", cls="wait-row")]
    return [row for row in rows if len(row) > 1], waiting


def test_the_schedule_holds_each_protocol_and_totals_the_seconds_they_hold(
    index: Node, project: Project
) -> None:
    rows, _ = schedule(index)
    body = rows[1:-1]
    assert [row[0] for row in body] == [one.title for one in project.protocols]
    assert [row[1] for row in body] == ["1", "2", "2"]
    assert [row[2] for row in body] == ["no sourced number", "1 h", "3 min"]
    assert sum(one.held_seconds[0] for one in project.protocols) == 3780.0
    assert rows[-1][2] == "1 h 3 min"


def test_the_schedule_counts_the_steps_holding_nothing_rather_than_timing_them(
    index: Node, project: Project
) -> None:
    rows, _ = schedule(index)
    assert [row[-1] for row in rows[1:-1]] == ["1", "1", "1"]
    assert rows[-1][-1] == "3"
    assert sum(one.held_seconds[1] for one in project.protocols) == 3


def test_the_schedule_draws_a_hole_where_nobody_stated_a_number(index: Node) -> None:
    rows, _ = schedule(index)
    assert rows[0][3:5] == ["Hands-on", "Unattended"]
    # Nothing states the hands-on share of the first two protocols, so neither it nor the
    # unattended share it would be subtracted from reads as a figure.
    assert rows[1][3] == "no sourced number"
    assert rows[1][4] == "no sourced number"
    assert "1 min" in rows[3][3]


def test_a_schedule_column_nothing_states_is_left_out_and_named_once_under_the_table() -> None:
    """A column that is a hole in every row says the same thing once per row, so it goes."""
    run = Project(
        "Untimed",
        protocols=(Protocol("Mix", steps=(Step("Pipette"),)),),
    )
    index = parse(render_index(run, folder_of(run)))
    [table] = index.find_all("table", cls="schedule")
    head = [cell.text for cell in table.find_all("tr")[0].find_all("th")]
    assert head == ["Protocol", "Steps", "Holding nothing"]
    [block] = index.find_all("section", cls="schedule")
    assert "No protocol here states held, hands-on or unattended time" in block.text


def test_the_schedule_gives_the_waiting_a_row_of_its_own_under_each_protocol(
    index: Node,
) -> None:
    _, waiting = schedule(index)
    assert waiting[0].startswith("Waiting on the pool to arrive")
    assert "10-15 working days" in waiting[0]
    assert waiting[1] == "No waiting recorded."
    assert "no sourced number" in waiting[2]


def test_the_index_carries_the_run_checks_and_the_holes_summed_over_its_protocols(
    index: Node,
) -> None:
    badges = index.find_all("li", cls="check")
    assert [badge.find_all("span", cls="check-name")[0].text for badge in badges] == [
        "coverage",
        "fidelity",
    ]
    assert [badge.find_all("span", cls="verdict")[0].text for badge in badges] == ["pass", "warn"]
    [banner] = index.find_all("p", cls="hole-count")
    assert banner.text.startswith("3 numbers in this run have no source")
    found = {
        hole.find_all("span", cls="hole-id")[0].text: hole.find_all("a")[0].attrs["href"]
        for hole in index.find_all("li", cls="hole")
    }
    assert found == {
        "H9": f"{REAGENTS_FILE}#hole-H9",
        "H3": "02-build-the-blocks.html#hole-H3",
        "H7": "03-pool-and-sequence.html#hole-H7",
    }


def test_the_index_hands_every_page_its_key_so_it_can_show_how_far_the_bench_got(
    index: Node, project: Project
) -> None:
    [listed] = index.find_all("ol", cls="chain-pages")
    items = listed.find_all("li")
    assert [item.attrs["data-page-key"] for item in items] == [
        page_key(one) for one in project.protocols
    ]
    assert [item.attrs["data-steps"] for item in items] == ["1", "2", "2"]
    assert [item.find_all("span", cls="page-progress")[0].text for item in items] == [
        "1 step",
        "2 steps",
        "2 steps",
    ]


def test_the_reagents_page_merges_a_material_two_protocols_buy_into_one_row(
    project: Project,
) -> None:
    page = parse(render_reagents(project, folder_of(project)))
    [block] = main_of(page).find_all("section", cls="materials")
    [table] = block.find_all("table")
    rows = [[cell.text for cell in row.find_all(("th", "td"))] for row in table.find_all("tr")]
    assert rows[0][0] == "Name"
    assert rows[0][-1] == "Used in"
    assert [row[0] for row in rows[1:]] == ["Nuclease-free water", "BsaI-HFv2"]
    assert rows[1][-1] == "Order the pool, Build the blocks"
    assert rows[2][-1] == "Build the blocks"


def test_the_reagents_page_lists_the_equipment_the_plasticware_and_the_run_bill(
    project: Project,
) -> None:
    main = main_of(parse(render_reagents(project, folder_of(project))))
    equipment = {
        item.find_all("strong")[0].text: item.find_all("span", cls="used-in")[0].text
        for item in main.find_all("li", cls="kit")
    }
    assert equipment["Thermocycler"] == "used in Order the pool, Build the blocks"
    assert equipment["Plate reader"] == "used in Build the blocks"
    assert equipment["pool tube"] == "used in Build the blocks"
    assert equipment["blocks"] == "used in Build the blocks"
    [bill] = main.find_all("section", cls="bill")
    assert "oligo pool" in bill.text
    assert "1200.00" in bill.text


def test_the_references_page_names_every_protocol_citing_each_document(project: Project) -> None:
    main = main_of(parse(render_references(project, folder_of(project))))
    [listed] = main.find_all("ol")
    cited = {
        item.text.split(".")[0]: item.find_all("span", cls="cited-by")[0].text
        for item in listed.find_all("li")
    }
    assert cited == {
        "Smith 2020": "cited by Order the pool, Build the blocks",
        "Jones 2019": "cited by Build the blocks",
    }
    sources = {
        item.find_all("strong")[0].text: item.find_all("span", cls="cited-by")[0].text
        for item in main.find_all("section", cls="sources")[0].find_all("li")
    }
    assert sources == {
        "NEB": "cited by Order the pool, Build the blocks",
        "M0491": "cited by Build the blocks",
    }


def test_a_step_says_what_it_waits_on_where_the_waiting_falls(project: Project) -> None:
    page = parse(render_html(project.protocols[0], folder=folder_of(project), here="01.html"))
    [wait] = page.find_all("li", cls="wait")
    assert wait.text.startswith("Waiting on the pool to arrive")
    assert "10-15 working days" in wait.text
    assert wait.find_all("a", cls="cite")[0].attrs["href"] == "#source-neb"


def test_a_step_carries_the_label_of_the_part_of_the_protocol_it_belongs_to(
    project: Project,
) -> None:
    page = parse(render_html(project.protocols[0], folder=folder_of(project), here="01.html"))
    assert [one.text for one in page.find_all("p", cls="step-section")] == ["Day 1"]


def test_every_link_a_filled_page_writes_resolves_inside_the_folder(tmp_path: Path) -> None:
    write_project_files(chain(), tmp_path)
    written = {path.name for path in tmp_path.iterdir()}
    for path in sorted(tmp_path.glob("*.html")):
        for anchor in parse(path.read_text(encoding="utf-8")).find_all("a"):
            href = anchor.attrs["href"]
            if href.startswith(("#", "http:", "https:", "mailto:")):
                continue
            assert href.split("#")[0] in written, f"{path.name} links to {href}"

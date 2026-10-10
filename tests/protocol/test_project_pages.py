"""What the three pages a run shares hold: the index, the reagents and the references.

The frame around them — the nav bar, the two columns and the links between pages — is
`test_folder.py`'s. This file reads what fills them.
"""

from dataclasses import replace
from pathlib import Path

import pytest

from mbio.protocol import (
    REAGENTS_FILE,
    REFERENCES_FILE,
    Bill,
    BillRow,
    Check,
    Citation,
    Folder,
    Hole,
    Incubation,
    Item,
    Material,
    Plate,
    Project,
    Protocol,
    Reference,
    Rule,
    Source,
    Stage,
    Step,
    ThermocyclerProgram,
    Timer,
    Topic,
    Vessel,
    Wait,
    page_key,
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
                BillRow(
                    "oligo pool",
                    1,
                    unit="pool",
                    key="S-1",
                    charge="1200.00",
                    citation=Citation("NEB", "price list"),
                ),
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
                    Material(
                        "BsaI-HFv2",
                        supplier="NEB",
                        catalog="R3733",
                        citation=Citation("M0491", "step 2"),
                        rules=(
                            Rule(
                                "forbids",
                                "PEG",
                                "PEG inhibits the ligase.",
                                citation=Citation("M0491"),
                            ),
                        ),
                    ),
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


@pytest.fixture(scope="module")
def project() -> Project:
    return chain()


@pytest.fixture(scope="module")
def index(project: Project) -> Node:
    return parse(render_index(project, Folder.of(project)))


def main_of(page: Node) -> Node:
    [main] = page.find_all("main", cls="page")
    return main


def test_the_index_explains_the_run_before_any_protocol_does(index: Node, project: Project) -> None:
    [topic] = index.find_all("section", cls="topic")
    assert topic.find_all("h2")[0].text == "Why one pool"
    assert [p.text for p in topic.find_all("p")] == list(project.background[0].body)


def test_the_index_links_a_file_the_run_writes_where_its_background_names_it() -> None:
    """The index's own text names sheets no protocol page does, and links them the same way."""
    run = Project(
        "AP-1 library",
        files=("../changes.tsv",),
        background=(Topic("What it costs", ("changes.tsv has wild type beside synthesised.",)),),
    )

    index = parse(render_index(run, Folder.of(run)))

    [link] = [a for a in index.find_all("a") if a.text.endswith(".tsv")]
    assert link.attrs["href"] == "../changes.tsv"


def test_a_protocol_page_shows_its_own_background_and_not_the_runs(project: Project) -> None:
    """The run says it once, on the index; a protocol page says only what it carries itself."""
    one = project.protocols[0]
    borrowed = parse(render_html(one, folder=Folder.of(project), here="01.html"))
    assert not borrowed.find_all("section", cls="topic")

    own = replace(one, background=(Topic("Why one tube", ("Two tubes would be mixed up.",)),))
    [topic] = parse(render_html(own, folder=Folder.of(project), here="01.html")).find_all(
        "section", cls="topic"
    )
    assert topic.find_all("h2")[0].text == "Why one tube"


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
    assert [box.attrs["href"] for box in boxes] == [page.href for page in Folder.of(project).pages]
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
    assert rows[-1][2].startswith("1 h 3 min")


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
    index = parse(render_index(run, Folder.of(run)))
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
        "handoffs",
        "sources",
    ]
    assert [badge.find_all("span", cls="verdict")[0].text for badge in badges] == [
        "pass",
        "warn",
        "fail",
        "pass",
    ]
    [banner] = index.find_all("p", cls="hole-count")
    assert banner.text.startswith("3 numbers in this run have no source")
    found = {
        hole.find_all("span", cls="hole-id")[0].text: hole.find_all("a")[0].attrs["href"]
        for hole in index.find_all("li", cls="hole")
    }
    assert found == {
        # The reagents page prints a price hole as the bill row whose money it stands in.
        "H9": f"{REAGENTS_FILE}#bill",
        "H3": "02-build-the-blocks.html#hole-H3",
        "H7": "03-pool-and-sequence.html#hole-H7",
    }


def test_a_hole_states_itself_in_sentences_rather_than_clauses_run_together(index: Node) -> None:
    """Four statements with nothing between them are re-parsed halfway through."""
    [listed] = index.find_all("section", cls="holes")
    [depth] = [one for one in listed.find_all("li", cls="hole") if one.attrs.get("id") == "hole-H7"]
    assert depth.text == (
        "H7 no sourced number — the depth one library needs. "
        "Waiting on a number nobody has published. Stands on Pool and sequence."
    )


@pytest.fixture(scope="module")
def unpriced() -> Node:
    """The index of a run whose bill leaves three rows unpriced for one reason, beside a hole."""
    fills = "a price record holding a row for this item"
    run = Project(
        "Unpriced",
        bill=Bill(
            tuple(
                BillRow(
                    f"reagent {n}",
                    1,
                    unit="tube",
                    key=f"R-{n}",
                    hole=Hole(
                        f"P{n}",
                        "nothing prices 1 tube",
                        "price",
                        where=f"reagent {n}, money",
                        filled_by=fills,
                    ),
                )
                for n in (1, 2, 3)
            )
        ),
        protocols=(
            Protocol(
                "Mix",
                holes=(Hole("H1", "how long the colonies grow", "lab", filled_by="the lab"),),
                steps=(Step("Pipette"),),
            ),
        ),
    )
    return parse(render_index(run, Folder.of(run)))


def test_holes_waiting_on_one_thing_say_it_once_and_every_one_of_them_stays_listed(
    unpriced: Node,
) -> None:
    """Three holes repeating one sentence is one fact told at three times the length."""
    index = unpriced
    [group] = index.find_all("details", cls="hole-group")
    [summary] = group.find_all("summary")
    assert "3 numbers, each waiting on the same thing" in summary.text
    assert summary.text.count("Waiting on a price record.") == 1
    assert summary.text.count("Filled by") == 1
    assert summary.find_all("span", cls="hole-id")[0].text == "P1, P2, P3"
    # None is dropped: each keeps its id, its own line and what it is missing, one click away.
    inside = group.find_all("li", cls="hole")
    assert [one.attrs["id"] for one in inside] == ["hole-P1", "hole-P2", "hole-P3"]
    assert "reagent 2, money: nothing prices 1 tube." in inside[1].text
    [banner] = index.find_all("p", cls="hole-count")
    assert banner.text.startswith("4 numbers in this run have no source")


def test_a_hole_standing_alone_is_not_put_behind_a_disclosure(unpriced: Node) -> None:
    [alone] = [
        one for one in unpriced.find_all("li", cls="hole") if one.attrs.get("id") == "hole-H1"
    ]
    assert not alone.find_all("details")
    assert "Waiting on the lab's own stock." in alone.text


def test_a_total_summed_over_part_of_the_run_says_so_beside_itself(index: Node) -> None:
    """A reader plans a week around this number, so it may not read as the whole run."""
    rows, _ = schedule(index)
    assert rows[-1][2].endswith("over 2 of 3 protocols")
    # The steps are counted on every protocol, so that total carries no such qualification.
    assert rows[-1][1] == "5"
    [block] = index.find_all("section", cls="schedule")
    assert "covers only part of the run says so beside it" in block.text


def test_a_total_over_the_whole_run_is_qualified_by_nothing() -> None:
    step = Step("Pipette", timers=(Timer("mix", 60),))
    run = Project("Timed", protocols=(Protocol("Mix", steps=(step,)),))
    index = parse(render_index(run, Folder.of(run)))
    [table] = index.find_all("table", cls="schedule")
    assert [cell.text for cell in table.find_all("tr")[-1].find_all(("th", "td"))] == [
        "Total",
        "1",
        "1 min",
        "0",
    ]
    assert "covers only part of the run" not in index.find_all("section", cls="schedule")[0].text


def test_the_index_hands_every_page_its_key_so_it_can_show_how_far_the_bench_got(
    index: Node, project: Project
) -> None:
    [listed] = index.find_all("ol", cls="chain-pages")
    items = listed.find_all("li")
    assert [item.attrs["data-page-key"] for item in items] == [
        page_key(one) for one in project.protocols
    ]
    assert [item.attrs["data-steps"] for item in items] == [
        "order-the-pool",
        "set-the-reaction-up plate-the-colonies",
        "run-the-thermocycler send-the-plate-away",
    ]
    assert [item.find_all("span", cls="page-progress")[0].text for item in items] == [
        "1 step",
        "2 steps",
        "2 steps",
    ]


JOB = "read every well back"
CALLS = Item("well calls", "one design name per well")


def choosing() -> Project:
    """A run whose middle place offers two ways of one job, and says how to pick between them."""
    return Project(
        "DMX",
        background=(Topic(JOB, ("Ligation scales; index PCR is quicker to set up.",)),),
        protocols=(
            Protocol("Pick the colonies", produces=(VECTOR,), steps=(Step("Pick"),)),
            Protocol(
                "Cargo validation: barcode ligation",
                choice=JOB,
                consumes=(VECTOR,),
                produces=(CALLS,),
                steps=(Step("Ligate"),),
            ),
            Protocol(
                "Cargo validation: index PCR",
                choice=JOB,
                consumes=(VECTOR,),
                produces=(CALLS,),
                steps=(Step("Amplify"), Step("Pool")),
            ),
            Protocol("Report", consumes=(CALLS,), steps=(Step("Write it up"),)),
        ),
    )


def test_two_ways_of_one_job_take_one_place_in_the_run_and_two_file_names(tmp_path: Path) -> None:
    written = write_project_files(choosing(), tmp_path)
    assert [path.name for path in written.protocols] == [
        "01-pick-the-colonies.html",
        "02-cargo-validation-barcode-ligation.html",
        "02-cargo-validation-index-pcr.html",
        "03-report.html",
    ]


def test_a_way_names_its_sibling_and_the_page_after_a_choice_names_neither(
    tmp_path: Path,
) -> None:
    write_project_files(choosing(), tmp_path)
    read = {path.name: parse(path.read_text(encoding="utf-8")) for path in tmp_path.glob("*.html")}
    [line] = read["02-cargo-validation-barcode-ligation.html"].find_all("p", cls="neighbours")
    assert line.text.startswith(
        "Protocol 2 of 3, and one of two ways to read every well back — do this one or "
        "Cargo validation: index PCR, not both."
    )
    [guide] = [a for a in line.find_all("a") if a.text == "How to choose"]
    assert guide.attrs["href"] == "index.html#topic-read-every-well-back"
    [after] = read["03-report.html"].find_all("p", cls="neighbours")
    assert "Comes after whichever of the two ways to read every well back you did." in after.text


def test_the_index_lists_a_choice_as_one_entry_with_the_ways_under_it() -> None:
    run = choosing()
    index = parse(render_index(run, Folder.of(run)))
    [listed] = index.find_all("ol", cls="chain-pages")
    # Three entries, because the run has three places and not four pages.
    assert [one.tag for one in listed.children if isinstance(one, Node)] == ["li", "li", "li"]
    [choice] = listed.find_all("li", cls="chain-choice")
    assert choice.find_all("span", cls="choice-job")[0].text == "Read every well back"
    assert "do one of these two" in choice.text
    [ways] = choice.find_all("ul", cls="chain-ways")
    assert [a.text for a in ways.find_all("a")] == [
        "Cargo validation: barcode ligation",
        "Cargo validation: index PCR",
    ]
    assert [one.attrs["data-steps"] for one in ways.find_all("li")] == ["ligate", "amplify pool"]


def test_the_flow_chart_draws_one_box_for_a_choice_and_hands_on_what_every_way_leaves() -> None:
    run = choosing()
    index = parse(render_index(run, Folder.of(run)))
    [box] = index.find_all("div", cls="flow-choice")
    assert box.find_all("span", cls="flow-title")[0].text == "Read every well back"
    assert [a.text for a in box.find_all("ul", cls="flow-ways")[0].find_all("a")] == [
        "Cargo validation: barcode ligation",
        "Cargo validation: index PCR",
    ]
    said = {
        item.find_all("span", cls="flow-item")[0].text: item.find_all("span", cls="flow-from")[
            0
        ].text
        for item in index.find_all("li", cls="flow-hand-item")
    }
    assert said["well calls"] == "from whichever way to read every well back you did"


def timed(title: str, held: int, hands: int | None = None, choice: str = "") -> Protocol:
    """One protocol of one step, holding `held` seconds and standing over `hands` of them."""
    return Protocol(
        title,
        choice=choice,
        steps=(Step("Do it", timers=(Timer("hold", held),), hands_on_seconds=hands),),
    )


def timed_choice() -> Project:
    """A run of three places whose middle one offers two ways, each taking its own time."""
    return Project(
        "DMX",
        protocols=(
            timed("Pick the colonies", 600, 60),
            timed("Barcode ligation", 3600, 300, JOB),
            timed("Index PCR", 7200, 120, JOB),
            timed("Report", 300),
        ),
    )


def test_the_schedule_gives_a_choice_one_row_spanning_the_ways_standing_under_it() -> None:
    """The bench does one way, so one figure across two would be a number nobody measured."""
    run = timed_choice()
    index = parse(render_index(run, Folder.of(run)))
    rows, _ = schedule(index)
    assert [row[0] for row in rows[1:-1]] == [
        "Pick the colonies",
        "Read every well back",
        "Barcode ligation",
        "Index PCR",
        "Report",
    ]
    [place] = index.find_all("tr", cls="schedule-choice")
    assert [cell.text for cell in place.find_all("td")][1:3] == ["1", "1 h to 2 h"]
    assert [row[2] for row in rows[3:5]] == ["1 h", "2 h"]
    [block] = index.find_all("section", cls="schedule")
    assert "the row above them gives each column's range" in block.text


def test_a_place_says_how_many_of_its_ways_stated_a_column_it_spans() -> None:
    """One way timed and the other not is a span over one way, never the place's own figure."""
    run = Project("DMX", protocols=(timed("A", 3600, 300, JOB), timed("B", 3600, choice=JOB)))
    index = parse(render_index(run, Folder.of(run)))
    [place] = index.find_all("tr", cls="schedule-choice")
    assert [cell.text for cell in place.find_all("td")][3] == "5 min over 1 of 2 ways"


def test_a_total_no_place_covers_in_full_keeps_the_zero_beside_what_was_measured() -> None:
    """The figure is a way's own measurement, so neither the hole mark nor a bare number fits."""
    run = Project("DMX", protocols=(timed("A", 3600, 300, JOB), timed("B", 3600, choice=JOB)))
    rows, _ = schedule(parse(render_index(run, Folder.of(run))))
    assert rows[-1][3] == "5 min over 0 of 1 protocol"


def test_the_schedule_total_spans_the_ways_and_counts_the_places_of_the_run() -> None:
    """The page says protocol 2 of 3 everywhere else, so a total over part of it says 3 too."""
    run = timed_choice()
    index = parse(render_index(run, Folder.of(run)))
    rows, _ = schedule(index)
    assert rows[-1][1:3] == ["3", "1 h 15 min to 2 h 15 min"]
    # Nothing states what the last protocol takes by hand, and that is one place of three.
    assert rows[-1][3] == "3 min to 6 min over 2 of 3 protocols"


def test_a_choice_whose_ways_take_the_same_time_prints_that_time_once() -> None:
    """A span printing one number twice reads as two, and the bench reads this to plan a week."""
    run = Project("DMX", protocols=(timed("A", 3600, choice=JOB), timed("B", 3600, choice=JOB)))
    index = parse(render_index(run, Folder.of(run)))
    rows, _ = schedule(index)
    [place] = index.find_all("tr", cls="schedule-choice")
    assert [cell.text for cell in place.find_all("td")][2] == "1 h"
    assert rows[-1][2] == "1 h"


def test_the_index_says_which_name_the_chain_hands_nobody(index: Node) -> None:
    """A verdict the audit computes reaches no reader unless the index draws its detail too."""
    [said] = [one for one in index.find_all("p", cls="check-detail") if "handoffs" in one.text]
    assert "consumes 'barcode map', which nothing hands it" in said.text


def test_the_reagents_page_merges_a_material_two_protocols_buy_into_one_row(
    project: Project,
) -> None:
    page = parse(render_reagents(project, Folder.of(project)))
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
    main = main_of(parse(render_reagents(project, Folder.of(project))))
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


def test_the_reagents_page_states_a_caution_two_protocols_both_bring_once() -> None:
    chilled = ("Keep the polymerase on ice.",)
    run = Project(
        "Two sittings",
        protocols=(
            Protocol("Amplify", materials=(Material("Q5", catalog="M0491", cautions=chilled),)),
            Protocol("Index", materials=(Material("Taq", catalog="M0267", cautions=chilled),)),
        ),
    )
    page = parse(render_reagents(run, Folder.of(run)))
    [block] = main_of(page).find_all("section", cls="materials")
    assert [p.text for p in block.find_all("p", cls="caution")] == [f"Caution: {chilled[0]}"]


def test_the_references_page_names_every_protocol_citing_each_document(project: Project) -> None:
    main = main_of(parse(render_references(project, Folder.of(project))))
    cited = {
        item.text.split(".")[0]: item.find_all("span", cls="cited-by")[0].text
        for item in main.find_all("section", cls="references")[0].find_all("li")
    }
    assert cited == {
        "Smith 2020": "cited by Order the pool, Build the blocks",
        "Jones 2019": "cited by Build the blocks",
    }
    sources = {
        item.attrs["id"]: item.find_all("span", cls="cited-by")[0].text
        for item in main.find_all("section", cls="sources")[0].find_all("li")
    }
    assert sources == {
        "source-neb": "cited by Order the pool, Build the blocks",
        "source-m0491": "cited by Build the blocks",
    }


def test_the_references_page_lists_the_record_the_run_bill_cites() -> None:
    """A price record prices the run's bill, which sits on a page no protocol owns."""
    run = Project(
        "Priced",
        sources={"prices": Source("prices.csv")},
        bill=Bill((BillRow("pool", 1, key="S-1", charge="9.00", citation=Citation("prices")),)),
        protocols=(Protocol("Order the pool", sources={"NEB": Source("NEB catalogue")}),),
    )
    folder = Folder.of(run)
    main = main_of(parse(render_references(run, folder)))
    sources = {
        item.attrs["id"]: item.find_all("span", cls="cited-by")[0].text
        for item in main.find_all("section", cls="sources")[0].find_all("li")
    }
    assert sources == {
        "source-prices": "cited by the bill",
        "source-neb": "cited by Order the pool",
    }
    [bill] = main_of(parse(render_reagents(run, folder))).find_all("section", cls="bill")
    assert [one.attrs["href"] for one in bill.find_all("a", cls="cite")] == [
        f"{REFERENCES_FILE}#source-prices"
    ]


def test_a_run_source_no_row_of_its_bill_cites_is_listed_with_no_citer() -> None:
    """A record pricing nothing leaves every row a hole, so the bill cites it nowhere."""
    run = Project(
        "Unpriced",
        sources={"prices": Source("prices.csv")},
        bill=Bill((BillRow("pool", 1, key="S-1", hole=Hole("H1", "what a pool costs", "price")),)),
    )
    main = main_of(parse(render_references(run, Folder.of(run))))
    [listed] = main.find_all("section", cls="sources")[0].find_all("li")
    assert listed.attrs["id"] == "source-prices"
    assert not listed.find_all("span", cls="cited-by")


def test_a_step_says_what_it_waits_on_where_the_waiting_falls(project: Project) -> None:
    page = parse(render_html(project.protocols[0], folder=Folder.of(project), here="01.html"))
    [wait] = page.find_all("li", cls="wait")
    assert wait.text.startswith("Waiting on the pool to arrive")
    assert "10-15 working days" in wait.text
    assert wait.find_all("a", cls="cite")[0].attrs["href"] == "#source-neb"


def test_a_step_carries_the_label_of_the_part_of_the_protocol_it_belongs_to(
    project: Project,
) -> None:
    page = parse(render_html(project.protocols[0], folder=Folder.of(project), here="01.html"))
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


def test_every_mark_a_filled_page_links_to_stands_on_the_page_it_names(tmp_path: Path) -> None:
    write_project_files(chain(), tmp_path)
    pages = {p.name: parse(p.read_text(encoding="utf-8")) for p in tmp_path.glob("*.html")}
    marks = {
        name: {node.attrs["id"] for node in page.iter() if "id" in node.attrs}
        for name, page in pages.items()
    }
    for name, page in pages.items():
        for anchor in page.find_all("a"):
            where, _, mark = anchor.attrs["href"].partition("#")
            if not mark or where.startswith(("http:", "https:", "mailto:")):
                continue
            found = marks.get(where or name, set())
            assert mark in found, f"{name} links to {anchor.attrs['href']}"


def test_a_citation_resolves_on_its_own_page_and_reaches_the_run_list_from_a_page_with_none(
    project: Project,
) -> None:
    folder = Folder.of(project)
    shared = parse(render_reagents(project, folder))
    assert {a.attrs["href"] for a in shared.find_all("a", cls="cite")} == {
        "references.html#source-m0491",
        "references.html#source-neb",
    }
    # A number on a page with no list of its own is the run list's, so it lands on that entry.
    run = main_of(parse(render_references(project, folder))).find_all("section", cls="sources")
    listed = [item.attrs["id"] for item in run[0].find_all("li")]
    for cited in shared.find_all("a", cls="cite"):
        assert cited.text == f"[{listed.index(cited.attrs['href'].partition('#')[2]) + 1}]"
    one = project.protocols[1]
    for page in (render_html(one), render_html(one, folder=folder, here=folder.pages[1].href)):
        assert {a.attrs["href"] for a in parse(page).find_all("a", cls="cite")} == {"#source-m0491"}

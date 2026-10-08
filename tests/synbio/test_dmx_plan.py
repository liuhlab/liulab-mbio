"""DMX's way in: the build file, the designs sheet, and the protocol folder a plan writes."""

import json
import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from liulab_synbio.cli import app
from liulab_synbio.dmx import ReadBackPlan, plan_dmx, read_build, read_designs
from liulab_synbio.dmx.method import INDEX_MARKS

SHEET = "name\tfragments\nshort\t2\nmiddle\t4\nlong\t8\n"

#: The published standalone run, whose page states every number asserted below.
DEMO = Path(__file__).parents[2] / "docs" / "examples" / "ap1-readback"

BUILD = {
    "name": "Shelf read-back",
    "designs": "designs.tsv",
    "archive": "cargo archive plate 1",
    "route": "index PCR",
    "validate_from": 4,
    "selection": "carbenicillin",
}


def written(directory: Path, **changed: object) -> Path:
    """Write a designs sheet and a build file beside it, and return the build file."""
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "designs.tsv").write_text(SHEET, encoding="utf-8")
    path = directory / "build.json"
    path.write_text(json.dumps(BUILD | changed), encoding="utf-8")
    return path


@pytest.fixture(scope="module")
def build(tmp_path_factory) -> Path:
    """The one build every plan test here reads: three designs, read from four fragments up."""
    return written(tmp_path_factory.mktemp("run"))


@pytest.fixture(scope="module")
def planned(build: Path) -> ReadBackPlan:
    """That build planned once, which every test that only reads the plan shares."""
    return plan_dmx(build)


def test_a_build_names_a_sheet_resolved_against_its_own_directory(build):
    """The sheet is two columns the lab types, and its path is relative to the build file."""
    one = read_build(build)
    assert one.designs == build.parent / "designs.tsv"
    assert [design.name for design in read_designs(one.designs)] == ["short", "middle", "long"]


@pytest.mark.parametrize(
    ("changed", "said"),
    [
        ({"archive": ""}, "archive plate"),
        ({"validate_from": -1}, "counts fragments"),
        ({"route": "ligation"}, "reads its wells back"),
        ({"route": "barcode ligation", "index_plate": "plate 1"}, "only the 'index PCR' route"),
        ({"validate_from": "4"}, r"^a build's validate_from is str, not a whole number$"),
        ({"archive": 1}, r"^a build's archive is int, not a string$"),
    ],
)
def test_a_build_is_refused_where_it_is_read(tmp_path, changed, said):
    """Each check says which field is wrong and what a build states instead."""
    with pytest.raises(ValueError, match=said):
        read_build(written(tmp_path, **changed))


def test_a_build_is_refused_by_its_keys_and_by_a_sheet_that_is_not_there(tmp_path):
    """A build is read by the keys it names, and the sheet it names is a file on disk."""
    thin = tmp_path / "thin.json"
    thin.write_text(json.dumps({"name": "thin"}), encoding="utf-8")
    with pytest.raises(
        ValueError, match=r"^a build is missing archive, designs, route, validate_from$"
    ):
        read_build(thin)
    with pytest.raises(ValueError, match=r"^a build carries unknown key\(s\) scheme$"):
        read_build(written(tmp_path, scheme="DMX"))
    said = f"a build's designs is 'nowhere.tsv', and {tmp_path / 'nowhere.tsv'} is no file"
    with pytest.raises(ValueError, match=f"^{re.escape(said)}$"):
        read_build(written(tmp_path, designs="nowhere.tsv"))


def test_a_sheet_that_is_not_two_named_columns_says_what_one_is(tmp_path):
    """A sheet is refused by its header, not read as whatever its first row happens to hold."""
    path = tmp_path / "designs.tsv"
    path.write_text("name\tsequence\nshort\tATGC\n", encoding="utf-8")
    with pytest.raises(ValueError, match="name, fragments"):
        read_designs(path)


def test_the_floor_chooses_the_designs_and_the_build_names_how_to_read_them(planned):
    """A design under the floor is never read, and the build's choices reach the read-back."""
    one = planned.validation
    assert [design.name for design in one.designs] == ["middle", "long"]
    assert len(planned.designs) == 3
    assert one.route.name == "index PCR"
    assert one.selection == "carbenicillin"


def test_a_floor_above_every_design_reads_nothing_and_is_refused(tmp_path):
    """A plan writes a page someone works through, and reading no design is not one."""
    with pytest.raises(ValueError, match="none of them back"):
        plan_dmx(written(tmp_path, validate_from=9))


def test_the_run_is_one_protocol_in_three_sections_whose_chain_resolves(planned):
    """Array and pick, mark, call: three labelled sections of one page, handed what it needs."""
    chain = planned.chain()
    assert chain.title == "Shelf read-back"
    assert len(chain.protocols) == 1
    assert [step.section for step in chain.protocols[0].steps] == [
        "Array and pick",
        "Array and pick",
        "Mark every well",
        "Mark every well",
        "Mark every well",
        "Call the wells",
    ]
    assert [check.status for check in chain.audit()] == ["pass", "pass"]
    assert chain.background


def test_every_plate_a_transfer_names_is_drawn_on_the_page(tmp_path, planned):
    """Both routes: the plate a well is moved into is one the page lists."""
    ligation = plan_dmx(written(tmp_path, route="barcode ligation", validate_from=0))
    for one in (planned.protocol(), ligation.protocol()):
        drawn = {plate.name for plate in one.plates}
        moved = {
            well.plate
            for step in one.steps
            for transfer in step.transfers
            for move in transfer.moves
            for well in (move.source, move.destination)
        }
        assert moved <= drawn
    assert ligation.validation.wells == 12


def test_naming_the_index_plate_closes_its_hole_and_leaves_it_open_otherwise(tmp_path, planned):
    """The index PCR marks are a plate the lab holds, so naming it is what fills the hole."""
    assert INDEX_MARKS in planned.protocol().all_holes
    named = plan_dmx(written(tmp_path, index_plate="index plate IDX-1")).protocol()
    assert named.all_holes == ()
    assert any(
        "index plate IDX-1" in instruction
        for step in named.steps
        for instruction in step.instructions
    )
    assert any("index plate IDX-1" in item.what for item in named.consumes)


def test_the_background_counts_the_plates_this_run_pools(planned):
    """The addressing topic follows the run, so no page claims two plates where one is pooled."""
    said = " ".join(line for topic in planned.chain().background for line in topic.body)
    assert "two plates" not in said
    assert "this run pools 1 plate, and the marks tell 96 plates apart" in said


def test_the_same_build_writes_the_same_bytes(tmp_path, build):
    """A plan is a function of its inputs, so a second run over them writes the file again."""
    first = plan_dmx(build).write(tmp_path / "one")
    second = plan_dmx(build).write(tmp_path / "two")
    for here, there in zip(first.paths, second.paths, strict=True):
        assert here.name == there.name
        assert here.read_bytes() == there.read_bytes()


def test_the_verb_writes_the_folder_and_turns_a_refusal_into_an_error_line(tmp_path, build):
    """One real run of the verb, and one build it cannot read."""
    out = tmp_path / "out"
    result = CliRunner().invoke(app, ["dmx", "plan", str(build), "--out", str(out)])
    assert result.exit_code == 0, result.output
    assert sorted(path.name for path in out.iterdir()) == [
        "01-design-read-back-index-pcr.html",
        "index.html",
        "project.json",
        "reagents.html",
        "references.html",
    ]
    refused = CliRunner().invoke(
        app, ["dmx", "plan", str(written(tmp_path / "bad", route="ligation")), "--out", str(out)]
    )
    assert refused.exit_code == 1
    assert "error:" in refused.output


def test_the_published_demo_plans_what_its_page_says():
    """The committed build and sheet still plan: 32 of 72 read back, and no hole left open."""
    made = plan_dmx(DEMO / "build.json")
    assert len(made.designs) == 72
    assert len(made.validation.designs) == 32
    assert made.validation.wells == 128
    assert [plate.wells for plate in made.validation.plates] == [384, 96, 96]
    assert made.protocol().all_holes == ()

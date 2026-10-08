"""DMX's way in: the build file, the designs sheet, and the protocol folder a plan writes."""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from liulab_synbio.cli import app
from liulab_synbio.dmx import plan_dmx, read_build, read_designs
from liulab_synbio.dmx.method import INDEX_MARKS

SHEET = "name\tfragments\nshort\t2\nmiddle\t4\nlong\t8\n"

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


def test_a_build_names_a_sheet_resolved_against_its_own_directory(tmp_path):
    """The sheet is two columns the lab types, and its path is relative to the build file."""
    one = read_build(written(tmp_path / "run"))
    assert one.designs == tmp_path / "run" / "designs.tsv"
    assert [(design.name, design.fragments) for design in read_designs(one.designs)] == [
        ("short", 2),
        ("middle", 4),
        ("long", 8),
    ]


@pytest.mark.parametrize(
    ("changed", "said"),
    [
        ({"archive": ""}, "archive plate"),
        ({"validate_from": -1}, "counts fragments"),
        ({"route": "ligation"}, "reads its wells back"),
        ({"route": "barcode ligation", "index_plate": "plate 1"}, "takes none of"),
        ({"validate_from": "4"}, "not a whole number"),
    ],
)
def test_a_build_is_refused_where_it_is_read(tmp_path, changed, said):
    """Each check says which field is wrong and what a build states instead."""
    with pytest.raises(ValueError, match=said):
        read_build(written(tmp_path, **changed))


def test_a_sheet_that_is_not_two_named_columns_says_what_one_is(tmp_path):
    """A sheet is refused by its header, not read as whatever its first row happens to hold."""
    path = tmp_path / "designs.tsv"
    path.write_text("name\tsequence\nshort\tATGC\n", encoding="utf-8")
    with pytest.raises(ValueError, match="name, fragments"):
        read_designs(path)


def test_the_floor_chooses_the_designs_and_the_bench_is_sized_from_them(tmp_path):
    """A design under the floor costs no well: the plates follow the designs actually read."""
    made = plan_dmx(written(tmp_path))
    one = made.validation
    assert [design.name for design in one.designs] == ["middle", "long"]
    assert len(made.designs) == 3
    assert one.wells == 8
    assert one.route.name == "index PCR"
    assert one.selection == "carbenicillin"


def test_a_floor_above_every_design_reads_nothing_and_is_refused(tmp_path):
    """A plan writes a page someone works through, and reading no design is not one."""
    with pytest.raises(ValueError, match="none of them back"):
        plan_dmx(written(tmp_path, validate_from=9))


def test_the_run_is_one_protocol_in_three_stages_whose_chain_resolves(tmp_path):
    """Array and pick, mark, call: three labelled stages of one page, handed what it consumes."""
    chain = plan_dmx(written(tmp_path)).chain()
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


def test_naming_the_index_plate_closes_its_hole_and_leaves_it_open_otherwise(tmp_path):
    """The index PCR marks are a plate the lab holds, so naming it is what fills the hole."""
    plain = plan_dmx(written(tmp_path / "plain")).protocol()
    assert INDEX_MARKS in plain.all_holes
    named = plan_dmx(written(tmp_path / "named", index_plate="index plate IDX-1")).protocol()
    assert named.all_holes == ()
    assert any(
        "index plate IDX-1" in instruction
        for step in named.steps
        for instruction in step.instructions
    )
    assert any(item.name == "index plate IDX-1" for item in named.consumes)


def test_the_same_build_writes_the_same_bytes(tmp_path):
    """A plan is a function of its inputs, so a second run over them writes the file again."""
    path = written(tmp_path)
    first = plan_dmx(path).write(tmp_path / "one")
    second = plan_dmx(path).write(tmp_path / "two")
    for here, there in zip(first.paths, second.paths, strict=True):
        assert here.name == there.name
        assert here.read_bytes() == there.read_bytes()


def test_the_verb_writes_the_folder_and_turns_a_refusal_into_an_error_line(tmp_path):
    """One real run of the verb, and one build it cannot read."""
    out = tmp_path / "out"
    result = CliRunner().invoke(app, ["dmx", "plan", str(written(tmp_path)), "--out", str(out)])
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

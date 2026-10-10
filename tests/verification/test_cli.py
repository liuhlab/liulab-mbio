"""`mbio sequence-verify`: its wiring, and what it prints for a planted change."""

import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from mbio.cli import app
from mbio.io import read_record
from mbio.sequence import Feature, Segment, SequenceRecord
from mbio.verification.cli import PAGE, read_result, report
from mbio.verification.judge import CONSENSUS_CAVEAT, verify
from mbio.verification.result import SequencingResult

from ..html import parse

#: The Golden Gate plan's product: its GFP insert between two tagged junctions.
PRODUCT = Path(__file__).parents[2] / "docs" / "examples" / "pUC19-GFP" / "product.dna"

#: A base inside the GFP insert, which spans 399 to 1112.
PLANTED = 700


@pytest.fixture(scope="module")
def bases() -> str:
    return read_record(PRODUCT).sequence


@pytest.fixture
def planted(bases, tmp_path) -> Path:
    """The product's consensus with one substitution in GFP."""
    swapped = "A" if bases[PLANTED] != "A" else "C"
    return _consensus(tmp_path / "planted.fasta", bases[:PLANTED] + swapped + bases[PLANTED + 1 :])


def _consensus(path: Path, bases: str) -> Path:
    path.write_text(f">consensus\n{bases}\n")
    return path


def _run(*arguments: str | Path) -> tuple[int, list[str]]:
    result = CliRunner().invoke(app, ["sequence-verify", *map(str, arguments)])
    return result.exit_code, re.sub(r"\x1b\[[0-9;]*m", "", result.output).splitlines()


def _verdicts(lines: list[str]) -> dict[str, str]:
    """Each check line's verdict by name, read off the columns the command pads."""
    columns = (re.split(r"\s{2,}", line, maxsplit=2) for line in lines)
    return {parts[0]: parts[1] for parts in columns if len(parts) == 3}


def test_a_substitution_in_the_insert_fails_it_and_both_junctions_pass(planted):
    code, lines = _run(PRODUCT, planted)
    verdicts = _verdicts(lines)
    assert code == 1
    assert [verdicts[one] for one in ("ATGA junction", "GFP insert", "TGGC junction")] == [
        "pass",
        "fail",
        "pass",
    ]
    assert f"substitution at {PLANTED + 1}" in next(one for one in lines if "GFP insert" in one)
    assert verdicts["planted.fasta"] == "pass"
    assert f"planted.fasta: {CONSENSUS_CAVEAT}" in lines
    assert lines[-1] == "not verified"


def test_the_same_consensus_unchanged_verifies(bases, tmp_path):
    code, lines = _run(PRODUCT, _consensus(tmp_path / "same.fasta", bases))
    assert code == 0
    assert set(_verdicts(lines).values()) == {"pass"}
    assert lines[-1] == "verified"


def test_a_feature_the_product_does_not_hold_is_refused():
    code, lines = _run(PRODUCT, PRODUCT, "--feature", "GFP", "--feature", "mCherry")
    assert code == 1
    assert lines == ["error: pUC19-GFP has no feature named 'mCherry'"]


def test_a_change_outside_every_region_is_listed_with_no_verdict():
    record = SequenceRecord(
        "GATTACAGGCATTCGACCTA", features=(Feature("tag", "misc_feature", (Segment(14, 18),)),)
    )
    read = SequencingResult("read.fasta", "GATTACAGGCATTCGACCTA".replace("CAT", "CCT"))
    region = Feature("insert", "misc_feature", (Segment(0, 8),))
    lines = report(verify(record, [read], [region]), [read], len(record))
    assert "substitution at 11, in no feature: outside every region, so no verdict" in lines
    assert lines[-1] == "verified"


def test_a_trace_that_is_not_this_product_prints_one_line_and_nothing_else(data_dir):
    code, lines = _run(PRODUCT, data_dir / "3730.ab1")
    assert code == 1
    [line, clone] = lines
    assert re.match(r"3730\s+fail\s+does not read as this product: \d+ of \d+ trusted bases", line)
    assert clone == "not verified"


def test_out_writes_the_clone_s_page_and_prints_its_path_last(planted, data_dir, tmp_path):
    out = tmp_path / "made"
    code, lines = _run(PRODUCT, planted, data_dir / "3730.ab1", "--out", out)
    assert code == 1
    assert lines[-2:] == ["not verified", str(out / PAGE)]
    page = parse((out / PAGE).read_text(encoding="utf-8"))
    rows = [row.find_all("th")[0].text for row in page.find_all("tbody")[0].find_all("tr")]
    assert rows[0] == "ATGA junction"
    assert "3730" in rows


@pytest.mark.parametrize("given", ["folder", "calls.tsv"])
def test_a_per_base_table_is_refused_until_it_is_read(tmp_path, given):
    path = tmp_path / given
    if given == "folder":
        path.mkdir()
    else:
        path.write_text("")
    with pytest.raises(NotImplementedError, match="per-base table is not read yet"):
        read_result(path)
    code, lines = _run(PRODUCT, path)
    assert code == 1
    assert lines == [
        f"error: {given}: the per-base table is not read yet, so give the consensus FASTA or "
        "GenBank file"
    ]

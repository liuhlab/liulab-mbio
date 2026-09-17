"""The `plot map` command: the file it writes, what it prints, and how it refuses bad input."""

import dataclasses
import re
from pathlib import Path

from typer.testing import CliRunner

from liulab_mbio.cli import app
from liulab_mbio.sequence import BindingSite, Feature, Primer, Segment, SequenceRecord, Strand
from liulab_mbio.snapgene import write_dna

from ..html import Node, parse


def _groups(page: Path) -> list[Node]:
    """What the page draws: each group a hover reads."""
    tree = parse(page.read_text(encoding="utf-8"))
    return [group for group in tree.find_all("g") if "data-kind" in group.attrs]


def _run(*arguments: str) -> tuple[int, list[str]]:
    result = CliRunner().invoke(app, ["plot", "map", *arguments])
    return result.exit_code, re.sub(r"\x1b\[[0-9;]*m", "", result.output).splitlines()


def test_the_command_writes_the_page_and_prints_the_file_written(
    puc19_file: Path, tmp_path: Path
) -> None:
    out = tmp_path / "pUC19.html"
    code, lines = _run(str(puc19_file), "-o", str(out))
    assert code == 0
    assert lines == [str(out)]
    assert "<svg" in out.read_text(encoding="utf-8")


def test_the_command_refuses_a_format_it_does_not_write(puc19_file: Path, tmp_path: Path) -> None:
    code, lines = _run(str(puc19_file), "--output", str(tmp_path / "pUC19.svg"))
    assert code == 1
    assert lines[0].startswith("error: cannot write a map as 'pUC19.svg'")


def test_the_command_refuses_a_file_it_cannot_read(tmp_path: Path) -> None:
    notes = tmp_path / "notes.txt"
    notes.write_text("ACGT")
    code, lines = _run(str(notes), "-o", str(tmp_path / "map.html"))
    assert code == 1
    assert lines[0].startswith("error: no reader for a 'txt' file")


def test_the_command_refuses_a_record_that_is_not_there(tmp_path: Path) -> None:
    code, _ = _run(str(tmp_path / "absent.dna"), "-o", str(tmp_path / "map.html"))
    assert code != 0
    assert not (tmp_path / "map.html").exists()


def test_the_command_switches_layers_and_feature_types_and_names_enzymes(
    puc19: SequenceRecord, tmp_path: Path
) -> None:
    record = tmp_path / "primed.dna"
    write_dna(
        dataclasses.replace(
            puc19,
            features=(*puc19.features, Feature("pUC19", "source", (Segment(0, len(puc19)),))),
            primers=(
                Primer(
                    "M13 fwd",
                    "GTAAAACGACGGCCAGT",
                    binding_sites=(BindingSite(378, 395, Strand.FORWARD),),
                ),
            ),
        ),
        record,
    )
    out = tmp_path / "map.html"
    code, _ = _run(
        str(record), "-o", str(out), "--enzyme", "BamHI", "--enzyme", "EcoRI", "--hide-type",
        "CDS", "--hide-type", "promoter", "--source", "--no-primers",
    )  # fmt: skip
    assert code == 0
    groups = _groups(out)
    types = {group.attrs["data-type"] for group in groups}
    assert "source" in types
    assert not types & {"CDS", "promoter", "primer"}
    assert [g.attrs["data-name"] for g in groups if g.attrs["data-kind"] == "cut_site"] == [
        "EcoRI",
        "BamHI",
    ]
    code, _ = _run(str(record), "-o", str(out), "--no-features", "--no-cut-sites")
    assert code == 0
    assert {group.attrs["data-kind"] for group in _groups(out)} == {"primer"}


def test_the_command_refuses_an_enzyme_no_shipped_enzyme_answers_to(
    puc19_file: Path, tmp_path: Path
) -> None:
    code, lines = _run(str(puc19_file), "-o", str(tmp_path / "map.html"), "--enzyme", "EcoRJ")
    assert code == 1
    assert lines[0].startswith("error: ")
    assert "EcoRJ" in lines[0]
    assert not (tmp_path / "map.html").exists()

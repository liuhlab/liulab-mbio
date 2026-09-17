"""The `plot map` command: the file it writes, what it prints, and how it refuses bad input."""

import dataclasses
import re
import struct
from pathlib import Path
from typing import Any

import pypdf
import pytest
from typer.testing import CliRunner

from liulab_mbio.cli import app
from liulab_mbio.plot import circular, draw_map
from liulab_mbio.sequence import BindingSite, Feature, Primer, Segment, SequenceRecord, Strand
from liulab_mbio.snapgene import write_dna

from ..html import Node, parse
from . import crowds


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


def test_the_command_writes_each_format_asked_for_laying_the_map_out_once(
    puc19_file: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    extent = draw_map(puc19_file).layout.extent
    lay_out, layouts = circular.layout, []

    def layout(*args: Any, **kwargs: Any) -> circular.CircularMap:
        layouts.append(args)
        return lay_out(*args, **kwargs)

    monkeypatch.setattr(circular, "layout", layout)
    outs = [tmp_path / name for name in ("pUC19.html", "pUC19.png", "pUC19.pdf")]
    code, lines = _run(str(puc19_file), *(f"--output={out}" for out in outs), "--dpi", "18")
    assert code == 0
    assert lines == [str(out) for out in outs]
    assert len(layouts) == 1
    assert "<svg" in outs[0].read_text(encoding="utf-8")
    width = struct.unpack(">I", outs[1].read_bytes()[16:20])[0]
    assert width == pytest.approx(extent.width * 18 / 72, abs=1)
    assert len(pypdf.PdfReader(outs[2]).pages) == 1


def test_the_command_prints_what_the_map_hid_after_the_files_written(tmp_path: Path) -> None:
    record, out = tmp_path / "crowd.dna", tmp_path / "crowd.html"
    write_dna(crowds.ecori_crowd(), record)
    code, lines = _run(str(record), "-o", str(out), "--enzyme", "EcoRI", "--enzyme", "HindIII")
    assert code == 0
    [written, notice] = lines
    assert written == str(out)
    assert re.fullmatch(r"\d+ enzyme sites are hidden", notice)


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


def _texts(page: Path) -> list[str]:
    return [text.text for text in parse(page.read_text(encoding="utf-8")).find_all("text")]


def test_the_command_draws_a_region_by_feature_or_one_based_span_and_opens_a_circle(
    puc19: SequenceRecord, puc19_file: Path, tmp_path: Path
) -> None:
    out = tmp_path / "map.html"
    for region, title in [
        ("MCS", "396 .. 452 (57 bp)"),
        ("396..452", "396 .. 452 (57 bp)"),
        ("396 .. 452", "396 .. 452 (57 bp)"),
        ("2680..10", "2680 .. 10 (17 bp)"),
        ("11..10", "11 .. 10 (2686 bp)"),
    ]:
        code, _ = _run(str(puc19_file), "-o", str(out), "--region", region)
        assert code == 0
        assert title in _texts(out)
    # Text of that form is a span, even where a feature is named so.
    named = tmp_path / "named.dna"
    write_dna(
        dataclasses.replace(
            puc19, features=(Feature("12..40", "misc_feature", (Segment(99, 200),)),)
        ),
        named,
    )
    code, _ = _run(str(named), "-o", str(out), "--region", "12..40")
    assert code == 0
    assert "12 .. 40 (29 bp)" in _texts(out)
    code, _ = _run(str(puc19_file), "-o", str(out), "--linear")
    assert code == 0
    assert "2686 bp" in _texts(out)
    assert len(parse(out.read_text(encoding="utf-8")).find_all("circle")) == 6


def test_the_command_refuses_a_region_off_the_record_or_naming_no_feature(
    puc19_file: Path, data_dir: Path, tmp_path: Path
) -> None:
    out = tmp_path / "map.html"
    for record, region, message in [
        (puc19_file, "0..10", "error: --region 0..10: a position runs from 1 to 2686"),
        (puc19_file, "10..2687", "error: --region 10..2687: a position runs from 1 to 2686"),
        (data_dir / "GFP.dna", "700..10", "error: --region 700..10 runs across the origin"),
        (puc19_file, "lacZ", "error: record 'pUC19' has no feature called 'lacZ'"),
    ]:
        code, lines = _run(str(record), "-o", str(out), "--region", region)
        assert code == 1
        assert lines[0].startswith(message)
        assert not out.exists()


def test_the_command_draws_the_sequence_view_in_rows_of_the_bases_asked_for_on_one_strand(
    gfp: SequenceRecord, gfp_file: Path, tmp_path: Path
) -> None:
    out = tmp_path / "map.html"
    code, lines = _run(
        str(gfp_file), "-o", str(out), "--sequence-view", "--bases-per-row", "100", "--one-strand"
    )
    assert (code, lines) == (0, [str(out)])
    [view] = parse(out.read_text(encoding="utf-8")).find_all("figure", cls="sequence-view")
    rows = view.find_all("g", cls="row")
    assert [row.find_all("g", cls="top")[0].text for row in rows] == [
        gfp.sequence[start : start + 100] for start in range(0, len(gfp), 100)
    ]
    assert not view.find_all("g", cls="bottom")
    code, lines = _run(str(gfp_file), "-o", str(out), "--sequence-view", "--bases-per-row", "0")
    assert code == 1
    assert lines[0] == "error: a row of the sequence view holds at least 1 base, not 0"

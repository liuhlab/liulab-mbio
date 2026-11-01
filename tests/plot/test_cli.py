"""The `plot map` command: what it hands molviz, the files it writes, and how it refuses."""

import dataclasses
import re
import struct
from pathlib import Path

import pytest
from typer.testing import CliRunner

from mbio.cli import app
from mbio.sequence import Feature, Segment, SequenceRecord
from mbio.snapgene import write_dna

from ..conftest import Handed


def _run(*arguments: str) -> tuple[int, list[str]]:
    result = CliRunner().invoke(app, ["plot", "map", *arguments])
    return result.exit_code, re.sub(r"\x1b\[[0-9;]*m", "", result.output).splitlines()


def test_one_run_hands_molviz_every_option_and_writes_each_file_asked_for(
    puc19: SequenceRecord, tmp_path: Path, handed: Handed
) -> None:
    # Text of a span's form is a span, even where a feature is named so.
    record = tmp_path / "named.dna"
    named = Feature("2680..10", "misc_feature", (Segment(99, 200),))
    write_dna(dataclasses.replace(puc19, features=(*puc19.features, named)), record)
    outs = [tmp_path / "maps" / name for name in ("pUC19.html", "pUC19.svg", "pUC19.png")]
    code, lines = _run(
        str(record), *(f"--output={out}" for out in outs), "--region", "2680..10",
        "--linear", "--sequence-view", "--no-features", "--no-primers", "--no-cut-sites",
        "--enzyme", "EcoRI-HF", "--enzyme", "BamHI", "--hide-type", "CDS", "--source",
        "--bases-per-row", "30", "--one-strand", "--dpi", "18",
    )  # fmt: skip
    assert (code, lines) == (0, [str(out) for out in outs])
    # A commercial name is handed over as the shipped enzyme it names, and the span 0-based,
    # ending past the length across the origin.
    [(_, options)] = handed
    assert options == {
        "kinds_off": ["features", "primers", "cuts"],
        "types_off": ["CDS"],
        "bottom_strand": False,
        "span": {"start": 2679, "end": 2696},
        "map": "line",
        "sequence_view": True,
        "enzymes": ["EcoRI", "BamHI"],
        "bases_per_row": 30,
    }
    width = re.search(r'<svg[^>]* width="([\d.]+)"', outs[1].read_text(encoding="utf-8"))
    assert width is not None
    # 18 dpi is a quarter of a pixel a point.
    png = struct.unpack(">I", outs[2].read_bytes()[16:20])[0]
    assert png == pytest.approx(float(width[1]) / 4, abs=1)


def test_the_command_prints_what_the_map_hid_after_the_files_written(tmp_path: Path) -> None:
    # Fifty EcoRI sites eight bases apart: more labels than the map has room for.
    bases = "A" * 100 + "GAATTCAA" * 50 + "A" * 2500
    record, out = tmp_path / "crowd.dna", tmp_path / "crowd.html"
    write_dna(SequenceRecord(bases, topology="circular", name="crowd"), record)
    code, lines = _run(str(record), "-o", str(out), "--enzyme", "EcoRI")
    assert code == 0
    [written, notice] = lines
    assert written == str(out)
    assert re.fullmatch(r"\d+ cut sites are hidden", notice)


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        (("-o", "map.txt"), "error: cannot write"),
        (("-o", "map.html", "--enzyme", "EcoRJ"), "error: \"no shipped enzyme is called 'EcoRJ'\""),
        (("-o", "map.html", "--region", "0..10"), "error: --region 0..10: a position runs from 1"),
        (("-o", "map.html", "--region", "10..2687"), "error: --region 10..2687: a position runs"),
        (("-o", "map.html", "--region", "lacZ"), "error: record 'pUC19' has no feature called"),
    ],
)
def test_the_command_refuses_a_format_an_enzyme_or_a_region_it_cannot_draw(
    puc19_file: Path, tmp_path: Path, arguments: tuple[str, ...], message: str
) -> None:
    flag, name, *rest = arguments
    code, lines = _run(str(puc19_file), flag, str(tmp_path / "maps" / name), *rest)
    assert code == 1
    assert lines[0].startswith(message)
    # Refused before anything is written, the directory named for the files included.
    assert not (tmp_path / "maps").exists()


def test_the_command_refuses_a_span_across_the_origin_of_a_linear_record(
    gfp_file: Path, tmp_path: Path
) -> None:
    code, lines = _run(str(gfp_file), "-o", str(tmp_path / "map.html"), "--region", "700..10")
    assert code == 1
    assert lines[0].startswith("error: --region 700..10 runs across the origin")


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

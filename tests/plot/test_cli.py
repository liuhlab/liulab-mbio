"""The `plot map` command: the file it writes, what it prints, and how it refuses bad input."""

import re
from pathlib import Path

from typer.testing import CliRunner

from liulab_mbio.cli import app


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
    code, lines = _run(str(puc19_file), "--output", str(tmp_path / "pUC19.png"))
    assert code == 1
    assert lines[0].startswith("error: cannot write a map as 'pUC19.png'")


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

"""The one verb over the pipeline, which is what the skill and the docs build call."""

from __future__ import annotations

import re
from pathlib import Path

from typer.testing import CliRunner

from liulab_mbio.cli import app
from liulab_mbio.io import read_record
from liulab_mbio.snapgene import write_dna

from .records import destination_vector, entry_clone


def test_the_cli_writes_the_outputs_into_the_directory_and_prints_their_paths(
    gfp_file: Path, tmp_path: Path
) -> None:
    entry = tmp_path / "entry.dna"
    destination = tmp_path / "destination.dna"
    write_dna(entry_clone(read_record(gfp_file).sequence), entry)
    write_dna(destination_vector(), destination)
    out = tmp_path / "run"

    result = CliRunner().invoke(
        app, ["cloning", "gateway", "plan", str(entry), str(destination), "--out", str(out)]
    )

    assert result.exit_code == 0, result.output
    summary, *paths = re.sub(r"\x1b\[[0-9;]*m", "", result.output).splitlines()
    names = ["product.dna", "protocol.json", "protocol.html"]
    assert paths == [str(out / name) for name in names]
    assert all((out / name).exists() for name in names)
    assert "attB1 at" in summary
    assert "checks pass" in summary

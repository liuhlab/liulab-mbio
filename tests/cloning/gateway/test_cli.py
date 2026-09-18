"""The one verb over the pipeline, which is what the skill and the docs build call."""

from __future__ import annotations

import re
from pathlib import Path

from typer.testing import CliRunner

from liulab_mbio.cli import app
from liulab_mbio.io import read_record
from liulab_mbio.snapgene import write_dna

from .records import attb_insert, destination_vector, donor_vector


def test_the_cli_writes_the_outputs_into_the_directory_and_prints_their_paths(
    gfp_file: Path, tmp_path: Path
) -> None:
    insert = tmp_path / "insert.dna"
    donor = tmp_path / "donor.dna"
    destination = tmp_path / "destination.dna"
    write_dna(attb_insert(read_record(gfp_file).sequence), insert)
    write_dna(donor_vector(), donor)
    write_dna(destination_vector(), destination)
    out = tmp_path / "run"

    verb = ["cloning", "gateway", "plan", str(insert), str(destination)]

    result = CliRunner().invoke(app, [*verb, "--donor", str(donor), "--out", str(out)])

    assert result.exit_code == 0, result.output
    summary, *paths = re.sub(r"\x1b\[[0-9;]*m", "", result.output).splitlines()
    names = ["entry-clone.dna", "product.dna", "protocol.json", "protocol.html"]
    assert paths == [str(out / name) for name in names]
    assert all((out / name).exists() for name in names)
    assert "BP then LR" in summary
    assert "attB1 at" in summary
    assert "checks pass" in summary

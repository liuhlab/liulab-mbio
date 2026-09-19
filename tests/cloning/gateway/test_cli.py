"""The one verb over the pipeline, which is what the skill and the docs build call."""

from __future__ import annotations

import re
from pathlib import Path

from typer.testing import CliRunner

from liulab_mbio.cli import app
from liulab_mbio.snapgene import write_dna

from .records import destination_vector, donor_vector

#: A strain the manuals name for plating a reaction, and not the default, so the page naming it
#: is the `--host` option arriving.
HOST = "OmniMAX 2 T1R chemically competent E. coli"


def test_the_cli_plans_every_option_it_was_given_and_prints_the_files_it_wrote(
    gfp_file: Path, tmp_path: Path
) -> None:
    donor = tmp_path / "donor.dna"
    destination = tmp_path / "destination.dna"
    write_dna(donor_vector(), donor)
    write_dna(destination_vector(), destination)
    out = tmp_path / "run"

    verb = ["cloning", "gateway", "plan", str(gfp_file), str(destination)]
    options = [
        "--out",
        str(out),
        "--donor",
        str(donor),
        "--amplify",
        "--fusion",
        "NONE",
        "--polymerase",
        "onetaq",
        "--host",
        HOST,
        "--name",
        "pEXP-GFP",
    ]

    result = CliRunner().invoke(app, [*verb, *options])

    assert result.exit_code == 0, result.output
    summary, *paths = re.sub(r"\x1b\[[0-9;]*m", "", result.output).splitlines()
    names = ["entry-clone.dna", "product.dna", "primers.tsv", "protocol.json", "protocol.html"]
    assert paths == [str(out / name) for name in names]
    assert all((out / name).exists() for name in names)
    assert summary.startswith("pEXP-GFP:")
    assert "BP then LR" in summary
    assert "attB1 at" in summary
    assert "checks" in summary
    sheet = (out / "primers.tsv").read_text(encoding="utf-8").splitlines()
    assert sheet[1].startswith("attB-GFP forward")
    page = (out / "protocol.json").read_text(encoding="utf-8")
    assert "OneTaq DNA Polymerase" in page
    assert HOST in page

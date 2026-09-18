"""The one verb over the pipeline, which is what the skill and the docs build call.

One test, because the spine it runs is already covered by the Golden Gate command tests and
these are the slowest tests the suite has.
"""

import re
from pathlib import Path

from typer.testing import CliRunner

from liulab_mbio.cli import app


def test_the_cli_writes_the_four_outputs_and_prints_what_it_planned(
    puc19_file: Path, gfp_file: Path, tmp_path: Path
) -> None:
    out = tmp_path / "run"
    result = CliRunner().invoke(
        app,
        ["cloning", "gibson", "plan", str(puc19_file), str(gfp_file), "--out", str(out)],
    )
    assert result.exit_code == 0, result.output
    summary, *paths = re.sub(r"\x1b\[[0-9;]*m", "", result.output).splitlines()
    names = ["product.dna", "primers.tsv", "protocol.json", "protocol.html"]
    assert paths == [str(out / name) for name in names]
    assert all((out / name).exists() for name in names)
    assert "3346 bp" in summary
    assert "overlaps 16 bp, 16 bp" in summary
    assert "NEBuilder HiFi" in summary

    # A refusal is one error line and exit 1, and nothing is written.
    refused = CliRunner().invoke(
        app,
        [
            "cloning",
            "gibson",
            "plan",
            str(puc19_file),
            str(gfp_file),
            "--out",
            str(tmp_path / "refused"),
            "--product",
            "KLD",
        ],
    )
    assert refused.exit_code == 1
    assert "no assembly product" in re.sub(r"\x1b\[[0-9;]*m", "", refused.output)
    assert not (tmp_path / "refused" / "product.dna").exists()

"""The one verb over the pipeline, which is what the skill and the docs build call.

One real run, whose single call wires every option through, and the refusal. The spine it runs
is already covered by the Golden Gate command tests, and a plan is the slowest thing the suite
does.
"""

import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from liulab_mbio.cli import app
from liulab_mbio.io import read_record
from liulab_mbio.sequence import reverse_complement

#: Where the fixture's own MCS feature sits, and where GFP is cut in two to go in as two inserts.
MCS = (395, 452)
SPLIT = 360


def plain(output: str) -> str:
    """The output with rich's styling taken off, which is how a substring is looked for."""
    return re.sub(r"\x1b\[[0-9;]*m", "", output)


@pytest.fixture
def halves(gfp_file: Path, tmp_path: Path) -> tuple[Path, Path]:
    """GFP cut in two as files, the second written on the strand that does not go in."""
    gfp = read_record(gfp_file).sequence
    written = []
    for name, bases in (("GFP5", gfp[:SPLIT]), ("GFP3", reverse_complement(gfp[SPLIT:]))):
        path = tmp_path / f"{name}.fa"
        path.write_text(f">{name}\n{bases}\n", encoding="utf-8")
        written.append(path)
    return written[0], written[1]


def test_the_cli_writes_the_four_outputs_and_wires_every_option_it_takes(
    puc19_file: Path, halves: tuple[Path, Path], tmp_path: Path
) -> None:
    out = tmp_path / "run"
    result = CliRunner().invoke(
        app,
        [
            "cloning",
            "gibson",
            "plan",
            str(puc19_file),
            *(str(path) for path in halves),
            "--out",
            str(out),
            "--site",
            f"{MCS[0]}-{MCS[1]}",
            "--orientation",
            "forward",
            "--orientation",
            "reverse",
            "--route",
            "amplify",
            "--route",
            "stitch",
            "--bridge",
            "GFP5:GFP3",
            "--product",
            "nebuilder",
            "--polymerase",
            "q5",
            "--host",
            "NEB 5-alpha",
            "--name",
            "pUC19-GFP halves",
        ],
    )
    assert result.exit_code == 0, result.output
    summary, *paths = plain(result.output).splitlines()
    names = ["product.dna", "primers.tsv", "protocol.json", "protocol.html"]
    assert paths == [str(out / name) for name in names]
    assert all((out / name).exists() for name in names)
    # The name and the site are the ones asked for, and the second insert went in turned round,
    # so the plasmid is the one the whole GFP makes.
    assert "pUC19-GFP halves: 3346 bp" in summary
    assert "3 fragments" in summary
    assert "NEBuilder HiFi" in summary
    assert "overlaps 16 bp, 20 bp, 16 bp" in summary
    product = read_record(out / "product.dna")
    assert read_record(puc19_file).sequence[: MCS[0]] in product.sequence
    # The stitched part's oligos and the bridging oligo are on the sheet, and the bridged
    # junction is drawn on the product.
    sheet = (out / "primers.tsv").read_text(encoding="utf-8")
    assert "GFP5-GFP3 bridge" in sheet
    assert "GFP3 oligo 1" in sheet
    # The polymerase is read whatever the case of its name, and it is the one the PCRs name.
    assert "Q5" in (out / "protocol.json").read_text(encoding="utf-8")


def test_an_option_naming_something_the_package_does_not_ship_is_refused(
    puc19_file: Path, gfp_file: Path, tmp_path: Path
) -> None:
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
    assert "no assembly product" in plain(refused.output)
    assert not (tmp_path / "refused" / "product.dna").exists()

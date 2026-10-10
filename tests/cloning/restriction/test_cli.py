"""The one verb over the pipeline, which is what the skill and the docs build call.

One real run wires every option through in a single call; the refusal fails before any planning.
"""

import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from mbio.cli import app
from mbio.io import read_record
from mbio.protocol import read_protocol
from mbio.snapgene import read_dna, write_dna

from .records import carrying

HOST = "NEB Stable Competent E. coli (C3040)"


@pytest.fixture(scope="module")
def source_file(puc19_file: Path, gfp_file: Path, tmp_path_factory) -> str:
    """The plasmid the insert is cut out of, as the command line takes it."""
    built = carrying(read_record(puc19_file), read_record(gfp_file), "pTrc-GFP")
    path = tmp_path_factory.mktemp("source") / "pTrc-GFP.dna"
    write_dna(built, path)
    return str(path)


def run(*arguments: str):
    """Invoke the command line with these arguments."""
    return CliRunner().invoke(app, ["cloning", "restriction", "plan", *arguments])


def plain(text: str) -> str:
    """The output with rich's styling taken off."""
    return re.sub(r"\x1b\[[0-9;]*m", "", text)


def test_the_cli_writes_the_outputs_with_every_option_it_was_handed(
    puc19_file: Path, gfp_file: Path, tmp_path: Path
) -> None:
    out = tmp_path / "run"
    result = run(
        str(puc19_file),
        str(gfp_file),
        "--out",
        str(out),
        "--enzyme",
        "EcoRI",
        "--enzyme",
        "BamHI",
        "--polymerase",
        "phusion",
        "--host",
        HOST,
        "--cleanup-kit",
        "D4003",
        "--name",
        "pUC19-GFP",
    )
    assert result.exit_code == 0, result.output
    summary, *paths = plain(result.output).splitlines()
    names = ["product.dna", "primers.tsv", "protocol.json", "product-map.html", "protocol.html"]
    assert paths == [str(out / name) for name in names]
    assert all((out / name).exists() for name in names)
    assert summary.startswith("pUC19-GFP: 3388 bp, EcoRI, BamHI")
    assert "GAATTC at 397" in summary
    assert read_dna(out / "product.dna").name == "pUC19-GFP"
    # A polymerase is read whatever the case of its name, and the strain is the one asked for.
    materials = [one.name for one in read_protocol(out / "protocol.json").materials]
    assert "Phusion DNA Polymerase and its reaction buffer" in materials
    assert HOST in materials
    # The clean-up kit the run named is the one listed; the gel kit is the method's own.
    assert "DNA Clean & Concentrator-5" in materials
    assert "Monarch Spin DNA Gel Extraction Kit" in materials
    # Nothing on the page quotes the NEB column's recovery, so its manual is not cited.
    cited = " ".join(one.text for one in read_protocol(out / "protocol.json").references)
    assert "PCR & DNA Cleanup Kit" not in cited


def test_an_enzyme_that_cuts_the_insert_twice_is_refused_before_anything_is_written(
    puc19_file: Path, source_file: str, tmp_path: Path
) -> None:
    out = tmp_path / "no"
    refused = run(
        str(puc19_file), source_file, "--out", str(out), "--enzyme", "BsaI", "--enzyme", "EcoRI"
    )
    assert refused.exit_code == 1
    assert "2 BsaI sites" in plain(refused.output)
    assert not (out / "product.dna").exists()

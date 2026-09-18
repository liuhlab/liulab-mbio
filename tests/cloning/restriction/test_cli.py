"""The one verb over the pipeline, which is what the skill and the docs build call."""

import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from liulab_mbio.cli import app
from liulab_mbio.io import read_record
from liulab_mbio.snapgene import write_dna

from .test_plan import carrying


@pytest.fixture(scope="module")
def vector_and_source(puc19_file: Path, gfp_file: Path, tmp_path_factory) -> tuple[str, str]:
    """The vector file, and the plasmid the insert is cut out of, as the command line takes them."""
    built = carrying(read_record(puc19_file), read_record(gfp_file), "pTrc-GFP")
    path = tmp_path_factory.mktemp("source") / "pTrc-GFP.dna"
    write_dna(built, path)
    return str(puc19_file), str(path)


def run(*arguments: str):
    """Invoke the command line with these arguments."""
    return CliRunner().invoke(app, ["cloning", "restriction", "plan", *arguments])


def test_the_cli_writes_the_four_outputs_and_refuses_an_enzyme_that_cuts_twice(
    vector_and_source, tmp_path: Path
) -> None:
    out = tmp_path / "run"
    result = run(*vector_and_source, "--out", str(out), "--enzyme", "EcoRI", "--enzyme", "BamHI")
    assert result.exit_code == 0, result.output
    summary, *paths = re.sub(r"\x1b\[[0-9;]*m", "", result.output).splitlines()
    names = ["product.dna", "primers.tsv", "protocol.json", "protocol.html"]
    assert paths == [str(out / name) for name in names]
    assert all((out / name).exists() for name in names)
    assert "3388 bp" in summary
    assert "GAATTC at 396" in summary

    refused = run(
        *vector_and_source, "--out", str(tmp_path / "no"), "--enzyme", "BsaI", "--enzyme", "EcoRI"
    )
    assert refused.exit_code == 1
    assert "2 BsaI sites" in re.sub(r"\x1b\[[0-9;]*m", "", refused.output)
    assert not (tmp_path / "no" / "product.dna").exists()

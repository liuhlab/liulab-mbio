"""The one verb over the pipeline, which is what the skill and the docs build call."""

from pathlib import Path

import pytest
from typer.testing import CliRunner

from liulab_mbio.cli import app

#: A second insert, written here rather than kept as a fixture: it carries no BbsI site and
#: begins on bases no other junction of the set can take.
LINKER = """\
>Linker
AACGGTTCAGGTGGATCTGGCGGTTCTGGAGGCAGCGGTTCAGGAGGTTCTGGCGGATCA
GGTGGTTCAGGAGGCTCAGGTTCTGGAGGATCTGGCGGTTCAGGAGGTTCTGGATCAGGT
TCTGGAGGCAGCGGTTCAGGAGGATCTGGT
"""


@pytest.fixture(scope="module")
def inputs(data_dir: Path) -> tuple[str, str]:
    """The vector file and the insert file, as the command line takes them."""
    return str(data_dir / "pUC19.dna"), str(data_dir / "GFP.dna")


def run(*arguments: str):
    """Invoke the command line with these arguments."""
    return CliRunner().invoke(app, ["goldengate", "plan", *arguments])


def test_the_cli_writes_the_three_outputs_into_the_directory(inputs, tmp_path: Path) -> None:
    out = tmp_path / "run"
    result = run(*inputs, "--out", str(out))
    assert result.exit_code == 0, result.output
    assert (out / "product.dna").exists()
    assert (out / "primers.tsv").exists()
    assert (out / "protocol.html").exists()
    assert "BbsI" in result.output
    assert "3347 bp" in result.output


def test_the_cli_takes_a_span_as_well_as_a_feature_name(inputs, tmp_path: Path) -> None:
    named = run(*inputs, "--out", str(tmp_path / "named"), "--site", "MCS")
    spanned = run(*inputs, "--out", str(tmp_path / "spanned"), "--site", "395-452")
    assert (named.exit_code, spanned.exit_code) == (0, 0)
    assert (tmp_path / "named" / "product.dna").read_bytes() == (
        tmp_path / "spanned" / "product.dna"
    ).read_bytes()


def test_the_cli_says_which_enzyme_it_refused(inputs, tmp_path: Path) -> None:
    result = run(*inputs, "--out", str(tmp_path / "run"), "--enzyme", "BsaI")
    assert result.exit_code == 1
    assert "BsaI" in result.output
    assert not (tmp_path / "run" / "product.dna").exists()


def test_the_cli_refuses_a_polymerase_the_package_does_not_ship(inputs, tmp_path: Path) -> None:
    result = run(*inputs, "--out", str(tmp_path / "run"), "--polymerase", "Pfu")
    assert result.exit_code == 1
    assert "this package ships Q5, Phusion, Taq, OneTaq" in result.output


def test_the_cli_reads_a_polymerase_name_in_any_case(inputs, tmp_path: Path) -> None:
    result = run(*inputs, "--out", str(tmp_path / "run"), "--polymerase", "phusion")
    assert result.exit_code == 0, result.output
    assert "Phusion DNA Polymerase" in (tmp_path / "run" / "protocol.html").read_text()


def test_the_cli_takes_more_than_one_insert_file(inputs, tmp_path: Path) -> None:
    linker = tmp_path / "linker.fasta"
    linker.write_text(LINKER, encoding="utf-8")
    out = tmp_path / "run"
    result = run(*inputs, str(linker), "--out", str(out))
    assert result.exit_code == 0, result.output
    assert "3497 bp" in result.output
    assert (out / "product.dna").exists()


def test_the_cli_refuses_a_sequence_file_that_is_not_there(inputs, tmp_path: Path) -> None:
    _, insert = inputs
    result = run(str(tmp_path / "absent.dna"), insert, "--out", str(tmp_path / "run"))
    assert result.exit_code != 0

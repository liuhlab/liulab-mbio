"""The one verb over the pipeline, which is what the skill and the docs build call."""

from pathlib import Path

from typer.testing import CliRunner

from liulab_mbio.cli import app

DATA = Path(__file__).parent / "data"
VECTOR = str(DATA / "pUC19.dna")
INSERT = str(DATA / "GFP.dna")


def run(*arguments: str):
    """Invoke the command line with these arguments."""
    return CliRunner().invoke(app, ["goldengate", "plan", *arguments])


def test_the_cli_writes_the_three_outputs_into_the_directory(tmp_path: Path) -> None:
    out = tmp_path / "run"
    result = run(VECTOR, INSERT, "--out", str(out))
    assert result.exit_code == 0, result.output
    assert (out / "product.dna").exists()
    assert (out / "primers.tsv").exists()
    assert (out / "protocol.html").exists()
    assert "BbsI" in result.output
    assert "3347 bp" in result.output


def test_the_cli_takes_a_span_as_well_as_a_feature_name(tmp_path: Path) -> None:
    named = run(VECTOR, INSERT, "--out", str(tmp_path / "named"), "--site", "MCS")
    spanned = run(VECTOR, INSERT, "--out", str(tmp_path / "spanned"), "--site", "395-452")
    assert (named.exit_code, spanned.exit_code) == (0, 0)
    assert (tmp_path / "named" / "product.dna").read_bytes() == (
        tmp_path / "spanned" / "product.dna"
    ).read_bytes()


def test_the_cli_says_which_enzyme_it_refused(tmp_path: Path) -> None:
    result = run(VECTOR, INSERT, "--out", str(tmp_path / "run"), "--enzyme", "BsaI")
    assert result.exit_code == 1
    assert "BsaI" in result.output
    assert not (tmp_path / "run" / "product.dna").exists()


def test_the_cli_refuses_a_polymerase_the_package_does_not_ship(tmp_path: Path) -> None:
    result = run(VECTOR, INSERT, "--out", str(tmp_path / "run"), "--polymerase", "Pfu")
    assert result.exit_code == 1
    assert "Q5" in result.output


def test_the_cli_refuses_a_sequence_file_that_is_not_there(tmp_path: Path) -> None:
    result = run(str(tmp_path / "absent.dna"), INSERT, "--out", str(tmp_path / "run"))
    assert result.exit_code != 0

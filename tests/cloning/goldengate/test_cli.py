"""The one verb over the pipeline, which is what the skill and the docs build call."""

import json
import re
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
def vector_and_insert(puc19_file: Path, gfp_file: Path) -> tuple[str, str]:
    """The vector file and the insert file, as the command line takes them."""
    return str(puc19_file), str(gfp_file)


def run(*arguments: str):
    """Invoke the command line with these arguments."""
    return CliRunner().invoke(app, ["cloning", "goldengate", "plan", *arguments])


def test_the_cli_writes_the_four_outputs_into_the_directory_and_prints_their_paths(
    vector_and_insert, tmp_path: Path
) -> None:
    out = tmp_path / "run"
    result = run(*vector_and_insert, "--out", str(out))
    assert result.exit_code == 0, result.output
    summary, *paths = re.sub(r"\x1b\[[0-9;]*m", "", result.output).splitlines()
    names = ["product.dna", "primers.tsv", "protocol.json", "protocol.html"]
    assert paths == [str(out / name) for name in names]
    assert all((out / name).exists() for name in names)
    assert "BbsI" in summary
    assert "3347 bp" in summary


def test_a_step_deleted_from_the_protocol_data_is_gone_from_the_page_rendered_again(
    vector_and_insert, tmp_path: Path
) -> None:
    out = tmp_path / "run"
    assert run(*vector_and_insert, "--out", str(out)).exit_code == 0
    data, page = out / "protocol.json", out / "protocol.html"
    dropped = "Digest the plasmid template with DpnI"
    assert dropped in page.read_text(encoding="utf-8")
    protocol = json.loads(data.read_text(encoding="utf-8"))
    kept = [step for step in protocol["steps"] if step["title"] != dropped]
    assert len(kept) == len(protocol["steps"]) - 1
    data.write_text(json.dumps({**protocol, "steps": kept}), encoding="utf-8")
    result = CliRunner().invoke(app, ["protocol", "render", str(data)])
    assert result.exit_code == 0, result.output
    assert re.sub(r"\x1b\[[0-9;]*m", "", result.output).splitlines() == [str(page)]
    rendered = page.read_text(encoding="utf-8")
    assert dropped not in rendered
    assert all(step["title"] in rendered for step in kept)


def test_the_cli_takes_a_span_as_well_as_a_feature_name(vector_and_insert, tmp_path: Path) -> None:
    named = run(*vector_and_insert, "--out", str(tmp_path / "named"), "--site", "MCS")
    spanned = run(*vector_and_insert, "--out", str(tmp_path / "spanned"), "--site", "395-452")
    assert (named.exit_code, spanned.exit_code) == (0, 0)
    assert (tmp_path / "named" / "product.dna").read_bytes() == (
        tmp_path / "spanned" / "product.dna"
    ).read_bytes()


def test_the_cli_says_which_enzyme_it_refused(vector_and_insert, tmp_path: Path) -> None:
    result = run(*vector_and_insert, "--out", str(tmp_path / "run"), "--enzyme", "BsaI")
    assert result.exit_code == 1
    assert "BsaI" in result.output
    assert not (tmp_path / "run" / "product.dna").exists()


def test_the_cli_refuses_a_polymerase_the_package_does_not_ship(
    vector_and_insert, tmp_path: Path
) -> None:
    result = run(*vector_and_insert, "--out", str(tmp_path / "run"), "--polymerase", "Pfu")
    assert result.exit_code == 1
    assert "this package ships Q5, Phusion, Taq, OneTaq" in result.output


def test_the_cli_refuses_a_codon_table_the_package_does_not_ship(
    vector_and_insert, tmp_path: Path
) -> None:
    result = run(*vector_and_insert, "--out", str(tmp_path / "run"), "--codon-table", "yeast")
    assert result.exit_code == 1
    assert "this package ships e-coli-k12, human, mouse" in result.output


def test_the_cli_reads_a_polymerase_name_in_any_case(vector_and_insert, tmp_path: Path) -> None:
    result = run(*vector_and_insert, "--out", str(tmp_path / "run"), "--polymerase", "phusion")
    assert result.exit_code == 0, result.output
    assert "Phusion DNA Polymerase" in (tmp_path / "run" / "protocol.html").read_text()


def test_the_cli_takes_more_than_one_insert_file(vector_and_insert, tmp_path: Path) -> None:
    linker = tmp_path / "linker.fasta"
    linker.write_text(LINKER, encoding="utf-8")
    out = tmp_path / "run"
    result = run(*vector_and_insert, str(linker), "--out", str(out))
    assert result.exit_code == 0, result.output
    assert "3497 bp" in result.output
    assert (out / "product.dna").exists()


def test_the_cli_refuses_a_sequence_file_that_is_not_there(
    vector_and_insert, tmp_path: Path
) -> None:
    _, insert = vector_and_insert
    result = run(str(tmp_path / "absent.dna"), insert, "--out", str(tmp_path / "run"))
    assert result.exit_code != 0

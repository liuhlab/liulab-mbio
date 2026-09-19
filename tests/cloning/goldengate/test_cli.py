"""The one verb over the pipeline, which is what the skill and the docs build call.

One real run exercises every option's wiring in one call, on the inputs `four` plans, so the
command line and the plan tests design on one product. The rest of the file refuses.
"""

import json
import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from liulab_mbio.cli import app
from liulab_mbio.sequence import reverse_complement
from liulab_mbio.snapgene import read_dna

from .inserts import LINKER, ORIENTATIONS

#: The strain and the product name the run asks for, neither of them the default.
HOST = "NEB Stable Competent E. coli (C3040)"
NAME = "Fusion clone"


def run(*arguments: str | Path):
    """Invoke the command line with these arguments."""
    return CliRunner().invoke(app, ["cloning", "goldengate", "plan", *map(str, arguments)])


def plain(text: str) -> str:
    """The output with rich's styling taken off, which is how a test reads it."""
    return re.sub(r"\x1b\[[0-9;]*m", "", text)


@pytest.fixture(scope="module")
def planned(puc19_file: Path, insert_files: tuple[Path, ...], tmp_path_factory):
    """One real run: the three inserts as files, the site as a span, and every other option."""
    out = tmp_path_factory.mktemp("run") / "out"
    result = run(
        puc19_file,
        *insert_files,
        "--out",
        out,
        "--site",
        "395-452",
        *(word for way in ORIENTATIONS for word in ("--orientation", way)),
        "--enzyme",
        "BbsI",
        "--polymerase",
        "phusion",
        "--host",
        HOST,
        "--name",
        NAME,
    )
    assert result.exit_code == 0, result.output
    return out, plain(result.output)


def test_one_run_writes_the_four_outputs_the_options_asked_for(planned) -> None:
    out, output = planned
    summary, *paths = output.splitlines()
    names = ["product.dna", "primers.tsv", "protocol.json", "protocol.html"]
    assert paths == [str(out / name) for name in names]
    assert all((out / name).exists() for name in names)
    # The span was read as coordinates, the enzyme was taken as given, and the product is the
    # one the plan tests read, under the name this run asked for.
    assert summary.startswith(f"{NAME}: 3645 bp, BbsI, 4 fragments")
    product = read_dna(out / "product.dna")
    assert (product.name, len(product)) == (NAME, 3645)
    # The linker went in the other way round, as the second --orientation asked.
    assert reverse_complement(LINKER) in product.sequence
    assert LINKER not in product.sequence
    # The polymerase was read whatever the case of its name, and the strain is the one named.
    page = (out / "protocol.html").read_text(encoding="utf-8")
    assert "Phusion DNA Polymerase" in page
    assert HOST in page


def test_a_step_deleted_from_the_protocol_data_is_gone_from_the_page_rendered_again(
    planned, tmp_path: Path
) -> None:
    out, _ = planned
    dropped = "Digest the plasmid template with DpnI"
    assert dropped in (out / "protocol.html").read_text(encoding="utf-8")
    protocol = json.loads((out / "protocol.json").read_text(encoding="utf-8"))
    kept = [step for step in protocol["steps"] if step["title"] != dropped]
    assert len(kept) == len(protocol["steps"]) - 1
    data = tmp_path / "protocol.json"
    data.write_text(json.dumps({**protocol, "steps": kept}), encoding="utf-8")
    result = CliRunner().invoke(app, ["protocol", "render", str(data)])
    assert result.exit_code == 0, result.output
    page = tmp_path / "protocol.html"
    assert plain(result.output).splitlines() == [str(page)]
    rendered = page.read_text(encoding="utf-8")
    assert dropped not in rendered
    assert all(step["title"] in rendered for step in kept)


def test_the_cli_says_which_enzyme_it_refused(
    puc19_file: Path, gfp_file: Path, tmp_path: Path
) -> None:
    result = run(puc19_file, gfp_file, "--out", tmp_path / "run", "--enzyme", "BsaI")
    assert result.exit_code == 1
    assert "BsaI" in result.output
    assert not (tmp_path / "run" / "product.dna").exists()


def test_the_cli_refuses_a_polymerase_the_package_does_not_ship(
    puc19_file: Path, gfp_file: Path, tmp_path: Path
) -> None:
    result = run(puc19_file, gfp_file, "--out", tmp_path / "run", "--polymerase", "Pfu")
    assert result.exit_code == 1
    assert "this package ships Q5, Phusion, Taq, OneTaq" in result.output


def test_the_cli_refuses_a_codon_table_the_package_does_not_ship(
    puc19_file: Path, gfp_file: Path, tmp_path: Path
) -> None:
    result = run(puc19_file, gfp_file, "--out", tmp_path / "run", "--codon-table", "yeast")
    assert result.exit_code == 1
    assert "this package ships e-coli-k12, human, mouse" in result.output


def test_the_cli_refuses_a_sequence_file_that_is_not_there(gfp_file: Path, tmp_path: Path) -> None:
    result = run(tmp_path / "absent.dna", gfp_file, "--out", tmp_path / "run")
    assert result.exit_code != 0

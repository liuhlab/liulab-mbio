"""The `library plan` verb: one command runs the whole thing and prints a summary and the paths."""

import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from liulab_mbio.cli import app
from liulab_mbio.library.scheme import read_scheme
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.snapgene import write_dna

EXAMPLE = Path(__file__).parents[2] / "docs" / "examples" / "protein-library" / "scheme.json"

LISTS = (
    {"N_a": "MKTAEK", "N_b": "MKTCEK"},
    {"bZIP_a": "WQAFAK", "bZIP_b": "WQAYAK"},
    {"C_a": "MKTHGK", "C_b": "MKTWGK"},
)


def plain(text: str) -> str:
    """The output with rich's styling taken off, which splits an option's first dash."""
    return re.sub(r"\x1b\[[0-9;]*m", "", text)


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    """A parts FASTA and a vector file, written once for every test here."""
    scheme = read_scheme(EXAMPLE)
    directory = tmp_path_factory.mktemp("inputs")
    fasta = directory / "parts.fasta"
    fasta.write_text(
        "".join(f">{name}\n{protein}\n" for one in LISTS for name, protein in one.items())
    )
    vector = directory / "vector.dna"
    write_dna(
        SequenceRecord(
            ("TA" * 40) + scheme.internal_stuffer(-1) + ("TA" * 40),
            topology="circular",
            name="carrier",
        ),
        vector,
    )
    return fasta, vector


def run(inputs, out: Path, *extra: str):
    """Invoke `library plan` over those inputs."""
    fasta, vector = inputs
    return CliRunner().invoke(
        app,
        [
            "library",
            "plan",
            str(fasta),
            "--scheme",
            str(EXAMPLE),
            "--vector",
            str(vector),
            "--host",
            "e-coli-k12",
            "--coverage",
            "10",
            "--out",
            str(out),
            *extra,
        ],
    )


def test_one_command_plans_the_library_and_prints_the_paths(inputs, tmp_path):
    out = tmp_path / "library"

    result = run(inputs, out)

    assert result.exit_code == 0, result.output
    lines = plain(result.output).splitlines()
    assert "8 construct(s)" in lines[0]
    assert "3 round(s)" in lines[0]
    assert "checks pass" in lines[0]
    written = [Path(line) for line in lines[1:] if line.strip()]
    assert [path.name for path in written] == [
        "parts.tsv",
        "barcodes.tsv",
        "changes.tsv",
        "round-1.dna",
        "round-2.dna",
        "product.dna",
        "protocol.json",
        "protocol.html",
    ]
    for path in written:
        assert path.is_file()
        assert path.read_bytes()


def test_a_kind_that_is_neither_protein_nor_dna_is_refused(inputs, tmp_path):
    result = run(inputs, tmp_path / "library", "--kind", "rna")

    assert result.exit_code == 1
    assert "--kind is 'protein' or 'dna'" in plain(result.output)


def test_a_name_that_says_no_position_is_refused_naming_it(inputs, tmp_path):
    _, vector = inputs
    fasta = tmp_path / "strange.fasta"
    fasta.write_text(">Q_zero\nMKTAEK\n")

    result = run((fasta, vector), tmp_path / "library")

    assert result.exit_code == 1
    assert "'Q_zero' says no position" in plain(result.output)

"""The `library plan` verb: one command runs the whole thing and prints a summary and the paths."""

import json
import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from liulab_mbio.protocol import read_protocol
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.snapgene import write_dna
from liulab_synbio.cli import app

LISTS = (
    {"N_a": "MKTAEK", "N_b": "MKTCEK"},
    {"bZIP_a": "WQAFAK", "bZIP_b": "WQAYAK"},
    {"C_a": "MKTHGK", "C_b": "MKTWGK"},
)

#: A price record as a user holds one, naming the synthesis order and nothing else.
PRICES = """key,item,bands,charge,basis,currency
synthesised blocks,gene fragments,count 1-100; length_nt 1-2000,1200.00,per order,USD
"""

#: The project every test here plans, as a user writes one.
PROJECT = {
    "name": "pool",
    "positions": ["N", "bZIP", "C"],
    "parts": "parts.fasta",
    "vector": "vector.dna",
    "host": "e-coli-k12",
    "oligo_length": 350,
    "batch_size": 96,
    "coverage": 10,
    "seed": 7,
}


def plain(text: str) -> str:
    """The output with rich's styling taken off, which splits an option's first dash."""
    return re.sub(r"\x1b\[[0-9;]*m", "", text)


def write_inputs(directory: Path, *, parts: str = "parts.fasta") -> Path:
    """Write a parts FASTA, a vector carrying no stuffer, and the project naming both."""
    (directory / parts).write_text(
        "".join(f">{name}\n{protein}\n" for one in LISTS for name, protein in one.items())
    )
    write_dna(
        SequenceRecord("TA" * 100, topology="circular", name="bare"), directory / "vector.dna"
    )
    project = directory / "project.json"
    project.write_text(json.dumps({**PROJECT, "parts": parts}), encoding="utf-8")
    return project


@pytest.fixture(scope="module")
def project(tmp_path_factory):
    """The project file, written once for every test here."""
    return write_inputs(tmp_path_factory.mktemp("inputs"))


def run(project: Path, out: Path, *extra: str):
    """Invoke `library plan` over that project."""
    return CliRunner().invoke(app, ["library", "plan", str(project), "--out", str(out), *extra])


def test_one_command_plans_the_library_and_prints_the_paths(project, tmp_path):
    out = tmp_path / "library"

    prices = tmp_path / "prices.csv"
    prices.write_text(PRICES, encoding="utf-8")

    # One run for the wiring of every option: what the FASTA holds, the span a stuffer goes at,
    # the pattern the part names are read with, and the price record the bill is costed against.
    # The rest is the project file's own.
    result = run(
        project,
        out,
        "--kind",
        "protein",
        "--site",
        "100-140",
        "--pattern",
        r"^{position}_",
        "--prices",
        str(prices),
    )

    assert result.exit_code == 0, result.output
    lines = plain(result.output).splitlines()
    assert lines[0].startswith("pool round 3:")
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
    bill = read_protocol(out / "protocol.json").bill
    assert bill is not None
    assert bill.rows[0].charge == "1200.00"
    assert bill.currency == "USD"


def test_a_kind_that_is_neither_protein_nor_dna_is_refused(project, tmp_path):
    result = run(project, tmp_path / "library", "--kind", "rna")

    assert result.exit_code == 1
    assert "--kind is 'protein' or 'dna'" in plain(result.output)


def test_a_name_that_says_no_position_is_refused_naming_it(tmp_path):
    strange = write_inputs(tmp_path, parts="strange.fasta")
    (tmp_path / "strange.fasta").write_text(">Q_zero\nMKTAEK\n")

    result = run(strange, tmp_path / "library")

    assert result.exit_code == 1
    assert "'Q_zero' says no position" in plain(result.output)


def test_a_project_the_barcode_frame_rule_refuses_fails_where_it_is_read(tmp_path):
    project = write_inputs(tmp_path)
    project.write_text(json.dumps({**PROJECT, "barcode": {"length": 12}}), encoding="utf-8")

    result = run(project, tmp_path / "library")

    assert result.exit_code == 1
    assert "barcode-frame" in plain(result.output)

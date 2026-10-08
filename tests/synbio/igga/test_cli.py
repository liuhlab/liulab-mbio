"""The `igga plan` verb: one command runs the whole thing and prints a summary and the paths."""

import io
import json
import re
import zipfile
from collections.abc import Sequence
from pathlib import Path

import pytest
from typer.testing import CliRunner

from liulab_mbio.ligase import RELATIONSHIP_NS, SPREADSHEET_NS
from liulab_mbio.protocol import read_project
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.snapgene import write_dna
from liulab_synbio.cli import app

from ...chains import whole

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
    "completeness": 0.99,
    "seed": 7,
}


def plain(text: str) -> str:
    """The output with rich's styling taken off, which splits an option's first dash."""
    return re.sub(r"\x1b\[[0-9;]*m", "", text)


def write_inputs(directory: Path, *, parts: str = "parts.fasta", working: bool = False) -> Path:
    """Write a parts FASTA, a vector carrying no stuffer, and the project naming both.

    With `working`, a second such vector is written and named as the one the library moves into.
    It carries no ccdB cassette either, so the run has to be told where to put one.
    """
    (directory / parts).write_text(
        "".join(f">{name}\n{protein}\n" for one in LISTS for name, protein in one.items())
    )
    write_dna(
        SequenceRecord("TA" * 100, topology="circular", name="bare"), directory / "vector.dna"
    )
    named = {}
    if working:
        write_dna(
            SequenceRecord("TA" * 100, topology="circular", name="stock"), directory / "stock.dna"
        )
        named = {"working_vector": "stock.dna"}
    project = directory / "project.json"
    project.write_text(json.dumps({**PROJECT, "parts": parts, **named}), encoding="utf-8")
    return project


@pytest.fixture(scope="module")
def project(tmp_path_factory):
    """The project file, written once for every test here."""
    return write_inputs(tmp_path_factory.mktemp("inputs"))


def run(project: Path, out: Path, *extra: str):
    """Invoke `igga plan` over that project."""
    return CliRunner().invoke(app, ["igga", "plan", str(project), "--out", str(out), *extra])


def one_run(out: Path):
    """Return every protocol the command wrote into `out`, read back as one."""
    return whole(read_project(out / "protocol" / "project.json"))


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
        "project.json",
        "index.html",
        "01-cargo-ordering-and-pool-preparation.html",
        "02-library-assembly-in-rounds.html",
        "03-final-cargo-ligation.html",
        "reagents.html",
        "references.html",
        "library-read-primers.tsv",
    ]
    for path in written:
        assert path.is_file()
        assert path.read_bytes()
    bill = one_run(out).bill
    assert bill is not None
    assert bill.rows[0].charge == "1200.00"
    assert bill.currency == "USD"


def test_a_working_vector_carrying_no_cassette_is_planned_from_the_site_named(tmp_path):
    """The project names the backbone; only the command line says where its cassette goes.

    A run of its own, because this project's bare vector donates no cargo and so cannot pass the
    gate the shared one passes.
    """
    out = tmp_path / "library"

    result = run(
        write_inputs(tmp_path, working=True),
        out,
        "--site",
        "100-140",
        "--working-site",
        "100-140",
    )

    assert result.exit_code == 0, result.output
    titles = [step.title for step in one_run(out).steps]
    assert any(one.startswith("Pick the working vector and confirm") for one in titles)


# A ligase matrix as a supplement serves one: several sheets, one per ligase and buffer, so a
# run that names none would answer a question about T7 with T4's numbers. Nothing of such a
# workbook is here -- the two sheet names below label a matrix this file writes.
SHEETS = ("File S6. T4 PEG", "File S7. T7 PEG")


@pytest.fixture(scope="module")
def ligase_matrix(tmp_path_factory) -> Path:
    """A matrix workbook as the user holds one, written once for every test here."""
    path = tmp_path_factory.mktemp("ligase") / "File S1_NAR.xlsx"
    path.write_bytes(_workbook(SHEETS))
    return path


def test_a_ligase_matrix_reports_each_round_on_the_sheet_it_names(project, tmp_path, ligase_matrix):
    out = tmp_path / "library"

    result = run(
        project,
        out,
        "--site",
        "100-140",
        "--ligase-matrix",
        str(ligase_matrix),
        "--ligase-sheet",
        SHEETS[1],
    )

    assert result.exit_code == 0, result.output
    checks = [one for one in one_run(out).checks if "on-target" in one.name]
    assert checks
    assert all(SHEETS[1] in one.detail for one in checks)
    assert SHEETS[0] not in checks[0].detail


def test_a_workbook_of_several_sheets_is_refused_until_one_is_named(
    project, tmp_path, ligase_matrix
):
    result = run(project, tmp_path / "library", "--ligase-matrix", str(ligase_matrix))

    assert result.exit_code == 1
    assert all(name in plain(result.output) for name in SHEETS)


def _workbook(names: Sequence[str]) -> bytes:
    """The three XML parts a workbook needs, one count matrix a sheet and nothing else.

    Each sheet names its own part through a relationship, which is how a supplement writes one.
    """
    parts = [f"xl/worksheets/sheet{number}.xml" for number, _ in enumerate(names, start=1)]
    package = "http://schemas.openxmlformats.org/package/2006/relationships"
    worksheet = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet"
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr(
            "xl/workbook.xml",
            f'<workbook xmlns="{SPREADSHEET_NS[1:-1]}" xmlns:r="{RELATIONSHIP_NS[1:-1]}"><sheets>'
            + "".join(
                f'<sheet name="{name}" sheetId="{number}" r:id="rId{number}"/>'
                for number, name in enumerate(names, start=1)
            )
            + "</sheets></workbook>",
        )
        archive.writestr(
            "xl/_rels/workbook.xml.rels",
            f'<Relationships xmlns="{package}">'
            + "".join(
                f'<Relationship Id="rId{number}" Type="{worksheet}" '
                f'Target="{part.removeprefix("xl/")}"/>'
                for number, part in enumerate(parts, start=1)
            )
            + "</Relationships>",
        )
        archive.writestr(
            "xl/sharedStrings.xml",
            f'<sst xmlns="{SPREADSHEET_NS[1:-1]}">'
            + "".join(f"<si><t>{label}</t></si>" for label in ("Overhang", "AAAA", "TTTT"))
            + "</sst>",
        )
        for number, part in enumerate(parts, start=1):
            archive.writestr(part, _sheet(number))
    return buffer.getvalue()


def _sheet(count: int) -> str:
    """One count matrix: two overhangs, each seen joining the other `count` times."""
    return (
        f'<worksheet xmlns="{SPREADSHEET_NS[1:-1]}"><sheetData>'
        '<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c>'
        '<c r="C1" t="s"><v>2</v></c></row>'
        f'<row r="2"><c r="A2" t="s"><v>1</v></c><c r="C2"><v>{count}</v></c></row>'
        f'<row r="3"><c r="A3" t="s"><v>2</v></c><c r="B3"><v>{count}</v></c></row>'
        "</sheetData></worksheet>"
    )


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

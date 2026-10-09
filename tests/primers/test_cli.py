"""The `primers design` command: the files it writes on a record and on a genome, and refusals."""

import random
import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from mbio.bench.oligos import SHEET_COLUMNS
from mbio.cli import app
from mbio.primers.cli import GENOME_COLUMNS, GENOME_FILE, SHEET_FILE

#: The genome region the designed pair amplifies, 1-based with both ends included.
REGION = "chrI:401..700"


def _run(*arguments: str) -> tuple[int, list[str]]:
    result = CliRunner().invoke(app, ["primers", "design", *arguments])
    return result.exit_code, re.sub(r"\x1b\[[0-9;]*m", "", result.output).splitlines()


@pytest.fixture(scope="module")
def genome(tmp_path_factory) -> Path:
    """One sequence of random bases, with the `.fai` index samtools would write beside it."""
    from ..fasta import write_fasta

    rng = random.Random(2686)
    path = tmp_path_factory.mktemp("genome") / "planted.fa"
    write_fasta(path, {"chrI": "".join(rng.choice("ACGT") for _ in range(1200))})
    return path


def test_the_command_writes_the_sheet_for_a_pair_on_a_record(
    data_dir: Path, tmp_path: Path
) -> None:
    out = tmp_path / "pair"
    code, lines = _run(
        str(data_dir / "primer-test.dna"),
        "--out",
        str(out),
        "--region",
        "20..340",
        "--forward-tail",
        "GGTCTCA",
        "--reverse-tail",
        "GGTCTCT",
        "--forward-name",
        "fwd1",
        "--reverse-name",
        "rev1",
        "--target-tm",
        "60",
        "--polymerase",
        "Taq",
    )
    assert code == 0, lines
    assert lines[0].startswith("primer-test: GGTCTCA")
    assert [Path(line).name for line in lines[1:]] == [SHEET_FILE]
    rows = (out / SHEET_FILE).read_text(encoding="utf-8").splitlines()
    assert rows[0] == "\t".join(SHEET_COLUMNS)
    assert [row.split("\t")[0] for row in rows[1:]] == ["fwd1", "rev1"]
    # The tails are ordered with the primer, and the region holds both annealing regions.
    assert rows[1].split("\t")[1].startswith("GGTCTCA")
    assert rows[2].split("\t")[1].startswith("GGTCTCT")


def test_the_command_writes_the_genome_report_beside_the_sheet(
    genome: Path, tmp_path: Path
) -> None:
    out = tmp_path / "genomic"
    code, lines = _run(
        str(genome),
        "--out",
        str(out),
        "--assembly",
        "planted",
        "--region",
        REGION,
        "--flank",
        "200",
    )
    assert code == 0, lines
    assert lines[0].startswith(f"planted {REGION}: ")
    assert lines[1].startswith("on planted: specific after 1 search(es), checks pass")
    assert [Path(line).name for line in lines[2:]] == [SHEET_FILE, GENOME_FILE]
    rows = (out / GENOME_FILE).read_text(encoding="utf-8").splitlines()
    assert rows[0] == "\t".join(GENOME_COLUMNS)
    # One amplicon, the intended one, and it covers the region asked for.
    [found] = rows[1:]
    location, _, made_by, _, intended = found.split("\t")
    name, _, span = location.partition(":")
    start, _, end = span.partition("..")
    assert (name, made_by, intended) == ("chrI", "pair", "yes")
    assert int(start) <= 401
    assert int(end) >= 700


@pytest.mark.parametrize(
    ("arguments", "said"),
    [
        (("--region", "chrI:1..10"), "--region is START..END on a record"),
        (("--region", "1..5000"), "--region 1..5000: a span runs from 1 to 360"),
        (("--polymerase", "Pfu"), "no polymerase called 'Pfu'"),
    ],
)
def test_the_command_refuses_what_it_cannot_read(
    data_dir: Path, tmp_path: Path, arguments: tuple[str, ...], said: str
) -> None:
    code, lines = _run(str(data_dir / "primer-test.dna"), "--out", str(tmp_path), *arguments)
    assert code == 1
    assert any(said in line for line in lines), lines


def test_the_command_refuses_a_genome_region_that_names_no_sequence(
    genome: Path, tmp_path: Path
) -> None:
    code, lines = _run(
        str(genome), "--out", str(tmp_path), "--assembly", "planted", "--region", "1..10"
    )
    assert code == 1
    assert any("--region is NAME:START..END on a genome" in line for line in lines), lines

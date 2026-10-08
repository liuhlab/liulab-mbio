"""The root command line: its own verbs, and that every sub-app mounts."""

import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from mbio import __version__
from mbio.barcodes import BarcodeRules, design_barcodes
from mbio.cli import (
    BARCODE_COLUMNS,
    BARCODE_FILE,
    CHANGE_COLUMNS,
    CHANGE_FILE,
    CHECK_COLUMNS,
    CHECK_FILE,
    CODING_FILE,
    app,
)
from mbio.io import read_record


def _plain(output: str) -> str:
    """Typer prints through rich, which styles what it prints; the words are what is asserted."""
    return re.sub(r"\x1b\[[0-9;]*m", "", output)


def test_the_version_verb_prints_the_installed_version() -> None:
    # hatch-vcs derives it from git, so the value moves. That it exists is the claim.
    assert __version__
    result = CliRunner().invoke(app, ["version"])
    assert result.exit_code == 0
    assert result.output.strip() == __version__


def test_a_bare_invocation_prints_help_rather_than_nothing() -> None:
    result = CliRunner().invoke(app, [])
    assert result.exit_code != 0
    assert result.output.strip()


def test_every_pipeline_is_mounted() -> None:
    for verb in ("cloning", "protocol", "plot", "primers"):
        result = CliRunner().invoke(app, [verb])
        # `no_args_is_help` on the sub-app, so it prints its own help and exits non-zero.
        assert result.exit_code != 0
        assert result.output.strip()


def test_the_cloning_group_lists_every_method_it_mounts() -> None:
    result = CliRunner().invoke(app, ["cloning", "--help"])
    listed = _plain(result.output)
    for method in ("gateway", "gibson", "goldengate", "restriction"):
        assert method in listed


def test_the_codon_optimize_verb_writes_a_protein_as_dna() -> None:
    result = CliRunner().invoke(
        app, ["codon-optimize", "MW*", "--kind", "protein", "--host", "e-coli-k12"]
    )
    assert result.exit_code == 0, result.output
    lines = _plain(result.output).splitlines()
    assert lines[0].startswith("coding sequence: 9 bp, 3 aa, host e-coli-k12, 0 codon change(s)")
    assert lines[-1] == "ATGTGGTAA"


def test_the_codon_optimize_verb_refuses_a_kind_that_is_neither() -> None:
    result = CliRunner().invoke(
        app, ["codon-optimize", "MW*", "--kind", "rna", "--host", "e-coli-k12"]
    )
    assert result.exit_code == 1
    assert "--kind is 'protein' or 'dna'" in result.output


def test_the_codon_optimize_verb_writes_the_sequence_and_the_codons_it_moved(
    tmp_path: Path,
) -> None:
    result = CliRunner().invoke(
        app,
        [
            "codon-optimize",
            "ATGGGTCTCGGCTAA",
            "--kind",
            "dna",
            "--host",
            "e-coli-k12",
            "--forbid",
            "BsaI",
            "--name",
            "part A",
            "--out",
            str(tmp_path / "gene"),
        ],
    )
    assert result.exit_code == 0, result.output
    lines = _plain(result.output).splitlines()
    assert lines[0].startswith("part A: 15 bp, 5 aa, host e-coli-k12, 1 codon change(s)")
    assert [Path(line).name for line in lines[1:]] == [CODING_FILE, CHANGE_FILE]
    assert read_record(tmp_path / "gene" / CODING_FILE).sequence == "ATGGGTCTGGGCTAA"
    assert (tmp_path / "gene" / CHANGE_FILE).read_text(encoding="utf-8").splitlines() == [
        "\t".join(CHANGE_COLUMNS),
        "2\tCTC\tCTG\tL\tBsaI",
    ]


def test_the_codon_optimize_verb_keeps_printing_the_sequence_when_no_out_is_given() -> None:
    # What anything piping it reads today: the bases on the last line, and no file written.
    result = CliRunner().invoke(
        app, ["codon-optimize", "ATGGGTCTCGGCTAA", "--kind", "dna", "--host", "e-coli-k12"]
    )
    assert result.exit_code == 0, result.output
    assert _plain(result.output).splitlines()[-1] == "ATGGGTCTCGGCTAA"


def test_the_barcode_design_verb_writes_the_set_and_its_checks(tmp_path: Path) -> None:
    result = CliRunner().invoke(
        app,
        [
            "barcode-design",
            "6",
            "--length",
            "11",
            "--scar",
            "AGCG",
            "--phase",
            "0",
            "--forbid",
            "BsaI",
            "--forbid",
            "BbsI",
            "--distance",
            "3",
            "--metric",
            "hamming",
            "--max-homopolymer",
            "4",
            "--gc-band",
            "0.3-0.7",
            "--seed",
            "7",
            "--out",
            str(tmp_path / "set"),
        ],
    )
    assert result.exit_code == 0, result.output
    lines = _plain(result.output).splitlines()
    assert lines[0] == (
        "6 barcode(s) of 11 bases, scar AGCG, 3 apart by hamming, free of BsaI, BbsI, "
        "seed 7, checks pass"
    )
    assert [Path(line).name for line in lines[1:]] == [BARCODE_FILE, CHECK_FILE]
    rows = (tmp_path / "set" / BARCODE_FILE).read_text(encoding="utf-8").splitlines()
    assert rows[0] == "\t".join(BARCODE_COLUMNS)
    barcodes = [row.split("\t")[1] for row in rows[1:]]
    assert barcodes == list(
        design_barcodes(
            6,
            BarcodeRules(
                11,
                scar="AGCG",
                forbidden=["BsaI", "BbsI"],
                max_homopolymer=4,
                gc_band=(0.3, 0.7),
                metric="hamming",
            ),
            seed=7,
        )
    )
    checks = (tmp_path / "set" / CHECK_FILE).read_text(encoding="utf-8").splitlines()
    assert checks[0] == "\t".join(CHECK_COLUMNS)
    judged = {row.split("\t")[0]: row.split("\t")[1] for row in checks[1:]}
    # Deletion ambiguity carries no verdict: no sourced threshold judges it.
    assert judged == {"rules": "pass", "separation": "pass", "deletion_ambiguity": ""}


def test_the_barcode_design_verb_takes_the_frame_rules_away_where_nothing_translates_it(
    tmp_path: Path,
) -> None:
    arguments = ["barcode-design", "4", "--length", "10", "--scar", "AGCG", "--out"]
    refused = CliRunner().invoke(app, [*arguments, str(tmp_path / "framed")])
    assert refused.exit_code == 1
    assert "not a whole number of codons" in refused.output
    result = CliRunner().invoke(app, [*arguments, str(tmp_path / "free"), "--untranslated"])
    assert result.exit_code == 0, result.output
    assert (tmp_path / "free" / BARCODE_FILE).exists()


@pytest.mark.parametrize(
    ("option", "value", "said"),
    [
        ("--metric", "mismatch", "--metric is 'hamming' or 'sequence-levenshtein'"),
        ("--gc-band", "0.4", "--gc-band is LOW-HIGH"),
    ],
)
def test_the_barcode_design_verb_refuses_a_dial_it_cannot_read(
    tmp_path: Path, option: str, value: str, said: str
) -> None:
    result = CliRunner().invoke(
        app, ["barcode-design", "4", "--length", "8", option, value, "--out", str(tmp_path)]
    )
    assert result.exit_code == 1
    assert said in _plain(result.output)

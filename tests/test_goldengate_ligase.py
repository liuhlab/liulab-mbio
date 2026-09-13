"""Scoring overhangs against a ligase fidelity matrix the user holds on their own disk.

The archive such a file comes from is CC BY-NC 4.0 and none of it is in this repository. Every
matrix here is four overhangs written by this file, with plausible integers in the cells.
"""

import io
import zipfile
from collections.abc import Mapping, Sequence
from pathlib import Path

import pytest
from typer.testing import CliRunner

from liulab_mbio.cli import app
from liulab_mbio.goldengate import plan_assembly
from liulab_mbio.goldengate.cli import LIGASE_MATRIX_ENV
from liulab_mbio.goldengate.design import Junction, design_overhangs, fidelity, ligation_matrix
from liulab_mbio.goldengate.ligase import SPREADSHEET_NS, read_profile

DATA = Path(__file__).parent / "data"

#: Two Watson-Crick pairs seen often, and one cross pair seen rarely. A row pairs with the
#: column spelling its reverse complement, so AAAA pairs with TTTT and GGAA with TTCC.
LABELS = ("AAAA", "TTTT", "GGAA", "TTCC")
COUNTS: Mapping[str, Mapping[str, int]] = {
    "AAAA": {"TTTT": 600, "GGAA": 3},
    "TTTT": {"AAAA": 600},
    "GGAA": {"AAAA": 3, "TTCC": 500},
    "TTCC": {"GGAA": 500},
}
OBSERVATIONS = 600 + 3 + 600 + 3 + 500 + 500

#: How the archive names a file: the ligase, the incubation and the temperature.
NAMED = "FileS03_T4_18h_25C"


def comma_separated(
    labels: Sequence[str] = LABELS,
    counts: Mapping[str, Mapping[str, int]] = COUNTS,
    header: str = "Overhang",
) -> str:
    """The matrix in the comma-separated shape, zeros written out as the archive writes them."""
    lines = [",".join((header, *labels))]
    lines += [
        ",".join((label, *(str(row.get(one, 0)) for one in labels)))
        for label, row in counts.items()
    ]
    return "\n".join(lines) + "\n"


def workbook(
    labels: Sequence[str] = LABELS,
    counts: Mapping[str, Mapping[str, int]] = COUNTS,
    header: str = "Overhang",
) -> bytes:
    """The same matrix as an `.xlsx`: the three XML parts a workbook needs and nothing else."""
    shared = (header, *labels, *(one for one in counts if one not in labels))
    number_of = {label: number for number, label in enumerate(shared)}
    cells = "".join(
        f'<c r="{_reference(column, 1)}" t="s"><v>{number_of[label]}</v></c>'
        for column, label in enumerate(shared, start=1)
    )
    rows = [f'<row r="1">{cells}</row>']
    for line, (label, row) in enumerate(counts.items(), start=2):
        body = f'<c r="{_reference(1, line)}" t="s"><v>{number_of[label]}</v></c>' + "".join(
            f'<c r="{_reference(column, line)}"><v>{row[one]}</v></c>'
            for column, one in enumerate(labels, start=2)
            if one in row
        )
        rows.append(f'<row r="{line}">{body}</row>')
    namespace = SPREADSHEET_NS[1:-1]
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr(
            "xl/workbook.xml",
            f'<workbook xmlns="{namespace}"><sheets><sheet name="Overhangs"/></sheets></workbook>',
        )
        archive.writestr(
            "xl/sharedStrings.xml",
            f'<sst xmlns="{namespace}">'
            + "".join(f"<si><t>{label}</t></si>" for label in shared)
            + "</sst>",
        )
        archive.writestr(
            "xl/worksheets/sheet1.xml",
            f'<worksheet xmlns="{namespace}"><sheetData>{"".join(rows)}</sheetData></worksheet>',
        )
    return buffer.getvalue()


def _reference(column: int, row: int) -> str:
    """A cell reference, for the few columns these fixtures have."""
    return f"{chr(ord('A') + column - 1)}{row}"


@pytest.fixture
def profile_path(tmp_path: Path) -> Path:
    """The fixture matrix on disk, named as the archive names one."""
    path = tmp_path / f"{NAMED}.csv"
    path.write_text(comma_separated(), encoding="utf-8")
    return path


def test_a_comma_separated_matrix_is_read_into_counts(profile_path: Path) -> None:
    profile = read_profile(profile_path)

    assert profile.count("AAAA", "TTTT") == 600
    assert profile.count("GGAA", "TTCC") == 500
    assert profile.overhang_length == 4
    assert profile.observations == OBSERVATIONS


def test_a_workbook_reads_into_the_same_counts_as_the_comma_separated_form(tmp_path: Path) -> None:
    path = tmp_path / f"{NAMED}.xlsx"
    path.write_bytes(workbook())

    assert (
        read_profile(path).counts
        == read_profile(_written(tmp_path / f"{NAMED}.csv", comma_separated())).counts
    )


def test_a_pair_never_observed_reads_as_no_ligations(profile_path: Path) -> None:
    assert read_profile(profile_path).count("TTTT", "GGAA") == 0


def test_counts_are_normalised_to_a_hundred_thousand_events(profile_path: Path) -> None:
    profile = read_profile(profile_path)

    assert profile.normalised("AAAA", "TTTT") == pytest.approx(100_000 * 600 / OBSERVATIONS)


def test_the_conditions_are_read_from_the_name_the_archive_gives_a_file(profile_path: Path) -> None:
    profile = read_profile(profile_path)

    assert profile.conditions == "T4 DNA ligase, 18 h at 25 °C"
    assert NAMED in profile.source


def test_a_file_whose_name_states_no_conditions_says_they_are_unstated(tmp_path: Path) -> None:
    profile = read_profile(_written(tmp_path / "mine.csv", comma_separated()))

    assert profile.conditions == ""
    assert "not stated" in profile.source


def test_the_caller_can_state_the_conditions_the_file_name_does_not(tmp_path: Path) -> None:
    path = _written(tmp_path / "mine.csv", comma_separated())

    profile = read_profile(path, conditions="T4 DNA ligase, 1 h at 37 °C")

    assert "1 h at 37 °C" in profile.source


def test_a_comma_separated_file_that_is_not_a_matrix_is_refused(tmp_path: Path) -> None:
    path = _written(
        tmp_path / "fragments.csv", "Fragment #,Sequence,Length (bp)\n1,ATGCATGC,8\n2,GGCC,4\n"
    )

    with pytest.raises(ValueError, match="header row of overhang labels"):
        read_profile(path)


def test_a_workbook_that_is_not_a_matrix_is_refused_rather_than_raising(tmp_path: Path) -> None:
    path = tmp_path / "assemblies.xlsx"
    path.write_bytes(workbook(labels=("Count", "Size"), counts={"A:B:C": {"Count": 9, "Size": 12}}))

    with pytest.raises(ValueError, match="not a ligation count matrix"):
        read_profile(path)


def test_a_file_of_a_kind_the_reader_does_not_know_is_refused(tmp_path: Path) -> None:
    path = _written(tmp_path / "matrix.txt", comma_separated())

    with pytest.raises(ValueError, match=r"\.xlsx"):
        read_profile(path)


def test_a_matrix_holding_no_counts_is_refused(tmp_path: Path) -> None:
    path = _written(tmp_path / "empty.csv", "Overhang,AAAA,TTTT\n")

    with pytest.raises(ValueError, match="no counts"):
        read_profile(path)


def test_an_enzyme_with_no_matrix_of_its_own_is_scored_on_the_profile(profile_path: Path) -> None:
    profile = read_profile(profile_path)

    report = fidelity(("AAAA", "GGAA"), "PaqCI", profile=profile)

    assert report.measured
    assert not report.enzyme_specific
    assert report.value > 0.99
    assert len(report.ligations) == 2


def test_the_report_says_a_profile_is_not_a_measurement_of_the_enzyme(profile_path: Path) -> None:
    report = fidelity(("AAAA", "GGAA"), "PaqCI", profile=read_profile(profile_path))

    assert NAMED in report.source
    assert "25 °C" in report.source
    assert "PaqCI" in report.label
    assert "not specific" in report.label


def test_without_a_profile_an_unmeasured_enzyme_is_still_scored_by_the_rules() -> None:
    report = fidelity(("AATG", "GCTT", "TACA"), "PaqCI")

    assert not report.measured
    assert report.enzyme_specific
    assert report.label == "rule-based estimate"
    assert "rule" in report.source


def test_an_enzyme_with_its_own_matrix_keeps_it_when_a_profile_is_there(profile_path: Path) -> None:
    profile = read_profile(profile_path)
    shipped = ligation_matrix("BsaI")
    assert shipped is not None

    report = fidelity(("AAAA", "GGAA"), "BsaI", profile=profile)

    assert report.enzyme_specific
    assert report.source == shipped.source


def test_a_caller_can_prefer_the_profile_over_the_shipped_matrix(profile_path: Path) -> None:
    profile = read_profile(profile_path)

    report = fidelity(("AAAA", "GGAA"), "BsaI", profile=profile, prefer_profile=True)

    assert not report.enzyme_specific
    assert NAMED in report.source


def test_a_profile_of_the_wrong_overhang_length_for_the_enzyme_is_refused(
    profile_path: Path,
) -> None:
    with pytest.raises(ValueError, match="3"):
        fidelity(("AAA", "GGA"), "BspQI", profile=read_profile(profile_path))


def test_a_design_ranks_free_candidates_by_the_profile_where_no_matrix_exists(
    profile_path: Path,
) -> None:
    designed = design_overhangs(
        [Junction("left"), Junction("right")], "PaqCI", profile=read_profile(profile_path)
    )

    assert designed.overhangs[0] == "GGAA"
    assert designed.fidelity.measured
    assert not designed.fidelity.enzyme_specific


def test_the_pipeline_scores_an_unmeasured_enzyme_on_the_profile_and_says_so(
    profile_path: Path,
) -> None:
    made = plan_assembly(DATA / "pUC19.dna", DATA / "GFP.dna", enzyme="PaqCI", profile=profile_path)

    scored = made.overhangs.fidelity
    assert scored.measured
    assert not scored.enzyme_specific
    protocol = made.protocol()
    assert "not specific to PaqCI" in protocol.overview["Fidelity"]
    assert NAMED in " ".join(one.text for one in protocol.references)


def test_without_a_matrix_the_pipeline_scores_that_enzyme_as_it_did_before() -> None:
    made = plan_assembly(DATA / "pUC19.dna", DATA / "GFP.dna", enzyme="PaqCI")

    assert not made.overhangs.fidelity.measured
    assert "rule-based estimate" in made.protocol().overview["Fidelity"]


def test_the_command_line_takes_the_matrix_as_an_option(profile_path: Path, tmp_path: Path) -> None:
    result = _run(tmp_path / "given", "--enzyme", "PaqCI", "--ligase-matrix", str(profile_path))

    assert result.exit_code == 0, result.output
    assert "not specific to PaqCI" in result.output


def test_the_command_line_takes_the_matrix_from_the_environment(
    profile_path: Path, tmp_path: Path
) -> None:
    result = _run(tmp_path / "env", "--enzyme", "PaqCI", env={LIGASE_MATRIX_ENV: str(profile_path)})

    assert result.exit_code == 0, result.output
    assert "not specific to PaqCI" in result.output


def test_the_command_line_refuses_a_file_that_is_not_a_matrix(tmp_path: Path) -> None:
    path = _written(tmp_path / "fragments.csv", "Fragment #,Sequence\n1,ATGC\n")

    result = _run(tmp_path / "run", "--ligase-matrix", str(path))

    assert result.exit_code == 1
    assert "not a ligation count matrix" in result.output


def _run(out: Path, *arguments: str, env: Mapping[str, str] | None = None):
    """Plan the fixture assembly on the command line, writing the three outputs into OUT."""
    return CliRunner().invoke(
        app,
        [
            "goldengate",
            "plan",
            str(DATA / "pUC19.dna"),
            str(DATA / "GFP.dna"),
            "--out",
            str(out),
            *arguments,
        ],
        env=dict(env or {}),
    )


def _written(path: Path, text: str) -> Path:
    """Write a fixture file and return where it went."""
    path.write_text(text, encoding="utf-8")
    return path

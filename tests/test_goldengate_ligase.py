"""Reading a ligase fidelity matrix the user holds on their own disk.

The archive such a file comes from is CC BY-NC 4.0 and none of it is in this repository. Every
matrix here is four overhangs written by this file, with plausible integers in the cells.
"""

import io
import zipfile
from collections.abc import Mapping, Sequence
from pathlib import Path

import pytest

from liulab_mbio.goldengate.ligase import SPREADSHEET_NS, read_profile

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


def _written(path: Path, text: str) -> Path:
    """Write a fixture file and return where it went."""
    path.write_text(text, encoding="utf-8")
    return path

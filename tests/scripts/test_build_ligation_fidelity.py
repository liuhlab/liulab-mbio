"""The ligation fidelity build script, on a workbook built here rather than fetched.

The workbook reader it calls is the package's, `mbio.ligase`, so the sheet a matrix of
several is read from is pinned here too, on a workbook built here the same way.
"""

import importlib.util
import io
import json
import sys
import zipfile
from collections.abc import Mapping, Sequence
from itertools import product
from pathlib import Path
from types import ModuleType

import pytest

from mbio.ligase import (
    RELATIONSHIP_NS,
    SPREADSHEET_NS,
    column_index,
    read_profile,
    read_workbook,
)
from mbio.sequence import reverse_complement

REPO = Path(__file__).resolve().parents[2]


def _script() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "build_ligation_fidelity", REPO / "scripts/build_ligation_fidelity.py"
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # Registered before it runs: `dataclasses` resolves annotations through `sys.modules`.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


build = _script()

SHEET = "S1 Table. BsaI-HFv2"

# Two rows of a matrix shaped as PLOS serves it: labels in the shared string table, counts as
# numbers, and the zeros a sparse file drops.
LABELS = ("Overhang", "AAAA", "AAAC", "TTTT")
CELLS: Mapping[str, tuple[str | None, int]] = {
    "A1": ("s", 0), "B1": ("s", 1), "C1": ("s", 2), "D1": ("s", 3),
    "A2": ("s", 3), "B2": (None, 635), "C2": (None, 0), "D2": (None, 4),
    "A3": ("s", 1), "B3": (None, 0), "D3": (None, 635),
}  # fmt: skip


def workbook(
    sheet: str = SHEET,
    labels: tuple[str, ...] = LABELS,
    cells: Mapping[str, tuple[str | None, int]] | None = None,
) -> bytes:
    """Write the three parts of an `.xlsx` the script reads, and nothing else."""
    cells = CELLS if cells is None else cells
    shared = "".join(f"<si><t>{label}</t></si>" for label in labels)
    rows = ""
    for number in sorted({int(reference[1:]) for reference in cells}):
        body = "".join(
            f'<c r="{reference}"{"" if kind is None else f' t="{kind}"'}><v>{value}</v></c>'
            for reference, (kind, value) in cells.items()
            if int(reference[1:]) == number
        )
        rows += f'<row r="{number}">{body}</row>'
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr(
            "xl/workbook.xml",
            f'<workbook xmlns="{SPREADSHEET_NS[1:-1]}"><sheets><sheet name="{sheet}"/>'
            "</sheets></workbook>",
        )
        archive.writestr(
            "xl/sharedStrings.xml", f'<sst xmlns="{SPREADSHEET_NS[1:-1]}">{shared}</sst>'
        )
        archive.writestr(
            "xl/worksheets/sheet1.xml",
            f'<worksheet xmlns="{SPREADSHEET_NS[1:-1]}"><sheetData>{rows}</sheetData></worksheet>',
        )
    return buffer.getvalue()


def test_a_workbook_is_read_into_rows_of_counts_without_a_spreadsheet_library() -> None:
    sheet, counts = read_workbook(workbook())

    assert sheet == SHEET
    assert counts == {"TTTT": {"AAAA": 635, "TTTT": 4}, "AAAA": {"TTTT": 635}}


def test_a_pair_never_observed_is_absent_rather_than_zero() -> None:
    _, counts = read_workbook(workbook())

    assert "AAAC" not in counts["TTTT"]
    assert "AAAA" not in counts["AAAA"]


def test_a_column_is_placed_by_its_letters_and_not_by_its_order_in_the_row() -> None:
    # A row whose middle cell is missing entirely, which is what a zero looks like in some rows.
    cells = {"A1": ("s", 0), "B1": ("s", 1), "C1": ("s", 2), "D1": ("s", 3),
             "A2": ("s", 3), "D2": (None, 7)}  # fmt: skip

    _, counts = read_workbook(workbook(cells=cells))

    assert counts == {"TTTT": {"TTTT": 7}}


def test_columns_past_the_alphabet_keep_their_place() -> None:
    assert (column_index("A1"), column_index("Z9"), column_index("AA1")) == (1, 26, 27)
    assert column_index("IW257") == 257


def test_a_cell_of_an_unread_type_is_refused_rather_than_guessed() -> None:
    cells = {"A1": ("s", 0), "B1": ("s", 1), "A2": ("s", 1), "B2": ("e", 0)}

    with pytest.raises(ValueError, match="unread type"):
        read_workbook(workbook(cells=cells))


def test_a_sheet_named_after_another_enzyme_is_refused() -> None:
    table = build.TABLES[0]

    with pytest.raises(ValueError, match="BsaI-HFv2"):
        build.matrix(table, workbook(sheet="Table S5. SapI"))


def test_a_label_that_is_not_an_overhang_is_refused() -> None:
    table = build.TABLES[0]

    with pytest.raises(ValueError, match="not an overhang"):
        build.matrix(table, workbook(labels=("Overhang", "AAAA", "AAAC", "TTT")))


def test_the_matrix_records_which_table_and_which_product_it_came_from() -> None:
    one = build.matrix(build.TABLES[0], workbook())

    assert (one["enzyme"], one["product"], one["table"]) == ("BsaI", "BsaI-HFv2", "S1 Table")
    assert one["file"] == "pone.0238592.s001.xlsx"
    assert (one["overhang_length"], one["cycling_celsius"]) == (4, [37, 16])
    assert one["observations"] == 635 + 4 + 635


def blobs() -> dict[str, bytes]:
    """One workbook per table, labelled with overhangs of the length that enzyme leaves."""
    return {
        table.enzyme: workbook(
            sheet=f"{table.table}. {table.product}",
            labels=tuple(
                label[: table.overhang_length] if index else label
                for index, label in enumerate(LABELS)
            ),
        )
        for table in build.TABLES
    }


def test_the_data_file_carries_the_licence_and_the_citation_it_ships_under() -> None:
    document = build.build(blobs())

    assert "CC BY 4.0" in document["source"]["licence"]
    assert "Pryor" in document["source"]["citation"]
    assert "reverse complement" in document["format"]


def test_one_row_per_line_is_still_json() -> None:
    document = build.build(blobs())

    rendered = build.dumps(document)

    assert json.loads(rendered) == document
    assert '"TTTT": {"AAAA":635,"TTTT":4}' in rendered


# Bilotti 2022's supplement is one workbook of eight matrices, one per ligase and buffer, so a
# reader taking the first sheet of a file answers a question about T7 with T4's numbers.
BILOTTI = (
    "File S1. T4",
    "File S2. T7",
    "File S3. hLig3",
    "File S4. T3",
    "File S5. PBCV-1",
    "File S6. T4 PEG",
    "File S7. T7 PEG",
    "File S8. hLig3 PEG",
)
T7 = BILOTTI[1]

#: The two numbers that sheet is recognised by: every four-base overhang, and its own total.
T7_OVERHANGS = 256
T7_OBSERVATIONS = 338_272


def bilotti() -> bytes:
    """A workbook shaped as that supplement is: eight sheets, one count matrix each.

    The T7 sheet carries every four-base overhang and the rest a two-overhang stub of its own
    size, so a sheet read in place of another is plain. The parts are numbered against the
    order the sheets are listed in, which an `.xlsx` allows, so only the relationships resolve
    them.
    """
    return sheets_workbook(
        {
            name: _overhang_matrix()
            if name == T7
            else {"AAAA": {"TTTT": number}, "TTTT": {"AAAA": number}}
            for number, name in enumerate(BILOTTI, start=1)
        }
    )


def sheets_workbook(sheets: Mapping[str, Mapping[str, Mapping[str, int]]]) -> bytes:
    """Write a workbook holding one count matrix per named sheet."""
    labels: dict[str, int] = {"Overhang": 0}
    for counts in sheets.values():
        for row, columns in counts.items():
            for label in (row, *columns):
                labels.setdefault(label, len(labels))
    parts = [f"xl/worksheets/sheet{len(sheets) - number}.xml" for number in range(len(sheets))]
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("xl/workbook.xml", _workbook_xml(tuple(sheets)))
        archive.writestr("xl/_rels/workbook.xml.rels", _rels_xml(parts))
        archive.writestr(
            "xl/sharedStrings.xml",
            f'<sst xmlns="{SPREADSHEET_NS[1:-1]}">'
            + "".join(f"<si><t>{label}</t></si>" for label in labels)
            + "</sst>",
        )
        for part, counts in zip(parts, sheets.values(), strict=True):
            archive.writestr(part, _worksheet_xml(counts, labels))
    return buffer.getvalue()


def _overhang_matrix() -> dict[str, dict[str, int]]:
    """Every four-base overhang against its own partner, totalling the T7 sheet's own count."""
    overhangs = ["".join(bases) for bases in product("ACGT", repeat=4)]
    each, extra = divmod(T7_OBSERVATIONS, len(overhangs))
    counts = {one: {reverse_complement(one): each} for one in overhangs}
    counts[overhangs[0]][reverse_complement(overhangs[0])] += extra
    return counts


def _workbook_xml(names: Sequence[str]) -> str:
    """The part listing the sheets, each pointing at its own part by relationship."""
    sheets = "".join(
        f'<sheet name="{name}" sheetId="{number}" r:id="rId{number}"/>'
        for number, name in enumerate(names, start=1)
    )
    return (
        f'<workbook xmlns="{SPREADSHEET_NS[1:-1]}" xmlns:r="{RELATIONSHIP_NS[1:-1]}">'
        f"<sheets>{sheets}</sheets></workbook>"
    )


def _rels_xml(parts: Sequence[str]) -> str:
    """The part saying which file each relationship id points at."""
    package = "http://schemas.openxmlformats.org/package/2006/relationships"
    kind = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet"
    entries = "".join(
        f'<Relationship Id="rId{number}" Type="{kind}" Target="{part.removeprefix("xl/")}"/>'
        for number, part in enumerate(parts, start=1)
    )
    return f'<Relationships xmlns="{package}">{entries}</Relationships>'


def _worksheet_xml(counts: Mapping[str, Mapping[str, int]], labels: Mapping[str, int]) -> str:
    """One sheet: the labels across the top and down the side, a count in each cell."""
    columns = sorted({label for row, one in counts.items() for label in (row, *one)})
    place = {label: number for number, label in enumerate(columns, start=2)}
    header = "".join(
        f'<c r="{_reference(number, 1)}" t="s"><v>{labels[label]}</v></c>'
        for number, label in enumerate(("Overhang", *columns), start=1)
    )
    rows = [f'<row r="1">{header}</row>']
    for line, (label, row) in enumerate(counts.items(), start=2):
        body = f'<c r="A{line}" t="s"><v>{labels[label]}</v></c>' + "".join(
            f'<c r="{_reference(place[one], line)}"><v>{value}</v></c>'
            for one, value in row.items()
        )
        rows.append(f'<row r="{line}">{body}</row>')
    return (
        f'<worksheet xmlns="{SPREADSHEET_NS[1:-1]}"><sheetData>{"".join(rows)}</sheetData>'
        "</worksheet>"
    )


def _reference(column: int, row: int) -> str:
    """A cell reference, for a matrix wider than the alphabet."""
    letters = ""
    while column:
        column, remainder = divmod(column - 1, 26)
        letters = chr(ord("A") + remainder) + letters
    return f"{letters}{row}"


def test_a_workbook_of_several_matrices_is_refused_until_a_sheet_is_named(tmp_path: Path) -> None:
    path = tmp_path / "File S1_NAR.xlsx"
    path.write_bytes(bilotti())

    with pytest.raises(ValueError, match="none was asked for") as refused:
        read_profile(path)

    assert all(name in str(refused.value) for name in BILOTTI)


def test_the_sheet_named_is_the_one_read_and_not_the_first(tmp_path: Path) -> None:
    path = tmp_path / "File S1_NAR.xlsx"
    path.write_bytes(bilotti())

    profile = read_profile(path, sheet=T7)

    assert len(profile.overhangs) == T7_OVERHANGS
    assert profile.observations == T7_OBSERVATIONS


def test_a_sheet_names_the_conditions_a_file_name_does_not(tmp_path: Path) -> None:
    path = tmp_path / "File S1_NAR.xlsx"
    path.write_bytes(bilotti())

    assert read_profile(path, sheet=T7).conditions == T7


def test_a_sheet_is_found_through_its_relationship_and_not_its_part_number() -> None:
    name, counts = read_workbook(bilotti(), sheet=0)

    assert name == BILOTTI[0]
    assert counts == {"AAAA": {"TTTT": 1}, "TTTT": {"AAAA": 1}}


def test_a_sheet_the_workbook_does_not_hold_is_refused_by_name() -> None:
    with pytest.raises(ValueError, match="no sheet named"):
        read_workbook(bilotti(), sheet="File S9. T7")

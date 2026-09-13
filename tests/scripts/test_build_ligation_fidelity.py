"""The ligation fidelity build script, on a workbook built here rather than fetched.

The workbook reader it calls is the package's, `liulab_mbio.goldengate.ligase`.
"""

import importlib.util
import io
import json
import sys
import zipfile
from collections.abc import Mapping
from pathlib import Path
from types import ModuleType

import pytest

from liulab_mbio.goldengate.ligase import SPREADSHEET_NS, column_index, read_workbook

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

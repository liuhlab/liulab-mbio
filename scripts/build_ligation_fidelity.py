#!/usr/bin/env python3
"""Rebuild the shipped ligation fidelity matrices from Pryor et al. 2020's S1-S5 Tables.

    pixi run python scripts/build_ligation_fidelity.py [--from DIR] [--out PATH]

The tables are the supplementary files of an open access PLOS One article under CC BY 4.0, so
they may be redistributed with attribution. `docs/research/ligation-fidelity.md` records the
licence, the axis convention and what each table measured.

Without `--from` the five workbooks are downloaded from `journals.plos.org`; with it they are
read from a directory holding files named as PLOS serves them. Each workbook is a zip of XML,
so `zipfile` and `xml.etree` read it and no spreadsheet library is needed.
"""

from __future__ import annotations

import argparse
import io
import json
import re
import sys
import urllib.request
import zipfile
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "src/liulab_mbio/data/ligation_fidelity.json"

DOI = "10.1371/journal.pone.0238592"

#: PLOS names a supplementary file after the article, without the `journal.` the DOI carries.
STEM = "pone.0238592"

#: PLOS serves a supplementary file by article DOI and file suffix, redirecting to storage.
SUPPLEMENTARY = (
    "https://journals.plos.org/plosone/article/file?id={doi}.{suffix}&type=supplementary"
)

#: The spreadsheet XML namespace, and the column letters of a cell reference.
_MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_COLUMN = re.compile(r"[A-Z]+")

#: What the file says about itself, so a reader need not find this script first.
FORMAT = (
    "Sparse integer count matrix. `counts` maps a top-strand overhang to the bottom-strand "
    "overhangs it was observed ligating to and how often; a pair never observed is absent "
    "rather than zero. Both axes are written 5' to 3', so the Watson-Crick entry of a row is "
    "the column spelling its reverse complement, not the column of the same name."
)

LICENCE = "CC BY 4.0, https://creativecommons.org/licenses/by/4.0/"

TERMS = (
    "This is an open access article distributed under the terms of the Creative Commons "
    "Attribution License, which permits unrestricted use, distribution, and reproduction in "
    "any medium, provided the original author and source are credited."
)

CITATION = (
    "Pryor, J.M., Potapov, V., Kucera, R.B., Bilotti, K., Cantor, E.J. and Lohman, G.J.S. "
    "(2020) Enabling one-pot Golden Gate assemblies of unprecedented complexity using "
    "data-optimized assembly design. PLoS One 15(9): e0238592."
)


@dataclass(frozen=True)
class Table:
    """One supplementary table: the enzyme it measured, and how that reaction was cycled."""

    enzyme: str
    product: str
    table: str
    suffix: str
    overhang_length: int
    cycling_celsius: tuple[int, int]


#: The five per-enzyme tables. `product` is what the sheet inside each workbook is named after,
#: which is the guard against a mis-mapped table.
TABLES: tuple[Table, ...] = (
    Table("BsaI", "BsaI-HFv2", "S1 Table", "s001", 4, (37, 16)),
    Table("BsmBI", "BsmBI-v2", "S2 Table", "s002", 4, (42, 16)),
    Table("Esp3I", "Esp3I", "S3 Table", "s003", 4, (37, 16)),
    Table("BbsI", "BbsI-HF", "S4 Table", "s004", 4, (37, 16)),
    Table("SapI", "SapI", "S5 Table", "s005", 3, (37, 16)),
)


def filename(table: Table) -> str:
    """Return the name PLOS serves this table's workbook under."""
    return f"{STEM}.{table.suffix}.xlsx"


def url(table: Table) -> str:
    """Return the download address of one table."""
    return SUPPLEMENTARY.format(doi=DOI, suffix=table.suffix)


def fetch(table: Table) -> bytes:
    """Download one workbook."""
    with urllib.request.urlopen(url(table), timeout=300) as response:
        return response.read()


def read_workbook(blob: bytes) -> tuple[str, dict[str, dict[str, int]]]:
    """Read one workbook into its sheet name and a sparse count matrix.

    The first row and the first column hold the overhang labels; every other cell is an
    observation count, and a zero is dropped.

    Raises
    ------
    ValueError
        If a cell is neither a shared string nor a number, which is not a shape these
        workbooks have and so is a sign the file is not the one expected.
    """
    with zipfile.ZipFile(io.BytesIO(blob)) as archive:
        sheet_name = _sheet_name(archive)
        shared = _shared_strings(archive)
        rows = _rows(archive, shared)
    header, *body = rows
    # Keyed by column number, because a row drops the cells it has no count for.
    columns = {_index(reference): str(label) for reference, label in header}
    counts: dict[str, dict[str, int]] = {}
    for cells in body:
        values = {_index(reference): value for reference, value in cells}
        counts[str(values[1])] = {
            columns[index]: value
            for index, value in values.items()
            if index != 1 and isinstance(value, int) and value
        }
    return sheet_name, counts


def _sheet_name(archive: zipfile.ZipFile) -> str:
    """Return the name of the workbook's first sheet."""
    workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
    sheets = workbook.find(f"{_MAIN}sheets")
    if sheets is None or len(sheets) == 0:
        raise ValueError("the workbook holds no sheet")
    return sheets[0].get("name", "")


def _shared_strings(archive: zipfile.ZipFile) -> list[str]:
    """Return the workbook's shared string table, which is where its labels live."""
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    table = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
    return ["".join(node.itertext()) for node in table]


def _rows(archive: zipfile.ZipFile, shared: list[str]) -> list[list[tuple[str, str | int]]]:
    """Return every non-empty cell of the first sheet, row by row, as reference and value."""
    sheet = ElementTree.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    data = sheet.find(f"{_MAIN}sheetData")
    if data is None:
        raise ValueError("the sheet holds no data")
    rows: list[list[tuple[str, str | int]]] = []
    for row in data:
        cells: list[tuple[str, str | int]] = []
        for cell in row:
            text = cell.findtext(f"{_MAIN}v")
            if text is None:
                continue
            kind = cell.get("t")
            if kind == "s":
                cells.append((cell.get("r", ""), shared[int(text)]))
            elif kind in (None, "n"):
                cells.append((cell.get("r", ""), int(float(text))))
            else:
                raise ValueError(f"cell {cell.get('r')} is of unread type {kind!r}")
        if cells:
            rows.append(cells)
    return rows


def _index(reference: str) -> int:
    """Return the 1-based column number of a cell reference such as ``"IW257"``."""
    letters = _COLUMN.match(reference)
    if letters is None:
        raise ValueError(f"{reference!r} is not a cell reference")
    number = 0
    for letter in letters.group():
        number = number * 26 + ord(letter) - ord("A") + 1
    return number


def matrix(table: Table, blob: bytes) -> dict[str, Any]:
    """Read one workbook into the record the data file carries for it.

    Raises
    ------
    ValueError
        If the sheet is not named after the product this table measured, or if a label is not
        an overhang of the length that enzyme leaves.
    """
    sheet_name, counts = read_workbook(blob)
    if table.product.lower() not in sheet_name.lower().replace(" ", ""):
        raise ValueError(
            f"{table.table} should measure {table.product}, and its sheet is named {sheet_name!r}"
        )
    for label in counts:
        if len(label) != table.overhang_length or set(label) - set("ACGT"):
            raise ValueError(f"{table.table} carries {label!r}, which is not an overhang")
    return {
        "enzyme": table.enzyme,
        "product": table.product,
        "table": table.table,
        "file": filename(table),
        "sheet": sheet_name,
        "overhang_length": table.overhang_length,
        "cycling_celsius": list(table.cycling_celsius),
        "observations": sum(sum(row.values()) for row in counts.values()),
        "counts": counts,
    }


def build(blobs: Mapping[str, bytes], *, built: str | None = None) -> dict[str, Any]:
    """Build the whole data file from one workbook per table, keyed by enzyme name."""
    return {
        "built": built or datetime.now(UTC).date().isoformat(),
        "source": {
            "citation": CITATION,
            "doi": DOI,
            "url": SUPPLEMENTARY.format(doi=DOI, suffix="sNNN"),
            "retrieved": built or datetime.now(UTC).date().isoformat(),
            "copyright": "Copyright (c) 2020 Pryor et al.",
            "licence": LICENCE,
            "terms": TERMS,
        },
        "format": FORMAT,
        "matrices": [matrix(table, blobs[table.enzyme]) for table in TABLES],
    }


def dumps(document: Mapping[str, Any]) -> str:
    """Render the data file with one overhang row per line, so a diff reads a row at a time."""
    head = [_field(key, value, "  ") for key, value in document.items() if key != "matrices"]
    blocks = []
    for one in document["matrices"]:
        fields = [_field(key, value, "    ") for key, value in one.items() if key != "counts"]
        rows = ",\n".join(
            f"      {json.dumps(row)}: {json.dumps(counts, separators=(',', ':'))}"
            for row, counts in one["counts"].items()
        )
        fields.append('    "counts": {\n' + rows + "\n    }")
        blocks.append("  {\n" + ",\n".join(fields) + "\n  }")
    head.append('  "matrices": [\n' + ",\n".join(blocks) + "\n  ]")
    return "{\n" + ",\n".join(head) + "\n}\n"


def _field(key: str, value: Any, indent: str) -> str:
    """Render one key and value, indenting the lines a nested value runs onto."""
    body = json.dumps(value, indent=2, ensure_ascii=False)
    return f"{indent}{json.dumps(key)}: " + body.replace("\n", "\n" + indent)


def main(argv: list[str] | None = None) -> int:
    """Fetch or read the five tables and write the data file."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--from",
        dest="source",
        type=Path,
        help="a directory of already downloaded workbooks, named as PLOS serves them",
    )
    parser.add_argument("--out", type=Path, default=DATA, help="where to write the data file")
    arguments = parser.parse_args(argv)

    blobs: dict[str, bytes] = {}
    for table in TABLES:
        if arguments.source is None:
            print(f"downloading {table.table} ({table.product})", file=sys.stderr)
            blobs[table.enzyme] = fetch(table)
        else:
            blobs[table.enzyme] = (arguments.source / filename(table)).read_bytes()

    document = build(blobs)
    arguments.out.write_text(dumps(document), encoding="utf-8")
    print(f"wrote {arguments.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

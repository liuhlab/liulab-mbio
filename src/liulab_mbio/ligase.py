"""A ligase's own overhang profile, read from a matrix file the user already holds.

Potapov et al. 2018 profiled T4 DNA ligase across all 256 four-base overhangs, which covers the
Type IIS enzymes nobody has published a matrix for. That archive is CC BY-NC 4.0, so **none of
it is in this package**: `read_profile` reads a copy from the user's own disk and nothing is
redistributed. `docs/research/ligation-fidelity.md` says where the archive is and what it holds.

A profile belongs to the ligase and the conditions it was measured under, not to the Type IIS
enzyme, and `liulab_mbio.overhangs.fidelity` says so in the report it returns.

The two shapes the archive uses are both read here with the standard library alone: an `.xlsx`
is a zip of XML, so `zipfile` and `xml.etree` are enough, and a `.csv` is read by `csv`.
"""

import csv
import io
import os
import re
import zipfile
from collections.abc import Mapping
from dataclasses import KW_ONLY, dataclass, field
from pathlib import Path
from xml.etree import ElementTree

#: The spreadsheet XML namespace, and the column letters of a cell reference.
SPREADSHEET_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_COLUMN = re.compile(r"[A-Z]+")

#: What the archive names a matrix file after: the ligase, the incubation, the temperature.
_CONDITIONS = re.compile(r"(T4|T7)[_-]0*(\d+)h[_-](\d+)C", re.IGNORECASE)

#: What a file has to hold to be one of these, said in the refusal so a caller need not guess.
EXPECTED = (
    "a header row of overhang labels, the same labels down the first column, and a count in "
    "each cell, as an .xlsx or a .csv file"
)


@dataclass(frozen=True, slots=True)
class LigaseProfile:
    """How often each overhang pair was seen ligating, measured with one ligase.

    Parameters
    ----------
    path
        The file it was read from.
    conditions
        The ligase, the incubation and the temperature. Empty where nothing states them.
    overhang_length
        How many bases the overhangs on both axes have.
    observations
        Every ligation event counted.
    counts
        Top-strand overhang to the bottom-strand overhangs it was seen ligating to. Both are
        written 5' to 3', so a row pairs with the column spelling its reverse complement.
    """

    path: Path
    _: KW_ONLY
    conditions: str
    overhang_length: int
    observations: int
    counts: Mapping[str, Mapping[str, int]] = field(hash=False)

    @property
    def overhangs(self) -> tuple[str, ...]:
        """Every overhang the measurement covers."""
        return tuple(self.counts)

    @property
    def source(self) -> str:
        """Where a report should say this number came from."""
        stated = self.conditions or "conditions not stated by the file name"
        return f"{stated}, read from {self.path.name}"

    def count(self, top: str, bottom: str) -> int:
        """How often a top-strand overhang was seen ligating to a bottom-strand one."""
        return self.counts.get(top.upper(), {}).get(bottom.upper(), 0)

    def normalised(self, top: str, bottom: str) -> float:
        """Return the count per 100,000 ligation events, the scale NEB's thresholds use."""
        return 100_000 * self.count(top, bottom) / self.observations


def read_profile(path: str | os.PathLike[str], *, conditions: str = "") -> LigaseProfile:
    """Read a ligation count matrix the user holds, as an `.xlsx` or a `.csv` file.

    Parameters
    ----------
    path
        The matrix. Its conditions are read from its name where that is named the way the
        archive names one, such as ``FileS03_T4_18h_25C.xlsx``.
    conditions
        What the reaction was, for a file whose name does not say, or to correct one that does.

    Raises
    ------
    ValueError
        If the file is not a count matrix, saying what one is.

    Examples
    --------
    >>> read_profile("FileS03_T4_18h_25C.xlsx").conditions  # doctest: +SKIP
    'T4 DNA ligase, 18 h at 25 °C'
    """
    one = Path(path)
    try:
        counts, length = _matrix(one)
    except (KeyError, ValueError, zipfile.BadZipFile, ElementTree.ParseError) as error:
        raise ValueError(
            f"{one.name} is not a ligation count matrix ({error}); expected {EXPECTED}"
        ) from error
    return LigaseProfile(
        one,
        conditions=conditions or _conditions(one.name),
        overhang_length=length,
        observations=sum(sum(row.values()) for row in counts.values()),
        counts=counts,
    )


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
    columns = {column_index(reference): str(label) for reference, label in header}
    counts: dict[str, dict[str, int]] = {}
    for cells in body:
        values = {column_index(reference): value for reference, value in cells}
        counts[str(values[1])] = {
            columns[index]: value
            for index, value in values.items()
            if index != 1 and isinstance(value, int) and value
        }
    return sheet_name, counts


def read_csv(text: str) -> dict[str, dict[str, int]]:
    """Read a comma-separated count matrix into the sparse form `read_workbook` returns.

    Raises
    ------
    ValueError
        If the file holds no row, or a cell holds something that is not a count.
    """
    rows = [row for row in csv.reader(io.StringIO(text)) if any(cell.strip() for cell in row)]
    if not rows:
        raise ValueError("it holds no rows")
    header, *body = rows
    columns = [cell.strip() for cell in header[1:]]
    counts: dict[str, dict[str, int]] = {}
    for cells in body:
        counts[cells[0].strip()] = {
            column: number
            for column, value in zip(columns, cells[1:], strict=False)
            if (number := _number(value))
        }
    return counts


def column_index(reference: str) -> int:
    """Return the 1-based column number of a cell reference such as ``"IW257"``.

    Raises
    ------
    ValueError
        If the reference carries no column letters.

    Examples
    --------
    >>> column_index("A1"), column_index("AA1")
    (1, 27)
    """
    letters = _COLUMN.match(reference)
    if letters is None:
        raise ValueError(f"{reference!r} is not a cell reference")
    number = 0
    for letter in letters.group():
        number = number * 26 + ord(letter) - ord("A") + 1
    return number


def _matrix(path: Path) -> tuple[dict[str, dict[str, int]], int]:
    """Read the counts of one file, and the length of the overhangs they are keyed by.

    Raises
    ------
    ValueError
        If the file is of a kind this does not read, holds no count, or is labelled with
        anything but overhangs of one length.
    """
    suffix = path.suffix.lower()
    if suffix == ".xlsx":
        _, counts = read_workbook(path.read_bytes())
    elif suffix == ".csv":
        counts = read_csv(path.read_text(encoding="utf-8-sig"))
    else:
        raise ValueError(f"a {suffix or 'suffixless'} file is neither .xlsx nor .csv")
    if not any(row for row in counts.values()):
        raise ValueError("it holds no counts at all")
    labels = {*counts, *(column for row in counts.values() for column in row)}
    lengths = {len(label) for label in labels}
    if len(lengths) != 1 or any(set(label) - set("ACGT") for label in labels):
        odd = sorted(label for label in labels if set(label) - set("ACGT")) or sorted(labels)
        raise ValueError(f"it is labelled {odd[0]!r}, which is not an overhang of one length")
    return counts, lengths.pop()


def _number(value: str) -> int:
    """Read one cell as a count, an empty cell counting as none.

    Raises
    ------
    ValueError
        If the cell holds something that is not a number.
    """
    text = value.strip()
    if not text:
        return 0
    try:
        return round(float(text))
    except ValueError:
        raise ValueError(f"{text!r} is not a count") from None


def _conditions(name: str) -> str:
    """Read the ligase and the reaction off a file named the way the archive names one."""
    found = _CONDITIONS.search(name)
    if found is None:
        return ""
    ligase, hours, celsius = found.groups()
    return f"{ligase.upper()} DNA ligase, {int(hours)} h at {celsius} °C"


def _sheet_name(archive: zipfile.ZipFile) -> str:
    """Return the name of the workbook's first sheet.

    Raises
    ------
    ValueError
        If the workbook holds no sheet.
    """
    workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
    sheets = workbook.find(f"{SPREADSHEET_NS}sheets")
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
    """Return every non-empty cell of the first sheet, row by row, as reference and value.

    Raises
    ------
    ValueError
        If the sheet holds no data, or a cell is of a type this does not read.
    """
    sheet = ElementTree.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    data = sheet.find(f"{SPREADSHEET_NS}sheetData")
    if data is None:
        raise ValueError("the sheet holds no data")
    rows: list[list[tuple[str, str | int]]] = []
    for row in data:
        cells: list[tuple[str, str | int]] = []
        for cell in row:
            text = cell.findtext(f"{SPREADSHEET_NS}v")
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

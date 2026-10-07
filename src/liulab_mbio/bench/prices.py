"""A banded price record the user holds, the charge it gives a quantity, and the bill.

Price is the one material parameter this package never ships: a tariff is the user's, it goes
stale, and a figure quoted here would be ours. `read_prices` reads a copy from the user's own
disk, exactly as `liulab_mbio.ligase.read_profile` reads a matrix, and nothing is redistributed.

A quantity always computes, because it comes from the design. Money comes only from a record,
and where no row prices a key the money cell is a hole carrying the key and the quantity that
went unpriced. No estimate is ever written.

**Price steers no design.** The lookup reports headroom — how far a quantity sits from the
nearest band edge — and that is the whole of what a cliff gets. A design that changed with
whose price list was loaded would not be reproducible between labs.
"""

import csv
import os
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import KW_ONLY, dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Literal

from liulab_mbio.protocol.model import Bill, BillRow, Citation, Hole, Source

#: What a basis says the charge is for: flat inside the band, or amortised over a pack.
type Basis = Literal["per order", "per unit"]

#: The columns a price record holds, in any order. `bands` is empty for an item that is not
#: banded, and otherwise holds ``quantity low-high`` separated by semicolons. One cell and not a
#: column pair a quantity: ``docs/adr/0011-bench-model.md`` says why.
COLUMNS = ("key", "item", "bands", "charge", "basis", "currency")

#: Where a price record is read from when a caller names none. The package ships none; this
#: points at a copy the user holds.
PRICES_ENV = "LIULAB_MBIO_PRICES"

#: What a file has to hold to be one of these, said in the refusal so a caller need not guess.
EXPECTED = (
    f"a CSV file with the columns {', '.join(COLUMNS)}, one currency throughout, a charge that "
    "is a number, a basis of 'per order' or 'per unit', and each band written "
    "'quantity low-high' with both ends inclusive, the high end left empty for a top tier with "
    "no top"
)


@dataclass(frozen=True, slots=True)
class Band:
    """One quantity a row is bounded by, inclusive at both ends.

    Both ends are inclusive because that is how every vendor writes a tier, and a transcription
    slip here is a wrong price. A `high` of ``None`` is a top tier with no top, so nobody writes
    a sentinel for one. This is a quantity, not a span in a sequence, so
    ``docs/adr/0001-coordinates.md`` does not reach it.
    """

    quantity: str
    low: Decimal
    high: Decimal | None = None

    def __post_init__(self) -> None:
        """Refuse a band whose ends are the wrong way round."""
        if self.high is not None and self.high < self.low:
            raise ValueError(f"band {self.quantity}: {self.low} to {self.high} is backwards")

    def holds(self, value: Decimal) -> bool:
        """Return whether `value` falls in the band."""
        return self.low <= value and (self.high is None or value <= self.high)

    def __str__(self) -> str:
        """Say the band as the record's own cell writes it, the high end empty for no top."""
        return f"{self.quantity} {self.low}-{self.high if self.high is not None else ''}"


@dataclass(frozen=True, slots=True)
class PriceRow:
    """One row of a price record: what it prices, in what band, for how much, on what basis."""

    key: str
    charge: Decimal
    basis: Basis
    _: KW_ONLY
    item: str = ""
    bands: tuple[Band, ...] = ()

    def holds(self, quantities: Mapping[str, float]) -> bool:
        """Return whether every band of the row holds the quantity it names."""
        return all(
            band.quantity in quantities and band.holds(Decimal(str(quantities[band.quantity])))
            for band in self.bands
        )

    @property
    def locator(self) -> str:
        """Where in the record this row is, as a `Citation` says it."""
        bands = "; ".join(str(band) for band in self.bands)
        return f"{self.key} {bands}".strip() if bands else self.key


@dataclass(frozen=True, slots=True)
class Headroom:
    """How far a quantity sits from the edges of the band it fell in."""

    quantity: str
    value: Decimal
    band: Band

    @property
    def above(self) -> Decimal | None:
        """How much more the quantity may grow before the band changes, or no top at all."""
        return None if self.band.high is None else self.band.high - self.value

    def __str__(self) -> str:
        """Say the fact the way a protocol shows it."""
        gap = self.above
        if gap is None:
            return f"{self.value:f} {self.quantity}, no band above it"
        if gap == 0:
            return f"{self.value:f} {self.quantity}, no slack above it at all"
        return f"{self.value:f} {self.quantity}, {gap:f} below the next band"


@dataclass(frozen=True, slots=True)
class PriceRecord:
    """A tariff the user holds: one currency, one date, and its rows.

    Parameters
    ----------
    path
        The file it was read from.
    name, date, currency
        What the record calls itself, when it was captured, and the one currency it holds.
    rows
        In the order the file lists them.
    """

    path: Path
    _: KW_ONLY
    name: str
    date: str
    currency: str
    rows: tuple[PriceRow, ...]

    @property
    def source(self) -> Source:
        """The record as one source, which is what a bill cites."""
        return Source(self.name, edition=self.date, read_as=f"read from {self.path.name}")

    def row(self, key: str, quantities: Mapping[str, float] | None = None) -> PriceRow | None:
        """Return the first row pricing `key` whose every band holds, or ``None`` for none."""
        asked = quantities or {}
        return next((row for row in self.rows if row.key == key and row.holds(asked)), None)

    def charge(
        self, key: str, quantities: Mapping[str, float] | None = None, *, units: float = 1
    ) -> Decimal | None:
        """Return what `key` costs at `quantities`, or ``None`` where no row prices it.

        A ``per order`` row is flat inside its band; a ``per unit`` row is multiplied by `units`.
        """
        row = self.row(key, quantities)
        if row is None:
            return None
        return row.charge if row.basis == "per order" else row.charge * Decimal(str(units))

    def headroom(
        self, key: str, quantities: Mapping[str, float] | None = None
    ) -> tuple[Headroom, ...]:
        """Return how far each banded quantity sits from its band's top edge."""
        row = self.row(key, quantities)
        asked = quantities or {}
        if row is None:
            return ()
        return tuple(
            Headroom(band.quantity, Decimal(str(asked[band.quantity])), band) for band in row.bands
        )


def read_prices(path: str | os.PathLike[str]) -> PriceRecord:
    """Read a banded price record the user holds, as a CSV file.

    The file's name and its newest modification are what the record calls itself, unless a
    leading ``# name:`` or ``# date:`` line says otherwise.

    Raises
    ------
    ValueError
        If the file is not a price record, saying what one is.

    Examples
    --------
    >>> read_prices("prices.csv").currency  # doctest: +SKIP
    'USD'
    """
    one = Path(path)
    text = one.read_text(encoding="utf-8-sig")
    named = {
        line[1:].split(":", 1)[0].strip().lower(): line.split(":", 1)[1].strip()
        for line in text.splitlines()
        if line.startswith("#") and ":" in line
    }
    body = [line for line in text.splitlines() if not line.startswith("#")]
    try:
        rows, currency = _rows(body)
    except (KeyError, TypeError, ValueError, InvalidOperation) as error:
        raise ValueError(f"{one.name} is not a price record ({error}); expected {EXPECTED}") from (
            error
        )
    return PriceRecord(
        one,
        name=named.get("name", one.name),
        date=named.get("date", ""),
        currency=currency,
        rows=rows,
    )


def _rows(lines: Sequence[str]) -> tuple[tuple[PriceRow, ...], str]:
    reader = csv.DictReader(lines)
    missing = sorted(set(COLUMNS) - set(reader.fieldnames or ()))
    if missing:
        raise ValueError(f"missing column(s) {', '.join(missing)}")
    rows: list[PriceRow] = []
    currencies: set[str] = set()
    for line, row in enumerate(reader, 2):
        basis = (row["basis"] or "").strip()
        if basis not in ("per order", "per unit"):
            raise ValueError(f"line {line}: basis is 'per order' or 'per unit', got {basis!r}")
        if not (row["key"] or "").strip():
            raise ValueError(f"line {line}: a row needs a key")
        currencies.add((row["currency"] or "").strip().upper())
        rows.append(
            PriceRow(
                row["key"].strip(),
                Decimal((row["charge"] or "").strip()),
                basis,  # pyright: ignore[reportArgumentType]
                item=(row["item"] or "").strip(),
                bands=_bands(row["bands"] or "", line),
            )
        )
    if not rows:
        raise ValueError("no row")
    if len(currencies) > 1:
        raise ValueError(f"more than one currency: {', '.join(sorted(currencies))}")
    return tuple(rows), currencies.pop()


def _bands(cell: str, line: int) -> tuple[Band, ...]:
    bands: list[Band] = []
    for one in (part.strip() for part in cell.split(";") if part.strip()):
        quantity, _, span = one.partition(" ")
        low, dash, high = span.partition("-")
        if not quantity or not dash or not low.strip():
            raise ValueError(f"line {line}: band {one!r} is not 'quantity low-high'")
        top = high.strip()
        bands.append(Band(quantity, Decimal(low.strip()), Decimal(top) if top else None))
    return tuple(bands)


@dataclass(frozen=True, slots=True)
class Item:
    """One thing this run consumes, before anything prices it.

    Parameters
    ----------
    item
        What it is.
    quantity, unit
        What the design says the run consumes.
    key
        What a price record prices it by, a catalogue number where it has one.
    quantities
        What the record's bands are read against, such as ``{"count": 153}``. The item's own
        quantity where it is not given.
    units
        What a ``per unit`` row is multiplied by. The item's own quantity where it is not given.
    headroom
        How far a quantity sits from its band's top edge, where the caller knows the bands
        without a tariff. A vendor's bands are a fact about the catalogue, so a quantity carries
        its band whether or not anyone holds a price for it. The record's own headroom is used
        instead wherever a row prices the item.
    """

    item: str
    quantity: float
    _: KW_ONLY
    unit: str = ""
    key: str = ""
    quantities: Mapping[str, float] | None = None
    units: float | None = None
    headroom: tuple[Headroom, ...] = ()


def bill(
    items: Iterable[Item],
    record: PriceRecord | None,
    *,
    title: str = "What this run consumes",
    source_key: str = "prices",
) -> Bill:
    """Return the bill for `items`, priced by `record` where it prices them.

    With no record, every quantity still computes and every money cell is a hole. A hole here
    names no issue: a price nobody loaded is a missing input of the user's, not a defect in what
    the package knows.

    Raises
    ------
    ValueError
        If there is nothing to bill.
    """
    rows: list[BillRow] = []
    total = Decimal(0)
    priced = False
    for n, one in enumerate(items, 1):
        quantities = one.quantities if one.quantities is not None else {}
        charge = (
            None
            if record is None
            else record.charge(
                one.key, quantities, units=one.units if one.units is not None else one.quantity
            )
        )
        if charge is None:
            rows.append(
                BillRow(
                    one.item,
                    one.quantity,
                    unit=one.unit,
                    key=one.key,
                    headroom="; ".join(str(gap) for gap in one.headroom),
                    hole=Hole(
                        f"P{n}",
                        f"no row prices {one.key or one.item!r} at "
                        f"{one.quantity:g} {one.unit}".strip(),
                        "price",
                        where=f"{one.item}, money",
                        filled_by="a price record holding a row for this key",
                    ),
                )
            )
            continue
        priced = True
        total += charge
        found = record.row(one.key, quantities) if record else None
        rows.append(
            BillRow(
                one.item,
                one.quantity,
                unit=one.unit,
                key=one.key,
                charge=f"{charge:f}",
                headroom="; ".join(
                    str(gap) for gap in (record.headroom(one.key, quantities) if record else ())
                ),
                citation=Citation(source_key, found.locator) if found else None,
            )
        )
    return Bill(
        tuple(rows),
        title=title,
        currency=record.currency if record else "",
        total=f"{total:f}" if priced else "",
        record=source_key if record else "",
    )

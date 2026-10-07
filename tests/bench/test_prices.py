"""The price record the user holds, the charge it gives a quantity, and the bill."""

from decimal import Decimal
from pathlib import Path

import pytest

from liulab_mbio.bench import prices

RECORD = """# name: Vendor list, captured 2026-10-06
# date: 2026-10-06
key,item,bands,charge,basis,currency
POOL,oligo pool,count 1-100; length_nt 1-200,5000.00,per order,USD
POOL,oligo pool,count 101-200; length_nt 1-200,12505.00,per order,USD
R3733,BsaI-HFv2,,0.42,per unit,USD
"""


@pytest.fixture
def record(tmp_path: Path) -> prices.PriceRecord:
    path = tmp_path / "prices.csv"
    path.write_text(RECORD, encoding="utf-8")
    return prices.read_prices(path)


def test_a_record_names_itself_and_holds_one_currency(record: prices.PriceRecord) -> None:
    assert record.currency == "USD"
    assert record.name.startswith("Vendor list")
    assert record.source.edition == "2026-10-06"


def test_bands_are_inclusive_at_both_ends(record: prices.PriceRecord) -> None:
    assert record.charge("POOL", {"count": 100, "length_nt": 200}) == Decimal("5000.00")
    assert record.charge("POOL", {"count": 101, "length_nt": 1}) == Decimal("12505.00")


def test_per_unit_amortises_and_per_order_is_flat_inside_the_band(
    record: prices.PriceRecord,
) -> None:
    assert record.charge("R3733", units=10) == Decimal("4.20")
    assert record.charge("POOL", {"count": 150, "length_nt": 50}) == Decimal("12505.00")


def test_headroom_says_how_far_a_quantity_sits_from_the_band_edge(
    record: prices.PriceRecord,
) -> None:
    count, length = record.headroom("POOL", {"count": 153, "length_nt": 200})
    assert str(count) == "153 count, 47 below the next band"
    assert str(length) == "200 length_nt, no slack above it at all"


def test_a_file_that_is_not_a_record_is_refused_and_says_what_one_is(tmp_path: Path) -> None:
    path = tmp_path / "wrong.csv"
    path.write_text("a,b\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match=r"is not a price record.*expected"):
        prices.read_prices(path)


def test_mixing_currencies_is_refused(tmp_path: Path) -> None:
    path = tmp_path / "mixed.csv"
    path.write_text(RECORD.replace("0.42,per unit,USD", "0.42,per unit,EUR"), encoding="utf-8")
    with pytest.raises(ValueError, match="more than one currency"):
        prices.read_prices(path)


def test_a_quantity_computes_with_no_record_and_the_money_is_a_hole() -> None:
    bill = prices.bill([prices.Item("oligo pool", 153, unit="oligos", key="POOL")], None)
    (row,) = bill.rows
    assert row.quantity == 153
    assert row.charge == ""
    assert row.hole is not None
    assert row.hole.kind == "price"
    # A missing price is a missing input of the user's, not a defect in what the package knows.
    assert row.hole.issue == ""
    assert bill.total == ""


def test_a_loaded_record_prices_the_bill_and_cites_the_row(record: prices.PriceRecord) -> None:
    bill = prices.bill(
        [
            prices.Item(
                "oligo pool",
                153,
                unit="oligos",
                key="POOL",
                quantities={"count": 153, "length_nt": 200},
            ),
            prices.Item("unlisted kit", 1, key="NOPE"),
        ],
        record,
    )
    priced, unpriced = bill.rows
    assert priced.charge == "12505.00"
    assert priced.citation is not None
    assert "47 below the next band" in priced.headroom
    assert unpriced.hole is not None
    assert bill.total == "12505.00"
    assert bill.currency == "USD"


def test_a_tier_with_no_top_leaves_its_high_end_empty(tmp_path: Path) -> None:
    path = tmp_path / "open.csv"
    path.write_text(RECORD.replace("count 101-200", "count 101-"), encoding="utf-8")

    record = prices.read_prices(path)

    assert record.charge("POOL", {"count": 5000, "length_nt": 200}) == Decimal("12505.00")
    count, _length = record.headroom("POOL", {"count": 5000, "length_nt": 200})
    assert str(count) == "5000 count, no band above it"


def test_a_band_that_names_no_span_is_refused(tmp_path: Path) -> None:
    path = tmp_path / "bad.csv"
    path.write_text(RECORD.replace("count 1-100", "count 100"), encoding="utf-8")

    with pytest.raises(ValueError, match="is not 'quantity low-high'"):
        prices.read_prices(path)

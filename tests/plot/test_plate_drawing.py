"""A plate drawn the way a record is: laid out as geometry, written as a page, a PNG or a PDF."""

from pathlib import Path

import pytest

from liulab_mbio.plot import draw_plate
from liulab_mbio.plot.plate import EMPTY
from liulab_mbio.plot.svg import Group, Rect
from liulab_mbio.protocol.model import FORMATS, Plate


def drawing(wells: int, **rest: object):
    one = Plate("plate", wells, **rest)  # pyright: ignore[reportArgumentType]
    return draw_plate(one.name, one.rows, one.columns, one.row_labels, seating=one.seating)


@pytest.mark.parametrize("wells", sorted(FORMATS))
def test_every_format_draws_one_shape_per_well_inside_the_page(wells: int) -> None:
    laid = drawing(wells).layout
    circles = [shape for shape in laid.shapes if isinstance(shape, Group)]
    assert len(circles) == wells
    for group in circles:
        (rect,) = group.shapes
        assert isinstance(rect, Rect)
        assert laid.extent.x <= rect.box.x
        assert rect.box.x + rect.box.width <= laid.extent.x + laid.extent.width


def test_a_well_holding_something_is_filled_and_an_empty_one_is_not() -> None:
    laid = drawing(96, seating={"A1": "UMI-1"}).layout
    groups = [shape for shape in laid.shapes if isinstance(shape, Group)]
    first, second = groups[0].shapes[0], groups[1].shapes[0]
    assert isinstance(first, Rect)
    assert isinstance(second, Rect)
    assert first.fill != EMPTY
    assert second.fill == EMPTY
    assert laid.legend == (("UMI-1", "#e48b8b"),)


def test_each_distinct_content_gets_a_fill_of_its_own() -> None:
    laid = drawing(96, seating={"A1": "a", "A2": "b", "A3": "a"}).layout
    assert len({fill for _, fill in laid.legend}) == 2


def test_a_drawing_is_laid_out_once_and_kept() -> None:
    one = drawing(1536)
    assert one.layout is one.layout


@pytest.mark.parametrize("suffix", [".html", ".png", ".pdf"])
def test_a_plate_writes_in_every_format_the_destination_names(tmp_path: Path, suffix: str) -> None:
    out = drawing(96, seating={"A1": "water"}).write(tmp_path / f"plate{suffix}")
    assert out.stat().st_size > 0


def test_a_suffix_nothing_writes_is_refused(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="the suffix must be"):
        drawing(96).write(tmp_path / "plate.jpg")


def test_a_row_label_short_of_its_rows_is_refused() -> None:
    with pytest.raises(ValueError, match="1 labels for 2 rows"):
        _ = draw_plate("plate", 2, 3, "A").layout

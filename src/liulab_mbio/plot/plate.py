"""A plate laid out as a grid of wells, in points, before anything draws it.

This module knows nothing of a protocol: it takes a row and column count, the row labels the
caller names them by, and what each well holds. A well's own name is its row label and its
1-based column, so the drawing numbers wells the way the bench does.

Laid out as geometry, like every other drawing here — `docs/adr/0006-drawings-as-geometry.md`.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from liulab_mbio.plot import svg
from liulab_mbio.plot.fonts import BOLD, SANS
from liulab_mbio.plot.labels import Box

#: The widest a plate is laid out, in points, before its wells shrink to fit. A 1536-well plate
#: is 48 wells across, so the pitch is what gives way, never the page.
WIDTH = 540.0

#: The widest one well's pitch is, in points, so a 12-well plate does not draw a wall of discs.
MAX_PITCH = 26.0

#: The margin round the grid, and the gap between a well's circle and its neighbour's.
MARGIN = 10.0
GAP = 0.18

#: An empty well, and the lines everything is drawn in.
EMPTY = "#f4f4f4"
INK = "#333333"
RULE = "#bbbbbb"

#: The fill each kind takes, in order, and so the most kinds a plate colours at all: the
#: petroff10 sequence, whose colours stay apart under every colour-vision deficiency. Matthew A.
#: Petroff, *Accessible Color Sequences for Data Visualization*, arXiv:2107.02270, as matplotlib
#: ships it (`_petroff10_data` in `_cm.py`).
PALETTE = (
    "#3f90da",
    "#ffa90e",
    "#bd1f01",
    "#94a4a2",
    "#832db6",
    "#a96b59",
    "#e76300",
    "#b9ac70",
    "#717581",
    "#92dadd",
)


@dataclass(frozen=True, slots=True)
class PlateMap:
    """Where every well went, and how big the drawing is.

    Parameters
    ----------
    shapes
        In the order they are drawn.
    extent
        The view box, in points.
    legend
        What each colour stands for, in the order the legend lists them. Empty where the plate
        holds more kinds than `PALETTE` has colours.
    kinds
        How many distinct contents the plate holds, coloured or not.
    """

    shapes: tuple[svg.Shape, ...]
    extent: Box
    legend: tuple[tuple[str, str], ...]
    kinds: int


def layout(
    rows: int,
    columns: int,
    row_labels: Sequence[str],
    *,
    title: str = "",
    seating: Mapping[str, str] | None = None,
) -> PlateMap:
    """Lay a plate out: its wells, its row and column labels, and a legend of what it holds.

    Parameters
    ----------
    rows, columns
        The array's shape.
    row_labels
        One per row, top to bottom, no two of them spelling one well's name.
    title
        Drawn above the grid.
    seating
        Well name to what sits there, named as this array names its wells: a row label and a
        1-based column. A well named nothing is drawn empty, and every distinct content gets a
        fill of its own while there are no more kinds than `PALETTE` has colours. Past that no
        fill tells two kinds apart and no legend is readable, so every well is drawn empty and
        the legend is dropped.

    Raises
    ------
    ValueError
        If the array has no wells, there is not one label per row, two wells would answer to
        one name, or a well is seated where the array has none.

    Examples
    --------
    >>> one = layout(2, 3, ("A", "B"), seating={"A1": "water"})
    >>> one.legend
    (('water', '#3f90da'),)
    """
    if rows < 1 or columns < 1:
        raise ValueError(f"a plate has at least one row and column, not {rows} by {columns}")
    if len(row_labels) != rows:
        raise ValueError(f"{len(row_labels)} labels for {rows} rows")
    held = seating or {}
    # Each well is drawn under the name this grid gives it, so a seating is checked against the
    # same names: one not in it would be dropped from the drawing and still colour the legend.
    grid = [[f"{label}{column + 1}" for column in range(columns)] for label in row_labels]
    named: set[str] = set()
    for name in (one for row in grid for one in row):
        if name in named:
            raise ValueError(f"two wells of a {rows} by {columns} plate are both named {name!r}")
        named.add(name)
    for well in held:
        if well not in named:
            raise ValueError(f"no well {well!r} on a {rows} by {columns} plate")
    kinds = tuple(dict.fromkeys(held.values()))
    fills = dict(zip(kinds, PALETTE, strict=False)) if len(kinds) <= len(PALETTE) else {}

    pitch = min(MAX_PITCH, (WIDTH - 2 * MARGIN) / (columns + 1))
    size = pitch * (1 - GAP)
    label_size = max(3.5, min(9.0, pitch * 0.5))
    left, top = MARGIN + pitch, MARGIN + pitch + (0 if not title else pitch)
    shapes: list[svg.Shape] = []
    if title:
        shapes.append(svg.Text(MARGIN, MARGIN + pitch * 0.7, title, BOLD, min(13.0, pitch * 0.8)))
    for column in range(columns):
        text = str(column + 1)
        shapes.append(
            svg.Text(
                left + column * pitch + (pitch - SANS.width(text, label_size)) / 2,
                top - pitch * 0.35,
                text,
                SANS,
                label_size,
            )
        )
    for row in range(rows):
        shapes.append(
            svg.Text(
                left - pitch * 0.35 - SANS.width(row_labels[row], label_size),
                top + row * pitch + pitch * 0.62,
                row_labels[row],
                SANS,
                label_size,
            )
        )
        for column in range(columns):
            name = grid[row][column]
            holds = held.get(name, "")
            box = Box(
                left + column * pitch + (pitch - size) / 2,
                top + row * pitch + (pitch - size) / 2,
                size,
                size,
            )
            shapes.append(
                svg.Group(
                    (svg.Rect(box, fills.get(holds, EMPTY), RULE, 0.4, size / 2),),
                    classes=("well",),
                    data={"well": name, "holds": holds} if holds else {"well": name},
                )
            )
    extent = Box(0.0, 0.0, left + columns * pitch + MARGIN, top + rows * pitch + MARGIN)
    return PlateMap(tuple(shapes), extent, tuple(fills.items()), len(kinds))

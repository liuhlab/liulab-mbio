"""Keep labels apart: the no-overlap rule, on bare boxes and the points their leaders start from.

Coordinates are in points, x to the right and y down, as an SVG draws them.

`columns` lays labels out round a circle as SnapGene's circular map does. Labels go in two
columns, right and left of the circle, in the order their sites run clockwise from the top. Each
column is spread to the nearest positions, in least squares, at which no two labels overlap.
Each label's inner edge touches an ellipse round the circle, and a straight leader runs from its
site to that edge. The ellipse grows until every column fits beside it and no leader enters the
circle.

No label can overlap the circle, since every label lies outside the ellipse and the circle inside
it. A leader running to the point where its label touches the ellipse stays inside the ellipse,
where no label is, so it crosses none; a leader meets its label nearer halfway up only where the
part of it outside the ellipse crosses no other label either.
"""

import itertools
import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import NamedTuple

#: How much the ellipse grows each time a column does not fit or a leader enters the circle.
GROWTH = 1.04

# How many heights from halfway up a label towards its touching point a leader tries first.
_STEPS = 4


class Point(NamedTuple):
    """A point, in points."""

    x: float
    y: float


@dataclass(frozen=True, slots=True)
class Box:
    """A label's box: its top-left corner, width and height."""

    x: float
    y: float
    width: float
    height: float


@dataclass(frozen=True, slots=True)
class Anchored:
    """A label to place: how big its box is, and the point its leader starts from."""

    width: float
    height: float
    anchor: Point


@dataclass(frozen=True, slots=True)
class Placed:
    """Where a label went: its box, and its leader from the anchor to the box."""

    box: Box
    leader: tuple[Point, Point]


def columns(
    labels: Sequence[Anchored], *, radius: float, gap: float, spacing: float, margin: float
) -> tuple[Placed, ...]:
    """Place labels in two columns round a circle of `radius` centred on (0, 0).

    Parameters
    ----------
    labels
        Each anchor lies outside the circle and less than `gap` beyond it.
    radius
        No box overlaps the circle, and no leader enters it.
    gap
        How far beyond the circle the ellipse starts.
    spacing
        The least distance between two boxes in one column.
    margin
        The least distance of each box from the vertical line through the centre, so the two
        columns never meet. Less than ``radius + gap``.

    Returns
    -------
    tuple[Placed, ...]
        One placement for each label, in the order given.
    """
    angles = [_clockwise(label.anchor) for label in labels]
    right = sorted((i for i, angle in enumerate(angles) if angle < math.pi), key=angles.__getitem__)
    left = sorted(
        (i for i, angle in enumerate(angles) if angle >= math.pi), key=lambda i: -angles[i]
    )
    start = radius + gap
    middles: dict[int, float] = {}
    for column in (right, left):
        targets = [-start * math.cos(angles[i]) for i in column]
        spread = _spread(targets, [labels[i].height for i in column], spacing)
        middles.update(zip(column, spread, strict=True))
    nearest = {i: _nearest(labels[i], middle) for i, middle in middles.items()}
    a = b = start
    while True:
        # A box level with the ellipse's top or bottom would reach the other column.
        if any(near**2 >= b * b * (1 - (margin / a) ** 2) for near in nearest.values()):
            b *= GROWTH
            continue
        placed: dict[int, Placed] = {}
        for column, side in ((right, 1.0), (left, -1.0)):
            boxes = {}
            for i in column:
                inner = side * a * math.sqrt(1 - (nearest[i] / b) ** 2)
                x = inner if side > 0 else inner - labels[i].width
                boxes[i] = Box(
                    x, middles[i] - labels[i].height / 2, labels[i].width, labels[i].height
                )
            for i in column:
                others = [box for j, box in boxes.items() if j != i]
                port = _port(labels[i].anchor, boxes[i], others, nearest[i], side, radius, (a, b))
                placed[i] = Placed(boxes[i], (labels[i].anchor, port))
        if not any(_enters(one.leader, radius) for one in placed.values()):
            return tuple(placed[i] for i in range(len(labels)))
        a *= GROWTH


def _clockwise(point: Point) -> float:
    """Return the angle of `point` about the centre, clockwise from the top, in [0, 2π)."""
    return math.atan2(point.x, -point.y) % math.tau


def _spread(targets: Sequence[float], heights: Sequence[float], spacing: float) -> list[float]:
    """Return the middles nearest `targets`, in least squares, that keep boxes `spacing` apart.

    Subtracting each box's least offset from the first turns this into isotonic regression, which
    pooling adjacent violators solves exactly.
    """
    if not targets:
        return []
    offsets = [0.0]
    for above, below in itertools.pairwise(heights):
        offsets.append(offsets[-1] + (above + below) / 2 + spacing)
    pools: list[list[float]] = []
    for target, offset in zip(targets, offsets, strict=True):
        pools.append([target - offset, 1.0])
        while len(pools) > 1 and pools[-2][0] / pools[-2][1] > pools[-1][0] / pools[-1][1]:
            total, count = pools.pop()
            pools[-1][0] += total
            pools[-1][1] += count
    fitted = [total / count for total, count in pools for _ in range(int(count))]
    return [value + offset for value, offset in zip(fitted, offsets, strict=True)]


def _nearest(label: Anchored, middle: float) -> float:
    """Return the height, within a label whose middle is at `middle`, nearest the centre's."""
    return min(max(0.0, middle - label.height / 2), middle + label.height / 2)


def _port(
    anchor: Point,
    box: Box,
    others: Sequence[Box],
    near: float,
    side: float,
    radius: float,
    ellipse: tuple[float, float],
) -> Point:
    """Return where a leader meets its label's inner edge: as near halfway up as keeps the rule.

    The point touching the ellipse always keeps it. One nearer the middle keeps it when the leader
    stays out of the circle and the part of it outside the ellipse crosses no other box.
    """
    edge = box.x if side > 0 else box.x + box.width
    middle = box.y + box.height / 2
    for step in range(_STEPS):
        port = Point(edge, middle + (near - middle) * step / _STEPS)
        low, high = sorted((_leaves(anchor, port, *ellipse), port.y))
        if not _enters((anchor, port), radius) and not any(
            _crosses((anchor, port), other)
            for other in others
            if other.y < high and low < other.y + other.height
        ):
            return port
    return Point(edge, near)


def _leaves(inside: Point, outside: Point, a: float, b: float) -> float:
    """Return the height at which a segment from a point inside the ellipse leaves it."""
    dx, dy = outside.x - inside.x, outside.y - inside.y
    qa = (dx / a) ** 2 + (dy / b) ** 2
    qb = 2 * (inside.x * dx / a**2 + inside.y * dy / b**2)
    qc = (inside.x / a) ** 2 + (inside.y / b) ** 2 - 1
    if qa == 0:
        return inside.y
    t = (-qb + math.sqrt(max(0.0, qb * qb - 4 * qa * qc))) / (2 * qa)
    return inside.y + min(t, 1.0) * dy


def _crosses(leader: tuple[Point, Point], box: Box) -> bool:
    """Return whether the segment passes through the inside of the box."""
    (x0, y0), (x1, y1) = leader
    low, high = 0.0, 1.0
    for p, q in (
        (x0 - x1, x0 - box.x),
        (x1 - x0, box.x + box.width - x0),
        (y0 - y1, y0 - box.y),
        (y1 - y0, box.y + box.height - y0),
    ):
        if p == 0:
            if q <= 0:
                return False
        elif p < 0:
            low = max(low, q / p)
        else:
            high = min(high, q / p)
    return low < high


def _enters(leader: tuple[Point, Point], radius: float) -> bool:
    """Return whether the segment comes nearer the centre than `radius`."""
    (x0, y0), (x1, y1) = leader
    dx, dy = x1 - x0, y1 - y0
    length = dx * dx + dy * dy
    t = 0.0 if length == 0 else min(1.0, max(0.0, -(x0 * dx + y0 * dy) / length))
    return math.hypot(x0 + t * dx, y0 + t * dy) < radius

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

`staircase` lays labels out in rows above a line, as SnapGene's linear map and sequence view do:
panorama labelling. Each box hangs from a vertical leader to its anchor, and sits above every
label whose anchor its own box reaches, so close neighbours rise in a staircase. No two boxes in a
row meet, and no leader crosses another label, since any box lying across a leader reaches its
anchor and so sits higher than the leader goes. `pack` puts bars in rows, each where it first fits
in order of start, which uses the fewest rows.

Neither caps how far labels reach. `hide_in_columns` and `hide_in_staircase` say which labels to
leave out so the rest stay within a cap: a crowd of labels that push one another past it hides its
labels in a given order until it fits, and a label in no such crowd never hides.
"""

import itertools
import math
from bisect import bisect_left, bisect_right
from collections.abc import Callable, Iterable, Sequence
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
    right, left = _columns(labels)
    start = radius + gap
    middles: dict[int, float] = {}
    for column in (right, left):
        spread, _ = _spread(column, labels, start, spacing)
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


def hide_in_columns(
    labels: Sequence[Anchored],
    order: Sequence[int],
    *,
    radius: float,
    gap: float,
    spacing: float,
    reach: float,
) -> tuple[int, ...]:
    """Return the labels to leave out so `columns` places the rest within `reach` of the centre.

    A crowd is a run of labels a column packs against one another. While a crowd reaches further
    than `reach` above or below the centre, the first of its labels in `order` hides, and the
    column is spread again. A crowd within `reach` hides nothing.

    Parameters
    ----------
    labels
        As `columns` takes them.
    order
        Every label's index once, the first to hide first.
    radius, gap, spacing
        As `columns` takes them.
    reach
        How far above or below the centre a box may reach.

    Returns
    -------
    tuple[int, ...]
        The labels hidden, in the order they hid.
    """
    rank = {index: position for position, index in enumerate(order)}
    hidden: list[int] = []
    for column in _columns(labels):
        shown = column
        while True:
            middles, crowds = _spread(shown, labels, radius + gap, spacing)
            at = dict(zip(shown, middles, strict=True))
            over = [
                i
                for crowd in crowds
                if any(abs(at[j]) + labels[j].height / 2 > reach for j in crowd)
                for i in crowd
            ]
            if not over:
                break
            first = min(over, key=rank.__getitem__)
            shown = [i for i in shown if i != first]
            hidden.append(first)
    return tuple(sorted(hidden, key=rank.__getitem__))


def staircase(labels: Sequence[Anchored], *, base: float, spacing: float) -> tuple[Placed, ...]:
    """Place labels in rows above a line, each box hanging from a vertical leader.

    Labels are placed from the right. A box hangs right from its anchor, just above every box
    already placed that it meets, so it sits above each label whose anchor it reaches, and labels
    at one anchor rise in the order given. It hangs left instead where that puts it lower and its
    box reaches no other anchor, so a label beside a crowd is not lifted over it.

    Parameters
    ----------
    labels
        Each anchor lies at or below `base`.
    base
        Where the lowest boxes end.
    spacing
        The least distance between two boxes, side by side or one above the other.

    Returns
    -------
    tuple[Placed, ...]
        One placement for each label, in the order given, its leader running up from the anchor
        to the bottom of its box.
    """
    xs = [label.anchor.x for label in labels]
    ordered = sorted(xs)
    placed: dict[int, tuple[float, float]] = {}
    # Each box placed: its left edge, its right edge and `spacing` beyond, and the lowest a box
    # above it may end.
    boxes: list[tuple[float, float, float]] = []

    def bottom(left: float, width: float) -> float:
        right = left + width + spacing
        return min((low for start, end, low in boxes if start < right and left < end), default=base)

    for i in sorted(range(len(labels)), key=lambda i: (-xs[i], i)):
        width, height = labels[i].width, labels[i].height
        right = bottom(xs[i], width)
        # Whether no other anchor lies under the box hung left, or within `spacing` before it.
        free = bisect_right(ordered, xs[i]) - bisect_left(ordered, xs[i] - width - spacing) <= 1
        left = bottom(xs[i] - width, width) if free else right
        placed[i] = (left, xs[i] - width) if left > right else (right, xs[i])
        low, start = placed[i]
        boxes.append((start, start + width + spacing, low - height - spacing))
    result = []
    for (low, x), label in zip((placed[i] for i in range(len(labels))), labels, strict=True):
        box = Box(x, low - label.height, label.width, label.height)
        result.append(Placed(box, (label.anchor, Point(label.anchor.x, low))))
    return tuple(result)


def pack(bars: Sequence[tuple[float, float]], *, spacing: float) -> tuple[int, ...]:
    """Return the row each bar goes in, from 0: the first it fits, taking bars by earliest start.

    Each bar is its start and end. A bar fits a row when it starts at least `spacing` past the end
    of every bar already there, and no other rule uses fewer rows.
    """
    ends: list[float] = []
    rows = [0] * len(bars)
    for i in sorted(range(len(bars)), key=lambda i: (bars[i][0], i)):
        start, end = bars[i]
        row = next((row for row, last in enumerate(ends) if last + spacing <= start), len(ends))
        if row == len(ends):
            ends.append(end)
        ends[row] = end
        rows[i] = row
    return tuple(rows)


def hide_in_staircase(
    labels: Sequence[Anchored],
    order: Sequence[int],
    *,
    base: float,
    spacing: float,
    ceiling: float,
) -> tuple[int, ...]:
    """Return the labels to leave out so `staircase` places the rest no higher than `ceiling`.

    Two labels can lift each other only when their anchors lie closer than their widths and
    `spacing` added up, and a crowd is a run of labels linked so. While a crowd rises above
    `ceiling`, its labels hide in `order`. A crowd below `ceiling` hides nothing.

    Parameters
    ----------
    labels
        As `staircase` takes them.
    order
        Every label's index once, the first to hide first.
    base, spacing
        As `staircase` takes them.
    ceiling
        The height no box's top may rise above.

    Returns
    -------
    tuple[int, ...]
        The labels hidden, in the order they hid.
    """
    rank = {index: position for position, index in enumerate(order)}
    hidden: list[int] = []

    def fits(indices: Sequence[int]) -> bool:
        placed = staircase([labels[i] for i in indices], base=base, spacing=spacing)
        return all(one.box.y >= ceiling for one in placed)

    pending = _crowds(labels, range(len(labels)), spacing)
    while pending:
        crowd = pending.pop()
        if fits(crowd):
            continue
        ranked = sorted(crowd, key=rank.__getitem__)

        def rest(count: int, crowd: list[int] = crowd, ranked: list[int] = ranked) -> list[int]:
            gone = set(ranked[:count])
            return [i for i in crowd if i not in gone]

        # Hiding more of a crowd almost never lifts the rest, so halving finds how many to hide;
        # the count it finds always fits.
        least = _least(0, len(ranked), lambda count: fits(rest(count)))
        if len(_crowds(labels, rest(least - 1), spacing)) <= 1:
            hidden.extend(ranked[:least])
            continue
        # The crowd breaks up before it fits: hide what it takes to break it, then each part alone.
        breaks = _least(0, least - 1, lambda count: len(_crowds(labels, rest(count), spacing)) > 1)
        hidden.extend(ranked[:breaks])
        pending.extend(_crowds(labels, rest(breaks), spacing))
    return tuple(sorted(hidden, key=rank.__getitem__))


def _columns(labels: Sequence[Anchored]) -> tuple[list[int], list[int]]:
    """Return the right column's labels and the left's, each from the top down."""
    angles = [_clockwise(label.anchor) for label in labels]
    right = sorted((i for i, angle in enumerate(angles) if angle < math.pi), key=angles.__getitem__)
    left = sorted(
        (i for i, angle in enumerate(angles) if angle >= math.pi), key=lambda i: -angles[i]
    )
    return right, left


def _clockwise(point: Point) -> float:
    """Return the angle of `point` about the centre, clockwise from the top, in [0, 2π)."""
    return math.atan2(point.x, -point.y) % math.tau


def _spread(
    column: Sequence[int], labels: Sequence[Anchored], start: float, spacing: float
) -> tuple[list[float], list[list[int]]]:
    """Return the middles nearest their anchors' heights, in least squares, keeping boxes apart.

    Each box's target is its anchor's angle carried out to `start`. Subtracting each box's least
    offset from the first turns this into isotonic regression, which pooling adjacent violators
    solves exactly. The pools come back too, as the crowds of labels packed against one another.
    """
    if not column:
        return [], []
    targets = [-start * math.cos(_clockwise(labels[i].anchor)) for i in column]
    offsets = [0.0]
    for above, below in itertools.pairwise(column):
        offsets.append(offsets[-1] + (labels[above].height + labels[below].height) / 2 + spacing)
    pools: list[list[float]] = []
    for target, offset in zip(targets, offsets, strict=True):
        pools.append([target - offset, 1.0])
        while len(pools) > 1 and pools[-2][0] / pools[-2][1] > pools[-1][0] / pools[-1][1]:
            total, count = pools.pop()
            pools[-1][0] += total
            pools[-1][1] += count
    fitted = [total / count for total, count in pools for _ in range(int(count))]
    crowds, first = [], 0
    for _, count in pools:
        crowds.append(list(column[first : first + int(count)]))
        first += int(count)
    return [value + offset for value, offset in zip(fitted, offsets, strict=True)], crowds


def _crowds(labels: Sequence[Anchored], indices: Iterable[int], spacing: float) -> list[list[int]]:
    """Return `indices` split into crowds along the line, each crowd's indices in increasing order.

    Each box hangs within its width of its anchor, so two labels whose anchors lie further apart
    than their widths and `spacing` added up never meet.
    """
    reaches = sorted((labels[i].anchor.x - labels[i].width - spacing / 2, i) for i in indices)
    crowds: list[list[int]] = []
    end = -math.inf
    for left, i in reaches:
        if left >= end:
            crowds.append([])
        crowds[-1].append(i)
        end = max(end, labels[i].anchor.x + labels[i].width + spacing / 2)
    return [sorted(crowd) for crowd in crowds]


def _least(low: int, high: int, holds: Callable[[int], bool]) -> int:
    """Return the least count above `low`, and at most `high`, for which `holds` does.

    `holds` fails at `low`, holds at `high`, and holds for every count above one it holds for.
    """
    while high - low > 1:
        middle = (low + high) // 2
        if holds(middle):
            high = middle
        else:
            low = middle
    return high


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

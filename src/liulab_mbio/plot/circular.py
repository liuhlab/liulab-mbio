"""The circular map: items in, shapes out, in points round a circle centred on (0, 0).

The backbone is two concentric lines, with a bp scale inside them and the record's name and length
in the centre. Each feature is drawn segment by segment as an arrow along its strand, or a box when
it has none, and a segment across the origin is one piece. Features that overlap go on rings
further in, the longest outermost. Each feature's name is boxed outside, in the columns
`labels.columns` lays out, so no label overlaps another label or the drawing.
"""

import itertools
import math
from collections.abc import Sequence
from dataclasses import dataclass

from liulab_mbio.plot import labels
from liulab_mbio.plot.fonts import BOLD, SANS, Font
from liulab_mbio.plot.labels import Box, Point
from liulab_mbio.plot.layers import Item, Span, outline_color, text_color, unwrapped
from liulab_mbio.plot.svg import Circle, Group, Line, Path, Rect, Shape, Text, number
from liulab_mbio.sequence import Strand

#: The backbone's radius, in points. Everything else on the circle is laid out from it.
RADIUS = 220.0

#: The size of a label's text, in points.
LABEL_SIZE = 11.5
#: The size of the record's name in the centre, in points.
NAME_SIZE = 13.5
#: The size of the record's length and of the scale's numbers, in points.
SMALL_SIZE = 10.5

_INK = "#252525"
_LEADER = "#7f7f7f"
_CONNECTOR = "#666666"

# The backbone's two lines, and the ticks inside them.
_BACKBONE = 2.6
_STRANDS_APART = 3.8
_TICK = 6.0
_ORIGIN_TICK = 12.0
# How far a scale number keeps from the origin's tick.
_CLEAR = 10.0

# A ring of features: the band's thickness, how far an arrowhead reaches beyond the band and
# along it, the distance between neighbouring rings, and the room left for the centre's text.
_BAND = 15.0
_OVERHANG = 4.0
_HEAD = 12.0
_RING_STEP = 25.0
_CENTRE = 0.4 * RADIUS
_RING_PADDING = 3.0

# The labels: padding inside a box, where leaders start and the ellipse begins beyond the
# drawing, the space between boxes in a column and between the columns, and the canvas's margin.
_PADDING = (3.0, 1.0)
_ANCHOR = 4.0
_GAP = 14.0
_SPACING = 2.0
_MARGIN = 3.0
_EDGE = 16.0


@dataclass(frozen=True, slots=True)
class Arrow:
    """One segment of a feature as drawn: a band round the centre, with a head at either end.

    Parameters
    ----------
    item, span
        What it draws.
    start, end
        Angles in radians, clockwise from the top; `end` passes a full turn when the segment
        crosses the origin.
    radius, width
        The band's middle, and its thickness.
    head_start, head_end
        Whether an arrowhead's tip ends the arrow at `start`, or at `end`.
    """

    item: Item
    span: Span
    start: float
    end: float
    radius: float
    width: float
    head_start: bool
    head_end: bool


@dataclass(frozen=True, slots=True)
class Label:
    """An item's name boxed outside the circle, and its leader from the circle to the box."""

    item: Item
    box: Box
    leader: tuple[Point, Point]


@dataclass(frozen=True, slots=True)
class CircularMap:
    """A circular map laid out: what went where, and the shapes that draw it.

    Parameters
    ----------
    extent
        The canvas.
    radius
        How far from the centre the drawing reaches. Only labels and their leaders lie further.
    arrows
        Every segment of every item, in the items' order.
    labels
        Every boxed name, in the items' order.
    shapes
        Everything, in the order it is drawn.
    """

    extent: Box
    radius: float
    arrows: tuple[Arrow, ...]
    labels: tuple[Label, ...]
    shapes: tuple[Shape, ...]


def layout(items: Sequence[Item], *, name: str, length: int) -> CircularMap:
    """Lay out `items` as a circular map of a record called `name`, `length` bases long."""
    drawing = RADIUS + _BACKBONE / 2
    arrows = _arrows(items, length)
    boxed = _labels(items, length, drawing)
    shapes = (
        *_backbone(length),
        *(_feature(item, [arrow for arrow in arrows if arrow.item is item]) for item in items),
        *(_label(label) for label in boxed),
        *_centre(name, length),
    )
    return CircularMap(_extent(drawing + _ANCHOR, boxed), drawing, arrows, boxed, shapes)


def _angle(position: float, length: int) -> float:
    return math.tau * position / length


def _point(radius: float, angle: float) -> Point:
    return Point(radius * math.sin(angle), -radius * math.cos(angle))


def _height(font: Font, size: float) -> float:
    return (font.ascender - font.descender) / font.units_per_em * size


def _baseline(font: Font, size: float, middle: float) -> float:
    """Return the baseline that centres a line of text on the height `middle`."""
    return middle + (font.ascender + font.descender) / 2 / font.units_per_em * size


def _hull(item: Item, length: int) -> tuple[int, int]:
    counted = unwrapped(item.spans, length)
    return counted[0][0], min(counted[-1][1], counted[0][0] + length)


def _rings(items: Sequence[Item], length: int) -> list[int]:
    """Return each item's ring, 0 outermost: the first it overlaps nothing on, longest first."""
    padding = length * _RING_PADDING / (math.tau * RADIUS)
    hulls = [_hull(item, length) for item in items]
    order = sorted(range(len(items)), key=lambda i: (hulls[i][0] - hulls[i][1], hulls[i][0]))
    rings: list[int] = [0] * len(items)
    taken: list[list[tuple[float, float]]] = []
    for i in order:
        start, end = hulls[i][0] - padding, hulls[i][1] + padding
        ring = next(
            (
                index
                for index, spans in enumerate(taken)
                if not any(
                    max(start, other_start + shift) < min(end, other_end + shift)
                    for other_start, other_end in spans
                    for shift in (-length, 0, length)
                )
            ),
            len(taken),
        )
        if ring == len(taken):
            taken.append([])
        taken[ring].append((start, end))
        rings[i] = ring
    return rings


def _arrows(items: Sequence[Item], length: int) -> tuple[Arrow, ...]:
    rings = _rings(items, length)
    count = max(rings, default=0) + 1
    outer = (
        RADIUS
        - _STRANDS_APART
        - _BACKBONE / 2
        - _TICK
        - 2 * _PADDING[1]
        - _height(SANS, SMALL_SIZE)
        - _OVERHANG
        - _BAND / 2
    )
    step = min(_RING_STEP, (outer - _CENTRE) / (count - 1)) if count > 1 else _RING_STEP
    width = _BAND * step / _RING_STEP
    arrows: list[Arrow] = []
    for item, ring in zip(items, rings, strict=True):
        last = len(item.spans) - 1
        for index, (span, (start, end)) in enumerate(
            zip(item.spans, unwrapped(item.spans, length), strict=True)
        ):
            arrows.append(
                Arrow(
                    item,
                    span,
                    _angle(start, length),
                    _angle(end, length),
                    outer - ring * step,
                    width,
                    head_start=index == 0 and item.strand in (Strand.REVERSE, Strand.BOTH),
                    head_end=index == last and item.strand in (Strand.FORWARD, Strand.BOTH),
                )
            )
    return tuple(arrows)


def _labels(items: Sequence[Item], length: int, drawing: float) -> tuple[Label, ...]:
    named = [item for item in items if SANS.drawn(item.label)]
    anchored = []
    for item in named:
        start, end = _hull(item, length)
        anchored.append(
            labels.Anchored(
                SANS.width(item.label, LABEL_SIZE) + 2 * _PADDING[0],
                _height(SANS, LABEL_SIZE) + 2 * _PADDING[1],
                _point(drawing + _ANCHOR, _angle((start + end) / 2, length)),
            )
        )
    placed = labels.columns(anchored, radius=drawing, gap=_GAP, spacing=_SPACING, margin=_MARGIN)
    return tuple(Label(item, one.box, one.leader) for item, one in zip(named, placed, strict=True))


def _extent(reach: float, boxed: Sequence[Label]) -> Box:
    left = min([-reach, *(label.box.x for label in boxed)]) - _EDGE
    right = max([reach, *(label.box.x + label.box.width for label in boxed)]) + _EDGE
    top = min([-reach, *(label.box.y for label in boxed)]) - _EDGE
    bottom = max([reach, *(label.box.y + label.box.height for label in boxed)]) + _EDGE
    return Box(left, top, right - left, bottom - top)


def _backbone(length: int) -> list[Shape]:
    inside = RADIUS - _STRANDS_APART - _BACKBONE / 2
    shapes: list[Shape] = [
        Circle(0, 0, RADIUS, _INK, _BACKBONE),
        Circle(0, 0, RADIUS - _STRANDS_APART, _INK, _BACKBONE),
        Line(0, -inside, 0, -inside + _ORIGIN_TICK, _INK, 2.4),
    ]
    middle = inside - _TICK - _PADDING[1] - _height(SANS, SMALL_SIZE) / 2
    scale: list[Shape] = []
    step = _step(length)
    for position in range(step, length, step):
        angle = _angle(position, length)
        scale.append(Line(*_point(inside, angle), *_point(inside - _TICK, angle), _INK, 1.2))
        text = str(position)
        width = SANS.width(text, SMALL_SIZE)
        if min(angle, math.tau - angle) * middle < width / 2 + _CLEAR:
            continue
        upright = math.cos(angle) >= 0
        centre = -middle if upright else middle
        number_text = Text(
            -width / 2, _baseline(SANS, SMALL_SIZE, centre), text, SANS, SMALL_SIZE, _INK
        )
        turn = math.degrees(angle) - (0 if upright else 180)
        scale.append(Group((number_text,), rotate=turn))
    shapes.append(Group(tuple(scale), classes=("scale",)))
    return shapes


def _step(length: int) -> int:
    """Return the least interval of 1, 2 or 5 times a power of ten marking at most nine numbers."""
    power = 1
    while True:
        for factor in (1, 2, 5):
            if factor * power * 10 >= length:
                return factor * power
        power *= 10


def _feature(item: Item, arrows: Sequence[Arrow]) -> Group:
    """Return an item's arrows, each segment joined to the next by a thin arc across any gap."""
    shapes: list[Shape] = []
    for before, after in itertools.pairwise(arrows):
        if after.start > before.end:
            joint = _move(before.radius, before.end) + _arc(before.radius, before.end, after.start)
            shapes.append(Path(joint, "none", _CONNECTOR, 1.5))
    shapes.extend(
        Path(_outline(arrow), arrow.span.color, outline_color(arrow.span.color), 0.8)
        for arrow in arrows
    )
    return Group(tuple(shapes), classes=(item.kind,), data={"kind": item.kind, **item.hover})


def _label(label: Label) -> Group:
    box, item = label.box, label.item
    (x1, y1), (x2, y2) = label.leader
    text = Text(
        box.x + _PADDING[0],
        _baseline(SANS, LABEL_SIZE, box.y + box.height / 2),
        item.label,
        SANS,
        LABEL_SIZE,
        text_color(item.color),
    )
    shapes = (
        Line(x1, y1, x2, y2, _LEADER, 0.8),
        Rect(box, item.color, outline_color(item.color), 0.8, corner=2.0),
        text,
    )
    return Group(shapes, classes=(item.kind, "label"), data={"kind": item.kind, **item.hover})


def _centre(name: str, length: int) -> list[Shape]:
    """Return the record's name in bold over its length, centred on the centre."""
    lines = [(name, BOLD, NAME_SIZE)] if BOLD.drawn(name) else []
    lines.append((f"{length} bp", SANS, SMALL_SIZE))
    top = -sum(_height(font, text_size) for _, font, text_size in lines) / 2
    shapes: list[Shape] = []
    for text, font, text_size in lines:
        height = _height(font, text_size)
        width = font.width(text, text_size)
        baseline = _baseline(font, text_size, top + height / 2)
        shapes.append(Text(-width / 2, baseline, text, font, text_size))
        top += height
    return shapes


def _outline(arrow: Arrow) -> str:
    """Return an arrow's outline as path data: its band, with a head at each end it points to."""
    radius, half = arrow.radius, arrow.width / 2
    reach = arrow.width * _OVERHANG / _BAND
    heads = arrow.head_start + arrow.head_end
    head = min(_HEAD / radius, (arrow.end - arrow.start) / heads) if heads else 0.0
    first = arrow.start + head * arrow.head_start
    last = arrow.end - head * arrow.head_end
    outer, inner = radius + half, radius - half
    if arrow.head_start:
        path = _move(radius, arrow.start) + _to(outer + reach, first) + _to(outer, first)
    else:
        path = _move(outer, first)
    path += _arc(outer, first, last)
    if arrow.head_end:
        path += _to(outer + reach, last) + _to(radius, arrow.end) + _to(inner - reach, last)
    path += _to(inner, last) + _arc(inner, last, first)
    if arrow.head_start:
        path += _to(inner - reach, first)
    return path + "Z"


def _move(radius: float, angle: float) -> str:
    x, y = _point(radius, angle)
    return f"M{number(x)} {number(y)}"


def _to(radius: float, angle: float) -> str:
    x, y = _point(radius, angle)
    return f"L{number(x)} {number(y)}"


def _arc(radius: float, start: float, end: float) -> str:
    """Return path data for an arc to `end` from `start`, in pieces of at most a half turn."""
    if end == start:
        return ""
    pieces = math.ceil(abs(end - start) / math.pi + 1e-9)
    sweep = 1 if end > start else 0
    path = ""
    for piece in range(1, pieces + 1):
        x, y = _point(radius, start + (end - start) * piece / pieces)
        path += f"A{number(radius)} {number(radius)} 0 0 {sweep} {number(x)} {number(y)}"
    return path

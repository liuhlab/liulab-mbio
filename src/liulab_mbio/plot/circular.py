"""The circular map: items in, shapes out, in points round a circle centred on (0, 0).

The backbone is two concentric lines, with a bp scale inside them and the record's name and length
in the centre. Each feature is drawn segment by segment as an arrow along its strand, or a box when
it has none, and a segment across the origin is one piece. Features that overlap go on rings
further in, the longest outermost. Each primer is a thin arrow outside the backbone at each binding
site, overlapping ones further out.

A feature's name sits on its arrow when it fits, curved along the band and upright on the lower
half. Every other label is outside, in the columns `labels.columns` lays out, so no label overlaps
another label or the drawing: a feature's name boxed in its colour, a primer's in purple, and a cut
site's enzymes in black, joined to the backbone at the cut.

The columns grow above and below the circle as far as `REACH`. Past it, labels hide in the order
`layers.hiding` gives, and a notice at the bottom right says how many. A name on an arrow never
hides.
"""

import itertools
import json
import math
from collections.abc import Sequence
from dataclasses import dataclass

from liulab_mbio.plot import labels
from liulab_mbio.plot.fonts import BOLD, SANS, Font
from liulab_mbio.plot.labels import Box, Point
from liulab_mbio.plot.layers import (
    Item,
    Span,
    hiding,
    notice,
    outline_color,
    text_color,
    unwrapped,
)
from liulab_mbio.plot.svg import (
    Circle,
    Group,
    Letters,
    Line,
    Path,
    Place,
    Rect,
    Shape,
    Text,
    number,
)
from liulab_mbio.sequence import Strand

#: The backbone's radius, in points. Everything else on the circle is laid out from it.
RADIUS = 220.0

#: How far above or below the centre a label may lie, in points, before labels hide.
REACH = 390.0

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

# A primer's arrow outside the backbone: its thickness, and the room either side of it.
_PRIMER_BAND = 4.0
_PRIMER_GAP = 2.0

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
    head
        The angle each arrowhead takes along the band.
    """

    item: Item
    span: Span
    start: float
    end: float
    radius: float
    width: float
    head_start: bool
    head_end: bool
    head: float


@dataclass(frozen=True, slots=True)
class Name:
    """An item's name set on one of its arrows, letter by letter along the band's middle."""

    arrow: Arrow
    letters: Letters


@dataclass(frozen=True, slots=True)
class Label:
    """An item's label outside the circle, and its leader from the circle to the label's box."""

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
        Every segment of every feature and every binding site of every primer, in the items'
        order.
    names
        Every feature's name set on its arrow, in the items' order.
    labels
        Every label, in the items' order.
    hidden
        The item of each label left out, in the order they hid.
    shapes
        Everything, in the order it is drawn.
    """

    extent: Box
    radius: float
    arrows: tuple[Arrow, ...]
    names: tuple[Name, ...]
    labels: tuple[Label, ...]
    hidden: tuple[Item, ...]
    shapes: tuple[Shape, ...]


def layout(items: Sequence[Item], *, name: str, length: int) -> CircularMap:
    """Lay out `items` as a circular map of a record called `name`, `length` bases long."""
    backbone = RADIUS + _BACKBONE / 2
    features = [item for item in items if item.kind == "feature"]
    primers = _primers([item for item in items if item.kind == "primer"], length, backbone)
    drawing = max(
        [backbone, *(arrow.radius + arrow.width / 2 + _reach(arrow.width) for arrow in primers)]
    )
    order = {id(item): index for index, item in enumerate(items)}
    arrows = tuple(
        sorted((*_arrows(features, length), *primers), key=lambda arrow: order[id(arrow.item)])
    )
    own = {id(item): [arrow for arrow in arrows if arrow.item is item] for item in items}
    names = {
        id(item): on_arrow
        for item in features
        if (on_arrow := _name(item, max(own[id(item)], key=_room)))
    }
    boxed, hidden = _labels([item for item in items if id(item) not in names], length, drawing)
    extent = _extent(drawing + _ANCHOR, boxed)
    said = _notice(hidden, extent)
    if said:
        below = _GAP + _height(SANS, SMALL_SIZE)
        extent = Box(extent.x, extent.y, extent.width, extent.height + below)
    shapes = (
        *_backbone(length),
        *(
            _feature(item, own[id(item)], names.get(id(item)))
            for item in items
            if item.kind != "cut_site"
        ),
        *(_label(label, backbone) for label in boxed),
        *_centre(name, length),
        *said,
    )
    return CircularMap(extent, drawing, arrows, tuple(names.values()), boxed, hidden, shapes)


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
        radius = outer - ring * step
        for index, (span, (start, end)) in enumerate(
            zip(item.spans, unwrapped(item.spans, length), strict=True)
        ):
            head_start = index == 0 and item.strand in (Strand.REVERSE, Strand.BOTH)
            head_end = index == last and item.strand in (Strand.FORWARD, Strand.BOTH)
            arrows.append(
                Arrow(
                    item,
                    span,
                    _angle(start, length),
                    _angle(end, length),
                    radius,
                    width,
                    head_start,
                    head_end,
                    _head(radius, _angle(end - start, length), head_start + head_end),
                )
            )
    return tuple(arrows)


def _primers(items: Sequence[Item], length: int, backbone: float) -> list[Arrow]:
    """Return each primer's arrow outside the backbone, overlapping ones on rings further out."""
    reach = _reach(_PRIMER_BAND)
    step = _PRIMER_BAND + 2 * reach + _PRIMER_GAP
    first = backbone + _PRIMER_GAP + reach + _PRIMER_BAND / 2
    arrows = []
    for item, ring in zip(items, _rings(items, length), strict=True):
        (span,) = item.spans
        radius = first + ring * step
        head_start, head_end = item.strand == Strand.REVERSE, item.strand == Strand.FORWARD
        arrows.append(
            Arrow(
                item,
                span,
                _angle(span.start, length),
                _angle(span.end, length),
                radius,
                _PRIMER_BAND,
                head_start,
                head_end,
                _head(radius, _angle(span.end - span.start, length), head_start + head_end),
            )
        )
    return arrows


def _head(radius: float, turn: float, heads: int) -> float:
    """Return the angle each of `heads` arrowheads takes along a band `turn` radians long."""
    return min(_HEAD / radius, turn / heads) if heads else 0.0


def _reach(width: float) -> float:
    """Return how far the heads of an arrow `width` thick reach beyond its band."""
    return width * _OVERHANG / _BAND


def _room(arrow: Arrow) -> float:
    """Return the angle an arrow's band runs, less its arrowheads."""
    return arrow.end - arrow.start - arrow.head * (arrow.head_start + arrow.head_end)


def _name(item: Item, arrow: Arrow) -> Name | None:
    """Return an item's name set along the middle of `arrow`'s band, or `None` if it does not fit.

    Each letter's box is as wide as its advance and as tall as the face, centred on the band's
    middle and turned square to it. The name fits when its width, padded at each end as a boxed
    label is, spans no more of the band's middle than the arrowheads leave, and no letter's corner
    reaches past the band's outer edge. The padding also holds the end letters' inner corners,
    which reach a little further round than their middles.
    """
    letters = SANS.letters(item.label, LABEL_SIZE)
    if not letters:
        return None
    radius = arrow.radius
    width = letters[-1].x + letters[-1].advance
    widest = max(letter.advance for letter in letters)
    if (
        width + 2 * _PADDING[0] > radius * _room(arrow)
        or math.hypot(radius + _height(SANS, LABEL_SIZE) / 2, widest / 2) > radius + arrow.width / 2
    ):
        return None
    middle = arrow.start + arrow.head * arrow.head_start + _room(arrow) / 2
    upright = math.cos(middle) >= 0
    rise = _baseline(SANS, LABEL_SIZE, 0.0)
    places = []
    for letter in letters:
        along = letter.x + letter.advance / 2 - width / 2
        angle = middle + along / radius if upright else middle - along / radius
        turn = angle if upright else angle + math.pi
        centre = _point(radius, angle)
        # From the letter's middle back half its advance along the baseline, and down to it.
        back = letter.advance / 2
        places.append(
            Place(
                centre.x - back * math.cos(turn) - rise * math.sin(turn),
                centre.y - back * math.sin(turn) + rise * math.cos(turn),
                (math.degrees(turn) + 180) % 360 - 180,
            )
        )
    color = text_color(arrow.span.color)
    return Name(arrow, Letters(tuple(places), item.label, SANS, LABEL_SIZE, color))


def _labels(
    items: Sequence[Item], length: int, drawing: float
) -> tuple[tuple[Label, ...], tuple[Item, ...]]:
    """Return the labels placed in the columns, and the items whose labels hid, in that order."""
    named = [item for item in items if SANS.drawn(item.label)]
    anchored = []
    for item in named:
        start, end = _hull(item, length)
        anchored.append(
            labels.Anchored(
                sum(_font(bold).width(text, LABEL_SIZE) for text, bold in item.runs)
                + 2 * _PADDING[0],
                _height(SANS, LABEL_SIZE) + 2 * _PADDING[1],
                _point(drawing + _ANCHOR, _angle((start + end) / 2, length)),
            )
        )
    order = sorted(range(len(named)), key=lambda i: hiding(named[i]))
    left_out = labels.hide_in_columns(
        anchored, order, radius=drawing, gap=_GAP, spacing=_SPACING, reach=REACH
    )
    gone = set(left_out)
    shown = [i for i in range(len(named)) if i not in gone]
    placed = labels.columns(
        [anchored[i] for i in shown], radius=drawing, gap=_GAP, spacing=_SPACING, margin=_MARGIN
    )
    boxed = tuple(
        Label(named[i], one.box, one.leader) for i, one in zip(shown, placed, strict=True)
    )
    return boxed, tuple(named[i] for i in left_out)


def _extent(reach: float, boxed: Sequence[Label]) -> Box:
    left = min([-reach, *(label.box.x for label in boxed)]) - _EDGE
    right = max([reach, *(label.box.x + label.box.width for label in boxed)]) + _EDGE
    top = min([-reach, *(label.box.y for label in boxed)]) - _EDGE
    bottom = max([reach, *(label.box.y + label.box.height for label in boxed)]) + _EDGE
    return Box(left, top, right - left, bottom - top)


def _notice(hidden: Sequence[Item], extent: Box) -> list[Shape]:
    """Return what says how many labels hid, under the bottom right of `extent`, if any did.

    It carries each hidden label, with its item's kind and type, for a page to list.
    """
    said = notice(hidden)
    if not said:
        return []
    width, height = SANS.width(said, SMALL_SIZE), _height(SANS, SMALL_SIZE)
    middle = extent.y + extent.height - _EDGE + _GAP + height / 2
    text = Text(
        extent.x + extent.width - _EDGE - width,
        _baseline(SANS, SMALL_SIZE, middle),
        said,
        SANS,
        SMALL_SIZE,
        _INK,
    )
    listed = json.dumps(
        [{"label": item.label, "kind": item.kind, "type": item.type} for item in hidden],
        ensure_ascii=False,
    )
    return [Group((text,), classes=("notice",), data={"hidden": listed})]


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


def _feature(item: Item, arrows: Sequence[Arrow], name: Name | None) -> Group:
    """Return an item's arrows, each segment joined to the next by a thin arc across any gap.

    A name set on one of them is drawn over it.
    """
    shapes: list[Shape] = []
    for before, after in itertools.pairwise(arrows):
        if after.start > before.end:
            joint = _move(before.radius, before.end) + _arc(before.radius, before.end, after.start)
            shapes.append(Path(joint, "none", _CONNECTOR, 1.5))
    shapes.extend(
        Path(_outline(arrow), arrow.span.color, outline_color(arrow.span.color), 0.8)
        for arrow in arrows
    )
    if name:
        shapes.append(name.letters)
    return Group(tuple(shapes), classes=(item.kind,), data={"kind": item.kind, **item.hover})


def _font(bold: bool) -> Font:
    return BOLD if bold else SANS


def _label(label: Label, backbone: float) -> Group:
    """Return a label: a feature's name boxed in its colour, any other label's text unboxed.

    A cut site's leader carries on to the backbone, so it points at the cut.
    """
    box, item = label.box, label.item
    (x1, y1), (x2, y2) = label.leader
    shapes: list[Shape] = [Line(x1, y1, x2, y2, _LEADER, 0.8)]
    if item.kind == "cut_site":
        inward = backbone / math.hypot(x1, y1)
        shapes.append(Line(x1 * inward, y1 * inward, x1, y1, _LEADER, 0.8))
    fill = item.color
    if item.kind == "feature":
        shapes.append(Rect(box, item.color, outline_color(item.color), 0.8, corner=2.0))
        fill = text_color(item.color)
    x = box.x + _PADDING[0]
    baseline = _baseline(SANS, LABEL_SIZE, box.y + box.height / 2)
    for text, bold in item.runs:
        shapes.append(Text(x, baseline, text, _font(bold), LABEL_SIZE, fill))
        x += _font(bold).width(text, LABEL_SIZE)
    return Group(
        tuple(shapes), classes=(item.kind, "label"), data={"kind": item.kind, **item.hover}
    )


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
    first = arrow.start + arrow.head * arrow.head_start
    last = arrow.end - arrow.head * arrow.head_end
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

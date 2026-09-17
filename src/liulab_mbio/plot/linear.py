"""The linear map: items in, shapes out, in points along a line running right from (0, 0).

The stretch drawn, a whole record or a region of it, spans the line and keeps the record's
numbering. The backbone is two lines with a bp scale under them, and `• • •` at each end the
molecule carries on past. Features lie under the scale as arrows along their strand, or boxes when
they have none, in rows packed by earliest start. Each primer is a thin arrow above the backbone at
each binding site, overlapping ones further up. Whatever the stretch cuts, it cuts cleanly, with no
arrowhead at the cut, so a feature across the origin of a circular record drawn opened lies at both
ends.

A feature's name goes inside its arrow when it fits, underneath when nothing else in its row lies
there, and otherwise in a box above. Boxed names, primers and cut sites are labelled in the rows
`labels.staircase` lays out above the drawing, each with a vertical leader down to where the item
lies, so no label overlaps another label or the drawing, and no leader crosses another label.

The rows rise above the drawing as far as `RISE`. Past it, labels hide in the order
`layers.hiding` gives, and a notice at the bottom right says how many. A name inside or underneath
a feature never hides.
"""

import dataclasses
import json
from bisect import bisect_left, bisect_right
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass

from liulab_mbio.plot import labels
from liulab_mbio.plot.fonts import BOLD, SANS, Font
from liulab_mbio.plot.labels import Box, Point
from liulab_mbio.plot.layers import (
    Item,
    Piece,
    Span,
    hiding,
    notice,
    outline_color,
    pieces,
    span_text,
    text_color,
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

#: The line's length, in points, however long the stretch it draws.
WIDTH = 720.0

#: How far above the drawing a label may reach, in points, before labels hide.
RISE = 730.0

#: The size of a label's text, in points.
LABEL_SIZE = 11.5
#: The size of the record's name, in points.
NAME_SIZE = 13.5
#: The size of the record's length and of the scale's numbers, in points.
SMALL_SIZE = 10.5

_INK = "#252525"
_LEADER = "#7f7f7f"
_CONNECTOR = "#666666"

# The backbone's two lines, the scale's ticks under them, and how far a number keeps from another
# and from the origin's tick.
_BACKBONE = 2.6
_STRANDS_APART = 3.8
_TICK = 6.0
_ORIGIN_TICK = 12.0
_CLEAR = 10.0

# A feature's arrow: the band's thickness, how far a head reaches beyond the band and along it, the
# room between two bars or names in a row, and between rows.
_BAND = 15.0
_OVERHANG = 4.0
_HEAD = 12.0
_BAR_GAP = 3.0
_ROW_GAP = 4.0

# A primer's arrow above the backbone: its thickness, and the room either side of it.
_PRIMER_BAND = 4.0
_PRIMER_GAP = 2.0

# The labels: padding inside a box, the room between the drawing and the lowest row, the space
# between boxes, and the canvas's margin.
_PADDING = (3.0, 1.0)
_GAP = 10.0
_SPACING = 2.0
_EDGE = 16.0

# The dots at an end: where the first sits beyond the line, how far apart they are, and their size.
_DOT_START = 7.0
_DOT_STEP = 6.0
_DOT = 0.9


@dataclass(frozen=True, slots=True)
class Arrow:
    """One piece of a feature's segment, or of a primer's binding site, as drawn along the line.

    Parameters
    ----------
    item, span
        What it draws.
    start, end
        Where it starts and ends along the line.
    middle
        The height of the band's middle.
    width
        The band's thickness.
    head_start, head_end
        Whether an arrowhead's tip ends the arrow at `start`, or at `end`; never where the
        stretch drawn cuts it.
    head
        How far along the band each arrowhead takes.
    """

    item: Item
    span: Span
    start: float
    end: float
    middle: float
    width: float
    head_start: bool
    head_end: bool
    head: float

    @property
    def body(self) -> Box:
        """The band less its arrowheads."""
        first = self.start + self.head * self.head_start
        last = self.end - self.head * self.head_end
        return Box(first, self.middle - self.width / 2, last - first, self.width)


@dataclass(frozen=True, slots=True)
class Name:
    """A feature's name set on the drawing: inside `arrow` when `inside`, otherwise underneath."""

    arrow: Arrow
    letters: Letters
    inside: bool


@dataclass(frozen=True, slots=True)
class Label:
    """An item's label above the drawing, and its vertical leader up from where the item lies."""

    item: Item
    box: Box
    leader: tuple[Point, Point]


@dataclass(frozen=True, slots=True)
class LinearMap:
    """A linear map laid out: what went where, and the shapes that draw it.

    Parameters
    ----------
    extent
        The canvas.
    top
        How high the drawing reaches. Only labels and their leaders lie higher.
    arrows
        Every piece drawn of every feature and primer, in the items' order.
    names
        Every feature's name set inside or underneath its arrows, in the items' order.
    labels
        Every label, in the items' order.
    hidden
        The item of each label left out, in the order they hid.
    shapes
        Everything, in the order it is drawn.
    """

    extent: Box
    top: float
    arrows: tuple[Arrow, ...]
    names: tuple[Name, ...]
    labels: tuple[Label, ...]
    hidden: tuple[Item, ...]
    shapes: tuple[Shape, ...]


@dataclass(frozen=True, slots=True)
class _Run:
    """Pieces of an item that follow on from one another: one bar, with an arrow for each span."""

    item: Item
    arrows: tuple[Arrow, ...]
    gaps: tuple[tuple[float, float], ...]
    start: float
    end: float


def layout(
    items: Sequence[Item],
    *,
    name: str,
    length: int,
    circular: bool = False,
    span: tuple[int, int] | None = None,
) -> LinearMap:
    """Lay out `items` along a line, as a map of a record called `name`, `length` bases long.

    Parameters
    ----------
    items
        What the record draws; only what lies in the stretch is drawn.
    name, length
        The record's.
    circular
        Whether the record is circular, so the line opens it and the molecule carries on past
        both ends.
    span
        The stretch drawn, 0-based and half-open, ending past `length` across the origin; the
        whole record when ``None``.
    """
    start, end = span if span is not None else (0, length)
    scale = WIDTH / (end - start)

    def at(position: float) -> float:
        return (position - start) * scale

    def runs(kind: str, width: float) -> list[_Run]:
        return [
            run
            for item in items
            if item.kind == kind
            for run in _runs(item, pieces(item, start, end, length, circular=circular), at, width)
        ]

    features = runs("feature", _BAND)
    rows = labels.pack([(run.start, run.end) for run in features], spacing=_BAR_GAP)
    inside, under, boxed = _names(features, rows)
    feature_arrows, names, bottom = _rows(features, rows, inside, under)

    primers = runs("primer", _PRIMER_BAND)
    primer_arrows, top = _primers(primers)

    arrows = {**feature_arrows, **primer_arrows}
    labelled = {id(run) for run in boxed}
    feet: dict[int, list[Point]] = {}
    for run in (*features, *primers):
        if id(run) in labelled or run.item.kind == "primer":
            foot = Point((run.start + run.end) / 2, _foot(run, arrows[id(run)]))
            feet.setdefault(id(run.item), []).append(foot)
    anchored: list[tuple[Item, labels.Anchored]] = []
    for item in items:
        if not SANS.drawn(item.label):
            continue
        if item.kind == "cut_site":
            feet[id(item)] = [
                Point(at(piece.start), -_BACKBONE / 2)
                for piece in pieces(item, start, end, length, circular=circular)
            ]
        anchored.extend((item, _anchored(item, foot)) for foot in feet.get(id(item), ()))
    hides = sorted(range(len(anchored)), key=lambda i: hiding(anchored[i][0]))
    left_out = labels.hide_in_staircase(
        [one for _, one in anchored],
        hides,
        base=top - _GAP,
        spacing=_SPACING,
        ceiling=top - RISE,
    )
    gone = set(left_out)
    shown = [pair for i, pair in enumerate(anchored) if i not in gone]
    placed = labels.staircase([one for _, one in shown], base=top - _GAP, spacing=_SPACING)
    boxes = tuple(
        Label(item, one.box, one.leader) for (item, _), one in zip(shown, placed, strict=True)
    )
    hidden = tuple({id(anchored[i][0]): anchored[i][0] for i in left_out}.values())

    order = {id(item): index for index, item in enumerate(items)}
    drawn = sorted((*features, *primers), key=lambda run: order[id(run.item)])
    scale_shapes, numbers = _scale(start, end, length, circular, at)
    ends = _ends(circular or start > 0, circular or end < length)
    title = _title(name, length, span, bottom)
    shapes = (
        *_backbone(),
        *ends,
        *scale_shapes,
        *_items(drawn, arrows, names),
        *(_label(label) for label in boxes),
        *title,
    )
    taken = [
        Box(0.0, -_BACKBONE / 2, WIDTH, _STRANDS_APART + _BACKBONE),
        *_boxes(ends),
        *numbers,
        *(_bounds(arrow) for run in drawn for arrow in arrows[id(run)]),
        *(_letters_box(name_.letters) for name_ in names.values()),
        *(label.box for label in boxes),
        *_boxes(title),
    ]
    said = _notice(hidden, _extent(taken))
    extent = _extent([*taken, *_boxes(said)])
    return LinearMap(
        extent,
        top,
        tuple(arrow for run in drawn for arrow in arrows[id(run)]),
        tuple(names[id(run)] for run in drawn if id(run) in names),
        boxes,
        hidden,
        (*shapes, *said),
    )


def _runs(
    item: Item, found: Sequence[Piece], at: Callable[[float], float], width: float
) -> list[_Run]:
    """Return an item's bars: its pieces that follow on from one another, each span an arrow."""
    groups: list[list[Piece]] = []
    for piece in found:
        if groups and piece.start == groups[-1][-1].end and not piece.starts:
            groups[-1].append(piece)
        else:
            groups.append([piece])
    bars = []
    for group in groups:
        arrows = tuple(
            _arrow(item, piece, span, at, width) for piece in group if (span := piece.span)
        )
        if arrows:
            gaps = tuple((at(piece.start), at(piece.end)) for piece in group if piece.span is None)
            bars.append(_Run(item, arrows, gaps, at(group[0].start), at(group[-1].end)))
    return bars


def _arrow(
    item: Item, piece: Piece, span: Span, at: Callable[[float], float], width: float
) -> Arrow:
    """Return a piece's arrow, its band's middle at 0, with a head only where the item ends."""
    start, end = at(piece.start), at(piece.end)
    head_start = piece.starts and item.strand in (Strand.REVERSE, Strand.BOTH)
    head_end = piece.ends and item.strand in (Strand.FORWARD, Strand.BOTH)
    heads = head_start + head_end
    head = min(_HEAD, (end - start) / heads) if heads else 0.0
    return Arrow(item, span, start, end, 0.0, width, head_start, head_end, head)


def _room(arrow: Arrow) -> float:
    return arrow.body.width


def _names(
    features: Sequence[_Run], rows: Sequence[int]
) -> tuple[dict[int, Arrow], dict[int, float], list[_Run]]:
    """Return where each bar's name goes: inside an arrow, underneath from a start, or boxed above.

    A name goes inside the arrow with most room when it fits there. Otherwise it goes underneath
    when, centred under its bar, it meets no other bar in its row and no name already set
    underneath there, and otherwise it is boxed.
    """
    inside: dict[int, Arrow] = {}
    under: dict[int, float] = {}
    boxed: list[_Run] = []
    height = _height(SANS, LABEL_SIZE)
    members: dict[int, list[_Run]] = {}
    for run, row in zip(features, rows, strict=True):
        members.setdefault(row, []).append(run)
    for row in sorted(members):
        # Bars in a row are apart, so their ends run in the order of their starts, and names set
        # underneath them in that order are apart in that order too.
        bars = sorted(members[row], key=lambda run: run.start)
        starts, ends = [run.start for run in bars], [run.end for run in bars]
        last = float("-inf")
        for run in bars:
            if not SANS.drawn(run.item.label):
                continue
            width = SANS.width(run.item.label, LABEL_SIZE)
            roomiest = max(run.arrows, key=_room)
            if width + 2 * _PADDING[0] <= _room(roomiest) and height <= roomiest.width:
                inside[id(run)] = roomiest
                continue
            left = (run.start + run.end - width) / 2
            near = bars[
                bisect_right(ends, left - _BAR_GAP) : bisect_left(starts, left + width + _BAR_GAP)
            ]
            if left < last + _BAR_GAP or any(other is not run for other in near):
                boxed.append(run)
            else:
                under[id(run)] = left
                last = left + width
    return inside, under, boxed


def _rows(
    features: Sequence[_Run],
    rows: Sequence[int],
    inside: dict[int, Arrow],
    under: dict[int, float],
) -> tuple[dict[int, tuple[Arrow, ...]], dict[int, Name], float]:
    """Return each bar's arrows in its row under the scale, their names, and where the rows end.

    A row holding a name underneath is as much taller as the name.
    """
    height = _height(SANS, LABEL_SIZE)
    y = _STRANDS_APART + _BACKBONE / 2 + _TICK + 2 * _PADDING[1] + _height(SANS, SMALL_SIZE)
    named = {row for run, row in zip(features, rows, strict=True) if id(run) in under}
    middles: dict[int, float] = {}
    for row in range(max(rows, default=-1) + 1):
        y += _ROW_GAP
        middles[row] = y + _OVERHANG + _BAND / 2
        y += 2 * _OVERHANG + _BAND
        if row in named:
            y += _PADDING[1] + height
    arrows: dict[int, tuple[Arrow, ...]] = {}
    names: dict[int, Name] = {}
    for run, row in zip(features, rows, strict=True):
        middle = middles[row]
        placed = tuple(dataclasses.replace(arrow, middle=middle) for arrow in run.arrows)
        arrows[id(run)] = placed
        label = run.item.label
        if id(run) in inside:
            arrow = placed[run.arrows.index(inside[id(run)])]
            body = arrow.body
            left = body.x + (body.width - SANS.width(label, LABEL_SIZE)) / 2
            letters = _set(label, left, middle, text_color(arrow.span.color))
            names[id(run)] = Name(arrow, letters, inside=True)
        elif id(run) in under:
            below = middle + _BAND / 2 + _OVERHANG + _PADDING[1] + height / 2
            letters = _set(label, under[id(run)], below, _INK)
            names[id(run)] = Name(max(placed, key=_room), letters, inside=False)
    return arrows, names, y


def _primers(primers: Sequence[_Run]) -> tuple[dict[int, tuple[Arrow, ...]], float]:
    """Return each primer's arrows in rows above the backbone, and how high the drawing reaches."""
    reach = _reach(_PRIMER_BAND)
    step = _PRIMER_BAND + 2 * reach + _PRIMER_GAP
    first = -_BACKBONE / 2 - _PRIMER_GAP - reach - _PRIMER_BAND / 2
    rows = labels.pack([(run.start, run.end) for run in primers], spacing=_BAR_GAP)
    arrows = {
        id(run): tuple(
            dataclasses.replace(arrow, middle=first - row * step) for arrow in run.arrows
        )
        for run, row in zip(primers, rows, strict=True)
    }
    top = min(
        [
            -_BACKBONE / 2,
            *(first - row * step - _PRIMER_BAND / 2 - reach for row in rows),
        ]
    )
    return arrows, top


def _reach(width: float) -> float:
    """Return how far the heads of an arrow `width` thick reach beyond its band."""
    return width * _OVERHANG / _BAND


def _height(font: Font, size: float) -> float:
    return (font.ascender - font.descender) / font.units_per_em * size


def _baseline(font: Font, size: float, middle: float) -> float:
    """Return the baseline that centres a line of text on the height `middle`."""
    return middle + (font.ascender + font.descender) / 2 / font.units_per_em * size


def _set(text: str, left: float, middle: float, fill: str) -> Letters:
    """Return `text` set letter by letter from `left`, centred on the height `middle`."""
    baseline = _baseline(SANS, LABEL_SIZE, middle)
    places = tuple(
        Place(left + letter.x, baseline, 0.0) for letter in SANS.letters(text, LABEL_SIZE)
    )
    return Letters(places, text, SANS, LABEL_SIZE, fill)


def _letters_box(letters: Letters) -> Box:
    font, size = letters.font, letters.size
    left, baseline = letters.places[0].x, letters.places[0].y
    top = baseline - font.ascender / font.units_per_em * size
    return Box(left, top, font.width(letters.text, size), _height(font, size))


def _foot(run: _Run, arrows: Sequence[Arrow]) -> float:
    """Return where a bar's leader starts: the top of a primer's arrow, or the backbone's top."""
    if run.item.kind == "primer":
        return min(arrow.middle - arrow.width / 2 - _reach(arrow.width) for arrow in arrows)
    return -_BACKBONE / 2


def _anchored(item: Item, anchor: Point) -> labels.Anchored:
    width = sum(_font(bold).width(text, LABEL_SIZE) for text, bold in item.runs)
    return labels.Anchored(
        width + 2 * _PADDING[0], _height(SANS, LABEL_SIZE) + 2 * _PADDING[1], anchor
    )


def _font(bold: bool) -> Font:
    return BOLD if bold else SANS


def _backbone() -> list[Shape]:
    return [
        Line(0.0, 0.0, WIDTH, 0.0, _INK, _BACKBONE),
        Line(0.0, _STRANDS_APART, WIDTH, _STRANDS_APART, _INK, _BACKBONE),
    ]


def _ends(left: bool, right: bool) -> list[Shape]:
    """Return `• • •` beyond each end of the line the molecule carries on past."""
    middle = _STRANDS_APART / 2
    offsets = [_DOT_START + index * _DOT_STEP for index in range(3)]
    dots = [
        *(Circle(-offset, middle, _DOT, _INK, 2 * _DOT) for offset in offsets if left),
        *(Circle(WIDTH + offset, middle, _DOT, _INK, 2 * _DOT) for offset in offsets if right),
    ]
    return [Group(tuple(dots), classes=("ends",))] if dots else []


def _scale(
    start: int, end: int, length: int, circular: bool, at: Callable[[float], float]
) -> tuple[list[Shape], list[Box]]:
    """Return the ticks and numbers under the backbone, and the boxes the numbers take.

    Numbers read as the record numbers its bases. A tick marks the end of each base whose number
    is a multiple of the step, and a longer one the origin, where the stretch runs across it.
    """
    below = _STRANDS_APART + _BACKBONE / 2
    middle = below + _TICK + _PADDING[1] + _height(SANS, SMALL_SIZE) / 2
    origin = at(length) if circular and start < length < end else None
    shapes: list[Shape] = []
    boxes: list[Box] = []
    if origin is not None:
        shapes.append(Line(origin, below, origin, below + _ORIGIN_TICK, _INK, 2.4))
    step = _step(end - start)
    right = float("-inf")
    for lap in (0, 1) if circular else (0,):
        offset = lap * length
        first = (max(start, offset) - offset) // step * step + step
        for position in range(offset + first, min(end, offset + length) + 1, step):
            x = at(position)
            shapes.append(Line(x, below, x, below + _TICK, _INK, 1.2))
            text = str(position - offset)
            width = SANS.width(text, SMALL_SIZE)
            if x - width / 2 < right + _CLEAR or (
                origin is not None and abs(x - origin) < width / 2 + _CLEAR
            ):
                continue
            baseline = _baseline(SANS, SMALL_SIZE, middle)
            shapes.append(Text(x - width / 2, baseline, text, SANS, SMALL_SIZE, _INK))
            boxes.extend(_boxes(shapes[-1:]))
            right = x + width / 2
    return [Group(tuple(shapes), classes=("scale",))], boxes


def _step(bases: int) -> int:
    """Return the least interval of 1, 2 or 5 times a power of ten marking at most ten numbers."""
    power = 1
    while True:
        for factor in (1, 2, 5):
            if factor * power * 10 >= bases:
                return factor * power
        power *= 10


def _items(
    runs: Sequence[_Run], arrows: dict[int, tuple[Arrow, ...]], names: dict[int, Name]
) -> list[Group]:
    """Return a group for each item: its arrows, the thin lines across its gaps, and its names."""
    groups: dict[int, tuple[Item, list[Shape]]] = {}
    for run in runs:
        _, shapes = groups.setdefault(id(run.item), (run.item, []))
        middle = arrows[id(run)][0].middle
        shapes.extend(Line(low, middle, high, middle, _CONNECTOR, 1.5) for low, high in run.gaps)
        shapes.extend(
            Path(_outline(arrow), arrow.span.color, outline_color(arrow.span.color), 0.8)
            for arrow in arrows[id(run)]
        )
        if id(run) in names:
            shapes.append(names[id(run)].letters)
    return [
        Group(tuple(shapes), classes=(item.kind,), data={"kind": item.kind, **item.hover})
        for item, shapes in groups.values()
    ]


def _label(label: Label) -> Group:
    """Return a label: a feature's name boxed in its colour, any other label's text unboxed."""
    box, item = label.box, label.item
    (x1, y1), (x2, y2) = label.leader
    shapes: list[Shape] = [Line(x1, y1, x2, y2, _LEADER, 0.8)]
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


def _title(name: str, length: int, span: tuple[int, int] | None, bottom: float) -> list[Shape]:
    """Return the record's name in bold over its length, or the region drawn, under the drawing."""
    lines = [(name, BOLD, NAME_SIZE)] if BOLD.drawn(name) else []
    if span is None:
        lines.append((f"{length} bp", SANS, SMALL_SIZE))
    else:
        lines.append((f"{span_text(*span, length)} ({span[1] - span[0]} bp)", SANS, SMALL_SIZE))
    top = bottom + _GAP
    shapes: list[Shape] = []
    for text, font, size in lines:
        height = _height(font, size)
        left = (WIDTH - font.width(text, size)) / 2
        shapes.append(Text(left, _baseline(font, size, top + height / 2), text, font, size))
        top += height
    return shapes


def _notice(hidden: Sequence[Item], extent: Box) -> list[Shape]:
    """Return what says how many labels hid, under the bottom right of `extent`, if any did.

    It carries the hidden labels for a page to list.
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
    listed = json.dumps([item.label for item in hidden], ensure_ascii=False)
    return [Group((text,), classes=("notice",), data={"hidden": listed})]


def _boxes(shapes: Iterable[Shape]) -> list[Box]:
    """Return the box each dot or line of text takes."""
    boxes = []
    for shape in shapes:
        match shape:
            case Group(inner):
                boxes.extend(_boxes(inner))
            case Circle(cx, cy, r, _, width):
                reach = r + width / 2
                boxes.append(Box(cx - reach, cy - reach, 2 * reach, 2 * reach))
            case Text(x, y, text, font, size):
                top = y - font.ascender / font.units_per_em * size
                boxes.append(Box(x, top, font.width(text, size), _height(font, size)))
            case _:
                pass
    return boxes


def _bounds(arrow: Arrow) -> Box:
    """Return the box an arrow takes, its heads included."""
    reach = _reach(arrow.width)
    return Box(
        arrow.start,
        arrow.middle - arrow.width / 2 - reach,
        arrow.end - arrow.start,
        arrow.width + 2 * reach,
    )


def _extent(boxes: Sequence[Box]) -> Box:
    left = min(box.x for box in boxes) - _EDGE
    right = max(box.x + box.width for box in boxes) + _EDGE
    top = min(box.y for box in boxes) - _EDGE
    bottom = max(box.y + box.height for box in boxes) + _EDGE
    return Box(left, top, right - left, bottom - top)


def _outline(arrow: Arrow) -> str:
    """Return an arrow's outline as path data: its band, with a head at each end it points to."""
    half, reach, y = arrow.width / 2, _reach(arrow.width), arrow.middle
    body = arrow.body
    first, last = body.x, body.x + body.width
    if arrow.head_start:
        path = _move(arrow.start, y) + _to(first, y - half - reach) + _to(first, y - half)
    else:
        path = _move(first, y - half)
    path += _to(last, y - half)
    if arrow.head_end:
        path += _to(last, y - half - reach) + _to(arrow.end, y) + _to(last, y + half + reach)
    path += _to(last, y + half) + _to(first, y + half)
    if arrow.head_start:
        path += _to(first, y + half + reach)
    return path + "Z"


def _move(x: float, y: float) -> str:
    return f"M{number(x)} {number(y)}"


def _to(x: float, y: float) -> str:
    return f"L{number(x)} {number(y)}"

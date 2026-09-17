"""The sequence view: items and bases in, shapes out, in points, one row block under another.

The stretch drawn, a whole record or a region of it, is set in rows of a fixed number of bases,
numbered as the record numbers them. Each row is a block of its own: a position ruler numbering
every tenth base, the top strand, a rail with a tick at every base, and the bottom strand, with
the position of the row's last base at its right.

Under the bases each feature lies as a bar, in tracks packed by earliest start: pointed where the
feature ends on its strand, cut flat where a row cuts it, and a thin line across the gap a joined
feature leaves. A CDS's translation runs above its bar, each amino acid's three letters centred on
its codon's middle base, a stop in red. A feature's name goes inside its bar when it fits,
underneath when nothing else in its track lies there, and otherwise in a box above the ruler, in
the rows `labels.staircase` lays out, as the linear map's are. A row grows to hold what it draws,
so nothing is hidden.
"""

import dataclasses
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import NamedTuple

from liulab_mbio.plot import labels
from liulab_mbio.plot.fonts import BOLD, MONO, SANS, Font
from liulab_mbio.plot.labels import Box, Point
from liulab_mbio.plot.layers import (
    Codon,
    Item,
    Piece,
    Span,
    outline_color,
    pieces,
    text_color,
)
from liulab_mbio.plot.svg import Group, Letters, Line, Path, Place, Rect, Shape, Text, number
from liulab_mbio.sequence import Strand, reverse_complement

#: How many bases a sequence view draws at most.
LIMIT = 100_000

#: How wide each base is, in points.
CELL = 11.0

#: The size of a base, in points.
BASE_SIZE = 11.0
#: The size of an amino acid's three letters, in points.
RESIDUE_SIZE = 10.0
#: The size of a label's text, in points.
LABEL_SIZE = 11.5
#: The size of the ruler's numbers and of each row's last position, in points.
SMALL_SIZE = 10.5

#: The colour of a stop codon.
STOP = "#CC3311"

_INK = "#252525"
_RAIL = "#7f7f7f"
_LEADER = "#7f7f7f"
_CONNECTOR = "#666666"

# The bases: the room under the ruler and either side of the rail, how far the rail's ticks reach
# either side of it at a base, a fifth and a tenth, and the room before a row's last position.
_RULER_GAP = 2.0
_RAIL_GAP = 2.0
_TICKS = (1.5, 2.5, 3.5)
_NUMBER_GAP = 14.0

# A feature's bar: its thickness, how far its point runs along it, the room between bars in a
# track, between tracks, under the bases and under a translation.
_BAR = 15.0
_POINT = 7.0
_BAR_GAP = 3.0
_TRACK_GAP = 3.0
_FEATURE_GAP = 4.0
_RESIDUE_GAP = 1.0

# The labels: padding inside a box, the room between the ruler and the lowest box, and between
# boxes; the room between rows, and the canvas's margin.
_PADDING = (3.0, 1.0)
_GAP = 6.0
_SPACING = 2.0
_ROW_GAP = 16.0
_EDGE = 16.0


@dataclass(frozen=True, slots=True)
class Bar:
    """One piece of a feature's segment in a row, as drawn under the bases.

    Parameters
    ----------
    item, span
        What it draws.
    start, end
        Where it starts and ends along the row.
    middle
        The height of the bar's middle.
    point_start, point_end
        Whether the bar comes to a point at `start`, or at `end`: only where the feature ends on
        its strand, never where the row cuts it.
    point
        How far along the bar each point runs.
    """

    item: Item
    span: Span
    start: float
    end: float
    middle: float
    point_start: bool
    point_end: bool
    point: float

    @property
    def body(self) -> Box:
        """The bar less its points."""
        first = self.start + self.point * self.point_start
        last = self.end - self.point * self.point_end
        return Box(first, self.middle - _BAR / 2, last - first, _BAR)


@dataclass(frozen=True, slots=True)
class Name:
    """A feature's name set in a row: inside `bar` when `inside`, otherwise underneath it."""

    bar: Bar
    letters: Letters
    inside: bool


@dataclass(frozen=True, slots=True)
class Label:
    """An item's label above a row's ruler, and its vertical leader up from where the item lies."""

    item: Item
    box: Box
    leader: tuple[Point, Point]


@dataclass(frozen=True, slots=True)
class Translation:
    """A CDS's amino acids in a row, three letters each: `letters` in ink, `stops` in red."""

    item: Item
    letters: Letters | None
    stops: Letters | None


@dataclass(frozen=True, slots=True)
class Row:
    """One row block of a sequence view: the bases it holds, and what went where in it.

    Parameters
    ----------
    start, end
        The bases it holds, 0-based and half-open, counted as the stretch drawn counts them.
    extent
        The box everything the row draws takes, its labels included.
    bases
        The box its ruler, strands and rail take.
    bars, names, translations, labels
        What it draws of each item, in order along the row.
    shapes
        Everything in the row, in the order it is drawn.
    """

    start: int
    end: int
    extent: Box
    bases: Box
    bars: tuple[Bar, ...]
    names: tuple[Name, ...]
    translations: tuple[Translation, ...]
    labels: tuple[Label, ...]
    shapes: tuple[Shape, ...]


@dataclass(frozen=True, slots=True)
class SequenceView:
    """A sequence view laid out: its rows, top to bottom, and the shapes that draw them.

    Parameters
    ----------
    extent
        The canvas.
    rows
        Each row block, in order down the view, none reaching into the next.
    shapes
        Each row's shapes, grouped by row.
    """

    extent: Box
    rows: tuple[Row, ...]
    shapes: tuple[Shape, ...]


@dataclass(frozen=True, slots=True)
class _Run:
    """Pieces of an item in a row that follow on from one another: one bar, gaps included."""

    item: Item
    pieces: tuple[Piece, ...]
    start: float
    end: float
    bars: tuple[Bar, ...]
    codons: tuple[tuple[int, Codon], ...]


def layout(
    items: Sequence[Item],
    *,
    bases: str,
    circular: bool = False,
    span: tuple[int, int] | None = None,
    bases_per_row: int = 60,
    both_strands: bool = True,
) -> SequenceView:
    """Lay out `items` over `bases`, a record's sequence, in rows of `bases_per_row` bases.

    Parameters
    ----------
    items
        What the record draws; the features that lie in the stretch are drawn.
    bases
        The record's whole sequence, whatever stretch is drawn.
    circular
        Whether the record is circular.
    span
        The stretch drawn, 0-based and half-open, ending past the record's length across the
        origin; the whole record when ``None``.
    bases_per_row
        How many bases a row holds, the last row perhaps fewer. At least 1.
    both_strands
        Whether the bottom strand is drawn under the top one.
    """
    length = len(bases)
    start, end = span if span is not None else (0, length)
    bounds = [
        (first, min(first + bases_per_row, end)) for first in range(start, end, bases_per_row)
    ]
    placed: list[list[tuple[Item, list[Piece]]]] = [[] for _ in bounds]
    read: list[dict[int, list[tuple[int, Codon]]]] = [{} for _ in bounds]
    for item in items:
        if item.kind != "feature":
            continue
        found = pieces(item, start, end, length, circular=circular)
        for index, cut in _by_row(found, start, bases_per_row).items():
            placed[index].append((item, cut))
        for codon in item.translation:
            # A codon's middle lies in the stretch once at most, since it is less than one turn.
            for at in (codon.middle, codon.middle + length):
                if start <= at < end:
                    codons = read[(at - start) // bases_per_row].setdefault(id(item), [])
                    codons.append((at, codon))
    widest = max(SANS.width(_last(last, length), SMALL_SIZE) for _, last in bounds)
    right = bases_per_row * CELL + _NUMBER_GAP + widest
    rows: list[Row] = []
    top = 0.0
    before: dict[int, int] = {}
    for (first, last), entries, codons in zip(bounds, placed, read, strict=True):
        row, before = _row(
            first,
            last,
            bases,
            entries,
            codons,
            before,
            top=top,
            right=right,
            width=bases_per_row * CELL,
            both_strands=both_strands,
        )
        rows.append(row)
        top = row.extent.y + row.extent.height + _ROW_GAP
    extent = _extent([row.extent for row in rows], margin=_EDGE)
    shapes = tuple(Group(row.shapes, classes=("row",)) for row in rows)
    return SequenceView(extent, tuple(rows), shapes)


def _by_row(found: Sequence[Piece], start: int, bases_per_row: int) -> dict[int, list[Piece]]:
    """Return what of each piece lies in each row, by the row's index, cut where rows meet."""
    rows: dict[int, list[Piece]] = {}
    for piece in found:
        if piece.start == piece.end:
            continue
        for index in range(
            (piece.start - start) // bases_per_row, (piece.end - 1 - start) // bases_per_row + 1
        ):
            low = max(piece.start, start + index * bases_per_row)
            high = min(piece.end, start + (index + 1) * bases_per_row)
            rows.setdefault(index, []).append(
                Piece(
                    low,
                    high,
                    piece.span,
                    piece.starts and low == piece.start,
                    piece.ends and high == piece.end,
                )
            )
    return rows


def _row(
    first: int,
    last: int,
    bases: str,
    entries: Sequence[tuple[Item, list[Piece]]],
    codons: dict[int, list[tuple[int, Codon]]],
    before: Mapping[int, int],
    *,
    top: float,
    right: float,
    width: float,
    both_strands: bool,
) -> tuple[Row, dict[int, int]]:
    """Lay out the row holding bases `first` to `last`, its top at `top`.

    `before` holds the track of each item the row before ran on into this one, by the item's id.
    Bars starting together are packed in the order of those tracks, so the items keep their order
    down the rows. Returns the row, and the same for the row after.
    """

    def at(position: float) -> float:
        return (position - first) * CELL

    runs = [
        run for item, found in entries for run in _runs(item, found, at, codons.get(id(item), ()))
    ]
    # `pack` takes bars starting together in the order given.
    runs.sort(key=lambda run: (run.start, before.get(id(run.item), len(before))))
    tracks = labels.pack([(run.start, run.end) for run in runs], spacing=_BAR_GAP)
    after = {
        id(run.item): track
        for run, track in zip(runs, tracks, strict=True)
        if run.pieces[-1].end == last and not run.pieces[-1].ends
    }
    inside, under, boxed = _names(runs, tracks, width)
    down = _down(runs, tracks, under, both_strands)
    anchored = [_anchored(run.item, Point((run.start + run.end) / 2, 0.0)) for run in boxed]
    stairs = labels.staircase(anchored, base=-_GAP, spacing=_SPACING)
    # Everything so far lies down from the ruler's top at 0; the row moves down to `top`.
    shift = top - min([0.0, *(one.box.y for one in stairs)])

    groups: dict[int, tuple[Item, list[Shape]]] = {}
    bars: list[Bar] = []
    names: list[Name] = []
    translations: list[Translation] = []
    for run, track in zip(runs, tracks, strict=True):
        middle = down.middles[track] + shift
        _, shapes = groups.setdefault(id(run.item), (run.item, []))
        mine = tuple(dataclasses.replace(bar, middle=middle) for bar in run.bars)
        bars.extend(mine)
        shapes.extend(
            Line(at(piece.start), middle, at(piece.end), middle, _CONNECTOR, 1.2)
            for piece in run.pieces
            if piece.span is None
        )
        shapes.extend(
            Path(_outline(bar), bar.span.color, outline_color(bar.span.color), 0.8) for bar in mine
        )
        if run.codons:
            translations.append(_translation(run, first, down.residues[track] + shift))
            lines = (translations[-1].letters, translations[-1].stops)
            shapes.append(Group(tuple(line for line in lines if line), classes=("translation",)))
        label = run.item.label
        if id(run) in inside:
            bar = mine[run.bars.index(inside[id(run)])]
            left = bar.body.x + (bar.body.width - SANS.width(label, LABEL_SIZE)) / 2
            names.append(Name(bar, _set(label, left, middle, text_color(bar.span.color)), True))
            shapes.append(names[-1].letters)
        elif id(run) in under:
            below = middle + _BAR / 2 + _PADDING[1] + _height(SANS, LABEL_SIZE) / 2
            bar = max(mine, key=lambda one: one.body.width)
            names.append(Name(bar, _set(label, under[id(run)], below, _INK), False))
            shapes.append(names[-1].letters)
    placed = tuple(
        Label(
            run.item,
            Box(one.box.x, one.box.y + shift, one.box.width, one.box.height),
            (
                Point(one.leader[0].x, one.leader[0].y + shift),
                Point(one.leader[1].x, one.leader[1].y + shift),
            ),
        )
        for run, one in zip(boxed, stairs, strict=True)
    )

    drawn, boxes = _bases_drawn(bases, first, last, down, shift, right, both_strands)
    boxes += [
        Box(0.0, shift, 0.0, down.bottom),
        *(Box(bar.start, bar.body.y, bar.end - bar.start, _BAR) for bar in bars),
        *(_letters_box(name.letters) for name in names),
        *(_letters_box(line) for one in translations for line in (one.letters, one.stops) if line),
        *(label.box for label in placed),
    ]
    shapes = (
        *drawn,
        *(
            Group(tuple(shapes), classes=(item.kind,), data={"kind": item.kind, **item.hover})
            for item, shapes in groups.values()
        ),
        *(_label(label) for label in placed),
    )
    row = Row(
        first,
        last,
        _extent(boxes, margin=0.0),
        Box(0.0, shift, at(last), down.bases),
        tuple(bars),
        tuple(names),
        tuple(translations),
        placed,
        shapes,
    )
    return row, after


class _Down(NamedTuple):
    """How far down a row, from the top of its ruler, each part of it lies.

    Parameters
    ----------
    strand, rail, complement
        The top strand's top, the rail, and the bottom strand's top.
    bases
        Where the bases end.
    residues, middles
        Each track's amino acids' middle, where it has any, and its bars' middle.
    bottom
        Where the row ends.
    """

    strand: float
    rail: float
    complement: float
    bases: float
    residues: dict[int, float]
    middles: dict[int, float]
    bottom: float


def _down(
    runs: Sequence[_Run], tracks: Sequence[int], under: dict[int, float], both_strands: bool
) -> _Down:
    """Return how far down a row each part lies.

    A track holding a translation, or a name underneath, is as much taller.
    """
    strand = _height(SANS, SMALL_SIZE) + _RULER_GAP
    rail = strand + _height(MONO, BASE_SIZE) + _RAIL_GAP + _TICKS[-1]
    complement = rail + _TICKS[-1] + _RAIL_GAP
    bases = complement + _height(MONO, BASE_SIZE) if both_strands else rail + _TICKS[-1]
    named = {track for run, track in zip(runs, tracks, strict=True) if id(run) in under}
    translated = {track for run, track in zip(runs, tracks, strict=True) if run.codons}
    residues: dict[int, float] = {}
    middles: dict[int, float] = {}
    y = bases + _FEATURE_GAP
    for track in range(max(tracks, default=-1) + 1):
        y += _TRACK_GAP
        if track in translated:
            residues[track] = y + _height(MONO, RESIDUE_SIZE) / 2
            y += _height(MONO, RESIDUE_SIZE) + _RESIDUE_GAP
        middles[track] = y + _BAR / 2
        y += _BAR
        if track in named:
            y += _PADDING[1] + _height(SANS, LABEL_SIZE)
    return _Down(strand, rail, complement, bases, residues, middles, y)


def _bases_drawn(
    bases: str,
    first: int,
    last: int,
    down: _Down,
    shift: float,
    right: float,
    both_strands: bool,
) -> tuple[list[Shape], list[Box]]:
    """Return a row's ruler, strands and last position, moved down by `shift`, and their boxes.

    The last position is set flush right at `right`, level with the rail.
    """
    length = len(bases)
    top = _bases(bases, first, last)
    ruler, boxes = _ruler(first, last, length, shift, down.rail + shift)
    strand = _strand(top, down.strand + shift)
    drawn: list[Shape] = [ruler, Group((strand,), classes=("strand", "top"))]
    boxes.append(_letters_box(strand))
    if both_strands:
        complement = _strand(reverse_complement(top)[::-1], down.complement + shift)
        drawn.append(Group((complement,), classes=("strand", "bottom")))
        boxes.append(_letters_box(complement))
    position = _last(last, length)
    left = right - SANS.width(position, SMALL_SIZE)
    baseline = _baseline(SANS, SMALL_SIZE, down.rail + shift)
    text = Text(left, baseline, position, SANS, SMALL_SIZE, _INK)
    drawn.append(Group((text,), classes=("position",)))
    boxes.append(_text_box(text))
    return drawn, boxes


def _runs(
    item: Item,
    found: Sequence[Piece],
    at: Callable[[float], float],
    codons: Sequence[tuple[int, Codon]],
) -> list[_Run]:
    """Return an item's runs in a row: its pieces that follow on, each span a bar."""
    groups: list[list[Piece]] = []
    for piece in found:
        if groups and piece.start == groups[-1][-1].end and not piece.starts:
            groups[-1].append(piece)
        else:
            groups.append([piece])
    runs = []
    for group in groups:
        low, high = group[0].start, group[-1].end
        runs.append(
            _Run(
                item,
                tuple(group),
                at(low),
                at(high),
                tuple(_bar(item, piece, span, at) for piece in group if (span := piece.span)),
                tuple(sorted((one for one in codons if low <= one[0] < high), key=_position)),
            )
        )
    return runs


def _position(codon: tuple[int, Codon]) -> int:
    return codon[0]


def _bar(item: Item, piece: Piece, span: Span, at: Callable[[float], float]) -> Bar:
    """Return a piece's bar, its middle at 0, pointed only where the item ends on its strand."""
    start, end = at(piece.start), at(piece.end)
    point_start = piece.starts and item.strand in (Strand.REVERSE, Strand.BOTH)
    point_end = piece.ends and item.strand in (Strand.FORWARD, Strand.BOTH)
    points = point_start + point_end
    point = min(_POINT, (end - start) / points) if points else 0.0
    return Bar(item, span, start, end, 0.0, point_start, point_end, point)


def _names(
    runs: Sequence[_Run], tracks: Sequence[int], width: float
) -> tuple[dict[int, Bar], dict[int, float], list[_Run]]:
    """Return where each run's name goes: inside a bar, underneath from a start, or boxed above.

    A name goes inside the bar with most room when it fits there. Otherwise it goes underneath,
    centred under its run and kept within the row, when it meets no other run in its track and no
    name already set underneath there; otherwise it is boxed.
    """
    inside: dict[int, Bar] = {}
    under: dict[int, float] = {}
    boxed: list[_Run] = []
    members: dict[int, list[_Run]] = {}
    for run, track in zip(runs, tracks, strict=True):
        members.setdefault(track, []).append(run)
    for track in sorted(members):
        # Runs in a track are apart, so a name set underneath ends past every one set before it.
        ordered = sorted(members[track], key=lambda run: run.start)
        last = float("-inf")
        for run in ordered:
            if not run.bars or not SANS.drawn(run.item.label):
                continue
            text = SANS.width(run.item.label, LABEL_SIZE)
            roomiest = max(run.bars, key=lambda bar: bar.body.width)
            if text + 2 * _PADDING[0] <= roomiest.body.width:
                inside[id(run)] = roomiest
                continue
            left = (run.start + run.end - text) / 2
            if text <= width:
                left = min(max(left, 0.0), width - text)
            if left < last + _BAR_GAP or any(
                other is not run
                and other.start < left + text + _BAR_GAP
                and left - _BAR_GAP < other.end
                for other in ordered
            ):
                boxed.append(run)
            else:
                under[id(run)] = left
                last = left + text
    return inside, under, boxed


def _translation(run: _Run, first: int, middle: float) -> Translation:
    """Return a run's amino acids, each centred on its codon's middle base, on the height `middle`."""
    placed: dict[bool, list[tuple[float, str]]] = {False: [], True: []}
    for position, codon in run.codons:
        placed[codon.stop].append(((position - first + 0.5) * CELL, codon.name))
    return Translation(
        run.item, _residues(placed[False], middle, _INK), _residues(placed[True], middle, STOP)
    )


def _residues(placed: Sequence[tuple[float, str]], middle: float, fill: str) -> Letters | None:
    """Return three-letter amino acids as one line, each centred where it is placed."""
    return _centred(placed, MONO, RESIDUE_SIZE, middle, fill) if placed else None


def _centred(
    placed: Sequence[tuple[float, str]], font: Font, size: float, middle: float, fill: str
) -> Letters:
    """Return words as one line of text, each centred where it is placed, a space between two."""
    baseline = _baseline(font, size, middle)
    places: list[Place] = []
    end = 0.0
    for index, (centre, word) in enumerate(placed):
        letters = font.letters(word, size)
        width = letters[-1].x + letters[-1].advance
        left = centre - width / 2
        if index:
            places.append(Place(end, baseline, 0.0))
        places.extend(Place(left + letter.x, baseline, 0.0) for letter in letters)
        end = left + width
    return Letters(tuple(places), " ".join(word for _, word in placed), font, size, fill)


def _ruler(first: int, last: int, length: int, top: float, rail: float) -> tuple[Group, list[Box]]:
    """Return a row's ruler and rail, and the box the ruler's numbers take.

    A number stands over every base the record numbers a multiple of ten. The rail ticks every
    base, longer at a multiple of five and longer still at one of ten.
    """
    numbers: list[tuple[float, str]] = []
    ticks = []
    reach = 0.0
    for position in range(first, last):
        shown = position % length + 1
        x = (position - first + 0.5) * CELL
        before, reach = reach, _TICKS[(shown % 5 == 0) + (shown % 10 == 0)]
        if position == first:
            ticks.append(f"M{number(x)} {number(rail - reach)}v{number(2 * reach)}")
        else:
            # On from the foot of the tick before.
            ticks.append(f"m{number(CELL)} {number(-before - reach)}v{number(2 * reach)}")
        if shown % 10 == 0:
            numbers.append((x, str(shown)))
    shapes: list[Shape] = [
        Line(0.0, rail, (last - first) * CELL, rail, _RAIL, 1.6),
        Path("".join(ticks), "none", _RAIL, 1.0),
    ]
    if not numbers:
        return Group(tuple(shapes), classes=("ruler",)), []
    text = _centred(numbers, SANS, SMALL_SIZE, top + _height(SANS, SMALL_SIZE) / 2, _INK)
    return Group((text, *shapes), classes=("ruler",)), [_letters_box(text)]


def _strand(text: str, top: float) -> Letters:
    """Return a strand's bases set one to a cell, its top at `top`."""
    baseline = _baseline(MONO, BASE_SIZE, top + _height(MONO, BASE_SIZE) / 2)
    offset = (CELL - MONO.width("A", BASE_SIZE)) / 2
    places = tuple(Place(index * CELL + offset, baseline, 0.0) for index in range(len(text)))
    return Letters(places, text, MONO, BASE_SIZE, _INK)


def _bases(bases: str, first: int, last: int) -> str:
    """Return the top strand from `first` to `last`, across the origin as needed."""
    length = len(bases)
    return "".join(bases[position % length] for position in range(first, last))


def _last(last: int, length: int) -> str:
    """Return the position of the base before `last`, as a person reads it."""
    return str((last - 1) % length + 1)


def _anchored(item: Item, anchor: Point) -> labels.Anchored:
    width = sum(_font(bold).width(text, LABEL_SIZE) for text, bold in item.runs)
    return labels.Anchored(
        width + 2 * _PADDING[0], _height(SANS, LABEL_SIZE) + 2 * _PADDING[1], anchor
    )


def _font(bold: bool) -> Font:
    return BOLD if bold else SANS


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


def _set(text: str, left: float, middle: float, fill: str) -> Letters:
    """Return `text` set letter by letter from `left`, centred on the height `middle`."""
    baseline = _baseline(SANS, LABEL_SIZE, middle)
    places = tuple(
        Place(left + letter.x, baseline, 0.0) for letter in SANS.letters(text, LABEL_SIZE)
    )
    return Letters(places, text, SANS, LABEL_SIZE, fill)


def _height(font: Font, size: float) -> float:
    return (font.ascender - font.descender) / font.units_per_em * size


def _baseline(font: Font, size: float, middle: float) -> float:
    """Return the baseline that centres a line of text on the height `middle`."""
    return middle + (font.ascender + font.descender) / 2 / font.units_per_em * size


def _text_box(text: Text) -> Box:
    font, size = text.font, text.size
    top = text.y - font.ascender / font.units_per_em * size
    return Box(text.x, top, font.width(text.text, size), _height(font, size))


def _letters_box(letters: Letters) -> Box:
    """Return the box a line of letters takes, from its first letter's start to its last's end."""
    font, size = letters.font, letters.size
    drawn = font.letters(letters.text, size)
    left, baseline = letters.places[0].x, letters.places[0].y
    top = baseline - font.ascender / font.units_per_em * size
    return Box(left, top, letters.places[-1].x + drawn[-1].advance - left, _height(font, size))


def _extent(boxes: Sequence[Box], *, margin: float) -> Box:
    left = min(box.x for box in boxes) - margin
    right = max(box.x + box.width for box in boxes) + margin
    top = min(box.y for box in boxes) - margin
    bottom = max(box.y + box.height for box in boxes) + margin
    return Box(left, top, right - left, bottom - top)


def _outline(bar: Bar) -> str:
    """Return a bar's outline as path data, coming to a point at each end it points to."""
    half, y = _BAR / 2, bar.middle
    body = bar.body
    first, last = body.x, body.x + body.width
    path = f"M{number(first)} {number(y - half)}H{number(last)}"
    if bar.point_end:
        path += f"L{number(bar.end)} {number(y)}L{number(last)} {number(y + half)}"
    else:
        path += f"V{number(y + half)}"
    path += f"H{number(first)}"
    if bar.point_start:
        path += f"L{number(bar.start)} {number(y)}"
    return path + "Z"

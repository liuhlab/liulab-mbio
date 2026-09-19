"""The sequence view: items and bases in, shapes out, in points, one row block under another.

The stretch drawn, a whole record or a region of it, is set in rows of a fixed number of bases,
numbered as the record numbers them. Each row is a block of its own: a position ruler numbering
every tenth base, the top strand, a rail with a tick at every base, and the bottom strand, with
the position of the row's last base at its right.

Each primer lies as an arrow along its binding site, beside the strand whose bases it spells: a
forward primer's between the ruler and the top strand, a reverse primer's under the bottom strand,
in tracks packed by earliest start. Its 5' tail runs on from the arrow's back a base to a cell,
bent away from the bases, and each base it does not pair with is marked in red. A cut site is drawn
through the top strand where its enzymes cut it, along the rail, and through the bottom strand
where each cuts that, so its overhang shows.

Under the bases each feature lies as a bar, in tracks packed by earliest start: pointed where the
feature ends on its strand, cut flat where a row cuts it, and a thin line across the gap a joined
feature leaves. A CDS's translation runs above its bar, each amino acid's three letters centred on
its codon's middle base, a stop in red. A feature's name goes inside its bar when it fits,
underneath when nothing else in its track lies there, and otherwise in a box above the ruler, in
the rows `labels.staircase` lays out, as the linear map's are. Primers and cut sites are labelled
there too, a cut site's enzymes one name a line. A row grows to hold what it draws, so nothing is
hidden.
"""

import dataclasses
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from functools import cache, lru_cache
from itertools import pairwise, repeat
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
#: The colour of the mark on a base a primer does not pair with.
MISMATCH = "#CC3311"

_INK = "#252525"
_RAIL = "#7f7f7f"
_LEADER = "#7f7f7f"
_CONNECTOR = "#666666"

# The bases: the room under the ruler and either side of the rail, how far the rail's ticks reach
# either side of it at a base, a fifth and a tenth, and the room before a row's last position.
_RULER_GAP = 2.0
_RAIL_GAP = 2.0
_TICKS = (1.5, 2.5, 3.5)
# Each tick after a row's first, on from the foot of the tick before, by how far the two reach.
_STEPS = {
    (before, reach): f"m{number(CELL)} {number(-before - reach)}v{number(2 * reach)}"
    for before in _TICKS
    for reach in _TICKS
}
_NUMBER_GAP = 14.0

# A feature's bar: its thickness, how far its point runs along it, the room between bars in a
# track, between tracks, under the bases and under a translation.
_BAR = 15.0
_POINT = 7.0
_BAR_GAP = 3.0
_TRACK_GAP = 3.0
_FEATURE_GAP = 4.0
_RESIDUE_GAP = 1.0

# A primer's arrow: its band's thickness, how far its head reaches beyond the band and along it,
# how far its tail's line lies from the band's middle and how thick it is, the room between tracks,
# and how far a mismatch's mark keeps within its cell. Then how thick a cut's lines are.
_PRIMER_BAND = 4.0
_PRIMER_REACH = 2.0
_PRIMER_HEAD = 6.0
_BEND = 6.0
_TAIL = 1.5
_PRIMER_GAP = 2.0
_MARK_INSET = 2.0
_CUT = 1.2

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
class Arrow:
    """A primer's binding site, or the part of it a row holds, as an arrow along the row.

    Parameters
    ----------
    item
        The primer at this binding site.
    start, end
        Where it starts and ends along the row.
    middle
        The height of the band's middle.
    head_start, head_end
        Whether a head's tip ends it at `start`, or at `end`: only at the primer's 3' end, never
        where a row cuts it.
    head
        How far along the band the head takes.
    """

    item: Item
    start: float
    end: float
    middle: float
    head_start: bool
    head_end: bool
    head: float

    @property
    def bounds(self) -> Box:
        """The box the arrow takes, reaching as far either side of its band as a head does."""
        reach = _PRIMER_BAND / 2 + _PRIMER_REACH
        return Box(self.start, self.middle - reach, self.end - self.start, 2 * reach)


@dataclass(frozen=True, slots=True)
class Tail:
    """What a row holds of a primer's 5' tail: a line, a base to a cell, bent away from the bases.

    Parameters
    ----------
    item
        The primer.
    points
        Along the line from its end nearer the binding site. Where the tail leaves the binding
        site in this row, the line starts on the arrow's middle and bends away within one cell.
    """

    item: Item
    points: tuple[Point, ...]


@dataclass(frozen=True, slots=True)
class Mismatch:
    """A mark, `box`, on a primer's arrow over a base the primer does not pair with."""

    item: Item
    box: Box


@dataclass(frozen=True, slots=True)
class Cut:
    """What a row holds of a cut site: lines through the strands, and along the rail between.

    Parameters
    ----------
    item
        The cut site.
    lines
        Each line's two ends: through the top strand where the site's enzymes cut it, along the
        rail to where each cuts the bottom strand, and through the bottom strand there, as much of
        each as lies in the row.
    """

    item: Item
    lines: tuple[tuple[Point, Point], ...]


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
        The box its ruler, strands, rail and primers take.
    strands
        The box its strands and the rail take.
    rail
        The height of its rail.
    bars, names, translations, arrows, tails, mismatches, cuts, labels
        What it draws of each item, in order along the row.
    shapes
        Everything in the row, in the order it is drawn.
    """

    start: int
    end: int
    extent: Box
    bases: Box
    strands: Box
    rail: float
    bars: tuple[Bar, ...]
    names: tuple[Name, ...]
    translations: tuple[Translation, ...]
    arrows: tuple[Arrow, ...]
    tails: tuple[Tail, ...]
    mismatches: tuple[Mismatch, ...]
    cuts: tuple[Cut, ...]
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
        One group of every row's shapes, grouped by row. What a page reads to find the base under
        a pointer is written on them: on the whole, the record's length and a base's width, and on
        each row, its `Row.start`, `Row.end` and `Row.rail`, and the top and bottom of its
        `Row.strands`.
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


@dataclass(frozen=True, slots=True)
class _Primed:
    """A primer's binding site and tail as a row holds them together: packed as one bar is.

    `pieces` is the one piece of the whole primer the row holds. `site` and `tail` are where along
    the row each lies, if it does; `bends` is whether the tail leaves the binding site in this row.
    """

    item: Item
    pieces: tuple[Piece, ...]
    start: float
    end: float
    site: tuple[float, float] | None
    head_start: bool
    head_end: bool
    tail: tuple[float, float] | None
    bends: bool


@dataclass(slots=True)
class _Found:
    """What lies in one row.

    Each feature's pieces and codons, each primer's pieces and mismatches, and each cut site's
    position, with whether the row labels it.
    """

    features: list[tuple[Item, list[Piece]]] = field(default_factory=list)
    codons: dict[int, list[tuple[int, Codon]]] = field(default_factory=dict)
    primers: list[tuple[Item, Piece]] = field(default_factory=list)
    mismatches: dict[int, list[int]] = field(default_factory=dict)
    cuts: list[tuple[Item, int, bool]] = field(default_factory=list)


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
        What the record draws; what lies in the stretch is drawn.
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
    found = [_Found() for _ in bounds]
    for item in items:
        match item.kind:
            case "feature":
                drawn = pieces(item, start, end, length, circular=circular)
                for index, cut in _by_row(drawn, start, bases_per_row).items():
                    found[index].features.append((item, cut))
                for codon in item.translation:
                    for at in _within(codon.middle, start, end, length):
                        index = (at - start) // bases_per_row
                        found[index].codons.setdefault(id(item), []).append((at, codon))
            case "primer":
                drawn = pieces(_whole(item), start, end, length, circular=circular)
                for index, cut in sorted(_by_row(drawn, start, bases_per_row).items()):
                    found[index].primers.extend((item, one) for one in cut)
                for position in item.mismatches:
                    for at in _within(position, start, end, length):
                        index = (at - start) // bases_per_row
                        found[index].mismatches.setdefault(id(item), []).append(at)
            case "cut_site":
                for piece in pieces(item, start, end, length, circular=circular):
                    for index, at, labelled in _cut_rows(
                        item, piece.start, start, end, length, circular, bases_per_row
                    ):
                        found[index].cuts.append((item, at, labelled))
    widest = max(SANS.width(_last(last, length), SMALL_SIZE) for _, last in bounds)
    right = bases_per_row * CELL + _NUMBER_GAP + widest
    rows: list[Row] = []
    top = 0.0
    before: dict[int, int] = {}
    for (first, last), here in zip(bounds, found, strict=True):
        row, before = _row(
            first,
            last,
            bases,
            here,
            before,
            top=top,
            right=right,
            width=bases_per_row * CELL,
            both_strands=both_strands,
        )
        rows.append(row)
        top = row.extent.y + row.extent.height + _ROW_GAP
    extent = _extent([row.extent for row in rows], margin=_EDGE)
    grouped = tuple(
        Group(
            row.shapes,
            classes=("row",),
            data={
                "start": str(row.start),
                "end": str(row.end),
                "top": number(row.strands.y),
                "rail": number(row.rail),
                "bottom": number(row.strands.y + row.strands.height),
            },
        )
        for row in rows
    )
    whole = Group(grouped, classes=("rows",), data={"length": str(length), "cell": number(CELL)})
    return SequenceView(extent, tuple(rows), (whole,))


def _within(position: int, start: int, end: int, length: int) -> list[int]:
    """Return where a position less than `length` lies in the stretch, counted as it counts them.

    It lies there once at most, since the stretch is at most one turn long.
    """
    return [at for at in (position, position + length) if start <= at < end]


def _whole(item: Item) -> Item:
    """Return a primer as one span from the 5' end of its tail to the 3' end of its binding site."""
    [(low, high)] = [(span.start, span.end) for span in item.spans]
    if item.strand == Strand.REVERSE:
        high += item.tail
    else:
        low -= item.tail
    return dataclasses.replace(item, spans=(Span(low, high, item.color),))


def _cut_rows(
    item: Item, at: int, start: int, end: int, length: int, circular: bool, bases_per_row: int
) -> list[tuple[int, int, bool]]:
    """Return each row a cut site is drawn in, where it cuts the top strand, and if it is labelled.

    `at` is where it cuts the top strand in the stretch. It is labelled in a row holding `at`: at
    the edge between two rows, the one its overhang lies in, the later one when it cuts blunt. It
    is drawn too in any other row its overhang reaches, a turn on or back across the origin of a
    circular record.
    """
    rows = len(range(start, end, bases_per_row))
    staggers = [cutter.stagger for cutter in item.cutters]
    low, high = min([0, *staggers]), max([0, *staggers])
    index, edge = divmod(at - start, bases_per_row)
    if index == rows or (index and not edge and low + high < 0):
        index -= 1
    drawn = {(index, at): True}
    for shift in (-length, 0, length) if circular else (0,):
        first, last = at + low + shift, at + high + shift
        for row in range(max(0, (first - start) // bases_per_row), rows):
            opens = start + row * bases_per_row
            if last <= opens:
                break
            if first < min(opens + bases_per_row, end):
                drawn.setdefault((row, at + shift), False)
    return [(row, position, labelled) for (row, position), labelled in drawn.items()]


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
    found: _Found,
    before: Mapping[int, int],
    *,
    top: float,
    right: float,
    width: float,
    both_strands: bool,
) -> tuple[Row, dict[int, int]]:
    """Lay out the row holding bases `first` to `last`, its top at `top`.

    `before` holds the track of each item the row before ran on into this one, by the item's id.
    Bars and primers starting together are packed in the order of those tracks, so the items keep
    their order down the rows. Returns the row, and the same for the row after.
    """

    def at(position: float) -> float:
        return (position - first) * CELL

    length = len(bases)
    runs = [
        run
        for item, cut in found.features
        for run in _runs(item, cut, at, found.codons.get(id(item), ()))
    ]
    tracks, after = _tracks(runs, before, last)
    primed = [_primed(item, piece, at, length) for item, piece in found.primers]
    forward = [run for run in primed if run.item.strand != Strand.REVERSE]
    reverse = [run for run in primed if run.item.strand == Strand.REVERSE]
    forward_tracks, forward_after = _tracks(forward, before, last)
    reverse_tracks, reverse_after = _tracks(reverse, before, last)
    after |= forward_after | reverse_after

    inside, under, boxed = _names(runs, tracks, width)
    down = _down(
        runs,
        tracks,
        under,
        both_strands,
        max(forward_tracks, default=-1) + 1,
        max(reverse_tracks, default=-1) + 1,
    )
    anchored = [(run.item, Point((run.start + run.end) / 2, 0.0)) for run in boxed]
    anchored += [
        (run.item, Point(sum(run.site) / 2, 0.0))
        for run in (*forward, *reverse)
        if run.site and SANS.drawn(run.item.label)
    ]
    anchored += [
        (item, Point(at(position), 0.0))
        for item, position, labelled in found.cuts
        if labelled and SANS.drawn(item.label)
    ]
    stairs = labels.staircase(
        [_anchored(item, anchor) for item, anchor in anchored], base=-_GAP, spacing=_SPACING
    )
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

    arrows: list[Arrow] = []
    tails: list[Tail] = []
    mismatches: list[Mismatch] = []
    for side, packed, middles, away in (
        (forward, forward_tracks, down.forwards, -1.0),
        (reverse, reverse_tracks, down.reverses, 1.0),
    ):
        for run, track in zip(side, packed, strict=True):
            middle = middles[track] + shift
            item = run.item
            _, shapes = groups.setdefault(id(item), (item, []))
            if run.tail:
                tails.append(Tail(item, _tail(run, middle, away)))
                points = " ".join(f"{number(x)} {number(y)}" for x, y in tails[-1].points)
                shapes.append(Path(f"M{points}", "none", item.color, _TAIL))
            if run.site:
                low, high = run.site
                head = min(_PRIMER_HEAD, high - low) if run.head_start or run.head_end else 0.0
                arrows.append(Arrow(item, low, high, middle, run.head_start, run.head_end, head))
                shapes.append(Path(_arrow(arrows[-1]), item.color, outline_color(item.color), 0.8))
                for position in found.mismatches.get(id(item), ()):
                    if low <= at(position) and at(position + 1) <= high:
                        reach = _PRIMER_BAND / 2 + _PRIMER_REACH
                        box = Box(
                            at(position) + _MARK_INSET,
                            middle - reach,
                            CELL - 2 * _MARK_INSET,
                            2 * reach,
                        )
                        mismatches.append(Mismatch(item, box))
                        shapes.append(Rect(box, MISMATCH, "none", 0.0))

    cuts: list[Cut] = []
    for item, position, _ in found.cuts:
        through, under = _cut(item, position, first, last, down, shift, both_strands)
        if not through + under:
            continue
        cuts.append(Cut(item, through + under))
        _, shapes = groups.setdefault(id(item), (item, []))
        shapes.extend(Line(x1, y1, x2, y2, _INK, _CUT) for (x1, y1), (x2, y2) in through)
        if under:
            # Grouped as the bottom strand is, for a page to hide with it.
            below = tuple(Line(x1, y1, x2, y2, _INK, _CUT) for (x1, y1), (x2, y2) in under)
            shapes.append(Group(below, classes=("bottom",)))

    placed = tuple(
        Label(
            item,
            Box(one.box.x, one.box.y + shift, one.box.width, one.box.height),
            (
                Point(one.leader[0].x, one.leader[0].y + shift),
                Point(one.leader[1].x, one.leader[1].y + shift),
            ),
        )
        for (item, _), one in zip(anchored, stairs, strict=True)
    )

    drawn, boxes = _bases_drawn(bases, first, last, down, shift, right, both_strands)
    boxes += [
        Box(0.0, shift, 0.0, down.bottom),
        *(Box(bar.start, bar.body.y, bar.end - bar.start, _BAR) for bar in bars),
        *(_letters_box(name.letters) for name in names),
        *(_letters_box(line) for one in translations for line in (one.letters, one.stops) if line),
        *(arrow.bounds for arrow in arrows),
        *(_line_box(tail.points) for tail in tails),
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
        Box(0.0, down.strand + shift, at(last), down.strands - down.strand),
        down.rail + shift,
        tuple(bars),
        tuple(names),
        tuple(translations),
        tuple(arrows),
        tuple(tails),
        tuple(mismatches),
        tuple(cuts),
        placed,
        shapes,
    )
    return row, after


def _tracks[R: (_Run, _Primed)](
    runs: list[R], before: Mapping[int, int], last: int
) -> tuple[tuple[int, ...], dict[int, int]]:
    """Sort `runs` as they are packed, and return each one's track and where items run on.

    Runs starting together are packed in the order of their items' tracks in `before`. What comes
    back with the tracks is the track of each item running on past `last`, by the item's id.
    """
    runs.sort(key=lambda run: (run.start, before.get(id(run.item), len(before))))
    tracks = labels.pack([(run.start, run.end) for run in runs], spacing=_BAR_GAP)
    after = {
        id(run.item): track
        for run, track in zip(runs, tracks, strict=True)
        if run.pieces[-1].end == last and not run.pieces[-1].ends
    }
    return tracks, after


class _Down(NamedTuple):
    """How far down a row, from the top of its ruler, each part of it lies.

    Parameters
    ----------
    forwards
        Each forward primer track's middle, the first nearest the top strand.
    strand, rail, complement, strands
        The top strand's top, the rail, the bottom strand's top, and where the strands end.
    reverses
        Each reverse primer track's middle, the first nearest the bottom strand.
    bases
        Where the bases, and the primers beside them, end.
    residues, middles
        Each feature track's amino acids' middle, where it has any, and its bars' middle.
    bottom
        Where the row ends.
    """

    forwards: dict[int, float]
    strand: float
    rail: float
    complement: float
    strands: float
    reverses: dict[int, float]
    bases: float
    residues: dict[int, float]
    middles: dict[int, float]
    bottom: float


def _down(
    runs: Sequence[_Run],
    tracks: Sequence[int],
    under: dict[int, float],
    both_strands: bool,
    forward: int,
    reverse: int,
) -> _Down:
    """Return how far down a row each part lies, with `forward` and `reverse` primer tracks.

    A feature track holding a translation, or a name underneath, is as much taller.
    """
    # A primer track reaches as far as its head towards the bases, and as its tail away from them.
    near = _PRIMER_BAND / 2 + _PRIMER_REACH
    away = max(_BEND + _TAIL / 2, near)
    y = _height(SANS, SMALL_SIZE) + _RULER_GAP
    forwards: dict[int, float] = {}
    for track in reversed(range(forward)):
        forwards[track] = y + away
        y += away + near + _PRIMER_GAP
    strand = y
    rail = strand + _height(MONO, BASE_SIZE) + _RAIL_GAP + _TICKS[-1]
    complement = rail + _TICKS[-1] + _RAIL_GAP
    y = strands = complement + _height(MONO, BASE_SIZE) if both_strands else rail + _TICKS[-1]
    reverses: dict[int, float] = {}
    for track in range(reverse):
        reverses[track] = y + _PRIMER_GAP + near
        y += _PRIMER_GAP + near + away
    bases = y
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
    return _Down(forwards, strand, rail, complement, strands, reverses, bases, residues, middles, y)


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


def _primed(item: Item, piece: Piece, at: Callable[[float], float], length: int) -> _Primed:
    """Return what of a primer a piece of its whole span holds: its binding site, its tail, both."""
    [(low, high)] = [(span.start, span.end) for span in item.spans]
    reverse = item.strand == Strand.REVERSE
    whole = (low, high + item.tail) if reverse else (low - item.tail, high)
    # The piece lies on the primer as it lies once in the stretch: here, or a turn on or back.
    shift = next(
        (
            shift
            for shift in (0, -length, length)
            if whole[0] + shift <= piece.start and piece.end <= whole[1] + shift
        ),
        0,
    )
    low, high = low + shift, high + shift
    first, last = max(piece.start, low), min(piece.end, high)
    site = (at(first), at(last)) if first < last else None
    if reverse:
        tail, bends = (max(piece.start, high), piece.end), piece.start <= high
    else:
        tail, bends = (piece.start, min(piece.end, low)), low <= piece.end
    return _Primed(
        item,
        (piece,),
        at(piece.start),
        at(piece.end),
        site,
        reverse and first == low,
        not reverse and last == high,
        (at(tail[0]), at(tail[1])) if tail[0] < tail[1] else None,
        bends,
    )


def _tail(run: _Primed, middle: float, away: float) -> tuple[Point, ...]:
    """Return a tail's line, from its end nearer the binding site, `away` up (-1) or down (1)."""
    assert run.tail is not None
    low, high = run.tail
    level = middle + away * _BEND
    near, far, step = (low, high, CELL) if away > 0 else (high, low, -CELL)
    if not run.bends:
        return (Point(near, level), Point(far, level))
    if abs(far - near) <= CELL:
        return (Point(near, middle), Point(far, level))
    return (Point(near, middle), Point(near + step, level), Point(far, level))


def _cut(
    item: Item, position: int, first: int, last: int, down: _Down, shift: float, both_strands: bool
) -> tuple[tuple[tuple[Point, Point], ...], tuple[tuple[Point, Point], ...]]:
    """Return the lines a cut site draws in a row, as `Cut.lines` gives them, in two parts.

    The first goes through the top strand, and the second along the rail and through the bottom
    strand.
    """
    width = (last - first) * CELL
    x = (position - first) * CELL
    rail = down.rail + shift
    through = []
    if first <= position <= last:
        through.append((Point(x, down.strand + shift), Point(x, rail)))
    if not both_strands:
        return tuple(through), ()
    under = []
    bottom = down.complement + _height(MONO, BASE_SIZE) + shift
    for stagger in sorted({cutter.stagger for cutter in item.cutters} or {0}):
        cut = (position + stagger - first) * CELL
        low, high = max(min(x, cut), 0.0), min(max(x, cut), width)
        if low < high:
            under.append((Point(low, rail), Point(high, rail)))
        if 0.0 <= cut <= width:
            under.append((Point(cut, rail), Point(cut, bottom)))
    return tuple(through), tuple(under)


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
    count = last - first
    shown = first % length + 1
    reach, steps, tens = _marks(shown % 10, count, min(count, length - shown + 1))
    numbers = [((index + 0.5) * CELL, str((first + index) % length + 1)) for index in tens]
    ticks = f"M{number(0.5 * CELL)} {number(rail - reach)}v{number(2 * reach)}{steps}"
    shapes: list[Shape] = [
        Line(0.0, rail, count * CELL, rail, _RAIL, 1.6),
        Path(ticks, "none", _RAIL, 1.0),
    ]
    if not numbers:
        return Group(tuple(shapes), classes=("ruler",)), []
    text = _centred(numbers, SANS, SMALL_SIZE, top + _height(SANS, SMALL_SIZE) / 2, _INK)
    return Group((text, *shapes), classes=("ruler",)), [_letters_box(text)]


@lru_cache(maxsize=1 << 8)
def _marks(phase: int, count: int, before: int) -> tuple[float, str, tuple[int, ...]]:
    """Return how far a row's first tick reaches, the ticks on from it, and which bases are tens.

    The row's first base is numbered `phase` modulo ten, and the numbering starts again at 1 after
    `before` bases, across the origin; rows alike in these share them.
    """
    shown = [*range(phase, phase + before), *range(1, count - before + 1)]
    reaches = [_TICKS[(one % 5 == 0) + (one % 10 == 0)] for one in shown]
    steps = "".join(map(_STEPS.__getitem__, pairwise(reaches)))
    return reaches[0], steps, tuple(index for index, one in enumerate(shown) if one % 10 == 0)


def _strand(text: str, top: float) -> Letters:
    """Return a strand's bases set one to a cell, its top at `top`."""
    baseline = _baseline(MONO, BASE_SIZE, top + _height(MONO, BASE_SIZE) / 2)
    # Made as the plain tuples places are, without a Python call for each base.
    places = zip(_cells(len(text)), repeat(baseline), repeat(0.0))
    return Letters(tuple(map(tuple.__new__, repeat(Place), places)), text, MONO, BASE_SIZE, _INK)


@cache
def _cells(count: int) -> tuple[float, ...]:
    """Return where each of `count` bases starts along a row, centred in its cell."""
    offset = (CELL - MONO.width("A", BASE_SIZE)) / 2
    return tuple(index * CELL + offset for index in range(count))


def _bases(bases: str, first: int, last: int) -> str:
    """Return the top strand from `first` to `last`, across the origin as needed."""
    start = first % len(bases)
    end = start + last - first
    return bases[start:end] if end <= len(bases) else bases[start:] + bases[: end - len(bases)]


def _last(last: int, length: int) -> str:
    """Return the position of the base before `last`, as a person reads it."""
    return str((last - 1) % length + 1)


def _lines(item: Item) -> tuple[tuple[tuple[str, bool], ...], ...]:
    """Return a label's lines, each in stretches of text marked bold or not.

    A cut site's enzymes go one a line, each bold where it cuts the record once.
    """
    if item.kind == "cut_site" and item.cutters:
        return tuple(((cutter.name, cutter.unique),) for cutter in item.cutters)
    return (item.runs,)


def _anchored(item: Item, anchor: Point) -> labels.Anchored:
    lines = _lines(item)
    width = max(sum(_font(bold).width(text, LABEL_SIZE) for text, bold in line) for line in lines)
    height = len(lines) * _height(SANS, LABEL_SIZE)
    return labels.Anchored(width + 2 * _PADDING[0], height + 2 * _PADDING[1], anchor)


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
    height = _height(SANS, LABEL_SIZE)
    for index, line in enumerate(_lines(item)):
        x = box.x + _PADDING[0]
        baseline = _baseline(SANS, LABEL_SIZE, box.y + _PADDING[1] + (index + 0.5) * height)
        for text, bold in line:
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


@cache
def _height(font: Font, size: float) -> float:
    return (font.ascender - font.descender) / font.units_per_em * size


def _baseline(font: Font, size: float, middle: float) -> float:
    """Return the baseline that centres a line of text on the height `middle`."""
    return middle + _drop(font, size)


@cache
def _drop(font: Font, size: float) -> float:
    return (font.ascender + font.descender) / 2 / font.units_per_em * size


def _text_box(text: Text) -> Box:
    font, size = text.font, text.size
    top = text.y - font.ascender / font.units_per_em * size
    return Box(text.x, top, font.width(text.text, size), _height(font, size))


def _letters_box(letters: Letters) -> Box:
    """Return the box a line of letters takes, from its first letter's start to its last's end."""
    font, size = letters.font, letters.size
    # A lone letter's width is its advance.
    last = font.width(font.drawn(letters.text)[-1], size)
    left, baseline = letters.places[0].x, letters.places[0].y
    top = baseline - font.ascender / font.units_per_em * size
    return Box(left, top, letters.places[-1].x + last - left, _height(font, size))


def _line_box(points: Sequence[Point]) -> Box:
    """Return the box a line through `points` takes, its thickness included."""
    xs, ys = [point.x for point in points], [point.y for point in points]
    half = _TAIL / 2
    return Box(min(xs), min(ys) - half, max(xs) - min(xs), max(ys) - min(ys) + 2 * half)


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


def _arrow(arrow: Arrow) -> str:
    """Return a primer's arrow as path data: its band, with its head where it has one."""
    half, reach, y = _PRIMER_BAND / 2, _PRIMER_BAND / 2 + _PRIMER_REACH, arrow.middle
    first = arrow.start + arrow.head * arrow.head_start
    last = arrow.end - arrow.head * arrow.head_end
    path = f"M{number(first)} {number(y - half)}"
    if arrow.head_end:
        path += f"H{number(last)}V{number(y - reach)}L{number(arrow.end)} {number(y)}"
        path += f"L{number(last)} {number(y + reach)}V{number(y + half)}"
    else:
        path += f"H{number(last)}V{number(y + half)}"
    path += f"H{number(first)}"
    if arrow.head_start:
        path += f"V{number(y + reach)}L{number(arrow.start)} {number(y)}"
        path += f"L{number(first)} {number(y - reach)}"
    return path + "Z"

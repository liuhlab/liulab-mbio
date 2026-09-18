"""Draw a record as a map, and write the drawing in the format its file name asks for."""

import dataclasses
import math
import os
from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import cast

from liulab_mbio.io import read_record
from liulab_mbio.plot import circular, convert, layers, page, svg
from liulab_mbio.plot import linear as line
from liulab_mbio.plot import sequence_view as view
from liulab_mbio.plot.labels import Box
from liulab_mbio.sequence import SequenceRecord

#: A region of a record: a feature's name, or a 0-based, half-open span ending past the record's
#: length across the origin.
type Region = str | tuple[int, int]

#: How many times its width a PDF's page of sequence view rows is tall: the proportions of A4, so
#: each page prints filling a sheet.
_PAGE = math.sqrt(2)

#: What a page's switch for each kind of item says, in the order the switches stand.
_KINDS: dict[layers.Kind, str] = {
    "feature": "Features",
    "primer": "Primers",
    "cut_site": "Cut sites",
}


@dataclass(frozen=True, slots=True)
class Switches:
    """Which items show: all a PNG or a PDF draws, and what a page shows first.

    Parameters
    ----------
    kinds
        The kinds of item switched on.
    types_off
        The feature types switched off.
    """

    kinds: frozenset[layers.Kind]
    types_off: frozenset[str]

    def show(self, kind: str, feature_type: str) -> bool:
        """Return whether an item of `kind` shows, and a feature also by its type."""
        return kind in self.kinds and (kind != "feature" or feature_type not in self.types_off)


@dataclass(frozen=True, slots=True, eq=False)
class Drawing:
    """A record to lay out as a map, to write once per format without laying it out again.

    A PNG or a PDF draws only the items that show, and the sequence view only when it is switched
    on. A page carries every item the record draws, laid out together, and the sequence view of
    any stretch up to `sequence_view.LIMIT` bases with both its strands; it shows first what is
    switched on, and switches the rest in place. Each distinct layout runs when it is first
    needed, and is kept.

    Parameters
    ----------
    record
        What is drawn.
    switches
        Which items show.
    enzymes
        The names of the enzymes whose every cut site is drawn; the shipped unique cutters when
        ``None``.
    span
        The stretch drawn, 0-based and half-open, ending past the record's length across the
        origin; the whole record when ``None``.
    linear
        Whether a circular record drawn whole is opened as a line.
    with_sequence_view
        Whether the sequence view is switched on.
    bases_per_row
        How many bases each row of the sequence view holds.
    both_strands
        Whether the sequence view's bottom strand, under the top one, is switched on.
    """

    record: SequenceRecord
    switches: Switches
    enzymes: str | tuple[str, ...] | None
    span: tuple[int, int] | None
    linear: bool
    with_sequence_view: bool
    bases_per_row: int
    both_strands: bool
    _kept: dict[tuple[object, ...], object] = field(default_factory=dict, init=False, repr=False)

    @property
    def layout(self) -> circular.CircularMap | line.LinearMap:
        """Where everything that shows went on the map, as a PNG or a PDF draws it."""
        shown = _shown(self)
        return (
            _line(self, shown) if self.linear or not _whole_circle(self) else _circle(self, shown)
        )

    @property
    def sequence_view(self) -> view.SequenceView | None:
        """Where everything that shows went in the sequence view, when it is switched on."""
        if not self.with_sequence_view:
            return None
        return _rows(self, _shown(self), both_strands=self.both_strands)

    @property
    def hidden(self) -> tuple[layers.Item, ...]:
        """The item of each label the map left out for want of room, in the order they hid.

        Each gives its kind, its name and the spans where it lies. `layers.notice` words them as
        the map does. This is the map a PNG or a PDF draws; a page counts what its own layout hid.
        """
        return self.layout.hidden

    def write(self, path: str | os.PathLike[str], *, dpi: float = 300) -> Path:
        """Write the drawing to `path`, in the format its suffix names, and return the path.

        A ``.html`` page keeps its text as text, and puts the sequence view beside the map, or
        under it on a narrow page. A ``.png`` at `dpi` and a ``.pdf`` draw every letter as its
        outline, so they look the same on any machine; `dpi` counts for the PNG alone. The PNG is
        one image, the sequence view under the map. The PDF has the map on its first page, and the
        sequence view's rows on the pages after, as many to a page as fit whole.

        A page carries every item behind its switches, a circular record drawn whole both as a
        circle and as a line, and the sequence view of a stretch up to `sequence_view.LIMIT`
        bases. A PNG and a PDF draw only what is switched on.

        Raises
        ------
        ValueError
            If the suffix names no format this writes, or a PNG at `dpi` would be too large to
            draw: past `convert.PNG_SIDE` pixels a side, or `convert.PNG_PIXELS` in all.
        """
        out = Path(path)
        writer = _WRITERS.get(out.suffix.lower())
        if writer is None:
            formats = ", ".join(_WRITERS)
            raise ValueError(f"cannot write a map as {out.name!r}: the suffix must be {formats}")
        writer(self, out, dpi)
        return out


def draw_map(
    record: SequenceRecord | str | os.PathLike[str],
    *,
    region: Region | None = None,
    linear: bool = False,
    sequence_view: bool = False,
    features: bool = True,
    primers: bool = True,
    cut_sites: bool = True,
    enzymes: str | Iterable[str] | None = None,
    hide_types: Iterable[str] = (),
    source: bool = False,
    bases_per_row: int = 60,
    both_strands: bool = True,
) -> Drawing:
    """Lay out a record as a map of its features, primers and cut sites: a circle or a line.

    A circular record is drawn as a circle, and opened as a line when `linear`. A linear record,
    and a region of any record, is always drawn as a line, numbered as the record is. The sequence
    view draws the same stretch base by base, as the map numbers it.

    What is switched off is left out of a PNG or a PDF, and a page keeps it behind its switches,
    shown first as asked: the layers, the feature types, the sequence view of a stretch up to
    `sequence_view.LIMIT` bases, and its bottom strand. Nothing is laid out until it is needed.

    Parameters
    ----------
    record
        A record, or a ``.dna``, GenBank or FASTA file to read one from.
    region
        The stretch to draw: the name of a feature, whose first segment's start to its last's end
        is drawn, or a ``(start, end)`` span, 0-based and half-open, ending past the record's
        length across the origin. A name matches whatever its case, and the first feature in the
        record's order answers when several do. The whole record when ``None``.
    linear
        Whether a circular record drawn whole is opened as a line.
    sequence_view
        Whether the sequence view is switched on, at most `sequence_view.LIMIT` bases.
    features, primers, cut_sites
        Whether each is switched on.
    enzymes
        The names of the enzymes whose every cut site is drawn; the shipped unique cutters when
        ``None``.
    hide_types
        The feature types switched off.
    source
        Whether a `source` feature is switched on.
    bases_per_row
        How many bases each row of the sequence view holds.
    both_strands
        Whether the sequence view's bottom strand, under the top one, is switched on.

    Raises
    ------
    KeyError
        If no shipped enzyme answers to a name in `enzymes`.
    ValueError
        If the file cannot be read as a record, the record has no bases, `region` names no
        feature or lies off the record, `bases_per_row` is less than 1, or a sequence view would
        draw more than `sequence_view.LIMIT` bases.

    Examples
    --------
    >>> draw_map("pUC19.dna").write("pUC19.html")  # doctest: +SKIP
    PosixPath('pUC19.html')
    """
    record = record if isinstance(record, SequenceRecord) else read_record(record)
    if not len(record):
        raise ValueError(f"record {record.name!r} has no bases to draw")
    if bases_per_row < 1:
        raise ValueError(f"a row of the sequence view holds at least 1 base, not {bases_per_row}")
    span = None if region is None else _span(record, region)
    start, end = span if span is not None else (0, len(record))
    if sequence_view and end - start > view.LIMIT:
        raise ValueError(
            f"a sequence view draws at most {view.LIMIT:,} bases, and {record.name!r} would "
            f"draw {end - start:,}: name a region to draw it"
        )
    named = enzymes if enzymes is None or isinstance(enzymes, str) else tuple(enzymes)
    # Reads the names alone, so an unknown one raises here rather than when a layout needs it.
    layers.items(record, features=False, primers=False, cut_sites=False, enzymes=named)
    kinds: dict[layers.Kind, bool] = {"feature": features, "primer": primers, "cut_site": cut_sites}
    switches = Switches(
        frozenset(kind for kind, on in kinds.items() if on),
        frozenset(hide_types) | (frozenset() if source else {"source"}),
    )
    return Drawing(
        record,
        switches,
        named,
        span,
        linear,
        sequence_view,
        bases_per_row,
        both_strands,
    )


def _span(record: SequenceRecord, region: Region) -> tuple[int, int]:
    """Return the span `region` names in `record`.

    Raises
    ------
    ValueError
        If `region` names no feature of the record, or is a span that does not lie on it.
    """
    length = len(record)
    if isinstance(region, str):
        feature = next(
            (one for one in record.features if one.name.casefold() == region.casefold()), None
        )
        if feature is None:
            raise ValueError(f"record {record.name!r} has no feature called {region!r} to draw")
        counted = layers.unwrapped(feature.segments, length)
        return counted[0][0], min(counted[-1][1], counted[0][0] + length)
    start, end = region
    if record.topology == "circular":
        if not 0 <= start < length or not start < end <= start + length:
            raise ValueError(
                f"region {region} does not lie on the circular record {record.name!r} of "
                f"{length} bases: it starts inside the record, and across the origin ends past "
                "its length, at most one turn on"
            )
    elif not 0 <= start < end <= length:
        raise ValueError(
            f"region {region} does not lie on the linear record {record.name!r} of {length} bases"
        )
    return start, end


def _whole_circle(drawing: Drawing) -> bool:
    return drawing.record.topology == "circular" and drawing.span is None


def _carries_sequence_view(drawing: Drawing) -> bool:
    """Whether a page carries the sequence view: of a stretch up to `sequence_view.LIMIT` bases."""
    start, end = drawing.span or (0, len(drawing.record))
    return end - start <= view.LIMIT


def _kept[T](drawing: Drawing, key: tuple[object, ...], make: Callable[[], T]) -> T:
    """Return what `make` gives, made the first time `key` is asked for and kept."""
    if key not in drawing._kept:
        drawing._kept[key] = make()
    return cast(T, drawing._kept[key])


def _uncut(drawing: Drawing) -> tuple[layers.Item, ...]:
    """Every feature and primer the record draws, whether it shows or not."""
    return _kept(
        drawing,
        ("uncut",),
        lambda: layers.items(
            drawing.record,
            cut_sites=False,
            source=True,
            translations=_carries_sequence_view(drawing),
        ),
    )


def _cuts(drawing: Drawing) -> tuple[layers.Item, ...]:
    """Every cut site the record draws, searched for only when a layout needs them."""
    return _kept(
        drawing,
        ("cuts",),
        lambda: layers.items(
            drawing.record, features=False, primers=False, enzymes=drawing.enzymes
        ),
    )


def _everything(drawing: Drawing) -> tuple[layers.Item, ...]:
    return _kept(drawing, ("everything",), lambda: (*_uncut(drawing), *_cuts(drawing)))


def _shown(drawing: Drawing) -> tuple[layers.Item, ...]:
    switches = drawing.switches
    found = _everything(drawing) if "cut_site" in switches.kinds else _uncut(drawing)
    return tuple(item for item in found if switches.show(item.kind, item.type))


def _circle(drawing: Drawing, items: tuple[layers.Item, ...]) -> circular.CircularMap:
    record = drawing.record
    return _kept(
        drawing,
        ("circle", items),
        lambda: circular.layout(items, name=record.name, length=len(record)),
    )


def _line(drawing: Drawing, items: tuple[layers.Item, ...]) -> line.LinearMap:
    record = drawing.record
    return _kept(
        drawing,
        ("line", items),
        lambda: line.layout(
            items,
            name=record.name,
            length=len(record),
            circular=record.topology == "circular",
            span=drawing.span,
        ),
    )


def _rows(
    drawing: Drawing, items: tuple[layers.Item, ...], *, both_strands: bool
) -> view.SequenceView:
    record = drawing.record
    return _kept(
        drawing,
        ("rows", items, both_strands),
        lambda: view.layout(
            items,
            bases=record.sequence,
            circular=record.topology == "circular",
            span=drawing.span,
            bases_per_row=drawing.bases_per_row,
            both_strands=both_strands,
        ),
    )


def _first_shown(
    shapes: Sequence[svg.Shape], switches: Switches, hidden: Sequence[layers.Item] = ()
) -> tuple[svg.Shape, ...]:
    """Return shapes as a page first shows them.

    Each group drawing an item that does not show is marked ``off``, and the notice of what a map
    hid counts only the `hidden` items that show.
    """
    return tuple(_first_shown_one(shape, switches, hidden) for shape in shapes)


def _first_shown_one(
    shape: svg.Shape, switches: Switches, hidden: Sequence[layers.Item]
) -> svg.Shape:
    if not isinstance(shape, svg.Group):
        return shape
    if "notice" in shape.classes:
        return _notice(shape, [item for item in hidden if switches.show(item.kind, item.type)])
    data = shape.data
    off = "kind" in data and not switches.show(data["kind"], data.get("type", ""))
    return dataclasses.replace(
        shape,
        shapes=_first_shown(shape.shapes, switches, hidden),
        classes=(*shape.classes, "off") if off else shape.classes,
    )


def _notice(notice: svg.Group, hidden: Sequence[layers.Item]) -> svg.Group:
    """Return a map's notice saying what hid among `hidden`, its right edge where it was.

    It is marked ``off`` when none of them hid.
    """
    said = layers.notice(hidden)
    if not said:
        return dataclasses.replace(notice, classes=(*notice.classes, "off"))
    lines = []
    for shape in notice.shapes:
        if isinstance(shape, svg.Text):
            right = shape.x + shape.font.width(shape.text, shape.size)
            shape = dataclasses.replace(
                shape, x=right - shape.font.width(said, shape.size), text=said
            )
        lines.append(shape)
    return dataclasses.replace(notice, shapes=tuple(lines))


def _switches(switches: Switches, drawn: line.LinearMap) -> list[page.Switch]:
    """Return a page's switch for each kind of item `drawn` holds, then for each feature type."""
    items = [
        *(arrow.item for arrow in drawn.arrows),
        *(label.item for label in drawn.labels),
        *drawn.hidden,
    ]
    kinds = {item.kind for item in items}
    types = sorted(
        {item.type for item in items if item.kind == "feature"},
        key=lambda one: (one.casefold(), one),
    )
    return [
        *(
            page.Switch("kind", kind, words, kind in switches.kinds)
            for kind, words in _KINDS.items()
            if kind in kinds
        ),
        *(page.Switch("type", one, one, one not in switches.types_off) for one in types),
    ]


def _html(drawing: Drawing, path: Path, dpi: float) -> None:
    everything = _everything(drawing)
    layouts: dict[str, circular.CircularMap | line.LinearMap] = {}
    if _whole_circle(drawing):
        layouts["circle"] = _circle(drawing, everything)
    layouts["line"] = _line(drawing, everything)
    switches = drawing.switches
    maps = {
        shape: svg.document(_first_shown(one.shapes, switches, one.hidden), one.extent)
        for shape, one in layouts.items()
    }
    beside = None
    if _carries_sequence_view(drawing):
        # Both strands, so hiding the bottom one moves nothing.
        rows = _rows(drawing, everything, both_strands=True)
        beside = svg.document(_first_shown(rows.shapes, switches), rows.extent)
    shown = "circle" if "circle" in layouts and not drawing.linear else "line"
    html = page.render(
        maps,
        title=drawing.record.name or path.stem,
        shown=shown,
        zooms=[shape for shape in layouts if shape == "circle"],
        switches=_switches(switches, layouts["line"]),
        sequence_view=beside,
        sequence_shown=drawing.with_sequence_view,
        both_strands=drawing.both_strands,
    )
    path.write_text(html, encoding="utf-8")


def _png(drawing: Drawing, path: Path, dpi: float) -> None:
    extent, shapes = _stacked(drawing)
    width, height = extent.width, extent.height
    least = convert.POINTS_PER_INCH / min(width, height)
    most = convert.POINTS_PER_INCH * min(
        convert.PNG_SIDE / max(width, height), math.sqrt(convert.PNG_PIXELS / (width * height))
    )
    if not least <= dpi <= most:
        drawn = f"can be drawn at {least:.3g} to {math.floor(most):,} dpi"
        if drawing.sequence_view is None:
            reason = f"a PNG of this map {drawn}, and a PDF at any size"
        else:
            reason = (
                f"a PNG of this map and its sequence view {drawn}; name a region to draw fewer "
                "bases, or write a PDF, which draws at any size"
            )
        raise ValueError(f"cannot draw {path.name!r} at {dpi:g} dpi: {reason}")
    path.write_bytes(convert.png(_outlined(shapes, extent), dpi=dpi))


def _pdf(drawing: Drawing, path: Path, dpi: float) -> None:
    path.write_bytes(convert.pdf(_pages(drawing)))


def _stacked(drawing: Drawing) -> tuple[Box, tuple[svg.Shape, ...]]:
    """Return the canvas and shapes of the map with the sequence view, if drawn, centred under it."""
    top, rows = drawing.layout.extent, drawing.sequence_view
    if rows is None:
        return top, drawing.layout.shapes
    under = rows.extent
    width = max(top.width, under.width)
    left = top.x - (width - top.width) / 2
    moved = svg.Group(
        rows.shapes,
        x=left + (width - under.width) / 2 - under.x,
        y=top.y + top.height - under.y,
    )
    return Box(left, top.y, width, top.height + under.height), (*drawing.layout.shapes, moved)


def _pages(drawing: Drawing) -> Iterator[str]:
    """Yield each page of a PDF: the map, then the sequence view's row blocks, none split.

    A page of rows is as wide as the view and `_PAGE` times as tall, and holds as many blocks as
    fit; a block too tall for a page has a page of its own, as tall as it needs.
    """
    yield _outlined(drawing.layout.shapes, drawing.layout.extent)
    rows = drawing.sequence_view
    if rows is None:
        return
    margin = rows.rows[0].extent.y - rows.extent.y
    tall = rows.extent.width * _PAGE
    page: list[view.Row] = []
    for row in rows.rows:
        if page and _bottom(row) + margin - (page[0].extent.y - margin) > tall:
            yield _page(rows.extent, page, margin, tall)
            page = []
        page.append(row)
    yield _page(rows.extent, page, margin, tall)


def _page(extent: Box, page: Sequence[view.Row], margin: float, tall: float) -> str:
    """Return a page of rows as wide as `extent`, `margin` clear round them and `tall` at least."""
    top = page[0].extent.y - margin
    bottom = max(top + tall, _bottom(page[-1]) + margin)
    box = Box(extent.x, top, extent.width, bottom - top)
    return _outlined([shape for row in page for shape in row.shapes], box)


def _bottom(row: view.Row) -> float:
    return row.extent.y + row.extent.height


def _outlined(shapes: Iterable[svg.Shape], extent: Box) -> str:
    """Return `shapes` as a PNG or a PDF draws them: on white, every letter an outline."""
    paper = svg.Rect(extent, "#ffffff", "none", 0)
    return svg.document((paper, *shapes), extent, outlines=True)


#: Each suffix a drawing is written as, and what writes it at a dpi.
_WRITERS: dict[str, Callable[[Drawing, Path, float], None]] = {
    ".html": _html,
    ".png": _png,
    ".pdf": _pdf,
}

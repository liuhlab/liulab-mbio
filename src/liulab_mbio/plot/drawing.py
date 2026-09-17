"""Draw a record as a map, and write the drawing in the format its file name asks for."""

import math
import os
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

from liulab_mbio.io import read_record
from liulab_mbio.plot import circular, convert, layers, page, svg
from liulab_mbio.plot import linear as line
from liulab_mbio.plot import sequence_view as view
from liulab_mbio.sequence import SequenceRecord

#: A region of a record: a feature's name, or a 0-based, half-open span ending past the record's
#: length across the origin.
type Region = str | tuple[int, int]


@dataclass(frozen=True, slots=True)
class Drawing:
    """A record laid out as a map, to write once per format without laying it out again.

    Parameters
    ----------
    record
        What is drawn.
    layout
        Where everything on the map went.
    sequence_view
        Where everything in the sequence view went, when it was asked for.
    """

    record: SequenceRecord
    layout: circular.CircularMap | line.LinearMap
    sequence_view: view.SequenceView | None = None

    @property
    def hidden(self) -> tuple[layers.Item, ...]:
        """The item of each label the map left out for want of room, in the order they hid.

        Each gives its kind, its name and the spans where it lies. `layers.notice` words them as
        the map does.
        """
        return self.layout.hidden

    def write(self, path: str | os.PathLike[str], *, dpi: float = 300) -> Path:
        """Write the drawing to `path`, in the format its suffix names, and return the path.

        A ``.html`` page keeps its text as text, and shows the sequence view beside the map. A
        ``.png`` at `dpi` and a ``.pdf`` draw the map alone, every letter as its outline, so they
        look the same on any machine; `dpi` counts for the PNG alone.

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
        Whether the sequence view is drawn beside the map, at most `sequence_view.LIMIT` bases.
    features, primers, cut_sites
        Whether each is drawn.
    enzymes
        The names of the enzymes whose every cut site is drawn; the shipped unique cutters when
        ``None``.
    hide_types
        The feature types left off.
    source
        Whether a `source` feature is drawn.
    bases_per_row
        How many bases each row of the sequence view holds.
    both_strands
        Whether the sequence view draws the bottom strand under the top one.

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
    items = layers.items(
        record,
        features=features,
        primers=primers,
        cut_sites=cut_sites,
        enzymes=enzymes,
        hide_types=hide_types,
        source=source,
        translations=sequence_view,
    )
    circle = record.topology == "circular"
    if circle and span is None and not linear:
        layout = circular.layout(items, name=record.name, length=len(record))
    else:
        layout = line.layout(
            items, name=record.name, length=len(record), circular=circle, span=span
        )
    rows = None
    if sequence_view:
        rows = view.layout(
            items,
            bases=record.sequence,
            circular=circle,
            span=span,
            bases_per_row=bases_per_row,
            both_strands=both_strands,
        )
    return Drawing(record, layout, rows)


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


def _html(drawing: Drawing, path: Path, dpi: float) -> None:
    image = svg.document(drawing.layout.shapes, drawing.layout.extent)
    rows = drawing.sequence_view
    beside = None if rows is None else svg.document(rows.shapes, rows.extent)
    title = drawing.record.name or path.stem
    path.write_text(page.render(image, title=title, sequence_view=beside), encoding="utf-8")


def _png(drawing: Drawing, path: Path, dpi: float) -> None:
    extent = drawing.layout.extent
    width, height = extent.width, extent.height
    least = convert.POINTS_PER_INCH / min(width, height)
    most = convert.POINTS_PER_INCH * min(
        convert.PNG_SIDE / max(width, height), math.sqrt(convert.PNG_PIXELS / (width * height))
    )
    if not least <= dpi <= most:
        raise ValueError(
            f"cannot draw {path.name!r} at {dpi:g} dpi: a PNG of this map can be drawn at "
            f"{least:.3g} to {math.floor(most):,} dpi, and a PDF at any size"
        )
    path.write_bytes(convert.png(_outlined(drawing), dpi=dpi))


def _pdf(drawing: Drawing, path: Path, dpi: float) -> None:
    path.write_bytes(convert.pdf([_outlined(drawing)]))


def _outlined(drawing: Drawing) -> str:
    """Return the drawing as a PNG or a PDF draws it: on white, every letter an outline."""
    extent = drawing.layout.extent
    paper = svg.Rect(extent, "#ffffff", "none", 0)
    return svg.document((paper, *drawing.layout.shapes), extent, outlines=True)


#: Each suffix a drawing is written as, and what writes it at a dpi.
_WRITERS: dict[str, Callable[[Drawing, Path, float], None]] = {
    ".html": _html,
    ".png": _png,
    ".pdf": _pdf,
}

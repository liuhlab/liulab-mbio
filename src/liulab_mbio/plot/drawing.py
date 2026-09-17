"""Draw a record as a map, and write the drawing in the format its file name asks for."""

import math
import os
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

from liulab_mbio.io import read_record
from liulab_mbio.plot import circular, convert, layers, page, svg
from liulab_mbio.sequence import SequenceRecord


@dataclass(frozen=True, slots=True)
class Drawing:
    """A record laid out as a map, to write once per format without laying it out again.

    Parameters
    ----------
    record
        What is drawn.
    layout
        Where everything went.
    """

    record: SequenceRecord
    layout: circular.CircularMap

    def write(self, path: str | os.PathLike[str], *, dpi: float = 300) -> Path:
        """Write the drawing to `path`, in the format its suffix names, and return the path.

        A ``.html`` page keeps its text as text. A ``.png`` at `dpi` and a ``.pdf`` draw every
        letter as its outline, so they look the same on any machine; `dpi` counts for the PNG
        alone.

        Raises
        ------
        ValueError
            If the suffix names no format this writes, or a PNG at `dpi` would be too large to
            draw.
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
    features: bool = True,
    primers: bool = True,
    cut_sites: bool = True,
    enzymes: str | Iterable[str] | None = None,
    hide_types: Iterable[str] = (),
    source: bool = False,
) -> Drawing:
    """Lay out a record as a circular map of its features, primers and cut sites.

    Parameters
    ----------
    record
        A record, or a ``.dna``, GenBank or FASTA file to read one from.
    features, primers, cut_sites
        Whether each is drawn.
    enzymes
        The names of the enzymes whose every cut site is drawn; the shipped unique cutters when
        ``None``.
    hide_types
        The feature types left off.
    source
        Whether a `source` feature is drawn.

    Raises
    ------
    KeyError
        If no shipped enzyme answers to a name in `enzymes`.
    ValueError
        If the file cannot be read as a record, or the record has no bases.

    Examples
    --------
    >>> draw_map("pUC19.dna").write("pUC19.html")  # doctest: +SKIP
    PosixPath('pUC19.html')
    """
    record = record if isinstance(record, SequenceRecord) else read_record(record)
    if not len(record):
        raise ValueError(f"record {record.name!r} has no bases to draw")
    items = layers.items(
        record,
        features=features,
        primers=primers,
        cut_sites=cut_sites,
        enzymes=enzymes,
        hide_types=hide_types,
        source=source,
    )
    layout = circular.layout(items, name=record.name, length=len(record))
    return Drawing(record, layout)


def _html(drawing: Drawing, path: Path, dpi: float) -> None:
    image = svg.document(drawing.layout.shapes, drawing.layout.extent)
    path.write_text(page.render(image, title=drawing.record.name or path.stem), encoding="utf-8")


def _png(drawing: Drawing, path: Path, dpi: float) -> None:
    extent = drawing.layout.extent
    sides = (extent.width, extent.height)
    least = convert.POINTS_PER_INCH / min(sides)
    most = convert.PNG_SIDE * convert.POINTS_PER_INCH / max(sides)
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

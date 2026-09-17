"""Draw a record as a map, and write the drawing in the format its file name asks for."""

import os
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from liulab_mbio.io import read_record
from liulab_mbio.plot import circular, layers, page, svg
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

    def write(self, path: str | os.PathLike[str]) -> Path:
        """Write the drawing to `path`, in the format its suffix names, and return the path.

        Raises
        ------
        ValueError
            If the suffix names no format this writes: ``.html`` is the one.
        """
        out = Path(path)
        writer = _WRITERS.get(out.suffix.lower())
        if writer is None:
            formats = ", ".join(_WRITERS)
            raise ValueError(f"cannot write a map as {out.name!r}: the suffix must be {formats}")
        writer(self, out)
        return out


def draw_map(record: SequenceRecord | str | os.PathLike[str]) -> Drawing:
    """Lay out a record as a circular map of its features.

    Parameters
    ----------
    record
        A record, or a ``.dna``, GenBank or FASTA file to read one from.

    Raises
    ------
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
    layout = circular.layout(layers.items(record), name=record.name, length=len(record))
    return Drawing(record, layout)


def _html(drawing: Drawing, path: Path) -> None:
    image = svg.document(drawing.layout.shapes, drawing.layout.extent)
    path.write_text(page.render(image, title=drawing.record.name or path.stem), encoding="utf-8")


#: Each suffix a drawing is written as, and what writes it.
_WRITERS: dict[str, Callable[[Drawing, Path], None]] = {".html": _html}

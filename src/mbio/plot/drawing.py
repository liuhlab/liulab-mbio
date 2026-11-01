"""Draw a record as a map through molviz.

`draw_map` hands molviz the record and the options a person chose, and returns molviz's own
drawing, which writes an HTML page, an SVG, a PNG or a PDF.
"""

import os
from collections.abc import Iterable
from typing import Any

import molviz

from mbio.enzymes import get_enzyme
from mbio.io import read_record
from mbio.plot.converter import molviz_record
from mbio.sequence import SequenceRecord, counted_round

#: A region of a record: a feature's name, or a 0-based, half-open span ending past the record's
#: length across the origin.
type Region = str | tuple[int, int]


def draw_map(
    record: SequenceRecord | str | os.PathLike[str],
    *,
    region: Region | None = None,
    highlight: str | Iterable[str] = (),
    linear: bool = False,
    sequence_view: bool = False,
    features: bool = True,
    primers: bool = True,
    cut_sites: bool = True,
    enzymes: str | Iterable[str] | None = None,
    hide_types: Iterable[str] = (),
    source: bool = False,
    bases_per_row: int | None = None,
    both_strands: bool = True,
) -> molviz.Drawing:
    """Draw a record as a map of its features, primers and cut sites: a circle or a line.

    A circular record is drawn as a circle, and opened as a line when `linear`. A linear record,
    and a region of any record, is always drawn as a line, numbered as the record is. The sequence
    view draws the same stretch base by base. molviz draws it, and the drawing it returns writes
    each format its suffix names; a page keeps what is switched off behind its switches.

    Parameters
    ----------
    record
        A record, or a ``.dna``, GenBank or FASTA file to read one from.
    region
        The stretch to draw: the name of a feature, whose first segment's start to its last's end
        is drawn, or a ``(start, end)`` span, 0-based and half-open, ending past the record's
        length across the origin. A name matches whatever its case, and the first feature in the
        record's order answers when several do. The whole record when ``None``.
    highlight
        The name, or names, of what the map points at: a feature, a primer, or an enzyme, which
        lights every cut site naming it. Every other item dims, and nothing moves.
    linear
        Whether a circular record drawn whole is opened as a line.
    sequence_view
        Whether the sequence view is switched on, at most `molviz.SEQUENCE_VIEW_LIMIT` bases.
    features, primers, cut_sites
        Whether each is switched on.
    enzymes
        The names of the enzymes whose every cut site is drawn; the shipped unique cutters when
        ``None``. A commercial or isoschizomer name answers for its shipped enzyme.
    hide_types
        The feature types switched off.
    source
        Whether a `source` feature is switched on.
    bases_per_row
        How many bases each row of the sequence view holds; a page fits its rows to its width,
        and a file holds 60, when ``None``.
    both_strands
        Whether the sequence view's bottom strand, under the top one, is switched on.

    Raises
    ------
    KeyError
        If no shipped enzyme answers to a name in `enzymes`.
    ValueError
        If the file cannot be read as a record, `region` names no feature or lies off the record,
        or molviz refuses the record or an option: a record with no bases, a `highlight` naming
        nothing the record draws, a `bases_per_row` under 1, or a sequence view past its limit.

    Examples
    --------
    >>> draw_map("pUC19.dna").write("pUC19.html")  # doctest: +SKIP
    """
    record = record if isinstance(record, SequenceRecord) else read_record(record)
    off = {"features": features, "primers": primers, "cuts": cut_sites}
    options: dict[str, Any] = {
        "kinds_off": [kind for kind, on in off.items() if not on],
        "types_off": sorted({*hide_types, *(() if source else ("source",))}),
        "bottom_strand": both_strands,
    }
    if region is not None:
        start, end = _span(record, region)
        options["span"] = {"start": start, "end": end}
    if linear:
        options["map"] = "line"
    if sequence_view:
        options["sequence_view"] = True
    if enzymes is not None:
        named = (enzymes,) if isinstance(enzymes, str) else enzymes
        options["enzymes"] = list(dict.fromkeys(get_enzyme(name).name for name in named))
    if lit := [highlight] if isinstance(highlight, str) else list(highlight):
        options["highlight"] = lit
    if bases_per_row is not None:
        options["bases_per_row"] = bases_per_row
    return molviz.draw(molviz_record(record), **options)


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
        starts = counted_round((one.start for one in feature.segments), length)
        last = feature.segments[-1]
        return starts[0], min(starts[-1] + last.end - last.start, starts[0] + length)
    start, end = region
    if start < end and record.fits(start, end):
        return start, end
    if record.topology == "circular":
        raise ValueError(
            f"region {region} does not lie on the circular record {record.name!r} of "
            f"{length} bases: it starts inside the record, and across the origin ends past "
            "its length, at most one turn on"
        )
    raise ValueError(
        f"region {region} does not lie on the linear record {record.name!r} of {length} bases"
    )

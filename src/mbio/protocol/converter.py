"""A protocol's figures, gels and plates in molviz's formats, as plain data for its draw calls.

A figure becomes its records, each read relative to the protocol's folder and converted by
`mbio.plot.converter.molviz_record`, and the options ``molviz.draw_figure`` takes, as
`mbio.plot.converter.molviz_options` writes them. A gel becomes molviz's gel and a plate molviz's plate, a plate's
labels and seating both drawn, the seating where both name a well. A transfer that is a stamp
becomes its two plates, every well it touches holding the transfer's title. The formats'
reference: <https://liuhlab.github.io/liulab-molviz/reference/record/>.
"""

import os
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from mbio.io import read_record
from mbio.plot.converter import molviz_options, molviz_record
from mbio.protocol.model import Figure, Gel, Plate, Transfer


def molviz_figure(figure: Figure, base: str | os.PathLike[str]) -> dict[str, Any]:
    """Return `figure` as ``{"records": [...], "options": {...}}`` for ``molviz.draw_figure``.

    Each record is read from its path relative to `base`, the folder the protocol was read from,
    and the options are `mbio.plot.converter.molviz_options`.

    Raises
    ------
    FileNotFoundError
        If a record is not there.
    KeyError
        If no shipped enzyme answers to a name the figure gives.
    """
    options = molviz_options(
        span=figure.span,
        linear=figure.linear,
        sequence_view=figure.sequence_view,
        enzymes=figure.enzymes,
        highlight=figure.highlight,
    )
    records = [molviz_record(read_record(Path(base) / path)) for path in figure.records]
    return {"records": records, "options": options}


def molviz_gel(gel: Gel) -> dict[str, Any]:
    """Return `gel` in molviz's gel format, its title left out where it has none.

    Examples
    --------
    >>> from mbio.protocol.model import Ladder, Lane
    >>> molviz_gel(Gel(Ladder("1 kb", (1000,)), (Lane("empty"),)))
    {'ladder': {'name': '1 kb', 'bands': [1000]}, 'lanes': [{'name': 'empty', 'bands': []}]}
    """
    out: dict[str, Any] = {"title": gel.title} if gel.title else {}
    out["ladder"] = {"name": gel.ladder.name, "bands": list(gel.ladder.bands_bp)}
    out["lanes"] = [{"name": lane.label, "bands": list(lane.bands_bp)} for lane in gel.lanes]
    return out


def molviz_plate(plate: Plate) -> dict[str, Any]:
    """Return `plate` in molviz's plate format, titled with its name, its wells in reading order."""
    return _plate(plate, plate.name, {**plate.labels, **plate.seating})


def molviz_stamp(transfer: Transfer, plates: Sequence[Plate]) -> list[dict[str, Any]]:
    """Return the plates a stamp takes from and fills, each touched well holding its title.

    They are titled ``From`` and ``Into`` each plate's name. A transfer that repeats no one
    pattern, names a plate `plates` does not declare, or moves through a well its plate has not
    got is no stamp a drawing can show, and returns no plate.
    """
    declared = {plate.name: plate for plate in plates}
    source = declared.get(transfer.plates[0])
    destination = declared.get(transfer.plates[-1])
    if transfer.stamp is None or source is None or destination is None:
        return []
    taken = {move.source.well for move in transfer.moves}
    filled = {move.destination.well for move in transfer.moves}
    if not (taken <= set(source.well_names) and filled <= set(destination.well_names)):
        return []
    return [
        _plate(source, f"From {source.name}", dict.fromkeys(taken, transfer.title)),
        _plate(destination, f"Into {destination.name}", dict.fromkeys(filled, transfer.title)),
    ]


def _plate(plate: Plate, title: str, wells: Mapping[str, str]) -> dict[str, Any]:
    out: dict[str, Any] = {"title": title, "rows": plate.rows, "columns": plate.columns}
    if wells:
        out["wells"] = {name: wells[name] for name in plate.well_names if name in wells}
    return out

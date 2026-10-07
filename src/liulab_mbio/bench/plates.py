"""Plates, where each thing sits in one, and the moves between them.

Format is one parameter — the well count — and `liulab_mbio.protocol.model.FORMATS` holds the
rows and columns of each, 1536 included. A plate says where a thing sits and nothing else; what
differs between the reactions stays on the reaction, so the two can never disagree.

Two functions cover every move the bench makes, because they are all one shape. `compact`
re-lays a set of wells out densely into a destination, which is a plain or an acoustic transfer
when the sources are a whole plate in order, and a compaction when the wells that failed are
left out of the sources. `pool` runs many wells into one, which is a pool and a per-plate donor
pool both.
"""

from collections.abc import Iterable, Mapping, Sequence

from liulab_mbio.protocol.model import (
    FORMATS,
    Citation,
    Move,
    Plate,
    Transfer,
    Vessel,
    Well,
    row_label,
)

__all__ = [
    "FORMATS",
    "Plate",
    "Transfer",
    "Vessel",
    "Well",
    "compact",
    "interleave",
    "plate",
    "pool",
    "row_label",
    "seat",
    "wells_of",
]


def plate(
    name: str,
    wells: int,
    *,
    catalog: str = "",
    holds: str = "",
    seating: Mapping[str, str] | None = None,
    labels: Mapping[str, str] | None = None,
    note: str = "",
) -> Plate:
    """Return a plate of `wells` wells, seating `seating` and labelled `labels`.

    `seating` names what sits in a well and resolves against the protocol; `labels` describes a
    well and resolves against nothing. `seat` writes either.

    Raises
    ------
    ValueError
        If `wells` is no format a plate comes in, or a seated or labelled well is off the array.

    Examples
    --------
    >>> plate("barcodes", 384, seating={"A1": "UMI-1"}).columns
    24
    """
    return Plate(
        name,
        wells,
        catalog=catalog,
        holds=holds,
        seating=dict(seating or {}),
        labels=dict(labels or {}),
        note=note,
    )


def wells_of(one: Plate, *, start: int = 0, count: int | None = None) -> tuple[Well, ...]:
    """Return `count` of the plate's wells in reading order, from `start`."""
    names = one.well_names[start : None if count is None else start + count]
    return tuple(Well(one.name, name) for name in names)


def seat(names: Iterable[str], wells: int, *, start: int = 0) -> dict[str, str]:
    """Return `names` seated in reading order from the `start`-th well of a `wells` format.

    Raises
    ------
    ValueError
        If `wells` is no format, or the names run past the last well.

    Examples
    --------
    >>> seat(["UMI-1", "UMI-2"], 96)
    {'A1': 'UMI-1', 'A2': 'UMI-2'}
    """
    if wells not in FORMATS:
        raise ValueError(f"{wells} wells is no format; one of {', '.join(str(n) for n in FORMATS)}")
    rows, columns = FORMATS[wells]
    names = tuple(names)
    if start + len(names) > wells:
        raise ValueError(
            f"{len(names)} names from well {start + 1} run past a {wells}-well plate's last well"
        )
    places = [
        f"{row_label(row)}{column}" for row in range(rows) for column in range(1, columns + 1)
    ]
    return dict(zip(places[start:], names, strict=False))


def interleave(wells: int, into: int) -> tuple[tuple[str, ...], ...]:
    """Return a `wells` plate's well names grouped into the passes an `into` plate covers it in.

    A denser plate's wells sit at a fraction of a smaller one's spacing, so one well in that many
    lines up under a head built for the smaller format and the plate is covered in that many
    passes. Each group is in the smaller format's own reading order, so group *i*'s *j*-th name
    is where the *j*-th well of pass *i* sits.

    Raises
    ------
    ValueError
        If either count is no format a plate comes in, or the smaller format's rows and columns
        do not divide the larger one's.

    Examples
    --------
    >>> passes = interleave(384, 96)
    >>> len(passes), len(passes[0]), passes[0][:2], passes[3][:2]
    (4, 96, ('A1', 'A3'), ('B2', 'B4'))
    """
    for count in (wells, into):
        if count not in FORMATS:
            raise ValueError(
                f"{count} wells is no format; one of {', '.join(str(n) for n in FORMATS)}"
            )
    rows, columns = FORMATS[wells]
    small_rows, small_columns = FORMATS[into]
    if rows % small_rows or columns % small_columns:
        raise ValueError(
            f"a {into}-well plate's {small_rows} x {small_columns} does not divide a "
            f"{wells}-well plate's {rows} x {columns}, so no head covers one in whole passes"
        )
    down, across = rows // small_rows, columns // small_columns
    groups: list[list[str]] = [[] for _ in range(down * across)]
    for row in range(rows):
        for column in range(columns):
            groups[(row % down) * across + column % across].append(f"{row_label(row)}{column + 1}")
    return tuple(tuple(one) for one in groups)


def compact(
    sources: Sequence[Well],
    destination: Plate,
    volume_ul: float,
    *,
    title: str,
    instrument: str = "",
    start: int = 0,
    note: str = "",
    citation: Citation | None = None,
) -> Transfer:
    """Return `sources` moved into `destination`'s wells densely, in reading order from `start`.

    The caller passes the wells that still hold something, so the ones that failed are simply
    absent and the destination is dense. Hand it a whole plate's wells in order and it is a
    plain or an acoustic transfer; hand it four plates' good wells and it is a compaction.

    Raises
    ------
    ValueError
        If there are no sources, or they run past the destination's last well.

    Examples
    --------
    >>> big, small = plate("picked", 384), plate("pooled", 96)
    >>> compact(wells_of(big, count=2), small, 2.0, title="Compact").moves[1].destination
    Well(plate='pooled', well='A2')
    """
    if not sources:
        raise ValueError(f"{title}: a transfer needs at least one source well")
    places = destination.well_names[start:]
    if len(sources) > len(places):
        raise ValueError(
            f"{title}: {len(sources)} wells do not fit a {destination.wells}-well plate "
            f"from well {start + 1}"
        )
    moves = tuple(
        Move(source, Well(destination.name, place), volume_ul)
        for source, place in zip(sources, places, strict=False)
    )
    return Transfer(title, moves, instrument=instrument, note=note, citation=citation)


def pool(
    sources: Sequence[Well],
    destination: Well,
    volume_ul: float,
    *,
    title: str,
    instrument: str = "",
    note: str = "",
    citation: Citation | None = None,
) -> Transfer:
    """Return every well of `sources` run into the one `destination`.

    Raises
    ------
    ValueError
        If there are no sources.

    Examples
    --------
    >>> one = plate("lysate", 1536)
    >>> pool(wells_of(one, count=3), Well("reservoir", "1"), 3.5, title="Pool").plates
    ('lysate', 'reservoir')
    """
    if not sources:
        raise ValueError(f"{title}: a pool needs at least one source well")
    moves = tuple(Move(source, destination, volume_ul) for source in sources)
    return Transfer(title, moves, instrument=instrument, note=note, citation=citation)

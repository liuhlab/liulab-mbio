"""The `Build`: what one DMX run chooses, read from its build file and checked as it is read.

A build holds what the method leaves open — which designs are read back, from which plate of
clonal stock, on which route, and how deep into the fragment counts to go — and nothing the
method's own numbers already carry. `docs/adr/0010-method-in-code.md` draws the line, so every
field here is named and documented rather than keyed by whatever a page is missing.

A design on disk is a name and a fragment count on a two-column sheet, never a sequence, so DMX
reads back anything a lab holds.
"""

import os
from dataclasses import KW_ONLY, dataclass
from pathlib import Path

from liulab_mbio import jsonfile
from liulab_synbio.dmx.method import ROUTE_INDEX_PCR, ROUTES, Design, Route


@dataclass(frozen=True, slots=True)
class Build:
    """What one DMX run chooses, checked on construction.

    Parameters
    ----------
    name
        What the run is called.
    designs
        The two-column sheet of designs to read back; `read_designs` reads one.
    archive
        What the lab calls the plate of clonal stock each design is spotted from. DMX reads
        clonal material, so a run with no plate to spot from has nothing to read.
    route
        One key of `liulab_synbio.dmx.ROUTES`: which marking route reads these wells back.
    validate_from
        The fragment-count floor at or above which a design is read back one well at a time.
        ``0`` reads every design on the sheet.
    selection
        The drug the archive's own vector selects on. Omitted, every plate says the vector's own
        antibiotic rather than naming one this run cannot know.
    index_plate
        What the lab calls its prepared plate of barcoded primer pairs. Only the index PCR route
        takes one, so only that route may name it.

    Raises
    ------
    ValueError
        If the run or the archive is unnamed, the floor is negative, the route is neither of the
        two, or an index plate is named on the route that takes none.
    """

    name: str
    _: KW_ONLY
    designs: Path
    archive: str
    route: str
    validate_from: int
    selection: str = ""
    index_plate: str = ""

    def __post_init__(self) -> None:
        """Normalise the path, then check everything a reader of the protocol will see."""
        object.__setattr__(self, "designs", Path(self.designs))
        if not self.name.strip():
            raise ValueError("a build names the run it plans")
        if not self.archive.strip():
            raise ValueError(
                "a build names the archive plate its designs are spotted from: DMX reads clonal "
                "material, and linear vendor DNA has no colony to pick"
            )
        if self.validate_from < 0:
            raise ValueError(
                f"validate_from is {self.validate_from}, and a fragment-count floor counts "
                "fragments; set 0 to read every design on the sheet"
            )
        if self.route not in ROUTES:
            raise ValueError(
                f"route is {self.route!r}, and a build reads its wells back on one of "
                f"{', '.join(repr(one) for one in ROUTES)}"
            )
        if self.index_plate and self.route != ROUTE_INDEX_PCR.name:
            raise ValueError(
                "index_plate names a plate of barcoded primer pairs, which only the "
                f"{ROUTE_INDEX_PCR.name!r} route takes"
            )

    @property
    def marking_route(self) -> Route:
        """The route this build named, which marks every well it reads back."""
        return ROUTES[self.route]


def read_build(path: str | os.PathLike[str]) -> Build:
    """Read a build from JSON, resolving the designs sheet against the file's own directory.

    Raises
    ------
    ValueError
        On a missing or unknown key, a value of another JSON type, a sheet naming no file, or
        any check `Build` makes.

    Examples
    --------
    >>> read_build("build.json").route  # doctest: +SKIP
    'index PCR'
    """
    file = Path(path)
    data = jsonfile.read_object(path)
    jsonfile.refuse_keys(data, _REQUIRED, _OPTIONAL, "a build")
    return Build(
        jsonfile.text(data, "name", "a build's"),
        designs=jsonfile.named_file(
            file, jsonfile.text(data, "designs", "a build's"), "designs", "a build's"
        ),
        archive=jsonfile.text(data, "archive", "a build's"),
        route=jsonfile.text(data, "route", "a build's"),
        validate_from=jsonfile.whole(data, "validate_from", "a build's"),
        selection=jsonfile.text(data, "selection", "a build's") if "selection" in data else "",
        index_plate=jsonfile.text(data, "index_plate", "a build's")
        if "index_plate" in data
        else "",
    )


def read_designs(path: str | os.PathLike[str]) -> tuple[Design, ...]:
    """Read the designs to read back from a two-column sheet the user holds.

    The columns are ``name`` and ``fragments``, tab separated, one design a row. A design's
    bases are never read, so any lab's own list of what it holds is a sheet anyone can type.

    Raises
    ------
    ValueError
        If the file is not that sheet, a row names no design, or a fragment count is not a
        whole number of at least one.
    """
    rows = Path(path).read_text(encoding="utf-8-sig").splitlines()
    if not rows or rows[0].split("\t")[:2] != ["name", "fragments"]:
        raise ValueError(
            f"{os.fspath(path)} is not a designs sheet; expected tab-separated columns "
            "name, fragments"
        )
    found: list[Design] = []
    for line, row in enumerate(rows[1:], 2):
        if not row.strip():
            continue
        cells = row.split("\t")
        if not cells[0].strip():
            raise ValueError(f"line {line} of {os.fspath(path)} leaves the design unnamed")
        if len(cells) < 2:
            raise ValueError(
                f"line {line} of {os.fspath(path)} names a design and no fragment count"
            )
        count = cells[1].strip()
        if not count.isdigit() or int(count) < 1:
            raise ValueError(
                f"line {line} of {os.fspath(path)} gives {count!r} fragments, and a design is "
                "built from at least one whole piece"
            )
        found.append(Design(cells[0].strip(), int(count)))
    if not found:
        raise ValueError(f"{os.fspath(path)} names no design to read back")
    return tuple(found)


#: The keys a build is written with, and the ones it may leave out.
_REQUIRED = frozenset({"name", "designs", "archive", "route", "validate_from"})
_OPTIONAL = frozenset({"selection", "index_plate"})

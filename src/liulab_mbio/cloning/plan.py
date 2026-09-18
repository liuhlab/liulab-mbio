"""What every cloning plan is: the files it writes, its status, and the records it was made from.

A method's own plan module holds its design, its `Plan` and the file set it writes. What every
one of them repeats is here: the shape a plan and its written files take, the names of the files
any plan writes, the protocol pair, the worst-of rule a plan's status follows, where a method
puts its inserts and which way round they go, how its oligos are judged as one, and taking a
record the caller already read.

The library pipeline is not a cloning method -- `docs/adr/0004-library-rounds.md` says why --
but its plan writes the same protocol pair, so it uses this module too.
"""

import os
import typing  # Spelled out: `Protocol` here is the bench protocol imported below.
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from liulab_mbio.checks import Check, Status, worst, worst_of
from liulab_mbio.io import read_record
from liulab_mbio.primers.evaluation import PrimerReport
from liulab_mbio.primers.thresholds import reading
from liulab_mbio.protocol.model import Protocol, read_protocol, write_protocol
from liulab_mbio.protocol.render import write_html
from liulab_mbio.sequence import Feature, SequenceRecord

#: What a plan calls the product it writes: the annotated plasmid the design makes.
PRODUCT_FILE = "product.dna"

#: What a plan calls the sheet it orders its oligos from.
PRIMER_FILE = "primers.tsv"

#: What a plan calls its protocol, as the data and as the page rendered from it.
PROTOCOL_DATA_FILE = "protocol.json"
PROTOCOL_FILE = "protocol.html"

#: Where a method puts its inserts: a feature name, a span, or `None` to look for `MCS_FEATURE`.
type Site = str | tuple[int, int] | None

#: The feature a vector names its cloning site with, looked for when the caller names none.
MCS_FEATURE = "MCS"

#: Which way round an insert goes into the vector.
type Orientation = Literal["forward", "reverse"]


class Written(typing.Protocol):
    """The files one plan wrote, as anything reporting them reads that record.

    A pipeline's own `Files` names each file it writes as a field. `paths` is the order those
    files are reported in, which is the order they were written.
    """

    @property
    def paths(self) -> tuple[Path, ...]:
        """Every file the plan wrote, the first written first."""
        ...


class Planned(typing.Protocol):
    """What every pipeline's plan promises: it writes its files into one directory."""

    def write(self, directory: str | os.PathLike[str], /) -> Written:
        """Write every file into `directory`, and say where each one went."""
        ...


@dataclass(frozen=True, slots=True)
class ProtocolFiles:
    """The protocol pair a plan writes, which the plan's own file set carries.

    Parameters
    ----------
    data
        The bench protocol as JSON, which ``protocol render`` turns back into a page.
    page
        The interactive bench protocol, as one self-contained HTML page rendered from `data`.
    """

    data: Path
    page: Path


def write_protocol_files(protocol: Protocol, directory: str | os.PathLike[str]) -> ProtocolFiles:
    """Write `protocol` into `directory` as `PROTOCOL_DATA_FILE` and `PROTOCOL_FILE`.

    The directory is made when it is not there, and the page is rendered from the data as
    written, so the two cannot disagree. The same protocol writes the same bytes.
    """
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    data = write_protocol(protocol, out / PROTOCOL_DATA_FILE)
    return ProtocolFiles(data, write_html(read_protocol(data), out / PROTOCOL_FILE))


def status(checks: Iterable[Check]) -> Status:
    """Return a plan's status: the worst verdict any of its checks carries."""
    return worst_of(checks)


def as_record(value: SequenceRecord | str | os.PathLike[str]) -> SequenceRecord:
    """Read a record from a path, or take one already read."""
    return value if isinstance(value, SequenceRecord) else read_record(value)


def insertion_span(vector: SequenceRecord, site: Site) -> tuple[int, int]:
    """Return the vector bases a method's inserts replace, named as `Site` allows.

    Every method reads its insertion site the same way, so a user moving between them learns
    nothing new. The span is the model's: 0-based, half-open, and inside the vector.

    Raises
    ------
    ValueError
        If no site is named and the vector annotates no `MCS_FEATURE` feature, if it annotates
        no feature of the name given, or if the span does not lie inside the vector.
    """
    if isinstance(site, tuple):
        start, end = site
    else:
        if site is None:
            found = _named(vector, MCS_FEATURE)
            if found is None:
                raise ValueError(
                    f"name the insertion site: this vector annotates no {MCS_FEATURE!r} "
                    "feature. Pass a feature name or a (start, end) span"
                )
        else:
            found = _named(vector, site)
            if found is None:
                raise ValueError(f"this vector annotates no feature called {site!r}")
        start, end = found.segments[0].start, found.segments[-1].end
    if not 0 <= start < end <= len(vector):
        raise ValueError(
            f"the insertion site {start}-{end} does not lie inside {len(vector)} bases"
        )
    return start, end


def orientations(
    orientation: Orientation | Sequence[Orientation], count: int
) -> tuple[Orientation, ...]:
    """Spread one orientation over every insert, or take the one given for each.

    Raises
    ------
    ValueError
        If a value is neither ``forward`` nor ``reverse``, or there is not one per insert.

    Examples
    --------
    >>> orientations("reverse", 2)
    ('reverse', 'reverse')
    """
    given = (orientation,) * count if isinstance(orientation, str) else tuple(orientation)
    if len(given) != count:
        raise ValueError(f"orientation has {len(given)} values for {count} insert(s)")
    for one in given:
        if one not in ("forward", "reverse"):
            raise ValueError(f"orientation is 'forward' or 'reverse', got {one!r}")
    return given


def _named(record: SequenceRecord, name: str) -> Feature | None:
    """Return the first feature of that name, whatever its case."""
    return next(
        (feature for feature in record.features if feature.name.lower() == name.lower()), None
    )


def primer_check(reports: Sequence[PrimerReport]) -> Check:
    """Return one verdict over every oligo a plan designed, the worst of theirs.

    A count alone cannot be acted on, so the detail names each kind of check that fired with
    the rows it covers, and names apart what no sourced threshold judged.
    """
    warned = sum(1 for report in reports if report.status == "warn")
    failed = sum(1 for report in reports if report.status == "fail")
    counted = Counter(
        reading(check).label
        for report in reports
        for check in report.checks
        if check.status not in (None, "pass")
    )
    unjudged = dict.fromkeys(
        reading(check).label
        for report in reports
        for check in report.checks
        if check.status is None
    )
    said = f"{len(reports)} designed, {warned} with a warning, {failed} failing"
    if counted:
        kinds = ", ".join(
            f"{label} on {rows}"
            for label, rows in sorted(counted.items(), key=lambda one: (-one[1], one[0]))
        )
        said += f": {kinds}"
    detail = f"{said}; not judged: {', '.join(unjudged)}" if unjudged else said
    return Check("primers", worst(report.status for report in reports), len(reports), detail)

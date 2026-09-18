"""What every cloning plan is: the files it writes, its status, and the records it was made from.

A method's own plan module holds its design, its `Plan` and the file set it writes. What every
one of them repeats is here: the shape a plan and its written files take, the names of the files
any plan writes, the protocol pair, the worst-of rule a plan's status follows, and taking a
record the caller already read.

The library pipeline is not a cloning method -- `docs/adr/0004-library-rounds.md` says why --
but its plan writes the same protocol pair, so it uses this module too.
"""

import os
import typing  # Spelled out: `Protocol` here is the bench protocol imported below.
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from liulab_mbio.checks import Check, Status, worst_of
from liulab_mbio.io import read_record
from liulab_mbio.protocol.model import Protocol, read_protocol, write_protocol
from liulab_mbio.protocol.render import write_html
from liulab_mbio.sequence import SequenceRecord

#: What a plan calls the product it writes: the annotated plasmid the design makes.
PRODUCT_FILE = "product.dna"

#: What a plan calls its protocol, as the data and as the page rendered from it.
PROTOCOL_DATA_FILE = "protocol.json"
PROTOCOL_FILE = "protocol.html"


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

"""The project: everything one library build chooses, read from JSON and checked as it is read.

A project holds what the method leaves open — which positions, which proteins, which vector, how
deep to sample — and nothing the method's own molecules already carry. A second project is a
second file and no change to this package. `docs/adr/0010-method-in-code.md` draws the line.

Nothing a project states replaces a method constant. Where the two meet they compose: the
enzymes a block is kept clear of are the method's unioned with `reserved_extra`, and the barcode
length is checked against the method's cloning scar rather than against a number stated here.
"""

import json
import os
from collections.abc import Mapping, Sequence
from dataclasses import KW_ONLY, dataclass, field
from pathlib import Path
from typing import Any

from liulab_mbio.barcodes import MIN_DISTANCE, SEED
from liulab_mbio.codons import codon_tables
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_synbio.igga.method import IGGA, Scheme, refuse

#: The method's own barcode length, which a project takes unless it states another.
BARCODE_LENGTH = 11


@dataclass(frozen=True, slots=True)
class Barcode:
    """What one part's barcode holds to.

    Parameters
    ----------
    length
        How many bases name one part.
    min_distance
        How far apart two barcodes of one part list must stand, by the metric
        `liulab_mbio.barcodes` measures with.
    """

    length: int = BARCODE_LENGTH
    min_distance: int = MIN_DISTANCE


@dataclass(frozen=True, slots=True)
class Project:
    """One library build's own choices, checked on construction.

    Parameters
    ----------
    name
        What the run and each round's product is called.
    positions
        The positions, in the order the rounds fill them, so the last is the terminal one. A
        position is a name: every one of them carries the method's own stuffers.
    parts
        The FASTA of every part list, each record named for the position it fills.
    vector
        The destination vector, circular: a ``.dna``, GenBank or FASTA file.
    host
        The codon usage table the coding bases are written for.
    oligo_length
        How long one synthesised oligo of the pool may be.
    primers
        The orthogonal primer set the pool is amplified by, as a two-column sheet the user
        holds. Without one no oligo can be written, because the sites are templated on it.
    bands
        The vendor's bands for each quantity the pool is banded by, written ``low-high`` with an
        empty top for a tier with no top. They are what headroom is measured against, and they
        are reported whether or not anyone holds a price.
    batch_size
        How many parts one batch of the bench work carries.
    coverage
        Colonies over products each round is sized for.
    seed
        The seed the barcodes are drawn with.
    reserved_extra
        Further enzymes this project needs a block kept clear of, added to the method's own.
    barcode
        What one part's barcode holds to.
    scheme
        The method the project is built by. There is one, and it is `IGGA`.

    Raises
    ------
    ValueError
        If a position is repeated or missing, a number is not positive, or the barcode and the
        method's cloning scar are not whole codons together — which names ``barcode-frame``.
    KeyError
        If `reserved_extra` names an enzyme this package does not ship, or `host` no shipped
        codon usage table.
    """

    name: str
    _: KW_ONLY
    positions: tuple[str, ...]
    parts: Path
    vector: Path
    host: str
    oligo_length: int
    batch_size: int
    coverage: float
    primers: Path | None = None
    bands: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    seed: int = SEED
    reserved_extra: tuple[str, ...] = ()
    barcode: Barcode = field(default_factory=Barcode)
    scheme: Scheme = IGGA

    def __post_init__(self) -> None:
        """Normalise the lists and paths, then check everything a designer will read."""
        object.__setattr__(self, "positions", tuple(self.positions))
        object.__setattr__(self, "reserved_extra", tuple(self.reserved_extra))
        object.__setattr__(self, "parts", Path(self.parts))
        object.__setattr__(self, "vector", Path(self.vector))
        object.__setattr__(self, "bands", dict(self.bands))
        self._check_positions()
        self._check_numbers()
        self._check_barcode_frame()
        if self.host not in codon_tables():
            raise KeyError(
                f"no shipped codon usage table is called {self.host!r}: {', '.join(codon_tables())}"
            )
        for named in self.reserved_extra:
            get_enzyme(named)

    @property
    def position_count(self) -> int:
        """How many positions one product joins."""
        return len(self.positions)

    @property
    def reserved(self) -> tuple[str, ...]:
        """Every enzyme a block is kept clear of: the method's, then this project's own.

        A project adds and never replaces, so an enzyme the method reserves stays reserved
        whatever a project says.
        """
        return tuple(dict.fromkeys((*self.scheme.reserved, *self.reserved_extra)))

    @property
    def reserved_enzymes(self) -> tuple[Enzyme, ...]:
        """Those enzymes, in that order."""
        return tuple(get_enzyme(name) for name in self.reserved)

    @property
    def barcode_block_length(self) -> int:
        """How long the finished barcode block is: one barcode a position, joined by the scar."""
        return self.scheme.barcode_block_length(self.barcode.length, self.position_count)

    @property
    def retained_length(self) -> int:
        """What the finished product keeps past its last part: stuffer and barcode block."""
        return self.scheme.retained_length(self.barcode.length, self.position_count)

    def _check_positions(self) -> None:
        """Refuse a project with no position, an unnamed one, or one named twice."""
        if not self.positions:
            raise ValueError("a project needs at least one position")
        if any(not one for one in self.positions):
            raise ValueError("every position of a project is named")
        if len(set(self.positions)) != len(self.positions):
            raise ValueError(f"the positions {', '.join(self.positions)} name one of them twice")

    def _check_numbers(self) -> None:
        """Refuse a project whose lengths, counts or coverage are not positive."""
        for named, value in (
            ("oligo_length", self.oligo_length),
            ("batch_size", self.batch_size),
            ("coverage", self.coverage),
            ("barcode.length", self.barcode.length),
            ("barcode.min_distance", self.barcode.min_distance),
        ):
            if value <= 0:
                raise ValueError(f"{named} is {value}, and a project states a positive one")

    def _check_barcode_frame(self) -> None:
        """Check a barcode and the scar joining it to the last make whole codons together."""
        scar = len(self.scheme.cloning_scar)
        unit = self.barcode.length + scar
        if unit % 3:
            refuse(
                "barcode-frame",
                f"a barcode of {self.barcode.length} bases and a cloning scar of {scar} make "
                f"{unit}, which is not a whole number of codons",
            )


def read_project(path: str | os.PathLike[str]) -> Project:
    """Read a project from JSON, resolving its two file paths against the file's own directory.

    Raises
    ------
    ValueError
        On a missing or unknown key, a value of another JSON type, a path naming no file, or any
        check `Project` makes.
    KeyError
        If it names an enzyme or a codon usage table this package does not ship.

    Examples
    --------
    >>> read_project("project.json").positions  # doctest: +SKIP
    ('N', 'DBD', 'C')
    """
    file = Path(path)
    data = json.loads(file.read_text(encoding="utf-8"))
    if not isinstance(data, Mapping):
        raise ValueError(f"{os.fspath(path)} holds {type(data).__name__}, not an object")
    _keys(data, _PROJECT_KEYS, _PROJECT_OPTIONAL, "a project")
    given = dict(data)
    return Project(
        _text(given, "name", "a project"),
        positions=tuple(
            _one_text(one, f"positions[{index}]")
            for index, one in enumerate(_sequence(given, "positions", "a project"))
        ),
        parts=_file(file, _text(given, "parts", "a project"), "parts"),
        vector=_file(file, _text(given, "vector", "a project"), "vector"),
        host=_text(given, "host", "a project"),
        oligo_length=_whole(given, "oligo_length", "a project"),
        batch_size=_whole(given, "batch_size", "a project"),
        coverage=_number(given, "coverage", "a project"),
        primers=(
            _file(file, _text(given, "primers", "a project"), "primers")
            if "primers" in given
            else None
        ),
        bands=_bands(given.get("bands")),
        seed=_whole(given, "seed", "a project") if "seed" in given else SEED,
        reserved_extra=tuple(
            _one_text(one, f"reserved_extra[{index}]")
            for index, one in enumerate(
                _sequence(given, "reserved_extra", "a project") if "reserved_extra" in given else ()
            )
        ),
        barcode=_barcode(given.get("barcode")),
    )


#: The keys a project is written with, and the ones it may leave out.
_PROJECT_KEYS = frozenset(
    {
        "name",
        "positions",
        "parts",
        "vector",
        "host",
        "oligo_length",
        "batch_size",
        "coverage",
    }
)
_PROJECT_OPTIONAL = frozenset({"seed", "reserved_extra", "barcode", "primers", "bands"})
_BARCODE_OPTIONAL = frozenset({"length", "min_distance"})


def _bands(entry: Any) -> Mapping[str, tuple[str, ...]]:
    """Read the vendor's bands a project states, as a quantity naming its tiers.

    Raises
    ------
    ValueError
        If it is not an object of lists of strings.
    """
    if entry is None:
        return {}
    if not isinstance(entry, Mapping):
        raise ValueError(f"a project's bands are {type(entry).__name__}, not an object")
    return {
        quantity: tuple(
            _one_text(one, f"bands {quantity}[{index}]")
            for index, one in enumerate(_sequence(entry, quantity, "a project's"))
        )
        for quantity in entry
    }


def _barcode(entry: Any) -> Barcode:
    """Build the barcode rules from parsed JSON, the method's own where it states none.

    Raises
    ------
    ValueError
        If it is not an object of the two keys a barcode is written with.
    """
    if entry is None:
        return Barcode()
    if not isinstance(entry, Mapping):
        raise ValueError(f"a project barcode is {type(entry).__name__}, not an object")
    _keys(entry, frozenset(), _BARCODE_OPTIONAL, "a project barcode")
    return Barcode(
        _whole(entry, "length", "a project barcode") if "length" in entry else BARCODE_LENGTH,
        _whole(entry, "min_distance", "a project barcode")
        if "min_distance" in entry
        else MIN_DISTANCE,
    )


def _file(project: Path, named: str, key: str) -> Path:
    """Resolve a path a project names against the project file's own directory.

    Raises
    ------
    ValueError
        If nothing is there to read.
    """
    found = Path(named)
    resolved = found if found.is_absolute() else project.parent / found
    if not resolved.is_file():
        raise ValueError(f"a project's {key} is {named!r}, and {os.fspath(resolved)} is no file")
    return resolved


def _keys(
    data: Mapping[str, Any], required: frozenset[str], optional: frozenset[str], where: str
) -> None:
    """Refuse a mapping that is missing a key or carries one this does not read.

    Raises
    ------
    ValueError
        Naming the keys and `where` they are.
    """
    if missing := sorted(required - set(data)):
        raise ValueError(f"{where} is missing {', '.join(missing)}")
    if unknown := sorted(set(data) - required - optional):
        raise ValueError(f"{where} carries unknown key(s) {', '.join(unknown)}")


def _text(data: Mapping[str, Any], key: str, where: str) -> str:
    """Return one string value.

    Raises
    ------
    ValueError
        If the value is of another JSON type.
    """
    return _one_text(data[key], f"{where} {key}")


def _one_text(value: Any, where: str) -> str:
    """Return one string.

    Raises
    ------
    ValueError
        If the value is of another JSON type.
    """
    if not isinstance(value, str):
        raise ValueError(f"{where} is {type(value).__name__}, not a string")
    return value


def _whole(data: Mapping[str, Any], key: str, where: str) -> int:
    """Return one whole number.

    Raises
    ------
    ValueError
        If the value is of another JSON type, true and false among them.
    """
    value = data[key]
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{where} {key} is {type(value).__name__}, not a whole number")
    return value


def _number(data: Mapping[str, Any], key: str, where: str) -> float:
    """Return one number, whole or not.

    Raises
    ------
    ValueError
        If the value is of another JSON type, true and false among them.
    """
    value = data[key]
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ValueError(f"{where} {key} is {type(value).__name__}, not a number")
    return float(value)


def _sequence(data: Mapping[str, Any], key: str, where: str) -> Sequence[Any]:
    """Return one list of values.

    Raises
    ------
    ValueError
        If the value is of another JSON type, a string among them.
    """
    value = data[key]
    if isinstance(value, str) or not isinstance(value, Sequence):
        raise ValueError(f"{where} {key} is {type(value).__name__}, not a list")
    return value

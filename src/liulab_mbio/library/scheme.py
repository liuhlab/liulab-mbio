"""The scheme: the architecture one library build is given, checked as it is read.

A scheme is data the user supplies rather than data this package ships, and
``docs/adr/0005-scheme-data.md`` says why. It is checked on construction, so an invalid scheme
fails where it is read and never at the bench.

Nothing here states an overhang. A part's entry overhang is read off its own 5' external stuffer
by cutting that stuffer with the external enzyme, and the cloning-scar overhang off its 3' one, so
the overhangs a design works from are the ones the ordered DNA will leave.
"""

import json
import os
from collections.abc import Mapping, Sequence
from dataclasses import KW_ONLY, dataclass
from pathlib import Path
from typing import Any, Literal, NoReturn

from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.sites import find_sites

#: What a scheme checks on construction. A refusal leads with the name of the invariant it broke.
type Invariant = Literal[
    "internal-stuffer-prefix",
    "internal-stuffer-cuts",
    "external-stuffer-5",
    "external-stuffer-3",
    "barcode-frame",
    "terminal-frame",
    "enzyme-regions",
]


def _refuse(invariant: Invariant, detail: str) -> NoReturn:
    """Refuse a scheme, naming the invariant first so a caller reports which one broke."""
    raise ValueError(f"{invariant}: {detail}")


def _cuts(bases: str, enzyme: Enzyme) -> bool:
    """Whether `enzyme` reads a site in `bases`."""
    return bool(find_sites(SequenceRecord(bases), enzyme))


def _overhangs(bases: str, enzyme: Enzyme) -> set[str]:
    """Return every overhang `enzyme` leaves inside `bases`, as the top strand spells it.

    A cut falling off the end leaves nothing to read and is not one of these, and neither is the
    blunt end of a chopper.
    """
    return {site.overhang for site in find_sites(SequenceRecord(bases), enzyme) if site.overhang}


def _one_overhang(bases: str, enzyme: Enzyme, invariant: Invariant, where: str) -> str:
    """Return the one overhang `enzyme` leaves inside `bases`, refusing unless there is one.

    Raises
    ------
    ValueError
        Naming `invariant`, if the stuffer carries no such cut or more than one.
    """
    found = _overhangs(bases, enzyme)
    if not found:
        _refuse(
            invariant,
            f"{where} is {bases!r}, in which {enzyme.name} leaves no overhang: a stuffer carries "
            "that enzyme's site and the bases its cut leaves single-stranded",
        )
    if len(found) > 1:
        _refuse(
            invariant,
            f"{where} is {bases!r}, in which {enzyme.name} leaves {len(found)} different "
            f"overhangs, {', '.join(sorted(found))}: a part enters on one",
        )
    return found.pop()


@dataclass(frozen=True, slots=True)
class Position:
    """One position of a scheme, and the stuffers a part filling it carries.

    Parameters
    ----------
    name
        What this position is called, such as ``"N"``. A part's name says which it fills.
    internal_stuffer_prefix
        The bases a part's internal stuffer opens with, before the shared core. It ends with the
        entry overhang of the position this one admits next, which is how a part carries its
        place.
    external_stuffer_5
        The 5' flank of the synthesised block, 5' to 3' as it is ordered, up to and including the
        bases the external enzyme's cut leaves single-stranded. Those bases are this position's
        entry overhang, which the scheme reads rather than states.
    external_stuffer_3
        The 3' flank, from the bases that cut leaves single-stranded — the cloning scar — to the
        end of the block.
    """

    name: str
    _: KW_ONLY
    internal_stuffer_prefix: str
    external_stuffer_5: str
    external_stuffer_3: str

    def __post_init__(self) -> None:
        """Upper-case the three sequences."""
        for field_name in ("internal_stuffer_prefix", "external_stuffer_5", "external_stuffer_3"):
            object.__setattr__(self, field_name, getattr(self, field_name).upper())


@dataclass(frozen=True, slots=True)
class Scheme:
    """The constant architecture of one library build, checked on construction.

    Parameters
    ----------
    name
        What this scheme is called, so a report can say what a design assumed.
    positions
        The positions, in the order the rounds fill them, so the last is the terminal one. Any
        number of them: the paper's three are one instance.
    internal_enzyme
        The name of the enzyme that opens the library built so far by excising its internal
        stuffer. Its sites are the ones the product keeps, which is what lets the next round open
        it.
    external_enzyme
        The name of the enzyme that releases a part from its synthesised block. Its sites leave on
        the discarded external stuffers.
    blunt_enzymes
        The names of the blunt choppers that cut the pieces a round throws away, so that neither
        can ligate back. Any number.
    internal_stuffer_core
        The bases every internal stuffer shares, following each position's own prefix.
    cloning_scar
        The bases every part's 3' end leaves at its junction, shared by every position.
    barcode_length
        How many bases name one part.
    source
        Where these values came from, so provenance travels with the scheme. Empty where nothing
        states it.

    Raises
    ------
    ValueError
        If the scheme holds no position, its barcode is not positive, or an invariant fails —
        which names the invariant: one of `Invariant`.
    KeyError
        If a named enzyme is not one this package ships.
    """

    name: str
    _: KW_ONLY
    positions: tuple[Position, ...]
    internal_enzyme: str
    external_enzyme: str
    blunt_enzymes: tuple[str, ...]
    internal_stuffer_core: str
    cloning_scar: str
    barcode_length: int
    source: str = ""

    def __post_init__(self) -> None:
        """Normalise the sequences, then check every invariant in turn."""
        object.__setattr__(self, "positions", tuple(self.positions))
        object.__setattr__(self, "blunt_enzymes", tuple(self.blunt_enzymes))
        object.__setattr__(self, "internal_stuffer_core", self.internal_stuffer_core.upper())
        object.__setattr__(self, "cloning_scar", self.cloning_scar.upper())
        if not self.positions:
            raise ValueError("a scheme needs at least one position")
        if self.barcode_length <= 0:
            raise ValueError(f"a barcode needs a positive length, got {self.barcode_length}")
        self._check_barcode_frame()
        self._check_external_stuffers()
        self._check_internal_stuffers()
        self._check_terminal_frame()
        self._check_enzyme_regions()

    @property
    def position_count(self) -> int:
        """How many positions one product joins."""
        return len(self.positions)

    @property
    def internal(self) -> Enzyme:
        """The enzyme that opens the library built so far."""
        return get_enzyme(self.internal_enzyme)

    @property
    def external(self) -> Enzyme:
        """The enzyme that releases a part from its block."""
        return get_enzyme(self.external_enzyme)

    @property
    def blunt(self) -> tuple[Enzyme, ...]:
        """The blunt choppers, in the order this scheme names them."""
        return tuple(get_enzyme(name) for name in self.blunt_enzymes)

    def internal_stuffer(self, index: int) -> str:
        """Return the internal stuffer a part at `index` carries: its prefix and the shared core."""
        return self.positions[index].internal_stuffer_prefix + self.internal_stuffer_core

    def entry_overhang(self, index: int) -> str:
        """Return the overhang admitting a part at `index`, read off its 5' external stuffer."""
        position = self.positions[index]
        return _one_overhang(
            position.external_stuffer_5,
            self.external,
            "external-stuffer-5",
            f"the 5' external stuffer of position {position.name!r}",
        )

    @property
    def entry_overhangs(self) -> tuple[str, ...]:
        """Every position's entry overhang, in position order."""
        return tuple(self.entry_overhang(index) for index in range(self.position_count))

    @property
    def scar_overhang(self) -> str:
        """The overhang every part's 3' external stuffer yields, shared by every position."""
        position = self.positions[0]
        return _one_overhang(
            position.external_stuffer_3,
            self.external,
            "external-stuffer-3",
            f"the 3' external stuffer of position {position.name!r}",
        )

    @property
    def barcode_block_length(self) -> int:
        """How long a finished barcode block is: one barcode a position, joined by the scar."""
        barcodes = self.position_count * self.barcode_length
        return barcodes + (self.position_count - 1) * len(self.cloning_scar)

    @property
    def retained_length(self) -> int:
        """What the finished product keeps past its last part: terminal stuffer and barcode block.

        The frame of everything downstream rides on this, which is why the terminal position's
        prefix is longer than an overhang where it has to be.
        """
        return len(self.internal_stuffer(-1)) + self.barcode_block_length

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Scheme":
        """Build a scheme from parsed JSON; keys are the field names and lists become tuples.

        Raises
        ------
        ValueError
            On a missing or unknown key, a value of another JSON type, or any invariant.
        """
        _keys(data, _SCHEME_KEYS, _SCHEME_OPTIONAL, "a scheme")
        positions = _sequence(data, "positions", "a scheme")
        return cls(
            _text(data, "name", "a scheme"),
            positions=tuple(_position(entry, index) for index, entry in enumerate(positions)),
            internal_enzyme=_text(data, "internal_enzyme", "a scheme"),
            external_enzyme=_text(data, "external_enzyme", "a scheme"),
            blunt_enzymes=tuple(
                _one_text(name, f"blunt_enzymes[{index}]")
                for index, name in enumerate(_sequence(data, "blunt_enzymes", "a scheme"))
            ),
            internal_stuffer_core=_text(data, "internal_stuffer_core", "a scheme"),
            cloning_scar=_text(data, "cloning_scar", "a scheme"),
            barcode_length=_whole(data, "barcode_length", "a scheme"),
            source=_one_text(data.get("source", ""), "a scheme source"),
        )

    def _check_barcode_frame(self) -> None:
        """Check a barcode and the scar joining it to the last make whole codons together."""
        unit = self.barcode_length + len(self.cloning_scar)
        if unit % 3:
            _refuse(
                "barcode-frame",
                f"a barcode of {self.barcode_length} bases and a cloning scar of "
                f"{len(self.cloning_scar)} make {unit}, which is not a whole number of codons",
            )

    def _check_external_stuffers(self) -> None:
        """Check every position enters on one overhang, and every 3' end leaves the same scar."""
        for index in range(self.position_count):
            self.entry_overhang(index)
        scars = {position.name: self._scar(index) for index, position in enumerate(self.positions)}
        distinct = set(scars.values())
        if len(distinct) > 1:
            listed = ", ".join(f"{name} {scar}" for name, scar in scars.items())
            _refuse(
                "external-stuffer-3",
                f"the 3' external stuffers yield {len(distinct)} different overhangs — {listed} — "
                "where every part shares one cloning scar",
            )
        yielded = distinct.pop()
        if yielded != self.cloning_scar:
            _refuse(
                "external-stuffer-3",
                f"the 3' external stuffers yield {yielded!r}, which is not this scheme's cloning "
                f"scar {self.cloning_scar!r}",
            )

    def _scar(self, index: int) -> str:
        """Return the overhang the 3' external stuffer at `index` yields."""
        position = self.positions[index]
        return _one_overhang(
            position.external_stuffer_3,
            self.external,
            "external-stuffer-3",
            f"the 3' external stuffer of position {position.name!r}",
        )

    def _check_internal_stuffers(self) -> None:
        """Check each prefix names the position it admits next, and that both cuts are there."""
        for index, position in enumerate(self.positions):
            following = self.entry_overhang((index + 1) % self.position_count)
            prefix = position.internal_stuffer_prefix
            if not prefix.endswith(following):
                _refuse(
                    "internal-stuffer-prefix",
                    f"the internal stuffer prefix of position {position.name!r} is {prefix!r}, "
                    f"which does not end with {following!r}, the entry overhang of the position "
                    "it admits next",
                )
            left = _overhangs(self.internal_stuffer(index), self.internal)
            if following not in left:
                spells = ", ".join(sorted(left)) or "no overhang at all"
                _refuse(
                    "internal-stuffer-prefix",
                    f"{self.internal.name} leaves {spells} in the internal stuffer of position "
                    f"{position.name!r}, and not the {following!r} the next part enters on",
                )
            if self.scar_overhang not in left:
                spells = ", ".join(sorted(left)) or "no overhang at all"
                _refuse(
                    "internal-stuffer-cuts",
                    f"{self.internal.name} leaves {spells} in the internal stuffer of position "
                    f"{position.name!r}, and not the cloning scar {self.scar_overhang!r}: a "
                    "stuffer carries both of the cuts that open it",
                )

    def _check_terminal_frame(self) -> None:
        """Check what the product keeps past its last part is a whole number of codons."""
        retained = self.retained_length
        if retained % 3:
            terminal = self.positions[-1]
            _refuse(
                "terminal-frame",
                f"the terminal position {terminal.name!r} leaves {retained} bases in the product "
                f"— a {len(self.internal_stuffer(-1))}-base internal stuffer and a "
                f"{self.barcode_block_length}-base barcode block — which is not a whole number "
                "of codons",
            )

    def _check_enzyme_regions(self) -> None:
        """Check no enzyme reads into a region another owns, and that every chopper chops."""
        internal, external = self.internal, self.external
        for index, position in enumerate(self.positions):
            for end, bases in (
                ("5'", position.external_stuffer_5),
                ("3'", position.external_stuffer_3),
            ):
                if _cuts(bases, internal):
                    _refuse(
                        "enzyme-regions",
                        f"{internal.name} reads a site in the {end} external stuffer of position "
                        f"{position.name!r}, which {external.name} owns: opening the library "
                        "would cut the donor block as well",
                    )
            if _cuts(self.internal_stuffer(index), external):
                _refuse(
                    "enzyme-regions",
                    f"{external.name} reads a site in the internal stuffer of position "
                    f"{position.name!r}, which {internal.name} owns: releasing the part would cut "
                    "it apart",
                )
        discarded = (
            self.internal_stuffer_core,
            *(
                bases
                for position in self.positions
                for bases in (position.external_stuffer_5, position.external_stuffer_3)
            ),
        )
        for chopper in self.blunt:
            if not any(_cuts(bases, chopper) for bases in discarded):
                _refuse(
                    "enzyme-regions",
                    f"the blunt enzyme {chopper.name} reads no site in the internal stuffer core "
                    "or any external stuffer, so it chops none of the pieces a round discards",
                )


def read_scheme(path: str | os.PathLike[str]) -> Scheme:
    """Read a scheme from a JSON file the user supplies; see `Scheme.from_dict`.

    Examples
    --------
    >>> read_scheme("scheme.json").entry_overhangs  # doctest: +SKIP
    ('CTCC', 'GGAG', 'CCGA')
    """
    return Scheme.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


#: The keys a scheme and one of its positions are written with, and the ones it may leave out.
_SCHEME_KEYS = frozenset(
    {
        "name",
        "positions",
        "internal_enzyme",
        "external_enzyme",
        "blunt_enzymes",
        "internal_stuffer_core",
        "cloning_scar",
        "barcode_length",
    }
)
_SCHEME_OPTIONAL = frozenset({"source"})
_POSITION_KEYS = frozenset(
    {"name", "internal_stuffer_prefix", "external_stuffer_5", "external_stuffer_3"}
)


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


def _position(entry: Any, index: int) -> Position:
    """Build one position from parsed JSON.

    Raises
    ------
    ValueError
        If it is not a mapping of the four keys a position is written with.
    """
    where = f"position {index}"
    if not isinstance(entry, Mapping):
        raise ValueError(f"{where} is {type(entry).__name__}, not an object")
    _keys(entry, _POSITION_KEYS, frozenset(), where)
    return Position(
        _text(entry, "name", where),
        internal_stuffer_prefix=_text(entry, "internal_stuffer_prefix", where),
        external_stuffer_5=_text(entry, "external_stuffer_5", where),
        external_stuffer_3=_text(entry, "external_stuffer_3", where),
    )


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

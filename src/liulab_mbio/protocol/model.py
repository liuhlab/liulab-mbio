"""The protocol model, and loading it from JSON.

Every class refuses, with `ValueError`, a value no bench could follow: an empty title, a
non-positive volume, time, cycle count or band size, or a link that is not http(s).
"""

import json
import math
import os
from collections.abc import Callable, Mapping
from dataclasses import KW_ONLY, MISSING, dataclass, field, fields, is_dataclass
from pathlib import Path
from types import NoneType, UnionType
from typing import Any, Literal, TypeAliasType, Union, get_args, get_origin, get_type_hints

from liulab_mbio.checks import STATUSES, Status

#: The most characters `Protocol.overview` gives one card. Anything longer is a sentence, which
#: reads badly in a grid and drags the cards to different heights; it belongs in
#: `Protocol.highlights`.
OVERVIEW_CHARS = 80


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


@dataclass(frozen=True, slots=True)
class Material:
    """A reagent, kit or consumable the protocol needs. An oligo is an `Oligo`.

    Parameters
    ----------
    name
        As it is labelled on the tube or shelf.
    supplier, catalog
        Who sells it and the number to order it by. Left empty where they are not known, and
        never guessed: a catalogue number is ordered as written.
    storage
        Such as ``"-20 °C"``.
    amount
        What one run takes, such as ``"1 µL per reaction"``.
    note
        Anything else the bench needs, such as a stock concentration.
    """

    name: str
    _: KW_ONLY
    supplier: str = ""
    catalog: str = ""
    storage: str = ""
    amount: str = ""
    note: str = ""


@dataclass(frozen=True, slots=True)
class Check:
    """One pass, warn or fail verdict on the work, shown in the header as a badge.

    Parameters
    ----------
    name
        What was judged, such as ``"junctions"``.
    status
        One of `STATUSES`.
    detail
        What a reader needs besides the verdict, shown only where the verdict is not a pass.
    """

    name: str
    status: Status
    detail: str = ""

    def __post_init__(self) -> None:
        """Refuse a verdict that is not one of the three."""
        _require(
            self.status in STATUSES,
            f"check {self.name!r}: status is one of {', '.join(STATUSES)}, got {self.status!r}",
        )


@dataclass(frozen=True, slots=True)
class Oligo:
    """One synthetic DNA to order: a row of the order sheet, kept apart from the reagents.

    Parameters
    ----------
    name
        What to order it under, and what its tube is labelled.
    sequence
        5' to 3'. Shown with a copy button, and its length is counted from it.
    purpose
        What it is for, such as the title of the step that uses it.
    tm_c
        Of the part that anneals, °C, or ``None`` where none was computed.
    stock
        The working dilution, such as ``"10 µM"``.
    note
        Anything else, such as a modification or a purification.
    status
        The row's own verdict, one of `STATUSES`, or ``None`` where nothing judged it. A row
        with no verdict says so rather than reading as a pass.
    checks
        Why the verdict is not a pass: the checks that fired, each in a few words. The sheet
        shows the verdict in the row and keeps these behind a toggle, so it stays a sheet.
    """

    name: str
    sequence: str
    _: KW_ONLY
    purpose: str = ""
    tm_c: float | None = None
    stock: str = ""
    note: str = ""
    status: Status | None = None
    checks: tuple[Check, ...] = ()

    def __post_init__(self) -> None:
        """Refuse an oligo with no sequence, or a verdict better than its own checks."""
        _require(bool(self.sequence.strip()), f"oligo {self.name!r} has no sequence")
        _require(
            self.status is None or self.status in STATUSES,
            f"oligo {self.name!r}: status is one of {', '.join(STATUSES)}, got {self.status!r}",
        )
        for check in self.checks:
            _require(
                self.status is not None
                and STATUSES.index(self.status) >= STATUSES.index(check.status),
                f"oligo {self.name!r}: {check.name} is {check.status!r} and the row says "
                f"{self.status!r}",
            )


@dataclass(frozen=True, slots=True)
class Component:
    """One line of a reaction table.

    Parameters
    ----------
    name
        What to pipette.
    volume_ul
        Microlitres per reaction.
    stock, final
        Free-text concentrations, such as ``"10x"`` and ``"1x"``.
    master_mix
        ``False`` for a component added to each tube separately, such as template.
    """

    name: str
    volume_ul: float
    _: KW_ONLY
    stock: str = ""
    final: str = ""
    master_mix: bool = True

    def __post_init__(self) -> None:
        """Refuse a volume that is not positive."""
        _require(self.volume_ul > 0, f"component {self.name!r}: volume_ul must be positive")


@dataclass(frozen=True, slots=True)
class ReactionTable:
    """Per-reaction volumes, scaled to a master mix for several reactions.

    Parameters
    ----------
    components
        In pipetting order; at least one.
    title
        Such as ``"PCR master mix"``.
    reactions
        The reaction count shown first; the reader can change it.
    overage
        Extra master mix as a fraction, so ``0.1`` makes enough for 10% more reactions.
    """

    components: tuple[Component, ...]
    _: KW_ONLY
    title: str = ""
    reactions: int = 1
    overage: float = 0.1

    def __post_init__(self) -> None:
        """Refuse an empty table, fewer than one reaction, or a negative overage."""
        _require(bool(self.components), f"reaction table {self.title!r} has no component")
        _require(self.reactions >= 1, "reactions must be at least 1")
        _require(self.overage >= 0, "overage must not be negative")

    def mix_volumes(self, reactions: int) -> tuple[float | None, ...]:
        """Return each component's master-mix volume in µL, to 0.01 µL.

        ``None`` marks a component added to each tube instead.

        Examples
        --------
        >>> ReactionTable((Component("Buffer", 2.5),), overage=0.1).mix_volumes(8)
        (22.0,)
        """
        scale = reactions * (1 + self.overage)
        return tuple(
            round(c.volume_ul * scale, 2) if c.master_mix else None for c in self.components
        )


@dataclass(frozen=True, slots=True)
class Incubation:
    """One temperature held for a time.

    Parameters
    ----------
    label
        Such as ``"Annealing"``.
    temperature_c
        Degrees Celsius.
    seconds
        ``None`` holds until the reader stops it.
    """

    label: str
    temperature_c: float
    seconds: float | None

    def __post_init__(self) -> None:
        """Refuse a time that is not positive."""
        _require(
            self.seconds is None or self.seconds > 0,
            f"incubation {self.label!r}: seconds must be positive or null",
        )


@dataclass(frozen=True, slots=True)
class Stage:
    """Incubations run in order, repeated `cycles` times."""

    incubations: tuple[Incubation, ...]
    _: KW_ONLY
    cycles: int = 1

    def __post_init__(self) -> None:
        """Refuse an empty stage or fewer than one cycle."""
        _require(bool(self.incubations), "a stage needs at least one incubation")
        _require(self.cycles >= 1, "cycles must be at least 1")


@dataclass(frozen=True, slots=True)
class ThermocyclerProgram:
    """Stages run in order.

    Parameters
    ----------
    stages
        In run order; at least one.
    title
        The name to save the program under.
    lid_temperature_c
        Heated lid, or ``None`` to leave it unstated.
    """

    stages: tuple[Stage, ...]
    _: KW_ONLY
    title: str = ""
    lid_temperature_c: float | None = None

    def __post_init__(self) -> None:
        """Refuse a program with no stage."""
        _require(bool(self.stages), f"program {self.title!r} has no stage")

    @property
    def duration_seconds(self) -> float:
        """Run time at the block temperatures, leaving out ramps and indefinite holds."""
        return sum(
            stage.cycles * sum(i.seconds or 0 for i in stage.incubations) for stage in self.stages
        )


def _check_bands(bands_bp: tuple[int, ...], owner: str) -> None:
    _require(all(bp > 0 for bp in bands_bp), f"{owner}: band sizes must be positive")


@dataclass(frozen=True, slots=True)
class Ladder:
    """A named DNA size marker and its band sizes in base pairs."""

    name: str
    bands_bp: tuple[int, ...]

    def __post_init__(self) -> None:
        """Refuse a ladder with no band or a non-positive band size."""
        _require(bool(self.bands_bp), f"ladder {self.name!r} has no band")
        _check_bands(self.bands_bp, f"ladder {self.name!r}")


@dataclass(frozen=True, slots=True)
class Lane:
    """One sample lane of a simulated gel; no band sizes draws an empty lane."""

    label: str
    bands_bp: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        """Refuse a non-positive band size."""
        _check_bands(self.bands_bp, f"lane {self.label!r}")


@dataclass(frozen=True, slots=True)
class Gel:
    """A simulated agarose gel: a ladder lane followed by sample lanes."""

    ladder: Ladder
    lanes: tuple[Lane, ...]
    _: KW_ONLY
    title: str = ""

    def migration(self, bp: float) -> float:
        """Return how far a band runs: 0 for the largest band on this gel, 1 for the smallest.

        Linear in the logarithm of size, the usual approximation for agarose.
        """
        sizes = [*self.ladder.bands_bp, *(bp for lane in self.lanes for bp in lane.bands_bp)]
        top, bottom = math.log10(max(sizes)), math.log10(min(sizes))
        if top == bottom:
            return 0.5
        return (top - math.log10(bp)) / (top - bottom)


@dataclass(frozen=True, slots=True)
class Timer:
    """A countdown the reader can start from the step."""

    label: str
    seconds: float

    def __post_init__(self) -> None:
        """Refuse a time that is not positive."""
        _require(self.seconds > 0, f"timer {self.label!r}: seconds must be positive")


@dataclass(frozen=True, slots=True)
class Troubleshooting:
    """A problem the reader may see at a step, and what to do about it."""

    problem: str
    solution: str


@dataclass(frozen=True, slots=True)
class Reference:
    """A citation, with an optional http(s) link."""

    text: str
    _: KW_ONLY
    url: str = ""

    def __post_init__(self) -> None:
        """Refuse a link that is not http or https."""
        _require(
            not self.url or self.url.startswith(("http://", "https://")),
            f"reference url must be http(s), got {self.url!r}",
        )


@dataclass(frozen=True, slots=True)
class Step:
    """One numbered step of a protocol.

    Parameters
    ----------
    title
        What the step achieves, such as ``"Run the thermocycler"``.
    instructions
        Ordered actions, one sentence each.
    cautions, notes
        Shown before and after the instructions.
    tables, programs, timers
        Reaction tables, thermocycler programs and countdowns the step uses.
    gels, expected
        What a successful step looks like.
    troubleshooting
        Problems the reader may see here.
    """

    title: str
    _: KW_ONLY
    instructions: tuple[str, ...] = ()
    cautions: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    tables: tuple[ReactionTable, ...] = ()
    programs: tuple[ThermocyclerProgram, ...] = ()
    timers: tuple[Timer, ...] = ()
    gels: tuple[Gel, ...] = ()
    expected: tuple[str, ...] = ()
    troubleshooting: tuple[Troubleshooting, ...] = ()

    def __post_init__(self) -> None:
        """Refuse an empty title."""
        _require(bool(self.title.strip()), "a step needs a title")


@dataclass(frozen=True, slots=True)
class Protocol:
    """A bench protocol.

    Parameters
    ----------
    title
        The page heading.
    summary
        One paragraph: what the protocol does.
    overview
        Short facts, label to value, shown as a grid of cards. A handful of words each, at most
        `OVERVIEW_CHARS` characters, so the row scans left to right and stays one height.
    highlights
        What a fact means, a sentence each, shown as prose under the cards. A statement a reader
        has to read rather than scan goes here and not in `overview`.
    checks
        Verdicts on the work, shown as a strip of badges, so a warning is seen and not read.
    materials, oligos, equipment
        The reagents, the oligos to order, and the hardware. Three lists and not one, because an
        order sheet and a reagent list want different columns.
    steps, references
        In the order they are shown.
    """

    title: str
    _: KW_ONLY
    summary: str = ""
    overview: Mapping[str, str] = field(default_factory=dict, hash=False)
    highlights: tuple[str, ...] = ()
    checks: tuple[Check, ...] = ()
    materials: tuple[Material, ...] = ()
    oligos: tuple[Oligo, ...] = ()
    equipment: tuple[str, ...] = ()
    steps: tuple[Step, ...] = ()
    references: tuple[Reference, ...] = ()

    def __post_init__(self) -> None:
        """Refuse an empty title, or an overview value too long to be a card."""
        _require(bool(self.title.strip()), "a protocol needs a title")
        for label, value in self.overview.items():
            _require(
                len(value) <= OVERVIEW_CHARS,
                f"overview {label!r} is {len(value)} characters, over {OVERVIEW_CHARS}: a card "
                "holds a few words, so put a sentence in highlights instead",
            )

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Protocol":
        """Build a protocol from parsed JSON; keys are the field names, lists become tuples.

        Each value takes its field's JSON type: a whole number for an ``int``, any number but
        not true or false for a ``float``, and null only where the field is optional.

        Raises
        ------
        ValueError
            On an unknown or missing key or a value of another JSON type, naming where it is,
            or on a value the model refuses.
        """
        return _PROTOCOL(data, "protocol")


def read_protocol(path: str | os.PathLike[str]) -> Protocol:
    """Read a protocol from a JSON file; see `Protocol.from_dict`."""
    return Protocol.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


type _Convert = Callable[[Any, str], Any]


def _refused(where: str, expected: str, value: Any) -> ValueError:
    match value:
        case str():
            got = "a string"
        case list():
            got = "a list"
        case Mapping():
            got = "an object"
        case None | bool():
            got = json.dumps(value)
        case int() | float():
            got = repr(value)
        case _:
            got = type(value).__name__
    return ValueError(f"{where}: expected {expected}, got {got}")


def _converter(hint: Any) -> _Convert:
    """Return a converter from parsed JSON to the annotation `hint`, refusing another JSON type."""
    if isinstance(hint, TypeAliasType):
        return _converter(hint.__value__)
    origin, args = get_origin(hint), get_args(hint)
    if origin is Literal:
        return _converter(type(args[0]))
    if origin in (Union, UnionType) and len(args) == 2 and NoneType in args:
        (inner,) = (arg for arg in args if arg is not NoneType)
        return _or_null(_converter(inner))
    if origin is tuple and args[1:] == (...,):
        return _list(_converter(args[0]))
    if origin is Mapping and args[0] is str:
        return _mapping(_converter(args[1]))
    if isinstance(hint, type) and is_dataclass(hint):
        return _object(hint)
    if hint in _SCALARS:
        return _scalar(*_SCALARS[hint])
    raise TypeError(f"a protocol field has no JSON form: {hint!r}")


def _object[T](cls: type[T]) -> Callable[[Any, str], T]:
    """Return a converter from a JSON object to the dataclass `cls`, checking its keys."""
    spec = fields(cls)  # pyright: ignore[reportArgumentType]
    hints = get_type_hints(cls)
    nested = {f.name: _converter(hints[f.name]) for f in spec}
    required = {f.name for f in spec if f.default is MISSING and f.default_factory is MISSING}

    def convert(data: Any, where: str) -> T:
        if not isinstance(data, Mapping):
            raise _refused(where, "an object", data)
        if unknown := sorted(set(data) - set(nested)):
            raise ValueError(f"{where}: unknown key(s) {', '.join(unknown)}")
        if missing := sorted(required - set(data)):
            raise ValueError(f"{where}: missing key(s) {', '.join(missing)}")
        return cls(**{key: nested[key](value, f"{where}.{key}") for key, value in data.items()})

    return convert


def _list(item: _Convert) -> _Convert:
    def convert(data: Any, where: str) -> tuple[Any, ...]:
        if not isinstance(data, list):
            raise _refused(where, "a list", data)
        return tuple(item(value, f"{where}[{i}]") for i, value in enumerate(data))

    return convert


def _mapping(value: _Convert) -> _Convert:
    def convert(data: Any, where: str) -> dict[str, Any]:
        if not isinstance(data, Mapping):
            raise _refused(where, "an object", data)
        return {
            key: value(item, f"{where}[{json.dumps(key, ensure_ascii=False)}]")
            for key, item in data.items()
        }

    return convert


def _or_null(convert: _Convert) -> _Convert:
    return lambda data, where: None if data is None else convert(data, where)


def _scalar(expected: str, accepts: Callable[[Any], bool]) -> _Convert:
    def convert(data: Any, where: str) -> Any:
        if not accepts(data):
            raise _refused(where, expected, data)
        return data

    return convert


# `bool` is a subclass of `int` in Python, and true is not a number in JSON.
_SCALARS: dict[type, tuple[str, Callable[[Any], bool]]] = {
    str: ("a string", lambda value: isinstance(value, str)),
    bool: ("true or false", lambda value: isinstance(value, bool)),
    int: ("a whole number", lambda value: type(value) is int),
    float: ("a number", lambda value: type(value) in (int, float)),
}
_PROTOCOL = _object(Protocol)

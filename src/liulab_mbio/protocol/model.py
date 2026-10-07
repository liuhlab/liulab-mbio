"""The protocol model, and reading and writing it as JSON.

Every class refuses, with `ValueError`, a value no bench could follow: an empty title, a
non-positive volume, time, cycle count or band size, or a link that is not http(s).
"""

import json
import math
import os
from collections.abc import Callable, Iterable, Mapping
from dataclasses import KW_ONLY, MISSING, asdict, dataclass, field, fields, is_dataclass
from pathlib import Path
from types import MappingProxyType, NoneType, UnionType
from typing import Any, Literal, TypeAliasType, Union, get_args, get_origin, get_type_hints

from liulab_mbio.checks import STATUSES, Status

#: The most characters `Protocol.overview` gives one card. Anything longer is a sentence, which
#: reads badly in a grid and drags the cards to different heights; it belongs in
#: `Protocol.highlights`.
OVERVIEW_CHARS = 80

#: The rows and columns of each well count a plate comes in. Format is one parameter, the well
#: count, and 1536 is here because four 384-well plates compress into one.
FORMATS: Mapping[int, tuple[int, int]] = MappingProxyType(
    {12: (3, 4), 24: (4, 6), 96: (8, 12), 384: (16, 24), 1536: (32, 48)}
)


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def row_label(row: int) -> str:
    """Return the label of a 0-based plate row: ``A`` to ``Z``, then ``AA`` on.

    Examples
    --------
    >>> row_label(0), row_label(25), row_label(31)
    ('A', 'Z', 'AF')
    """
    letters = ""
    while True:
        row, rest = divmod(row, 26)
        letters = chr(ord("A") + rest) + letters
        if row == 0:
            return letters
        row -= 1


@dataclass(frozen=True, slots=True)
class Source:
    """A document a number was read from, named once and cited by key.

    Parameters
    ----------
    document
        What it is, as a reader would cite it.
    edition
        The revision, version or year it carries.
    url
        Where it was read, http(s).
    read_as
        How it was read, such as ``"plain curl"``.
    date
        When it was read.
    """

    document: str
    _: KW_ONLY
    edition: str = ""
    url: str = ""
    read_as: str = ""
    date: str = ""

    def __post_init__(self) -> None:
        """Refuse a source with no document, or a link that is not http or https."""
        _require(bool(self.document.strip()), "a source needs a document")
        _require(
            not self.url or self.url.startswith(("http://", "https://")),
            f"source url must be http(s), got {self.url!r}",
        )


@dataclass(frozen=True, slots=True)
class Citation:
    """Where in a source one row's number stands.

    Provenance is per row, not per number: a citation hangs on the component, the incubation,
    the material or the bill row that carries the number.

    Parameters
    ----------
    source
        The key of a `Protocol.sources` entry.
    locator
        Where in it, such as ``"p. 11"``, ``"Day 1.1"`` or ``"FAQ 14"``.
    """

    source: str
    locator: str = ""

    def __post_init__(self) -> None:
        """Refuse a citation naming no source."""
        _require(bool(self.source.strip()), "a citation needs a source")


#: What a rule does with what it names: keeps it out of the tube, or keeps it in.
type RuleKind = Literal["forbids", "requires"]


@dataclass(frozen=True, slots=True)
class Rule:
    """A prohibition or a requirement a material carries wherever it is used.

    A rule hangs on the material, not on a step, so it follows the material into every step that
    uses it and no edit to a step can drop it. A rule that computes a number is not one of these;
    that stays a function in `liulab_mbio.bench`.

    Parameters
    ----------
    kind
        ``"forbids"`` keeps `subject` out of a step using the material; ``"requires"`` keeps it in.
    subject
        What is forbidden or required, matched against what a step names: what it pipettes, the
        programs it runs, its title and its instructions.
    detail
        Why, in a sentence. Shown wherever the material is.
    when
        The tube contents the rule holds under, matched against every component's name and what
        each component's material brings with it. Empty holds always.
    citation
        The source that states it.
    """

    kind: RuleKind
    subject: str
    detail: str
    _: KW_ONLY
    when: str = ""
    citation: Citation | None = None

    def __post_init__(self) -> None:
        """Refuse a rule with no subject, no reason, or a kind that is neither."""
        _require(
            self.kind in ("forbids", "requires"),
            f"rule kind is forbids or requires, got {self.kind!r}",
        )
        _require(bool(self.subject.strip()), "a rule needs a subject")
        _require(bool(self.detail.strip()), f"rule {self.subject!r} needs a reason")

    def holds(self, named: Iterable[str], contents: Iterable[str]) -> bool:
        """Return whether the rule is kept by a step naming `named` with `contents` in the tube."""
        if self.when and not _names(self.when, contents):
            return True
        return _names(self.subject, named) == (self.kind == "requires")


def _names(subject: str, among: Iterable[str]) -> bool:
    """Return whether `subject` is named among `among`, whatever the case."""
    wanted = subject.casefold()
    return any(wanted in one.casefold() for one in among)


#: Why a number is missing. A ``"price"`` hole names no issue: a price nobody loaded is a missing
#: input of the user's, not a defect in what the package knows.
type HoleKind = Literal["undecided", "unpublished", "lab", "unread", "price"]


@dataclass(frozen=True, slots=True)
class Hole:
    """A number nobody sourced, standing where the number would be.

    A hole is never a value and is never judged. The field it belongs to stays empty and the
    hole stands beside it, so no loader has to read a union and no reader can mistake one for a
    figure. Filling one with a guess is the defect a hole exists to prevent.

    Parameters
    ----------
    id
        Stable between runs, such as ``"H24"``.
    missing
        What is not known, in a few words.
    kind
        Why: the method has not decided, nobody published it, it is the lab's own stock, a
        source was not read, or no price record prices it.
    where
        What the number belongs to, such as ``"ligase units per reaction"``.
    filled_by
        What would close it.
    issue
        The ticket it is routed to. A ``"price"`` hole names none.
    """

    id: str
    missing: str
    kind: HoleKind
    _: KW_ONLY
    where: str = ""
    filled_by: str = ""
    issue: str = ""

    def __post_init__(self) -> None:
        """Refuse an unnamed hole, an unknown kind, or a price hole routed to an issue."""
        _require(bool(self.id.strip()), "a hole needs an id")
        _require(bool(self.missing.strip()), f"hole {self.id!r}: say what is missing")
        _require(
            self.kind in get_args(HoleKind.__value__),
            f"hole {self.id!r}: unknown kind {self.kind!r}",
        )
        _require(
            self.kind != "price" or not self.issue,
            f"hole {self.id!r}: a missing price is a missing input, so it names no issue",
        )


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
    contains
        What it brings into the tube besides itself, such as a crowding agent in a buffer. A
        rule conditional on the tube's contents is matched against these.
    rules
        What must not, or must, happen where this material is used.
    citation
        Where its parameters were read.
    """

    name: str
    _: KW_ONLY
    supplier: str = ""
    catalog: str = ""
    storage: str = ""
    amount: str = ""
    note: str = ""
    contains: tuple[str, ...] = ()
    rules: tuple[Rule, ...] = ()
    citation: Citation | None = None


@dataclass(frozen=True, slots=True)
class Check:
    """One verdict on the work, shown in the header as a badge.

    Parameters
    ----------
    name
        What was judged, such as ``"junctions"``.
    status
        One of `STATUSES`, or ``None`` where no sourced threshold judges it. A badge with no
        verdict says so rather than reading as a pass.
    detail
        What a reader needs besides the verdict, shown wherever the verdict is not a pass.
    """

    name: str
    status: Status | None
    detail: str = ""

    def __post_init__(self) -> None:
        """Refuse a verdict that is neither one of the three nor none at all."""
        _require(
            self.status is None or self.status in STATUSES,
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
                check.status is None
                or (
                    self.status is not None
                    and STATUSES.index(self.status) >= STATUSES.index(check.status)
                ),
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
    citation
        Where the volume was read.
    """

    name: str
    volume_ul: float
    _: KW_ONLY
    stock: str = ""
    final: str = ""
    master_mix: bool = True
    citation: Citation | None = None

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
    citation
        Where the temperature and the time were read.
    """

    label: str
    temperature_c: float
    seconds: float | None
    citation: Citation | None = None

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
class Vessel:
    """Something the bench holds material in whose contents have no positions.

    A tube, a flask, a reservoir, a 25 cm plate of colonies, a flow cell. A vessel whose
    positions form an array is a `Plate`.

    Parameters
    ----------
    name
        What a step calls it.
    kind, catalog, holds, note
        What sort of vessel it is, the number to order it by, what is in it, and anything else.
    """

    name: str
    _: KW_ONLY
    kind: str = ""
    catalog: str = ""
    holds: str = ""
    note: str = ""

    def __post_init__(self) -> None:
        """Refuse a vessel with no name."""
        _require(bool(self.name.strip()), "a vessel needs a name")


@dataclass(frozen=True, slots=True)
class Plate:
    """A vessel whose positions form an array, and where each thing sits in it.

    A plate is a seating and not a second count of the reactions: it says where each one sits
    and nothing else, so it can never disagree with what the reaction says differs between them.
    Plates sit on the protocol rather than on a step, since they outlive the step that fills one.

    Parameters
    ----------
    name
        What a step calls it, and what a `Well` names.
    wells
        The format: how many wells it has, one of `FORMATS`.
    catalog, note
        The number to order it by, and anything else.
    holds
        What its wells hold where they hold nothing named, such as ``"one picked colony each"``.
    seating
        Well name to the name of what sits there. Each name resolves against the protocol's
        materials, oligos, vessels and other plates, and `Protocol.audit` reports one that does
        not; a dangling name is a defect the protocol can report.
    """

    name: str
    wells: int
    _: KW_ONLY
    catalog: str = ""
    holds: str = ""
    seating: Mapping[str, str] = field(default_factory=dict, hash=False)
    note: str = ""

    def __post_init__(self) -> None:
        """Refuse a nameless plate, a format no plate comes in, or a well off the array."""
        _require(bool(self.name.strip()), "a plate needs a name")
        _require(
            self.wells in FORMATS,
            f"plate {self.name!r}: {self.wells} wells is no format; "
            f"one of {', '.join(str(n) for n in FORMATS)}",
        )
        known = set(self.well_names)
        for well in self.seating:
            _require(well in known, f"plate {self.name!r} has no well {well!r}")

    @property
    def rows(self) -> int:
        """How many rows the format has."""
        return FORMATS[self.wells][0]

    @property
    def columns(self) -> int:
        """How many columns the format has."""
        return FORMATS[self.wells][1]

    @property
    def row_labels(self) -> tuple[str, ...]:
        """Each row's label, top to bottom."""
        return tuple(row_label(row) for row in range(self.rows))

    @property
    def well_names(self) -> tuple[str, ...]:
        """Every well, in reading order: row by row, each row left to right.

        Examples
        --------
        >>> Plate("culture", 384).well_names[:2]
        ('A1', 'A2')
        """
        return tuple(
            f"{label}{column}" for label in self.row_labels for column in range(1, self.columns + 1)
        )


@dataclass(frozen=True, slots=True)
class Well:
    """One position of a plate: the plate's name, and the well's, such as ``A1``.

    This is how a step says where a thing is instead of describing it.
    """

    plate: str
    well: str

    def __post_init__(self) -> None:
        """Refuse a well naming no plate or no position."""
        _require(bool(self.plate.strip()), "a well needs a plate")
        _require(bool(self.well.strip()), f"a well of {self.plate!r} needs a position")


@dataclass(frozen=True, slots=True)
class Move:
    """One well's worth moved into another well."""

    source: Well
    destination: Well
    volume_ul: float

    def __post_init__(self) -> None:
        """Refuse a volume that is not positive."""
        _require(self.volume_ul > 0, "a move's volume_ul must be positive")


@dataclass(frozen=True, slots=True)
class Transfer:
    """Material moved between wells, shown as a table beside a step's prose.

    One shape covers every move the bench makes: one source to one destination is a plain or an
    acoustic transfer, many sources to one destination is a pool, and a dense re-layout that
    leaves out the wells that failed is a compaction. Four names would be four things to keep
    in step.

    Parameters
    ----------
    title
        What the move achieves.
    moves
        At least one.
    instrument
        What performs it, such as an acoustic liquid handler.
    note, citation
        Anything else, and where the volume was read.
    """

    title: str
    moves: tuple[Move, ...]
    _: KW_ONLY
    instrument: str = ""
    note: str = ""
    citation: Citation | None = None

    def __post_init__(self) -> None:
        """Refuse a transfer with no title or no move."""
        _require(bool(self.title.strip()), "a transfer needs a title")
        _require(bool(self.moves), f"transfer {self.title!r} has no move")

    @property
    def plates(self) -> tuple[str, ...]:
        """Every plate the transfer touches, source before destination, each named once."""
        seen = dict.fromkeys(
            name for move in self.moves for name in (move.source.plate, move.destination.plate)
        )
        return tuple(seen)


@dataclass(frozen=True, slots=True)
class BillRow:
    """What one item of this run costs, or the hole where no record prices it.

    The quantity is computed from the design and owes nothing to a price, so it is here whether
    or not a price record was loaded. The money is a hole where nothing priced the key, and no
    figure is ever estimated.

    Parameters
    ----------
    item
        What it is.
    quantity, unit
        What this run consumes.
    key
        What a price record prices it by: a catalogue number where it has one.
    charge
        The money, written out, such as ``"12505.00"``. Empty where `hole` stands instead.
    headroom
        How far the quantity sits from the nearest band edge, in words.
    citation
        The record row that priced it.
    hole
        Stands where `charge` would, when no row prices the key.
    """

    item: str
    quantity: float
    _: KW_ONLY
    unit: str = ""
    key: str = ""
    charge: str = ""
    headroom: str = ""
    citation: Citation | None = None
    hole: Hole | None = None

    def __post_init__(self) -> None:
        """Refuse a row that is both priced and holed, and a hole of another kind."""
        _require(bool(self.item.strip()), "a bill row needs an item")
        _require(
            not (self.charge and self.hole),
            f"bill row {self.item!r} is both priced and holed",
        )
        _require(
            self.hole is None or self.hole.kind == "price",
            f"bill row {self.item!r}: an unpriced row holds a price hole",
        )


@dataclass(frozen=True, slots=True)
class Bill:
    """What this run consumes, and what it costs where a record prices it.

    Parameters
    ----------
    rows
        In the order they are shown.
    title, currency
        What to call it, and the one currency its record holds.
    total
        The sum over the rows that priced, written out. Empty where none did.
    record
        The key of the `Protocol.sources` entry naming the price record.
    """

    rows: tuple[BillRow, ...]
    _: KW_ONLY
    title: str = ""
    currency: str = ""
    total: str = ""
    record: str = ""

    def __post_init__(self) -> None:
        """Refuse an empty bill."""
        _require(bool(self.rows), f"bill {self.title!r} has no row")


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
    transfers
        What the step moves between wells. A step references a well by holding the transfer,
        never by describing where a thing is.
    gels, expected
        What a successful step looks like.
    troubleshooting
        Problems the reader may see here.
    holes
        The numbers this step would otherwise have to invent.
    """

    title: str
    _: KW_ONLY
    instructions: tuple[str, ...] = ()
    cautions: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    tables: tuple[ReactionTable, ...] = ()
    programs: tuple[ThermocyclerProgram, ...] = ()
    timers: tuple[Timer, ...] = ()
    transfers: tuple[Transfer, ...] = ()
    gels: tuple[Gel, ...] = ()
    expected: tuple[str, ...] = ()
    troubleshooting: tuple[Troubleshooting, ...] = ()
    holes: tuple[Hole, ...] = ()

    def __post_init__(self) -> None:
        """Refuse an empty title."""
        _require(bool(self.title.strip()), "a step needs a title")

    @property
    def named(self) -> tuple[str, ...]:
        """Everything the step names: its title, its actions, what it pipettes and what it runs.

        A `Rule`'s subject is matched against these.
        """
        return (
            self.title,
            *self.instructions,
            *(c.name for table in self.tables for c in table.components),
            *(p.title for p in self.programs),
            *(i.label for p in self.programs for s in p.stages for i in s.incubations),
            *(t.title for t in self.transfers),
            *(t.label for t in self.timers),
        )

    @property
    def pipetted(self) -> tuple[str, ...]:
        """What the step puts in a tube: every reaction table's component names."""
        return tuple(c.name for table in self.tables for c in table.components)


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
    vessels, plates
        What the run holds material in, and where each thing sits.
    steps, references
        In the order they are shown.
    sources
        Every document a number was read from, keyed by what a `Citation` names it.
    holes
        Numbers missing from the run as a whole. One belonging to a step sits on that step.
    bill
        What the run consumes, and what it costs where a price record prices it.
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
    vessels: tuple[Vessel, ...] = ()
    plates: tuple[Plate, ...] = ()
    steps: tuple[Step, ...] = ()
    references: tuple[Reference, ...] = ()
    sources: Mapping[str, Source] = field(default_factory=dict, hash=False)
    holes: tuple[Hole, ...] = ()
    bill: Bill | None = None

    def __post_init__(self) -> None:
        """Refuse an empty title, or an overview value too long to be a card."""
        _require(bool(self.title.strip()), "a protocol needs a title")
        for label, value in self.overview.items():
            _require(
                len(value) <= OVERVIEW_CHARS,
                f"overview {label!r} is {len(value)} characters, over {OVERVIEW_CHARS}: a card "
                "holds a few words, so put a sentence in highlights instead",
            )

    @property
    def all_holes(self) -> tuple[Hole, ...]:
        """Every hole the protocol carries, the run's own first, then each step's, then the bill's.

        Collected with their stable ids, in the order they are shown.
        """
        return (
            *self.holes,
            *(hole for step in self.steps for hole in step.holes),
            *(row.hole for row in (self.bill.rows if self.bill else ()) if row.hole),
        )

    def rules_for(self, step: Step) -> tuple[tuple[Material, Rule], ...]:
        """Return each rule that bears on `step`, with the material carrying it.

        A material bears on a step that names it, and its rules come with it: a rule cannot be
        edited out of a step because it was never written into one.
        """
        contents = self.contents_of(step)
        found: list[tuple[Material, Rule]] = []
        for material in self.materials:
            if not _names(material.name, step.named):
                continue
            found += [
                (material, rule)
                for rule in material.rules
                if not rule.when or _names(rule.when, contents)
            ]
        return tuple(found)

    def contents_of(self, step: Step) -> tuple[str, ...]:
        """Return what is in the step's tubes: what it pipettes, and what each of those brings."""
        by_name = {material.name: material for material in self.materials}
        brought = [
            one
            for pipetted in step.pipetted
            for material in (by_name.get(pipetted),)
            if material
            for one in material.contains
        ]
        return (*step.pipetted, *brought)

    def audit(self) -> tuple[Check, ...]:
        """Judge the protocol itself: its citations, its wells, its rules and its holes.

        The protocol is its own gate. A check it carries is paid by the protocol; a gate added
        to the repo's own is paid by everyone who works here afterwards.

        Examples
        --------
        >>> [check.name for check in Protocol("Demo").audit()]
        ['sources', 'wells', 'rules', 'holes']
        """
        return (self._sources(), self._wells(), self._rules(), self._holes())

    def _sources(self) -> Check:
        cited = {
            citation.source
            for citation in (
                *(m.citation for m in self.materials),
                *(r.citation for m in self.materials for r in m.rules),
                *(c.citation for s in self.steps for t in s.tables for c in t.components),
                *(
                    i.citation
                    for s in self.steps
                    for p in s.programs
                    for g in p.stages
                    for i in g.incubations
                ),
                *(t.citation for s in self.steps for t in s.transfers),
                *(row.citation for row in (self.bill.rows if self.bill else ())),
            )
            if citation
        }
        dangling = sorted(cited - set(self.sources))
        if dangling:
            return Check("sources", "fail", f"cited but not named: {', '.join(dangling)}")
        return Check("sources", "pass", f"{len(cited)} of {len(cited)} citations resolve")

    def _wells(self) -> Check:
        known = {
            *(m.name for m in self.materials),
            *(o.name for o in self.oligos),
            *(v.name for v in self.vessels),
            *(p.name for p in self.plates),
        }
        dangling = sorted(
            f"{plate.name} {well} holds {held!r}"
            for plate in self.plates
            for well, held in plate.seating.items()
            if held not in known
        )
        seated = sum(len(plate.seating) for plate in self.plates)
        if dangling:
            return Check("wells", "fail", "; ".join(dangling))
        return Check("wells", "pass", f"{seated} seated wells resolve")

    def _rules(self) -> Check:
        broken = [
            f"{step.title}: {material.name} {rule.kind} {rule.subject}"
            for step in self.steps
            for material, rule in self.rules_for(step)
            if not rule.holds(step.named, self.contents_of(step))
        ]
        kept = sum(len(self.rules_for(step)) for step in self.steps)
        if broken:
            return Check("rules", "fail", "; ".join(broken))
        return Check("rules", "pass", f"{kept} rules hold where their material is used")

    def _holes(self) -> Check:
        holes = self.all_holes
        if not holes:
            return Check("holes", "pass", "no number is missing")
        ids = ", ".join(hole.id for hole in holes)
        count = f"{len(holes)} numbers have" if len(holes) > 1 else "1 number has"
        return Check("holes", None, f"{count} no source: {ids}")

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


def write_protocol(protocol: Protocol, path: str | os.PathLike[str]) -> Path:
    """Write `protocol` to `path` as JSON that `read_protocol` reads back equal; return the path.

    Every field is written, empty and default ones too, in the order its class declares them,
    indented two spaces and as UTF-8 ending in a newline, so one protocol always writes the same
    bytes.
    """
    out = Path(path)
    text = json.dumps(asdict(protocol), ensure_ascii=False, indent=2)
    out.write_text(text + "\n", encoding="utf-8")
    return out


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

"""The protocol model, and reading and writing it as JSON.

Every class refuses, with `ValueError`, a value no bench could follow: an empty title, a
non-positive volume, time, cycle count or band size, or a link that is not http(s).
"""

import json
import math
import os
import re
from collections.abc import Callable, Iterable, Mapping
from dataclasses import KW_ONLY, dataclass, field, fields, is_dataclass, replace
from decimal import ROUND_HALF_UP, Decimal, localcontext
from pathlib import Path
from types import MappingProxyType
from typing import Any, Literal, cast, get_args

from mbio import jsonfile
from mbio.checks import STATUSES, Status

#: The most characters `Protocol.overview` gives one card. Anything longer is a sentence, which
#: reads badly in a grid and drags the cards to different heights; it belongs in
#: `Protocol.highlights`.
OVERVIEW_CHARS = 80

#: The rows and columns of each well count a plate comes in. Format is one parameter, the well
#: count, and 1536 is here because four 384-well plates compress into one.
FORMATS: Mapping[int, tuple[int, int]] = MappingProxyType(
    {12: (3, 4), 24: (4, 6), 96: (8, 12), 384: (16, 24), 1536: (32, 48)}
)

#: The longest a slug runs: a step's key, and a protocol's title in the file its page is written
#: to. Long enough to stay readable, short enough for a file name on any filesystem.
NAME_CHARS = 48


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def slug(text: str) -> str:
    """Return `text` as a handle: lowercase, one dash where it is not a letter or a digit.

    Capped at `NAME_CHARS`, so it fits a file name and reads as one word. Slugging a slug
    returns it unchanged, so a key a builder already wrote this way survives a round trip.

    Examples
    --------
    >>> slug("Set up the Golden Gate reaction")
    'set-up-the-golden-gate-reaction'
    >>> slug("Digest, then ligate")
    'digest-then-ligate'
    """
    dashed = "".join(char if char.isalnum() else "-" for char in text.casefold())
    return re.sub("-+", "-", dashed).strip("-")[:NAME_CHARS].strip("-")


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


def well_at(well: str) -> tuple[int, int] | None:
    """Return a well name's 0-based row and column, or `None` where it names neither.

    The inverse of `row_label` with the column read off the digits, so a position a plate
    could hold reads back as the place it names and anything else says it is not one.

    Examples
    --------
    >>> well_at("A1"), well_at("AF48"), well_at("reservoir")
    ((0, 0), (31, 47), None)
    """
    name = well.strip().upper()
    if not name.isascii():
        return None
    cut = len(name) - len(name.lstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))
    letters, digits = name[:cut], name[cut:]
    if not letters or not digits.isdigit():
        return None
    row = 0
    for letter in letters:
        row = row * 26 + ord(letter) - ord("A") + 1
    return row - 1, int(digits) - 1


_SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def number(value: float) -> str:
    """Return a number as a page prints it: plain text, to three significant figures.

    Fixed, with thousands separators, from 0.001 to 999,999, and never rounding away a digit
    before the point; outside that, a power of ten. Only zero prints as ``0``.

    Examples
    --------
    >>> number(1161.6), number(0.05), number(0.0009), number(2.5e-9), number(1e7), number(0)
    ('1,162', '0.05', '9 × 10⁻⁴', '2.5 × 10⁻⁹', '1 × 10⁷', '0')
    """
    if value == 0:
        return "0"
    # The float's exact value, rounded half up as `protocol.js` rounds it when the count changes.
    exact = Decimal(value)
    with localcontext(rounding=ROUND_HALF_UP):
        mantissa, power = f"{exact:.2e}".split("e")
        if 0.001 <= abs(value) < 1e6:
            fixed = f"{exact:,.{max(0, 2 - int(power))}f}"
            return fixed.rstrip("0").rstrip(".") if "." in fixed else fixed
    return f"{mantissa.rstrip('0').rstrip('.')} × 10{str(int(power)).translate(_SUPERSCRIPT)}"


@dataclass(frozen=True, slots=True)
class Source:
    """A document a claim on the page rests on, named once and cited by key.

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
    note
        The research note it was read into, as a repo-relative path. Empty where the document
        has no note of its own.
    """

    document: str
    _: KW_ONLY
    edition: str = ""
    url: str = ""
    read_as: str = ""
    date: str = ""
    note: str = ""

    def __post_init__(self) -> None:
        """Refuse a source with no document, or a link that is not http or https."""
        _require(bool(self.document.strip()), "a source needs a document")
        _require(
            not self.url or self.url.startswith(("http://", "https://")),
            f"source url must be http(s), got {self.url!r}",
        )


@dataclass(frozen=True, slots=True)
class Citation:
    """Where in a source the claim one row or sentence makes stands.

    Provenance is per row, not per number: a citation hangs on the component, the incubation,
    the cycling stage, the material or the bill row that carries the number, and on the note,
    caution or troubleshooting entry whose sentence it backs.

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
    that stays a function in `mbio.bench`.

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
        if self.when and not names(self.when, contents):
            return True
        return names(self.subject, named) == (self.kind == "requires")


def names(subject: str, among: Iterable[str]) -> bool:
    """Return whether `subject` is named among `among`, whatever the case.

    It is how a `Rule` finds its material and how a protocol split across several pages finds
    which of them a reagent belongs on.

    Examples
    --------
    >>> names("BsaI", ("Digest with BsaI-HFv2",))
    True
    """
    wanted = subject.casefold()
    return any(wanted in one.casefold() for one in among)


#: Why a number is missing. A ``"price"`` hole is a missing input of the user's, not a defect in
#: what the package knows. Only ``"unread"`` fails `Protocol.audit`: the other four name a gap no
#: source closes, which is what a finished plan keeps, while a source nobody read is work left
#: undone.
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
    """

    id: str
    missing: str
    kind: HoleKind
    _: KW_ONLY
    where: str = ""
    filled_by: str = ""

    def __post_init__(self) -> None:
        """Refuse an unnamed hole or an unknown kind."""
        _require(bool(self.id.strip()), "a hole needs an id")
        _require(bool(self.missing.strip()), f"hole {self.id!r}: say what is missing")
        _require(
            self.kind in get_args(HoleKind.__value__),
            f"hole {self.id!r}: unknown kind {self.kind!r}",
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
    cautions
        What to watch out for wherever it is used, a sentence each. It hangs here for the
        reason a rule does, and shows on every step naming the material.
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
    cautions: tuple[str, ...] = ()
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


def _sources_check(cited: frozenset[str], named: frozenset[str]) -> Check:
    """Return the `sources` verdict: every key `cited` is one of those `named`.

    A protocol is judged against the sources it names itself; a project's bill against the
    sources it and its protocols name, which are the ones the run's pages list.
    """
    dangling = sorted(cited - named)
    if dangling:
        return Check("sources", "fail", f"cited but not named: {', '.join(dangling)}")
    return Check("sources", "pass", f"{len(cited)} of {len(cited)} citations resolve")


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
class AmountToVolume:
    """A calculator on a reaction row: the volume that carries its DNA at the bench's concentration.

    The page opens it at the concentration the row's `volume_ul` assumes.

    Parameters
    ----------
    nanograms
        What the row has to carry.
    made_up_by
        The row of the same table that gives way, such as the water.
    too_dilute
        The `Troubleshooting.problem` the page shows once the DNA no longer fits, read from the
        row's own step or else another step of the protocol.
    least_ul, too_concentrated
        The least volume a method lets its row take, and the problem shown below it. Unset,
        nothing warns of a volume too small.
    """

    nanograms: float
    _: KW_ONLY
    made_up_by: str
    too_dilute: str = ""
    least_ul: float | None = None
    too_concentrated: str = ""

    def __post_init__(self) -> None:
        """Refuse an amount that is not positive, or a least volume without its problem."""
        _require(self.nanograms > 0, "a calculator's nanograms must be positive")
        _require(self.least_ul is None or self.least_ul > 0, "a calculator's least_ul is positive")
        _require(
            (self.least_ul is None) == (not self.too_concentrated),
            "a calculator's least_ul and too_concentrated come together",
        )


@dataclass(frozen=True, slots=True)
class CountToNet:
    """A calculator on an expected result: what a plated count comes to, against its floor.

    The reader counts the plate and the control beside it and says what the dilution was; the
    page takes the control off, scales back up and reads the net against `floor`. It opens at
    the floor itself, undiluted and with nothing on the control, so an untouched page states
    what the plan asked for.

    Parameters
    ----------
    floor
        The net count the design asked for. Nothing here works it out: a method computes it and
        the page only compares.
    counted
        What the plate is called, in the step's own words.
    counting
        What is counted, such as colonies, for the line that states the net.
    control
        What the control plate is called, whose count comes off. Empty where a step plates none.
    below_floor
        The `Troubleshooting.problem` the page shows once the net is short of the floor, read
        from the step's own entries or else another step's.
    """

    floor: int
    _: KW_ONLY
    counted: str
    counting: str
    control: str = ""
    below_floor: str = ""

    def __post_init__(self) -> None:
        """Refuse a floor that is not positive, or a count of nothing named."""
        _require(self.floor > 0, "a calculator's floor must be positive")
        _require(bool(self.counted.strip()), "a calculator says what was counted")
        _require(bool(self.counting.strip()), "a calculator says what it counts")


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
    calculator
        Where the bench measures this row's concentration, which the page then takes in place
        of `stock`.
    """

    name: str
    volume_ul: float
    _: KW_ONLY
    stock: str = ""
    final: str = ""
    master_mix: bool = True
    citation: Citation | None = None
    calculator: AmountToVolume | None = None

    def __post_init__(self) -> None:
        """Refuse a volume that is not positive, or a stock beside the one the bench measures."""
        _require(self.volume_ul > 0, f"component {self.name!r}: volume_ul must be positive")
        _require(
            self.calculator is None or not self.stock,
            f"component {self.name!r}: the bench measures its stock, so it states none",
        )


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
        """Refuse an empty table, no reaction, a negative overage, or a calculator's stray row."""
        _require(bool(self.components), f"reaction table {self.title!r} has no component")
        _require(self.reactions >= 1, "reactions must be at least 1")
        _require(self.overage >= 0, "overage must not be negative")
        plain = {c.name for c in self.components if c.calculator is None}
        for component in self.components:
            if component.calculator is not None:
                _require(
                    component.calculator.made_up_by in plain,
                    f"reaction table {self.title!r}: {component.name!r} is made up by "
                    f"{component.calculator.made_up_by!r}, which is no other row of it",
                )

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
    """One temperature held for a time, stepped each cycle where it is a touchdown.

    Parameters
    ----------
    label
        Such as ``"Annealing"``.
    temperature_c
        Degrees Celsius, the first cycle's where `delta_c` steps it.
    seconds
        ``None`` holds until the reader stops it.
    delta_c
        Degrees Celsius added each later cycle, negative for a touchdown; ``None`` holds the
        temperature. `last_c` reads the end off it, so a touchdown is one cycled stage.
    citation
        Where the temperature and the time were read.
    """

    label: str
    temperature_c: float
    seconds: float | None
    _: KW_ONLY
    delta_c: float | None = None
    citation: Citation | None = None

    def __post_init__(self) -> None:
        """Refuse a time that is not positive, or a step that steps nowhere."""
        _require(
            self.seconds is None or self.seconds > 0,
            f"incubation {self.label!r}: seconds must be positive or null",
        )
        _require(
            self.delta_c != 0,
            f"incubation {self.label!r}: delta_c must be non-zero or null",
        )

    def last_c(self, cycles: int) -> float:
        """Return the temperature of the `cycles`-th cycle.

        Examples
        --------
        >>> Incubation("Anneal", 68.0, 20, delta_c=-0.5).last_c(10)
        63.5
        """
        return self.temperature_c + (self.delta_c or 0) * (cycles - 1)


@dataclass(frozen=True, slots=True)
class Stage:
    """Incubations run in order, repeated `cycles` times.

    Parameters
    ----------
    incubations
        In run order; at least one.
    cycles
        How many times the stage runs. ``None`` leaves the count blank, for a stage nothing
        sources: the reader is given the rule that stops it instead of an invented number.
    citation
        Where the count was read.
    """

    incubations: tuple[Incubation, ...]
    _: KW_ONLY
    cycles: int | None = 1
    citation: Citation | None = None

    def __post_init__(self) -> None:
        """Refuse an empty stage or fewer than one cycle."""
        _require(bool(self.incubations), "a stage needs at least one incubation")
        _require(self.cycles is None or self.cycles >= 1, "cycles must be at least 1, or null")


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
    def duration_seconds(self) -> float | None:
        """Run time at the block temperatures, leaving out ramps and indefinite holds.

        ``None`` where a stage's cycle count is blank, since nothing then bounds the run.
        """
        total = 0.0
        for stage in self.stages:
            if stage.cycles is None:
                return None
            total += stage.cycles * sum(i.seconds or 0 for i in stage.incubations)
        return total

    @property
    def timer(self) -> "Timer | None":
        """The countdown the page gives this run, so no step holds a `Timer` for a program.

        It counts `duration_seconds`, and is ``None`` where that bounds nothing: a blank cycle
        count, or a program that is only an open hold.

        Examples
        --------
        >>> hold = Stage((Incubation("Hold", 4.0, None),))
        >>> ligate = Stage((Incubation("Ligate", 25.0, 600),))
        >>> ThermocyclerProgram((ligate, hold), title="Ligation").timer
        Timer(label='Ligation', seconds=600.0)
        >>> ThermocyclerProgram((hold,)).timer is None
        True
        >>> ThermocyclerProgram((Stage((ligate.incubations[0],), cycles=None),)).timer is None
        True
        """
        seconds = self.duration_seconds
        return Timer(self.title or "Thermocycler program", seconds) if seconds else None


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
class Figure:
    """A record drawn as a map beside a step, so what the step makes is shown and not described.

    A spec and never a shape, as a `Gel` and a `Plate` are: the figure says what to draw and the
    page draws it, so the figure cannot disagree with the design and the JSON stays readable.

    Parameters
    ----------
    records
        What is drawn: a path to a sequence file, relative to the directory the protocol was read
        from. Several are drawn as stacked rows, in the order given, each to the same span and
        the same switches, which is how a figure shows one molecule becoming the next.
    caption
        What the figure shows, in the words a step uses.
    span
        The stretch drawn, ``(start, end)``, 0-based and half-open, ending past the record's
        length across the origin. Null draws the whole record.
    linear
        Whether a circular record drawn whole is opened as a line.
    sequence_view
        Whether the bases are drawn under the map.
    enzymes
        The enzymes whose every cut site is drawn; null draws the shipped unique cutters.
    highlight
        What the map points at: a feature, a primer or an enzyme. Every other item dims.
    citation
        Where a figure taken from a published map was read. One computed from the design cites
        nothing.
    """

    records: tuple[str, ...]
    caption: str
    _: KW_ONLY
    span: tuple[int, int] | None = None
    linear: bool = False
    sequence_view: bool = False
    enzymes: tuple[str, ...] | None = None
    highlight: tuple[str, ...] = ()
    citation: Citation | None = None

    def __post_init__(self) -> None:
        """Refuse a figure with no caption, with no record, or with a span that is empty."""
        _require(bool(self.caption.strip()), "a figure needs a caption")
        _require(bool(self.records), f"figure {self.caption!r}: name a record to draw")
        if self.span is not None:
            start, end = self.span
            _require(
                0 <= start < end,
                f"figure {self.caption!r}: span {self.span} is 0-based and half-open, so it "
                "starts at zero or more and ends past its start",
            )


@dataclass(frozen=True, slots=True)
class Timer:
    """A countdown the reader can start from the step.

    `seconds` is the plan's time. The page lets the reader type their own, and keeps it as it
    keeps a check mark, so the protocol never holds it.
    """

    label: str
    seconds: float

    def __post_init__(self) -> None:
        """Refuse a time that is not positive."""
        _require(self.seconds > 0, f"timer {self.label!r}: seconds must be positive")


@dataclass(frozen=True, slots=True)
class Wait:
    """Time the run spends waiting on someone else, which nobody attends.

    A vendor's turnaround: an oligo pool ordered, a plate sent away to be sequenced. An
    `Incubation` is a thermocycler stage and a `Timer` is a countdown someone starts, so neither
    of them is this, and most of a real project's calendar is this.

    Parameters
    ----------
    what
        What is waited on, such as ``"the oligo pool to arrive"``.
    duration
        How long, as whoever states it does: ``"10-15 working days"``. Free text, because that
        is the honest shape, and empty is an admitted unknown rather than a zero.
    citation
        Where the turnaround was read.
    """

    what: str
    duration: str = ""
    citation: Citation | None = None

    def __post_init__(self) -> None:
        """Refuse a wait that does not say what it waits on."""
        _require(bool(self.what.strip()), "a wait needs what it waits on")


@dataclass(frozen=True, slots=True)
class Note:
    """A *why* a step shows after its instructions, and where the *why* came from.

    Parameters
    ----------
    text
        The sentence the reader sees.
    citation
        The document it was read from. Absent where the note is what this run computed or
        chose, which no document states, or what the reader checks against the page itself.
    """

    text: str
    _: KW_ONLY
    citation: Citation | None = None


@dataclass(frozen=True, slots=True)
class Caution:
    """What a step warns of before its instructions, and where the warning came from.

    Parameters
    ----------
    text
        The sentence the reader sees.
    citation
        The document it was read from. Absent where the hazard is general bench practice no
        document states, or a limit this run chose.
    """

    text: str
    _: KW_ONLY
    citation: Citation | None = None


@dataclass(frozen=True, slots=True)
class Troubleshooting:
    """A problem the reader may see at a step, and what to do about it.

    Parameters
    ----------
    problem, solution
        What the reader sees, and what to do.
    citation
        The source that gives the solution.
    """

    problem: str
    solution: str
    _: KW_ONLY
    citation: Citation | None = None


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
    labels
        Well name to text describing that well, such as ``"quarter 3"`` or its printed address.
        A label resolves against nothing and is drawn on the layout as it stands: it is what
        `holds` is for a plate, said a well at a time. A well a step has to name occupies
        `seating`, never this.
    rules
        What handling this plate forbids or requires on this page, for the same reason a
        `Material` carries its own: the rule reaches every step naming the plate, and no edit
        to a step can drop it. A plate stands on the protocol that lays it out, so a rule a
        later protocol has to keep hangs on the `Item` handed over instead, which travels.
    """

    name: str
    wells: int
    _: KW_ONLY
    catalog: str = ""
    holds: str = ""
    seating: Mapping[str, str] = field(default_factory=dict, hash=False)
    labels: Mapping[str, str] = field(default_factory=dict, hash=False)
    rules: tuple[Rule, ...] = ()
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
        for well in (*self.seating, *self.labels):
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


#: A stride in words, one entry per stride a stamp can have. The widest a format allows is 4:
#: a 1536-well plate sampled into a 96-well one takes every fourth row and column.
_STRIDES = ("", "", "every other", "every third", "every fourth")


@dataclass(frozen=True, slots=True)
class Stamp:
    """The one pattern a transfer repeats, in place of the moves that spell it out.

    Destination (row, column) comes from source (``stride`` · row + `row`, ``stride`` ·
    column + `column`), so three numbers stand for the whole move list and a page can draw
    the pattern rather than print a row each.

    Parameters
    ----------
    stride
        How far apart the source wells stand: 1 well for well, 2 every other row and column.
    row, column
        Where on the source the pattern starts, 0-based.
    """

    stride: int
    row: int
    column: int

    def __post_init__(self) -> None:
        """Refuse a stride no plate format allows, or a start off the plate."""
        _require(
            1 <= self.stride < len(_STRIDES),
            f"a stamp's stride is 1 to {len(_STRIDES) - 1}, not {self.stride}",
        )
        _require(min(self.row, self.column) >= 0, "a stamp starts on the plate")

    @property
    def start(self) -> str:
        """The source well the pattern starts at, such as ``A2``."""
        return f"{row_label(self.row)}{self.column + 1}"

    @property
    def words(self) -> str:
        """The pattern as the bench follows it: which source wells, and where it starts.

        Examples
        --------
        >>> Stamp(1, 0, 0).words, Stamp(2, 0, 1).words
        ('each well into the same well', 'every other row and column, starting A2')
        """
        if self.stride == 1:
            if not self.row and not self.column:
                return "each well into the same well"
            return f"the same layout, starting {self.start}"
        return f"{_STRIDES[self.stride]} row and column, starting {self.start}"


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

    @property
    def stamp(self) -> Stamp | None:
        """The one pattern every move repeats, or `None` where the moves are not one pattern.

        A stamp runs one volume from one plate into another and takes every destination well
        from the source at one stride. That is what the bench does by hand, so a page draws it
        once instead of printing a row per move; anything else is a list and stays one.

        Examples
        --------
        >>> picked = (Well("picked", "A1"), Well("picked", "A3"))
        >>> moves = tuple(Move(one, Well("index", f"A{n}"), 1.0) for n, one in enumerate(picked, 1))
        >>> Transfer("Sample a quarter", moves).stamp
        Stamp(stride=2, row=0, column=0)
        """
        if len(self.plates) != 2 or len({move.volume_ul for move in self.moves}) != 1:
            return None
        source, destination = self.plates
        places: list[tuple[tuple[int, int], tuple[int, int]]] = []
        for move in self.moves:
            if (move.source.plate, move.destination.plate) != (source, destination):
                return None
            at, to = well_at(move.source.well), well_at(move.destination.well)
            if at is None or to is None:
                return None
            places.append((at, to))
        for stride in range(1, len(_STRIDES)):
            offsets = {(at[0] - stride * to[0], at[1] - stride * to[1]) for at, to in places}
            if len(offsets) != 1:
                continue
            row, column = offsets.pop()
            if min(row, column) >= 0:
                return Stamp(stride, row, column)
        return None


@dataclass(frozen=True, slots=True)
class BillRow:
    """What one item of this run costs, or the hole where no record prices it.

    The quantity is computed from the design and owes nothing to a price, so it is here whether
    or not a price record was loaded. The money is a hole where nothing priced the key, and no
    figure is ever estimated. **A row carries a charge or a hole**, so an item this run names
    and nobody prices cannot reach a page as an empty cell.

    Parameters
    ----------
    item
        What it is.
    quantity, unit
        What this run consumes. A `quantity` of ``None`` is an item the design does not size,
        and `unit` then carries what the protocol states the run takes, in its own words.
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
    quantity: float | None = None
    _: KW_ONLY
    unit: str = ""
    key: str = ""
    charge: str = ""
    headroom: str = ""
    citation: Citation | None = None
    hole: Hole | None = None

    def __post_init__(self) -> None:
        """Refuse a row that is neither priced nor holed, both, or holed for another reason."""
        _require(bool(self.item.strip()), "a bill row needs an item")
        _require(
            bool(self.charge) != (self.hole is not None),
            f"bill row {self.item!r} carries a charge or a price hole, and not both",
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
        The key of the `sources` entry naming the price record, on whichever of `Protocol` and
        `Project` holds the bill.
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

    @property
    def citations(self) -> tuple[Citation, ...]:
        """Every citation this bill's rows carry, row by row."""
        return tuple(row.citation for row in self.rows if row.citation)

    @property
    def cited(self) -> frozenset[str]:
        """Every source key this bill's rows name.

        Examples
        --------
        >>> Bill((BillRow("cells", 1, charge="9.00", citation=Citation("NEB")),)).cited
        frozenset({'NEB'})
        """
        return frozenset(citation.source for citation in self.citations)


@dataclass(frozen=True, slots=True)
class Item:
    """One thing a protocol consumes or produces: the name is the contract, the rest is prose.

    A `Project` chains protocols by these names, matched as a `Plate.seating` name is. Nothing
    parses `spec`: the protocol producing an item states what it makes, the one consuming it
    restates what it needs, and the reader compares the two.

    Parameters
    ----------
    name
        What a project matches it by, so two protocols handing it over spell it alike.
    what
        What it is, for the bench.
    spec
        What it has to meet, a phrase each, such as ``"≥100 ng/µL"``.
    storage
        Where it waits until the protocol consuming it takes it, such as ``"-20 °C"``.
    rules
        What handling it forbids or requires. The rule hangs on the item, so it travels with the
        item: the protocol producing it and every protocol handed it show the same rule, and the
        page that keeps the rule is never the only page that states it.
    """

    name: str
    what: str
    _: KW_ONLY
    spec: tuple[str, ...] = ()
    storage: str = ""
    rules: tuple[Rule, ...] = ()

    def __post_init__(self) -> None:
        """Refuse an item with no name, or one nothing is said about."""
        _require(bool(self.name.strip()), "an item needs a name")
        _require(bool(self.what.strip()), f"item {self.name!r}: say what it is")


@dataclass(frozen=True, slots=True)
class Topic:
    """One thing a run's reader is told before any of its protocols: a heading and paragraphs.

    It explains what no one protocol can — why the run is shaped as it is — so a step never has
    to stop and explain a decision.

    Parameters
    ----------
    title
        The heading.
    body
        A paragraph each.
    """

    title: str
    body: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Refuse a topic with no title."""
        _require(bool(self.title.strip()), "a topic needs a title")


@dataclass(frozen=True, slots=True)
class Step:
    """One numbered step of a protocol.

    Parameters
    ----------
    title
        What the step achieves, such as ``"Run the thermocycler"``.
    key
        What this step is, whatever its title is reworded to: the handle its anchor and the
        bench's check mark are kept under. The builder assigns it, for the step's job and not
        its wording, and it is slugged as a title is. Empty, it falls back to the title's slug.
    section
        What stage of the protocol the step belongs to, such as ``"Day 1"``. A label and not a
        container: the steps stay one list and the numbering runs through it.
    instructions
        Ordered actions, one sentence each.
    cautions
        Shown before the instructions. A bare string is a caution citing nothing.
    notes
        Shown after them. A bare string is a note citing nothing.
    tables, programs, timers
        Reaction tables, thermocycler programs and countdowns the step uses.
    waits
        Time the step spends waiting on someone else, which nobody attends.
    hands_on_seconds
        How much of the step someone stands over, which `waits`, `timers` and `programs` say
        nothing about. ``None`` is an admitted unknown and never a zero, and the unattended
        share is the time held less this wherever both are known. Nature Protocols publishes
        both numbers, so a protocol is expected to know them:
        ``docs/research/gantt-conventions.md`` section 5.
    transfers
        What the step moves between wells. A step references a well by holding the transfer,
        never by describing where a thing is.
    figures
        The records the step draws, shown under its instructions.
    gels, expected
        What a successful step looks like.
    calculator
        Where the reader has a count of their own to work out and read against what `expected`
        states, shown under it.
    troubleshooting
        Problems the reader may see here.
    holes
        The numbers this step would otherwise have to invent.
    """

    title: str
    _: KW_ONLY
    key: str = ""
    section: str = ""
    instructions: tuple[str, ...] = ()
    cautions: tuple[Caution | str, ...] = ()
    notes: tuple[Note | str, ...] = ()
    tables: tuple[ReactionTable, ...] = ()
    programs: tuple[ThermocyclerProgram, ...] = ()
    timers: tuple[Timer, ...] = ()
    waits: tuple[Wait, ...] = ()
    hands_on_seconds: int | None = None
    transfers: tuple[Transfer, ...] = ()
    figures: tuple[Figure, ...] = ()
    gels: tuple[Gel, ...] = ()
    expected: tuple[str, ...] = ()
    calculator: CountToNet | None = None
    troubleshooting: tuple[Troubleshooting, ...] = ()
    holes: tuple[Hole, ...] = ()

    def __post_init__(self) -> None:
        """Refuse an empty title or a negative hands-on time; slug the key; make each text one."""
        _require(bool(self.title.strip()), "a step needs a title")
        _require(
            self.hands_on_seconds is None or self.hands_on_seconds >= 0,
            f"step {self.title!r}: hands_on_seconds is zero or more, or null where nobody "
            "stated it",
        )
        object.__setattr__(self, "key", slug(self.key or self.title))
        object.__setattr__(
            self, "notes", tuple(Note(n) if isinstance(n, str) else n for n in self.notes)
        )
        object.__setattr__(
            self,
            "cautions",
            tuple(Caution(c) if isinstance(c, str) else c for c in self.cautions),
        )

    @property
    def cautioned(self) -> tuple[Caution, ...]:
        """Every caution, each a `Caution`: `__post_init__` makes one of any bare string.

        Examples
        --------
        >>> Step("Thaw the mix", cautions=("Keep it on ice.",)).cautioned
        (Caution(text='Keep it on ice.', citation=None),)
        """
        return cast("tuple[Caution, ...]", self.cautions)

    @property
    def noted(self) -> tuple[Note, ...]:
        """Every note, each a `Note`: `__post_init__` makes one of any bare string.

        Examples
        --------
        >>> Step("Rest the tube", notes=("It settles.",)).noted
        (Note(text='It settles.', citation=None),)
        """
        return cast("tuple[Note, ...]", self.notes)

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
    def held_seconds(self) -> float | None:
        """How long this step holds the bench, where anything in it bounds the time.

        Summed over its timers and its thermocycler programs. ``None`` where nothing bounds
        one — a step with neither, or one whose only program runs to an open end — because a
        step holding nothing is a question and never a zero.

        Examples
        --------
        >>> Step("Rest the tube", key="rest", timers=(Timer("rest", 600),)).held_seconds
        600.0
        >>> Step("Mix the reaction", key="mix").held_seconds is None
        True
        """
        bounded = [float(timer.seconds) for timer in self.timers]
        bounded += [p.timer.seconds for p in self.programs if p.timer is not None]
        return sum(bounded) if bounded else None

    @property
    def pipetted(self) -> tuple[str, ...]:
        """What the step puts in a tube: every reaction table's component names."""
        return tuple(c.name for table in self.tables for c in table.components)


def sectioned(section: str, *steps: Step) -> tuple[Step, ...]:
    """Return `steps` labelled as one stretch of a protocol.

    A section labels a run of consecutive steps, so it is named once where the list is
    assembled. Sentence case, in the voice a step title has, naming what the bench achieves
    over that stretch.

    Examples
    --------
    >>> made = sectioned("Day 1", Step("Thaw the cells"), Step("Plate them"))
    >>> [step.section for step in made]
    ['Day 1', 'Day 1']
    """
    return tuple(replace(step, section=section) for step in steps)


def figured(step: Step, figure: Figure | None) -> Step:
    """Return `step` showing `figure`, or unchanged where there is no record to draw.

    What a step draws is the plan's to choose, so a builder that knows the chemistry hands the
    step over without one.

    Examples
    --------
    >>> drawn = figured(Step("Ligate"), Figure(("product.dna",), "the join"))
    >>> drawn.figures[0].caption
    'the join'
    >>> figured(Step("Ligate"), None).figures
    ()
    """
    return step if figure is None else replace(step, figures=(figure,))


@dataclass(frozen=True, slots=True)
class Protocol:
    """A bench protocol.

    Parameters
    ----------
    title
        The page heading.
    key
        What the rendered page remembers the bench's check marks under, minted from the content
        where a pipeline writes the JSON and left alone after, so an agent's edit keeps the ticks
        already made. Empty where nobody has minted one.
    summary
        One paragraph: what the protocol does.
    overview
        Short facts, label to value, shown as a grid of cards. A handful of words each, at most
        `OVERVIEW_CHARS` characters, so the row scans left to right and stays one height.
    highlights
        What a fact means, a sentence each, shown as prose under the cards. A statement a reader
        has to read rather than scan goes here and not in `overview`.
    background
        Why the work is shaped as it is, a topic at a time, read before the first step. A run of
        several protocols says it once, on `Project.background`; a protocol written as a page of
        its own says it here, so the rationale reaches its reader either way.
    checks
        Verdicts on the work, shown as a strip of badges, so a warning is seen and not read.
    choice
        The job this protocol is one way of doing, where several ways are offered. Protocols
        naming the same job are the ways: they stand together in `Project.protocols`, take one
        place in the run, and the bench does one of them. It is read at the bench, so it is a
        lowercase job that sits inside "one of two ways to read every well back". Empty where
        the protocol is simply a step of the chain.
    consumes, produces
        What the bench is handed before this protocol, and what it leaves for the next one. A
        `Project` chains protocols by these names; a protocol rendered alone states them for its
        reader. Everything else a protocol needs it declares for itself, in `materials`,
        `equipment` and `oligos`.
    materials, oligos, equipment
        The reagents, the oligos to order, and the hardware. Three lists and not one, because an
        order sheet and a reagent list want different columns.
    order_sheet
        The file holding `oligos` as a sheet to send a supplier, as a path from the page. The
        page links it rather than being the thing that is ordered from.
    files
        The files the run writes, each as a path from the page. Rendered text naming one by its
        file name links it, so a bare filename is never one the reader cannot open. A file the
        page never names costs nothing, so a run declares what it writes once rather than each
        page declaring what it mentions.
    vessels, plates
        What the run holds material in, and where each thing sits.
    steps
        In the order they are shown.
    sources
        Every document a claim on the page rests on, keyed by what a `Citation` names it. A
        page lists only the ones its citations name.
    holes
        Numbers missing from the run as a whole. One belonging to a step sits on that step.
    bill
        What the run consumes, and what it costs where a price record prices it.
    """

    title: str
    _: KW_ONLY
    key: str = ""
    summary: str = ""
    overview: Mapping[str, str] = field(default_factory=dict, hash=False)
    highlights: tuple[str, ...] = ()
    background: tuple[Topic, ...] = ()
    checks: tuple[Check, ...] = ()
    choice: str = ""
    consumes: tuple[Item, ...] = ()
    produces: tuple[Item, ...] = ()
    materials: tuple[Material, ...] = ()
    oligos: tuple[Oligo, ...] = ()
    order_sheet: str = ""
    files: tuple[str, ...] = ()
    equipment: tuple[str, ...] = ()
    vessels: tuple[Vessel, ...] = ()
    plates: tuple[Plate, ...] = ()
    steps: tuple[Step, ...] = ()
    sources: Mapping[str, Source] = field(default_factory=dict, hash=False)
    holes: tuple[Hole, ...] = ()
    bill: Bill | None = None

    def __post_init__(self) -> None:
        """Refuse an empty title, or an overview value too long to be a card; slug the key."""
        _require(bool(self.title.strip()), "a protocol needs a title")
        object.__setattr__(self, "key", slug(self.key))
        for label, value in self.overview.items():
            _require(
                len(value) <= OVERVIEW_CHARS,
                f"overview {label!r} is {len(value)} characters, over {OVERVIEW_CHARS}: a card "
                "holds a few words, so put a sentence in highlights instead",
            )

    @property
    def all_holes(self) -> tuple[Hole, ...]:
        """Every hole the protocol carries, the run's own first, then each step's, then the bill's.

        One id is one hole: a number the run and a step both miss is the same number, so it is
        collected once, where it is first shown.
        """
        found: dict[str, Hole] = {}
        for hole in (
            *self.holes,
            *(hole for step in self.steps for hole in step.holes),
            *(row.hole for row in (self.bill.rows if self.bill else ()) if row.hole),
        ):
            found.setdefault(hole.id, hole)
        return tuple(found.values())

    @property
    def held_seconds(self) -> tuple[float, int]:
        """How long this protocol holds the bench, and how many of its steps hold nothing.

        One place computes it, so a protocol's own page and the schedule a project draws can
        never disagree. A step nothing times is counted, never summed as a zero.

        Examples
        --------
        >>> steps = (Step("Mix", key="mix"), Step("Rest", key="rest", timers=(Timer("rest", 60),)))
        >>> Protocol("Demo", steps=steps).held_seconds
        (60.0, 1)
        """
        held = [step.held_seconds for step in self.steps]
        return sum((one for one in held if one is not None), 0.0), held.count(None)

    def materials_for(self, step: Step) -> tuple[Material, ...]:
        """Return the materials `step` names, which are the ones whose facts reach it.

        What a material carries reaches a step through this and nowhere else.
        """
        return tuple(one for one in self.materials if names(one.name, step.named))

    def rules_for(self, step: Step) -> tuple[tuple[str, Rule], ...]:
        """Return each rule that bears on `step`, with the name of whatever carries it.

        A material, a plate, or anything the protocol is handed or leaves bears on a step that
        names it, and its rules come with it: a rule cannot be edited out of a step because it
        was never written into one. Only the carrier's name is handed back, which is all a
        reader and the audit need.
        """
        contents = self.contents_of(step)
        handled: tuple[Material | Plate | Item, ...] = (
            *self.plates,
            *self.consumes,
            *self.produces,
        )
        carriers: tuple[Material | Plate | Item, ...] = (
            *self.materials_for(step),
            *(one for one in handled if names(one.name, step.named)),
        )
        return tuple(
            (carrier.name, rule)
            for carrier in carriers
            for rule in carrier.rules
            if not rule.when or names(rule.when, contents)
        )

    def cautions_for(self, step: Step) -> tuple[Caution, ...]:
        """Return every caution `step` shows: its materials' first, then its own, each once.

        A step's own stands where no material carries one. Where one sentence is both, the
        one naming a document wins, so deduplication never costs the page a citation.
        """
        carried = (
            Caution(one) for material in self.materials_for(step) for one in material.cautions
        )
        each: dict[str, Caution] = {}
        for one in (*carried, *step.cautioned):
            if one.text not in each or (one.citation and not each[one.text].citation):
                each[one.text] = one
        return tuple(each.values())

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

    @property
    def citations(self) -> tuple[Citation, ...]:
        """Every citation this protocol's own rows and its bill carry, kind by kind.

        Examples
        --------
        >>> Protocol("Demo").citations
        ()
        """
        return tuple(
            citation
            for citation in (
                *(m.citation for m in self.materials),
                *(r.citation for m in self.materials for r in m.rules),
                *(c.citation for s in self.steps for t in s.tables for c in t.components),
                *(g.citation for s in self.steps for p in s.programs for g in p.stages),
                *(
                    i.citation
                    for s in self.steps
                    for p in s.programs
                    for g in p.stages
                    for i in g.incubations
                ),
                *(t.citation for s in self.steps for t in s.transfers),
                *(f.citation for s in self.steps for f in s.figures),
                *(t.citation for s in self.steps for t in s.troubleshooting),
                *(c.citation for s in self.steps for c in s.cautioned),
                *(n.citation for s in self.steps for n in s.noted),
                *(w.citation for s in self.steps for w in s.waits),
            )
            if citation
        ) + (self.bill.citations if self.bill else ())

    @property
    def cited(self) -> frozenset[str]:
        """Every source key this protocol's own citations name.

        Examples
        --------
        >>> Protocol("Demo").cited
        frozenset()
        """
        return frozenset(citation.source for citation in self.citations)

    def _sources(self) -> Check:
        return _sources_check(self.cited, frozenset(self.sources))

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
            f"{step.title}: {carrier} {rule.kind} {rule.subject}"
            for step in self.steps
            for carrier, rule in self.rules_for(step)
            if not rule.holds(step.named, self.contents_of(step))
        ]
        kept = sum(len(self.rules_for(step)) for step in self.steps)
        if broken:
            return Check("rules", "fail", "; ".join(broken))
        return Check("rules", "pass", f"{kept} rules hold where their carrier is used")

    def _holes(self) -> Check:
        holes = self.all_holes
        if not holes:
            return Check("holes", "pass", "no number is missing")
        unread = [hole.id for hole in holes if hole.kind == "unread"]
        if unread:
            return Check("holes", "fail", f"a source nobody read would fill: {', '.join(unread)}")
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


def by_place[T](items: Iterable[T], choice: Callable[[T], str]) -> tuple[tuple[T, ...], ...]:
    """Return `items` one place of the run at a time, `choice` naming the job each is a way of.

    The ways of one job take a single place, so they are numbered alike and listed as one entry;
    everything else stands alone. Nothing stores an ordinal.

    Examples
    --------
    >>> by_place(("", "read every well back", "read every well back", ""), lambda one: one)
    (('',), ('read every well back', 'read every well back'), ('',))
    """
    groups: list[list[T]] = []
    last = ""
    for item in items:
        here = choice(item)
        if not here or here != last:
            groups.append([])
        groups[-1].append(item)
        last = here
    return tuple(tuple(group) for group in groups)


def left_by_every(protocols: Iterable[Protocol]) -> set[str]:
    """Return the names every one of `protocols` produces, which is what one place hands on.

    Where they are the ways of one job the bench did only one of them, so a name one way makes
    and another does not is not certainly there.

    Examples
    --------
    >>> left_by_every((Protocol("A", produces=(Item("calls", "per well"),)),))
    {'calls'}
    """
    return set.intersection(*({item.name for item in one.produces} for one in protocols))


def _leaving(every: set[str]) -> str:
    """Return what every way of one job leaves behind, as the `choices` check words it."""
    if not every:
        return "nothing"
    return "the same thing" if len(every) == 1 else f"the same {len(every)} things"


@dataclass(frozen=True, slots=True)
class Project:
    """Protocols run in order, each handed what the ones before it produced.

    A protocol is a document someone follows in one sitting; a project is the run they are part
    of. The chain is by name alone: a protocol consuming ``"entry clone"`` is handed whatever
    the project was given or an earlier protocol produced under that name, and `audit` reports a
    name nothing hands over as a badge rather than refusing to build the project. A chain with a
    dangling input is still a document someone can read.

    A protocol may name the job it is one way of doing. The ways of one job stand together and
    take one place in the run, the bench does one of them, and what passes out of that place is
    what every way leaves behind.

    Parameters
    ----------
    title
        What the run is called.
    key
        What the run's own pages remember under, as `Protocol.key` is. Each protocol carries its
        own and no two may share one, since a key names one page's store.
    summary
        One paragraph: what the run achieves.
    background
        What the reader is told before the first protocol: why the run is shaped as it is, a
        topic at a time. Explanation a step would otherwise carry belongs here.
    files
        The files the run writes, as `Protocol.files` holds them, for the run's own pages.
    inputs
        What the bench already holds before the first protocol.
    protocols
        In the order they are run.
    checks
        Verdicts on the design, which one protocol of the run cannot judge alone. A verdict on
        one protocol's own work stays on that protocol.
    sources
        Every document the run's own pages cite, as `Protocol.sources` holds a page's. The
        record pricing `bill` is one.
    bill
        What the run consumes, and what it costs where a price record prices it. It is the
        run's, not each protocol's, because two protocols buying the same cells would otherwise
        be counted twice.
    """

    title: str
    _: KW_ONLY
    key: str = ""
    summary: str = ""
    background: tuple[Topic, ...] = ()
    files: tuple[str, ...] = ()
    inputs: tuple[Item, ...] = ()
    protocols: tuple[Protocol, ...] = ()
    checks: tuple[Check, ...] = ()
    sources: Mapping[str, Source] = field(default_factory=dict, hash=False)
    bill: Bill | None = None

    def __post_init__(self) -> None:
        """Refuse a project with no title, two protocols keyed alike, or ways standing apart.

        Slug the key. The ways of one job take one place in the run, and a place is contiguous,
        so ways with another protocol between them are a builder defect rather than a badge.
        """
        _require(bool(self.title.strip()), "a project needs a title")
        object.__setattr__(self, "key", slug(self.key))
        keys = [one.key for one in self.protocols if one.key]
        shared = next((key for key in keys if keys.count(key) > 1), "")
        _require(
            not shared,
            f"two protocols are keyed {shared!r}: a key names one page's store, so a protocol "
            "copied from another needs its own key or none",
        )
        places = by_place(self.protocols, lambda one: one.choice)
        jobs = [group[0].choice for group in places if group[0].choice]
        apart = next((job for job in jobs if jobs.count(job) > 1), "")
        _require(
            not apart,
            f"ways to {apart!r} stand apart in protocols: the ways of one job take one place in "
            "the run, and a place is contiguous, so they go next to each other",
        )

    def audit(self) -> tuple[Check, ...]:
        """Judge the chain, the choices it offers, and the sources its own bill cites.

        A consumed name resolves to an input or an earlier protocol's output, and a bill row's
        citation to a source the run or one of its protocols names. Each protocol judges its own
        citations. A run offering no choice is judged by two checks, as it always was.

        Examples
        --------
        >>> [check.name for check in Project("Demo").audit()]
        ['handoffs', 'sources']
        """
        offered = (self._choices(),) if any(one.choice for one in self.protocols) else ()
        return (self._handoffs(), *offered, self._sources())

    def _handoffs(self) -> Check:
        handed = {item.name for item in self.inputs}
        dangling: list[str] = []
        consumed = 0
        for group in by_place(self.protocols, lambda one: one.choice):
            for protocol in group:
                consumed += len(protocol.consumes)
                dangling += [
                    f"{protocol.title} consumes {item.name!r}, which nothing hands it"
                    for item in protocol.consumes
                    if item.name not in handed
                ]
            handed |= left_by_every(group)
        if dangling:
            return Check("handoffs", "fail", "; ".join(dangling))
        counted = (
            f"{consumed} consumed items resolve" if consumed != 1 else "1 consumed item resolves"
        )
        return Check("handoffs", "pass", counted)

    def _choices(self) -> Check:
        wrong: list[str] = []
        unhelped: list[str] = []
        said: list[str] = []
        # A topic is linked by the slug of its title, so that is what names it here: a
        # title opening in capitals still carries the guidance a lowercase job asks for.
        titled = {slug(topic.title) for topic in self.background}
        for group in by_place(self.protocols, lambda one: one.choice):
            job = group[0].choice
            if not job:
                continue
            every = left_by_every(group)
            odd = sorted({item.name for one in group for item in one.produces} - every)
            if len(group) < 2:
                wrong.append(f"Only one way to {job} is written, so there is nothing to choose.")
            elif odd:
                listed = ", ".join(f'"{name}"' for name in odd)
                comes = "comes" if len(odd) == 1 else "come"
                wrong.append(
                    f"The ways to {job} do not leave the same things: "
                    f"{listed} {comes} from only one of them."
                )
            else:
                said.append(f"{len(group)} ways to {job}, each leaving {_leaving(every)}.")
            if slug(job) not in titled:
                unhelped.append(f"Nothing on the overview says how to pick a way to {job}.")
        if wrong:
            return Check("choices", "fail", " ".join(wrong))
        if unhelped:
            return Check("choices", "warn", " ".join(unhelped))
        return Check("choices", "pass", " ".join(said))

    def _sources(self) -> Check:
        named = frozenset(self.sources).union(
            key for protocol in self.protocols for key in protocol.sources
        )
        return _sources_check(self.bill.cited if self.bill else frozenset(), named)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Project":
        """Build a project from parsed JSON, by the rules `Protocol.from_dict` follows."""
        return _PROJECT(data, "project")


def citing(protocol: Protocol) -> Protocol:
    """Return `protocol` with the sources its own citations name, and no others.

    A builder hands in every document the method might read from; which of them this run cited
    depends on the steps it built. Dropping the rest keeps the data to documents the reader can
    follow back to a sentence on the page.

    Examples
    --------
    >>> one = citing(Protocol("Demo", sources={"M0491": Source("NEB M0491")}))
    >>> dict(one.sources)
    {}
    """
    kept = {key: source for key, source in protocol.sources.items() if key in protocol.cited}
    return replace(protocol, sources=kept)


def read_protocol(path: str | os.PathLike[str]) -> Protocol:
    """Read a protocol from a JSON file; see `Protocol.from_dict`."""
    return Protocol.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


def write_protocol(protocol: Protocol, path: str | os.PathLike[str]) -> Path:
    """Write `protocol` to `path` as JSON that `read_protocol` reads back equal; return the path.

    Every field is written, empty and default ones too, in the order its class declares them,
    indented two spaces and as UTF-8 ending in a newline, so one protocol always writes the same
    bytes.
    """
    return _write(protocol, path)


def read_project(path: str | os.PathLike[str]) -> Project:
    """Read a project from a JSON file; see `Project.from_dict`."""
    return Project.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


def write_project(project: Project, path: str | os.PathLike[str]) -> Path:
    """Write `project` to `path` as JSON that `read_project` reads back equal; return the path.

    Its protocols are written where they stand, by the rule `write_protocol` follows, so the
    chain is one file an agent edits and renders again.
    """
    return _write(project, path)


def _write(what: Protocol | Project, path: str | os.PathLike[str]) -> Path:
    out = Path(path)
    out.write_text(json.dumps(_plain(what), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return out


def _plain(value: Any) -> Any:
    """Return `value` as JSON data, writing a note or caution citing nothing as its text alone.

    So the common sentence stays the bare string a hand-edited file writes, and only one
    carrying a citation grows an object.
    """
    match value:
        case Note(citation=None) | Caution(citation=None):
            return value.text
        case _ if is_dataclass(value) and not isinstance(value, type):
            return {f.name: _plain(getattr(value, f.name)) for f in fields(value)}
        case Mapping():
            return {key: _plain(item) for key, item in value.items()}
        case tuple() | list():
            return [_plain(item) for item in value]
        case _:
            return value


_PROTOCOL = jsonfile.reader(Protocol)
_PROJECT = jsonfile.reader(Project)

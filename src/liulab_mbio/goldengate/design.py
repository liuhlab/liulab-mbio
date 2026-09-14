"""Choosing the Type IIS enzyme an assembly uses, and the overhangs its junctions are cut to.

Three decisions, in the order a design makes them:

- **Which enzyme.** Rank the candidates by how many sites they read in the parts that end up in
  the product. An enzyme with no site needs nothing done to the parts, so domestication is
  proposed only when no candidate is free, and then as a report of what it would change.
- **Which overhangs.** Each junction takes an overhang of the length the enzyme leaves. A
  palindrome would let a fragment ligate to itself, a repeat would let two junctions swap, and
  a near-duplicate is what mis-ligates; each refusal is reported with the rule that made it.
- **How well the set should ligate.** Pryor et al. 2020 measured every overhang pair for five
  enzymes, and that data ships here: fidelity is the product over the junctions of correct
  ligations over all ligations. An enzyme they did not measure is scored against a ligase
  profile the caller holds (`liulab_mbio.goldengate.ligase`) where there is one and by the rules
  where there is not, and the report says which of the three scored it.

The fidelity data is `src/liulab_mbio/data/ligation_fidelity.json`, from the supplementary
tables of Pryor, J.M., Potapov, V., Kucera, R.B., Bilotti, K., Cantor, E.J. and Lohman, G.J.S.
(2020) *PLoS One* 15(9): e0238592, used under CC BY 4.0.
`docs/research/ligation-fidelity.md` records the licence, the axis convention and the rules.
"""

import json
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import KW_ONLY, dataclass, field
from functools import cache
from importlib.resources import files
from itertools import product
from typing import Literal

from liulab_mbio.codons import CodonUsage
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.goldengate.ligase import LigaseProfile
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand, reverse_complement
from liulab_mbio.sites import (
    CutSite,
    Domestication,
    EnzymeLike,
    domesticate,
    primer_tail,
    site_counts,
)

#: The Type IIS enzymes a design ranks by default: the ones NEB's Ligase Master Mix table
#: covers, plus BtgZI, which that table does not.
GOLDEN_GATE_ENZYMES: tuple[str, ...] = (
    "BbsI",
    "BsaI",
    "BsmBI",
    "BspQI",
    "Esp3I",
    "PaqCI",
    "SapI",
    "BtgZI",
)

#: Ranked behind every other enzyme whatever its site count: NEB publishes no Golden Gate
#: protocol for BtgZI, and its cut ends re-ligate poorly.
LAST_RESORT = frozenset({"BtgZI"})

#: How many bases two overhangs in one set must differ by. The two-mismatch rule is the one
#: modular cloning standards state; Potapov 2018 measured it as stricter than it needs to be,
#: which is why it is an argument and not a constant.
MIN_DISTANCE = 2

#: Ligations per 100,000 events at or above which NEB's Ligase Fidelity Viewer calls a
#: Watson-Crick pair strong. Potapov 2018 normalises to that many events and the Viewer's
#: thresholds are stated in those units, not in one matrix's own totals.
STRONG_LIGATION = 100

#: Ligations per 100,000 events at or above which it calls a mismatch modest, one to avoid.
MODEST_MISMATCH = 10

#: What the rule-based fallback takes off a junction for a near-duplicate partner and for an
#: overhang of one base kind. Ranking choices, not measurements: a set scored this way is
#: comparable with another scored this way and with nothing else.
_NEAR_PENALTY = 0.05
_UNIFORM_PENALTY = 0.02

_RULE_SOURCE = "rule-based estimate: no published ligation data covers this enzyme"

#: Why a candidate overhang was refused. `_refuse` returns every one but ``"stop"``, which belongs
#: to a caller reading a candidate in frame, as `liulab_mbio.library.standard` does.
type RejectionRule = Literal[
    "length", "palindrome", "uniform", "repeat", "near-duplicate", "site", "stop"
]

#: What a set of overhangs is scored against: the enzyme's own matrix, or a ligase's profile.
type Scoring = LigationMatrix | LigaseProfile


@dataclass(frozen=True, slots=True)
class LigationMatrix:
    """How often each overhang pair was seen ligating, measured with one enzyme.

    Parameters
    ----------
    enzyme, product
        The enzyme, and the supplier's product the measurement used.
    table
        Which supplementary table of the paper this is.
    overhang_length
        How many bases the overhangs on both axes have.
    cycling_celsius
        The two temperatures the measured reaction was cycled between.
    observations
        Every ligation event counted.
    counts
        Top-strand overhang to the bottom-strand overhangs it was seen ligating to. Both are
        written 5' to 3', so a row pairs with the column spelling its reverse complement.
    citation
        The paper the data comes from.
    """

    enzyme: str
    product: str
    table: str
    _: KW_ONLY
    overhang_length: int
    cycling_celsius: tuple[int, int]
    observations: int
    counts: Mapping[str, Mapping[str, int]] = field(hash=False)
    citation: str = ""

    @property
    def overhangs(self) -> tuple[str, ...]:
        """Every overhang the measurement covers."""
        return tuple(self.counts)

    @property
    def source(self) -> str:
        """Where a report should say this number came from."""
        return f"{self.citation} {self.table}, measured with {self.product}"

    def count(self, top: str, bottom: str) -> int:
        """How often a top-strand overhang was seen ligating to a bottom-strand one."""
        return self.counts.get(top.upper(), {}).get(bottom.upper(), 0)

    def normalised(self, top: str, bottom: str) -> float:
        """Return the count per 100,000 ligation events, the scale NEB's thresholds use."""
        return 100_000 * self.count(top, bottom) / self.observations


@dataclass(frozen=True, slots=True)
class Ligation:
    """What the data says about one overhang of a set.

    Parameters
    ----------
    overhang
        The overhang, written on the top strand.
    correct
        Ligations of its two ends to each other, which is what the junction is for.
    total
        Ligations of its two ends to any end the reaction holds, correct and mismatched.
    """

    overhang: str
    correct: int
    total: int

    @property
    def value(self) -> float:
        """The share of this junction's ligations that were correct."""
        return self.correct / self.total if self.total else 0.0


@dataclass(frozen=True, slots=True)
class FidelityReport:
    """How well a set of overhangs should ligate, and where the number came from.

    Parameters
    ----------
    enzyme
        The enzyme whose data scored the set.
    source
        The measurement, or a statement that the rules scored it instead.
    measured
        ``False`` when no published data covers this enzyme and the rules stood in.
    value
        The probability that every junction ligates to its own partner.
    ligations
        One entry per overhang, empty when the rules scored the set.
    weak
        Overhangs whose Watson-Crick pair was seen fewer than `STRONG_LIGATION` times per
        100,000 ligation events.
    mismatches
        Top overhang, bottom overhang and normalised count, for every cross pair seen at least
        `MODEST_MISMATCH` times per 100,000 ligation events, the worst first.
    enzyme_specific
        ``False`` when a ligase profile scored the set: a measurement of the ligase and the
        conditions, and not of this enzyme.
    """

    enzyme: str
    source: str
    _: KW_ONLY
    measured: bool
    value: float
    ligations: tuple[Ligation, ...] = ()
    weak: tuple[str, ...] = ()
    mismatches: tuple[tuple[str, str, float], ...] = ()
    enzyme_specific: bool = True

    @property
    def label(self) -> str:
        """What kind of number this is, for a report printing it beside the value."""
        if not self.measured:
            return "rule-based estimate"
        if not self.enzyme_specific:
            return f"measured ligase profile, not specific to {self.enzyme}"
        return "measured"


@dataclass(frozen=True, slots=True)
class EnzymeChoice:
    """One enzyme a design could use, and what using it would cost.

    Parameters
    ----------
    enzyme
        The enzyme.
    sites
        How many sites it reads across the parts.
    measured
        Whether shipped ligation data can score its overhangs.
    changes
        The synonymous codon changes that would take its sites away. Empty while another
        candidate is free, because domestication is proposed only when none is.
    outside_cds
        Sites lying in no coding sequence, which a human has to decide about.
    unchanged
        Sites in a coding sequence that no synonymous change could take away.
    """

    enzyme: Enzyme
    sites: int
    _: KW_ONLY
    measured: bool = False
    changes: tuple[Domestication, ...] = ()
    outside_cds: tuple[CutSite, ...] = ()
    unchanged: tuple[CutSite, ...] = ()

    @property
    def free(self) -> bool:
        """Whether the parts hold no site of this enzyme already."""
        return self.sites == 0

    @property
    def clean(self) -> bool:
        """Whether the parts can be made free of it without anyone deciding anything."""
        return not (self.outside_cds or self.unchanged)


@dataclass(frozen=True, slots=True)
class Junction:
    """Where two parts meet in the product, and how free its overhang is.

    Parameters
    ----------
    name
        What a report calls this junction.
    record
        The part the junction is read from, for a scarless or in-frame junction.
    position
        Where the two parts separate: the 0-based top-strand index the downstream part
        begins at.
    scarless
        Take the overhang from `record`, so the product gains no bases at this junction.
    in_frame
        A scarless junction inside a coding sequence that has to read through it: the junction
        sits on a codon boundary and moves only by whole codons.
    window
        How many bases the junction may move by to get past a rule. A scarless junction that
        moves does not change what the product spells, only where the cut falls.
    overhang
        Fixed by the caller. Still held to every rule, and the refusal says which one.

    Raises
    ------
    ValueError
        If a fixed overhang is asked to be scarless too, or a scarless junction names no
        record, or the window is negative.
    """

    name: str
    _: KW_ONLY
    record: SequenceRecord | None = None
    position: int = 0
    scarless: bool = False
    in_frame: bool = False
    window: int = 0
    overhang: str | None = None

    def __post_init__(self) -> None:
        """Refuse a junction whose fields ask for two different things."""
        if self.overhang is not None and (self.scarless or self.in_frame):
            raise ValueError(
                f"junction {self.name!r} both fixes an overhang and reads one from a record"
            )
        if (self.scarless or self.in_frame) and self.record is None:
            raise ValueError(f"junction {self.name!r} needs a record to read its overhang from")
        if self.window < 0:
            raise ValueError(f"junction {self.name!r} has a negative window")

    @property
    def fixed(self) -> bool:
        """Whether the caller already chose this junction's overhang."""
        return self.overhang is not None


@dataclass(frozen=True, slots=True)
class Rejection:
    """One candidate overhang a rule refused.

    Parameters
    ----------
    overhang
        The candidate.
    rule
        Which rule refused it.
    detail
        Why that rule refused this candidate.
    """

    overhang: str
    rule: RejectionRule
    detail: str


@dataclass(frozen=True, slots=True)
class Choice:
    """The overhang one junction took.

    Parameters
    ----------
    junction
        The junction this is for.
    overhang
        The overhang it takes, written on the top strand.
    offset
        How far the junction moved from the position it was asked for.
    rejected
        Every candidate refused before this one, in the order they were tried.
    """

    junction: Junction
    overhang: str
    offset: int = 0
    rejected: tuple[Rejection, ...] = ()


@dataclass(frozen=True, slots=True)
class OverhangSet:
    """The overhangs a design chose, one per junction.

    Parameters
    ----------
    enzyme
        The enzyme that will cut them.
    choices
        One per junction, in the order the junctions were given.
    fidelity
        How well the set should ligate.
    """

    enzyme: Enzyme
    choices: tuple[Choice, ...]
    fidelity: FidelityReport

    @property
    def overhangs(self) -> tuple[str, ...]:
        """The chosen overhangs, in the order the junctions were given."""
        return tuple(choice.overhang for choice in self.choices)


@cache
def _shipped() -> Mapping[str, LigationMatrix]:
    """Read the shipped ligation data, keyed by the enzyme it was measured with."""
    text = (files("liulab_mbio") / "data" / "ligation_fidelity.json").read_text(encoding="utf-8")
    document = json.loads(text)
    citation = document["source"]["citation"]
    return {
        one["enzyme"]: LigationMatrix(
            one["enzyme"],
            one["product"],
            one["table"],
            overhang_length=one["overhang_length"],
            cycling_celsius=tuple(one["cycling_celsius"]),
            observations=one["observations"],
            counts=one["counts"],
            citation=citation,
        )
        for one in document["matrices"]
    }


def ligation_matrix(enzyme: EnzymeLike) -> LigationMatrix | None:
    """Return the shipped ligation data measured with this enzyme, or ``None`` where there is none.

    Matched by name. The measurement is of one supplier's product cycled at one pair of
    temperatures, so an isoschizomer sold by somebody else does not inherit it.

    Examples
    --------
    >>> ligation_matrix("BsaI").count("TTTT", "AAAA")
    635
    >>> ligation_matrix("PaqCI") is None
    True
    """
    return _shipped().get(_one(enzyme).name)


def choose_enzyme(
    parts: Iterable[SequenceRecord],
    *,
    enzymes: Iterable[EnzymeLike] = GOLDEN_GATE_ENZYMES,
    usage: CodonUsage | None = None,
    avoid: Iterable[EnzymeLike] = (),
) -> tuple[EnzymeChoice, ...]:
    """Rank Type IIS enzymes for assembling these parts, best first.

    The parts are the pieces that end up in the product, so a site anywhere in them is a site
    the enzyme would cut during assembly. An enzyme with no site is preferred; only when none
    is free is domestication proposed, and then every candidate carries what it would change.

    A tie between free enzymes goes to the one whose overhangs shipped data can score, then to
    the one leaving the longer overhang. `LAST_RESORT` enzymes rank behind everything.

    Parameters
    ----------
    parts
        The records going into the assembly.
    enzymes
        The candidates. Every one must be Type IIS.
    usage
        The host's codon usage, for a proposed domestication.
    avoid
        Further enzymes whose sites a domestication may not create.

    Raises
    ------
    ValueError
        If a candidate cuts inside its own recognition site, leaving no overhang to design.
    """
    records = tuple(parts)
    candidates = _resolve(enzymes)
    for one in candidates:
        _type_iis(one)
    counts = site_counts(records, candidates)
    anyone_free = any(counts[one.name] == 0 for one in candidates)
    choices: list[EnzymeChoice] = []
    for one in candidates:
        changes, outside, unchanged = _proposal(records, one, usage, avoid, anyone_free)
        choices.append(
            EnzymeChoice(
                one,
                counts[one.name],
                measured=ligation_matrix(one) is not None,
                changes=changes,
                outside_cds=outside,
                unchanged=unchanged,
            )
        )
    return tuple(sorted(choices, key=_rank))


def _proposal(
    records: Sequence[SequenceRecord],
    enzyme: Enzyme,
    usage: CodonUsage | None,
    avoid: Iterable[EnzymeLike],
    anyone_free: bool,
) -> tuple[tuple[Domestication, ...], tuple[CutSite, ...], tuple[CutSite, ...]]:
    """Run domestication over the parts without keeping the edits, and gather what it reports."""
    if anyone_free:
        return (), (), ()
    changes: list[Domestication] = []
    outside: list[CutSite] = []
    unchanged: list[CutSite] = []
    for record in records:
        _, report = domesticate(record, enzyme, usage=usage, avoid=avoid)
        changes.extend(report.changes)
        outside.extend(report.outside_cds)
        unchanged.extend(report.unchanged)
    return tuple(changes), tuple(outside), tuple(unchanged)


def _rank(choice: EnzymeChoice) -> tuple[bool, int, bool, bool, int, str]:
    """Order the candidates: usable protocol, fewest sites, no decision to make, then data."""
    return (
        choice.enzyme.name in LAST_RESORT,
        choice.sites,
        not choice.clean,
        not choice.measured,
        -choice.enzyme.overhang_length,
        choice.enzyme.name,
    )


def design_overhangs(
    junctions: Iterable[Junction],
    enzyme: EnzymeLike,
    *,
    avoid: Iterable[EnzymeLike] = (),
    min_distance: int = MIN_DISTANCE,
    allow_uniform: bool = False,
    profile: LigaseProfile | None = None,
    prefer_profile: bool = False,
) -> OverhangSet:
    """Choose one overhang per junction, and say what was refused on the way.

    Junctions whose overhang is already fixed are settled first, then the ones read from a
    record, then the free ones, so the constrained junctions are never blocked by a free one.
    The result lists them in the order they were given.

    Parameters
    ----------
    junctions
        The junctions to design for.
    enzyme
        The enzyme that will cut them, which sets the overhang length.
    avoid
        Enzymes besides this one whose sites a primer tail carrying the overhang must not spell.
    min_distance
        How many bases two overhangs in the set must differ by.
    allow_uniform
        Accept an overhang of one base kind. Refused by default: an all-GC junction truncates.
    profile
        A ligase's own matrix, read by `liulab_mbio.goldengate.ligase.read_profile`. It ranks
        the free candidates and scores the set where no shipped matrix covers the enzyme.
    prefer_profile
        Use `profile` even where a shipped matrix covers the enzyme.

    Raises
    ------
    ValueError
        If a junction has no candidate left, naming the rule that refused the last one, or if
        a profile covers overhangs of another length than the enzyme leaves.

    Examples
    --------
    >>> record = SequenceRecord("AAAACCTGAGGGTTTT")
    >>> junction = Junction("insert", record=record, position=4, scarless=True)
    >>> design_overhangs([junction], "BsaI").overhangs
    ('CCTG',)
    """
    one = _type_iis(_one(enzyme))
    others = _resolve(avoid)
    table, _ = _scoring(one, profile, prefer_profile)
    wanted = tuple(junctions)
    order = sorted(range(len(wanted)), key=lambda index: _freedom(wanted[index]))
    chosen: dict[int, Choice] = {}
    taken: list[str] = []
    for index in order:
        junction = wanted[index]
        rejected: list[Rejection] = []
        for offset, candidate in _candidates(junction, one, table):
            refusal = _refuse(candidate, one, taken, others, min_distance, allow_uniform)
            if refusal is None:
                chosen[index] = Choice(junction, candidate, offset, tuple(rejected))
                taken.append(candidate)
                break
            rejected.append(refusal)
        else:
            raise _stuck(junction, rejected)
    choices = tuple(chosen[index] for index in range(len(wanted)))
    scored = fidelity(
        [choice.overhang for choice in choices],
        one,
        profile=profile,
        prefer_profile=prefer_profile,
    )
    return OverhangSet(one, choices, scored)


def _freedom(junction: Junction) -> int:
    """How much choice a junction leaves, least first, which is the order to settle them in."""
    if junction.fixed:
        return 0
    if junction.in_frame:
        return 1
    if junction.scarless:
        return 2
    return 3


def _stuck(junction: Junction, rejected: Sequence[Rejection]) -> ValueError:
    """Return the error for a junction every candidate was refused for."""
    if not rejected:
        return ValueError(f"junction {junction.name!r} has no candidate overhang at all")
    last = rejected[-1]
    rules = ", ".join(sorted({one.rule for one in rejected}))
    return ValueError(
        f"junction {junction.name!r} has no overhang left: {rules} refused every candidate, "
        f"the last being {last.overhang!r} because {last.detail}"
    )


def _candidates(
    junction: Junction, enzyme: Enzyme, table: Scoring | None
) -> Iterator[tuple[int, str]]:
    """Every overhang this junction could take, best first, with how far it moved to get it."""
    length = enzyme.overhang_length
    if junction.overhang is not None:
        yield 0, junction.overhang.upper()
        return
    if junction.scarless or junction.in_frame:
        # The constructor refuses a scarless junction naming no record.
        record = junction.record
        assert record is not None
        step = 3 if junction.in_frame else 1
        if junction.in_frame:
            _codon_boundary(record, junction.position)
        for offset in _offsets(junction.window, step):
            bases = _read(record, junction.position + offset, length)
            if bases is not None:
                yield offset, bases
        return
    free = ("".join(bases) for bases in product("ACGT", repeat=length))
    if table is None:
        yield from ((0, candidate) for candidate in free)
        return
    ranked = sorted(free, key=lambda one: (-table.count(one, reverse_complement(one)), one))
    yield from ((0, candidate) for candidate in ranked)


def _offsets(window: int, step: int) -> Iterator[int]:
    """Where the junction may sit, nearest to the asked position first."""
    yield 0
    for distance in range(step, window + 1, step):
        yield -distance
        yield distance


def _read(record: SequenceRecord, position: int, length: int) -> str | None:
    """Return the bases a junction here would leave, or ``None`` when they run off the record."""
    if position < 0:
        return None
    try:
        return record.extract(Segment(position, position + length))
    except ValueError:
        return None


def _codon_boundary(record: SequenceRecord, position: int) -> None:
    """Refuse a junction that would fall inside a codon of the coding sequence it lies in.

    Raises
    ------
    ValueError
        If no coding sequence covers the position, or the junction is not on a codon boundary.
    """
    for feature in record.features:
        if feature.type != "CDS":
            continue
        offset = _frame_offset(feature, position)
        if offset is None:
            continue
        if offset % 3:
            raise ValueError(
                f"position {position} is {offset % 3} base(s) past a codon boundary of "
                f"{feature.name!r}, so a junction there would not read in frame"
            )
        return
    raise ValueError(f"no coding sequence covers position {position}, so nothing holds a frame")


def _frame_offset(feature: Feature, position: int) -> int | None:
    """How far into a coding sequence a top-strand position is, or ``None`` when outside it."""
    segments = feature.segments
    if feature.strand == Strand.REVERSE:
        segments = tuple(reversed(segments))
    seen = 0
    for segment in segments:
        if segment.start <= position <= segment.end:
            if feature.strand == Strand.REVERSE:
                return seen + segment.end - position
            return seen + position - segment.start
        seen += segment.end - segment.start
    return None


def _refuse(
    candidate: str,
    enzyme: Enzyme,
    taken: Sequence[str],
    avoid: Sequence[Enzyme],
    min_distance: int,
    allow_uniform: bool,
) -> Rejection | None:
    """Why this candidate will not do, or ``None`` when it will."""
    length = enzyme.overhang_length
    if len(candidate) != length or set(candidate) - set("ACGT"):
        return Rejection(
            candidate,
            "length",
            f"{enzyme.name} leaves a {length}-base overhang of A, C, G and T",
        )
    if candidate == reverse_complement(candidate):
        return Rejection(
            candidate,
            "palindrome",
            "it reads the same on both strands, so a fragment carrying it ligates to itself",
        )
    if not allow_uniform and (set(candidate) <= set("GC") or set(candidate) <= set("AT")):
        kind = "G and C" if set(candidate) <= set("GC") else "A and T"
        return Rejection(
            candidate,
            "uniform",
            f"it is {kind} throughout, and such a junction is prone to truncation",
        )
    for other in taken:
        if candidate in (other, reverse_complement(other)):
            return Rejection(
                candidate,
                "repeat",
                f"{other} is already a junction, and two junctions sharing an overhang swap",
            )
    for other in taken:
        for partner in (other, reverse_complement(other)):
            distance = _distance(candidate, partner)
            if distance < min_distance:
                return Rejection(
                    candidate,
                    "near-duplicate",
                    f"it differs from {partner} in {distance} base(s), fewer than {min_distance}",
                )
    try:
        primer_tail(enzyme, candidate, avoid=avoid)
    except ValueError as error:
        return Rejection(candidate, "site", f"no primer tail carries it: {error}")
    return None


def _distance(one: str, other: str) -> int:
    """How many positions two overhangs of one length differ in."""
    return sum(a != b for a, b in zip(one, other, strict=True))


def refusal(
    candidate: str,
    enzyme: EnzymeLike,
    *,
    taken: Iterable[str] = (),
    avoid: Iterable[EnzymeLike] = (),
    min_distance: int = MIN_DISTANCE,
    allow_uniform: bool = False,
) -> Rejection | None:
    """Why `candidate` will not join `taken`, or ``None`` when it will.

    The rules `design_overhangs` chooses by, for a caller searching for a set of its own and
    needing to weigh one candidate against a partial set rather than design a whole one.

    Examples
    --------
    >>> refusal("AGGT", "BsaI") is None
    True
    >>> refusal("AGGT", "BsaI", taken=["AGGT"]).rule
    'repeat'
    """
    return _refuse(
        candidate.upper(),
        _type_iis(_one(enzyme)),
        tuple(one.upper() for one in taken),
        _resolve(avoid),
        min_distance,
        allow_uniform,
    )


def fidelity(
    overhangs: Iterable[str],
    enzyme: EnzymeLike,
    *,
    profile: LigaseProfile | None = None,
    prefer_profile: bool = False,
) -> FidelityReport:
    """Score how well a set of overhangs should ligate to their own partners and nothing else.

    With shipped data, this is Pryor 2020's definition: the product over the junctions of
    correct ligations divided by every ligation the junction was seen making with an end the
    reaction holds. A junction has two ends, one presenting the top-strand overhang and one
    presenting its reverse complement, and both are counted — which is what reproduces the
    fidelity the paper reports for its own worked examples.

    A `profile` stands in where no shipped matrix covers the enzyme, and `prefer_profile` uses
    it even where one does. A profile is a measurement of the ligase and the conditions and not
    of the enzyme, so `FidelityReport.enzyme_specific` is then ``False`` and the source says so.

    With neither, the rules score the set instead and `FidelityReport.measured` is ``False``.
    Those numbers are a ranking, not a prediction: they compare one candidate set with another
    scored the same way and with nothing else.

    Raises
    ------
    ValueError
        If an overhang is not one this enzyme leaves, or a profile covers overhangs of another
        length than the enzyme leaves.

    Examples
    --------
    >>> round(fidelity(["AAAA"], "BsaI").value, 3)
    1.0
    """
    one = _type_iis(_one(enzyme))
    chosen = tuple(overhang.upper() for overhang in overhangs)
    for overhang in chosen:
        if len(overhang) != one.overhang_length or set(overhang) - set("ACGT"):
            raise ValueError(
                f"{one.name} leaves a {one.overhang_length}-base overhang, "
                f"and {overhang!r} is {len(overhang)}"
            )
    table, specific = _scoring(one, profile, prefer_profile)
    if table is None:
        return _by_rule(one, chosen)
    query = sorted({*chosen, *(reverse_complement(overhang) for overhang in chosen)})
    ligations: list[Ligation] = []
    weak: list[str] = []
    value = 1.0
    for overhang in chosen:
        partner = reverse_complement(overhang)
        ends = (overhang, partner)
        ligation = Ligation(
            overhang,
            sum(table.count(end, other) for end, other in (ends, ends[::-1])),
            sum(table.count(end, column) for end in ends for column in query),
        )
        ligations.append(ligation)
        value *= ligation.value
        if table.normalised(overhang, partner) < STRONG_LIGATION:
            weak.append(overhang)
    mismatches = sorted(
        (
            (row, column, table.normalised(row, column))
            for row in query
            for column in query
            if column != reverse_complement(row)
            and table.normalised(row, column) >= MODEST_MISMATCH
        ),
        key=lambda entry: (-entry[2], entry[0], entry[1]),
    )
    source = table.source
    if not specific:
        source = f"{source}, a ligase profile and not a measurement of {one.name}"
    return FidelityReport(
        one.name,
        source,
        measured=True,
        value=value,
        ligations=tuple(ligations),
        weak=tuple(weak),
        mismatches=tuple(mismatches),
        enzyme_specific=specific,
    )


def _by_rule(enzyme: Enzyme, chosen: Sequence[str]) -> FidelityReport:
    """Score a set by the rules, for an enzyme nobody has measured."""
    value = 1.0
    for overhang in chosen:
        uniform = set(overhang) <= set("GC") or set(overhang) <= set("AT")
        penalty = _UNIFORM_PENALTY if uniform else 0.0
        partners = [
            partner
            for other in chosen
            if other != overhang
            for partner in (other, reverse_complement(other))
        ]
        penalty += _NEAR_PENALTY * sum(
            1 for partner in partners if _distance(overhang, partner) < MIN_DISTANCE
        )
        value *= max(0.0, 1.0 - penalty)
    return FidelityReport(enzyme.name, _RULE_SOURCE, measured=False, value=value)


def _scoring(
    enzyme: Enzyme, profile: LigaseProfile | None, prefer_profile: bool
) -> tuple[Scoring | None, bool]:
    """Return what scores this enzyme's overhangs, and whether the enzyme itself was measured.

    The enzyme's own matrix wins unless the caller asks for the profile, because a ligase
    profile stands in for a measurement nobody has made rather than replacing one they have.

    Raises
    ------
    ValueError
        If the profile covers overhangs of another length than the enzyme leaves.
    """
    matrix = ligation_matrix(enzyme)
    if profile is not None and (prefer_profile or matrix is None):
        if profile.overhang_length != enzyme.overhang_length:
            raise ValueError(
                f"{profile.path.name} covers {profile.overhang_length}-base overhangs and "
                f"{enzyme.name} leaves {enzyme.overhang_length}"
            )
        return profile, False
    return matrix, matrix is not None


def _one(enzyme: EnzymeLike) -> Enzyme:
    """Read one enzyme, by name or by record."""
    return get_enzyme(enzyme) if isinstance(enzyme, str) else enzyme


def _resolve(enzymes: Iterable[EnzymeLike]) -> tuple[Enzyme, ...]:
    """Read several enzymes, by name or by record."""
    return tuple(_one(one) for one in enzymes)


def _type_iis(enzyme: Enzyme) -> Enzyme:
    """Refuse an enzyme that cuts inside its own site, which leaves no overhang to design.

    Raises
    ------
    ValueError
        If the enzyme is not Type IIS.
    """
    if enzyme.type != "IIS":
        raise ValueError(
            f"{enzyme.name} cuts inside its own site, so its overhang is not one to design"
        )
    return enzyme

"""The rules a method's cut ends are held to, and how well a set of overhangs should ligate.

Every method that cuts DNA and joins it again asks the first question, and every method cutting
with a Type IIS enzyme asks the other two, so all three sit here below the pipelines rather than
inside one of them -- `docs/adr/0007-cloning-methods.md` says why:

- **Whether two ends anneal.** An end is its overhang and its end type together, which is what
  `End` holds and `compatible` weighs. The bases alone do not decide it: a 5' overhang and a 3'
  overhang spelling the same bases run the wrong way for each other.
- **Which overhangs a set may hold.** A palindrome would let a fragment ligate to itself, a
  repeat would let two junctions swap, and a near-duplicate is what mis-ligates. An overhang a
  step outside the set already spends is reserved, and held out by name. `refusal` weighs one
  candidate against the set so far and names the rule that refuses it.
- **How well the set should ligate.** Pryor et al. 2020 measured every overhang pair for five
  enzymes, and that data ships here: fidelity is the product over the junctions of correct
  ligations over all ligations. An enzyme they did not measure is scored on a shipped matrix of
  the same overhang length standing in, then on a ligase profile the caller holds
  (`liulab_mbio.ligase`), then by the rules, and the report says which of the four scored it
  and whether it is specific to the enzyme.

The fidelity data is `src/liulab_mbio/data/ligation_fidelity.json`, from the supplementary
tables of Pryor et al. 2020, used under CC BY 4.0. It carries its own citation, which
`ligation_source` reads for anything printing the paper in full.
`docs/research/ligation-fidelity.md` records the licence, the axis convention and the rules.
"""

import json
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import KW_ONLY, dataclass, field
from functools import cache
from importlib.resources import files
from typing import Any, Literal

from liulab_mbio.enzymes import EndType, Enzyme, get_enzyme
from liulab_mbio.ligase import LigaseProfile
from liulab_mbio.sequence import SequenceRecord, reverse_complement
from liulab_mbio.sites import EnzymeLike, primer_tail

#: How many bases two overhangs in one set must differ by, counting each one's reverse complement
#: too. Two is read off the shipped matrices rather than taken from a standard: it is the smallest
#: separation at which every one of them holds its measured cross-ligations below
#: `MODEST_MISMATCH`, where a one-base separation leaves a third to a half of its pairs at or
#: above it. It is an argument because a caller may hold a set to more.
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

#: Why a candidate overhang was refused. `refusal` returns every one but ``"stop"``, which belongs
#: to a caller reading a candidate in frame, as the library pipeline's overhang choice does.
type RejectionRule = Literal[
    "length", "palindrome", "uniform", "repeat", "reserved", "near-duplicate", "site", "stop"
]

#: What a set of overhangs is scored against: the enzyme's own matrix, or a ligase's profile.
type Scoring = LigationMatrix | LigaseProfile


@dataclass(frozen=True, slots=True)
class End:
    """One cut end, as the rule joining two of them sees it.

    Parameters
    ----------
    overhang
        The bases the cut left single-stranded, written as the top strand reads them 5' to 3',
        and ``""`` for a blunt end. `liulab_mbio.sites.Fragment` writes both of its ends that
        way, so two ends that meet are compared as they are written. Stored upper-case.
    end_type
        Which strand carries those bases, or ``"blunt"``. `cut_by` reads it off the enzyme.

    Raises
    ------
    ValueError
        If a blunt end carries bases, or an overhanging one carries none.
    """

    overhang: str
    end_type: EndType

    def __post_init__(self) -> None:
        """Upper-case the overhang and refuse an end type its own bases contradict."""
        object.__setattr__(self, "overhang", self.overhang.upper())
        if bool(self.overhang) == (self.end_type == "blunt"):
            raise ValueError(
                f"a {self.end_type} end cannot leave {self.overhang!r} single-stranded"
            )

    @classmethod
    def cut_by(cls, enzyme: EnzymeLike, overhang: str) -> "End":
        """Return the end `enzyme` leaves where its cut left `overhang` single-stranded.

        Raises
        ------
        ValueError
            If the overhang is not as long as the one this enzyme leaves.

        Examples
        --------
        >>> End.cut_by("EcoRI", "AATT")
        End(overhang='AATT', end_type="5'")
        """
        one = _one(enzyme)
        if len(overhang) != one.overhang_length:
            raise ValueError(
                f"{one.name} leaves a {one.overhang_length}-base overhang, "
                f"and {overhang!r} is {len(overhang)}"
            )
        return cls(overhang, one.end)


def compatible(one: End, other: End) -> bool:
    """Whether two cut ends anneal, so a ligase can seal them into one molecule.

    Two blunt ends anneal. Two overhangs anneal when they spell the same bases and sit on the
    same strand, which is why an end is its overhang and its end type together.

    Examples
    --------
    >>> compatible(End.cut_by("SalI", "TCGA"), End.cut_by("XhoI", "TCGA"))
    True
    >>> compatible(End.cut_by("HindIII", "AGCT"), End.cut_by("SacI", "AGCT"))
    False
    """
    return one.end_type == other.end_type and one.overhang == other.overhang


@dataclass(frozen=True, slots=True)
class LigationSource:
    """The paper the shipped ligation data was measured in.

    Parameters
    ----------
    citation
        The paper in full, as a reference list prints it.
    doi
        Its DOI, bare.
    """

    citation: str
    doi: str

    @property
    def url(self) -> str:
        """Where the paper resolves."""
        return f"https://doi.org/{self.doi}"


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
    def cited(self) -> str:
        """The paper by first author and year, or the whole citation where it gives neither."""
        # The data ships the citation whole and not in parts, so the short form is read off it.
        surname, _, rest = self.citation.partition(",")
        year = re.search(r"\((\d{4})\)", rest)
        return f"{surname} {year.group(1)}" if surname and year else self.citation

    @property
    def source(self) -> str:
        """The paper in passing and which of its tables, as a check cites it.

        The reference in full belongs under References, not in a sentence read at the bench.
        """
        return f"{self.cited} {self.table}, measured with {self.product}"

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
        ``False`` when the number is not a measurement of this enzyme: a ligase profile, which
        measures the ligase and the conditions, or another Type IIS enzyme's matrix standing in.
    stand_in
        The product whose matrix stood in, where one did. Empty otherwise, which is what keeps
        a stand-in apart from a ligase profile in `label`.
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
    stand_in: str = ""

    @property
    def label(self) -> str:
        """What kind of number this is, for a report printing it beside the value."""
        if not self.measured:
            return "rule-based estimate"
        if self.enzyme_specific:
            return "measured"
        if self.stand_in:
            return f"measured with {self.stand_in}, not specific to {self.enzyme}"
        return f"measured ligase profile, not specific to {self.enzyme}"


@dataclass(frozen=True, slots=True)
class OnTarget:
    """How often one ligase was seen joining one overhang to its own partner.

    Parameters
    ----------
    overhang
        The overhang, written on the top strand.
    rate
        Its correct Watson-Crick pair, per 100,000 ligation events.
    floor
        The rate it is held to.
    """

    overhang: str
    rate: float
    floor: float

    @property
    def weak(self) -> bool:
        """Whether this ligase joins it below the floor."""
        return self.rate < self.floor


@dataclass(frozen=True, slots=True)
class OnTargetReport:
    """How often one ligase joins each overhang of a set, and where the numbers came from.

    Fidelity asks whether a set's junctions can be told apart. This asks how often each one is
    made at all, which is a different measurement of a different thing: a rate below the floor
    costs correct joins, so it costs colonies rather than giving the wrong product. Nothing here
    is a score, and nothing ranks on it.

    Parameters
    ----------
    source
        The conditions the profile was measured under and the file it was read from.
    floor
        The rate each overhang was held to.
    rates
        One entry per overhang, in the order given.
    """

    source: str
    _: KW_ONLY
    floor: float
    rates: tuple[OnTarget, ...] = ()

    @property
    def weak(self) -> tuple[str, ...]:
        """The overhangs this ligase joins below the floor, the worst first."""
        return tuple(
            one.overhang for one in sorted(self.rates, key=lambda one: one.rate) if one.weak
        )

    @property
    def label(self) -> str:
        """Which measurement this is, for a report printing it beside a fidelity score."""
        return f"measured on {self.source}"


@dataclass(frozen=True, slots=True)
class Junction:
    """A junction as it is asked for, before anything is cut: where two parts are to meet.

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


@cache
def _document() -> Mapping[str, Any]:
    """Read the shipped ligation data file whole."""
    text = (files("liulab_mbio") / "data" / "ligation_fidelity.json").read_text(encoding="utf-8")
    return json.loads(text)


def ligation_source() -> LigationSource:
    """Return the paper the shipped ligation data was measured in.

    The data file is the citation's one home, so anything printing the paper -- a reference
    list, a check naming it in passing -- reads it from here.

    Examples
    --------
    >>> ligation_source().url
    'https://doi.org/10.1371/journal.pone.0238592'
    """
    source = _document()["source"]
    return LigationSource(source["citation"], source["doi"])


@cache
def _shipped() -> Mapping[str, LigationMatrix]:
    """Read the shipped ligation data, keyed by the enzyme it was measured with."""
    document = _document()
    citation = ligation_source().citation
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


def refusal(
    candidate: str,
    enzyme: EnzymeLike,
    *,
    taken: Iterable[str] = (),
    reserved: Iterable[str] = (),
    avoid: Iterable[EnzymeLike] = (),
    min_distance: int = MIN_DISTANCE,
    allow_uniform: bool = False,
) -> Rejection | None:
    """Why `candidate` will not join `taken`, or ``None`` when it will.

    Every rule a set of overhangs is held to, weighed one candidate at a time: its length, a
    palindrome, one base class, a repeat, a reserved overhang, a near-duplicate, and a primer
    tail spelling a further site.

    `taken` is what this set already holds; `reserved` is what a step outside it spends, which
    the set is to leave alone. Both bind the candidate the same way — it may neither take one
    nor sit nearer than `min_distance` to one — and taking a reserved one is refused as
    ``"reserved"`` rather than as a repeat, so a caller is told what actually happened.

    Examples
    --------
    >>> refusal("AGGT", "BsaI") is None
    True
    >>> refusal("AGGT", "BsaI", taken=["AGGT"]).rule
    'repeat'
    >>> refusal("AGGT", "BsaI", reserved=["AGGT"]).rule
    'reserved'
    """
    return _refuse(
        candidate.upper(),
        _type_iis(_one(enzyme)),
        tuple(one.upper() for one in taken),
        tuple(one.upper() for one in reserved),
        _resolve(avoid),
        min_distance,
        allow_uniform,
    )


def _refuse(
    candidate: str,
    enzyme: Enzyme,
    taken: Sequence[str],
    reserved: Sequence[str],
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
    for other in reserved:
        if candidate in (other, reverse_complement(other)):
            return Rejection(
                candidate,
                "reserved",
                f"{other} is reserved, and a junction taking it would ligate to whatever "
                "reserved it",
            )
    for other in (*taken, *reserved):
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

    Where no shipped matrix covers the enzyme, `stand_in_matrix` scores it instead, and a
    `profile` the caller holds comes after that; `prefer_profile` lifts the profile above
    everything. Neither is a measurement of this enzyme, so `FidelityReport.enzyme_specific` is
    then ``False`` and both the source and the label say which stood in.

    With none of them, the rules score the set instead and `FidelityReport.measured` is ``False``.
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
    table, specific = scoring(one, profile=profile, prefer_profile=prefer_profile)
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
    stood_in = ""
    if not specific and isinstance(table, LigationMatrix):
        stood_in = table.product
        source = f"{source}, standing in for {one.name}, which nobody has measured"
    elif not specific:
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
        stand_in=stood_in,
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


def on_target(
    overhangs: Iterable[str], profile: LigaseProfile, *, floor: float = STRONG_LIGATION
) -> OnTargetReport:
    """Report how often the ligase `profile` measured joined each overhang to its own partner.

    A set is scored for fidelity on the enzyme's own matrix, which is T4's chemistry. A method
    ligating with another ligase joins the same overhangs at its own rates, and the two disagree
    most on the A/T-rich ones. This reads one cell an overhang and says which fall below `floor`,
    so a caller holding a profile sees what the score could not. It returns no score, because
    what is at risk is how many correct joins are made and not which product they make.

    `floor` defaults to `STRONG_LIGATION`, the rate NEB's Viewer calls a Watson-Crick pair strong
    at; `docs/research/ligation-fidelity.md` section 8 holds what it warns on.

    Raises
    ------
    ValueError
        If an overhang is not as long as the ones the profile covers.
    """
    chosen = tuple(overhang.upper() for overhang in overhangs)
    for overhang in chosen:
        if len(overhang) != profile.overhang_length or set(overhang) - set("ACGT"):
            raise ValueError(
                f"{profile.source} covers {profile.overhang_length}-base overhangs, "
                f"and {overhang!r} is {len(overhang)}"
            )
    return OnTargetReport(
        profile.source,
        floor=floor,
        rates=tuple(
            OnTarget(overhang, profile.normalised(overhang, reverse_complement(overhang)), floor)
            for overhang in chosen
        ),
    )


def best_overhang(
    junction: Junction,
    candidates: Iterable[tuple[int, str]],
    enzyme: EnzymeLike,
    *,
    taken: Sequence[str] = (),
    reserved: Iterable[str] = (),
    avoid: Iterable[EnzymeLike] = (),
    min_distance: int = MIN_DISTANCE,
    allow_uniform: bool = False,
    profile: LigaseProfile | None = None,
    prefer_profile: bool = False,
) -> tuple[Choice | None, tuple[Rejection, ...]]:
    """Take the candidate the whole set scores best with, and say what every rule refused.

    A candidate is ranked by `fidelity` over `taken` and itself together, not by its own
    on-target count: what a junction costs is what it does to the junctions beside it, and a
    candidate that scores well alone can be the one that mis-ligates. A tie goes to the
    candidate nearest the position asked for, then to the bases, so the answer does not depend
    on the order the candidates arrive in.

    Parameters
    ----------
    junction
        What the choice is for, which the answer carries.
    candidates
        How far each candidate moves the junction, and the overhang it would leave.
    enzyme
        The enzyme cutting it, which sets the overhang length.
    taken
        The overhangs the set already holds.
    reserved, avoid, min_distance, allow_uniform
        What `refusal` holds each candidate to.
    profile, prefer_profile
        A ligase's own matrix, and whether to score with it where shipped data covers the
        enzyme.

    Returns
    -------
    The choice, or ``None`` where every candidate was refused, and the refusals either way.

    Examples
    --------
    >>> choice, refused = best_overhang(Junction("one"), [(0, "AGGT"), (0, "ACGT")], "BsaI")
    >>> choice.overhang, refused[0].rule
    ('AGGT', 'palindrome')
    """
    one = _type_iis(_one(enzyme))
    others = _resolve(avoid)
    held = tuple(taken)
    rejected: list[Rejection] = []
    best: tuple[float, int, str] | None = None
    offset = 0
    for moved, candidate in candidates:
        refused = refusal(
            candidate,
            one,
            taken=held,
            reserved=reserved,
            avoid=others,
            min_distance=min_distance,
            allow_uniform=allow_uniform,
        )
        if refused is not None:
            rejected.append(refused)
            continue
        scored = fidelity(
            [*held, candidate.upper()], one, profile=profile, prefer_profile=prefer_profile
        )
        key = (-scored.value, abs(moved), candidate.upper())
        if best is None or key < best:
            best, offset = key, moved
    if best is None:
        return None, tuple(rejected)
    return Choice(junction, best[2], offset, tuple(rejected)), tuple(rejected)


def stand_in_matrix(enzyme: EnzymeLike) -> LigationMatrix | None:
    """Return the shipped matrix that stands in for an enzyme nobody measured, or ``None``.

    An isoschizomer first: an enzyme reading the same site cuts it the same way, so a matrix
    measured with one is a measurement of the other under another name. Otherwise the shipped
    matrix of the same overhang length with the most ligations behind it, which Pryor 2020's
    Discussion is what licenses: the predicted fidelity "is unlikely to be significantly
    impacted by the choice of Type IIS restriction enzyme". Both are read off the shipped data
    rather than listed, so a matrix added later stands in without a table to edit.

    An enzyme with its own matrix needs no stand-in and gets ``None``.

    Examples
    --------
    >>> stand_in_matrix("PaqCI").enzyme
    'Esp3I'
    >>> stand_in_matrix("BspQI").enzyme
    'SapI'
    >>> stand_in_matrix("BsaI") is None
    True
    """
    one = _one(enzyme)
    if ligation_matrix(one) is not None:
        return None
    shipped = [(matrix, get_enzyme(matrix.enzyme)) for matrix in _shipped().values()]
    same_site = [matrix for matrix, other in shipped if other.site == one.site]
    if same_site:
        return min(same_site, key=lambda matrix: matrix.enzyme)
    length = [matrix for matrix, other in shipped if other.overhang_length == one.overhang_length]
    if not length:
        return None
    return max(length, key=lambda matrix: (matrix.observations, matrix.enzyme))


def scoring(
    enzyme: EnzymeLike,
    *,
    profile: LigaseProfile | None = None,
    prefer_profile: bool = False,
) -> tuple[Scoring | None, bool]:
    """Return what scores this enzyme's overhangs, and whether the enzyme itself was measured.

    The enzyme's own matrix wins unless the caller asks for the profile. Where nobody measured
    the enzyme, `stand_in_matrix` is next: a matrix of the same reaction run with another Type
    IIS enzyme is a smaller error than one of a different reaction, measured either way. A
    profile is the third fallback, for an enzyme no shipped matrix shares an overhang length
    with, and `prefer_profile` lifts it above everything. None of the three: ``None``, and the
    rules score the set.

    Raises
    ------
    ValueError
        If the profile covers overhangs of another length than the enzyme leaves.

    Examples
    --------
    >>> table, specific = scoring("BsaI")
    >>> table.enzyme, specific
    ('BsaI', True)
    >>> table, specific = scoring("PaqCI")
    >>> table.enzyme, specific
    ('Esp3I', False)
    """
    one = _one(enzyme)
    if profile is not None and prefer_profile:
        return _fitting(profile, one), False
    matrix = ligation_matrix(one)
    if matrix is not None:
        return matrix, True
    stand_in = stand_in_matrix(one)
    if stand_in is not None:
        return stand_in, False
    if profile is not None:
        return _fitting(profile, one), False
    return None, False


def _fitting(profile: LigaseProfile, enzyme: Enzyme) -> LigaseProfile:
    """Return the profile, refusing one whose overhangs this enzyme could not leave.

    Raises
    ------
    ValueError
        If the profile covers overhangs of another length than the enzyme leaves.
    """
    if profile.overhang_length != enzyme.overhang_length:
        raise ValueError(
            f"{profile.path.name} covers {profile.overhang_length}-base overhangs and "
            f"{enzyme.name} leaves {enzyme.overhang_length}"
        )
    return profile


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

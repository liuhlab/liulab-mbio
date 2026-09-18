"""The rules a set of Type IIS overhangs is held to, and how well the set should ligate.

Any method that cuts with a Type IIS enzyme asks the same two questions, so they sit here below
the pipelines rather than inside one of them:

- **Which overhangs a set may hold.** A palindrome would let a fragment ligate to itself, a
  repeat would let two junctions swap, and a near-duplicate is what mis-ligates. `refusal`
  weighs one candidate against the set so far and names the rule that refuses it.
- **How well the set should ligate.** Pryor et al. 2020 measured every overhang pair for five
  enzymes, and that data ships here: fidelity is the product over the junctions of correct
  ligations over all ligations. An enzyme they did not measure is scored against a ligase
  profile the caller holds (`liulab_mbio.ligase`) where there is one and by the rules where
  there is not, and the report says which of the three scored it.

The fidelity data is `src/liulab_mbio/data/ligation_fidelity.json`, from the supplementary
tables of Pryor, J.M., Potapov, V., Kucera, R.B., Bilotti, K., Cantor, E.J. and Lohman, G.J.S.
(2020) *PLoS One* 15(9): e0238592, used under CC BY 4.0.
`docs/research/ligation-fidelity.md` records the licence, the axis convention and the rules.
"""

import json
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import KW_ONLY, dataclass, field
from functools import cache
from importlib.resources import files
from typing import Literal

from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.ligase import LigaseProfile
from liulab_mbio.sequence import SequenceRecord, reverse_complement
from liulab_mbio.sites import EnzymeLike, primer_tail

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

#: Why a candidate overhang was refused. `refusal` returns every one but ``"stop"``, which belongs
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

    Every rule a set of overhangs is held to, weighed one candidate at a time: its length, a
    palindrome, one base class, a repeat, a near-duplicate, and a primer tail spelling a further
    site.

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


def scoring(
    enzyme: EnzymeLike,
    *,
    profile: LigaseProfile | None = None,
    prefer_profile: bool = False,
) -> tuple[Scoring | None, bool]:
    """Return what scores this enzyme's overhangs, and whether the enzyme itself was measured.

    The enzyme's own matrix wins unless the caller asks for the profile, because a ligase
    profile stands in for a measurement nobody has made rather than replacing one they have.
    Neither: ``None``, and the rules score the set.

    Raises
    ------
    ValueError
        If the profile covers overhangs of another length than the enzyme leaves.

    Examples
    --------
    >>> table, specific = scoring("BsaI")
    >>> table.enzyme, specific
    ('BsaI', True)
    """
    one = _one(enzyme)
    matrix = ligation_matrix(one)
    if profile is not None and (prefer_profile or matrix is None):
        if profile.overhang_length != one.overhang_length:
            raise ValueError(
                f"{profile.path.name} covers {profile.overhang_length}-base overhangs and "
                f"{one.name} leaves {one.overhang_length}"
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

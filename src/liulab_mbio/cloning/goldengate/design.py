"""Choosing the Type IIS enzyme an assembly uses, and the overhangs its junctions are cut to.

Two decisions, in the order a design makes them:

- **Which enzyme.** Rank the candidates by how many sites they read in the parts that end up in
  the product. An enzyme with no site needs nothing done to the parts, so domestication is
  proposed only when no candidate is free, and then as a report of what it would change.
- **Which overhangs.** Each junction takes an overhang of the length the enzyme leaves: the
  first candidate no rule of `liulab_mbio.overhangs` refuses, and a junction free to move
  slides within its window to reach one. The chosen set carries the fidelity scored there.
"""

from collections.abc import Iterable, Iterator, Sequence
from dataclasses import KW_ONLY, dataclass
from itertools import product

from liulab_mbio.codons import CodonUsage
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.ligase import LigaseProfile
from liulab_mbio.overhangs import (
    MIN_DISTANCE,
    Choice,
    FidelityReport,
    Junction,
    Rejection,
    Scoring,
    fidelity,
    ligation_matrix,
    refusal,
    scoring,
)
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand, reverse_complement
from liulab_mbio.sites import (
    CutSite,
    Domestication,
    EnzymeLike,
    domesticate,
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
        A ligase's own matrix, read by `liulab_mbio.ligase.read_profile`. It ranks the free
        candidates and scores the set where no shipped matrix covers the enzyme.
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
    table, _ = scoring(one, profile=profile, prefer_profile=prefer_profile)
    wanted = tuple(junctions)
    order = sorted(range(len(wanted)), key=lambda index: _freedom(wanted[index]))
    chosen: dict[int, Choice] = {}
    taken: list[str] = []
    for index in order:
        junction = wanted[index]
        rejected: list[Rejection] = []
        for offset, candidate in _candidates(junction, one, table):
            refused = refusal(
                candidate,
                one,
                taken=taken,
                avoid=others,
                min_distance=min_distance,
                allow_uniform=allow_uniform,
            )
            if refused is None:
                chosen[index] = Choice(junction, candidate, offset, tuple(rejected))
                taken.append(candidate)
                break
            rejected.append(refused)
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

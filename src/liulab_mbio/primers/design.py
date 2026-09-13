"""Design a primer, or a pair, by choosing the annealing region a placement allows."""

import heapq
import math
from collections.abc import Iterator
from dataclasses import dataclass
from itertools import groupby

from liulab_mbio.checks import STATUSES, Check, Status, worst
from liulab_mbio.primers.evaluation import PrimerReport, evaluate_primer, pair_checks
from liulab_mbio.primers.placement import Placement
from liulab_mbio.primers.polymerase import Q5, Polymerase
from liulab_mbio.primers.thresholds import TARGET_TM, THRESHOLDS, Thresholds
from liulab_mbio.sequence import (
    BindingSite,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)

#: How a candidate ranks: its worst status, how many checks did not pass, whether its length
#: lies outside the band a length passes on, how far its Tm sits from the target, how far its
#: 5' end sits from the position asked for, its length, and that end, which settles a tie.
type _Rank = tuple[int, int, int, float, int, int, int]

#: A pair's rank, or the best one a pair could still reach: the same keys, summed.
type _PairRank = tuple[float, ...]


def design_primer(
    template: SequenceRecord,
    position: int,
    strand: Strand,
    *,
    placement: Placement | None = None,
    tail: str = "",
    name: str = "",
    target_tm: float = TARGET_TM,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
) -> Primer:
    """Design a primer annealing at a position on a template, or anywhere a `Placement` allows.

    A forward primer's annealing region starts at `position` and reads towards higher
    coordinates; a reverse primer's ends there and reads back. Either crosses the origin of a
    circular template. `placement` widens the search to every binding site it allows, and
    `position` stays the 5' end asked for. Each candidate, at every length `Thresholds.length`
    does not fail, is judged whole by `evaluate_primer` on the template, structures, runs and
    off-target sites included, so a designed primer is worse than a pass only where nothing the
    placement allows passes. Candidates rank by worst status, then by how many checks did not
    pass, then by the band a length passes on, then by the Tm nearest `target_tm`, then by
    nearness to `position`, and the shorter primer breaks a tie. `tail` joins its 5' end and
    stays out of the binding site.

    Raises
    ------
    ValueError
        If no annealing region fits the template there, or fits the placement.

    Examples
    --------
    >>> template = SequenceRecord(
    ...     "GGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCT"
    ...     "TACAGGATCCACTAGTAACGGCCGCCAGTGTGCTGGAATTCGCCCTTA"
    ... )
    >>> design_primer(template, 0, Strand.FORWARD, name="MCS fwd").sequence
    'GGCGTAATCATGGTCATAGC'

    Anchored, where a tail has to join the 5' end at one place and only the length may vary:

    >>> anchored = Placement(five_prime=Segment(48, 49))
    >>> design_primer(template, 48, Strand.FORWARD, placement=anchored, tail="TTGGTCTCA").sequence
    'TTGGTCTCATACAGGATCCACTAGTAACGG'

    Near a target, where the 3' end has to land 10 to 21 bases before a junction:

    >>> junction = 48
    >>> near = Placement(three_prime=Segment(junction - 21, junction - 9))
    >>> primer = design_primer(template, junction - 40, Strand.FORWARD, placement=near)
    >>> primer.sequence, junction - primer.binding_sites[0].end
    ('ATCATGGTCATAGCTGTTTCC', 21)

    In a region, where both ends are free inside a flank:

    >>> region = Placement(five_prime=Segment(0, 6), three_prime=Segment(20, 40))
    >>> design_primer(template, 0, Strand.FORWARD, placement=region).sequence
    'CGTAATCATGGTCATAGCTGTT'
    """
    options = _options(
        template, position, strand, placement, tail, name, target_tm, polymerase, thresholds
    )
    return min(_fitted(options, position, placement), key=lambda option: option.rank).primer


def design_pair(
    template: SequenceRecord,
    start: int,
    end: int,
    *,
    forward_placement: Placement | None = None,
    reverse_placement: Placement | None = None,
    forward_tail: str = "",
    reverse_tail: str = "",
    forward_name: str = "",
    reverse_name: str = "",
    target_tm: float = TARGET_TM,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
) -> tuple[Primer, Primer]:
    """Design a pair amplifying `start` to `end`, judged on what a pair adds as well.

    `end` passes the length of a circular template when the amplicon crosses the origin. Each
    candidate is judged as in `design_primer`, and every pair of them adds `pair_checks`: the
    gap between the two Tms, the heterodimers and the products. The pair chosen ranks best on
    its worst status, then on the same keys summed over both primers, so a pair worse than a
    pass means nothing here passes. A placement widens either primer's search; the two bound
    the amplicon between them, which is why nothing takes an amplicon size.

    Raises
    ------
    ValueError
        If no annealing region fits the template at either end, or fits that end's placement.

    Examples
    --------
    >>> template = SequenceRecord(
    ...     "GGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCT"
    ...     "TACAGGATCCACTAGTAACGGCCGCCAGTGTGCTGGAATTCGCCCTTA"
    ... )
    >>> [primer.sequence for primer in design_pair(template, 0, len(template))]
    ['GGCGTAATCATGGTCATAGC', 'TAAGGGCGAATTCCAGCA']

    In a region, where each primer is free inside its own flank and the two bound the amplicon:

    >>> forward = Placement(five_prime=Segment(0, 6), three_prime=Segment(15, 35))
    >>> reverse = Placement(five_prime=Segment(90, 96), three_prime=Segment(60, 80))
    >>> pair = design_pair(
    ...     template,
    ...     0,
    ...     len(template),
    ...     forward_placement=forward,
    ...     reverse_placement=reverse,
    ... )
    >>> [primer.sequence for primer in pair]
    ['CGTAATCATGGTCATAGCTGTT', 'CGAATTCCAGCACACTGG']
    """
    forwards = _options(
        template,
        start,
        Strand.FORWARD,
        forward_placement,
        forward_tail,
        forward_name,
        target_tm,
        polymerase,
        thresholds,
    )
    reverses = _options(
        template,
        end,
        Strand.REVERSE,
        reverse_placement,
        reverse_tail,
        reverse_name,
        target_tm,
        polymerase,
        thresholds,
    )
    return _best_pair(
        _fitted(forwards, start, forward_placement),
        _fitted(reverses, end, reverse_placement),
        template,
        thresholds,
    )


@dataclass(frozen=True, slots=True)
class _Option:
    """One candidate binding site: the primer it makes, what it scored, how it ranks."""

    primer: Primer
    report: PrimerReport
    rank: _Rank


def _options(
    template: SequenceRecord,
    position: int,
    strand: Strand,
    placement: Placement | None,
    tail: str,
    name: str,
    target_tm: float,
    polymerase: Polymerase,
    thresholds: Thresholds,
) -> list[_Option]:
    """Return every binding site the placement allows, each judged whole and ranked."""
    options = []
    for site in _sites(template, position, strand, placement, thresholds):
        bases = template.extract(Segment(site.start, site.end))
        sequence = bases if strand is Strand.FORWARD else reverse_complement(bases)
        primer = Primer(name, tail + sequence, binding_sites=(site,))
        report = evaluate_primer(primer, template, polymerase=polymerase, thresholds=thresholds)
        five = site.start if strand is Strand.FORWARD else site.end
        options.append(_Option(primer, report, _rank(report, target_tm, five, position, template)))
    return options


def _sites(
    template: SequenceRecord,
    position: int,
    strand: Strand,
    placement: Placement | None,
    thresholds: Thresholds,
) -> Iterator[BindingSite]:
    """Yield every binding site a placement allows, at every length the band does not fail."""
    length = len(template)
    circular = template.topology == "circular"
    shortest = int(_finite(thresholds.length.warn_low, thresholds.length.low))
    longest = min(int(_finite(thresholds.length.warn_high, thresholds.length.high)), length)
    anchors, five = _anchors(template, position, placement)
    for anchor in anchors:
        for size in range(shortest, longest + 1):
            # An anchor is one end: the start of the site where it is the 5' end of a forward
            # primer or the 3' end of a reverse one, and a length back from it otherwise.
            start = anchor if (strand is Strand.FORWARD) == five else anchor - size
            if circular:
                start %= length
            elif not 0 <= start <= length - size:
                continue
            site = BindingSite(start, start + size, strand)
            if placement is None or placement.allows(site, template):
                yield site


def _anchors(
    template: SequenceRecord, position: int, placement: Placement | None
) -> tuple[list[int], bool]:
    """Return the positions the end a placement bounds may take, and whether that end is 5'."""
    if placement is None:
        return [position], True
    five = _positions(placement.five_prime, template)
    return (five, True) if five else (_positions(placement.three_prime, template), False)


def _positions(span: Segment | None, template: SequenceRecord) -> list[int]:
    """Return every position a span of them holds, wrapping one across the origin."""
    if span is None:
        return []
    length = len(template)
    circular = template.topology == "circular"
    return [
        (span.start + step) % length if circular else span.start + step
        for step in range(span.end - span.start)
    ]


def _fitted(options: list[_Option], position: int, placement: Placement | None) -> list[_Option]:
    """Return the options, refusing a position or a placement holding no annealing region."""
    if not options:
        where = (
            f"this template at {position}" if placement is None else f"this placement: {placement}"
        )
        raise ValueError(f"no annealing region fits {where}")
    return options


def _rank(
    report: PrimerReport, target_tm: float, five: int, position: int, template: SequenceRecord
) -> _Rank:
    return (
        STATUSES.index(report.status),
        sum(check.status not in (None, "pass") for check in report.checks),
        int(report["length"].status != "pass"),
        abs(report["tm"].value - target_tm),
        _apart(five, position, template),
        int(report["length"].value),
        five % len(template),
    )


def _apart(position: int, wanted: int, template: SequenceRecord) -> int:
    """Return how far a 5' end lies from the position asked for, the short way round a circle."""
    if template.topology != "circular":
        return abs(position - wanted)
    gap = (position - wanted) % len(template)
    return min(gap, len(template) - gap)


def _best_pair(
    forwards: list[_Option],
    reverses: list[_Option],
    template: SequenceRecord,
    thresholds: Thresholds,
) -> tuple[Primer, Primer]:
    """Return the pair ranking best, judging only pairs that can still reach the best rank."""
    pairs = _by_bound(_ranked(forwards), _ranked(reverses))
    bound, forward, reverse = next(pairs)
    best = (_judged(bound, forward, reverse, template, thresholds), forward, reverse)
    for bound, forward, reverse in pairs:
        if bound >= best[0]:
            break
        rank = _judged(bound, forward, reverse, template, thresholds)
        if rank < best[0]:
            best = (rank, forward, reverse)
    return best[1].primer, best[2].primer


def _judged(
    bound: _PairRank,
    forward: _Option,
    reverse: _Option,
    template: SequenceRecord,
    thresholds: Thresholds,
) -> _PairRank:
    """Return how a pair ranks once what only a pair can be judged on has been judged."""
    checks = pair_checks(forward.report, reverse.report, template, thresholds=thresholds)
    return _pair_rank(bound, forward, reverse, checks)


def _pair_rank(
    bound: _PairRank, forward: _Option, reverse: _Option, checks: tuple[Check, ...]
) -> _PairRank:
    verdicts: list[Status] = [check.status for check in checks if check.status is not None]
    status = worst((forward.report.status, reverse.report.status, *verdicts))
    return (
        STATUSES.index(status),
        bound[1] + sum(one != "pass" for one in verdicts),
        *bound[2:],
        forward.rank[-1],
    )


def _ranked(options: list[_Option]) -> list[_Option]:
    """Return the options in ranked order, which is the order a pair search walks them in."""
    return sorted(options, key=lambda option: option.rank)


def _by_bound(
    forwards: list[_Option], reverses: list[_Option]
) -> Iterator[tuple[_PairRank, _Option, _Option]]:
    """Yield every pair of ranked options, lowest bound first, judging none of them."""
    groups = [
        _summed(ones, others) for ones in _by_status(forwards) for others in _by_status(reverses)
    ]
    return heapq.merge(*groups, key=lambda entry: entry[0])


def _by_status(options: list[_Option]) -> list[list[_Option]]:
    """Split ranked options into runs of one status, inside which the rest of a rank rises."""
    return [list(group) for _, group in groupby(options, key=lambda option: option.rank[0])]


def _summed(
    forwards: list[_Option], reverses: list[_Option]
) -> Iterator[tuple[_PairRank, _Option, _Option]]:
    """Yield pairs of one status each, in rising order of their bound.

    The bounds of a run rise with either index, so the next pair is always a step from one
    already taken: the walk of a sorted matrix, and it reaches the far corner only if asked.
    """
    heap = [(_bound(forwards[0], reverses[0]), 0, 0)]
    seen = {(0, 0)}
    while heap:
        bound, one, other = heapq.heappop(heap)
        yield bound, forwards[one], reverses[other]
        for step in ((one + 1, other), (one, other + 1)):
            if step[0] < len(forwards) and step[1] < len(reverses) and step not in seen:
                seen.add(step)
                heapq.heappush(heap, (_bound(forwards[step[0]], reverses[step[1]]), *step))


def _bound(forward: _Option, reverse: _Option) -> _PairRank:
    """Return the best rank a pair could reach before what only a pair adds is judged.

    Those checks can only raise the worst status or add to the checks that did not pass, so no
    pair ranks better than this.
    """
    return (
        max(forward.rank[0], reverse.rank[0]),
        *(one + other for one, other in zip(forward.rank[1:], reverse.rank[1:], strict=True)),
    )


def _finite(value: float, fallback: float) -> float:
    return value if math.isfinite(value) else fallback

"""Design a primer, or a pair, by choosing the annealing region at a position on a template."""

import math
from dataclasses import dataclass

from liulab_mbio.checks import STATUSES, Status, worst
from liulab_mbio.primers.polymerase import Q5, Polymerase, melting_temperature
from liulab_mbio.primers.thresholds import TARGET_TM, THRESHOLDS, Thresholds
from liulab_mbio.sequence import (
    BindingSite,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)


def design_primer(
    template: SequenceRecord,
    position: int,
    strand: Strand,
    *,
    tail: str = "",
    name: str = "",
    target_tm: float = TARGET_TM,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
) -> Primer:
    """Design a primer annealing at a position on a template.

    A forward primer's annealing region starts at `position` and reads towards higher
    coordinates; a reverse primer's ends there and reads back. Either crosses the origin of a
    circular template. Of the lengths `Thresholds.length` allows, the chosen one grades best
    and then lands closest to `target_tm`. `tail` joins its 5' end and stays out of the
    binding site.

    Raises
    ------
    ValueError
        If no annealing region fits the template there.

    Examples
    --------
    >>> template = SequenceRecord("GGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCT")
    >>> design_primer(template, 0, Strand.FORWARD, name="MCS fwd").sequence
    'GGCGTAATCATGGTCATAGC'
    """
    options = _annealing_options(template, position, strand, polymerase, thresholds)
    if not options:
        raise ValueError(f"no annealing region fits this template at {position}")
    best = min(options, key=lambda option: _option_score(option, target_tm))
    return Primer(name, tail + best.sequence, binding_sites=(best.site,))


def design_pair(
    template: SequenceRecord,
    start: int,
    end: int,
    *,
    forward_tail: str = "",
    reverse_tail: str = "",
    forward_name: str = "",
    reverse_name: str = "",
    target_tm: float = TARGET_TM,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
) -> tuple[Primer, Primer]:
    """Design a pair amplifying `start` to `end`, with their Tms as near each other as they go.

    `end` passes the length of a circular template when the amplicon crosses the origin. The
    two lengths are chosen together: the pair grading best, then the Tms sitting closest to
    each other and to `target_tm`.

    Raises
    ------
    ValueError
        If no annealing region fits the template at either end.
    """
    forwards = _annealing_options(template, start, Strand.FORWARD, polymerase, thresholds)
    reverses = _annealing_options(template, end, Strand.REVERSE, polymerase, thresholds)
    if not forwards or not reverses:
        raise ValueError(f"no annealing region fits this template at {start} or at {end}")
    forward, reverse = min(
        ((one, other) for one in forwards for other in reverses),
        key=lambda pair: _pair_score(pair[0], pair[1], target_tm, thresholds),
    )
    return (
        Primer(forward_name, forward_tail + forward.sequence, binding_sites=(forward.site,)),
        Primer(reverse_name, reverse_tail + reverse.sequence, binding_sites=(reverse.site,)),
    )


@dataclass(frozen=True, slots=True)
class _Option:
    site: BindingSite
    sequence: str
    tm: float
    grade: Status


def _annealing_options(
    template: SequenceRecord,
    position: int,
    strand: Strand,
    polymerase: Polymerase,
    thresholds: Thresholds,
) -> list[_Option]:
    length = len(template)
    circular = template.topology == "circular"
    if not circular and not 0 <= position <= length:
        return []
    shortest = int(_finite(thresholds.length.warn_low, thresholds.length.low))
    longest = min(int(_finite(thresholds.length.warn_high, thresholds.length.high)), length)
    options = []
    for size in range(shortest, longest + 1):
        start = position if strand is Strand.FORWARD else position - size
        if circular:
            start %= length
        elif not 0 <= start <= length - size:
            continue
        bases = template.extract(Segment(start, start + size))
        sequence = bases if strand is Strand.FORWARD else reverse_complement(bases)
        tm = melting_temperature(sequence, polymerase)
        options.append(
            _Option(
                BindingSite(start, start + size, strand),
                sequence,
                tm,
                worst((thresholds.length.grade(size), thresholds.tm.grade(tm))),
            )
        )
    return options


def _option_score(option: _Option, target_tm: float) -> tuple[int, float, int]:
    return STATUSES.index(option.grade), abs(option.tm - target_tm), len(option.sequence)


def _pair_score(
    forward: _Option, reverse: _Option, target_tm: float, thresholds: Thresholds
) -> tuple[int, float, int]:
    difference = abs(forward.tm - reverse.tm)
    grade = worst((forward.grade, reverse.grade, thresholds.tm_difference.grade(difference)))
    drift = max(abs(forward.tm - target_tm), abs(reverse.tm - target_tm))
    return STATUSES.index(grade), drift + difference, len(forward.sequence) + len(reverse.sequence)


def _finite(value: float, fallback: float) -> float:
    return value if math.isfinite(value) else fallback

"""A throwaway splitter for issue #261: cut a cargo into oligo-sized Golden Gate fragments.

Two searches over the same inputs, so they can be measured against each other.

`simple_split` is what `cloning.goldengate.design.design_overhangs` does today, lifted to
positions it chooses itself: cuts spaced evenly, then one greedy first-fit pass picking each
overhang by its own on-target count, with no backtracking and no set score in the loop.

`clever_split` is the shape `docs/research/dad-split.md` §3.4 recommends: an interval dynamic
program pins the fewest fragments and the exact window each cut may sit in, then depth-first
branch and bound over the per-cut candidate lists maximises Pryor 2020 fidelity.

The bound. Adding an overhang to a set can only raise the denominators of the overhangs
already in it, and every further factor is at most one, so a partial set's fidelity is an
upper bound on any set extending it. That makes the bound admissible, which is what lets the
search prune and still be exact over the candidate lists it was given.

`allow_edits` is GoldenHinges' third answer, built here because that package will not install:
a candidate list is widened by the 4-mers one synonymous codon change can reach.

Nothing here is imported by `liulab_mbio`; the package is read and never modified.
"""

from __future__ import annotations

import time
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from functools import cache

from liulab_mbio import overhangs
from liulab_mbio.codons import amino_acid
from liulab_mbio.ligase import LigaseProfile
from liulab_mbio.sequence import reverse_complement

OVERHANG = 4

# A candidate: where the cut sits, the overhang it spells, and the codon edit that spells it.
Candidate = tuple[int, str, tuple[int, str] | None]


@dataclass(frozen=True, slots=True)
class Budget:
    """What one oligo leaves for cargo, and the bounds a fragment is held to.

    Parameters
    ----------
    oligo
        Total synthesised length of one oligo.
    overhead
        Everything on the oligo that is not cargo: both primer sites, both Type IIS sites and
        their spacers. The README says what was assumed; §6.1 is not decided here.
    floor
        Shortest cargo stretch a fragment may carry.
    """

    oligo: int = 300
    overhead: int = 54
    floor: int = 40

    @property
    def ceiling(self) -> int:
        """Longest cargo stretch one oligo can carry."""
        return self.oligo - self.overhead


# A default every search shares, so a caller changes only the knob it names.
DEFAULT = Budget()


@dataclass(frozen=True, slots=True)
class Split:
    """One answer: where the cuts fall, what they spell, and what it cost to find.

    Parameters
    ----------
    cuts
        The 0-based top-strand index each downstream fragment begins at, ascending. The cargo's
        two ends are not among them.
    overhangs
        The 4-mer at each cut, in the same order.
    value
        Pryor 2020 fidelity of the internal set, 0.0 where no set was found.
    fragments
        How many pieces the cargo was cut into.
    seconds
        Wall time the search took.
    scored
        How many times `fidelity` was called.
    feasible
        Whether every cut found a legal overhang.
    edits
        Synonymous codon changes the chosen set spends.
    sequence
        The cargo as the answer leaves it, which differs from the input only where it edited.
    note
        What happened, where the numbers do not say it.
    """

    cuts: tuple[int, ...]
    overhangs: tuple[str, ...]
    value: float
    fragments: int
    seconds: float
    scored: int = 0
    feasible: bool = True
    edits: int = 0
    sequence: str = ""
    note: str = ""


@dataclass(slots=True)
class _Scorer:
    """Pryor 2020 fidelity over one enzyme and one ligase profile, with a call count."""

    enzyme: str
    profile: LigaseProfile | None
    calls: int = 0
    _seen: dict[tuple[str, ...], float] = field(default_factory=dict)

    def __call__(self, chosen: Sequence[str]) -> float:
        """Score a set, remembering what it already scored."""
        if not chosen:
            return 1.0
        key = tuple(chosen)
        if key in self._seen:
            return self._seen[key]
        self.calls += 1
        value = overhangs.fidelity(
            key, self.enzyme, profile=self.profile, prefer_profile=self.profile is not None
        ).value
        self._seen[key] = value
        return value


def _masks(length: int, budget: Budget, count: int) -> tuple[list[int], list[int]] | None:
    """Forward and backward reachability bitmasks over cut positions, or ``None`` where none.

    `reach[t]` holds the positions a prefix of `t` fragments can end at; `back[t]` the positions
    whose suffix `t` fragments can cover.
    """
    span = length - OVERHANG
    if span < budget.floor:
        return None
    steps = range(budget.floor, budget.ceiling + 1)
    mask = (1 << (span + 1)) - 1
    reach = [1] + [0] * count
    for depth in range(1, count + 1):
        total = 0
        for step in steps:
            total |= reach[depth - 1] << step
        reach[depth] = total & mask
    if not reach[count] >> span & 1:
        return None
    back = [1 << span] + [0] * count
    for depth in range(1, count + 1):
        total = 0
        for step in steps:
            total |= back[depth - 1] >> step
        back[depth] = total & mask
    return reach, back


def fragment_windows(length: int, budget: Budget, count: int) -> tuple[tuple[int, ...], ...] | None:
    """Every position each cut may take, for a cargo cut into `count` fragments.

    The interval dynamic program of `dad-split.md` §5, run forward and backward as bitmasks: a
    cut at depth `t` is legal when the prefix is coverable in `t` fragments and the suffix in
    `count - t`. The intersection is exact — no position is admitted that no whole partition
    contains, and none that one does is left out.

    A window is a marginal: it holds no position that no whole partition contains, but picking
    one position from each window independently can still break the chain, so the search below
    also holds each step to the length bound.

    Returns ``None`` where no partition into `count` fragments exists.
    """
    found = _masks(length, budget, count)
    if found is None:
        return None
    reach, back = found
    span = length - OVERHANG
    windows = []
    for depth in range(1, count):
        live = reach[depth] & back[count - depth]
        windows.append(tuple(at for at in range(span + 1) if live >> at & 1))
        if not windows[-1]:
            return None
    return tuple(windows)


def fewest_fragments(length: int, budget: Budget) -> int:
    """How few oligos this cargo needs, by the same dynamic program.

    Raises
    ------
    ValueError
        If no partition exists at any count.
    """
    for count in range(1, length // budget.floor + 2):
        if fragment_windows(length, budget, count) is not None:
            return count
    raise ValueError(f"no partition of {length} bases fits {budget.floor}..{budget.ceiling}")


@cache
def _synonyms(codon: str) -> tuple[str, ...]:
    """Every other codon spelling the same amino acid."""
    if len(codon) < 3 or set(codon) - set("ACGT"):
        return ()
    wanted = amino_acid(codon)
    return tuple(
        other
        for first in "ACGT"
        for second in "ACGT"
        for third in "ACGT"
        if (other := first + second + third) != codon and amino_acid(other) == wanted
    )


def _reachable(sequence: str, at: int, frame: int) -> list[tuple[str, tuple[int, str]]]:
    """Every 4-mer one synonymous codon change can spell at this cut, with the change."""
    made = []
    for base in range(at, at + OVERHANG):
        if base < frame:
            continue
        head = frame + (base - frame) // 3 * 3
        if head + 3 > len(sequence):
            continue
        codon = sequence[head : head + 3]
        for other in _synonyms(codon):
            edited = sequence[:head] + other + sequence[head + 3 :]
            made.append((edited[at : at + OVERHANG], (head, other)))
    return made


def _options(
    sequence: str,
    positions: Iterable[int],
    enzyme: str,
    reserved: Sequence[str],
    min_distance: int,
    frame: int | None,
) -> list[Candidate]:
    """One cut's candidates: every legal 4-mer the window spells, then what an edit could spell.

    Unedited candidates come first, so the search reaches a set that spends no edit before one
    that does.
    """
    made: list[Candidate] = []
    spots = list(positions)
    for at in spots:
        bases = sequence[at : at + OVERHANG]
        if len(bases) == OVERHANG and _allowed(bases, enzyme, reserved, min_distance):
            made.append((at, bases, None))
    if frame is not None:
        for at in spots:
            for bases, edit in _reachable(sequence, at, frame):
                if _allowed(bases, enzyme, reserved, min_distance):
                    made.append((at, bases, edit))
    return made


def _allowed(bases: str, enzyme: str, reserved: Sequence[str], min_distance: int) -> bool:
    """Whether a candidate passes every rule that looks at it alone."""
    return overhangs.refusal(bases, enzyme, reserved=reserved, min_distance=min_distance) is None


def _joins(
    candidate: str, taken: Sequence[str], enzyme: str, reserved: Sequence[str], min_distance: int
) -> bool:
    """Whether a candidate may join a partial set, by the pairwise rules alone."""
    return (
        overhangs.refusal(
            candidate, enzyme, taken=taken, reserved=reserved, min_distance=min_distance
        )
        is None
    )


def clever_split(
    sequence: str,
    *,
    budget: Budget = DEFAULT,
    enzyme: str = "BsaI",
    profile: LigaseProfile | None = None,
    reserved: Sequence[str] = ("AGGA", "TTCC"),
    min_distance: int = overhangs.MIN_DISTANCE,
    count: int | None = None,
    cap: int | None = None,
    frame: int | None = None,
    seconds_limit: float = 10.0,
) -> Split:
    """Pin the fewest cuts by the dynamic program, then branch and bound over their overhangs.

    `cap` keeps at most that many candidates per cut, nearest an evenly spaced anchor first,
    and is the one knob that trades runtime for optimality. Nearest-first matters: thinning a
    window evenly instead leaves picks too far apart to satisfy the length chain, and the search
    then reports an infeasibility that is the thinning's and not the cargo's. `frame` allows synonymous edits, widening every
    candidate list. `seconds_limit` stops a search that will not finish and says so in the
    note, so a timeout is reported and never hidden.
    """
    start = time.perf_counter()
    fragments = count or fewest_fragments(len(sequence), budget)
    found = _masks(len(sequence), budget, fragments)
    windows = fragment_windows(len(sequence), budget, fragments)
    if found is None or windows is None:
        return _nothing(fragments, start, f"no partition into {fragments} fragments")
    _, back = found
    span = len(sequence) - OVERHANG
    options = [
        _options(sequence, window, enzyme, reserved, min_distance, frame) for window in windows
    ]
    for at, legal in enumerate(options):
        if not legal:
            return _nothing(
                fragments, start, f"cut {at + 1} spells no legal overhang anywhere in its window"
            )
        anchor = round(span * (at + 1) / fragments)
        ranked = sorted(legal, key=lambda one: (one[2] is not None, abs(one[0] - anchor)))
        options[at] = ranked[:cap] if cap is not None else ranked
    scorer = _Scorer(enzyme, profile)
    best: list[Candidate] | None = None
    highest = 0.0
    nodes = 0
    stopped = False

    def walk(depth: int, previous: int, taken: list[Candidate]) -> None:
        nonlocal best, highest, nodes, stopped
        nodes += 1
        if nodes % 256 == 0 and time.perf_counter() - start > seconds_limit:
            stopped = True
            return
        chosen = [bases for _, bases, _ in taken]
        bound = scorer(sorted(chosen))
        if best is not None and bound <= highest:
            return
        if depth == len(options):
            best, highest = list(taken), bound
            return
        for candidate in options[depth]:
            at = candidate[0]
            if not budget.floor <= at - previous <= budget.ceiling:
                continue
            if not back[fragments - depth - 1] >> at & 1:
                continue
            if _joins(candidate[1], chosen, enzyme, reserved, min_distance):
                taken.append(candidate)
                walk(depth + 1, at, taken)
                taken.pop()
                if stopped:
                    return

    walk(0, 0, [])
    seconds = time.perf_counter() - start
    note = f"stopped at {nodes} nodes: best seen, not proven" if stopped else ""
    if best is None:
        return _nothing(
            fragments,
            start,
            note or "no set of overhangs satisfies every rule and the length bound at once",
            scorer.calls,
        )
    settled = list(best)
    built = sequence
    for _, _, edit in settled:
        if edit is not None:
            head, codon = edit
            built = built[:head] + codon + built[head + 3 :]
    return Split(
        tuple(at for at, _, _ in settled),
        tuple(bases for _, bases, _ in settled),
        highest,
        fragments,
        seconds,
        scorer.calls,
        edits=sum(edit is not None for _, _, edit in settled),
        sequence=built,
        note=note,
    )


def _nothing(fragments: int, start: float, note: str, scored: int = 0) -> Split:
    """Say there is no answer, which is a real outcome and names why."""
    return Split(
        (), (), 0.0, fragments, time.perf_counter() - start, scored, feasible=False, note=note
    )


def simple_split(
    sequence: str,
    *,
    budget: Budget = DEFAULT,
    enzyme: str = "BsaI",
    profile: LigaseProfile | None = None,
    reserved: Sequence[str] = ("AGGA", "TTCC"),
    min_distance: int = overhangs.MIN_DISTANCE,
    count: int | None = None,
    window: int = 16,
) -> Split:
    """Space the cuts evenly, then take each overhang first-fit, as `design_overhangs` does.

    Candidates are ranked by their own on-target count, the per-overhang proxy the shipped
    greedy uses. The set score is computed once at the end, to report and not to choose.
    """
    start = time.perf_counter()
    fragments = count or fewest_fragments(len(sequence), budget)
    span = len(sequence) - OVERHANG
    anchors = [round(span * at / fragments) for at in range(1, fragments)]
    table, _ = overhangs.scoring(enzyme, profile=profile, prefer_profile=profile is not None)
    taken: list[str] = []
    places: list[int] = []
    for anchor in anchors:
        previous = places[-1] if places else 0
        offsets = [0, *(sign * step for step in range(1, window + 1) for sign in (-1, 1))]
        found = [
            (at, sequence[at : at + OVERHANG])
            for offset in offsets
            if 0 <= (at := anchor + offset) <= span
            and budget.floor <= at - previous <= budget.ceiling
        ]
        ranked = sorted(
            found,
            key=lambda pair: -(table.count(pair[1], reverse_complement(pair[1])) if table else 0),
        )
        for at, bases in ranked:
            if _joins(bases, taken, enzyme, reserved, min_distance):
                taken.append(bases)
                places.append(at)
                break
        else:
            return _nothing(
                fragments, start, f"no overhang within {window} bases of {anchor} joins the set"
            )
    value = overhangs.fidelity(
        taken, enzyme, profile=profile, prefer_profile=profile is not None
    ).value
    return Split(
        tuple(places),
        tuple(taken),
        value,
        fragments,
        time.perf_counter() - start,
        1,
        sequence=sequence,
    )


def assembles(sequence: str, split: Split) -> bool:
    """Reassemble the fragments overhang by overhang and check the cargo comes back.

    OMEGA's round trip, kept: a design that does not rebuild its own input is not a design.
    """
    if not split.feasible:
        return False
    built = split.sequence or sequence
    edges = (0, *split.cuts, len(built) - OVERHANG)
    pieces = [built[edges[at] : edges[at + 1] + OVERHANG] for at in range(len(edges) - 1)]
    for piece, overhang in zip(pieces, split.overhangs, strict=False):
        if piece[-OVERHANG:] != overhang:
            return False
    whole = pieces[0]
    for piece in pieces[1:]:
        if whole[-OVERHANG:] != piece[:OVERHANG]:
            return False
        whole += piece[OVERHANG:]
    return whole == built


def cargo_lengths(sequence: str, split: Split) -> tuple[int, ...]:
    """How many bases of cargo each oligo carries, overhang included."""
    if not split.feasible:
        return ()
    edges = (0, *split.cuts, len(sequence) - OVERHANG)
    return tuple(edges[at + 1] + OVERHANG - edges[at] for at in range(len(edges) - 1))

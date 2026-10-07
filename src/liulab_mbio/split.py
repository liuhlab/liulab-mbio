"""Split a cargo too long to synthesise into fragments that assemble back into it.

A synthesised oligo holds only so many bases, so a gene longer than one oligo is ordered as
several and joined in one pot. Two questions follow, and this answers both:

- **Where the cuts fall.** An interval dynamic program pins the fewest fragments the budget
  allows and, for each junction, every position it could still sit at with a legal completion
  behind and ahead of it. That window is exact: a position outside it spells no design.
- **What each cut leaves.** One pass over the junctions takes, at each, the overhang scoring
  the whole set best by `liulab_mbio.overhangs.fidelity`. The candidates are the positions the
  window allows, so the search never leaves the set of designs the first question proved.

Fragment count is the axis that decides whether an assembly works, so the dynamic program
minimises it and nothing trades it away. Fidelity ranks only the choices left at that count.

Every span here is 0-based and half-open, and a fragment's span holds the overhang at each of
its ends: an internal overhang belongs to the two fragments it joins, so the fragments overlap
by that much and `sum(len(fragment)) == len(cargo) + overhang * (pieces - 1)`.
`docs/adr/0001-coordinates.md` has the rule.
"""

from collections.abc import Iterable, Sequence
from dataclasses import KW_ONLY, dataclass
from math import ceil

from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.ligase import LigaseProfile
from liulab_mbio.overhangs import (
    MIN_DISTANCE,
    Choice,
    FidelityReport,
    Junction,
    Rejection,
    best_overhang,
    fidelity,
)
from liulab_mbio.sequence import Segment, SequenceRecord, Strand
from liulab_mbio.sites import EnzymeLike, find_sites


def shortest_fragment(overhang: int) -> int:
    """Return the shortest fragment an enzyme leaving `overhang` bases cuts out of an oligo.

    A fragment carries an overhang at each end, so a shorter one is overhang throughout and has
    no double-stranded core between them. Nothing published sets a floor on a Golden Gate
    fragment, so this is the package's own arithmetic and `Budget.minimum` is where a caller
    sets a longer one.

    Examples
    --------
    >>> shortest_fragment(4)
    9
    """
    return 2 * overhang + 1


@dataclass(frozen=True, slots=True)
class Budget:
    """How much of a cargo one oligo may carry.

    Parameters
    ----------
    oligo
        The whole oligo, every base the vendor synthesises.
    overhead
        What that oligo spends on anything but the cargo: its primer sites, the recognition
        sites that cut the fragment out, and their spacers.
    minimum
        The shortest fragment to design, or ``None`` to take the enzyme's own floor.
    """

    oligo: int
    overhead: int
    minimum: int | None = None

    @property
    def span(self) -> int:
        """The templated span: the most of the cargo one oligo carries, overhangs included."""
        return self.oligo - self.overhead

    def floor(self, overhang: int) -> int:
        """Return the shortest fragment to design, for an enzyme leaving `overhang` bases.

        Raises
        ------
        ValueError
            If the overhead leaves no room for a fragment that long.

        Examples
        --------
        >>> Budget(350, 74).floor(4)
        9
        """
        minimum = self.minimum if self.minimum is not None else shortest_fragment(overhang)
        if self.span < minimum:
            raise ValueError(
                f"an oligo of {self.oligo} nt spending {self.overhead} nt on overhead carries "
                f"{self.span} nt of cargo, short of the {minimum} nt minimum fragment"
            )
        return minimum


@dataclass(frozen=True, slots=True)
class Fragment:
    """One piece of the cargo, as one oligo templates it.

    Parameters
    ----------
    index
        Which piece this is, counting from zero in cargo order.
    span
        Where it lies in the cargo, its overhang at each end included.
    overhangs
        What its two cuts leave, 5' then 3', as the cargo's top strand spells them.
    """

    index: int
    span: Segment
    overhangs: tuple[str, str]

    @property
    def length(self) -> int:
        """How many bases the oligo templates for this fragment."""
        return self.span.end - self.span.start


@dataclass(frozen=True, slots=True)
class CargoSplit:
    """A cargo and the fragments that assemble back into it.

    Parameters
    ----------
    cargo
        What was split.
    enzyme
        The Type IIS enzyme cutting each fragment out of its oligo.
    budget
        What one oligo could carry.
    fragments
        In cargo order.
    junctions
        One per internal cut, each carrying the window the dynamic program proved.
    choices
        What each junction took, and what was refused on the way.
    fidelity
        How well the whole set should ligate.
    """

    cargo: SequenceRecord
    _: KW_ONLY
    enzyme: Enzyme
    budget: Budget
    fragments: tuple[Fragment, ...]
    junctions: tuple[Junction, ...]
    choices: tuple[Choice, ...]
    fidelity: FidelityReport

    @property
    def pieces(self) -> int:
        """How many oligos this cargo takes."""
        return len(self.fragments)

    @property
    def overhangs(self) -> tuple[str, ...]:
        """The internal overhangs, in cargo order. The terminal two are the cargo's own ends."""
        return tuple(choice.overhang for choice in self.choices)

    def sequence(self, fragment: Fragment) -> str:
        """Return the bases one fragment templates."""
        return self.cargo.extract(fragment.span)

    def reassembled(self) -> str:
        """Join the fragments back on their shared overhangs, which must spell the cargo."""
        bases = self.sequence(self.fragments[0])
        for fragment in self.fragments[1:]:
            overlap = len(fragment.overhangs[0])
            bases += self.sequence(fragment)[overlap:]
        return bases


def fewest_pieces(length: int, budget: Budget) -> int:
    """How few oligos a cargo of `length` could take, counting the overhang each cut shares.

    This is the arithmetic floor and not yet a design: a junction also has to spell a legal
    overhang, which `split_cargo` is what decides.

    Examples
    --------
    >>> fewest_pieces(1149, Budget(350, 74))
    5
    >>> fewest_pieces(276, Budget(350, 74))
    1
    """
    overhang = 4
    span = budget.span
    if length <= span:
        return 1
    return ceil((length - overhang) / (span - overhang))


def split_cargo(
    cargo: SequenceRecord,
    enzyme: EnzymeLike,
    *,
    budget: Budget,
    reserved: Iterable[str] = (),
    avoid: Iterable[EnzymeLike] = (),
    min_distance: int = MIN_DISTANCE,
    allow_uniform: bool = False,
    profile: LigaseProfile | None = None,
    prefer_profile: bool = False,
) -> CargoSplit:
    """Split `cargo` into the fewest oligo-sized fragments that assemble back into it.

    The fragments overlap by one overhang at each internal junction, and the two outermost
    overhangs are the cargo's own first and last bases, which are what it joins the rest of the
    construct on. `reserved` holds overhangs out by name, so an internal junction never spells
    what a step outside this assembly already spends.

    Parameters
    ----------
    cargo
        The bases to split. Linear: a span across an origin has no fragment order.
    enzyme
        The Type IIS enzyme cutting each fragment out of its oligo, which sets overhang length.
    budget
        What one oligo carries.
    reserved
        Overhangs a step outside this assembly spends, which no internal junction may take.
    avoid
        Enzymes besides this one whose sites the oligo must not spell.
    min_distance
        How many bases two overhangs of the set must differ by.
    allow_uniform
        Accept an overhang of one base kind.
    profile, prefer_profile
        A ligase's own matrix, and whether to score with it where shipped data covers the
        enzyme. `liulab_mbio.ligase.read_profile` reads one.

    Raises
    ------
    ValueError
        If the cargo is circular, carries a site of `enzyme` itself, is shorter than one
        fragment, or spells no legal overhang set at any fragment count the budget allows --
        which says how far the budget does reach.

    Examples
    --------
    >>> from liulab_mbio.sequence import SequenceRecord
    >>> split = split_cargo(SequenceRecord("ACGT" * 60), "BsaI", budget=Budget(200, 74))
    >>> split.pieces, split.reassembled() == "ACGT" * 60
    (2, True)
    """
    one = get_enzyme(enzyme) if isinstance(enzyme, str) else enzyme
    if cargo.topology == "circular":
        raise ValueError("a circular record has no first fragment, so it cannot be split")
    _check_cutter(cargo, one)
    bases = str(cargo.sequence)
    length = len(bases)
    overhang = one.overhang_length
    minimum = budget.floor(overhang)
    if length < shortest_fragment(overhang):
        raise ValueError(
            f"a cargo of {length} bases is shorter than the {shortest_fragment(overhang)} bases "
            f"a {one.name} fragment needs between its two overhangs"
        )
    if length <= budget.span:
        return _whole(cargo, one, budget, bases, profile, prefer_profile)
    held = tuple(str(one).upper() for one in reserved)
    refused: list[str] = []
    most = min(length // minimum, _supply(overhang, min_distance))
    for pieces in range(fewest_pieces(length, budget), most + 1):
        windows = _windows(length, pieces, budget, minimum, overhang)
        if windows is None:
            continue
        found = _choose(
            cargo,
            bases,
            one,
            windows,
            pieces,
            budget,
            minimum,
            held,
            avoid,
            min_distance,
            allow_uniform,
            profile,
            prefer_profile,
        )
        if isinstance(found, str):
            refused.append(f"at {pieces} fragments, {found}")
            continue
        return _split(cargo, one, budget, bases, found, length, overhang, profile, prefer_profile)
    raise ValueError(_ceiling(length, budget, overhang, min_distance, refused))


def _check_cutter(cargo: SequenceRecord, enzyme: Enzyme) -> None:
    """Refuse a cargo spelling the site of the enzyme that cuts its fragments out.

    The cuts this design places are the two the oligo carries at a fragment's ends, so a site
    the cargo spells itself is a third one, inside a fragment. No split clears it: it is the
    cargo's bases, and domestication is what takes it out.

    Raises
    ------
    ValueError
        Naming how many sites, and where the first one reads.
    """
    carried = find_sites(cargo, enzyme)
    if not carried:
        return
    first = carried[0]
    strand = "forward" if first.strand == Strand.FORWARD else "reverse"
    raise ValueError(
        f"the cargo spells {len(carried)} {enzyme.name} site(s), the first at {first.start} on "
        f"the {strand} strand, and {enzyme.name} is what cuts each fragment out of its oligo, "
        "so that site would cut a fragment apart"
    )


def _supply(overhang: int, min_distance: int) -> int:
    """Return the most fragments the overhang supply holds, whatever the cargo.

    The internal overhangs of a split into `pieces` fragments, their reverse complements with
    them, are ``2 * (pieces - 1)`` words of a distance-`min_distance` code of length `overhang`
    over four bases. The Singleton bound holds such a code to
    ``4 ** (overhang - min_distance + 1)`` words.
    """
    return 4 ** max(overhang - min_distance + 1, 0) // 2 + 1


def _ceiling(
    length: int, budget: Budget, overhang: int, min_distance: int, refused: Sequence[str]
) -> str:
    """Say how far this budget reaches and why this cargo is past it."""
    reason = "; ".join(refused) if refused else "no fragment count fits the length budget"
    return (
        f"a cargo of {length} bases needs at least {fewest_pieces(length, budget)} fragments at "
        f"{budget.span} bases an oligo, and no legal overhang set exists: {reason}. The supply of "
        f"{overhang}-base overhangs standing {min_distance} bases from one another and from every "
        f"reverse complement runs out at {_supply(overhang, min_distance)} fragments, so a longer "
        "oligo is what lifts this and a longer search is not"
    )


def _whole(
    cargo: SequenceRecord,
    enzyme: Enzyme,
    budget: Budget,
    bases: str,
    profile: LigaseProfile | None,
    prefer_profile: bool,
) -> CargoSplit:
    """Return the split of a cargo that already fits one oligo: one fragment, no junction."""
    overhang = enzyme.overhang_length
    ends = (bases[:overhang], bases[-overhang:])
    return CargoSplit(
        cargo,
        enzyme=enzyme,
        budget=budget,
        fragments=(Fragment(0, Segment(0, len(bases)), ends),),
        junctions=(),
        choices=(),
        fidelity=fidelity((), enzyme, profile=profile, prefer_profile=prefer_profile),
    )


def _spread(mask: int, width: int) -> int:
    """Return the OR of `mask` shifted left by every distance below `width`."""
    out, done = mask, 1
    while done < width:
        step = min(done, width - done)
        out |= out << step
        done += step
    return out


def _windows(
    length: int, pieces: int, budget: Budget, minimum: int, overhang: int
) -> tuple[int, ...] | None:
    """Every position each internal cut could take, as a bitmask, or ``None`` for none at all.

    A cut at `p` opens an overhang at ``p`` to ``p + overhang``, so the fragment before it runs
    from the cut before to ``p + overhang`` and the last fragment runs from the last cut to the
    end of the cargo. Both are held to the budget, which is what keeps the final fragment inside
    the oligo as well as every other.
    """
    low = minimum - overhang
    high = budget.span - overhang
    if low < 1 or high < low:
        return None
    whole = (1 << (length + 1)) - 1
    forward = [1]
    for _ in range(pieces - 1):
        forward.append((_spread(forward[-1], high - low + 1) << low) & whole)
    last = ((1 << (length - minimum + 1)) - 1) ^ ((1 << max(length - budget.span, 0)) - 1)
    backward = [last]
    for _ in range(pieces - 2):
        backward.append(_spread(backward[-1] >> high, high - low + 1) & whole)
    backward.reverse()
    windows = tuple(forward[index] & backward[index - 1] for index in range(1, pieces))
    return None if any(window == 0 for window in windows) else windows


def _choose(
    cargo: SequenceRecord,
    bases: str,
    enzyme: Enzyme,
    windows: Sequence[int],
    pieces: int,
    budget: Budget,
    minimum: int,
    reserved: Sequence[str],
    avoid: Iterable[EnzymeLike],
    min_distance: int,
    allow_uniform: bool,
    profile: LigaseProfile | None,
    prefer_profile: bool,
) -> tuple[tuple[Junction, Choice], ...] | str:
    """Take one overhang per junction, the whole set's best first, or say what refused one.

    A junction's candidates are the positions one legal fragment away from the cut before it
    that can still reach the end of the cargo, so every choice leaves a legal design behind it
    and the last fragment is inside the oligo as surely as the first.
    """
    overhang = enzyme.overhang_length
    low = minimum - overhang
    high = budget.span - overhang
    others = tuple(avoid)
    taken: list[str] = []
    settled: list[tuple[Junction, Choice]] = []
    previous = 0
    for index, window in enumerate(windows, 1):
        positions = _positions(window & _step(previous, low, high, len(bases)))
        if not positions:
            return f"junction {index}/{pieces} has no position a legal fragment away"
        anchor = round(index * (len(bases) - overhang) / pieces)
        junction = Junction(
            f"junction {index}/{pieces}",
            record=cargo,
            position=anchor,
            scarless=True,
            window=max(abs(one - anchor) for one in positions),
        )
        choice, rejected = best_overhang(
            junction,
            tuple((one - anchor, bases[one : one + overhang]) for one in positions),
            enzyme,
            taken=taken,
            reserved=reserved,
            avoid=others,
            min_distance=min_distance,
            allow_uniform=allow_uniform,
            profile=profile,
            prefer_profile=prefer_profile,
        )
        if choice is None:
            return _stuck(index, pieces, rejected)
        settled.append((junction, choice))
        taken.append(choice.overhang)
        previous = anchor + choice.offset
    return tuple(settled)


def _stuck(index: int, pieces: int, rejected: Sequence[Rejection]) -> str:
    """Say which rules refused every candidate of one junction."""
    rules = ", ".join(sorted({one.rule for one in rejected})) or "no candidate at all"
    return f"junction {index}/{pieces} has no overhang left ({rules})"


def _step(previous: int, low: int, high: int, length: int) -> int:
    """Every position one legal fragment past `previous`, as a bitmask."""
    start = previous + low
    stop = min(previous + high, length)
    if stop < start:
        return 0
    return ((1 << (stop + 1)) - 1) ^ ((1 << start) - 1)


def _positions(mask: int) -> tuple[int, ...]:
    """Return the positions a bitmask holds, lowest first."""
    found: list[int] = []
    while mask:
        low = mask & -mask
        found.append(low.bit_length() - 1)
        mask ^= low
    return tuple(found)


def _split(
    cargo: SequenceRecord,
    enzyme: Enzyme,
    budget: Budget,
    bases: str,
    settled: Sequence[tuple[Junction, Choice]],
    length: int,
    overhang: int,
    profile: LigaseProfile | None,
    prefer_profile: bool,
) -> CargoSplit:
    """Build the split from the junctions that were settled."""
    cuts = [junction.position + choice.offset for junction, choice in settled]
    edges = [0, *cuts, length]
    fragments: list[Fragment] = []
    for index in range(len(edges) - 1):
        start = edges[index]
        end = edges[index + 1] + (overhang if index < len(cuts) else 0)
        head = bases[start : start + overhang]
        tail = bases[end - overhang : end]
        fragments.append(Fragment(index, Segment(start, end), (head, tail)))
    choices = tuple(choice for _, choice in settled)
    return CargoSplit(
        cargo,
        enzyme=enzyme,
        budget=budget,
        fragments=tuple(fragments),
        junctions=tuple(junction for junction, _ in settled),
        choices=choices,
        fidelity=fidelity(
            [choice.overhang for choice in choices],
            enzyme,
            profile=profile,
            prefer_profile=prefer_profile,
        ),
    )

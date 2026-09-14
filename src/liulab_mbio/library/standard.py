"""Choosing the overhang standard a library build works to, and what it costs the proteins.

One standard serves the whole library, so a position's entry overhang forces terminal codons on
every member of the part lists either side of it. The standard is scored across every member
rather than taken from the first set that passes, and it reports each part's wild-type residues
beside the synthesised ones, which is the bookkeeping a change table is written from.

A junction spells whole codons. The coding sequence upstream donates what the overhang needs to
finish the first of them, and the overhang's remaining bases spell codons of their own, which the
part either side may own: that split is chosen with the overhang, for the fewest changes. The
first position has the vector upstream rather than a part list, and a vector's stuffer is padding,
so it absorbs those codons and charges no protein for them.

The overhang rules are `liulab_mbio.goldengate.design`'s, reached through its `refusal`: length,
palindrome, one base class, repeat, near-duplicate, and a tail spelling no further site. That is
the one edge `library/` has to that pipeline. An overhang spelling a stop where a coding sequence
reads through it is no candidate at all, whichever part would own the codon.
"""

from collections.abc import Mapping, Sequence
from dataclasses import KW_ONLY, dataclass
from functools import cache
from itertools import product
from typing import Literal

from liulab_mbio.codons import CodonUsage, amino_acid, codon_usage
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.goldengate.design import (
    MIN_DISTANCE,
    Choice,
    Junction,
    Rejection,
    ligation_matrix,
    refusal,
)
from liulab_mbio.library.scheme import Scheme
from liulab_mbio.sequence import reverse_complement

#: One part list: the name each member is ordered under, and the protein it codes for.
type PartList = Mapping[str, str]

#: Which end of a part a junction spells.
type End = Literal["5'", "3'"]

#: What a report calls the junction every part's 3' end leaves.
SCAR = "cloning scar"


def junction_residues(overhang_length: int) -> tuple[int, int]:
    """Bases the upstream coding sequence donates to a junction, and codons the overhang spells.

    A junction spells whole codons: the sequence upstream donates what the overhang needs to
    finish the first, and the rest of the overhang spells codons of its own.

    Raises
    ------
    ValueError
        If `overhang_length` is not positive.

    Examples
    --------
    >>> junction_residues(4)
    (2, 1)
    >>> junction_residues(3)
    (0, 1)
    """
    if overhang_length < 1:
        raise ValueError(f"an overhang has at least one base, got {overhang_length}")
    donated = (-overhang_length) % 3
    return donated, (overhang_length - (3 - donated) % 3) // 3


@dataclass(frozen=True, slots=True)
class Terminus:
    """What one junction spells at one end of one part, wild type beside synthesised.

    Parameters
    ----------
    part
        The member of the part list, as it is named and ordered.
    position
        The position that part fills.
    end
        Which end of the part the junction spells.
    wild_type
        The residues the protein has there, in protein order.
    synthesised
        The residues the standard spells there, one for one against `wild_type`.
    """

    part: str
    position: str
    end: End
    wild_type: str
    synthesised: str

    @property
    def changed_residues(self) -> int:
        """How many of these residues the standard changes."""
        return sum(
            one != other for one, other in zip(self.wild_type, self.synthesised, strict=True)
        )

    @property
    def changed(self) -> bool:
        """Whether the standard changes this part at this end at all."""
        return bool(self.changed_residues)


@dataclass(frozen=True, slots=True)
class Standard:
    """The overhang standard one library build works to, and what it costs the proteins.

    Parameters
    ----------
    entry_overhangs
        One per position, in the order the rounds fill them.
    scar_overhang
        The overhang every part's 3' end leaves, shared by every position.
    choices
        What each junction took and what was refused before it, the positions in order and the
        cloning scar last.
    termini
        Every part's wild-type and synthesised residues at a junction, so a change table is
        written from this alone.
    forced
        The junctions whose cheapest candidate still changes a residue, reported rather than
        changed in silence. For a junction left free that is every overhang the rules allow.
    """

    entry_overhangs: tuple[str, ...]
    _: KW_ONLY
    scar_overhang: str
    choices: tuple[Choice, ...] = ()
    termini: tuple[Terminus, ...] = ()
    forced: tuple[str, ...] = ()

    @property
    def cost(self) -> int:
        """How many residues the standard changes, over every member of every part list."""
        return sum(one.changed_residues for one in self.termini)

    @property
    def changes(self) -> tuple[Terminus, ...]:
        """Only the parts the standard changes, in the order the termini are reported."""
        return tuple(one for one in self.termini if one.changed)


def design_standard(
    scheme: Scheme,
    part_lists: Sequence[PartList],
    *,
    pinned: Mapping[str, str] | None = None,
    usage: CodonUsage | None = None,
    min_distance: int = MIN_DISTANCE,
    allow_uniform: bool = False,
) -> Standard:
    """Choose the overhang standard costing these part lists the fewest amino acids.

    Every junction is scored across every member of the part lists either side of it, and the
    set with the lowest total is the one chosen, so the standard is the cheapest for the actual
    sequences rather than the first that passes. The cloning scar is the scheme's own and is held
    to the same rules as the overhangs designed around it.

    Parameters
    ----------
    scheme
        The architecture the standard is written into: its enzymes, its positions and the
        overhang length they leave. Its stuffers are not read for overhangs, which is what this
        chooses.
    part_lists
        One per position, in the scheme's order: each member's name and the protein it codes for,
        as one-letter amino acids.
    pinned
        Entry overhangs fixed by hand, keyed by position name. A pinned overhang is held to every
        rule, and a refusal names the rule that refused it.
    usage
        The host's codon usage, which settles which residue stands in where the wild-type one
        cannot be spelled. Defaults to the shipped table.
    min_distance
        How many bases two overhangs of the standard must differ by.
    allow_uniform
        Accept an overhang of one base kind, as `design_overhangs` does.

    Raises
    ------
    ValueError
        If the part lists do not match the scheme's positions, a protein is too short or holds a
        letter that is not an amino acid, a pinned position is not one of the scheme's, or a
        junction has no overhang left — which names the rules and the last candidate refused.
    KeyError
        If the scheme names an enzyme this package does not ship.
    """
    table = usage if usage is not None else codon_usage()
    enzyme = scheme.internal
    avoid = (scheme.external, *scheme.blunt)
    codons = junction_residues(enzyme.overhang_length)[1]
    sites = _sites(scheme, part_lists, pinned or {}, codons)
    pool = _pool(enzyme, avoid, min_distance, allow_uniform)
    readings = [
        _readings(site, pool, table, enzyme, avoid, min_distance, allow_uniform) for site in sites
    ]
    order = sorted(range(len(sites)), key=lambda index: sites[index].pinned is None)
    chosen = _settle(readings, order, sites, enzyme, avoid, min_distance, allow_uniform)
    return _standard(sites, readings, order, chosen, enzyme, avoid, min_distance, allow_uniform)


@dataclass(frozen=True, slots=True)
class _Site:
    """One junction a standard chooses an overhang for, and the part lists either side of it."""

    name: str
    _: KW_ONLY
    upstream: PartList | None
    upstream_position: str
    downstream: PartList | None
    downstream_position: str
    pinned: str | None


def _sites(
    scheme: Scheme, part_lists: Sequence[PartList], pinned: Mapping[str, str], codons: int
) -> tuple[_Site, ...]:
    """Build one junction per position, plus the cloning scar, refusing part lists that do not fit.

    Raises
    ------
    ValueError
        On a part list count, an empty list, an unreadable protein, or a pinned name that is not
        a position of this scheme.
    """
    names = [position.name for position in scheme.positions]
    if len(part_lists) != len(names):
        raise ValueError(
            f"this scheme has {len(names)} position(s) — {', '.join(names)} — and "
            f"{len(part_lists)} part list(s) were given"
        )
    if strange := sorted(set(pinned) - set(names)):
        raise ValueError(f"pinned name(s) {', '.join(strange)} are not positions of this scheme")
    for name, parts in zip(names, part_lists, strict=True):
        _checked(name, parts, codons)
    made = [
        _Site(
            name,
            upstream=part_lists[index - 1] if index else None,
            upstream_position=names[index - 1] if index else "",
            downstream=part_lists[index],
            downstream_position=name,
            pinned=pinned.get(name),
        )
        for index, name in enumerate(names)
    ]
    made.append(
        _Site(
            SCAR,
            upstream=None,
            upstream_position="",
            downstream=None,
            downstream_position="",
            pinned=scheme.cloning_scar,
        )
    )
    return tuple(made)


def _checked(position: str, parts: PartList, codons: int) -> None:
    """Refuse a part list a junction could not be spelled across.

    Raises
    ------
    ValueError
        If the list is empty, a protein holds a letter that is not an amino acid, or a protein is
        shorter than the residues a junction spells.
    """
    if not parts:
        raise ValueError(f"the part list for position {position!r} holds no part")
    shortest = codons + 1
    residues = set(_codons_for()) - {"*"}
    for name, protein in parts.items():
        upper = protein.upper()
        if bad := sorted(set(upper) - residues):
            raise ValueError(f"part {name!r} is not one-letter amino acids: {''.join(bad)}")
        if len(upper) < shortest:
            raise ValueError(
                f"part {name!r} is {len(upper)} amino acid(s), and a junction of this scheme "
                f"spells {shortest} at one end"
            )


def _pool(
    enzyme: Enzyme, avoid: Sequence[Enzyme], min_distance: int, allow_uniform: bool
) -> tuple[str, ...]:
    """Every overhang the rules allow on its own, the best ligating first.

    Ranked by the enzyme's own measured ligation where data covers it, so a tie on amino acids
    is broken the way `design_overhangs` breaks one rather than by a rule of this module's own.
    """
    matrix = ligation_matrix(enzyme)
    every = ("".join(bases) for bases in product("ACGT", repeat=enzyme.overhang_length))
    ranked = sorted(
        every,
        key=lambda one: (-(matrix.count(one, reverse_complement(one)) if matrix else 0), one),
    )
    return tuple(
        one
        for one in ranked
        if refusal(one, enzyme, avoid=avoid, min_distance=min_distance, allow_uniform=allow_uniform)
        is None
    )


@dataclass(frozen=True, slots=True)
class _Reading:
    """One candidate overhang at one junction: what it spells, and what that costs."""

    overhang: str
    termini: tuple[Terminus, ...]
    cost: int


def _readings(
    site: _Site,
    pool: Sequence[str],
    usage: CodonUsage,
    enzyme: Enzyme,
    avoid: Sequence[Enzyme],
    min_distance: int,
    allow_uniform: bool,
) -> tuple[_Reading, ...]:
    """Every overhang this junction could take, cheapest first.

    A pinned overhang is the only candidate, and is held to the rules here so that a refusal
    names the rule rather than reporting an empty search.

    Raises
    ------
    ValueError
        If a pinned overhang breaks a rule or spells a stop, or nothing at all can be spelled.
    """
    if site.pinned is not None:
        broken = refusal(
            site.pinned,
            enzyme,
            avoid=avoid,
            min_distance=min_distance,
            allow_uniform=allow_uniform,
        )
        if broken is not None:
            raise ValueError(
                f"the overhang {site.pinned.upper()!r} pinned at {site.name!r} is refused: "
                f"{broken.rule}, because {broken.detail}"
            )
        candidates: Sequence[str] = (site.pinned.upper(),)
    else:
        candidates = pool
    made = [
        reading
        for reading in (_reading(site, one, usage) for one in candidates)
        if reading is not None
    ]
    if not made:
        raise ValueError(
            f"junction {site.name!r} has no overhang left: every candidate spells a stop where "
            "the product reads through it"
        )
    return tuple(sorted(made, key=lambda one: one.cost))


def _reading(site: _Site, overhang: str, usage: CodonUsage) -> _Reading | None:
    """Return the cheapest way to read this overhang here, or ``None`` where only a stop fits.

    The whole codons an overhang spells may be owned by the part either side of it, and that
    split is what is chosen: the same DNA, charged where it costs least.
    """
    if site.upstream is None and site.downstream is None:
        return _Reading(overhang, (), 0)
    suffix, codons = _spelling(overhang)
    if any(amino_acid(one) == "*" for one in codons):
        return None
    best: _Reading | None = None
    for split in range(len(codons) + 1):
        before = _termini(
            site.upstream, site.upstream_position, "3'", suffix, codons[:split], usage
        )
        after = _termini(site.downstream, site.downstream_position, "5'", "", codons[split:], usage)
        if before is None or after is None:
            continue
        termini = before + after
        cost = sum(one.changed_residues for one in termini)
        if best is None or cost < best.cost:
            best = _Reading(overhang, termini, cost)
    return best


def _spelling(overhang: str) -> tuple[str, tuple[str, ...]]:
    """Return the bases finishing the upstream part's last codon, and the codons after them."""
    donated, _ = junction_residues(len(overhang))
    start = (3 - donated) % 3
    return overhang[:start], tuple(overhang[at : at + 3] for at in range(start, len(overhang), 3))


def _termini(
    parts: PartList | None,
    position: str,
    end: End,
    suffix: str,
    codons: Sequence[str],
    usage: CodonUsage,
) -> tuple[Terminus, ...] | None:
    """Return what this junction spells at each member's end, or ``None`` where only a stop fits."""
    spelled = "".join(amino_acid(one) for one in codons)
    count = len(spelled) + (1 if end == "3'" and suffix else 0)
    if parts is None or not count:
        return ()
    made: list[Terminus] = []
    for name, protein in parts.items():
        upper = protein.upper()
        if end == "3'":
            wild = upper[-count:]
            shared = _shared(wild[0], suffix, usage) if suffix else ""
            if shared is None:
                return None
            synthesised = shared + spelled
        else:
            wild = upper[:count]
            synthesised = spelled
        made.append(Terminus(name, position, end, wild, synthesised))
    return tuple(made)


def _shared(wild: str, suffix: str, usage: CodonUsage) -> str | None:
    """Return the residue a codon ending in `suffix` spells for `wild`, or ``None`` for a stop.

    The wild-type residue stands where any of its codons ends in those bases. Where none does,
    the codon this host counts most often for it is overwritten at the bases the overhang fixes,
    which is the smallest change to the DNA that spells the junction.
    """
    if any(codon.endswith(suffix) for codon in _codons_for()[wild]):
        return wild
    preferred = usage.synonymous(_codons_for()[wild][0])[0]
    fitting = [
        codon for codon in _every_codon() if codon.endswith(suffix) and amino_acid(codon) != "*"
    ]
    if not fitting:
        return None
    best = min(
        fitting,
        key=lambda codon: (-_agreement(codon, preferred), -usage.counts[codon], codon),
    )
    return amino_acid(best)


def _agreement(one: str, other: str) -> int:
    """How many positions two codons spell alike."""
    return sum(a == b for a, b in zip(one, other, strict=True))


@cache
def _every_codon() -> tuple[str, ...]:
    """All 64 codons."""
    return tuple("".join(bases) for bases in product("ACGT", repeat=3))


@cache
def _codons_for() -> Mapping[str, tuple[str, ...]]:
    """Each amino acid, and every codon that spells it."""
    families: dict[str, list[str]] = {}
    for codon in _every_codon():
        families.setdefault(amino_acid(codon), []).append(codon)
    return {one: tuple(codons) for one, codons in families.items()}


def _settle(
    readings: Sequence[Sequence[_Reading]],
    order: Sequence[int],
    sites: Sequence[_Site],
    enzyme: Enzyme,
    avoid: Sequence[Enzyme],
    min_distance: int,
    allow_uniform: bool,
) -> tuple[str, ...]:
    """Return the overhang each junction takes, for the lowest total over every part list.

    The search is exhaustive with a bound: every rule is either about one overhang or about a
    pair of them, so a set is allowed exactly when each of its pairs is, and a branch whose
    running total plus the cheapest possible remainder already matches the best set found can
    hold nothing better.

    Raises
    ------
    ValueError
        If no set satisfies the rules at once.
    """
    joins = _joiner(enzyme, avoid, min_distance, allow_uniform)
    floor = [0] * (len(order) + 1)
    for at in reversed(range(len(order))):
        floor[at] = floor[at + 1] + readings[order[at]][0].cost
    best: list[str] | None = None
    lowest = 0

    def walk(at: int, taken: list[str], running: int) -> None:
        nonlocal best, lowest
        if best is not None and running + floor[at] >= lowest:
            return
        if at == len(order):
            best, lowest = list(taken), running
            return
        for reading in readings[order[at]]:
            if joins(reading.overhang, taken):
                taken.append(reading.overhang)
                walk(at + 1, taken, running + reading.cost)
                taken.pop()

    walk(0, [], 0)
    if best is None:
        raise ValueError(
            "no overhang standard satisfies the rules at every junction at once: "
            f"{', '.join(repr(sites[index].name) for index in order)} cannot be settled together"
        )
    chosen = [""] * len(readings)
    for at, index in enumerate(order):
        chosen[index] = best[at]
    return tuple(chosen)


def _joiner(enzyme: Enzyme, avoid: Sequence[Enzyme], min_distance: int, allow_uniform: bool):
    """Return a test of whether a candidate may join a partial set, remembering each pair."""
    seen: dict[tuple[str, str], bool] = {}

    def joins(candidate: str, taken: Sequence[str]) -> bool:
        for other in taken:
            pair = (candidate, other)
            if pair not in seen:
                seen[pair] = (
                    refusal(
                        candidate,
                        enzyme,
                        taken=(other,),
                        avoid=avoid,
                        min_distance=min_distance,
                        allow_uniform=allow_uniform,
                    )
                    is None
                )
            if not seen[pair]:
                return False
        return True

    return joins


def _standard(
    sites: Sequence[_Site],
    readings: Sequence[Sequence[_Reading]],
    order: Sequence[int],
    chosen: Sequence[str],
    enzyme: Enzyme,
    avoid: Sequence[Enzyme],
    min_distance: int,
    allow_uniform: bool,
) -> Standard:
    """Gather the chosen overhangs, the trail of what each junction refused, and the termini."""
    taken: list[str] = []
    trails: dict[int, tuple[Rejection, ...]] = {}
    for index in order:
        rejected: list[Rejection] = []
        for reading in readings[index]:
            if reading.overhang == chosen[index]:
                break
            broken = refusal(
                reading.overhang,
                enzyme,
                taken=taken,
                avoid=avoid,
                min_distance=min_distance,
                allow_uniform=allow_uniform,
            )
            if broken is not None:
                rejected.append(broken)
        trails[index] = tuple(rejected)
        taken.append(chosen[index])
    choices: list[Choice] = []
    termini: list[Terminus] = []
    forced: list[str] = []
    for index, site in enumerate(sites):
        junction = (
            Junction(site.name, overhang=site.pinned)
            if site.pinned is not None
            else Junction(site.name)
        )
        choices.append(Choice(junction, chosen[index], 0, trails[index]))
        taking = next(one for one in readings[index] if one.overhang == chosen[index])
        termini.extend(taking.termini)
        if readings[index][0].cost:
            forced.append(site.name)
    return Standard(
        tuple(chosen[: len(sites) - 1]),
        scar_overhang=chosen[-1],
        choices=tuple(choices),
        termini=tuple(termini),
        forced=tuple(forced),
    )

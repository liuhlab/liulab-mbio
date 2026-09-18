"""The att sites, which one pairs with which, and what a recombination writes at a junction.

Gateway cuts nothing and ligates nothing: two att sites recombine and the reaction rewrites the
sites themselves. So what this module owns is the sites' own bases and one piece of arithmetic.
Strand exchange crosses a 7 bp overlap, so the first `CROSSOVER` bases of a 25 bp recombination
region come from one partner and the rest from the other; BP and LR are that one rule run in
the two directions.

The eight sequences ship, but a record's sites are **found** rather than matched: the vendor
states that the sequences drift between products, so `identify` holds a window to its overlap
and lets the flanks differ. Every value here cites `docs/research/gateway-cloning.md`.

`liulab_mbio.sites` is enzyme cut sites and has nothing to do with this module. Coordinates are
the model's, 0-based and half-open, and a site across the origin of a circular record ends past
the record's length.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

from liulab_mbio.sequence import Segment, SequenceRecord, Strand, reverse_complement

#: The two reactions there are: BP joins attB with attP, LR joins attL with attR back again
#: (Hartley 2000, and MAN0000470 page 5; note §4).
type Reaction = Literal["BP", "LR"]

#: The recombination region every att site shares, and the parts it is made of: a 5 bp flank, a
#: 15 bp core and a 5 bp flank, with a 7 bp overlap inside the core (MAN0000470 page 5 and
#: US 7,670,823 B1; note §3).
REGION_BP = 25
FLANK_BP = 5
CORE_BP = 15
OVERLAP_BP = 7

#: The first base of the 7 bp overlap, 0-based: where strand exchange begins, so where one
#: partner's bases give way to the other's. The overlap is positions 9-15 of the region written
#: in the orientation used here (note §3, §5, §6).
CROSSOVER = 8

#: How far a found site may stand from the sequence shipped for it, outside the overlap. Two att
#: sites of one number stand as few as two bases apart -- attB1 and attL1 differ only in their
#: 5' flank -- so one mismatch is as far as a window may drift and still be nearer one site than
#: any other. This package's choice, measured on `REGIONS`; what the note settles is that a plan
#: tolerates mismatch outside the core and none inside the overlap (note §18).
MISMATCHES = 1

#: The 25 bp recombination region of each att site, 5' to 3' in the site's own orientation, so
#: that site 1 and site 2 line up base for base. From US 7,670,823 B1 FIG. 9, checked against
#: the vendor's own diagrams and against phage lambda J02459 (note §2). The patent carries no
#: copyright notice and every string occurs in deposited GenBank records, which is what lets
#: them ship; no vector carrying them does.
REGIONS: Mapping[str, str] = {
    "attB1": "ACAAGTTTGTACAAAAAAGCAGGCT",
    "attB2": "ACCACTTTGTACAAGAAAGCTGGGT",
    "attP1": "CCAACTTTGTACAAAAAAGCTGAAC",
    "attP2": "CCAACTTTGTACAAGAAAGCTGAAC",
    "attL1": "CCAACTTTGTACAAAAAAGCAGGCT",
    "attL2": "CCAACTTTGTACAAGAAAGCTGGGT",
    "attR1": "ACAAGTTTGTACAAAAAAGCTGAAC",
    "attR2": "ACCACTTTGTACAAGAAAGCTGAAC",
}

#: Which site reacts with which: a site reacts only with its partner of the same number, which
#: is what holds the moved segment's direction (MAN0000470 page 5; note §4).
PARTNERS: Mapping[str, str] = {
    "attB1": "attP1",
    "attP1": "attB1",
    "attB2": "attP2",
    "attP2": "attB2",
    "attL1": "attR1",
    "attR1": "attL1",
    "attL2": "attR2",
    "attR2": "attL2",
}

#: What each reaction takes: the kind of site on the record whose segment moves, then the kind
#: on the record whose backbone takes it. BP moves an attB substrate into a donor vector; LR
#: moves an entry clone into a destination vector (note §4, §5, §6).
SUBSTRATES: Mapping[Reaction, tuple[str, str]] = {"BP": ("attB", "attP"), "LR": ("attL", "attR")}

#: The two overlaps there are, one per site number. Site 1 and site 2 differ inside the overlap,
#: which is why it is matched exactly and the flanks are not.
_OVERLAPS = frozenset(bases[CROSSOVER : CROSSOVER + OVERLAP_BP] for bases in REGIONS.values())
_REVERSE_OVERLAPS = frozenset(reverse_complement(bases) for bases in _OVERLAPS)

#: Where a reverse-strand site's overlap falls, read on the top strand.
_REVERSE_AT = REGION_BP - CROSSOVER - OVERLAP_BP


@dataclass(frozen=True, slots=True)
class AttSite:
    """One att site found in a record: which site it is, where it stands, and what it spells.

    Parameters
    ----------
    name
        Which of `REGIONS` its bases read as.
    start, end
        The 25 bp recombination region, on the top strand. `end` passes the record's length
        where the region runs across the origin.
    strand
        `Strand.FORWARD` where the region reads along the top strand.
    bases
        The 25 bases in the site's own orientation, as the record spells them, which is what a
        drifted site is carried into a product by.
    mismatches
        How far those bases stand from `REGIONS[name]`, all of it outside the 7 bp overlap.
    """

    name: str
    start: int
    end: int
    strand: Strand
    bases: str
    mismatches: int

    @property
    def kind(self) -> str:
        """``attB``, ``attP``, ``attL`` or ``attR``."""
        return self.name[:-1]

    @property
    def number(self) -> int:
        """1 or 2. A site reacts only with its own number, which fixes the segment's direction."""
        return int(self.name[-1])

    @property
    def span(self) -> Segment:
        """The recombination region, as a segment of the record it was found in."""
        return Segment(self.start, self.end)


def core(bases: str) -> str:
    """Return the 15 bp core of a recombination region, which every site of a number shares.

    Examples
    --------
    >>> core(REGIONS["attB1"]) == core(REGIONS["attP1"])
    True
    """
    return bases[FLANK_BP : FLANK_BP + CORE_BP]


def overlap(bases: str) -> str:
    """Return the 7 bp overlap of a recombination region, where strand exchange takes place.

    Examples
    --------
    >>> overlap(REGIONS["attB1"]), overlap(REGIONS["attB2"])
    ('GTACAAA', 'GTACAAG')
    """
    return bases[CROSSOVER : CROSSOVER + OVERLAP_BP]


def joined(before: str, after: str) -> str:
    """Return what a junction spells where `before`'s molecule gives way to `after`'s.

    The whole of both reactions: the first `CROSSOVER` bases of the region come from `before`
    and the rest from `after`, because that is where the 7 bp overlap begins and where the two
    strands are exchanged (note §3, §5, §6). Both regions are written in their own orientation.

    Examples
    --------
    >>> joined(REGIONS["attR1"], REGIONS["attL1"]) == REGIONS["attB1"]
    True
    """
    return before[:CROSSOVER] + after[CROSSOVER:]


def writes(before: str, after: str) -> str:
    """Return the name of the site the junction between these two att sites spells.

    LR run at site 1 is ``writes("attR1", "attL1") == "attB1"`` on the expression clone and
    ``writes("attL1", "attR1") == "attP1"`` on the by-product. That by-product's attP is
    shorter than the donor vector's: attR was built with 33 bp of the arm deleted, which is
    what makes the reaction irreversible (note §3, §6). Only the 25 bp region is named here.

    Raises
    ------
    KeyError
        If either name is not an att site, if the two do not react, or if they spell a region
        no att site carries.

    Examples
    --------
    >>> writes("attP1", "attB1")
    'attL1'
    """
    if PARTNERS[before] != after:
        raise KeyError(f"{before} reacts with {PARTNERS[before]}, not with {after}")
    bases = joined(REGIONS[before], REGIONS[after])
    for name, region in REGIONS.items():
        if region == bases:
            return name
    raise KeyError(f"{before} and {after} spell no att site: {bases}")


def identify(bases: str, *, names: Sequence[str] = (), mismatches: int = MISMATCHES) -> str | None:
    """Return which att site these 25 bases read as, or ``None`` where none is nearer than another.

    The 7 bp overlap has to match exactly: it is what tells site 1 from site 2, and the patent
    reports that a substitution inside it changes which site the region recombines with. The
    flanks and the rest of the core may differ by up to `mismatches`, because the vendor states
    that the sequences drift between products (note §2, §3, §18). Two att sites can stand as few
    as two bases apart, so a window as near one as another is named by neither.

    Parameters
    ----------
    bases
        A 25 bp recombination region, in the site's own orientation.
    names
        Which of `REGIONS` to consider; all of them when empty.
    mismatches
        How far the bases outside the overlap may differ.

    Raises
    ------
    ValueError
        If `bases` is not `REGION_BP` long.

    Examples
    --------
    >>> identify(REGIONS["attL2"])
    'attL2'
    """
    if len(bases) != REGION_BP:
        raise ValueError(f"an att recombination region is {REGION_BP} bases, got {len(bases)}")
    found = overlap(bases)
    scored = sorted(
        (_differences(bases, REGIONS[name]), name)
        for name in (names or tuple(REGIONS))
        if overlap(REGIONS[name]) == found
    )
    if not scored or scored[0][0] > mismatches:
        return None
    if len(scored) > 1 and scored[1][0] == scored[0][0]:
        return None
    return scored[0][1]


def find_att_sites(
    record: SequenceRecord, *, names: Sequence[str] = (), mismatches: int = MISMATCHES
) -> tuple[AttSite, ...]:
    """Return every att site in `record`, read on either strand, the first one first.

    A circular record is read across its origin as well, so a site there is found once, with
    its `end` past the record's length.

    Parameters
    ----------
    record
        The record to read.
    names, mismatches
        As `identify`.
    """
    wanted = tuple(names) if names else tuple(REGIONS)
    haystack = record.sequence
    if record.topology == "circular":
        haystack += record.sequence[: REGION_BP - 1]
    found = []
    for start in range(len(haystack) - REGION_BP + 1):
        window = haystack[start : start + REGION_BP]
        candidates = []
        if overlap(window) in _OVERLAPS:
            candidates.append((Strand.FORWARD, window))
        if window[_REVERSE_AT : _REVERSE_AT + OVERLAP_BP] in _REVERSE_OVERLAPS:
            candidates.append((Strand.REVERSE, reverse_complement(window)))
        for strand, bases in candidates:
            name = identify(bases, names=wanted, mismatches=mismatches)
            if name is not None:
                found.append(
                    AttSite(
                        name,
                        start,
                        start + REGION_BP,
                        strand,
                        bases,
                        _differences(bases, REGIONS[name]),
                    )
                )
    return tuple(found)


def att_pair(
    record: SequenceRecord, kind: str, *, mismatches: int = MISMATCHES
) -> tuple[AttSite, AttSite]:
    """Return `record`'s site 1 and site 2 of this kind, in that order.

    A Gateway substrate carries two sites of one kind pointing at each other across the segment
    that moves, so a pair is what a reaction needs and one site is not. The two must lie on
    opposite strands: a site reacts only with its own number, and no source read for the note
    reports what a site the other way round does (note §4, and its open gaps).

    Raises
    ------
    KeyError
        If `kind` is not one this module holds a pair for.
    ValueError
        If the record does not carry one site of each number, or if the two do not point at
        each other.
    """
    wanted = tuple(f"{kind}{number}" for number in (1, 2))
    for name in wanted:
        if name not in REGIONS:
            raise KeyError(f"no att site called {name!r}; this module holds {', '.join(REGIONS)}")
    found = find_att_sites(record, names=wanted, mismatches=mismatches)
    pair = []
    for name in wanted:
        sites = [site for site in found if site.name == name]
        if len(sites) != 1:
            raise ValueError(
                f"{_named(record)} carries {len(sites)} {name} sites, needing one: looked for "
                f"{' and '.join(wanted)} on either strand, matching the {OVERLAP_BP} bp overlap "
                f"exactly and the rest of the {REGION_BP} bp region within {mismatches} base(s)"
            )
        pair.append(sites[0])
    first, second = pair
    if first.strand == second.strand:
        raise ValueError(
            f"{_named(record)} carries {first.name} and {second.name} on the same strand, so "
            "they do not point at each other across a segment that could move"
        )
    return first, second


def _named(record: SequenceRecord) -> str:
    """Return what to call a record in a refusal."""
    return record.name or "this record"


def _differences(one: str, other: str) -> int:
    """Count where two strings of one length disagree."""
    return sum(first != second for first, second in zip(one, other, strict=True))

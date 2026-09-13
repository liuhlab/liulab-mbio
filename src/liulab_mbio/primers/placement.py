"""Where a primer may anneal on a template, where it does, where else it primes, and amplicons.

A site or an amplicon crosses the origin of a circular template, ending past its length.
"""

from dataclasses import dataclass

from liulab_mbio.primers.thresholds import THRESHOLDS, Thresholds
from liulab_mbio.sequence import (
    BindingSite,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)


@dataclass(frozen=True, slots=True)
class Placement:
    """Where a primer's binding site may lie: the positions each of its two ends may take.

    A tail fixes the 5' end, a target to read across bounds each end's distance from it, and a
    region bounds the amplicon. A position is a boundary, 0-based and half-open as everywhere:
    a forward site's 5' end is its start and its 3' end its end, and a reverse site's the other
    way round, which is the position `design_primer` reads. A span runs past the length of a
    circular template when it crosses the origin.

    Parameters
    ----------
    five_prime
        The positions its 5' end may take. `Segment(442, 463)` allows 442 to 462. ``None``
        leaves it wherever a length puts it.
    three_prime
        The positions its 3' end may take, the same way.

    Raises
    ------
    ValueError
        If neither end is bounded, which would say nothing at all.

    Examples
    --------
    >>> template = SequenceRecord("ACGT" * 12)
    >>> anchored = Placement(five_prime=Segment(0, 1))
    >>> anchored.allows(BindingSite(0, 20, Strand.FORWARD), template)
    True
    >>> anchored.allows(BindingSite(1, 21, Strand.FORWARD), template)
    False
    """

    five_prime: Segment | None = None
    three_prime: Segment | None = None

    def __post_init__(self) -> None:
        """Refuse a placement that bounds neither end."""
        if self.five_prime is None and self.three_prime is None:
            raise ValueError("a placement bounds at least one end of a binding site")

    def allows(self, site: BindingSite, template: SequenceRecord) -> bool:
        """Return whether both ends of a binding site lie where this placement puts them."""
        ends = (site.start, site.end) if site.strand is Strand.FORWARD else (site.end, site.start)
        return all(
            span is None or _holds(span, position, template)
            for span, position in zip((self.five_prime, self.three_prime), ends, strict=True)
        )


@dataclass(frozen=True, slots=True)
class PrimingSite:
    """Where a primer's 3' end can anneal on a template, and how strongly.

    Parameters
    ----------
    site
        What the annealing region covers, running past the end of a circular template when it
        crosses the origin.
    tm
        Of the primer annealed there, °C at primer3's default conditions.
    mismatches
        Bases of the annealing region that do not pair.
    """

    site: BindingSite
    tm: float
    mismatches: int


def find_binding_sites(
    sequence: str, template: SequenceRecord, *, thresholds: Thresholds = THRESHOLDS
) -> tuple[BindingSite, ...]:
    """Return every place a primer's 3' end matches a template exactly, on either strand.

    The match runs from the 3' end back towards the 5' end, so a tail hangs off it; it must
    reach `Thresholds.binding_min_length`. A site crosses the origin of a circular template.

    Examples
    --------
    >>> template = SequenceRecord("CGGCGTAATCATGGTCATAGCTGTTTCC")
    >>> sites = find_binding_sites("TTGGTCTCAGGCGTAATCATGGTCATAGC", template)
    >>> [(site.start, site.end, site.strand.name) for site in sites]
    [(1, 21, 'FORWARD')]
    """
    dna = sequence.upper()
    length = len(template)
    reverse = reverse_complement(dna)
    sites = []
    for index in range(length):
        matched = _matched(template, index, dna[::-1], -1)
        if matched >= thresholds.binding_min_length:
            start = (index - matched + 1) % length
            sites.append(BindingSite(start, start + matched, Strand.FORWARD))
        matched = _matched(template, index, reverse, 1)
        if matched >= thresholds.binding_min_length:
            sites.append(BindingSite(index, index + matched, Strand.REVERSE))
    return tuple(sorted(sites, key=lambda site: (site.start, site.strand)))


def find_priming_sites(
    sequence: str, template: SequenceRecord, *, thresholds: Thresholds = THRESHOLDS
) -> tuple[PrimingSite, ...]:
    """Return every place an annealing region can prime a template, on either strand.

    A place counts when few enough of its bases mismatch, fewest of all at the 3' end, and
    when the primer annealed there melts within `Thresholds.off_target_margin` of a perfect
    match. Sites cross the origin of a circular template.
    """
    import primer3

    dna = sequence.upper()
    size = len(dna)
    length = len(template)
    if size > length:
        return ()
    circular = template.topology == "circular"
    top = template.sequence + (template.sequence[: size - 1] if circular else "")
    reverse = reverse_complement(dna)
    edge = min(size, thresholds.off_target_3prime_window)
    floor = primer3.calc_end_stability(dna, reverse).tm - thresholds.off_target_margin
    ends = tuple(
        _near(window, set(top), thresholds.off_target_3prime_mismatches)
        for window in (dna[-edge:], reverse[:edge])
    )
    found = []
    for start in range(length if circular else length - size + 1):
        forward_end = ends[0] is None or top[start + size - edge : start + size] in ends[0]
        reverse_end = ends[1] is None or top[start : start + edge] in ends[1]
        if not (forward_end or reverse_end):
            continue
        here = top[start : start + size]
        for strand, probe, anchored, annealed in (
            (Strand.FORWARD, dna, forward_end, None),
            (Strand.REVERSE, reverse, reverse_end, here),
        ):
            if not anchored:
                continue
            mismatches = _mismatches(probe, here)
            if mismatches > thresholds.off_target_mismatches:
                continue
            tm = primer3.calc_end_stability(dna, annealed or reverse_complement(here)).tm
            if tm >= floor:
                found.append(PrimingSite(BindingSite(start, start + size, strand), tm, mismatches))
    return tuple(found)


def amplicon_sizes(
    forward: Primer,
    reverse: Primer,
    template: SequenceRecord,
    *,
    thresholds: Thresholds = THRESHOLDS,
) -> tuple[int, ...]:
    """Return the size of every amplicon a pair can make on a template, smallest first.

    Tails count, as they do in `PairReport.amplicon_length`. A primer carrying binding sites is
    taken to bind where they say, so asking about a template other than the one it was placed on
    means clearing them first.
    """
    placed: list[tuple[BindingSite, int]] = []
    for primer in (forward, reverse):
        sites = primer.binding_sites or find_binding_sites(
            primer.sequence, template, thresholds=thresholds
        )
        annealing = max((site.end - site.start for site in sites), default=len(primer.sequence))
        placed.extend((site, len(primer.sequence) - annealing) for site in sites)
    products = []
    for site, tail in placed:
        if site.strand is not Strand.FORWARD:
            continue
        for other, other_tail in placed:
            if other.strand is not Strand.REVERSE:
                continue
            span = _span(site, other, template)
            if span is not None:
                products.append(span + tail + other_tail)
    return tuple(sorted(products))


def _holds(span: Segment, position: int, template: SequenceRecord) -> bool:
    """Return whether a span of positions holds one, counting round a circular template."""
    if template.topology == "circular":
        return (position - span.start) % len(template) < span.end - span.start
    return span.start <= position < span.end


def _matched(template: SequenceRecord, index: int, probe: str, step: int) -> int:
    """Return how many bases of `probe` match the template from `index`, walking by `step`."""
    length = len(template)
    circular = template.topology == "circular"
    matched = 0
    while matched < min(len(probe), length):
        position = index + step * matched
        if not circular and not 0 <= position < length:
            break
        if template.sequence[position % length] != probe[matched]:
            break
        matched += 1
    return matched


def _near(window: str, alphabet: set[str], mismatches: int) -> frozenset[str] | None:
    """Return every string `mismatches` substitutions or fewer from `window`, over `alphabet`.

    ``None`` stands for every string of that length, which is no filter at all. A template is
    scanned once per base, so testing its 3' window against this set is what keeps a design's
    search of every length affordable.
    """
    if mismatches >= len(window):
        return None
    near = {window}
    for _ in range(mismatches):
        near |= {
            one[:index] + base + one[index + 1 :]
            for one in near
            for index in range(len(one))
            for base in alphabet
        }
    return frozenset(near)


def _mismatches(one: str, other: str) -> int:
    return sum(base != base_here for base, base_here in zip(one, other, strict=True))


def _span(forward: BindingSite, reverse: BindingSite, template: SequenceRecord) -> int | None:
    length = len(template)
    if template.topology == "circular":
        return (reverse.end - forward.start) % length or length
    span = reverse.end - forward.start
    return span if span > 0 else None

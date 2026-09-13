"""Where a primer anneals on a template, where else it can prime, and what a pair amplifies.

A site or an amplicon crosses the origin of a circular template, ending past its length.
"""

from dataclasses import dataclass

from liulab_mbio.primers.thresholds import THRESHOLDS, Thresholds
from liulab_mbio.sequence import BindingSite, Primer, SequenceRecord, Strand, reverse_complement


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
    window = thresholds.off_target_3prime_window
    floor = primer3.calc_end_stability(dna, reverse).tm - thresholds.off_target_margin
    found = []
    for start in range(length if circular else length - size + 1):
        here = top[start : start + size]
        for strand, probe, anchor, annealed in (
            (Strand.FORWARD, dna, _mismatches(dna[-window:], here[-window:]), None),
            (Strand.REVERSE, reverse, _mismatches(reverse[:window], here[:window]), here),
        ):
            if anchor > thresholds.off_target_3prime_mismatches:
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


def _mismatches(one: str, other: str) -> int:
    return sum(base != base_here for base, base_here in zip(one, other, strict=True))


def _span(forward: BindingSite, reverse: BindingSite, template: SequenceRecord) -> int | None:
    length = len(template)
    if template.topology == "circular":
        return (reverse.end - forward.start) % length or length
    span = reverse.end - forward.start
    return span if span > 0 else None

"""Where a primer may anneal on a template, where it does, where else it primes, and amplicons.

A site or an amplicon crosses the origin of a circular template, ending past its length.
"""

from dataclasses import dataclass
from functools import lru_cache

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

    A plan places the same primers on the same records again and again, so where a sequence
    binds is remembered: it depends on nothing but the sequence, the template and the thresholds.
    """
    return _binding_sites(sequence.upper(), template.sequence, template.topology, thresholds)


@lru_cache(maxsize=1 << 12)
def _binding_sites(
    dna: str, sequence: str, topology: str, thresholds: Thresholds
) -> tuple[BindingSite, ...]:
    length = len(sequence)
    circular = topology == "circular"
    reverse = reverse_complement(dna)
    sites = []
    for index in range(length):
        matched = _matched(sequence, circular, index, dna[::-1], -1)
        if matched >= thresholds.binding_min_length:
            start = (index - matched + 1) % length
            sites.append(BindingSite(start, start + matched, Strand.FORWARD))
        matched = _matched(sequence, circular, index, reverse, 1)
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

    A design tests thousands of candidates against one template, so what a sequence primes there
    is remembered: it depends on nothing but the sequence, the template and the thresholds.
    """
    return _priming_sites(sequence.upper(), template.sequence, template.topology, thresholds)


def find_off_target_sites(
    sequence: str,
    template: SequenceRecord,
    placed: tuple[BindingSite, ...],
    *,
    thresholds: Thresholds = THRESHOLDS,
) -> tuple[PrimingSite, ...]:
    """Return where an annealing region primes a template other than where it is placed.

    What `find_priming_sites` returns, less every site sharing its 3' end with one of `placed`.
    Dropping those first is what a primer is judged on, and the thermodynamics of a place no
    one asked about is never worked out.

    Examples
    --------
    >>> template = SequenceRecord("CGGCGTAATCATGGTCATAGCTGTTTCC")
    >>> intended = find_binding_sites("GGCGTAATCATGGTCATAGC", template)
    >>> find_off_target_sites("GGCGTAATCATGGTCATAGC", template, intended)
    ()
    """
    dna = sequence.upper()
    bases, topology = template.sequence, template.topology
    intended = {_three_prime_end(site, len(bases)) for site in placed}
    return _primed(
        dna,
        bases,
        topology,
        thresholds,
        tuple(
            annealing
            for annealing in _annealings(dna, bases, topology, thresholds)
            if _three_prime_end(annealing[0], len(bases)) not in intended
        ),
    )


def _three_prime_end(site: BindingSite, length: int) -> tuple[int, Strand]:
    """Return the base a site reads from and its strand, which is where a primer extends."""
    end = site.start if site.strand is Strand.REVERSE else site.end - 1
    return end % length, site.strand


@lru_cache(maxsize=1 << 15)
def _perfect_stability(dna: str) -> float:
    """Return the end stability of a sequence annealed to its own complement.

    It depends on the sequence alone, so the templates a sequence is searched against share
    one value.
    """
    import primer3

    return primer3.calc_end_stability(dna, reverse_complement(dna)).tm


@lru_cache(maxsize=1 << 12)
def _priming_sites(
    dna: str, sequence: str, topology: str, thresholds: Thresholds
) -> tuple[PrimingSite, ...]:
    annealings = _annealings(dna, sequence, topology, thresholds)
    return _primed(dna, sequence, topology, thresholds, annealings)


@lru_cache(maxsize=1 << 15)
def _annealings(
    dna: str, sequence: str, topology: str, thresholds: Thresholds
) -> tuple[tuple[BindingSite, int], ...]:
    """Return every site whose 3' window anchors and whose bases mismatch within the cap.

    Each with how many of its bases mismatch, by start. What a site melts at is settled apart
    from this, so a caller dropping a site pays primer3 nothing for it.
    """
    size = len(dna)
    length = len(sequence)
    if size > length:
        return ()
    circular = topology == "circular"
    top = sequence + (sequence[: size - 1] if circular else "")
    reverse = reverse_complement(dna)
    edge = min(size, thresholds.off_target_3prime_window)
    starts = range(length if circular else length - size + 1)
    ends = tuple(
        _anchored(
            window,
            sequence,
            circular,
            edge,
            offset,
            thresholds.off_target_3prime_mismatches,
            starts,
        )
        for window, offset in ((dna[-edge:], size - edge), (reverse[:edge], 0))
    )
    found = []
    for start in sorted(ends[0] | ends[1]):
        here = top[start : start + size]
        for strand, probe, anchored in (
            (Strand.FORWARD, dna, start in ends[0]),
            (Strand.REVERSE, reverse, start in ends[1]),
        ):
            if not anchored:
                continue
            mismatches = _mismatches(probe, here, thresholds.off_target_mismatches)
            if mismatches <= thresholds.off_target_mismatches:
                found.append((BindingSite(start, start + size, strand), mismatches))
    return tuple(found)


def _primed(
    dna: str,
    sequence: str,
    topology: str,
    thresholds: Thresholds,
    annealings: tuple[tuple[BindingSite, int], ...],
) -> tuple[PrimingSite, ...]:
    """Return those annealing places melting within the margin of a perfect match."""
    if not annealings:
        return ()
    import primer3

    perfect = _perfect_stability(dna)
    floor = perfect - thresholds.off_target_margin
    found = []
    for site, mismatches in annealings:
        # A perfect match anneals exactly as the floor's duplex does.
        if mismatches:
            here = _covered(sequence, topology, site)
            annealed = here if site.strand is Strand.REVERSE else reverse_complement(here)
            tm = primer3.calc_end_stability(dna, annealed).tm
        else:
            tm = perfect
        if tm >= floor:
            found.append(PrimingSite(site, tm, mismatches))
    return tuple(found)


def _covered(sequence: str, topology: str, site: BindingSite) -> str:
    """Return the bases of a template a site covers, wrapping once round a circular origin."""
    if topology == "circular" and site.end > len(sequence):
        return sequence[site.start :] + sequence[: site.end - len(sequence)]
    return sequence[site.start : site.end]


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


def _matched(sequence: str, circular: bool, index: int, probe: str, step: int) -> int:
    """Return how many bases of `probe` match a template from `index`, walking by `step`."""
    length = len(sequence)
    matched = 0
    while matched < min(len(probe), length):
        position = index + step * matched
        if not circular and not 0 <= position < length:
            break
        if sequence[position % length] != probe[matched]:
            break
        matched += 1
    return matched


@lru_cache(maxsize=8)
def _alphabet(sequence: str) -> frozenset[str]:
    """Return the letters a template holds, which every candidate tested against it asks for."""
    return frozenset(sequence)


def _near(window: str, alphabet: frozenset[str], mismatches: int) -> frozenset[str] | None:
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


def _anchored(
    window: str,
    sequence: str,
    circular: bool,
    edge: int,
    offset: int,
    mismatches: int,
    starts: range,
) -> set[int]:
    """Return every start whose 3' window can match the template, `offset` bases into the site.

    A design judges thousands of candidates against one template, so the template is indexed
    once and each candidate looks up only the few windows its own can match, rather than walking
    every position itself.
    """
    hits = _hits(window, sequence, edge, circular, mismatches)
    if hits is None:
        return set(starts)
    length = len(sequence)
    if circular:
        return {(position - offset) % length for position in hits}
    return {position - offset for position in hits if position - offset in starts}


@lru_cache(maxsize=1 << 12)
def _hits(
    window: str, sequence: str, edge: int, circular: bool, mismatches: int
) -> tuple[int, ...] | None:
    """Return where a template's own windows lie within `mismatches` of this one.

    ``None`` stands for everywhere. Candidates of every length at one 3' end share a window,
    so the lookup is remembered.
    """
    near = _near(window, _alphabet(sequence), mismatches)
    if near is None:
        return None
    index = _windows(sequence, edge, circular)
    return tuple(position for probe in near for position in index.get(probe, ()))


@lru_cache(maxsize=8)
def _windows(sequence: str, edge: int, circular: bool) -> dict[str, tuple[int, ...]]:
    """Return where every window of `edge` bases starts, wrapping one round the origin."""
    scan = sequence + (sequence[: edge - 1] if circular else "")
    index: dict[str, list[int]] = {}
    for position in range(len(scan) - edge + 1):
        index.setdefault(scan[position : position + edge], []).append(position)
    return {window: tuple(positions) for window, positions in index.items()}


def _mismatches(one: str, other: str, cap: int) -> int:
    """Return how many bases differ, counting no further than one past `cap`.

    A count over the cap is only ever compared with it, and most places a primer is tested
    against differ within the first few bases.
    """
    count = 0
    for base, base_here in zip(one, other, strict=True):
        if base != base_here:
            count += 1
            if count > cap:
                break
    return count


def _span(forward: BindingSite, reverse: BindingSite, template: SequenceRecord) -> int | None:
    length = len(template)
    if template.topology == "circular":
        return (reverse.end - forward.start) % length or length
    span = reverse.end - forward.start
    return span if span > 0 else None

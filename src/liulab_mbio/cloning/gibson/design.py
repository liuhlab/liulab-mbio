"""The overlap a Gibson junction takes: how long it is, and which part it is taken from.

An **overlap** is the bases two fragments both spell where they meet. One of them already spells
them; the other carries them as a 5' tail on its primer, so the junction adds nothing and reads
exactly what the two parts read.

Two rules, both the supplier's, through ``docs/research/gibson-assembly.md``:

- **Length and melting temperature are one rule, not two** (note §3 and §4). The product
  documents a band, and NEB's worked example lengthens the overlap inside it until the Wallace
  melting temperature reaches the floor: "The length of the overlap sequence is determined by
  the number of nucleotides needed to reach a Tm >= 48 C". So `overlap_before` and
  `overlap_after` return the shortest overlap in the band that reaches it, and the longest the
  band allows where nothing in it does.
- **Where it is taken from** (note §3). Both NEB manuals let either fragment carry it, and both
  ask for the whole overlap to come from the vector where the vector is PCR-generated and
  reused: "the entire overlap sequence must originate from the vector sequence and must be
  added to primers that will be used to amplify the insert". So a junction with the vector on
  one side takes its bases from the vector, and the insert tails them.

Coordinates are the model's, 0-based and half-open, and an overlap read off a circular vector
may run across the origin.

`overlap_checks` then judges what was chosen. Which of those carry a verdict is the note's
answer and not this module's: length and melting temperature are the supplier's own bands, a
palindrome and a repeat at an overlap's end are rules both NEB manuals state, and GC content and
two overlaps being alike are measured and left unjudged because nobody quantifies either (§10,
§18).

Two jobs deserve neither a synthesis order nor a PCR, and this module lays out the oligos for
both. `stitch_oligos` tiles a short part with overlapping oligos on both strands, and
`bridging_oligo` joins two parts sharing no homology at all. Each takes single-stranded oligos
into the reaction, which only some assembly products document (§14, §15), so each refuses a
product whose supplier does not.
"""

from collections.abc import Iterator, Sequence

from liulab_mbio.checks import Check
from liulab_mbio.cloning.gibson.assembly import Bridge, Part, StitchOligo
from liulab_mbio.cloning.gibson.bench import ASSEMBLY_PRODUCTS, AssemblyProduct, OverlapRule
from liulab_mbio.sequence import Segment, SequenceRecord, Strand, reverse_complement

#: The Wallace rule's degrees per base pair, which is the only way either NEB manual names for
#: an overlap's melting temperature: "AT pair = 2 C and GC pair = 4 C" (note §4).
WALLACE_AT = 2.0
WALLACE_GC = 4.0

#: Bases of self-complementary sequence that make an overlap "highly palindromic", which NEB
#: asks a design to avoid and puts at up to a tenfold loss of colonies (§10). NEB quantifies the
#: loss and not the trigger, so the floor is this package's: five pairs of stem. Lower fires on
#: the single restriction site an overlap may legitimately carry.
PALINDROME_BASES = 10

#: What counts as a repeating sequence at an overlap's end, which both NEB manuals refuse: bases
#: covered by a tandem repeat of a unit this long. NEB's own example is a His-tag, three copies
#: of a codon, and the floor is that example in bases (§10); the unit band is this package's.
REPEAT_UNIT_BP = (2, 6)
REPEAT_BASES = 9

#: Bases NEB asks to be left outside such a repeat at each end: "you must flank the His-tag
#: sequence on both sides with at least 2 nucleotides that are not part of the His-tag repeating
#: sequence" (§10).
REPEAT_FLANK = 2

#: Gibson 2011's own stitching oligos: 60 bases each, neighbours overlapping by 20 bp (§14). An
#: oligo is shortened where the part cannot be tiled evenly by full-length ones; the overlap is
#: not.
STITCH_OLIGO_BASES = 60
STITCH_OVERLAP_BP = 20

#: The most oligos Gibson 2011 puts in one reaction: "only eight to twelve 60-base oligos are
#: assembled at one time", for the error rate of chemical synthesis rather than for chemistry
#: (§14). NEB's own FAQ allows the same twelve.
STITCH_OLIGOS = 12

#: The part sizes Addgene puts this route between: "too long to include on overlapping PCR
#: primers (>60 bp) but too short to make its own part (<150 bp)" (§14). The only size window
#: any source states, and a recommendation rather than a limit, so it warns and does not refuse.
STITCH_WINDOW_BP = (60, 150)

#: What NEB asks for of each oligo in the reaction: "45 nM of each oligonucleotide that is less
#: than or equal to twelve 60-base oligonucleotides containing 30-base overlaps" (§14).
STITCH_OLIGO_NM = 45.0

#: Bases a bridging oligo takes from each side. De Saeger measured 20 bp a side at 80-100%
#: correct clones and the Addgene protocol uses the same; NEB's own protocol uses 25 and
#: Thermo's 15, and nothing below 15 is recommended by anyone (§15).
BRIDGE_HOMOLOGY_BP = 20

#: What NEB's own bridging protocol puts in one reaction: 5 µL of the oligo at 0.2 µM, which is
#: this many picomoles against 30 ng of linearised vector (§15).
BRIDGE_OLIGO_PMOL = 1.0

#: What either route's oligos are held to. Both suppliers say the same: "Standard, desalted
#: primers may be used" (NEB), "Gel or HPLC purification of oligonucleotides is not required"
#: (Takara) (§14, §18).
OLIGO_PURITY = "Standard desalted; no gel or HPLC purification is needed."

#: Why such an oligo's sheet row carries no verdict: no source gives one a melting
#: temperature, a hairpin or a dimer threshold, and the three roles
#: `liulab_mbio.primers.thresholds` knows are a primer's (§18).
NO_VERDICT = (
    "It primes nothing, so none of the thresholds a primer is judged by measures it and this "
    "row carries no verdict."
)


def wallace_tm(bases: str) -> float:
    """Return the melting temperature of `bases` by the Wallace rule, °C.

    Two degrees for each A or T and four for each G or C. A base that is neither adds nothing,
    because the rule has nothing to say about one.

    Examples
    --------
    >>> wallace_tm("GCGCGCGCGCGCGCG")
    60.0
    >>> wallace_tm("ATATATATATATATA")
    30.0
    """
    upper = bases.upper()
    return WALLACE_AT * sum(upper.count(base) for base in "AT") + WALLACE_GC * sum(
        upper.count(base) for base in "GC"
    )


def overlap_before(record: SequenceRecord, at: int, rule: OverlapRule) -> str:
    """Return the overlap ending at `at`, taken from `record`'s own bases before that point.

    Raises
    ------
    ValueError
        If `record` holds fewer bases than the rule's longest overlap.

    Examples
    --------
    >>> rule = OverlapRule(4, 6, 18.0)
    >>> overlap_before(SequenceRecord("AAAAGGCCTTTT"), 8, rule)
    'AGGCC'
    """
    return _chosen((_window(record, at - size, size) for size in _sizes(record, rule)), rule)


def overlap_after(record: SequenceRecord, at: int, rule: OverlapRule) -> str:
    """Return the overlap beginning at `at`, taken from `record`'s own bases from that point.

    Raises
    ------
    ValueError
        If `record` holds fewer bases than the rule's longest overlap.

    Examples
    --------
    >>> rule = OverlapRule(4, 6, 18.0)
    >>> overlap_after(SequenceRecord("AAAAGGCCTTTT"), 3, rule)
    'AGGCC'
    """
    return _chosen((_window(record, at, size) for size in _sizes(record, rule)), rule)


def _sizes(record: SequenceRecord, rule: OverlapRule) -> range:
    """Return the overlap lengths the rule allows, shortest first.

    Raises
    ------
    ValueError
        If the record cannot spell the longest of them.
    """
    if len(record) < rule.longest:
        raise ValueError(
            f"{record.name or 'this record'} is {len(record)} bases, too short for the "
            f"{rule.shortest}-{rule.longest} bp overlap this assembly product documents"
        )
    return range(rule.shortest, rule.longest + 1)


def _window(record: SequenceRecord, start: int, size: int) -> str:
    """Return `size` bases of `record` from `start`, reading across the origin where it must.

    Raises
    ------
    ValueError
        If a linear record runs out before the overlap does.
    """
    if record.topology == "circular":
        first = start % len(record)
    elif start >= 0 and start + size <= len(record):
        first = start
    else:
        raise ValueError(
            f"an overlap of {size} bases does not fit {record.name or 'this record'} at {start}"
        )
    return record.extract(Segment(first, first + size))


def _chosen(candidates: Iterator[str], rule: OverlapRule) -> str:
    """Return the shortest candidate reaching the rule's floor, or the longest one offered.

    A product stating no melting temperature takes the shortest its band allows, which is what
    its manual asks for.
    """
    taken = ""
    for bases in candidates:
        taken = bases
        if rule.tm_floor is None or wallace_tm(bases) >= rule.tm_floor:
            return bases
    return taken


def stitch_limit_bp(
    *, oligo_bases: int = STITCH_OLIGO_BASES, overlap_bp: int = STITCH_OVERLAP_BP
) -> int:
    """Return the longest part `STITCH_OLIGOS` oligos of that shape tile.

    *n* oligos of length *L* overlapping by *v* tile *n*(*L*-*v*) + *v* bases, the arithmetic
    that reproduces Gibson's worked 340 bp from eight 60-mers at 20 bp (§14).

    Examples
    --------
    >>> stitch_limit_bp()
    500
    """
    return STITCH_OLIGOS * (oligo_bases - overlap_bp) + overlap_bp


def stitch_oligos(
    bases: str,
    *,
    name: str = "",
    product: AssemblyProduct,
    oligo_bases: int = STITCH_OLIGO_BASES,
    overlap_bp: int = STITCH_OVERLAP_BP,
) -> tuple[StitchOligo, ...]:
    """Lay `bases` out as overlapping oligos tiling both strands, to assemble rather than amplify.

    Neighbours alternate strands and overlap by exactly `overlap_bp`, so the set anneals into
    the whole molecule in the same isothermal reaction as everything else -- Gibson 2011 anneals
    them in no separate step (§14). As few oligos as will do it are used, and none is longer
    than `oligo_bases`; where they cannot all be that long the length is spread evenly, because
    the overlap is the number the sources state and the length is the number they round.

    Parameters
    ----------
    bases
        The whole molecule the oligos make, the junction overlaps at its two ends included.
    name
        What the part is called; each oligo is numbered after it.
    product
        The assembly product on the bench, which has to take single-stranded oligos at all.
    oligo_bases, overlap_bp
        The note's oligo shape, for a caller holding a supplier that states another.

    Returns
    -------
    tuple[StitchOligo, ...]
        Two or more, in the order they run along the molecule.

    Raises
    ------
    ValueError
        If `product` documents no single-stranded oligo in its reaction, if `bases` is no longer
        than one overlap, or if tiling it would take more than `STITCH_OLIGOS` oligos.

    Examples
    --------
    >>> from liulab_mbio.cloning.gibson.bench import NEBUILDER_HIFI
    >>> [one.strand.name for one in stitch_oligos("AT" * 75, product=NEBUILDER_HIFI)]
    ['FORWARD', 'REVERSE', 'FORWARD', 'REVERSE']
    """
    requires_oligos(product, "stitch a part out of overlapping oligos")
    if len(bases) <= overlap_bp:
        raise ValueError(
            f"{name or 'a stitched part'} is {len(bases)} bases, which is no longer than the "
            f"{overlap_bp} bp its oligos overlap by"
        )
    count = max(2, -(-(len(bases) - overlap_bp) // (oligo_bases - overlap_bp)))
    if count > STITCH_OLIGOS:
        limit = stitch_limit_bp(oligo_bases=oligo_bases, overlap_bp=overlap_bp)
        raise ValueError(
            f"{name or 'this part'} is {len(bases)} bases, too long to stitch: Gibson 2011 "
            f"assembles at most {STITCH_OLIGOS} {oligo_bases}-base oligos at one time, which "
            f"tile {limit} bases at a {overlap_bp} bp overlap. Amplify it instead"
        )
    spread = len(bases) + (count - 1) * overlap_bp
    lengths = [spread // count + (1 if index < spread % count else 0) for index in range(count)]
    oligos = []
    at = 0
    for index, length in enumerate(lengths):
        piece = bases[at : at + length]
        forward = index % 2 == 0
        oligos.append(
            StitchOligo(
                f"{name} oligo {index + 1}".strip(),
                piece if forward else reverse_complement(piece),
                Strand.FORWARD if forward else Strand.REVERSE,
                at,
            )
        )
        at += length - overlap_bp
    return tuple(oligos)


def stitch_checks(parts: Sequence[Part]) -> tuple[Check, ...]:
    """Judge every stitched part's size against the window Addgene puts this route between.

    Nothing is returned where no part is stitched. Addgene states the window as a
    recommendation, so a part outside it warns and is not refused: below it a primer tail
    carries the bases, above it a PCR or a synthesis order is cheaper (§14).
    """
    sized = [(part.name, part.fragment_length) for part in parts if part.stitched]
    if not sized:
        return ()
    low, high = STITCH_WINDOW_BP
    outside = [label for label, size in sized if not low <= size <= high]
    return (
        Check(
            "stitched size",
            "warn" if outside else "pass",
            len(sized),
            _said(
                sized,
                " bp",
                f"Addgene uses this route between {low} and {high} bp, above a primer tail and "
                "below a part of its own",
            ),
        ),
    )


def bridging_oligo(
    before: str,
    after: str,
    *,
    name: str = "",
    product: AssemblyProduct,
    homology_bp: int = BRIDGE_HOMOLOGY_BP,
) -> Bridge:
    """Design one oligo joining two fragments that share no homology, so neither needs a tail.

    It takes `homology_bp` bases from the end of `before` and as many from the start of `after`
    and spells them end to end, so the exonuclease leaves each fragment an end the oligo anneals
    to and the two are joined with nothing between them. Which strand it is written on does not
    matter: De Saeger measured both and found no difference (§15).

    Parameters
    ----------
    before, after
        The bases of the two fragments, in the order they go round the product.
    name
        What to order the oligo under.
    product
        The assembly product on the bench, which has to take single-stranded oligos at all.
    homology_bp
        Bases taken from each side.

    Raises
    ------
    ValueError
        If `product` documents no single-stranded oligo in its reaction, or if either fragment
        holds fewer bases than the homology asked for.

    Examples
    --------
    >>> from liulab_mbio.cloning.gibson.bench import NEBUILDER_HIFI
    >>> bridging_oligo("AAAACCCC", "GGGGTTTT", product=NEBUILDER_HIFI, homology_bp=4).sequence
    'CCCCGGGG'
    """
    requires_oligos(product, "join two fragments with a bridging oligo")
    for label, bases in (("before", before), ("after", after)):
        if len(bases) < homology_bp:
            raise ValueError(
                f"the fragment {label} this junction is {len(bases)} bases, too short for the "
                f"{homology_bp} bp of homology a bridging oligo takes from each side"
            )
    return Bridge(name or "bridge", before[-homology_bp:] + after[:homology_bp], homology_bp)


def requires_oligos(product: AssemblyProduct, route: str) -> None:
    """Refuse a route on a product whose supplier documents no single-stranded oligo in it.

    Both oligo routes need one. Where a supplier documents none the answer is the note's --
    the route is not supported -- rather than another supplier's numbers under this label.

    Raises
    ------
    ValueError
        Whenever it does not, naming the products that do.
    """
    if product.single_stranded_oligos:
        return
    supported = ", ".join(one.name for one in ASSEMBLY_PRODUCTS if one.single_stranded_oligos)
    raise ValueError(
        f"{product.supplier} documents no single-stranded oligo in the {product.name} reaction, "
        f"so it cannot {route}; {supported} can"
    )


def gc_percent(bases: str) -> float:
    """Return what share of `bases` is G or C, as a percentage.

    Examples
    --------
    >>> gc_percent("GGAT")
    50.0
    """
    upper = bases.upper()
    return 100.0 * sum(upper.count(base) for base in "GC") / len(upper) if upper else 0.0


def longest_palindrome(bases: str) -> int:
    """Return the longest stretch of `bases` that reads the same on both strands.

    Such a stretch pairs with itself rather than with the fragment it is meant to anneal to,
    which is the hairpin NEB asks a junction to be kept clear of.

    Examples
    --------
    >>> longest_palindrome("TTGAATTCAA")
    10
    >>> longest_palindrome("AAAAA")
    0
    """
    upper = bases.upper()
    longest = 0
    for centre in range(len(upper) - 1):
        arm = 0
        while (
            centre - arm >= 0
            and centre + 1 + arm < len(upper)
            and upper[centre - arm] == reverse_complement(upper[centre + 1 + arm])
        ):
            arm += 1
        longest = max(longest, 2 * arm)
    return longest


def end_repeat(bases: str) -> int:
    """Return the bases a tandem repeat covers at either end of `bases`, or 0 where none does.

    A repeat reaching an end is what NEB's His-tag rule refuses: two fragments whose ends spell
    the same repeating unit anneal to each other as readily as to the partner they are meant
    for. A repeat with `REPEAT_FLANK` bases of anything else outside it is left alone, which is
    the fix NEB gives.

    Examples
    --------
    >>> end_repeat("CACCACCACGTTA")
    9
    >>> end_repeat("GTCACCACCACTA")
    0
    """
    upper = bases.upper()
    shortest, longest = REPEAT_UNIT_BP
    covered = 0
    for unit_bp in range(shortest, longest + 1):
        for start in range(len(upper) - unit_bp + 1):
            unit = upper[start : start + unit_bp]
            copies = 1
            while upper[start + copies * unit_bp : start + (copies + 1) * unit_bp] == unit:
                copies += 1
            size = copies * unit_bp
            reaches = start < REPEAT_FLANK or len(upper) - (start + size) < REPEAT_FLANK
            if copies > 1 and size >= REPEAT_BASES and reaches:
                covered = max(covered, size)
    return covered


def most_alike(overlaps: Sequence[str]) -> tuple[int, int, int]:
    """Return the two overlaps with the longest stretch in common, and how long it is.

    Two indices and a count of bases, on either strand. Nothing says how alike is too alike --
    see the module docstring -- so this measures and does not judge.

    Raises
    ------
    ValueError
        If fewer than two overlaps are given.

    Examples
    --------
    >>> most_alike(["GGGGCATTAC", "TTTTTTTTTT", "AAAACATTAC"])
    (0, 2, 6)
    """
    if len(overlaps) < 2:
        raise ValueError(f"two overlaps are needed to compare, got {len(overlaps)}")
    pairs = (
        (first, second)
        for first in range(len(overlaps))
        for second in range(first + 1, len(overlaps))
    )
    return max(
        ((first, second, _in_common(overlaps[first], overlaps[second])) for first, second in pairs),
        key=lambda one: one[2],
    )


def overlap_checks(named: Sequence[tuple[str, str]], product: AssemblyProduct) -> tuple[Check, ...]:
    """Judge every overlap of one reaction, one check for each thing measured.

    `named` is each junction's name and the bases it spells, in the product's own order; there
    is one junction for each fragment, so the tier is read off how many there are. Each check
    covers every junction and names them all, as one verdict over many oligos does.
    """
    rule = product.tier(len(named)).overlap
    band = (
        f"{rule.shortest} bp"
        if rule.shortest == rule.longest
        else f"{rule.shortest}-{rule.longest} bp"
    )
    lengths = [(label, len(bases)) for label, bases in named]
    below = [label for label, size in lengths if size < product.shortest_overlap_bp]
    outside = [label for label, size in lengths if not rule.shortest <= size <= rule.longest]
    return (
        Check(
            "overlap length",
            "fail" if below else "warn" if outside else "pass",
            len(named),
            _said(
                lengths,
                " bp",
                f"{product.name} documents {band} here and calls "
                f"{product.shortest_overlap_bp} bp the shortest that works at all",
            ),
        ),
        _tm_check(named, rule),
        Check(
            "overlap gc",
            None,
            len(named),
            _said(
                [(label, round(gc_percent(bases))) for label, bases in named],
                "%",
                "no source quantifies a GC band for an overlap, so nothing judges this",
            ),
        ),
        _structure_check(
            "overlap hairpin",
            [(label, longest_palindrome(bases)) for label, bases in named],
            PALINDROME_BASES,
            "a highly palindromic overlap costs NEB up to tenfold fewer colonies; "
            f"{PALINDROME_BASES} bases reading the same on both strands is this package's floor",
        ),
        _structure_check(
            "overlap repeat",
            [(label, end_repeat(bases)) for label, bases in named],
            REPEAT_BASES,
            f"both NEB manuals refuse a repeating sequence at an overlap's end and ask for "
            f"{REPEAT_FLANK} bases of anything else outside it",
        ),
        _alike_check(named),
    )


def _tm_check(named: Sequence[tuple[str, str]], rule: OverlapRule) -> Check:
    """Judge every overlap's melting temperature, where the product states a floor for one."""
    read = [(label, round(wallace_tm(bases))) for label, bases in named]
    if rule.tm_floor is None:
        return _unjudged("overlap tm", read, " °C", "this product states no floor for one")
    cold = [label for label, tm in read if tm < rule.tm_floor]
    return Check(
        "overlap tm",
        "warn" if cold else "pass",
        len(named),
        _said(read, " °C", f"by the Wallace rule, against a floor of {rule.tm_floor:g} °C"),
    )


def _structure_check(name: str, measured: Sequence[tuple[str, int]], floor: int, why: str) -> Check:
    """Judge a structure that costs an assembly colonies, naming what each junction carries."""
    over = [label for label, size in measured if size >= floor]
    return Check(
        name,
        "warn" if over else "pass",
        len(measured),
        _said(measured, " bp", why),
    )


def _alike_check(named: Sequence[tuple[str, str]]) -> Check:
    """Report the two overlaps of one reaction with the most in common, and judge neither."""
    if len(named) < 2:
        return _unjudged("overlap similarity", [], "", "one junction has nothing to resemble")
    first, second, shared = most_alike([bases for _, bases in named])
    return Check(
        "overlap similarity",
        None,
        shared,
        f"{named[first][0]} and {named[second][0]} share {shared} bases, the most any two here "
        "do; NEB asks for a unique overlap at each junction but nobody quantifies how alike is "
        "too alike, so nothing judges this",
    )


def _unjudged(name: str, measured: Sequence[tuple[str, int]], unit: str, why: str) -> Check:
    """Return a check carrying no verdict, saying where it is shown that nothing judged it."""
    return Check(name, None, len(measured), _said(measured, unit, why))


def _said(measured: Sequence[tuple[str, float]], unit: str, why: str) -> str:
    """Return what each junction measured, then the band or the reason there is none."""
    read = ", ".join(f"{label} {value:g}{unit}" for label, value in measured)
    return f"{read}; {why}" if read else why


def _in_common(one: str, other: str) -> int:
    """Return the longest stretch two overlaps both spell, on either strand."""
    return max(_longest_shared(one, reading) for reading in (other, reverse_complement(other)))


def _longest_shared(one: str, other: str) -> int:
    """Return the longest substring of `one` that `other` also spells."""
    row = [0] * (len(other) + 1)
    longest = 0
    for base in one:
        carried = 0
        for index, against in enumerate(other, start=1):
            carried, row[index] = row[index], carried + 1 if base == against else 0
            longest = max(longest, row[index])
    return longest

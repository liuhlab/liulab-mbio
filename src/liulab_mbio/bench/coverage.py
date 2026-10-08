"""How many constructs a library design yields, and how many colonies each round needs.

Completeness is the chance that no product of a round is missing from its colonies. Each round
is sized for it rather than the library sized once at the end: a round short of its floor loses
members no later round can put back, and every later round multiplies what survives. Coverage,
colonies over products, is what the floor works out at.

`plan_coverage` is the way in. The completeness rule is Clarke & Carbon's, in `REFERENCES`, and
it holds only where every member of a part list is equally represented.
"""

import math
from collections.abc import Sequence
from dataclasses import KW_ONLY, dataclass

from liulab_mbio.protocol.model import Reference


def _positive(products: int) -> None:
    """Refuse a round with nothing to cover."""
    if products < 1:
        raise ValueError(f"products must be positive, got {products}")


def constructs(part_list_sizes: Sequence[int]) -> int:
    """Return how many distinct constructs part lists of these sizes yield.

    One round appends one part list to every member of the library, so the count is the product
    of the sizes.

    Raises
    ------
    ValueError
        If no part list is given, or one of them is empty.

    Examples
    --------
    >>> constructs((24, 24, 24))
    13824
    """
    if not part_list_sizes:
        raise ValueError("a library needs at least one part list")
    for size in part_list_sizes:
        if size < 1:
            raise ValueError(f"a part list holds at least one part, got {size}")
    return math.prod(part_list_sizes)


def colonies_for_coverage(products: int, coverage: float) -> int:
    """Return the colonies a round needs to cover `products` that many times over.

    Coverage is colonies over distinct products, so this is that definition read backwards and
    rounded up to a whole colony.

    Raises
    ------
    ValueError
        If `products` or `coverage` is not positive.

    Examples
    --------
    >>> colonies_for_coverage(576, 10)
    5760
    """
    _positive(products)
    if coverage <= 0:
        raise ValueError(f"coverage must be positive, got {coverage}")
    return math.ceil(products * coverage)


def absent_probability(products: int, colonies: int) -> float:
    """Return the chance one named product is missing from this many colonies.

    Each colony is one draw from `products` equally represented products, so the chance of
    missing a given one is ``(1 - 1 / products) ** colonies``.

    Raises
    ------
    ValueError
        If `products` is not positive, or `colonies` is negative.

    Examples
    --------
    >>> round(absent_probability(100, 500), 4)
    0.0066
    """
    _positive(products)
    if colonies < 0:
        raise ValueError(f"colonies cannot be negative, got {colonies}")
    return (1.0 - 1.0 / products) ** colonies


def colonies_for_completeness(products: int, probability: float) -> int:
    """Return the colonies a round needs for `probability` that no product is missing.

    Clarke & Carbon's ``N = ln(1 - P) / ln(1 - f)``, with ``f = 1 / products`` for a library of
    equally represented members, asked of every member at once by giving each
    ``probability ** (1 / products)``.

    Equal representation is the assumption, and a real part list breaks it: where synthesis or
    an earlier round under-delivers a member, this is a floor rather than an answer.

    Raises
    ------
    ValueError
        If `products` is not positive, or `probability` does not lie between 0 and 1.

    Examples
    --------
    >>> colonies_for_completeness(100, 0.95)
    754
    """
    _positive(products)
    if not 0.0 < probability < 1.0:
        raise ValueError(f"probability must lie between 0 and 1, got {probability}")
    if products == 1:
        return 1
    each = probability ** (1.0 / products)
    return math.ceil(math.log(1.0 - each) / math.log(1.0 - 1.0 / products))


@dataclass(frozen=True, slots=True)
class RoundCoverage:
    """What one round has to cover, and what the colonies asked for leave out.

    Parameters
    ----------
    number
        Which round it is, counting from one.
    part_list_size
        How many parts that round's part list holds.
    products
        The distinct products the round can make: every part list up to and including its own.
    completeness
        The chance the round is sized for that none of its products is missing.
    colonies
        The floor that completeness takes.
    coverage
        What that floor works out at, colonies over products. It is derived, not asked for.
    absent_probability
        The chance one named product is missing from that many colonies.
    """

    number: int
    part_list_size: int
    _: KW_ONLY
    products: int
    completeness: float
    colonies: int
    coverage: float
    absent_probability: float


def plan_coverage(
    part_list_sizes: Sequence[int], *, completeness: float
) -> tuple[RoundCoverage, ...]:
    """Return what each round has to cover, one row per round, in the order they run.

    Round *n* joins the *n*-th part list to the library the rounds before it built, so its
    products are the sizes up to and including its own, multiplied. `completeness` has no
    default: how much of a library a round may lose is the user's call, and no source sets it.

    Raises
    ------
    ValueError
        If no part list is given, one of them is empty, or `completeness` does not lie between
        0 and 1.

    Examples
    --------
    >>> [row.colonies for row in plan_coverage((4, 6), completeness=0.99)]
    [21, 183]
    """
    constructs(part_list_sizes)  # refuses the whole list before a single row is built
    rows: list[RoundCoverage] = []
    for number, size in enumerate(part_list_sizes, start=1):
        products = constructs(part_list_sizes[:number])
        colonies = colonies_for_completeness(products, completeness)
        rows.append(
            RoundCoverage(
                number,
                size,
                products=products,
                completeness=completeness,
                colonies=colonies,
                coverage=colonies / products,
                absent_probability=absent_probability(products, colonies),
            )
        )
    return tuple(rows)


def skew_ratio(counts: Sequence[float]) -> float:
    """Return a count distribution's skew: its 90th percentile over its 10th.

    Joung's measure of how evenly a pooled library is represented, read off the counts one
    member at a time. A perfectly even library answers 1.

    Raises
    ------
    ValueError
        If no count is given, one of them is negative, or the 10th percentile is zero — a
        library a tenth of whose members went unread carries no ratio, only a share seen.

    Examples
    --------
    >>> skew_ratio([5, 5, 5, 5])
    1.0
    >>> skew_ratio([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    5.0
    """
    if not counts:
        raise ValueError("a skew ratio is measured over at least one count")
    for one in counts:
        if one < 0:
            raise ValueError(f"a member cannot be read {one} times")
    ranked = sorted(float(one) for one in counts)
    low = _percentile(ranked, 0.10)
    if low == 0.0:
        raise ValueError(
            "the 10th percentile of these counts is zero, so they carry no skew ratio: report "
            "the share of members seen instead"
        )
    return _percentile(ranked, 0.90) / low


def _percentile(ranked: Sequence[float], fraction: float) -> float:
    """Return a percentile of already sorted counts, interpolating between neighbours."""
    if len(ranked) == 1:
        return ranked[0]
    at = fraction * (len(ranked) - 1)
    below = math.floor(at)
    above = math.ceil(at)
    return ranked[below] + (ranked[above] - ranked[below]) * (at - below)


def reads_for_representation(constructs: int, per_member: int = 100) -> int:
    """Return the reads a representation read takes: one member's share, over every member.

    Joung judges a pooled library at over 100 reads a member, so the depth follows from the
    library's width and is computed rather than stated.

    Raises
    ------
    ValueError
        If either number is not positive.

    Examples
    --------
    >>> reads_for_representation(13824)
    1382400
    """
    if constructs < 1:
        raise ValueError(f"a library holds at least one construct, got {constructs}")
    if per_member < 1:
        raise ValueError(f"a member is read at least once, got {per_member}")
    return constructs * per_member


@dataclass(frozen=True, slots=True)
class RepresentationMarks:
    """What a pooled library's representation read has to reach.

    Every one is Joung's, measured on an sgRNA library counted by a barcode amplicon before a
    screen. They transfer because the iGGA representation read amplifies the barcode block
    alone, so every member's counted amplicon is the same length. A project may tighten each
    and may not loosen it.

    Parameters
    ----------
    seen
        The share of the library's combinations that have to be read at all.
    skew
        The 90th/10th percentile ratio the counts have to stay under.
    reads_per_member
        The depth the other two are judged at.
    """

    seen: float
    skew: float
    reads_per_member: int


#: Joung's acceptance bar for a pooled library, quoted in `docs/research/vector-qc-panel.md`
#: section 3.5: fewer than 0.5% of members undetected, a skew ratio under 10, judged at over 100
#: reads a member.
REPRESENTATION_MARKS = RepresentationMarks(seen=0.995, skew=10.0, reads_per_member=100)


#: Where the completeness rule and the representation marks come from, ready for a protocol's
#: reference list.
REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Clarke, L. and Carbon, J. (1976) A colony bank containing synthetic ColE1 hybrid "
        "plasmids representative of the entire E. coli genome. Cell 9, 91-99, for the clones a "
        "library needs to represent every member",
        url="https://doi.org/10.1016/0092-8674(76)90055-6",
    ),
    Reference(
        "Joung, J. et al. (2017) Genome-scale CRISPR-Cas9 knockout and transcriptional "
        "activation screening. Nat. Protoc. 12, 828-863, for a pooled library's acceptance "
        "bar: under 0.5% of members undetected, a 90th/10th percentile skew ratio under 10, "
        "judged at over 100 reads a member",
        url="https://doi.org/10.1038/nprot.2017.016",
    ),
    Reference(
        "Imkeller, K. et al. (2020) gscreend: modelling asymmetric count ratios in CRISPR "
        "screens to decrease experiment size and improve phenotype detection. Genome Biol. 21, "
        "53, Table 2, for the screen coverage a library's own skew demands: p90/p10 2.5 at "
        "200x, 5 at 300x, 10 at 400x",
        url="https://doi.org/10.1186/s13059-020-1939-1",
    ),
    Reference(
        "Heo, S.-J. et al. (2024) Compact CRISPR genetic screens enabled by improved guide RNA "
        "library cloning. Genome Biol. 25, 25, for a library skewed under 2 matching a "
        "1,000-fold screen at 100-fold coverage",
        url="https://doi.org/10.1186/s13059-023-03132-3",
    ),
)

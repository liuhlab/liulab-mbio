"""How many constructs a library design yields, and how many colonies each round needs.

Library coverage is one round's colonies counted against the distinct products that round could
make. It is counted per round rather than once at the end: a round short of the coverage asked
for loses members no later round can put back, and every later round multiplies what survives.

`plan_coverage` is the way in. The completeness rule is Clarke & Carbon's, in `REFERENCES`, and
it holds only where every member of a part list is equally represented.
"""

import math
from collections.abc import Sequence
from dataclasses import KW_ONLY, dataclass

from liulab_mbio.protocol import Reference


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
    coverage
        The library coverage asked for, colonies over products.
    colonies
        What that coverage takes.
    absent_probability
        The chance one named product is missing from that many colonies.
    """

    number: int
    part_list_size: int
    _: KW_ONLY
    products: int
    coverage: float
    colonies: int
    absent_probability: float


def plan_coverage(part_list_sizes: Sequence[int], *, coverage: float) -> tuple[RoundCoverage, ...]:
    """Return what each round has to cover, one row per round, in the order they run.

    Round *n* joins the *n*-th part list to the library the rounds before it built, so its
    products are the sizes up to and including its own, multiplied. `coverage` has no default:
    how much of a library a round may lose is the user's call, and no source sets it.

    Raises
    ------
    ValueError
        If no part list is given, one of them is empty, or `coverage` is not positive.

    Examples
    --------
    >>> [row.colonies for row in plan_coverage((4, 6), coverage=10)]
    [40, 240]
    """
    constructs(part_list_sizes)  # refuses the whole list before a single row is built
    rows: list[RoundCoverage] = []
    for number, size in enumerate(part_list_sizes, start=1):
        products = constructs(part_list_sizes[:number])
        colonies = colonies_for_coverage(products, coverage)
        rows.append(
            RoundCoverage(
                number,
                size,
                products=products,
                coverage=coverage,
                colonies=colonies,
                absent_probability=absent_probability(products, colonies),
            )
        )
    return tuple(rows)


#: Where the completeness rule comes from, ready for a protocol's reference list.
REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Clarke, L. and Carbon, J. (1976) A colony bank containing synthetic ColE1 hybrid "
        "plasmids representative of the entire E. coli genome. Cell 9, 91-99, for the clones a "
        "library needs to represent every member",
        url="https://doi.org/10.1016/0092-8674(76)90055-6",
    ),
)

"""What a Gibson assembly takes at the bench: its product, its reaction and its incubation.

An **assembly product** is the kit on the bench, and it is one row of data: the overlap band it
documents, how long it is incubated, how much DNA it takes and at what molar ratio, and how many
inserts it is documented for. Those vary by fragment count, so a product carries one `Tier` per
band its manual gives and `AssemblyProduct.tier` picks the one a reaction falls in. A fourth
product is another row here, never another module.

Every number is the supplier's, through ``docs/research/gibson-assembly.md``; the section each
comes from is named beside it. What any cloning pipeline shares -- DNA amounts, PCR, gels and
validation -- is `liulab_mbio.bench`.
"""

from collections.abc import Sequence
from dataclasses import dataclass

from liulab_mbio.bench.amounts import Amount, dna_amount, to_pmol
from liulab_mbio.bench.reactions import reaction_table
from liulab_mbio.protocol.model import (
    Component,
    Incubation,
    ReactionTable,
    Reference,
    Stage,
    ThermocyclerProgram,
)


@dataclass(frozen=True, slots=True)
class OverlapRule:
    """The overlap one assembly product documents for one fragment count.

    Parameters
    ----------
    shortest, longest
        The band the supplier asks for, in base pairs.
    tm_floor
        The melting temperature the overlap has to reach, °C by the Wallace rule, or ``None``
        where the supplier states none. A design takes the shortest overlap in the band that
        reaches it.
    """

    shortest: int
    longest: int
    tm_floor: float | None


@dataclass(frozen=True, slots=True)
class Tier:
    """What an assembly product documents for reactions of up to `fragments` fragments.

    The supplier's bands move together: a reaction of more fragments takes a longer overlap, a
    longer incubation, more DNA and a different molar ratio.

    Parameters
    ----------
    fragments
        The most fragments this tier covers, the vector counted among them.
    overlap
        The overlap rule for it.
    incubation_seconds
        How long the reaction is held at the product's temperature.
    total_pmol
        The total picomoles of DNA the supplier's table takes, low and high.
    insert_ratio
        Moles of each insert for each mole of vector.
    """

    fragments: int
    overlap: OverlapRule
    incubation_seconds: int
    total_pmol: tuple[float, float]
    insert_ratio: float


@dataclass(frozen=True, slots=True)
class AssemblyProduct:
    """One assembly product a lab holds, and everything its manual documents about it.

    Parameters
    ----------
    name, supplier, catalog
        What to order, and from whom.
    master_mix_fold
        The master mix's concentration, as its label reads.
    master_mix_ul, reaction_ul
        Microlitres of master mix, and what the reaction is made up to.
    celsius
        The incubation temperature.
    tiers
        One per band the manual gives, fewest fragments first.
    shortest_overlap_bp
        The shortest overlap the manual documents at all, below which a design refuses.
    inserts_limit
        The most inserts the supplier recommends in one reaction.
    unpurified_fraction
        How much of the reaction unpurified PCR product may be.
    references
        Where the numbers above come from.
    """

    name: str
    supplier: str
    catalog: str
    master_mix_fold: str
    master_mix_ul: float
    reaction_ul: float
    celsius: float
    tiers: tuple[Tier, ...]
    shortest_overlap_bp: int
    inserts_limit: int
    unpurified_fraction: float
    references: tuple[Reference, ...]

    def tier(self, fragments: int) -> Tier:
        """Return what this product documents for a reaction of `fragments` fragments.

        A reaction of more fragments than any tier covers takes the last one, which is the most
        the manual says.

        Raises
        ------
        ValueError
            If fewer than two fragments are given.
        """
        if fragments < 2:
            raise ValueError(f"an assembly joins at least two fragments, got {fragments}")
        return next((tier for tier in self.tiers if fragments <= tier.fragments), self.tiers[-1])

    @property
    def unpurified_ul(self) -> float:
        """The most unpurified PCR product one reaction takes, µL."""
        return round(self.reaction_ul * self.unpurified_fraction, 2)


#: The NEBuilder HiFi manual and NEB's own posting of its reaction, which is the one document
#: here whose licence lets its table be reproduced. Note §1 records the verdict for each.
NEBUILDER_REFERENCES: tuple[Reference, ...] = (
    Reference(
        "NEB, NEBuilder HiFi DNA Assembly Master Mix / Cloning Kit instruction manual, "
        "NEB #E2621S/L/X and #E5520S, version 6.0_1/26",
        url="https://www.neb.com/-/media/nebus/files/manuals/manuale2621_e5520.pdf",
    ),
    Reference(
        "New England Biolabs (2022) NEBuilder HiFi DNA Assembly Reaction (E2621), protocols.io, "
        "under CC BY",
        url="https://dx.doi.org/10.17504/protocols.io.bfhrjj56",
    ),
)

#: NEBuilder HiFi, the product a plan assumes. Its two tiers are the two columns of the manual's
#: reaction table: overlap bands from note §3, incubations from §7, totals and ratios from §8,
#: the 48 °C Wallace floor from §4, the five-insert ceiling from §9, and the 20% cap on
#: unpurified PCR product from §11.
NEBUILDER_HIFI = AssemblyProduct(
    "NEBuilder HiFi DNA Assembly Master Mix",
    "New England Biolabs",
    "E2621",
    "2X",
    10.0,
    20.0,
    50.0,
    (
        Tier(3, OverlapRule(15, 20, 48.0), 900, (0.03, 0.2), 2.0),
        Tier(6, OverlapRule(20, 30, 48.0), 3600, (0.2, 0.5), 1.0),
    ),
    12,
    5,
    0.2,
    NEBUILDER_REFERENCES,
)

#: Nanograms of linearised vector NEB's own footnote calls optimal, the floor of its 50-100 ng
#: band (§8). Each insert is then `Tier.insert_ratio` times its moles.
VECTOR_NG = 50.0

#: NEB asks for a fivefold molar excess of any insert shorter than this, whatever the tier (§8).
SHORT_INSERT_BP = 200
SHORT_INSERT_RATIO = 5.0


def assembly_amounts(
    vector: tuple[str, int],
    inserts: Sequence[tuple[str, int]],
    *,
    tier: Tier,
    vector_ng: float = VECTOR_NG,
) -> tuple[Amount, ...]:
    """Return what to put in the assembly, vector first.

    The vector and each insert are a name and a length in base pairs. The vector goes in at
    `vector_ng`, and each insert at `Tier.insert_ratio` times its moles -- or at
    `SHORT_INSERT_RATIO` where the insert is shorter than `SHORT_INSERT_BP`, which is NEB's own
    exception.

    Raises
    ------
    ValueError
        If a length is not positive.

    Examples
    --------
    >>> [one.pmol for one in assembly_amounts(("pUC19", 2629), [("GFP", 717)], tier=NEBUILDER_HIFI.tiers[0])]
    [0.0309, 0.0618]
    """
    name, length_bp = vector
    if length_bp <= 0:
        raise ValueError(f"amount {name!r}: length_bp must be positive")
    pmol = round(to_pmol(vector_ng, length_bp), 4)
    return (
        dna_amount(name, length_bp, pmol=pmol),
        *(
            dna_amount(
                insert,
                bases,
                pmol=round(pmol * _ratio(bases, tier), 4),
            )
            for insert, bases in inserts
        ),
    )


def _ratio(length_bp: int, tier: Tier) -> float:
    """Return the molar excess one insert of this length goes in at."""
    return SHORT_INSERT_RATIO if length_bp < SHORT_INSERT_BP else tier.insert_ratio


def assembly_reaction(
    product: AssemblyProduct, amounts: Sequence[Amount], *, reactions: int = 1
) -> ReactionTable:
    """Return the one-tube assembly reaction this product's manual sets up.

    The DNA goes in each tube and the master mix is half the reaction, which is the order and
    the proportion NEB's table gives.

    Raises
    ------
    ValueError
        If fewer than two fragments are given, or the DNA does not fit the reaction volume.
    """
    if len(amounts) < 2:
        raise ValueError("an assembly joins at least two fragments")
    return reaction_table(
        amounts,
        (
            Component(
                product.name,
                product.master_mix_ul,
                stock=product.master_mix_fold,
                final="1X",
            ),
        ),
        volume_ul=product.reaction_ul,
        title=f"{product.name} reaction",
        reactions=reactions,
    )


def assembly_program(product: AssemblyProduct, *, fragments: int) -> ThermocyclerProgram:
    """Return the incubation this product asks for, at the fragment count's own tier.

    `fragments` counts the vector. There is no cycling: the reaction is isothermal, and it goes
    on ice or to -20 °C afterwards.

    Raises
    ------
    ValueError
        If fewer than two fragments are given.

    Examples
    --------
    >>> assembly_program(NEBUILDER_HIFI, fragments=2).stages[0].incubations[0].seconds
    900
    """
    tier = product.tier(fragments)
    return ThermocyclerProgram(
        (
            Stage((Incubation("Assembly", product.celsius, tier.incubation_seconds),)),
            Stage((Incubation("Hold", 4.0, None),)),
        ),
        title="Assembly incubation",
    )

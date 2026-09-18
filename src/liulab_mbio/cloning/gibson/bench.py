"""What a Gibson assembly takes at the bench: its product, its reaction and its incubation.

An **assembly product** is the kit on the bench, and it is one row of data: the overlap band it
documents, how long it is incubated, how much DNA it takes and at what molar ratio, and how many
inserts it is documented for. Those vary by fragment count, so a product carries one `Tier` per
band its manual gives and `AssemblyProduct.tier` picks the one a reaction falls in. A fourth
product is another row here, never another module.

The three shipped are `ASSEMBLY_PRODUCTS`. Two are NEB's and share a mechanism; In-Fusion's is
its own, so it inherits none of their numbers and states none of its own where the note records
that Takara publishes none. A field is ``None`` exactly there, and the check reading it then
carries no verdict.

Every number is the supplier's, through ``docs/research/gibson-assembly.md``; the section each
comes from is named beside it. What any cloning pipeline shares -- DNA amounts, PCR, gels and
validation -- is `liulab_mbio.bench`.
"""

from collections.abc import Sequence
from dataclasses import dataclass

from liulab_mbio.bench.amounts import Amount, dna_amount, to_pmol
from liulab_mbio.bench.reactions import reaction_table
from liulab_mbio.checks import Check
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
        The total picomoles of DNA the supplier's table takes, low and high, or ``None`` where
        the supplier states a mass instead and no picomole band is sourced.
    insert_ratio
        Moles of each insert for each mole of vector.
    """

    fragments: int
    overlap: OverlapRule
    incubation_seconds: int
    total_pmol: tuple[float, float] | None
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
    incubation_note
        What the supplier says about incubating longer than its table asks for.
    tiers
        One per band the manual gives, fewest fragments first.
    shortest_overlap_bp
        The shortest overlap the manual documents at all, below which a design refuses.
    inserts_limit
        The most inserts the supplier recommends in one reaction.
    short_insert_bp, short_insert_ratio
        An insert shorter than `short_insert_bp` goes in at `short_insert_ratio` times the
        vector's moles instead of its tier's, which is the supplier's own exception.
    unpurified_fraction
        How much of the reaction unpurified PCR product may be, or ``None`` where the supplier
        asks for purified DNA and documents no such allowance.
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
    incubation_note: str
    tiers: tuple[Tier, ...]
    shortest_overlap_bp: int
    inserts_limit: int
    short_insert_bp: int
    short_insert_ratio: float
    unpurified_fraction: float | None
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
    def unpurified_ul(self) -> float | None:
        """The most unpurified PCR product one reaction takes, µL, where the supplier allows it."""
        if self.unpurified_fraction is None:
            return None
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
    "Incubating longer than 15 minutes improves some assemblies, up to four hours; overnight "
    "is refused.",
    (
        Tier(3, OverlapRule(15, 20, 48.0), 900, (0.03, 0.2), 2.0),
        Tier(6, OverlapRule(20, 30, 48.0), 3600, (0.2, 0.5), 1.0),
    ),
    12,
    5,
    200,
    5.0,
    0.2,
    NEBUILDER_REFERENCES,
)

#: The Gibson Assembly Master Mix manual, and NEB's own CC BY posting of its reaction (§1).
GIBSON_REFERENCES: tuple[Reference, ...] = (
    Reference(
        "NEB, Gibson Assembly Master Mix / Gibson Assembly Cloning Kit instruction manual, "
        "NEB #E2611S/L and #E5510S, version 3.0_1/26",
        url="https://www.neb.com/-/media/nebus/files/manuals/manuale2611_e5510.pdf",
    ),
    Reference(
        "New England Biolabs (2022) Gibson Assembly Master Mix - Assembly (E2611), "
        "protocols.io, under CC BY",
        url="https://www.protocols.io/view/gibson-assembly-master-mix-assembly-e2611-bdd8i29w",
    ),
)

#: The Gibson Assembly Master Mix, which takes more insert than NEBuilder at low fragment counts
#: and a longer overlap at high ones: bands from §3, incubations from §7, totals and ratios from
#: §8. Its 40 bp ceiling is FAQ 6's, the figure tied to the exonuclease quantity; the same
#: manual's FAQ 23 says 80 nt, and §3 says not to average the two.
GIBSON_MASTER_MIX = AssemblyProduct(
    "Gibson Assembly Master Mix",
    "New England Biolabs",
    "E2611",
    "2X",
    10.0,
    20.0,
    50.0,
    "Incubating longer than 15 minutes improves some assemblies, up to four hours; overnight "
    "is refused.",
    (
        Tier(3, OverlapRule(15, 25, 48.0), 900, (0.02, 0.5), 2.0),
        Tier(6, OverlapRule(20, 40, 48.0), 3600, (0.2, 1.0), 1.0),
    ),
    12,
    5,
    200,
    5.0,
    0.2,
    GIBSON_REFERENCES,
)

#: Takara's manuals for In-Fusion, whose rules are its own throughout (§1, §2).
IN_FUSION_REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Takara Bio, In-Fusion Snap Assembly User Manual (060822)",
        url="https://www.takarabio.com/documents/User%20Manual/In/In-Fusion%20Snap%20Assembly%20User%20Manual.pdf",
    ),
    Reference(
        "Takara Bio, In-Fusion Cloning FAQs, retrieved 2026-09-18",
        url="https://www.takarabio.com/learning-centers/cloning/in-fusion-cloning-faqs",
    ),
)

#: In-Fusion, whose mechanism is not NEB's and which inherits none of NEB's numbers (§2). Takara
#: states one overlap length for one insert and another above two fragments rather than a band
#: (§3), no melting temperature for it at all (§4), a 10 µL reaction of which 2 µL is a 5X mix
#: (§6), one incubation whatever the fragment count (§7), a 2:1 molar ratio with its own
#: small-fragment exception (§8), five inserts (§9) and purified PCR product only (§6).
IN_FUSION = AssemblyProduct(
    "In-Fusion Snap Assembly Master Mix",
    "Takara Bio",
    "638947",
    "5X",
    2.0,
    10.0,
    50.0,
    "The reaction finishes inside its 15 minutes; incubating longer does not help and can "
    "leave ends that anneal worse.",
    (
        Tier(2, OverlapRule(15, 15, None), 900, None, 2.0),
        Tier(6, OverlapRule(20, 20, None), 900, None, 2.0),
    ),
    12,
    5,
    350,
    3.0,
    None,
    IN_FUSION_REFERENCES,
)

#: Every assembly product a plan can be run with, the default first.
ASSEMBLY_PRODUCTS: tuple[AssemblyProduct, ...] = (NEBUILDER_HIFI, GIBSON_MASTER_MIX, IN_FUSION)

#: Nanograms of linearised vector both suppliers call optimal, the floor of NEB's 50-100 ng band
#: and of Takara's own (§8). Each insert is then its tier's ratio times the vector's moles.
VECTOR_NG = 50.0


def assembly_product(name: str) -> AssemblyProduct:
    """Return the assembly product of that name, or the one it is the beginning of.

    Raises
    ------
    ValueError
        If no shipped product answers to it, or more than one does.

    Examples
    --------
    >>> assembly_product("in-fusion").supplier
    'Takara Bio'
    """
    wanted = name.strip().lower()
    found = [one for one in ASSEMBLY_PRODUCTS if one.name.lower().startswith(wanted)]
    if len(found) != 1:
        shipped = ", ".join(one.name for one in ASSEMBLY_PRODUCTS)
        said = "no assembly product" if not found else "more than one assembly product"
        raise ValueError(f"{said} called {name!r}; this package ships {shipped}")
    return found[0]


def assembly_amounts(
    vector: tuple[str, int],
    inserts: Sequence[tuple[str, int]],
    *,
    product: AssemblyProduct,
    vector_ng: float = VECTOR_NG,
) -> tuple[Amount, ...]:
    """Return what to put in the assembly, vector first.

    The vector and each insert are a name and a length in base pairs. The vector goes in at
    `vector_ng`, and each insert at the molar ratio `product` documents for this many
    fragments -- or at `AssemblyProduct.short_insert_ratio` where the insert is shorter than
    the product's own small-fragment exception.

    Raises
    ------
    ValueError
        If a length is not positive, or fewer than two fragments are given.

    Examples
    --------
    >>> [one.pmol for one in assembly_amounts(("pUC19", 2629), [("GFP", 717)], product=NEBUILDER_HIFI)]
    [0.0309, 0.0618]
    """
    name, length_bp = vector
    if length_bp <= 0:
        raise ValueError(f"amount {name!r}: length_bp must be positive")
    tier = product.tier(len(inserts) + 1)
    pmol = round(to_pmol(vector_ng, length_bp), 4)
    return (
        dna_amount(name, length_bp, pmol=pmol),
        *(
            dna_amount(
                insert,
                bases,
                pmol=round(pmol * _ratio(bases, product, tier), 4),
            )
            for insert, bases in inserts
        ),
    )


def _ratio(length_bp: int, product: AssemblyProduct, tier: Tier) -> float:
    """Return the molar excess one insert of this length goes in at."""
    return product.short_insert_ratio if length_bp < product.short_insert_bp else tier.insert_ratio


def fragment_check(product: AssemblyProduct, inserts: int) -> Check:
    """Judge how many inserts one reaction joins against what `product` is documented for.

    Above the supplier's own ceiling the assembly still runs; what falls is the share of clones
    that carry the right one, so this warns rather than refuses.

    Examples
    --------
    >>> fragment_check(NEBUILDER_HIFI, 6).status
    'warn'
    """
    limit = product.inserts_limit
    within = inserts <= limit
    return Check(
        "fragment count",
        "pass" if within else "warn",
        inserts,
        f"{inserts} insert{'' if inserts == 1 else 's'} in one reaction; {product.supplier} "
        f"recommends {limit} or fewer for {product.name}"
        + ("" if within else ", and assembling the rest in a second round"),
    )


def assembly_dna_check(product: AssemblyProduct, amounts: Sequence[Amount]) -> Check:
    """Judge the picomoles in the reaction against the band `product` documents for them.

    A supplier stating a mass and no picomole band judges nothing, so the check carries no
    verdict and says which supplier states what.
    """
    total = sum(amount.pmol for amount in amounts)
    band = product.tier(len(amounts)).total_pmol
    if band is None:
        return Check(
            "assembly DNA",
            None,
            round(total, 4),
            f"{total:g} pmol in {product.reaction_ul:g} µL; {product.supplier} states a mass "
            f"for {product.name} and no picomole band, so nothing judges this",
        )
    low, high = band
    return Check(
        "assembly DNA",
        "pass" if low <= total <= high else "warn",
        round(total, 4),
        f"{total:g} pmol in {product.reaction_ul:g} µL; {product.name} takes {low:g}-{high:g} "
        "pmol at this fragment count",
    )


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

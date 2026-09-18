"""The Gateway reactions: what each takes, how long it runs, and how it is stopped.

What any bench shares is `liulab_mbio.bench`. Here is only what no shared table covers: the LR
reaction's own components, its incubation, the proteinase K that ends it, and the
transformation Invitrogen's manual prescribes, which is not the one NEB's does. Every number
cites `docs/research/gateway-cloning.md`.
"""

from collections.abc import Sequence

from liulab_mbio.bench.amounts import DNA_VOLUME_UL, Amount, to_pmol
from liulab_mbio.bench.reactions import dna_components, fits
from liulab_mbio.bench.steps import Transformation
from liulab_mbio.protocol.model import Component, ReactionTable, Reference

#: The LR reaction with LR Clonase II, from MAN0001032 revision A.0 page 3 and MAN0000470
#: page 32 (note §11): microlitres, then what the DNA is made up to before the enzyme goes in.
LR_VOLUME_UL = 10.0
LR_TE_UL = 8.0
LR_CLONASE_UL = 2.0

#: What the reaction takes of each plasmid, ng. The manual's band for the entry clone is
#: 50-150 ng, and it names a symptom for going over: colonies carrying several molecules, with
#: a small-colony phenotype. Highest yields come from 150 ng of each (MAN0001032 page 2, §11).
ENTRY_NG = 150.0
ENTRY_MIN_NG = 50.0
DESTINATION_NG = 150.0

#: The incubation, and how far it may be extended for a plasmid the manual calls large
#: (MAN0000470 page 32; note §11): degrees Celsius, seconds and base pairs.
LR_CELSIUS = 25.0
LR_SECONDS = 3600
LR_LONG_SECONDS = 64800
LR_LONG_BP = 10000

#: What ends it: proteinase K, supplied with the enzyme mix (MAN0001032 page 3 and MAN0000470
#: page 32; note §11, §12). Microlitres, micrograms per microlitre, Celsius and seconds.
PROTEINASE_K_UL = 1.0
PROTEINASE_K_UG_UL = 2.0
STOP_CELSIUS = 37.0
STOP_SECONDS = 600

#: What the vendor calls an efficient LR reaction, with the whole of it transformed and plated,
#: and the competent cells that count assumes, cfu/µg (MAN0001032 page 3 and MAN0000470
#: page 32; note §14). What fraction of them is correct is not published: see the note's gaps.
LR_COLONIES = 5000
CELL_EFFICIENCY_CFU_UG = 100_000_000

#: What the buffer is called on the bench, and what the enzyme mix is ordered as. The catalogue
#: number is the 20-reaction pack (note §17).
TE_BUFFER = "TE buffer, pH 8.0"
LR_CLONASE = "Gateway LR Clonase II enzyme mix"
LR_CLONASE_CATALOG = "11791020"
SUPPLIER = "Invitrogen / Thermo Fisher"

#: The strain to select the clone in. A general cloning strain is what the vendor's own support
#: page names for plating an LR reaction, and CcdB kills it, so unreacted vector and by-product
#: do not grow. Never one carrying F', whose ccdA cancels that (MAN0000470 pages 24 and 31;
#: note §9).
DEFAULT_HOST = "TOP10 chemically competent E. coli"

#: Invitrogen's transformation, which is not NEB's (MAN0000470 page 26; note §14, §17). No
#: manual in the set states a thaw time, so the step prints none rather than another kit's.
GATEWAY_TRANSFORMATION = Transformation(
    cells_ul=50.0,
    reaction_ul=1.0,
    source="the LR reaction",
    thaw_seconds=None,
    ice_seconds=1800,
    heat_shock_celsius=42.0,
    heat_shock_seconds=30,
    recover_seconds=120,
    outgrowth_ul=250.0,
    outgrowth_celsius=37.0,
    outgrowth_seconds=3600,
    plate_ul=100.0,
    dilution=10,
)

#: Where the numbers above come from, ready for a protocol's reference list.
REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Invitrogen, Gateway LR Clonase II Enzyme Mix, MAN0001032 revision A.0, for the "
        "reaction, its incubation and how it is stopped"
    ),
    Reference(
        "Invitrogen, Gateway Technology with Clonase II, MAN0000470, revised 2 April 2012, for "
        "the transformation, the plating and the strains"
    ),
    Reference(
        "Hartley, J.L., Temple, G.F. and Brasch, M.A. (2000) DNA cloning using in vitro "
        "site-specific recombination. Genome Res. 10, 1788-1795",
        url="https://doi.org/10.1101/gr.143000",
    ),
    Reference(
        "Brasch, M., Cheo, D., Hartley, J. and Temple, G., US 7,670,823 B1, FIG. 9, for the "
        "att site sequences"
    ),
)


def lr_amounts(entry: tuple[str, int], destination: tuple[str, int]) -> tuple[Amount, Amount]:
    """Return what the LR reaction takes of the entry clone and the destination vector.

    The manual gives both in nanograms, so the picomoles are worked out from each plasmid's own
    length rather than assumed.

    Parameters
    ----------
    entry, destination
        Each plasmid's name and length in base pairs.
    """
    return (
        _amount(*entry, ENTRY_NG),
        _amount(*destination, DESTINATION_NG),
    )


def _amount(name: str, length_bp: int, nanograms: float) -> Amount:
    """Return one plasmid's row: the manual's nanograms, and what they are at that length.

    The nanograms are the number the manual states, so they are kept rather than worked back
    from a rounded picomole figure.
    """
    return Amount(
        name,
        length_bp,
        pmol=round(to_pmol(nanograms, length_bp), 3),
        nanograms=nanograms,
        volume_ul=DNA_VOLUME_UL,
    )


def lr_reaction(amounts: Sequence[Amount], *, reactions: int = 1) -> ReactionTable:
    """Return the LR reaction table: the two plasmids, TE to `LR_TE_UL`, then the enzyme mix.

    The enzyme goes in after the volume is made up, which is why this is not
    `liulab_mbio.bench.reactions.reaction_table`: there the filler is the last line.

    Raises
    ------
    ValueError
        If the DNA does not leave room for the TE buffer.
    """
    fits(amounts, volume_ul=LR_TE_UL, what="the LR reaction")
    used = sum(amount.volume_ul for amount in amounts)
    return ReactionTable(
        (
            *dna_components(amounts),
            Component(TE_BUFFER, round(LR_TE_UL - used, 2), final=f"to {LR_TE_UL:g} µL"),
            Component(LR_CLONASE, LR_CLONASE_UL, final=f"to {LR_VOLUME_UL:g} µL"),
        ),
        title="LR reaction",
        reactions=reactions,
    )

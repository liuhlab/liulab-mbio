"""The Gateway reactions: what each takes, how long it runs, and how it is stopped.

What any bench shares is `liulab_mbio.bench`. Here is only what no shared table covers: each
reaction's own components, its incubation, the proteinase K that ends it, and the
transformation Invitrogen's manual prescribes, which is not the one NEB's does. BP and LR are
the same table filled from two manuals. Every number cites `docs/research/gateway-cloning.md`.
"""

from collections.abc import Mapping, Sequence
from dataclasses import replace

from liulab_mbio.bench.amounts import DNA_VOLUME_UL, Amount, dna_amount, to_pmol
from liulab_mbio.bench.reactions import dna_components, fits
from liulab_mbio.bench.steps import Transformation
from liulab_mbio.protocol.model import Component, ReactionTable, Reference

#: The BP reaction with BP Clonase II, from product sheet 11789.II.pps revision 31 October 2010
#: page 3 and MAN0000470 page 23 (note §10): microlitres, then what the DNA is made up to
#: before the enzyme goes in.
BP_VOLUME_UL = 10.0
BP_TE_UL = 8.0
BP_CLONASE_UL = 2.0

#: What the reaction takes of each DNA, picomoles: an equimolar 50 fmol of each, with the attB
#: substrate allowed down to 20 fmol. The manual prints 150 ng only as what 50 fmol of its own
#: 4.3-4.8 kb donor weighs, so a weight is worked out from each record's length rather than
#: quoted (MAN0000470 pages 21-22; note §10).
BP_PMOL = 0.05
BP_SUBSTRATE_MIN_PMOL = 0.02

#: What the reaction may not go over, ng: of donor vector, and of DNA altogether, "as excess
#: DNA will inhibit the reaction" (MAN0000470 page 21; note §10).
BP_DONOR_MAX_NG = 250.0
BP_TOTAL_MAX_NG = 500.0

#: The incubation, and how far it is extended for an attB substrate the manual calls large
#: (MAN0000470 page 23 and its troubleshooting table, pages 43-44; note §10, §16): degrees
#: Celsius, seconds and base pairs.
BP_CELSIUS = 25.0
BP_SECONDS = 3600
BP_LONG_SECONDS = 64800
BP_LONG_BP = 5000

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

#: What ends either reaction: proteinase K, supplied with the enzyme mix (MAN0001032 page 3 and
#: MAN0000470 pages 23 and 32; note §10, §11, §12). Microlitres, micrograms per microlitre,
#: Celsius and seconds.
PROTEINASE_K_UL = 1.0
PROTEINASE_K_UG_UL = 2.0
STOP_CELSIUS = 37.0
STOP_SECONDS = 600

#: What the vendor calls an efficient reaction of each kind, with the whole of it transformed
#: and plated, and the competent cells those counts assume, cfu/µg (11789.II.pps page 3,
#: MAN0001032 page 3 and MAN0000470 pages 26 and 32; note §14). What fraction of them is
#: correct is not published: see the note's gaps.
BP_COLONIES = 1500
LR_COLONIES = 5000
CELL_EFFICIENCY_CFU_UG = 100_000_000

#: What the buffer is called on the bench, and what each enzyme mix is ordered as. The catalogue
#: numbers are the 20-reaction packs (note §17).
TE_BUFFER = "TE buffer, pH 8.0"
BP_CLONASE = "Gateway BP Clonase II enzyme mix"
BP_CLONASE_CATALOG = "11789020"
LR_CLONASE = "Gateway LR Clonase II enzyme mix"
LR_CLONASE_CATALOG = "11791020"
SUPPLIER = "Invitrogen / Thermo Fisher"

#: What the vendor's one-tube protocol yields, quoted, which is why this plan stages the two
#: reactions instead: it runs BP for 4 hours and LR for 2 in the same tube, on its own timings,
#: and says to sequence every clone it gives (MAN0000470 pages 45-46; note §10, §17).
ONE_TUBE_YIELD = "at least 10-20% of the expression clones the staged route gives"

#: The strain to select the clone in. A general cloning strain is what the vendor's own support
#: page names for plating an LR reaction, and CcdB kills it, so unreacted vector and by-product
#: do not grow. Never one carrying F', whose ccdA cancels that (MAN0000470 pages 24 and 31;
#: note §9).
DEFAULT_HOST = "TOP10 chemically competent E. coli"

#: The strain a donor or destination vector is grown in, because it is the only kind that grows
#: one: "To propagate and maintain your destination vector, you must use ccdB Survival T1R
#: E. coli" (MAN0000470 page 36). What is sold for it now is ccdB Survival 2 T1R, A10460
#: (MAN0000761 revision 2.0; note §9).
PROPAGATION_HOST = "One Shot ccdB Survival 2 T1R"

#: What the manuals say about each strain note §9 names, keyed by the name to look for in the
#: host a plan is given. ``"f-prime"`` is the one prohibition the sources carry: "Do not use
#: E. coli strains that contain the F' episome (e.g. TOP10F') for transformation. These strains
#: contain the ccdA gene and will prevent negative selection with the ccdB gene" (MAN0000470
#: pages 24 and 31). ``"sensitive"`` is a strain the vendor names for plating a BP or LR
#: reaction; ``"resistant"`` is one a ccdB vector is propagated in, which no manual forbids
#: plating on -- that the selection would be lost there is an inference the note marks as one.
HOSTS: Mapping[str, str] = {
    "TOP10F'": "f-prime",
    "ccdB Survival": "resistant",
    "DB3.1": "resistant",
    "TOP10": "sensitive",
    "DH5": "sensitive",
    "OmniMAX 2": "sensitive",
}

#: The counter-screen the lost cassette allows, µg/mL: the expression clone no longer carries
#: CmR, and "A true expression clone will not grow in the presence of chloramphenicol"
#: (MAN0000223 page 5; note §13).
CHLORAMPHENICOL_UG_ML = 30

#: Invitrogen's transformation, which is not NEB's, and which is one protocol for both
#: reactions (MAN0000470 page 26; note §14, §17). No manual in the set states a thaw time, so
#: the step prints none rather than another kit's.
LR_TRANSFORMATION = Transformation(
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

#: The same transformation after BP; only what goes into the cells is named differently.
BP_TRANSFORMATION = replace(LR_TRANSFORMATION, source="the BP reaction")

#: Where the numbers above come from, ready for a protocol's reference list.
REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Invitrogen, Gateway BP Clonase II Enzyme Mix, product sheet 11789.II.pps revision 31 "
        "October 2010, for the BP reaction, its incubation and how it is stopped"
    ),
    Reference(
        "Invitrogen, Gateway LR Clonase II Enzyme Mix, MAN0001032 revision A.0, for the "
        "reaction, its incubation and how it is stopped"
    ),
    Reference(
        "Invitrogen, Gateway pDONR Vectors, MAN0000291, revised 29 March 2012, for what an "
        "entry clone is selected on and grown up from"
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


def bp_amounts(substrate: tuple[str, int], donor: tuple[str, int]) -> tuple[Amount, Amount]:
    """Return what the BP reaction takes of the attB substrate and the donor vector.

    The manual asks for an equimolar amount of each, so the picomoles are the number and the
    weights follow from the two records' own lengths.

    Parameters
    ----------
    substrate, donor
        The attB-flanked DNA's name and length in base pairs, then the donor vector's.
    """
    return (
        dna_amount(substrate[0], substrate[1], pmol=BP_PMOL),
        dna_amount(donor[0], donor[1], pmol=BP_PMOL),
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


def bp_reaction(amounts: Sequence[Amount], *, reactions: int = 1) -> ReactionTable:
    """Return the BP reaction table: the two DNAs, TE to `BP_TE_UL`, then the enzyme mix.

    Raises
    ------
    ValueError
        If the DNA does not leave room for the TE buffer.
    """
    return _clonase_reaction(
        amounts,
        title="BP reaction",
        clonase=BP_CLONASE,
        clonase_ul=BP_CLONASE_UL,
        te_ul=BP_TE_UL,
        volume_ul=BP_VOLUME_UL,
        reactions=reactions,
    )


def lr_reaction(amounts: Sequence[Amount], *, reactions: int = 1) -> ReactionTable:
    """Return the LR reaction table: the two plasmids, TE to `LR_TE_UL`, then the enzyme mix.

    Raises
    ------
    ValueError
        If the DNA does not leave room for the TE buffer.
    """
    return _clonase_reaction(
        amounts,
        title="LR reaction",
        clonase=LR_CLONASE,
        clonase_ul=LR_CLONASE_UL,
        te_ul=LR_TE_UL,
        volume_ul=LR_VOLUME_UL,
        reactions=reactions,
    )


def _clonase_reaction(
    amounts: Sequence[Amount],
    *,
    title: str,
    clonase: str,
    clonase_ul: float,
    te_ul: float,
    volume_ul: float,
    reactions: int,
) -> ReactionTable:
    """Return one Clonase II reaction: the DNA, TE to `te_ul`, then the enzyme mix on top.

    The enzyme goes in after the volume is made up, which is why this is not
    `liulab_mbio.bench.reactions.reaction_table`: there the filler is the last line.
    """
    fits(amounts, volume_ul=te_ul, what=f"the {title}")
    used = sum(amount.volume_ul for amount in amounts)
    return ReactionTable(
        (
            *dna_components(amounts),
            Component(TE_BUFFER, round(te_ul - used, 2), final=f"to {te_ul:g} µL"),
            Component(clonase, clonase_ul, final=f"to {volume_ul:g} µL"),
        ),
        title=title,
        reactions=reactions,
    )

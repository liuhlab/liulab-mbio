"""What one round of a library assembly takes at the bench, computed from its own lengths.

The volumes, times, temperatures, masses and the molar ratio below are the method's own, each
sourced in ``docs/research/protein-library-assembly.md``. Every amount is computed from them: a
mass becomes picomoles at the length of the DNA actually being cut or ligated, through
`liulab_mbio.bench.amounts`, so no quantity here is copied from the paper's own example.

Two sourced numbers are easy to get wrong, so both are named rather than written into a step:
both growth steps run at `GROWTH_CELSIUS` and not at 37 °C, and the two SPRI ratios differ.

The reactions and the cycling a round runs are here too, because every protocol of the chain
that cuts, joins or grows reads them from one place. The method publishes neither the ligase's
units nor its volume, so the last line of a ligation is the supplier's own, and the enzymes a
round uses belong to the scheme.
"""

from collections.abc import Sequence
from dataclasses import KW_ONLY, dataclass

from liulab_mbio.bench.amounts import Amount, dna_amount, to_nanograms, to_pmol
from liulab_mbio.bench.coverage import RoundCoverage
from liulab_mbio.bench.materials import Electroporation, electroporation
from liulab_mbio.bench.reactions import fits, floor_ng_ul, reaction_table
from liulab_mbio.bench.steps import listed
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.protocol.model import (
    Citation,
    Component,
    Incubation,
    ReactionTable,
    Reference,
    Source,
    Stage,
    ThermocyclerProgram,
)
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.sites import find_sites
from liulab_synbio.igga.method import Scheme

#: DNA into one digest, ng, and the volume it is cut in, with CutSmart buffer, µL.
#: METHOD DETAILS p. e4, which runs the last transfer's digest on the same numbers.
DIGEST_NG = 1000.0
DIGEST_VOLUME_UL = 50.0

#: Each of the two enzymes, µL: the first alone for `DIGEST_SECONDS` at `DIGEST_CELSIUS`, then
#: the second for a second hour. The destination digest is "the same protocol".
#: METHOD DETAILS p. e4.
ENZYME_UL = 2.5
DIGEST_SECONDS = 3600
DIGEST_CELSIUS = 37.0

#: What one µL of a product is worth, U/µL, keyed by the catalogue number that determines it —
#: size letter and all, never the enzyme name. BbsI-HF ships as R3539S/L at 20 U/µL and as
#: R3539M at 50 U/µL, where the same volume is 125 units and not 50. Each is that product's own
#: NEB specification. ``docs/research/bench-numbers.md``, "H23 re-checked".
DIGEST_STOCK_UNITS_UL: dict[str, float] = {
    "R3733L": 20.0,
    "R0629L": 20.0,
    "R3539L": 20.0,
    "R0560L": 10.0,
}

#: Which of those products the split digest runs on, against the catalogue number the shipped
#: enzyme record carries, which has no size. The paper's Key Resources Table names each one, so
#: this half of the unit count is the paper's and the concentration above is NEB's.
DIGEST_PRODUCTS: dict[str, str] = {
    "R3733": "R3733L",
    "R0629": "R0629L",
    "R3539": "R3539L",
    "R0560": "R0560L",
}

#: The documents those concentrations were read from, and the key a digest line cites them by.
DIGEST_SOURCE_KEY = "NEB specifications"
DIGEST_SOURCE = Source(
    "New England Biolabs product specification sheets",
    edition="#R3733S/L, #R3539S/L, #R0629S/L, #R0560S/L",
    read_as="plain curl",
    date="2026-10-06",
    note="docs/research/bench-numbers.md",
)

#: SPRI bead volume ratios, both eluted in water. They differ: twice the volume after the
#: digest, once after the ligation. METHOD DETAILS pp. e4-e5, which cleans the last transfer up
#: the same way, immediately before the electroporation.
SPRI_AFTER_DIGEST = 2.0
SPRI_AFTER_LIGATION = 1.0

#: Digested destination per ligation, ng, and the volume it is ligated in, µL: 20 ng per 200 µL.
#: METHOD DETAILS p. e4.
LIGATION_DESTINATION_NG = 20.0
LIGATION_VOLUME_UL = 200.0

#: Donor to destination molar ratio, 1:1, with T7 ligase in StickTogether buffer for
#: `LIGATION_SECONDS` at room temperature. METHOD DETAILS p. e4, which names neither the two
#: species the ratio is between nor the ligase's units.
MOLAR_RATIO = 1.0
LIGATION_SECONDS = 1800

#: The most purified ligation product one electroporation into Endura cells takes, ng.
#: METHOD DETAILS p. e5.
TRANSFORMATION_NG = 100.0

#: Recovery with shaking, then outgrowth from 12 to 16 hours, seconds. METHOD DETAILS p. e5.
RECOVERY_SECONDS = 3600
OUTGROWTH_SECONDS: tuple[int, int] = (43200, 57600)

#: The temperature both growth steps run at, °C. It is 30 and not 37, which is a library
#: precaution rather than an oversight. METHOD DETAILS p. e5.
GROWTH_CELSIUS = 30.0


def digest_amount(
    dna: tuple[str, int],
    *,
    nanograms: float = DIGEST_NG,
    volume_ul: float = DIGEST_VOLUME_UL,
    concentration_ng_ul: float | None = None,
) -> Amount:
    """Return what one digest takes, as picomoles of the DNA's own length.

    `dna` is the name and the length in base pairs of what is being cut, which for either
    digest is a whole plasmid. Without `concentration_ng_ul` the volume is a placeholder, so
    nothing can be shown not to fit.

    Raises
    ------
    ValueError
        If the length or the concentration is not positive, or if the DNA and the two enzymes
        do not fit `volume_ul`.

    Examples
    --------
    >>> round(digest_amount(("N part list", 5000)).pmol, 3)
    0.325
    """
    name, length_bp = dna
    amount = dna_amount(
        name,
        length_bp,
        pmol=to_pmol(nanograms, length_bp),
        concentration_ng_ul=concentration_ng_ul,
    )
    fits((amount,), volume_ul=volume_ul, taken_ul=2 * ENZYME_UL, what=f"digest of {name}")
    return amount


def pool_floor_ng_ul(
    *,
    nanograms: float = DIGEST_NG,
    volume_ul: float = DIGEST_VOLUME_UL,
    taken_ul: float = 2 * ENZYME_UL,
) -> float:
    """Return the least a pooled part list may be concentrated at, ng/µL.

    The round's digest takes `nanograms` of the pool in `volume_ul`, and its two enzymes take
    `taken_ul` of that: `liulab_mbio.bench.reactions.floor_ng_ul` on the round's own numbers.
    No source states this floor, which is why the vendors' 10 ng/µL resuspension floor does not
    settle it.

    Raises
    ------
    ValueError
        If the enzymes leave no volume for the DNA.

    Examples
    --------
    >>> round(pool_floor_ng_ul(), 2)
    22.22
    """
    return floor_ng_ul(nanograms, volume_ul=volume_ul, taken_ul=taken_ul)


def ligation_amounts(
    destination: tuple[str, int],
    donor: tuple[str, int],
    *,
    ratio: float = MOLAR_RATIO,
    nanograms: float = LIGATION_DESTINATION_NG,
    volume_ul: float = LIGATION_VOLUME_UL,
    destination_ng_ul: float | None = None,
    donor_ng_ul: float | None = None,
) -> tuple[Amount, Amount]:
    """Return what one round's ligation takes, the opened destination first.

    `destination` and `donor` are the digest fragments the round joins, each a name and a length
    in base pairs. The destination's `nanograms` become picomoles of its own length, and the
    donor gets `ratio` times as many picomoles, weighed at the donor's length: a donor shorter
    than the destination weighs less for the same number of molecules, so matching the two
    masses would not match their molecules.

    Only the DNA is measured against `volume_ul`, the ligase's own volume not being published.

    Raises
    ------
    ValueError
        If a length, a concentration or `ratio` is not positive, or if the DNA does not fit
        `volume_ul`.

    Examples
    --------
    >>> opened, released = ligation_amounts(("library", 5000), ("C part list", 1200))
    >>> opened.nanograms, released.nanograms
    (20.0, 4.8)
    """
    if ratio <= 0:
        raise ValueError(f"molar ratio must be positive, got {ratio}")
    destination_name, destination_bp = destination
    donor_name, donor_bp = donor
    pmol = to_pmol(nanograms, destination_bp)
    opened = dna_amount(
        destination_name, destination_bp, pmol=pmol, concentration_ng_ul=destination_ng_ul
    )
    released = dna_amount(donor_name, donor_bp, pmol=pmol * ratio, concentration_ng_ul=donor_ng_ul)
    fits(
        (opened, released),
        volume_ul=volume_ul,
        what=f"ligation of {donor_name} into {destination_name}",
    )
    return opened, released


def final_assembly_amounts(
    cargo: tuple[str, int],
    vector: tuple[str, int],
    *,
    vector_ng: float,
    ratio: float,
    volume_ul: float = LIGATION_VOLUME_UL,
    release_ul: float = DIGEST_VOLUME_UL,
    vector_ng_ul: float | None = None,
) -> tuple[Amount, Amount]:
    """Return what the final one-pot assembly holds, the cargo in its release tube first.

    `vector_ng` of the opened working vector becomes picomoles of its own length, and the cargo
    gets `ratio` times as many, weighed at the cargo's. The cargo is never pipetted: it arrives
    in the whole `release_ul` the release ran in, which is the volume of its row.

    Raises
    ------
    ValueError
        If a length, a concentration or `ratio` is not positive, or if the DNA and the cargo
        enzyme do not fit `volume_ul`.

    Examples
    --------
    >>> freed, opened = final_assembly_amounts(("cargo", 1076), ("vector", 9899),
    ...                                        vector_ng=75.0, ratio=2.0)
    >>> freed.nanograms, opened.nanograms
    (16.31, 75.0)
    """
    if ratio <= 0:
        raise ValueError(f"molar ratio must be positive, got {ratio}")
    cargo_name, cargo_bp = cargo
    vector_name, vector_bp = vector
    pmol = to_pmol(vector_ng, vector_bp)
    opened = dna_amount(vector_name, vector_bp, pmol=pmol, concentration_ng_ul=vector_ng_ul)
    freed = dna_amount(
        cargo_name,
        cargo_bp,
        pmol=pmol * ratio,
        concentration_ng_ul=to_nanograms(pmol * ratio, cargo_bp) / release_ul,
    )
    fits(
        (freed, opened),
        volume_ul=volume_ul,
        taken_ul=ENZYME_UL,
        what=f"final assembly into {vector_name}",
    )
    return freed, opened


def transformation_amount(
    product: tuple[str, int],
    *,
    nanograms: float = TRANSFORMATION_NG,
    concentration_ng_ul: float | None = None,
) -> Amount:
    """Return what one electroporation takes, at most `nanograms` of the ligation product.

    `product` is the name and the length in base pairs of what the round ligated, and its
    picomoles are read off that length.

    Raises
    ------
    ValueError
        If the length or the concentration is not positive.

    Examples
    --------
    >>> round(transformation_amount(("round 1 library", 6200)).pmol, 4)
    0.0262
    """
    name, length_bp = product
    return dna_amount(
        name,
        length_bp,
        pmol=to_pmol(nanograms, length_bp),
        concentration_ng_ul=concentration_ng_ul,
    )


#: Where the numbers above come from, ready for a protocol's reference list.
REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Takacsi-Nagy, O. et al. (2026) Synthetic transcription factors designed by domain "
        "recombination enhance CAR T cell antitumor function. Cell 189, 1-20, STAR Methods "
        "METHOD DETAILS pp. e4-e5. CC BY 4.0",
        url="https://doi.org/10.1016/j.cell.2026.07.054",
    ),
)


#: The buffer both digests run in, and the ligase and buffer the ligation runs in. The method
#: names all three and publishes neither the ligase's units nor either buffer's strength.
CUTSMART = "CutSmart Buffer"
LIGASE = "T7 DNA Ligase"
LIGASE_BUFFER = "StickTogether DNA Ligase Buffer"

#: What each clean-up is done with. The method gives the two volume ratios and not the product.
SPRI_BEADS = "SPRI paramagnetic beads"

#: The electrocompetent strain the method names, and who sells it.
STRAIN = "Endura ElectroCompetent Cells"
STRAIN_SUPPLIER = "Lucigen"
STRAIN_CATALOG = "60242-2"


def pulse() -> Electroporation:
    """Return the strain's own program, which belongs to the cells and not to the step.

    Raises
    ------
    LookupError
        If the package stops shipping a program for these cells, rather than printing none.
    """
    found = electroporation(STRAIN_CATALOG)
    if found is None:
        raise LookupError(f"no electroporation program ships for {STRAIN_CATALOG}")
    return found


@dataclass(frozen=True, slots=True)
class RoundBench:
    """What one round takes at the bench, computed from that round's own lengths.

    Parameters
    ----------
    number
        Which round it is, counting from one.
    position
        The position it fills, which names its part list.
    destination_digest, donor_digest
        What goes into each of the two digests: the library built so far, and the part list.
    ligation
        What the ligation takes, the opened destination first and the released part list second.
    transformation
        The most of the purified ligation one electroporation takes.
    coverage
        What this round has to cover, and what the colonies asked for leave out.
    """

    number: int
    position: str
    _: KW_ONLY
    destination_digest: Amount
    donor_digest: Amount
    ligation: tuple[Amount, Amount]
    transformation: Amount
    coverage: RoundCoverage


def choppers(scheme: Scheme) -> tuple[tuple[Enzyme, ...], tuple[Enzyme, ...]]:
    """Return the blunt enzymes each digest carries, the internal digest's first.

    A chopper belongs to the digest whose discarded piece it cuts: the internal stuffer core the
    internal digest excises, or the external stuffers the external digest leaves behind. One that
    reads a site in both is named by both, and the method has already refused one that reads a
    site in neither.
    """
    core = SequenceRecord(scheme.internal_stuffer_core)
    flanks = [
        SequenceRecord(bases) for bases in (scheme.external_stuffer_5, scheme.external_stuffer_3)
    ]
    inside = tuple(one for one in scheme.blunt if find_sites(core, one))
    outside = tuple(one for one in scheme.blunt if any(find_sites(flank, one) for flank in flanks))
    return inside, outside


def digest_reaction(
    dna: Amount,
    enzymes: Sequence[Enzyme],
    *,
    volume_ul: float = DIGEST_VOLUME_UL,
    reactions: int = 1,
) -> ReactionTable:
    """Return one digest: the DNA, each enzyme in turn, and the buffer and water that fill it.

    The buffer is one line with the water because the method names `CUTSMART` and not the
    strength it is supplied at, so its own volume is the supplier's to set.

    An enzyme whose product `DIGEST_PRODUCTS` names carries that product and its
    concentration, and what the volume is worth in units; one it does not name carries the
    volume alone.

    Raises
    ------
    ValueError
        If the DNA and the enzymes do not fit `volume_ul`.
    """
    return reaction_table(
        (dna,),
        tuple(_digest_component(one) for one in enzymes),
        volume_ul=volume_ul,
        title=f"Digest with {listed([one.name for one in enzymes])}",
        filler=f"{CUTSMART} and nuclease-free water",
        reactions=reactions,
    )


def _digest_component(enzyme: Enzyme) -> Component:
    """Return one enzyme's line of a digest, with its unit count where its product is known."""
    catalog = DIGEST_PRODUCTS.get(enzyme.catalog_number or "")
    units_ul = DIGEST_STOCK_UNITS_UL.get(catalog or "")
    if catalog is None or units_ul is None:
        return Component(enzyme.supplier_label, ENZYME_UL)
    return Component(
        enzyme.supplier_label,
        ENZYME_UL,
        stock=f"{catalog}, {units_ul:g} U/µL",
        final=f"{units_ul * ENZYME_UL:g} units",
        citation=Citation(DIGEST_SOURCE_KEY, f"#{catalog}, concentration"),
    )


def digest_program(enzymes: Sequence[Enzyme]) -> ThermocyclerProgram:
    """Return the two hours a digest runs: the first enzyme alone, then the rest added to it."""
    first, rest = enzymes[0], tuple(enzymes[1:])
    stages = [Stage((Incubation(first.name, DIGEST_CELSIUS, DIGEST_SECONDS),))]
    if rest:
        added = listed([one.name for one in rest])
        stages.append(
            Stage((Incubation(f"{first.name} and {added}", DIGEST_CELSIUS, DIGEST_SECONDS),))
        )
    return ThermocyclerProgram(
        tuple(stages), title=f"Digest with {listed([one.name for one in enzymes])}"
    )


def ligation_reaction(
    amounts: tuple[Amount, Amount],
    *,
    volume_ul: float = LIGATION_VOLUME_UL,
    reactions: int = 1,
) -> ReactionTable:
    """Return one round's ligation: the two digest fragments, then the ligase, buffer and water.

    The last line carries all three together: the method publishes neither the ligase's units nor
    its volume, so the split between them is the supplier's to set and is not invented here.

    Raises
    ------
    ValueError
        If the two fragments do not fit `volume_ul`.
    """
    return reaction_table(
        amounts,
        volume_ul=volume_ul,
        title="Ligation",
        filler=f"{LIGASE} in {LIGASE_BUFFER}, and nuclease-free water",
        reactions=reactions,
    )


def final_assembly_reaction(
    amounts: tuple[Amount, Amount],
    enzyme: Enzyme,
    *,
    volume_ul: float = LIGATION_VOLUME_UL,
) -> ReactionTable:
    """Return the final one-pot assembly: the release, the working vector, `enzyme`, then the rest.

    `amounts` is `final_assembly_amounts`' pair, so the first row is the release the cargo
    arrives in. The last line carries the ligase with its buffer and the water, as a round's
    ligation does, neither the ligase's units nor its volume being published.

    Raises
    ------
    ValueError
        If the DNA and the enzyme do not fit `volume_ul`.
    """
    return reaction_table(
        amounts,
        (_digest_component(enzyme),),
        volume_ul=volume_ul,
        title="One-pot assembly",
        filler=f"{LIGASE} in {LIGASE_BUFFER}, and nuclease-free water",
    )


def growth_program() -> ThermocyclerProgram:
    """Return the recovery and the outgrowth, both at `GROWTH_CELSIUS` and not at 37 °C.

    The outgrowth carries the shorter of the two times the method gives; the step says the range.
    """
    return ThermocyclerProgram(
        (
            Stage((Incubation("Recovery, shaking", GROWTH_CELSIUS, RECOVERY_SECONDS),)),
            Stage((Incubation("Outgrowth", GROWTH_CELSIUS, OUTGROWTH_SECONDS[0]),)),
        ),
        title=f"Recovery and outgrowth at {GROWTH_CELSIUS:g} °C",
    )

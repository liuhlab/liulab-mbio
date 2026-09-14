"""What one round of a library assembly takes at the bench, computed from its own lengths.

The volumes, times, temperatures, masses and the molar ratio below are the method's own, each
sourced in ``docs/research/protein-library-assembly.md``. Every amount is computed from them: a
mass becomes picomoles at the length of the DNA actually being cut or ligated, through
`liulab_mbio.bench.amounts`, so no quantity here is copied from the paper's own example.

Two sourced numbers are easy to get wrong, so both are named rather than written into a step:
both growth steps run at `GROWTH_CELSIUS` and not at 37 °C, and the two SPRI ratios differ.

No reaction table is built here. The method publishes neither the ligase's units nor its volume,
and the enzymes a round uses belong to the scheme.
"""

from collections.abc import Sequence

from liulab_mbio.bench.amounts import Amount, dna_amount, to_pmol
from liulab_mbio.bench.steps import listed
from liulab_mbio.protocol import Reference

#: DNA into one digest, ng, and the volume it is cut in, with CutSmart buffer, µL.
#: METHOD DETAILS p. e4.
DIGEST_NG = 1000.0
DIGEST_VOLUME_UL = 50.0

#: Each of the two enzymes, µL: the first alone for `DIGEST_SECONDS` at `DIGEST_CELSIUS`, then
#: the second for a second hour. The destination digest is "the same protocol".
#: METHOD DETAILS p. e4.
ENZYME_UL = 2.5
DIGEST_SECONDS = 3600
DIGEST_CELSIUS = 37.0

#: SPRI bead volume ratios, both eluted in water. They differ: twice the volume after the
#: digest, once after the ligation. METHOD DETAILS pp. e4-e5.
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


def _fit(
    what: str,
    amounts: Sequence[Amount],
    *,
    volume_ul: float,
    taken_ul: float = 0.0,
    taken_by: str = "",
) -> None:
    """Refuse DNA that does not fit its reaction, naming what to concentrate or scale."""
    dna_ul = sum(amount.volume_ul for amount in amounts)
    if dna_ul + taken_ul <= volume_ul:
        return
    also = f" and {taken_ul:g} µL of {taken_by}" if taken_ul else ""
    fix = "scale the reaction up"
    room = volume_ul - taken_ul
    if room > 0:
        needed = sum(amount.nanograms for amount in amounts) / room
        names = listed([amount.name for amount in amounts])
        fix = f"concentrate {names} to {needed:.3g} ng/µL or more, or {fix}"
    raise ValueError(
        f"{what}: {dna_ul:g} µL of DNA{also} exceeds the {volume_ul:g} µL reaction; {fix}"
    )


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
    _fit(
        f"digest of {name}",
        (amount,),
        volume_ul=volume_ul,
        taken_ul=2 * ENZYME_UL,
        taken_by="the two enzymes",
    )
    return amount


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
    _fit(
        f"ligation of {donor_name} into {destination_name}",
        (opened, released),
        volume_ul=volume_ul,
    )
    return opened, released


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

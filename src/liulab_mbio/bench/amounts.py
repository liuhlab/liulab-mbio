"""How much DNA a reaction takes: its weight, picomoles from nanograms, and what to pipette.

The conversion is NEBioCalculator's, through ``docs/research/primer-design-and-pcr.md``.
"""

from dataclasses import KW_ONLY, dataclass

#: NEBioCalculator's double-stranded DNA weight, g/mol: `_DUPLEX_ENDS + bp * _BASE_PAIR`. NEB's
#: manuals use 650 Da per base pair instead, which differs by about 5%, so a protocol says which.
_DUPLEX_ENDS = 36.04
_BASE_PAIR = 615.94

#: What to pipette of a fragment whose concentration is not known yet, µL.
DNA_VOLUME_UL = 1.0


def molecular_weight(length_bp: int) -> float:
    """Return the weight of a double-stranded DNA molecule, g/mol."""
    return _DUPLEX_ENDS + length_bp * _BASE_PAIR


def to_pmol(nanograms: float, length_bp: int) -> float:
    """Return how many picomoles a mass of double-stranded DNA of this length is.

    Examples
    --------
    >>> round(to_pmol(1000, 2686), 3)
    0.604
    """
    return 1000.0 * nanograms / molecular_weight(length_bp)


def to_nanograms(pmol: float, length_bp: int) -> float:
    """Return what those picomoles of double-stranded DNA of this length weigh, ng."""
    return pmol * molecular_weight(length_bp) / 1000.0


@dataclass(frozen=True, slots=True)
class Amount:
    """How much of one fragment a reaction takes, in the three units a bench needs.

    Parameters
    ----------
    name, length_bp
        The fragment's.
    pmol
        What the protocol asks for.
    nanograms
        The same amount weighed.
    volume_ul
        What to pipette at the fragment's concentration, or `DNA_VOLUME_UL` without one.

    Raises
    ------
    ValueError
        If any of the three is not positive.
    """

    name: str
    length_bp: int
    _: KW_ONLY
    pmol: float
    nanograms: float
    volume_ul: float

    def __post_init__(self) -> None:
        """Refuse an amount no one can pipette."""
        for field, value in (
            ("pmol", self.pmol),
            ("nanograms", self.nanograms),
            ("volume_ul", self.volume_ul),
        ):
            if value <= 0:
                raise ValueError(f"amount {self.name!r}: {field} must be positive")

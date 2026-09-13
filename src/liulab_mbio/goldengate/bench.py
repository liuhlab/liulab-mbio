"""What a Golden Gate experiment takes at the bench: quantities, reactions and programs.

Every number is NEB's, through ``docs/research/golden-gate-assembly.md`` for the assembly and
``docs/research/primer-design-and-pcr.md`` for the conversions, the PCR programs and the
ladders. Functions return `liulab_mbio.protocol` values, so a protocol prints them unchanged.
"""

from dataclasses import KW_ONLY, dataclass

#: NEBioCalculator's double-stranded DNA weight, g/mol: `_DUPLEX_ENDS + bp * _BASE_PAIR`. NEB's
#: manuals use 650 Da per base pair instead, which differs by about 5%, so a protocol says which.
_DUPLEX_ENDS = 36.04
_BASE_PAIR = 615.94

#: Picomoles of each fragment NEB's Golden Gate reactions take, whichever system is used.
FRAGMENT_PMOL = 0.05

#: Insert to vector molar ratio, from the E1601 manual revision 5.0_6/26 for amplicon inserts.
#: An earlier revision of the same page asked for 2:1.
INSERT_RATIO = 1.0

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
class Fragment:
    """One piece of DNA going into an assembly.

    Parameters
    ----------
    name
        What its tube is labelled.
    length_bp
        Base pairs.
    concentration_ng_ul
        Of that tube, or ``None`` when it is not measured yet.

    Raises
    ------
    ValueError
        If either number is not positive.
    """

    name: str
    length_bp: int
    _: KW_ONLY
    concentration_ng_ul: float | None = None

    def __post_init__(self) -> None:
        """Refuse a length or a concentration that is not positive."""
        if self.length_bp <= 0:
            raise ValueError(f"fragment {self.name!r}: length_bp must be positive")
        if self.concentration_ng_ul is not None and self.concentration_ng_ul <= 0:
            raise ValueError(f"fragment {self.name!r}: concentration_ng_ul must be positive")


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


def assembly_amounts(
    vector: Fragment,
    inserts: tuple[Fragment, ...],
    *,
    vector_pmol: float = FRAGMENT_PMOL,
    insert_ratio: float = INSERT_RATIO,
) -> tuple[Amount, ...]:
    """Return what to put in the assembly, vector first.

    NEB asks for `FRAGMENT_PMOL` of the destination plasmid and the same of each precloned
    insert; `insert_ratio` is the insert to vector molar ratio for amplicon inserts.
    """
    return tuple(
        _amount(fragment, pmol)
        for fragment, pmol in (
            (vector, vector_pmol),
            *((insert, vector_pmol * insert_ratio) for insert in inserts),
        )
    )


def _amount(fragment: Fragment, pmol: float) -> Amount:
    nanograms = to_nanograms(pmol, fragment.length_bp)
    volume = (
        DNA_VOLUME_UL
        if fragment.concentration_ng_ul is None
        else nanograms / fragment.concentration_ng_ul
    )
    return Amount(
        fragment.name,
        fragment.length_bp,
        pmol=pmol,
        nanograms=round(nanograms, 2),
        volume_ul=round(volume, 2),
    )

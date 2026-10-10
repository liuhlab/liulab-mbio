"""The ladder and the agarose percentage a range of bands takes.

Both are NEB's, through ``docs/research/primer-design-and-pcr.md``.
"""

from collections.abc import Mapping
from types import MappingProxyType

from mbio.protocol.model import Citation, Ladder, Note, Source

#: NEB's two ladders, each with 500 and 517 counted as the one band NEB counts them as.
LADDER_100_BP = Ladder(
    "NEB 100 bp DNA Ladder (N3231)",
    (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1517),
)
LADDER_1_KB_PLUS = Ladder(
    "NEB 1 kb Plus DNA Ladder (N3200)",
    (
        100,
        200,
        300,
        400,
        500,
        600,
        700,
        800,
        900,
        1000,
        1200,
        1517,
        2017,
        3001,
        4001,
        5001,
        6001,
        8001,
        10002,
    ),
)

#: Where the 100 bp ladder at 2% agarose gives way to the 1 kb Plus ladder at 1%, bp.
_LADDER_LIMIT = 1000


def choose_ladder(bands_bp: tuple[int, ...]) -> Ladder:
    """Return the ladder covering these bands.

    Raises
    ------
    ValueError
        If no band is given.
    """
    return LADDER_100_BP if _largest(bands_bp) < _LADDER_LIMIT else LADDER_1_KB_PLUS


def agarose_percent(bands_bp: tuple[int, ...]) -> float:
    """Return the agarose percentage these bands resolve on.

    Raises
    ------
    ValueError
        If no band is given.
    """
    return 2.0 if _largest(bands_bp) < _LADDER_LIMIT else 1.0


def _largest(bands_bp: tuple[int, ...]) -> int:
    if not bands_bp:
        raise ValueError("a gel needs at least one expected band")
    return max(bands_bp)


#: The documents the ladders and the percentages are read from, keyed as a `Citation` names
#: them. A ladder's page is keyed by its catalogue number, which is how its material row cites it.
SOURCES: Mapping[str, Source] = MappingProxyType(
    {
        "N3200": Source(
            "New England Biolabs #N3200 1 kb Plus DNA Ladder product page",
            edition="capture 2026-03-07",
            read_as="Wayback Machine",
            note="docs/research/primer-design-and-pcr.md",
        ),
        "N3231": Source(
            "New England Biolabs #N3231 100 bp DNA Ladder product page",
            edition="capture 2025-10-31",
            read_as="Wayback Machine",
            note="docs/research/primer-design-and-pcr.md",
        ),
        "agarose-resolution": Source(
            "New England Biolabs, Agarose Gel Resolution",
            edition="capture 2024-09-18",
            read_as="Wayback Machine",
            note="docs/research/primer-design-and-pcr.md",
        ),
    }
)

#: Where `agarose_percent` is read.
RESOLUTION_CITATION = Citation("agarose-resolution", "optimum resolution for linear DNA")


def resolution_note(percent: float) -> Note:
    """Return the note saying this percentage resolves the bands the step expects."""
    return Note(
        f"{percent:g}% agarose resolves bands of these sizes.", citation=RESOLUTION_CITATION
    )


#: Where each ladder's bands are read, keyed by its catalogue number.
LADDER_CITATIONS: Mapping[str, Citation] = MappingProxyType(
    {catalog: Citation(catalog, "band sizes") for catalog in ("N3200", "N3231")}
)

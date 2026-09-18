"""The ladder and the agarose percentage a range of bands takes.

Both are NEB's, through ``docs/research/primer-design-and-pcr.md``.
"""

from liulab_mbio.protocol.model import Ladder, Reference

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


#: Where the numbers above come from, ready for a protocol's reference list.
REFERENCES: tuple[Reference, ...] = (
    Reference("NEB product pages: 1 kb Plus DNA Ladder (N3200), 100 bp DNA Ladder (N3231)"),
    Reference("NEB, Agarose Gel Resolution, for the percentage a band range resolves on"),
)

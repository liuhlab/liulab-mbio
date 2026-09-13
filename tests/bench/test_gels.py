"""Gel ladders and agarose pinned against NEB's pages.

The source is `docs/research/primer-design-and-pcr.md` for the ladders.
"""

from liulab_mbio.bench.gels import LADDER_1_KB_PLUS, LADDER_100_BP, agarose_percent, choose_ladder


def test_a_band_over_a_kilobase_takes_the_wider_ladder() -> None:
    assert choose_ladder((137, 797)) is LADDER_100_BP
    assert choose_ladder((137, 3000)) is LADDER_1_KB_PLUS
    assert agarose_percent((137, 797)) == 2.0
    assert agarose_percent((137, 3000)) == 1.0

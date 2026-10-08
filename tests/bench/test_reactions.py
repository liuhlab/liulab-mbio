"""The one fill-to-volume rule every pipeline's reaction follows, and what it refuses."""

import re

import pytest

from liulab_mbio.bench.amounts import dna_amount, to_pmol
from liulab_mbio.bench.reactions import reaction_table

#: The one refusal, which names both ways out: concentrate the DNA, or scale the reaction.
REFUSAL = (
    "Assembly: 21.56 µL of components exceeds the 15 µL reaction; "
    "concentrate insert to 1.44 ng/µL or more, or scale the reaction up"
)


def test_dna_that_does_not_fit_is_refused_with_both_ways_out_named() -> None:
    dilute = dna_amount("insert", 700, pmol=0.05, concentration_ng_ul=1.0)

    with pytest.raises(ValueError, match=re.escape(REFUSAL)):
        reaction_table((dilute,), volume_ul=15.0, title="Assembly")


def test_a_rows_picomoles_print_three_figures_and_not_six() -> None:
    """A weighed picomole is unrounded, so the row states it as the page does."""
    weighed = dna_amount("pUC19", 2686, pmol=to_pmol(1000.0, 2686))
    table = reaction_table((weighed,), volume_ul=15.0)

    assert table.components[0].final == "0.604 pmol (1000 ng)"

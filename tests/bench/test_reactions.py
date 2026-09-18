"""The one fill-to-volume rule every pipeline's reaction follows, and what it refuses."""

import re

import pytest

from liulab_mbio.bench.amounts import dna_amount
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

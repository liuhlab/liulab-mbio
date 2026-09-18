"""Where a junction lands in the product, and the site the join puts back there.

A fragment writes both of its overhangs on the top strand, so which strand overhangs is what
says where the bases the two ends paired on lie: the following piece's first bases for a 5'
overhang, the last bases of the piece before for a 3' one, and none at all for a blunt end. One
enzyme cuts both ends of these, so the site every join restores is that enzyme's own.
"""

import pytest

from liulab_mbio.cloning.restriction.digest import excised, opened, resolve
from liulab_mbio.cloning.restriction.ligation import ligate
from liulab_mbio.sequence import SequenceRecord

#: Bases spelling no site of the enzymes named below: the vector's, and what it is given.
FILLER = "TAAGGTCA" * 2
INSERT = "GTTCAGTA" * 2

#: One shipped enzyme of each end type, its site, and where its two junctions land.
END_TYPES = [
    pytest.param("SphI", "GCATGC", (17, 39), id="3'"),
    pytest.param("EcoRI", "GAATTC", (17, 39), id="5'"),
    pytest.param("SmaI", "CCCGGG", (19, 41), id="blunt"),
]


@pytest.mark.parametrize(("enzyme", "site", "positions"), END_TYPES)
def test_a_junction_sits_on_the_bases_its_ends_paired_on_and_spells_the_site_they_restore(
    enzyme: str, site: str, positions: tuple[int, int]
) -> None:
    chosen = resolve([enzyme])
    vector = SequenceRecord(FILLER + site + FILLER * 4, topology="circular", name="pVector")
    holder = SequenceRecord(
        FILLER + site + INSERT + site + FILLER, topology="circular", name="pHolder"
    )
    backbone = opened(vector, chosen)[0]
    made = ligate(backbone, excised(holder, chosen, into=backbone)[0])
    # The junctions check reads the product across each join for the bases it paired on.
    assert made.status == "pass"
    assert made.junction_positions == positions
    for one in made.junctions:
        assert (one.spells, one.enzyme) == (site, enzyme)
        assert one.site is not None
        assert made.product.extract(one.site) == site
        # The paired bases lie inside the site, wherever in it the cut left them.
        assert one.site.start <= one.start
        assert one.end <= one.site.end

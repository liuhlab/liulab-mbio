"""Where a junction lands in the product, and the site the join puts back there.

A fragment writes both of its overhangs on the top strand, so which strand overhangs is what
says where the bases the two ends paired on lie: the following piece's first bases for a 5'
overhang, the last bases of the piece before for a 3' one, and none at all for a blunt end. One
enzyme cuts both ends of the records written in code, so the site every join restores is that
enzyme's own; the pUC19 and GFP fixtures are where a two-enzyme product is pinned.
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


def test_the_product_carries_the_insert_once_with_everything_both_pieces_annotate(
    puc19, source, gfp
):
    chosen = resolve(["EcoRI", "BamHI"])
    backbone = opened(puc19, chosen)[0]
    insert = excised(source, chosen, into=backbone)[0]
    made = ligate(backbone, insert, name="pUC19-pTrc-GFP")
    product = made.product
    assert len(product) == len(backbone.bases) + len(insert.bases)
    assert product.sequence.count(insert.bases) == 1
    coding = next(one for one in product.features if one.name == gfp.name)
    assert product.extract(coding.segments[0]) == gfp.sequence
    # What each piece annotated came with it: the vector draws its M13 sites as features.
    assert {"M13 fwd", "M13 rev"} <= {one.name for one in product.features}
    assert made.status == "pass"


def test_each_junction_is_marked_on_the_product_and_spells_the_site_its_ends_came_from(
    puc19, source
):
    chosen = resolve(["EcoRI", "BamHI"])
    backbone = opened(puc19, chosen)[0]
    made = ligate(backbone, excised(source, chosen, into=backbone)[0])
    assert [(one.start, one.spells, one.enzyme) for one in made.junctions] == [
        (396, "GAATTC", "EcoRI"),
        (1119, "GGATCC", "BamHI"),
    ]
    drawn = {one.name: one for one in made.product.features if one.name.endswith("junction")}
    assert set(drawn) == {"GAATTC junction", "GGATCC junction"}
    for one in made.junctions:
        feature = drawn[f"{one.spells} junction"]
        assert made.product.extract(feature.segments[0]) == one.spells
        assert one.span is not None
        assert made.product.extract(one.span) == one.overhang

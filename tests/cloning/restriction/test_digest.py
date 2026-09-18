"""What a digest leaves, over records small enough to read.

Two ends this method leaves may anneal to each other, and neither the backbone closing on itself
nor an insert that goes in either way round is a refusal. What is refused is a digest that leaves
no one piece to take.
"""

import pytest

from liulab_mbio.cloning.restriction.digest import (
    cut,
    excised,
    opened,
    resolve,
    reversible,
    self_closing,
    turned,
)
from liulab_mbio.cloning.restriction.ligation import ligate
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.sites import find_sites

#: Bases spelling no site of any enzyme named below, for padding a record written in code.
FILLER = "ACGT" * 6


def plasmid(sequence: str, name: str) -> SequenceRecord:
    """A circular record written in code."""
    return SequenceRecord(sequence, topology="circular", name=name)


#: One enzyme cutting both ends, and two leaving the same 5' TCGA: the vector and the plasmid the
#: insert is cut out of, for each.
ONE_ENZYME = (
    ["EcoRI"],
    plasmid(FILLER * 2 + "GAATTC" + FILLER * 8, "pOneSite"),
    plasmid(FILLER * 2 + "GAATTC" + FILLER * 2 + "GAATTC" + FILLER * 2, "pTwoSites"),
)
COMPATIBLE_PAIR = (
    ["SalI", "XhoI"],
    plasmid(FILLER * 2 + "GTCGAC" + FILLER * 8 + "CTCGAG" + FILLER * 2, "pSalXho"),
    plasmid(FILLER * 2 + "GTCGAC" + FILLER * 2 + "CTCGAG" + FILLER * 2, "pSalXhoInsert"),
)


@pytest.mark.parametrize(("enzymes", "vector", "holder"), [ONE_ENZYME, COMPATIBLE_PAIR])
def test_a_backbone_whose_two_ends_anneal_is_opened_and_its_insert_goes_in_either_way(
    enzymes: list[str], vector: SequenceRecord, holder: SequenceRecord
) -> None:
    chosen = resolve(enzymes)
    backbone = opened(vector, chosen)[0]
    insert = excised(holder, chosen, into=backbone)[0]
    assert self_closing(backbone)
    assert reversible(backbone, insert)


def test_a_directional_pair_neither_closes_on_itself_nor_reverses():
    chosen = resolve(["EcoRI", "BamHI"])
    vector = plasmid(FILLER * 2 + "GAATTC" + FILLER * 8 + "GGATCC" + FILLER * 2, "pEcoBam")
    holder = plasmid(FILLER * 2 + "GAATTC" + FILLER * 2 + "GGATCC" + FILLER * 2, "pEcoBamInsert")
    backbone = opened(vector, chosen)[0]
    insert = excised(holder, chosen, into=backbone)[0]
    assert not self_closing(backbone)
    assert not reversible(backbone, insert)


def test_an_insert_turned_over_is_cut_again_and_ligates_in_with_both_sites_put_back():
    enzymes, vector, _ = ONE_ENZYME
    chosen = resolve(enzymes)
    backbone = opened(vector, chosen)[0]
    body = "AAACCCTTT"
    holder = plasmid(FILLER * 2 + "GAATTC" + body + "GAATTC" + FILLER * 2, "pLopsided")
    insert = excised(holder, chosen, into=backbone)[0]
    over = turned(insert)
    assert insert.bases == "AATTC" + body + "G"
    # The bottom strand read 5' to 3' keeps its own AATT at the front, where reversing the top
    # strand in place would put it at the back and leave neither site whole.
    assert over.bases == "AATTC" + "AAAGGGTTT" + "G"
    assert (over.name, over.left_end, over.right_end) == (
        insert.name,
        insert.left_end,
        insert.right_end,
    )
    back = ligate(backbone, over).product
    assert [site.start for site in find_sites(back, chosen)] == [48, 63]


def test_one_enzyme_cuts_the_insert_out_on_its_own_and_a_further_cut_is_refused():
    enzymes, _, holder = ONE_ENZYME
    chosen = resolve(enzymes)
    # Two sites of the one enzyme are the two cuts that release the piece between them.
    pieces = excised(holder, chosen, into=opened(holder, chosen)[0])
    insert, rest = pieces
    assert insert.length < rest.length
    assert insert.length + rest.length == len(holder)
    assert (insert.left_enzyme.name, insert.right_enzyme.name) == ("EcoRI", "EcoRI")
    # A third cut falls inside that piece, which is what the rule is there for. The vector is
    # held to less, so this is the insert's own digest.
    third = plasmid(holder.sequence + "GAATTC" + FILLER, "pThreeSites")
    with pytest.raises(ValueError, match=r"EcoRI cuts pThreeSites 3 time\(s\), at "):
        cut(third, chosen)

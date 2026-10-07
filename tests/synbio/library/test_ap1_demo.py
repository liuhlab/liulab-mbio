"""The AP-1 demo planned end to end, from the amino acid sequences and nothing else.

The inputs are `docs/examples/ap1-library`, specified by `docs/research/ap1-demo-project.md`:
72 proteins, a project file and a destination. Every assertion here is on the planned product and
the blocks that make it, not on how the design reached them.
"""

from collections import Counter
from pathlib import Path

import pytest

from liulab_mbio.enzymes import get_enzyme
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.sites import find_sites
from liulab_mbio.translate import translate
from liulab_synbio.library import plan_library

DEMO = Path(__file__).parents[3] / "docs" / "examples" / "ap1-library"

#: The enzymes the method reserves: the three a round uses, and BsmBI, which seats a part in its
#: carrier and so must not cut a designed block either.
RESERVED = ("BsaI", "BsmBI", "BbsI", "SrfI", "PmeI")

#: What each block's own stuffers put there, and nothing else: BsaI releases the part, BbsI
#: opens the destination it becomes, SrfI and PmeI shred what a round throws away.
EXPECTED_SITES = {"BsaI": 2, "BbsI": 2, "SrfI": 1, "PmeI": 2}


@pytest.fixture(scope="module")
def plan():
    return plan_library(DEMO / "project.json")


def test_the_demo_plans_every_part_and_the_whole_library(plan):
    assert [len(one) for one in plan.part_lists] == [24, 24, 24]
    assert len(plan.parts) == 72
    assert plan.constructs == 24**3 == 13824
    assert [one.position for one in plan.rounds] == ["N", "DBD", "C"]
    assert plan.status == "pass"


def test_no_reserved_enzyme_reads_a_part_outside_its_stuffers(plan):
    enzymes = tuple(get_enzyme(one) for one in RESERVED)
    for part in plan.parts:
        assert not find_sites(SequenceRecord(part.coding_sequence), enzymes), part.name
        found = find_sites(SequenceRecord(part.sequence), enzymes)
        assert Counter(site.enzyme.name for site in found) == EXPECTED_SITES, part.name


def test_the_product_reads_in_frame_from_its_first_part_to_its_last_barcode(plan):
    bases = str(plan.product.sequence)
    start = bases.index(plan.representative_parts[0].coding_sequence)
    stretch = bases[start : plan.rounds[-1].retained.end]
    assert len(stretch) % 3 == 0
    assert "*" not in translate(stretch)


def test_the_product_keeps_one_barcode_a_round_in_reverse_order(plan):
    bases = str(plan.product.sequence)
    scar = plan.scheme.cloning_scar
    codes = [one.barcode for one in plan.representative_parts]
    assert scar.join(reversed(codes)) in bases
    retained = bases[plan.rounds[-1].retained.start : plan.rounds[-1].retained.end]
    assert len(retained) == plan.project.retained_length == 75
    assert "*" not in translate(retained)

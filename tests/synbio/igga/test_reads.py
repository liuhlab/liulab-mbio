"""The reads of a finished library: what each pair spans, how long it runs, and on what.

The plan is the small three-by-two build `test_plan` writes, so a whole library is simulated
and both pairs are designed against real records rather than against a stub.
"""

import pytest

from liulab_mbio.primers.placement import amplicon_sizes
from liulab_mbio.sequence import SequenceRecord
from liulab_synbio.igga.plan import plan_igga
from liulab_synbio.igga.reads import read_pairs, read_sheet
from liulab_synbio.igga.vector import working_vector

from .test_plan import LISTS, inputs, project, scheme  # noqa: F401


@pytest.fixture(scope="module")
def plan(inputs):  # noqa: F811
    return plan_igga(project(inputs), parts=LISTS)


@pytest.fixture(scope="module")
def pairs(plan):
    return read_pairs(plan.scheme, plan.rounds, plan.working)


def test_linkage_spans_the_whole_cargo_and_representation_spans_the_block(pairs, plan):
    """The two reads differ in what they have to carry, which is what sets each one's platform."""
    assert pairs.linkage.amplicon_length > pairs.representation.amplicon_length
    assert pairs.linkage.platform == "long read"
    assert pairs.representation.platform == "short read"
    block = plan.project.barcode_block_length
    assert pairs.representation.amplicon_length > block


def test_both_pairs_amplify_one_product_on_the_record_they_were_designed_against(pairs, plan):
    """A pair that amplified twice would count one member as two."""
    product = plan.rounds[-1].product
    for one in (pairs.linkage, pairs.representation):
        assert amplicon_sizes(one.forward, one.reverse, product) == (one.amplicon_length,)


def test_the_representation_forward_primer_lies_wholly_in_the_retained_stuffer(pairs, plan):
    """It has to read the same for every member, so it may not touch a member's coding bases."""
    stuffer = plan.rounds[-1].stuffer
    site = pairs.representation.forward.binding_sites[0]
    assert site.start >= stuffer.start
    assert site.end <= stuffer.end


def test_a_build_with_no_working_vector_designs_no_final_pair(pairs):
    """There is nowhere to design it against, which is what H31 already says."""
    assert pairs.final_representation is None
    assert len(pairs.designed) == 2


def test_the_final_pair_keeps_the_forward_anchor_and_moves_the_reverse_one(plan):
    """The reverse anchor was in the library backbone, which the move leaves behind."""
    stock = SequenceRecord("ACGATCGTTA" * 20, topology="circular", name="pWORK")
    working = working_vector(stock, [], site=(0, 1))

    moved = read_pairs(plan.scheme, plan.rounds, working).final_representation

    assert moved is not None
    before = read_pairs(plan.scheme, plan.rounds).representation
    stuffer = plan.rounds[-1].stuffer
    retained = plan.rounds[-1].product.bases(stuffer.start, stuffer.end)
    assert moved.forward.sequence in retained
    assert before.forward.sequence in retained
    assert moved.reverse.sequence != before.reverse.sequence
    assert moved.platform == "short read"


def test_the_sheet_carries_every_designed_primer_with_its_amplicon(pairs):
    """What the plan writes beside the other sheets, one primer a row."""
    rows = read_sheet(pairs).splitlines()

    assert rows[0].split("\t") == [
        "read",
        "name",
        "sequence",
        "length",
        "amplicon_bp",
        "platform",
    ]
    assert len(rows) == 1 + 2 * len(pairs.designed)
    assert rows[1].split("\t")[2] == pairs.linkage.forward.sequence


def test_a_library_with_no_rounds_is_read_by_nothing(plan):
    """A read is designed against a record, and there is none."""
    with pytest.raises(ValueError, match="rounds that built it"):
        read_pairs(plan.scheme, ())

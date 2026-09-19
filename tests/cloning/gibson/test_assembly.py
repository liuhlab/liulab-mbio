"""The parts and the product, as units: what a tail does, and what joining the parts leaves.

These are the plan's own parts, built here from the same records and thresholds, so a failure
names the function that decided the value rather than the pipeline that read it.
"""

import pytest

from liulab_mbio.cloning.gibson.assembly import (
    Assembly,
    Part,
    amplify,
    assemble,
    given_vector,
    open_vector,
)
from liulab_mbio.cloning.gibson.bench import NEBUILDER_HIFI
from liulab_mbio.cloning.gibson.design import overlap_after, overlap_before
from liulab_mbio.primers.thresholds import THRESHOLDS_FOR
from liulab_mbio.sequence import SequenceRecord, reverse_complement

#: Where the fixture's own MCS feature sits, which is the span the insert replaces.
MCS = (395, 452)

#: What the plan designs its amplifications by.
AMPLIFICATION = THRESHOLDS_FOR["amplification"]


@pytest.fixture(scope="module")
def overlaps(puc19: SequenceRecord) -> tuple[str, str]:
    """The two overlaps the vector spells either side of that span."""
    rule = NEBUILDER_HIFI.tiers[0].overlap
    return overlap_before(puc19, MCS[0], rule), overlap_after(puc19, MCS[1], rule)


@pytest.fixture(scope="module")
def parts(
    puc19: SequenceRecord, gfp: SequenceRecord, overlaps: tuple[str, str]
) -> tuple[Part, Part]:
    """The pUC19 backbone opened across its MCS, and GFP tailed with the bases it meets it on."""
    left, right = overlaps
    return (
        open_vector(puc19, *MCS, name="pUC19 backbone", thresholds=AMPLIFICATION),
        amplify(
            gfp,
            0,
            len(gfp),
            left_tail=left,
            right_tail=right,
            name="GFP",
            thresholds=AMPLIFICATION,
        ),
    )


@pytest.fixture(scope="module")
def built(parts: tuple[Part, Part]) -> Assembly:
    """Those two parts joined, which is the plasmid the common path makes."""
    return assemble(parts, name="pUC19-GFP")


def test_a_part_carries_its_neighbours_bases_and_none_of_its_own(parts, overlaps, gfp):
    backbone, insert = parts
    left, right = overlaps
    assert (insert.left_tail, insert.right_tail) == (left, right)
    assert insert.forward.sequence.startswith(left)
    assert insert.reverse.sequence.startswith(reverse_complement(right))
    assert insert.amplicon.sequence == left + gfp.sequence + right
    # The amplicon is longer than what the part puts into the product, by the two tails.
    assert insert.length == insert.fragment_length + len(left) + len(right)
    # A junction with the vector on either side takes the vector's own bases, so it needs none.
    assert (backbone.left_tail, backbone.right_tail) == ("", "")


def test_the_tm_of_a_tailed_primer_is_read_from_its_annealing_region_alone(parts):
    report = parts[1].report.forward
    annealing = report.primer.binding_sites[0]
    assert len(report.primer.sequence) == 42
    assert annealing.end - annealing.start == 26
    assert report["length"].value == 26
    # The whole oligo's melting temperature is reported beside it, and carries no verdict.
    assert report["tm_full"].value > report["tm"].value
    assert report["tm_full"].status is None
    # The structures are judged on the whole oligo, tail included.
    assert report["hairpin"].status is not None


def test_a_vector_already_linear_is_taken_as_the_opened_part(gfp, puc19):
    backbone = given_vector(gfp, name="GFP backbone")
    # No PCR opens it, so it has no primers, nothing to judge and no template to cut.
    assert not backbone.amplified
    assert (backbone.forward, backbone.reverse, backbone.report) == (None, None, None)
    assert not backbone.dpni
    assert (backbone.span, backbone.bases) == ((0, len(gfp)), gfp.sequence)
    with pytest.raises(ValueError, match="circular"):
        given_vector(puc19)


def test_the_product_carries_every_feature_at_its_new_coordinates(built, puc19, gfp):
    product = built.product
    assert len(product) == len(puc19) - (MCS[1] - MCS[0]) + len(gfp) == 3346
    features = {
        one.name: [(bit.start, bit.end) for bit in one.segments] for one in product.features
    }
    # The vector's reporter meets the replaced span twice, so it keeps a segment either side.
    reporter = next(name for name in features if name.startswith("lacZ"))
    assert features[reporter] == [(145, 395), (1112, 1129)]
    assert features["GFP"] == [(395, 1112)]
    # A feature lying wholly inside the span the insert replaced is gone.
    assert "MCS" not in features
    # One further round the plasmid, shifted by what the insert added and the span removed.
    assert features["M13 rev"] == [(1124, 1141)]
    # The vector keeps its own origin, so no junction sits at base zero.
    assert product.sequence[: MCS[0]] == puc19.sequence[: MCS[0]]
    assert product.topology == "circular"
    assert all(one.start > 0 for one in built.junctions)


def test_a_part_gives_way_where_the_shared_bases_end_and_not_where_they_begin(built, gfp):
    # Both overlaps are the vector's own bases: the first sits at the end of the backbone's
    # share and the last at the start of it. Reading the junction starts as boundaries instead
    # would put a screening primer 16 bases out and turn vector bases round with the insert.
    first, last = built.junctions
    assert built.junction_positions == (first.start, last.start) == (379, 1112)
    assert built.boundaries == (first.end, last.start)
    assert built.boundaries == built.insert_span == (MCS[0], MCS[0] + len(gfp))


def test_each_junction_is_marked_and_says_which_part_the_bases_came_from(built):
    marks = {one.name: one for one in built.product.features if one.name.endswith("overlap")}
    assert set(marks) == {"pUC19 backbone-GFP overlap", "GFP-pUC19 backbone overlap"}
    for junction in built.junctions:
        mark = marks[f"{junction.before}-{junction.after} overlap"]
        assert (mark.segments[0].start, mark.segments[0].end) == (junction.start, junction.end)
        assert junction.taken_from in mark.qualifiers["note"][0]

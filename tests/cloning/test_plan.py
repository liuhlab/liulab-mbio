"""Where every cloning method reads its insertion site from, on records built here."""

import pytest

from liulab_mbio.cloning.plan import insertion_span
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Topology


def vector(*segments: Segment, topology: Topology = "circular") -> SequenceRecord:
    """A 100 bp vector whose `MCS` feature lies over `segments`."""
    return SequenceRecord(
        "ACGT" * 25,
        topology=topology,
        name="pSmall",
        features=(Feature("MCS", "misc_feature", segments),),
    )


@pytest.mark.parametrize(
    "segments",
    [(Segment(90, 108),), (Segment(90, 98), Segment(1, 8))],
    ids=["one segment", "two segments"],
)
def test_a_named_site_across_the_origin_is_refused_and_the_refusal_says_so(segments):
    with pytest.raises(ValueError, match="'MCS' runs across the origin of pSmall"):
        insertion_span(vector(*segments), None)


@pytest.mark.parametrize(
    ("segments", "span"),
    [
        ((Segment(90, 100),), (90, 100)),
        ((Segment(10, 20), Segment(30, 40)), (10, 40)),
        ((Segment(10, 20), Segment(15, 30)), (10, 30)),
        ((Segment(10, 20), Segment(12, 15)), (10, 20)),
    ],
    ids=[
        "ends at the last base",
        "two segments in order",
        "two segments overlapping",
        "one segment inside another",
    ],
)
def test_a_named_site_on_one_side_of_the_origin_runs_from_its_first_base_to_its_last(
    segments, span
):
    assert insertion_span(vector(*segments), "MCS") == span


def test_a_span_given_across_the_origin_is_refused_as_lying_outside_the_vector():
    with pytest.raises(ValueError, match="the insertion site 90-108 does not lie inside 100 bases"):
        insertion_span(vector(Segment(10, 20)), (90, 108))


def test_a_named_site_at_both_ends_of_a_linear_vector_is_refused_without_naming_an_origin():
    linear = vector(Segment(90, 98), Segment(1, 8), topology="linear")
    with pytest.raises(
        ValueError, match="pSmall is linear, and the insertion site 'MCS' lies in two pieces"
    ) as refused:
        insertion_span(linear, "MCS")
    assert "origin" not in str(refused.value)

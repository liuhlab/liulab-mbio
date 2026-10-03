"""What the map reads off an item's spans: where it starts and ends, and the pieces it lies in."""

from liulab_mbio.plot import layers
from liulab_mbio.sequence import Feature, Segment, SequenceRecord


def test_a_join_whose_segments_overlap_spans_them_and_not_the_whole_circle() -> None:
    joined = Feature("join", "CDS", (Segment(99, 200), Segment(199, 300)))
    [item] = layers.items(
        SequenceRecord("ACGT" * 250, topology="circular", features=(joined,)),
        primers=False,
        cut_sites=False,
    )

    laid = layers.pieces(item, 0, 1000, 1000, circular=True)

    assert [(one.start, one.end) for one in laid] == [(99, 200), (199, 300)]
    assert layers.hull(item.spans, 1000) == (99, 300)

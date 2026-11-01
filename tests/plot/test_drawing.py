"""A record drawn through `draw_map`: what it hands molviz that the command does not."""

from mbio.plot import draw_map
from mbio.sequence import SequenceRecord

from ..conftest import Handed


def test_a_feature_across_the_origin_is_drawn_from_its_first_segment_to_its_last(
    across_origin: SequenceRecord, handed: Handed
) -> None:
    # Feature g reads 91..98 then 2..8 on 100 bases, so it runs from index 90 to 108; its name
    # matches whatever its case.
    draw_map(across_origin, region="G", highlight="g")
    [(_, options)] = handed
    assert (options["span"], options["highlight"]) == ({"start": 90, "end": 108}, ["g"])

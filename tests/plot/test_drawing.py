"""A record drawn through `draw_map`: what it hands molviz that the command does not."""

from pathlib import Path
from typing import Any

import molviz

from mbio.plot import draw_map
from mbio.sequence import Feature, Segment, SequenceRecord

from ..html import blocks


def _options(drawing: molviz.Drawing, page: Path) -> dict[str, Any]:
    drawing.write(page)
    return blocks(page.read_text(encoding="utf-8"))["options"]


def test_a_feature_across_the_origin_is_drawn_from_its_first_segment_to_its_last(
    tmp_path: Path,
) -> None:
    # On a 100 bp record, one reading 91..98 then 2..8 runs from index 90 to 108.
    ori = Feature("ori", "rep_origin", (Segment(90, 98), Segment(101, 108)))
    record = SequenceRecord("ACGT" * 25, topology="circular", name="small", features=(ori,))
    drawn = draw_map(record, region="ORI", highlight="ori")
    options = _options(drawn, tmp_path / "map.html")
    assert (options["span"], options["highlight"]) == ({"start": 90, "end": 108}, ["ori"])

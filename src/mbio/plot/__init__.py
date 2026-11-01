"""Draw a sequence record as a map through molviz.

`draw_map` is the way in: it takes a record, or a file to read one from, hands it to molviz with
the options chosen, and returns molviz's `Drawing`, whose `write` writes the format a file's
suffix names. `converter` writes a record in molviz's own format, as `molviz_record`.

The other submodules are mbio's own layout and drawing, which no way in uses. `layers`
resolves what a record draws, with its names, colours and hover details; `circular` lays those
items out as shapes round a circle, `linear` along a line, and `sequence_view` base by base in
rows; `labels` keeps labels apart, on bare boxes; `fonts` measures text in the pinned faces; `svg`
holds the shapes and writes them; `page` wraps the SVG in one offline HTML page; `convert` turns
it into a PNG or a PDF.

Only the ways in are re-exported. Everything else is imported by module.
"""

from mbio.plot.converter import molviz_record
from mbio.plot.drawing import draw_map

__all__ = ["draw_map", "molviz_record"]

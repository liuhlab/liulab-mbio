"""Draw a sequence record as a map.

`draw_map` is the way in: it takes a record, or a file to read one from, and lays it out once;
`Drawing.write` writes the result in the format a file's suffix names.

The submodules are the stages it is made of. `layers` resolves what a record draws, with its
names, colours and hover details; `circular` lays those items out as shapes round a circle, and
`linear` along a line; `labels` keeps labels apart, on bare boxes; `fonts` measures text in the
pinned faces; `svg` holds the shapes and writes them; `page` wraps the SVG in one offline HTML
page; `convert` turns it into a PNG or a PDF. A drawing is laid out as geometry before anything
draws, per `docs/adr/0006-drawings-as-geometry.md`.

Only the way in is re-exported. Everything else is imported by module.
"""

from liulab_mbio.plot.drawing import Drawing, draw_map

__all__ = ["Drawing", "draw_map"]

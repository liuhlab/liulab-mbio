"""Draw a sequence record as a map, and a plate as a map of its wells.

`draw_map` is the way in: it takes a record, or a file to read one from; `Drawing.write` writes it
in the format a file's suffix names, laying out each distinct layout once.

The submodules are the stages it is made of. `layers` resolves what a record draws, with its
names, colours and hover details; `circular` lays those items out as shapes round a circle,
`linear` along a line, and `sequence_view` base by base in rows; `labels` keeps labels apart, on
bare boxes; `fonts` measures text in the pinned faces; `svg` holds the shapes and writes them;
`page` wraps the SVG in one offline HTML page; `convert` turns it into a PNG or a PDF; `plate`
lays a plate's wells out on the same substrate, which `draw_plate` is the way in to. A drawing
is laid out as geometry before anything draws, per `docs/adr/0006-drawings-as-geometry.md`.

Only the way in is re-exported. Everything else is imported by module.
"""

from liulab_mbio.plot.drawing import Drawing, PlateDrawing, draw_map, draw_plate

__all__ = ["Drawing", "PlateDrawing", "draw_map", "draw_plate"]

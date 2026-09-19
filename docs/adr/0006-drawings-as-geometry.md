---
search:
  exclude: true
---

# A drawing is laid out as geometry in our own SVG before anything draws

`liulab_mbio.plot` computes every shape of a map in Python, in points: arrows, ticks, label boxes
and their leaders, with each label as wide as the pinned font the package ships measures it. Only
then does `svg` write those shapes, and the page, the PNG and the PDF all draw that one SVG. No
plotting library draws a map, and nothing is laid out in the browser.

The reason is the promise that no label overlaps another label or the drawing, and no leader
crosses a label, on any record. A rule keeps that only if it knows each label's true size and
owns every shape the label must clear, and a test checks it only on numbers it can read. As
geometry, the label rule holds by construction, and the tests assert it on boxes and leaders in
milliseconds rather than on pixels.

## Considered options

- **A plotting library**, such as pyCirclize, pyGenomeViz, DNA Features Viewer or matplotlib.
  Measured on crowded plasmids, they overlapped labels or ran leaders through them. None splits a
  feature across the origin or knows a primer or a cut site, and the HTML page would still have to
  draw the map a second time.
- **Laying out in the browser**, measuring text as it draws. A PNG or PDF would then need a headless
  browser, the layout could not be tested in Python, and one record would lay out differently
  wherever the fonts differ.

## Consequences

- The package draws what a library would have given: arcs, arrowheads, ticks and text placement.
- A width holds only while every output draws with the faces it was measured in. The page embeds
  them and pins each line of text to its width with `textLength`; a PNG or PDF draws each letter
  as its outline, so its text cannot be searched or copied.
- A switch in the page hides shapes where they are. Moving labels when a layer is switched off
  would need a layout per combination of switches, so each distinct layout runs once and is kept.
- The page zooms the line by carrying it laid out along lines two, four and eight times as long,
  and fits a narrow page by carrying the sequence view at half its row width too. Each is a
  layout of its own, run once and kept, and a PNG or PDF draws only the first. The circle zooms by
  scaling its one drawing.

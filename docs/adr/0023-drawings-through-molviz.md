---
search:
  exclude: true
---

# mbio hands its drawings to a pinned liulab-molviz release as data, and draws nothing itself

Every drawing mbio makes, a map, a sequence view, a protocol's figure, gel and plate, and a
clone's verification, is drawn by `liulab-molviz`, pinned to one exact release. mbio converts
what it knows into molviz's input format and calls molviz; it lays out no shape, measures no
text and writes no SVG. This replaces the decision that mbio lay out every drawing as its own
geometry in Python.

What mbio hands over is its biology: every shipped enzyme's cuts, a CDS's translation read from
the bases, the bases a primer does not pair with, each result's reads and disagreements. molviz
picks the unique cutters, places the labels and hides what will not fit. A drawing need molviz
lacks, or a look that differs from what mbio drew, is a molviz issue, never mbio code.

## Why

A page carried every layout it might show, each laid out in Python in advance: one map page held
the circle, the line at four zoom levels and the sequence view at two row widths, so every base
appeared four times. A viewer lays out in the browser from data a few kilobytes long, and the
same code writes the PNG and PDF. The reasons for laying out in Python no longer held: a
JavaScript viewer exports PNG and PDF, its label rule is tested in JavaScript, and identical
pixels across browsers is not a goal. Other lab packages draw sequences too, so the viewer lives
in its own repository and is written once.

## Considered options

- **Keep mbio's own geometry.** It kept the label promise testable in Python, at the cost of
  pages several times the size of their data and a second, JavaScript copy of the drawing for the
  page's switches.
- **A plotting library or an open-source viewer.** None splits a feature across the origin,
  draws a trace, or keeps labels clear on a crowded plasmid.

## Consequences

- molviz's format is the one definition. The converters, `molviz_record`, `molviz_verification`
  and those in `mbio.protocol.converter`, return plain dicts and lists, and no mbio type mirrors
  the format. Each layer converts what it owns.
- The pin is exact (`==`), because a drawing's bytes are that release's. Moving it regenerates
  every docs example that holds a drawing.
- A test holds a converter's output to `molviz.check`, `check_gel` or `check_plate`. A drawing
  starts a JavaScript runtime, so a test draws only where the drawn page is the promise.
- mbio ships no fonts and depends on neither vl-convert nor pypdf; molviz brings both.

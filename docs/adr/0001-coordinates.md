---
search:
  exclude: true
---

# Half-open coordinates, and segments across the origin end past the length

Every module shares one coordinate convention: 0-based and half-open, as in Python slicing
and Biopython's `SimpleLocation`. `Segment(start, end)` holds bases `start` to `end - 1`.
SnapGene's 1-based inclusive `"start-end"` ranges are converted where the file is read and
written, and nowhere else.

A segment across the origin of a circular record keeps `start < end` and lets `end` pass the
record's length. On a 2686 bp plasmid, a 6 bp site starting at index 2683 is
`Segment(2683, 2689)`. Length is always `end - start`, and no segment is empty.

A position a person reads, such as a map's label or hover text, is 1-based and inclusive, as
SnapGene and GenBank print it, and runs across the origin as it reads: that site is `2684 .. 3`.
It converts where the text is written, as a file format converts where it is read.

A position a person types takes the same form. `plot map --region 2684..3` names that site as
the map prints it, `END` before `START` running across the origin, and the command converts it
to `(2683, 2689)` where it reads the text. Text of that form is always a span, even where a
feature is named so. `goldengate plan --site START-END` differs: it passes its numbers on
unconverted, 0-based and half-open, and cannot run across the origin.

## Considered options

- **`end < start` for a wrapping span**, as SnapGene writes it. Every length and containment
  calculation needs a branch, and `start == end` cannot say whether it means nothing or the
  whole circle.
- **Split at the origin into two segments**, as Biopython does with a join. One contiguous
  span becomes two, so a round trip cannot tell it from a genuine two-segment feature, and a
  site or primer across the origin stops being one span.

## Consequences

Slicing `sequence` directly is wrong across the origin; `SequenceRecord.extract` is the
wrap-aware read. Edits that shift coordinates or rotate the origin reduce positions modulo the
length.

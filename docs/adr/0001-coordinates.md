---
search:
  exclude: true
---

# Half-open coordinates, and spans across the origin count past the length

Every module shares one coordinate convention: 0-based and half-open, as in Python slicing
and Biopython's `SimpleLocation`. `Segment(start, end)` holds bases `start` to `end - 1`.
SnapGene's 1-based inclusive `"start-end"` ranges are converted where the file is read and
written, and nowhere else.

A segment across the origin of a circular record keeps `start < end` and lets `end` pass the
record's length. On a 2686 bp plasmid, a 6 bp site starting at index 2683 is
`Segment(2683, 2689)`. Length is always `end - start`, and no segment is empty.

A feature extends the idiom: each segment after the origin starts past the length. On a 100 bp
record, one reading 91..98 then 2..8 is `(90, 98), (101, 108)`. A circular record refuses
descending starts and a feature longer than one turn. GenBank and SnapGene number such a
segment from the first base again, converted where the file is read and written.

A position a person reads, such as a map's label, is 1-based and inclusive, as SnapGene and
GenBank print it, and runs across the origin as it reads: that site is `2684 .. 3`. It converts
where the text is written.

A position a person types takes the same form. `plot map --region 2684..3` names that site,
converted to `(2683, 2689)` where the text is read. Text of that form is always a span, even
where a feature is named so. `cloning goldengate plan --site START-END` passes its numbers on
unconverted, 0-based and half-open, and cannot run across the origin.

## Considered options

- **`end < start` across the origin**, as SnapGene writes it. Every length and containment
  calculation needs a branch, and `start == end` cannot say whether it means nothing or the
  whole circle.
- **Split at the origin into two segments**, as Biopython does with a join. A round trip
  cannot tell one span from a genuine two-segment feature, and a site or primer stops being
  one span.
- **Later segments numbered within the record, guarded by a validator.** Sorted, a feature
  across the origin is a well-formed one that does not cross it, leaving a validator no
  evidence to check.

## Consequences

Slicing `sequence` is wrong across the origin; `SequenceRecord.extract` reads across it. Edits that
shift coordinates or rotate the origin reduce positions modulo the length, then count features round
with `Feature.counted_round`, and positions read in order, such as junctions, with `counted_round`.
A feature lists its segments in reading order; on a circular record the representation carries that
rule, not prose: starts ascend, so a sort changes nothing. A feature cut apart at a linear record's
two ends still lists the higher start first, which a sort breaks.

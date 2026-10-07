---
search:
  exclude: true
---

# A plate is a seating, a number attaches to what changes it, and a hole is never a value

**A plate is a seating.** A reaction already says how many vessels it runs and what differs
between them; a `Plate` says only where each sits, so the two cannot disagree. A `Vessel` is the
general thing; a `Plate` is one whose positions form an array. Format is one parameter, the well
count, reaching 1536 because four 384-well plates compress into one. A well's name resolves
against the protocol's own items, so a dangling one is reported. One `Transfer` covers every
move, from one to one to a compaction. Plates sit on the protocol, outliving the step that
fills one.

**Provenance is per row.** A `Source` lives once in the protocol's map, and a `Citation` hangs
on each component, incubation, material and bill row. Per scalar would rewrite every number as
a value object.

**A number attaches to the thing that changes it.** NEB #C3020 and Endura disagree on every
electroporation setting, because the program belongs to the cells. A material's parameters are
therefore keyed by catalogue number in `bench`: swap the cells and the number changes while the
step does not.

**A rule is `forbids` or `requires` on the material.** It travels with the material into every
step using it, and a `when` matches the tube's contents, so *never add PEG to the T7 ligase
reaction* cannot be edited out of a step that never held it. A rule computing a number stays a
function.

**A hole is not judged and names its ticket.** The field stays empty and the hole stands beside
it, so no loader reads a union and no reader sees a guess. A price hole names no issue: a
missing input is not a defect.

**The bill.** The price record is the user's, read by `bench.prices`, never shipped. Quantities
compute with no record, and the money is then a hole. Headroom is reported beside each banded
row, and **price steers no design**: one changing with whose list was loaded would not be
reproducible between labs.

## A band is a cell, not a column

One `bands` cell holds them all, semicolons between: `count 1-100; length_nt 1-200`. Columns per
quantity would widen the record to every quantity any row bands, and re-shape it the first time
a vendor bands by a new one. The cell keeps the shape fixed for one grammar: a name, a low, a
high, both ends inclusive, and an empty high for a tier with no top, so nobody writes a
sentinel.

## Consequences

`Protocol.audit` judges the protocol's own citations, wells, rules and holes. Nothing is added
to `pixi run check`, which everyone afterwards pays for.

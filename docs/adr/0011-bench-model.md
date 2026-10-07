---
search:
  exclude: true
---

# A plate is a seating, a number attaches to what changes it, and a hole is never a value

**A plate is a seating.** A reaction already says how many vessels it runs and what differs
between them; a `Plate` says only where each sits, so the two cannot disagree. A `Vessel` is the
general thing — a tube, a BioAssay plate — and a `Plate` is a vessel whose positions form an
array. Format is one parameter, the well count, reaching 1536 because four 384-well
plates compress into one. A well holds a name resolving against the protocol's own materials,
oligos, vessels and plates, so a dangling one is reported. One `Transfer` covers every move: one
to one is acoustic, many to one is a pool, and a dense re-layout dropping the failed wells is a
compaction. A step references a well by holding the transfer, and plates sit on the protocol,
outliving the step that fills one.

**Provenance is per row.** A `Source` lives once in the protocol's map, and a `Citation` of
source plus locator hangs on each component, incubation, material and bill row. Per scalar would
rewrite every number as a value object.

**A number attaches to the thing that changes it.** NEB #C3020 pulses at 2.0 kV, 200 Ω and
25 µF; Endura at 1800 V, 600 Ω and 10 µF, disagreeing on every setting because the program
belongs to the cells. A material's parameters are therefore keyed by catalogue number in
`bench`: swap the cells and the number changes while the step does not.

**A rule is `forbids` or `requires` on the material.** It travels with the material into every
step using it, and a `when` matches the tube's contents, so *never add PEG to the T7 ligase
reaction* and *do not heat-inactivate it in PEG* cannot be edited out of a step that never held
them. A rule computing a number stays a function.

**A hole is not judged and names its ticket.** The field stays empty and the hole stands beside
it, so no loader reads a union and no reader sees a guess. A price hole names no issue: it is a
missing input, not a defect in what we know.

**The bill.** The price record is the user's, read by `bench.prices`, never shipped. Quantities
compute with no record loaded, and the money is then a hole. Headroom is reported beside each
banded row, and **price steers no design**: one changing with whose list was loaded would not be
reproducible between labs.

## Consequences

`Protocol.audit` judges the protocol's own citations, wells, rules and holes. Nothing is added
to `pixi run check`, which everyone afterwards pays for.

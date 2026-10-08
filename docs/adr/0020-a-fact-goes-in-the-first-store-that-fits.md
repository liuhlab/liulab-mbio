---
search:
  exclude: true
---

# A fact goes in the first of four stores that fits

Four stores for a fact exist and each has a working precedent, but nothing said which one a new
fact goes in, so each author chose, and the fact ended up where its builder sat. `src/` carries
846 `#:` source lines, 37 in `igga/steps.py` and 56 in `dmx.py`: the choice is made often, and
made wrongly the fact dies with the next edit to the step holding it.

Ask these four in order. The first yes is the store. Four noes, and the next section has it.

1. **Does a licence forbid shipping it, or is it the lab's own and going stale?** A file the user
   holds, read at run time. Where the file is absent a `Hole` stands in its place, never an
   estimate (`ligase.py`, `bench/prices.py`).
2. **Does a script rebuild it from a tracked input?** `src/liulab_mbio/data/`, built by a script
   in `scripts/` and sourced in `docs/research/` (`enzymes.json`, `ligation_fidelity.json`).
3. **Does it belong to a thing the bench names by catalogue number, and hold wherever that thing
   is used?** A table in `bench/materials.py` keyed by that number, read through `material()`:
   what the material `contains`, the `rules` it travels with, the electroporation settings its
   cells decide. It then reaches every step naming the material, so no step edit drops it, and a
   caution is the next fact to land here.
4. **Is it true of this run only, or does it interpolate a computed number?** The protocol JSON
   the builder writes (ADR 0002), where prose true of one run belongs (ADR 0019).

This extends ADR 0011's rule that a number attaches to the thing that changes it. Question 3 is
that rule, third because a licence and a rebuild script settle a fact's file before its subject
does.

## ADR 0010 asks a different question, and the two get confused

ADR 0010 settles whose number it is — the method's or the project's. This one settles where a
fact is stored. A number the method chooses answers no question above: it is a `#:`-sourced
constant in the method's own bench module, in the shape of `igga/bench.py`, whose 19 `#:` lines
sit beside the four functions that spend them. ADR 0010 is not reopened.

## Considered options

- **A `Note` type with a `kind` set** — `forbids`, `requires`, `caution`, `expects`, `hole`, with
  a resolver reading what each subject carries. Nothing measured repeats on `expects`, a hole is
  an absent value and keeps its own type, and the rename changes the `rules` JSON key and both
  construction sites without removing one duplicate.
- **An addressed `knowhow.toml` sidecar**, resolved into a protocol as it renders. It wants an
  address grammar and a step key the model does not have, and ADR 0019 already declined a
  hand-edited file of prose that de-duplicates nothing.
- **A shared module constant.** It removes the copies in the source, but an edit to the protocol
  JSON still deletes the fact.

## Consequences

Nothing is built and nothing is gated: all four stores already ship. An author with a new fact
answers four questions rather than choosing, and a reviewer reads the same four back.

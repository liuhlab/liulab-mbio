---
search:
  exclude: true
---

# A fact goes in the first of four stores that fits

Four stores for a fact exist and each has a working precedent, but nothing said which one a new
fact goes in, so it landed in the nearest builder function: `src/` carries 846 `#:` source lines,
37 of them in `igga/steps.py` and 56 in `dmx.py`. A fact written into a builder is a fact the next
edit to that builder deletes.

Ask these four in order. The first yes is the store.

1. **Does a licence forbid shipping it, or is it the lab's own and going stale?** A file the user
   holds, read at run time. Where the file is absent a `Hole` stands in its place, never an
   estimate (`ligase.py`, `bench/prices.py`).
2. **Does a script rebuild it from a tracked input?** `src/liulab_mbio/data/`, built by a script
   in `scripts/` and sourced in `docs/research/` (`enzymes.json`, `ligation_fidelity.json`).
3. **Does it belong to a thing the bench names by catalogue number, and hold wherever that thing
   is used?** A table in `bench/materials.py` keyed by that number, read through `material()` and
   carried on the `Material` as `contains`, `rules` or `cautions`. It reaches every step that
   names the material, so no step edit drops it.
4. **Is it true of this run only, or does it interpolate a computed number?** The protocol JSON,
   written by the builder (ADR 0019).

This extends ADR 0011's rule that a number attaches to the thing that changes it. Question 3 is
that rule; the other three say where a fact goes when no material owns it.

## ADR 0010 asks a different question, and the two get confused

ADR 0010 settles whose number it is — the method's or the project's. This one settles where a
fact is stored. A number the method chooses is on neither list above: it is a `#:`-sourced
constant in the method's own bench module, in the shape of `igga/bench.py`, whose 16 constants
carry 19 `#:` lines across 4 functions. ADR 0010 is not reopened.

## Considered options

- **A `Note` type with a `kind` set** — `forbids`, `requires`, `caution`, `expects`, `hole`, with
  a resolver reading what each subject carries. Nothing measured repeats on `expects`, a hole
  keeps its own type under ADR 0011, and the rename changes the `rules` JSON key and both
  construction sites without removing one duplicate.
- **An addressed `knowhow.toml` sidecar**, resolved into a protocol as it renders. It wants an
  address grammar and a step key the model does not have, and it is a hand-edited data file,
  which `CLAUDE.md` and ADR 0019 both refuse.
- **A shared module constant.** It removes the copies in the source, but an edit to the protocol
  JSON still deletes the fact, which is the failure being fixed.

## Consequences

Nothing is built and nothing is gated: all four stores already ship. An author with a new fact
answers four questions instead of reaching for the nearest builder, and a reviewer reads the same
four back.

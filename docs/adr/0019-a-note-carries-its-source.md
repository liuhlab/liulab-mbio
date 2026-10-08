---
search:
  exclude: true
---

# A troubleshooting answer carries its source, or it is deleted

`Troubleshooting` held a problem and a solution and nothing saying where the solution came from,
so a reader at the bench could not tell advice a vendor published from advice someone remembered.
It now carries an optional `citation`, the `Citation` every other row already carries. An entry
no document supports is deleted rather than kept with the field empty.

## No library ships

A set of notes keyed by problem, shipped once and reused by every method, is refused on
measurement. Across all 69 call sites there is not one duplicate `(problem, solution)` pair.
Eight problems repeat — "No band", "Few colonies later", a reaction that outgrows its volume —
and each carries a different solution, because a solution names that method's enzyme, template or
vessel. Five entries interpolate a number the run computed, which no static file holds. A library
would remove no duplication and could not carry the entries that matter most.

`CLAUDE.md` also holds package data to a script in `scripts/` rebuilding it from a tracked input.
Prose has no input to rebuild from, so a note library would be the first hand-edited file in a
directory whose rule is never to hand-edit one.

## Nothing new is built, and nothing new is gated

`Source`, `Citation`, `Protocol.sources` and `Protocol.cited` already carry provenance per row,
and `citing` already drops the sources a run never cited. One line in `cited` reads the
troubleshooting entries, so the existing `sources` audit reports a key naming no source — the
only failure a machine can see. The page renders the citation beside the solution as an anchor
into its own Sources section, so the reader never leaves the page to follow one.

There is no audit check counting uncited entries. "Nobody sourced this" is a defect in the
author, not in the protocol's data, and once the entries are swept the count is zero for good, so
the check would cost every protocol something and report nothing. If uncited entries return, that
check is the escalation.

A figure cites the same map through the same `Citation`. The artefact differs — a drawing ships
as package data a script rebuilds, a note does not — but provenance is one mechanism.

## Considered options

- **A shipped library of notes.** Refused above: nothing to de-duplicate, and the entries
  carrying a computed number could not live in it.
- **A required citation.** Every call site would gain a field before anyone had read a source
  for it, and the cheapest way through is a citation that does not support what it sits on.
- **A fifth check on `Protocol.audit`.** A permanently silent check is paid by every protocol
  that renders.

## Consequences

The sweep of the existing entries is author work, one entry at a time: cite it from what
`reference_docs/` holds, or delete it. A step whose entries all go shows no troubleshooting
section, which already renders correctly.

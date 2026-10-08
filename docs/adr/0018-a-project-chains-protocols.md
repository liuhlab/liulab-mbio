---
search:
  exclude: true
---

# A protocol declares what it consumes and produces, and a project chains them by name

A run is several protocols, each a document someone follows in one sitting. What joins them is
what the bench carries from one to the next, so a protocol declares two things: `consumes` and
`produces`, each a tuple of an `Item`. A `Project` runs protocols in order, and `Project.audit`
returns a `handoffs` check: every consumed name is a project input or an earlier protocol's
output.

**The name is the contract.** A project matches by name, as a `Plate.seating` name already
resolves against the protocol's own items. `Item.what`, `spec` and `storage` are prose for the
bench, and nothing parses them.

**A mismatch is a badge.** The chain is judged like any other check, never refused: a run with a
dangling input is still a document someone can read, and the page says which name nothing hands
over.

**The project holds what no single protocol owns.** The bill is the run's, since two protocols
buying the same cells would be counted twice. Checks that judge the design rather than the bench
are the project's; a protocol rendered alone keeps its own. The explanation of why the run is
shaped as it is becomes `background`, a `Topic` at a time, so a step never stops to explain one.
Materials, equipment, references, sources and holes stay on the protocol that owns them.

`Step.section` is a label on the step, not a container: the steps stay one list and the numbering
runs through it, so `step-14` is one place in one document.

## Considered options

- **`requires`, beside what a protocol consumes.** Reagents are `materials`, hardware is
  `equipment`, and what the protocol itself promises is `Protocol.audit` and `checks`. A third
  list would restate one of those.
- **A `config` field.** A build's choices are the project file, and a rendered protocol has every
  number already resolved. ADR 0010 names a method-as-configuration-language the thing to guard
  against, and this would be one.
- **Typed inputs — a plasmid, a plate, a tube.** A type nothing reads is a word with a check
  behind it. A name and a sentence say more to the reader and bind no later protocol.
- **Matching a spec against the producer's**, so `≥100 ng/µL` is read and compared. Parsing a
  quantity out of prose invents a grammar for every unit the bench writes. The producer states
  what it makes, the consumer restates what it needs, and the reader compares them.
- **A declared item list in place of `Rule`.** A rule governs what is in a tube — its `when`
  matches the tube's contents — not what passes between protocols. The two do not overlap, and
  the T7 ligase rules are the only thing keeping PEG out of a reaction.

## Consequences

Every one-protocol plan writes a project of one, through `cloning.plan.as_project`, so whatever
reads a run reads one shape. How a project is written and rendered as a folder of pages is its
own decision, and this one fixes only what a protocol declares.

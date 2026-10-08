---
search:
  exclude: true
---

# A protocol declares what it consumes and produces, and a project chains them by name

A run is several protocols, each a document someone follows in one sitting. What joins them is
what the bench carries from one to the next, so a protocol declares `consumes` and
`produces`, tuples of an `Item`. A `Project` runs protocols a place at a time, and
`Project.audit` returns a `handoffs` check: every consumed name is a project input, or a name
every way of an earlier place produced.

**A protocol may name the job it is one way of doing.** Protocols sharing that `choice` are the
ways: one place in the run, and the bench does one of them. So `handoffs` reads
across the branch, not along it. A second check, `choices`, fails a job one protocol names alone,
or ways leaving different things, and warns where nothing says how to pick. Nothing picks a way;
what to weigh is the `background` topic the job titles.

**The name is the contract.** A project matches by name, as a `Plate.seating` name already
resolves against the protocol's own items. `Item.what`, `spec` and `storage` are prose nothing parses.

**A mismatch is a badge**, never a refusal: a run with a dangling input is still a document
someone can read, and the page says which name nothing hands over.

**The project holds what no single protocol owns.** The bill is the run's: two protocols buying
the same cells would be counted twice. Checks judging the design are the project's; a protocol
rendered alone keeps its own. Why the run is shaped as it is becomes `background`, a `Topic` at a
time. Materials, equipment, references and holes stay on the protocol that owns them.

`Step.section` labels a step, not a container: the steps stay one list and the numbering runs
through it.

## Considered options

- **`requires`, beside what a protocol consumes.** Reagents are `materials`, hardware is
  `equipment`, and what the protocol promises is `Protocol.audit` and `checks`.
- **A `config` field.** A build's choices are the project file, and a rendered protocol has every
  number resolved; ADR 0010 guards against method-as-configuration-language.
- **Typed inputs — a plasmid, a plate, a tube.** A type nothing reads is a word with a check
  behind it; a name and a sentence say more.
- **Matching a spec against the producer's.** Parsing a quantity out of prose invents a grammar
  for every unit the bench writes; each side states its own and the reader compares.
- **A declared item list in place of `Rule`.** A rule governs what is in a tube, not what passes
  between protocols; the two do not overlap.
- **A condition on a protocol, or a project holding a graph.** Either makes the run a program,
  ADR 0010's own warning. Guidance held on the ways went with them: two copies of one comparison
  can disagree.

## Consequences

Every one-protocol plan writes a project of one, through `cloning.plan.as_project`, so whatever
reads a run reads one shape.

A project is written as a folder: `project.json`, one page per protocol named by its place in
the chain, an index, and the two shared pages. Routing inside a page cannot work on disk, so each
page is self-contained, links by relative address, and names its neighbours in a printed line.

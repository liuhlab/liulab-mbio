---
search:
  exclude: true
---

# The docs are arranged by what a reader is doing, and the reference is a hand-kept list

The site was a tree of what the package holds. It becomes five sections named for what an
experimental biologist came to do — Home, Methods, Protocols, Projects, Reference — and the
import and CLI names become `mbio` and `synbio`, the distribution staying `liulab-mbio`. The page
list is #490, the section specs #488 to #496, the map #487.

The axis costs the other reader. Someone who thinks in packages has no branch to walk down and no
page per module; the source is their reference, CLAUDE.md's architecture table the only
tree. That reader is a contributor, agent-facing anyway.

## The reference is 61 symbols, and no gate keeps the list honest

A module directive documents what a package holds, not what anyone calls: `api.md`'s 132
rendered over a thousand headings nobody read. `reference/python.md` names 61 symbols, grouped
by the reader's task.

So a symbol someone wants may be documented nowhere, and the list goes stale the first time a
signature changes and nobody opens the page. **No conformance rule polices it, deliberately.** A
rule would have to encode which symbols are user-facing, and that judgement *is* the page: it
would either restate the list or demand back everything the cut removed. A list someone fixes on
noticing is cheaper than a gate everyone pays for.

`reference/cli.md` is the exception: a test reads the apps' `registered_commands` and
`registered_groups`, not `--help`, and asserts the page's verbs match. A verb coming or going is
a fact the app knows about itself, not a judgement.

## A human page restates what an agent file says

`docs/adr/`, `docs/agents/` and `docs/research/` are out of the nav and excluded from search, so
a link into one is a dead end for its intended reader. A published page restates it in
its own words instead — two copies that can drift with nothing comparing them, taken because the
link fails today and the drift only might.

## Considered options

- **Nav by package.** The shape of the repo, not of a job; every tool in #489's survey that lost
  this reader had it.
- **Module directives, filtered.** `__all__` or a private-name rule still answers "what is in
  this module", not the reference page's question.
- **Linking the research notes and the ADRs.** Publishing them puts uncapped,
  jargon-exempt prose in front of a bench reader.
- **Sharding `examples-check`.** Measured worth real time, and deferred: the serial run sits
  inside the docs job's budget. Reopen above thirty seconds warm.
- **A shim for the old names.** A second spelling to keep correct, for no
  downstream this repo cannot edit.

## Consequences

Every page showing what a command printed is backed by a generated, committed corpus, so no page
can claim a number the package stopped printing; six show no output and carry none. The axis
changed the package, not the reverse: barcodes and primers gain a design verb and codon-optimize
an `--out`, three Methods pages otherwise showing no command that writes a file. The rename lands
as one pull request, no alias period. Its trap: conformance rule 14 reads its import
names from `[tool.liulab.must-not-import]` and goes silently vacuous, not red, if the key and the
directory move apart. #494 holds the sweep and its guards.

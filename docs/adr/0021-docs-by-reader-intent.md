---
search:
  exclude: true
---

# The docs are arranged by what a reader is doing, and the reference is a hand-kept list

The site was a tree of what the package holds. It becomes five sections named for what an
experimental biologist came to do — Home, Methods, Protocols, Projects, Reference — and the
import and CLI names become `mbio` and `synbio`, the distribution staying `liulab-mbio`. The page
list is #490, the section specs #488 to #496, and the map they hang from is #487. None of that is
restated here.

The axis costs the other reader. Someone who thinks in packages has no branch to walk down and no
page indexing a module; the source is their reference, and CLAUDE.md's architecture table is the
only tree. That reader is a contributor, who is agent-facing anyway.

## The reference is 59 symbols, and no gate keeps it honest

A module directive documents what a package holds, and nothing narrows it to what anyone calls:
the 96 on `api.md` rendered over a thousand headings nobody read. `reference/python.md` names 59
symbols one at a time, against 1,191 public top-level names, grouped by the task the reader has.

So a symbol someone wants may be documented nowhere, and the hand-kept list goes stale the first
time a signature changes and nobody opens the page. **No conformance rule polices it,
deliberately.** A rule would have to encode which symbols are user-facing, and that judgement
*is* the page, so it either restates the list or demands back everything the cut removed. A list
someone fixes on noticing is cheaper than a gate everyone pays for.

## A human page restates what an agent file says

`docs/adr/`, `docs/agents/` and `docs/research/` are out of the nav and excluded from search, so
a link into one is a dead end for the reader it was written for. A published page restates the
content in its own words instead. The cost is two copies that can drift, with nothing comparing
them — taken because the link fails for that reader today and the drift only might.

## Considered options

- **Nav by package**, one branch each. The shape of the repo rather than of a job, and every tool
  in #489's survey that lost this reader had it.
- **Module directives, filtered.** `__all__` or a private-name rule still answers "what is in
  this module", which is not the question a reference page is asked.
- **Linking the research notes and the ADRs.** Publishing them to make the links work puts
  uncapped, jargon-exempt prose in front of a bench reader.
- **Sharding `examples-check`.** Measured worth real time, and deferred: the serial run sits
  inside the docs job's budget. Reopen when the warm run passes thirty seconds.
- **A shim for the old import and CLI names.** A second spelling to keep correct, for no
  downstream this repo cannot edit.

## Consequences

Every page showing what a command printed is backed by a generated, committed corpus, so no page
can claim a number the package stopped printing; seven pages show no output and carry none. The
rename lands as one pull request with no alias period. Its trap: conformance rule 14 reads its
import names from `[tool.liulab.must-not-import]` and goes silently vacuous rather than red if
the key and the directory move apart. #494 holds the sweep and its guards.

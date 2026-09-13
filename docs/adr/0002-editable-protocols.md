---
search:
  exclude: true
---

# A pipeline hands back editable protocol data, and the agent invents no number

The package supplies building blocks and a few whole pipelines; the agent combines them into
what the user asks for. So a pipeline writes its protocol as data: `Plan.write` puts
`protocol.json` in the directory and renders `protocol.html` from it. Asked for something the
pipeline did not anticipate — skip a step, add one, reword one — the agent edits that file and
renders it again, as `build-protocol` says.

Two rules hold the split. Numbers and verdicts come only from the package: the agent never
produces from its own knowledge one the package computes or checks, and where the package
cannot compute what a request needs, it says so instead. Changes to what a pipeline computes —
the enzyme, the polymerase, the inserts, their orientation — go back through the pipeline, so
everything depending on them is computed again.

The loader refuses a value of the wrong JSON type and names where it is, so an editing mistake
fails when the page is rendered rather than printing on the page.

## Considered options

- **A flag per step on every pipeline**, such as one that skips the DpnI digest. Each
  unanticipated request would add an option, so the options grow without bound while the
  pipeline learns nothing general.
- **A change log or an edit history in the protocol.** Each run is one session's job, so a
  record of the edits is another format to write and keep true, for a reader who has the final
  page in front of them.
- **Keeping the pipeline's own copy beside the final protocol.** Two files that differ and
  neither named as the current one, so the user has to work out which is theirs. Keeping copies
  is the user's call.

## Consequences

The user's folder holds one final `protocol.json` and one `protocol.html`. Dropping a step
leaves stale mentions — the materials, another step's troubleshooting, an oligo row — which no
code checks, so `build-protocol` tells the agent to fix them. An assembly plan's `primers.tsv`
stays the whole order sheet, which is the record of what was designed. A pipeline gains an
option only when the same request keeps coming back.

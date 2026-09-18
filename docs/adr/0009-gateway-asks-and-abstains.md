---
search:
  exclude: true
---

# Where a record or a source runs out, the plan asks or abstains — it never infers

Two Gateway checks needed something neither the user's records nor the manuals supply. Both are
answered the same way, and the way is worth writing down because inferring would have been
easier and would have looked right.

**A fusion is an input.** The reading frame across an att junction is judged only at an end
`--fusion` names. Nothing in a plasmid map says which coding sequence is a tag: the one nearest
the junction is as likely to be the resistance gene, so detection would fail correct work, which
`AGENTS.md` counts as evidence against a rule. The vendor states both frame rules as intentions
— *if you wish to fuse your product to a tag* — so an input matches the rule as written.

**The ccdB-resistant-host rule carries no verdict.** The manuals require that strain for growing
a donor or destination vector, and forbid an F′ strain, whose `ccdA` cancels the selection. Both
of those are judged. Whether a reaction plated on a resistant strain still counter-selects is
stated nowhere, and the vendor's support page leans the other way, so the check reaches the page
as guidance with `status=None` and reads "not judged".

## Considered options

- **Detect the fusion from the record**, by reading the nearest coding sequence across each
  junction. It fires on a vector with no tag, where that sequence is the marker, and the cost
  lands on a user who did everything right.
- **Pass the ccdB host check.** A pass reads as sourced, and nothing sources it.
- **Fail it.** No manual forbids what it would fail.
- **Leave it off the page.** The guidance is what a bench most needs here, and dropping it reads
  as though nobody asked.

## Consequences

A user who forgets `--fusion` gets no frame verdict rather than a wrong one, so the option has
to be named in the method file and on the page. One badge on every Gateway protocol reads "not
judged"; a reader has to be told that means unmeasured and never fine. Both are cheap to reverse
the day a source states either rule.

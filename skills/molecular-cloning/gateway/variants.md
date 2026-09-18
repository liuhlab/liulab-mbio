# Gateway variants and kit quirks

Read this only when `METHOD.md` has sent you here. Each one below is a published route this
command does not plan, and none of their numbers may be borrowed by a plan that ran the staged
route.

## The one-tube BP and LR

The vendor publishes a route that runs both reactions in one tube, on its own timings, with one
stop and one transformation. It is not the staged protocol with a step dropped: it has its own
incubations and its own reaction sizes, and it costs yield and certainty. The manual puts it at
a fraction of the expression clones the staged route gives, and tells you to sequence every
clone it does give.

So a user asking to skip the miniprep is asking for a different protocol. Say that, and plan the
staged route, which is what the numbers on the page belong to.

## The two-step adapter PCR

Above about 70 bases of primer the manual switches to a two-step PCR: a first round with short
gene-specific primers, then a second with universal attB adapter primers that finish the tail.
The plan prints this as a note on the amplify step rather than designing it.

Two things people get backwards. The short figure the sources give is a **gene-specific overlap**
in that second route, not a shorter tail — the product still carries the whole attB tail either
way. And it is conditional on the two primers' stoichiometry, not a flat number. Do not shorten
a one-step tail to it.

## Entry vectors whose own sites drift

A vendor entry vector can carry a mutated att site. That is why the plan searches for sites
instead of matching a shipped sequence, and why it designs its own sequencing primers rather
than naming the vendor's: `GW1` and `GW2` suit one kit's vector alone, and an M13 primer crosses
well over 100 bases of vector before it reaches the insert.

## MultiSite Gateway

More than one fragment, and every att site beyond attB1, attB2 and their partners. It stays on
the chooser's unsupported list. Say so and stop, rather than planning it out of this method's
arithmetic.

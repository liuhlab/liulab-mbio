---
search:
  exclude: true
---

# The gate reads reactions, and mbio holds the predicates it calls

`liulab_mbio.reaction` holds one tube: the molecules in it, each bound to a role, the enzymes
acting there, and how many vessels differ only in which species fills a role. It carries no
verdict. `liulab_synbio.library.gate` composes this method's reactions from a finished design and
judges each one, calling `liulab_mbio`'s predicates — `sites.find_sites`, `sites.digest`,
`translate.in_frame`, `translate.stop_codons`, `barcodes.separation`, `barcodes.check_barcodes`,
`overhangs.fidelity` — with this method's own parameters, and reporting in this method's words.

A predicate returns a finding, and the gate builds the `Check`. A finding is whatever domain
object mbio already returns, so there is no finding type of its own: `sites.CutSite` carries the
enzyme, the span, the strand and both cuts, and a stop codon is a `Segment`. A `Judgement` pairs
the check with the findings behind it, so a map can draw a failure on the record it occurred in.
`liulab_mbio.checks.Check` does not change, so the cloning pipelines do not move.

A failing check names what is wrong and need not name a remedy. Requiring one is a cost paid by
everyone who writes a check afterwards; a hole stays a hole and goes to an issue.

## The pipeline judges nothing of its own

`plan_library` reports the gate's verdict, and `rounds.py`'s three checks are retired. Two were
already bounded, by `donor releases` with `destination_vector`, and by `terminal block` with
`terminal stop`. The third found a gap: no round opens the last product, so `check_product`
gains `product opens`, the two-cut rule every destination is held to. What a round also checked,
the stuffer's exact span, is a plan's coordinate, which a gate may not read.

## Considered options

- **The gate reads bare records.** Nothing then says which tube a failure belongs to, and a
  molecule clean in one tube and fatal in another — the destination's own sites, inert where no
  enzyme of that tube reads them — has no place to be either.
- **The gate reads the plan.** It would then judge the designer's reasoning rather than the
  molecules, and an agent's design could not be judged at all.
- **A finding type of its own.** Every predicate would have to wrap what it already returns, and
  a drawing would have to unwrap it again.

## Consequences

A design an agent composed and one `plan_library` wrote are judged identically, because both
arrive as finished records bound to roles. What must be true inside one tube is the method's, so
nothing in `liulab_mbio` knows one. The regression corpus is the proof: `docs/examples/ap1-library`
passes, and one deliberate break per rule fails on the check that should catch it.

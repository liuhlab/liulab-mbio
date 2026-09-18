---
search:
  exclude: true
---

# A library is assembled in rounds, and `cloning/goldengate/` does not generalise to one

A combinatorial library is built one round at a time. Each round is its own digest, ligation,
transformation and prep: the round's part list is cut externally to release each part, the library
built so far is cut internally to open it, the two ligate, and the product becomes the next round's
destination.

One pot here is not a shortcut but a different mechanism. The product deliberately keeps the
internal enzyme's sites, and that is what lets the next round open it, so one reaction holding
every enzyme would cut the product it just made. `docs/research/protein-library-assembly.md`
carries the scheme, its enzymes and the frame arithmetic.

This needs recording because `cloning/goldengate/` sits next door and a reader will assume it
generalises. It does not. `Part` gives none of its fields a default, so `forward`, `reverse` and `report` are
all required, and a synthesised part — which has no primers and no pair report — cannot be one.
Nor will `_chain` build a product unless exactly one part begins with the overhang the part before
ends with, where a part list is many parts sharing one entry overhang by design. So `library/` is
a second pipeline beside it, and the round boundary shapes everything above it rather than being a
flag to flip later.

## Considered options

- **One pot, every enzyme in one reaction.** It cuts its own product, for the reason above. A
  method that avoided that would be a different one with its own failure modes, which #55 puts
  out of scope.
- **A different enzyme pair each round**, so no product need keep a site. Every round then needs
  its own pair and its own overhang standard, and the rounds are capped by the enzymes free in
  every part.
- **Reusing `cloning.goldengate.assemble` once per round.** A part list is many parts sharing
  one entry overhang and carrying no primers, which its chain rule refuses and its `Part` cannot hold, so
  admitting one means loosening both for the single-reaction pipeline as well.

## Consequences

Rounds are what the method costs: a transformation, a growth and a prep each, and each round's
loss carries into the next, so library coverage is counted per round and not once at the end. A
plan records every round's intermediate product, and the protocol covers every round. The product
reports cut sites for the scheme's own enzymes, which is the design working rather than a defect.
`library/` imports no cloning method: the overhang rules both need sit below every pipeline, per
`docs/adr/0007-cloning-methods.md`.

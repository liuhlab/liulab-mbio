---
search:
  exclude: true
---

# A record knows which tube it is digested in

`synbio.igga.vector` held one barred set: a destination had to be free of every enzyme
the method names, outside the piece it gives up. That refused the method's own vector. A round's
destination is a plasmid, transformed and prepped every round, and its backbone is what frees the
cargo at the end — so it carries the external enzyme's sites, and a blunt chopper's outboard of
those, on purpose. The shipped predicate named one of them and turned the vector away.

The barred set is now the enzymes of **one** digest. A `Cassette` says which piece, to which
enzyme, and in which tube, and a refusal names that tube. A round's destination bars the internal
enzyme and the chopper that shreds the stuffer it gives up. The working vector bars the cargo
enzyme and what the method reserves. A donor backbone bars nothing.

## A donor bars nothing, because it keeps nothing

"No site outside the cassette" assumes the record keeps its backbone. Two tubes do. The donor's
does not: the digest frees the cargo and throws the backbone away, so the releasing cuts lie
outboard of the cassette by construction, and so does the chopper aimed at what they leave.
Measured on the rebuilt DMX vector, barring either refuses the molecule the method is built on.
A third releasing site is caught where it does harm — the gate counts the cuts in the tube — and
where a chopper may sit is `gate.check_dmx_vector`'s question.

One stuffer, two tubes: the same 34 bases are what the internal enzyme frees from a round's
destination and what the external enzyme frees from a donor. The cassettes differ only in the
enzyme and the tube.

## Considered options

- **Keep one barred set and relax it by enzyme.** The set would then mean nothing in particular,
  and the next method would inherit a list with no rule behind it.
- **Keep the one-tube, linear-donor model and say so on the method page.** It answers the donor
  side and leaves the destination side broken: nothing could then release the finished cargo.
- **Give `Cassette` the tube and bar every enzyme of it.** Closest to right, and wrong for the
  donor, whose own enzymes read outboard by design.

## Consequences

`vector.round_cassette` stops barring the external enzyme and the donor's chopper, and
`vector.donor_cassette` is its sibling. `gate._cargo` takes the released piece by the end it
leaves on — the cloning scar, which every position shares — rather than by being the only piece,
because a circular donor gives up two. A rebuilt DMX destination is accepted with its outboard
sites, and `vector.released_cargo` starts answering on a library built in one.

A method enzyme is now judged where it acts and nowhere else, so a record clean in one tube and
fatal in another has a place to be both.

---
search:
  exclude: true
---

# A well may be labelled as well as seated, and only the seating has to resolve

`Protocol.audit` resolves every seated well against a declared material, oligo, vessel or plate,
and reports one that does not. Three plates could not meet that and say what they know. A picked
plate's wells carry `quarter 3`, which groups wells rather than naming an occupant. A Route B
index plate's carry `forward 2, reverse 2`, the well's printed address, whose occupant is
already declared as one material. A PCR2 plate's carry a block's name, and a block is none of
the four kinds. Declaring the first two made `audit` report 384 dangling entries.

So `Plate` gains `labels`, a second well-level channel beside `seating`: well to text, drawn on
the layout and resolving against nothing. `seating` keeps its contract word for word, and the
check keeps the one defect it catches — a well naming a material nobody declared. The three
plates move their strings across and seat nothing; `bench.plates.plate` takes `labels=` and
`seat` writes either; the renderer draws both, the seating winning a well that carries both.

The model already had this distinction at plate level, in `Plate.holds`: what a plate's wells
hold where they hold nothing named. A label is `holds` said a well at a time.

With the validation plates now declarable, `igga.steps.protocol` passes the picked plates and
Route B's index plates to the protocol, which it never did, so every well those transfers name
belongs to a plate the page draws.

## Considered options

- **Loosen the check to accept any string.** It throws away the one defect the check finds, for
  every plate in the package, to describe three. ADR 0011 put the check there deliberately.
- **Declare each block, quarter and address as a material.** Seventy-two blocks become 72
  materials restating what one order sheet already lists, which is the second count ADR 0011
  keeps a plate from being.
- **Leave the wells blank.** The model knows which block sits where. Printing nothing because
  the only channel is the wrong one is the model's failure, not the reader's.

## Consequences

Every plate declares what it knows. A reviewer reads `seating` as a reference and `labels` as
prose, and `wells` stays a real check rather than a formality. The protocol JSON grows about
2.4x on the demo, because a label does not replace a transfer: 288 moves are still 288 moves.
Nothing measures that size as a problem, only that it grew, so it stands. A whole-plate pattern
on `Transfer` would be its own decision, with its own evidence.

---
search:
  exclude: true
---

# A plate's drug is read off the vector, never quoted from the protocol it came from

Every number this method's read-back prints was quoted from the paper the step came from. The
drug cannot be, because the method changed the plasmid. Departure D11 rebuilds the DMX vector
AmpR to KanR before first use, so the last transfer moves cargo between two vectors that do not
share a marker. The source protocol plates on carbenicillin. Carbenicillin on a KanR plasmid
selects nothing, and a reader following the page as shipped reads an empty plate or a lawn.

So the marker decides. `liulab_mbio.bench.phenotype.selection_marker` reads it off the record and
`SELECTION` names its drug, exactly as the four cloning methods already do.
`liulab_synbio.igga.stages.selection_for` adds the concentration this method has a source for and
returns one phrase; `dmx.Validation` carries that phrase as `selection` and prints it wherever a
plate, a well or a broth is named. Where the record annotates no marker, the page says *the
vector's own antibiotic* and hole **H22** stands. Where it annotates one, H22 is answered and
`stages.holes_for` leaves it out.

A temperature, a colony density and a well volume are properties of the step, not of the
plasmid, so those citations stay where they were. Only the drug moves.

## Considered options

- **Quote the source's carbenicillin and note the departure beside it.** The note is prose
  beside a number, and the number is what a reader at the bench follows.
- **Make the drug a project field.** It asks the user to restate something their own file
  already spells, and lets the two disagree.
- **Name the drug with no concentration.** A plate needs one, and both this method's drugs have
  a cited concentration here: kanamycin at 50 µg/mL from the guide to its own part carrier, and
  carbenicillin at 100 µg/mL from the source protocol, on the one stage whose vector is AmpR.

## Consequences

A vector that names its marker gets a page naming the drug, the concentration and the marker it
was read from. A vector that names none gets a hole rather than a guess, which is what the demo
shows, its toy vector annotating nothing. A drug this package has no cited concentration for is
named without one rather than given a number nobody published.

Each stage now takes the drug its own vector carries: the cargo read-back and every round on the
KanR destination, the final transfer on the working vector's own. Adding a marker to
`SELECTION`, or a concentration to `stages.SELECTION_PLATE`, changes every page at once.

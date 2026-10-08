# Example: reading the AP-1 cargo back

A run of its own, over the cargo the [AP-1 library](../ap1-library/index.md) leaves on a plate.
It is not a step of that run, and nothing in that run waits for it. The library has its own
read-back pages; this one does the same job on its own, from a plate of frozen stock and a list
of names.

Every file in `protocol/` was written by the command below. Nothing here is edited by hand.

```bash
pixi run liulab_synbio dmx plan docs/examples/ap1-readback/build.json \
  --out docs/examples/ap1-readback/protocol
```

To change what these pages say, change the code and run that command again. The same inputs
write the same bytes.

## What goes in

| File | What it is |
| --- | --- |
| [build.json](build.json) | what this run chose: the sheet, the plate it spots from, the way it marks a well, and how far down the list it goes |
| [designs.tsv](designs.tsv) | 72 designs, one a row: a name, and how many pieces it was built from |

A design here is a name and a count, never a sequence. So any frozen stock a lab holds can be
read back, whoever built it and however. Two columns anyone can type are the whole of it.

`archive` names the plate each design is spotted from. `selection` is the drug that plate's own
vector carries, kanamycin here. `index_plate` is what this lab calls its own plate of barcoded
primer pairs, and naming it is what stops the pages saying that set has no source. **Replace all
three with your own before you run this.**

`route` is one of the two ways a well is marked. `validate_from` is a floor on the piece count,
and a design built from fewer pieces is left out. This run marks by index PCR and sets the floor
at 2, so the 32 designs stitched from two pieces or more are read back and the 40 made in one
piece are not. Set the floor to `0` to read every design on the sheet.

## What comes out

| File | What it is |
| --- | --- |
| [protocol/index.html](protocol/index.html) | the way in: what is read, what the bench is handed, and how long it is held |
| [protocol/01-design-read-back-index-pcr.html](protocol/01-design-read-back-index-pcr.html) | the bench page: array and pick, mark every well, call the wells |
| [protocol/reagents.html](protocol/reagents.html) | what to order |
| [protocol/references.html](protocol/references.html) | where each number was read from |
| [protocol/project.json](protocol/project.json) | the run as data, which `liulab_mbio protocol render` turns back into these pages |

This run picks 128 wells into one 384-well plate, and marks them over two 96-well plates. No DNA
is designed, so there is no sequence file and no primer sheet here.

[Highly parallel DNA synthesis and assembly](../../synthesis-and-assembly.md) explains the
method these pages follow.

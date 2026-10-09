# DMX: reading back many wells

A plate of picked colonies tells you nothing until each well is named. DMX marks every well so
that sequencing says which well a read came from, then calls each well a pass, a fail, or no
verdict at all. One command chooses the plates, the reaction and the cycling for the route you
name, and writes the bench protocol.

## When to use it

- You hold a plate of frozen clones and need to know which design is in each well.
- The wells came from an assembly that is not clean every time, so identity must be read.
- You are not running a pooled screen, which takes its identity from the screen itself.

The [iGGA library](igga.md) reads its archived parts this way; here the same job stands on
its own, over any frozen stock.

## What you need

Two files and one plate of lab stock. The build names the route, the archive plate and the
floor; the sheet lists the designs. A design here is **a name and a fragment count, never a
sequence**, so any frozen stock a lab holds can be read back, whoever built it. Two columns
anyone can type. The run below uses files that ship with the repo.

| Input | What it is |
| --- | --- |
| [build.json](../examples/ap1-readback/build.json) | the sheet, the plate it spots from, the drug that plate's vector carries, the route, the floor, and what this lab calls its own plate of barcoded primer pairs. Replace all of those before you run it |
| [designs.tsv](../examples/ap1-readback/designs.tsv) | 72 designs, one a row: a name, and how many fragments it was built from |

The package ships no barcode sequences. The kit is 96 plasmids in four groups of 24, chaining
on `AGGA`, `GTTC`, `CCTT`, `TCAG` and `TTCC`, and its sequences are read from a copy your lab
holds, named by `LIULAB_SYNBIO_DMX_BARCODES`.
[The DMX barcode kit](../reference/dmx-barcode-kit.md) gives its structure and the columns a
copy must hold.

## Run it

```bash
pixi run synbio dmx plan docs/examples/ap1-readback/build.json --out readback/
```

```text
AP-1 cargo read-back: 32 of 72 designs read back on index PCR, 128 wells over 1 picked plate
readback/project.json
readback/index.html
readback/01-design-read-back-index-pcr.html
readback/reagents.html
readback/references.html
```

The floor is a fragment count, and this build set it at 2: the 32 designs stitched from two
pieces or more are read, and the 40 made in one are not. Set it to zero to read every design
on the sheet. Four colonies a design is 128 wells, which fits one 384-well plate and fills two
index plates.

Nothing here is designed, so the command writes a protocol folder, with no sequence file and
no primer sheet.

## What it wrote

| File | What it is |
| --- | --- |
| [index.html](../examples/ap1-readback/protocol/index.html) | the way in: what is read, what the bench is handed, and how long it is held |
| [01-design-read-back-index-pcr.html](../examples/ap1-readback/protocol/01-design-read-back-index-pcr.html) | the bench page: array and pick, mark every well, call every well |
| [reagents.html](../examples/ap1-readback/protocol/reagents.html) | what to order |
| [references.html](../examples/ap1-readback/protocol/references.html) | where each number was read from |
| [project.json](../examples/ap1-readback/protocol/project.json) | the run as data, which `mbio protocol render` turns back into these pages |

Each name links to what that run wrote, published here unedited; the same two inputs always
write the same bytes.

## What the bench works through

| Step | What it does |
| --- | --- |
| Array and grow | one spot a design on a BioAssay plate |
| Pick | four colonies a design into a 384-well plate, a quarter of the plate at a time |
| Sample | a quarter of the picked plate into each index plate |
| Mark | amplify every well with the pair its own address names |
| Pool and sequence | one pool a plate |
| Call | a consensus a well, read against that well's design |

A well's address is **worked out from where the well is, not looked up in a file**, so a
demultiplexer can check an address instead of trusting one. One axis of it carries the plate,
which is how two plates are told apart on one flow cell. The run above holds the bench
1 h 35 min across six steps.

## Choosing a route

| | Barcode ligation | Index PCR |
| --- | --- | --- |
| What marks a well | four kit barcodes, ligated into the construct in lysate | one barcoded primer pair |
| What addresses what | three groups the well, the fourth the plate | 96 forward marks the well, 96 reverse the plate |
| Where it runs | one 1536-well plate, compressed from four 384-well picked plates | one 96-well plate for each quarter of a picked plate |
| Wells it reaches | over 330,000 | 9,216, on 192 primers already held |
| Called above | 150 reads | 20 reads, and a well at 10 or more still carries a verdict, with a warning |
| What it asks of the bench | one incubator | one thermocycled reaction in every well |

**Do one, never both.** Barcode ligation pays $695 of sequencing and five days whatever the
library's size, so it costs less for each design the more designs share it: $19.52 a design
over 100, and $3.38 over 2,000. Nobody has published the same table for index PCR, so none is
stated here. How far each reaches does not decide it: both reach far more wells than most runs
need. A **reaction** is one tube of a mix, so one amplification across a 96-well array is one
reaction, not 96.

## What counts as a pass

Two questions, in order.

**Is the read deep enough to call?** If not, the well has **no verdict**. It is read again or
picked again, and that is not a failure. Reformatting compacts out the wells that failed and
leaves a well with no verdict in place. A build may raise the depth it wants; it cannot lower
the depth the route tolerates. The two routes' floors are not a strict and a lenient setting
of one scale, so neither carries over to the other.

**Does the call match the design exactly?** Both entry overhangs, the fragment, the stuffer
and the barcode. Anything less fails, a silent change included, because a barcode that no
longer names its member cannot be put right later. A well with more than one consensus is
mixed, and mixed fails.

## Before you order

- Read the way in and the bench page once, top to bottom.
- Name one route and hold to it.
- Check the prepared marking plate is on the shelf: a run calls for one and never builds it.
- Set the fragment-count floor deliberately, since a design in more pieces needs more picks.
- Spot from a copy of the archive plate, never the only one you hold.

## Tips and troubleshooting

**A well reads as the right protein and the wrong barcode.** It fails. A barcode-kit member
marks a well, where a barcode names a library part, and a barcode that no longer names its
member cannot be put right by a later read.

**A design comes back with no passing well.** Pick it again from the same archive spot before
re-synthesising it. Four colonies a design is the only measured anchor, and the chance of a
clean one falls as a design is built from more pieces.

**Nothing reads deeply enough.** Those wells carry no verdict rather than a failure. Read or
pick them again.

## Where the numbers come from

Each depth floor, reaction and cycling program is the published route's own, and every number
carries its **source** inside the protocol: the document, its edition, and the date it was
read.

Two departures from the paper change what the bench does. The vector was rebuilt to carry
kanamycin, because the last transfer moves cargo into a backbone carrying the other marker, so
every plate, well and broth here takes your own vector's drug rather than the published one.
And barcode ligation runs on a throwaway heat-lysed aliquot, so the culture in the well stays
the stock.

Two things nothing published fills. The package holds no index marks: the annealing regions
bind the vector as published, but which 192 sequences sit on their 5' ends is the plate your
lab holds. And the route's second criterion, a mean error below 10%, is not judged here,
because a well reaches this package as a read count and a call.

`synbio dmx plan` takes one option, `--out`. Everything else a run chooses is in the build
file.

## References

| Source | What it gives |
| --- | --- |
| Qian, Z. et al. (2026) Accelerating protein design by scaling experimental characterization. *Nat. Commun.* [doi:10.1038/s41467-026-76740-9](https://doi.org/10.1038/s41467-026-76740-9) | The barcode kit, barcoding 1536 wells in lysate, and a consensus called above 150 reads |
| Long, Y., Mora, A., Li, F.-Z., Gürsoy, E., Johnston, K. E. and Arnold, F. H. (2025) LevSeq: rapid generation of sequence-function data for directed evolution and machine learning. *ACS Synth. Biol.* [doi:10.1021/acssynbio.4c00625](https://doi.org/10.1021/acssynbio.4c00625) | Index PCR on one barcoded primer pair a well, its reaction and touchdown cycling, and twenty reads wanted with ten tolerated |
| Lund, S., Potapov, V., Johnson, S. R., Buss, J. and Tanner, N. A. (2024) Highly parallelized construction of DNA from low-cost oligonucleotide mixtures using Data-optimized Assembly Design and Golden Gate. *ACS Synth. Biol.* 13, 745-751 | Four colonies giving a clean copy of 343 of 458 genes, and how that rate falls with the number of fragments |
| New England Biolabs, Product Specification: Taq DNA Polymerase with ThermoPol Buffer, PS-M0267S/L/X/E v2.0, effective 12 February 2020 | The Taq stock concentration the index PCR is pipetted against, which the route's own paper never states |

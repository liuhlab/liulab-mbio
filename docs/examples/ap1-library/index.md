# Example: an AP-1 domain library

A combinatorial library over three positions, built from 72 protein sequences and nothing else.
Every output file below was written by the command, from the input files beside it.
Nothing here is edited by hand.

Read this page to see what a whole run writes. Its [project.json](project.json) holds the
choices this library made, bar the two the command below adds; copy it to start one of your
own. The method it is built by is the package's own and is not in the file.

```bash
pixi run liulab_synbio igga plan docs/examples/ap1-library/project.json \
  --out docs/examples/ap1-library \
  --working-site EGFP --prices docs/examples/ap1-library/prices.csv
```

`--working-site EGFP` says where the backbone's ccdB cassette goes, since this one carries
none. `--prices` points at the price list this lab holds.

To change what these files say, change the code and run that command again. The same inputs
write the same bytes, so a run that changes nothing leaves them alone.

## What goes in

| File | What it is |
| --- | --- |
| [parts.fasta](parts.fasta) | 72 proteins: 24 for each of the N, DBD and C positions |
| [project.json](project.json) | what this library chose: its three positions, its host, its completeness, its barcode rules, how its primers are plated, that every design is read back on both routes, and the numbers this lab set for itself |
| [primers.tsv](primers.tsv) | the orthogonal primer set the oligo pool is amplified by |
| [vector.gb](vector.gb) | the destination the first round opens |
| [working-vector.gb](working-vector.gb) | the backbone the finished library ends in: pLVX-TetOne-Puro-GFP, Addgene 171123, with its eight BsaI and BsmBI sites taken out. Swap in your own backbone and name a site in it |
| [prices.csv](prices.csv) | what this lab pays for each thing the run buys. Replace every row with your own quote |

`primer_plates` says how this lab lays the primers out: what the vendor delivers in a well,
what it is resuspended to, and what a working well holds. Leave it out and the primers are
ordered with the pool, and the run has one sitting fewer. The package states none of the three:
nothing publishes them, so they are the project author's to give.

Two more keys are this lab's own and not the method's. `linkage_fidelity` is the share of
reads whose barcode must still name its part, which this project puts at 0.9. `final_assembly`
sizes the last tube: 75 ng of working vector, with the cargo meeting it at 2:1. Leave either
out and the pages say the number has no source instead of printing one. **Replace both with
your own before you run this.** The pages label them as this run's own figures, so a reader who
copies the file and does not change them is carrying someone else's settings to the bench.

`routes` names the ways the picked wells are read back, and `validate_from` says which designs
are read at all: a fragment count, at or above which a design is read back. This project names
both routes and sets the floor to `0`, so all 72 designs are read — 288 wells and one 384-well
pick plate either way, and three 96-well index plates on index PCR. Leave the two keys out and
nothing is read, which is what a library headed for a pooled screen wants.

The two routes are two ways of doing one job. They are not two steps, and nobody does both. A
real run names the one it will do; this example names both so you can read each one, and every
page that offers them says to pick one. What to weigh is on
[the way in](protocol/index.html#topic-read-every-well-back). Barcode ligation pays a fixed
price once, so it costs less per design the more designs share it. Index PCR runs one
thermocycled reaction in every well, so a lot of wells need a lot of thermocyclers. How many
wells each route reaches does not decide it: both reach far more than this library needs.

Reading wells back is a job on its own too. [Reading the AP-1 cargo back](../ap1-readback/index.md)
takes the archive plate this run leaves, and a list of design names, and reads those wells back
with one command. It is a separate run, not a step of this chain.

The names say the position: `N_JUN` fills N, `DBD_JUN` fills DBD, `C_JUN` fills C. No overhang,
stuffer, barcode or codon is given. The planner chooses all four.

## What comes out

| File | What it is |
| --- | --- |
| [parts.tsv](parts.tsv) | the synthesis order sheet: 72 blocks, 133 to 1,149 bases |
| [pool.tsv](pool.tsv) | the oligo pool those blocks are synthesised as: 131 oligos, every one 350 nt |
| [pool-primers.tsv](pool-primers.tsv) | the 74 primers that amplify the pool, and how many oligos each one pulls out |
| [oligo.dna](oligo.dna) | one oligo of that pool, carrying the three primers that amplify it |
| [barcodes.tsv](barcodes.tsv) | which barcode names which part, and where it sits |
| [library-read-primers.tsv](library-read-primers.tsv) | the pairs that read linkage and representation back, with each amplicon |
| [changes.tsv](changes.tsv) | every amino acid the overhang standard moved, wild type beside synthesised |
| [round-1.dna](round-1.dna), [round-2.dna](round-2.dna) | one annotated record a round |
| [product.dna](product.dna) | one member of the finished library, 6,435 bases |
| [block-vector-1.dna](block-vector-1.dna), [block-vector-2.dna](block-vector-2.dna), [block-vector-3.dna](block-vector-3.dna) | the vector each position's blocks are built in, one a position |
| [working-vector-ccdb.dna](working-vector-ccdb.dna) | the backbone above with the ccdB cassette put in at `--working-site`, 10,284 bases: the vector the last sitting opens |
| [protocol/](protocol/index.html) | the bench protocols, as a folder of pages to work from |

## What the bench works through

The run is not one sitting. The `protocol/` folder holds one page for each, in the order
someone does them, and each page says what it is handed and what it leaves behind. Reading every
well back is one place in the run with a page for each route, and the bench works through one of
them:

| Page | What it is handed | What it leaves |
| --- | --- | --- |
| [Primer plates](protocol/01-primer-plates.html) | nothing yet | a primer stock plate, and a working copy |
| [Cargo ordering and pool preparation](protocol/02-cargo-ordering-and-pool-preparation.html) | nothing yet | the oligo pool |
| [Cargo creation](protocol/03-cargo-creation.html) | the pool, the working plate, and a block vector a position | one archived well a design |
| [Cargo validation: barcode ligation](protocol/04-cargo-validation-barcode-ligation.html) | the archive plate, and the DMX barcode kit | clonal wells, and a call for each |
| [Cargo validation: index PCR](protocol/04-cargo-validation-index-pcr.html) | the archive plate, and the index primer plate | clonal wells, and a call for each |
| [Library assembly in rounds](protocol/05-library-assembly-in-rounds.html) | the clonal wells and their calls | the library after round 3 |
| [Final cargo ligation](protocol/06-final-cargo-ligation.html) | the library after round 3 | the library in a working vector |

[The way in](protocol/index.html) explains the design, draws the chain, and lists the run's own
checks, how long it holds the bench and every number it has no source for.
[The reagents](protocol/reagents.html), which is what to order, and
[the references](protocol/references.html) are shared by every page. The whole chain is [protocol/project.json](protocol/project.json), which
`liulab_mbio protocol render` turns back into these pages after an edit.

A page names the one before it and the one after it, so each is complete on its own: mail one
page to whoever runs that sitting, or zip the folder and send the run.

The three part lists make 24 x 24 x 24 = 13,824 members. The product file holds one of them,
with the rest differing only in which protein and which barcode sits at each position.

The blocks of `parts.tsv` are not bought. The protocol orders the pool and its primers, pulls
each batch out of the pool, pulls each block out of its batch, and clones that block's cargo
into the vector of the position it fills. Only the cargo is synthesised: the stuffers either
side of it are the vector's own bases. The bill buys the oligos and the primers; it never buys
the blocks as well.

A part enters on the overhang the round before it leaves behind, so each position needs a vector
offering that overhang. The three `block-vector` files are one vector written three times, four
bases apart. The first is `vector.gb` itself, which is also what round one opens.

The protocol ends by moving the finished library into a working vector, which is where an
application gets it. The `working_vector` key names the backbone this library ends in, and the
steps then carry its enzyme and that enzyme's cycling. The plan reads the cargo enzyme off that
backbone first, then keeps every block clear of it: an enzyme a block spells cannot be the one
that admits the library. Leave the key out and those steps say what a vector would have fixed
instead. A backbone carrying no ccdB cassette needs `--working-site` as well, which says where
one goes; here that is EGFP, the feature the cassette goes in at.

## Where the working vector comes from

[working-vector.gb](working-vector.gb) is not the file Addgene sends. Addgene 171123 reads six
BsaI and two BsmBI sites, and a round of this method would cut every one of them. Taking them
out is a job of its own, with [its own bench page](working-vector-domestication.html) and
[the protocol behind it](working-vector-domestication.json), which lists every base that
moved. The library run opens with that record already in hand and never touches the backbone
again.

**Two of the eight are the hard ones, and the demo says so rather than hiding it.** They sit in
TAR, the stem-loop at each end of the viral genome. Five of the eight sit outside a gene too,
but TAR has to fold to work, so the change alters the shape the virus needs rather than a
spelling nothing reads. Whether the virus tolerates that is still open: the one published
backbone carrying the change reports no titre for it. The bench page carries that as a hole
rather than printing a number nobody measured. This demo keeps a hard backbone on purpose.
Domesticating one is worth seeing.

## What a price costs

[prices.csv](prices.csv) is what this lab pays. It is a CSV of `key`, `item`, `bands`, `charge`,
`basis` and `currency`, found the way the ligase matrix is, and a row nobody priced leaves a
hole rather than a guess. The nine keys are the bill's own: `oligo-pool`, `pool-primers`, each
enzyme's catalogue number, `60242-2`, `cuvettes` and `plasmid prep`.

Every row says what kind of price it is and the day it was read. Seven are the maker's own list
price. One, the oligo pool, is this lab's account pricing, because Twist publishes no pool
prices in public. One, the cells, is a distributor price that may carry a surcharge. Each
enzyme row divides a pack price by what the pack holds, so the money beside 7.5 µL is what
those 7.5 µL cost rather than what a vial costs. **Prices move.** Check every row against the
seller before you order, and swap in your own account's numbers.

## The figures

Each record is drawn by the same command, from the file beside it. A map labels only the
enzymes that cut a record once, so the four this method uses are named: all but SrfI cut twice,
and the maps would otherwise show none of them.

```bash
D=docs/examples/ap1-library
pixi run liulab_mbio plot map $D/product.dna \
  --enzyme BbsI --enzyme BsaI --enzyme PmeI --enzyme SrfI -o $D/product-map.pdf
pixi run liulab_mbio plot map $D/product.dna \
  --region 1368..1442 --sequence-view -o $D/barcode-block.pdf
```

| File | What it is |
| --- | --- |
| [vector-map.pdf](vector-map.pdf) | the destination the first round opens |
| [round-1-map.pdf](round-1-map.pdf), [round-2-map.pdf](round-2-map.pdf) | what each round leaves |
| [product-map.pdf](product-map.pdf) | one finished member |
| [barcode-block.pdf](barcode-block.pdf) | the 41 bp barcode block and the stuffer ahead of it: 75 bases, every one shown on the second page |

The gel and the plate layouts are in the protocol pages, drawn where the step that uses them is.

## Where the proteins come from

The 72 sequences are the AP-1 domains of Takacsi-Nagy, O. et al. (2026) Synthetic transcription
factors designed by domain recombination enhance CAR T cell antitumor function. *Cell* 189,
1-20, [doi:10.1016/j.cell.2026.07.054](https://doi.org/10.1016/j.cell.2026.07.054), used under
[CC BY 4.0](http://creativecommons.org/licenses/by/4.0/). They were read from that paper's
Table S1 and translated back to protein; the DNA here is this package's own.

[Highly parallel DNA synthesis and assembly](../../synthesis-and-assembly.md) explains the
method these files follow.

# Example: an AP-1 domain library

A combinatorial library over three positions, built from 72 protein sequences and nothing else.
Every output file below was written by the command, from the input files beside it.
Nothing here is edited by hand.

Read this page to see what a whole run writes. Its [project.json](project.json) holds every
choice this library made; copy it to start one of your own. The method it is built by is the
package's own and is not in the file.

```bash
pixi run liulab_synbio igga plan docs/examples/ap1-library/project.json \
  --out docs/examples/ap1-library
```

To change what these files say, change the code and run that command again. The same inputs
write the same bytes, so a run that changes nothing leaves them alone. Point the command at
`project-barcode-ligation.json` and an `--out` of your own to see the other read-back route;
the design is the same and only the protocol differs.

## What goes in

| File | What it is |
| --- | --- |
| [parts.fasta](parts.fasta) | 72 proteins: 24 for each of the N, DBD and C positions |
| [project.json](project.json) | what this library chose: its three positions, its host, its completeness, its barcode rules, and that every design is read back by index PCR |
| [project-barcode-ligation.json](project-barcode-ligation.json) | the same library read back by barcode ligation instead |
| [primers.tsv](primers.tsv) | the orthogonal primer set the oligo pool is amplified by |
| [vector.gb](vector.gb) | the destination the first round opens |

The two project files differ in one key. `route` says which of the two read-back routes reads
the picked wells, and `validate_from` says which designs are read at all: a fragment count, at
or above which a design is read back. Both files set it to `0`, so all 72 designs are read — 288
wells, one 384-well pick plate and three 96-well index plates. Leave the two keys out and
nothing is read, which is what a library headed for a pooled screen wants.

The names say the position: `N_JUN` fills N, `DBD_JUN` fills DBD, `C_JUN` fills C. No overhang,
stuffer, barcode or codon is given. The planner chooses all four.

## What comes out

| File | What it is |
| --- | --- |
| [parts.tsv](parts.tsv) | the synthesis order sheet: 72 blocks, 133 to 1,149 bases |
| [pool.tsv](pool.tsv) | the oligo pool those blocks are synthesised as: 131 oligos, every one 350 nt |
| [pool-primers.tsv](pool-primers.tsv) | the 74 primers that amplify the pool, and how many oligos each one pulls out |
| [barcodes.tsv](barcodes.tsv) | which barcode names which part, and where it sits |
| [library-read-primers.tsv](library-read-primers.tsv) | the pairs that read linkage and representation back, with each amplicon |
| [changes.tsv](changes.tsv) | every amino acid the overhang standard moved, wild type beside synthesised |
| [round-1.dna](round-1.dna), [round-2.dna](round-2.dna) | one annotated record a round |
| [product.dna](product.dna) | one member of the finished library, 6,435 bases |
| [block-vector-1.dna](block-vector-1.dna), [block-vector-2.dna](block-vector-2.dna), [block-vector-3.dna](block-vector-3.dna) | the vector each position's blocks are built in, one a position |
| [protocol/](protocol/index.html) | the bench protocols, as a folder of pages to work from |

## What the bench works through

The run is not one sitting. The `protocol/` folder holds one page for each, in the order
someone does them, and each page says what it is handed and what it leaves behind:

| Page | What it is handed | What it leaves |
| --- | --- | --- |
| [Cargo ordering and pool preparation](protocol/01-cargo-ordering-and-pool-preparation.html) | nothing yet | the oligo pool |
| [Cargo creation](protocol/02-cargo-creation.html) | the pool, and a block vector a position | one archived well a design |
| [Cargo validation](protocol/03-cargo-validation-index-pcr.html) | the archive plate, and the index primer plate | clonal wells, and a call for each |
| [Library assembly in rounds](protocol/04-library-assembly-in-rounds.html) | the clonal wells and their calls | the library after round 3 |
| [Final cargo ligation](protocol/05-final-cargo-ligation.html) | the library after round 3 | the library in a working vector |

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
application gets it. Add a `working_vector` key to name the backbone yours ends in and the
steps carry its enzyme and that enzyme's cycling; this project names none, so those steps say
what a vector would have fixed instead. Name one and the plan reads its cargo enzyme off that
backbone first, then keeps every block clear of it: an enzyme a block spells cannot be the one
that admits the library. A backbone carrying no ccdB cassette yet needs `--working-site` as
well, which says where one goes.

## What a price costs

Add a price record of your own and the bill carries money beside every quantity:

```bash
pixi run liulab_synbio igga plan docs/examples/ap1-library/project.json \
  --out library/ --prices prices.csv
```

The package ships no prices, so none are here. A record is a CSV of `key`, `item`, `bands`,
`charge`, `basis` and `currency`, found the way the ligase matrix is, and a row nobody priced
leaves a hole rather than a guess. The keys this plan asks for are the bill's own: `oligo-pool`,
`pool-primers`, each enzyme's catalogue number, `60242-2`, `cuvettes` and `plasmid prep`.

## The figures

Each record is drawn by the same command, from the file beside it:

```bash
D=docs/examples/ap1-library
pixi run liulab_mbio plot map $D/product.dna -o $D/product-map.pdf
pixi run liulab_mbio plot map $D/product.dna \
  --region 1368..1442 --sequence-view -o $D/barcode-block.pdf
```

| File | What it is |
| --- | --- |
| [vector-map.pdf](vector-map.pdf) | the destination the first round opens |
| [round-1-map.pdf](round-1-map.pdf), [round-2-map.pdf](round-2-map.pdf) | what each round leaves |
| [product-map.pdf](product-map.pdf) | one finished member |
| [barcode-block.pdf](barcode-block.pdf) | the 75 bases the barcode block reads, with the bases shown |

The gel and the plate layouts are in the protocol pages, drawn where the step that uses them is.

## Where the proteins come from

The 72 sequences are the AP-1 domains of Takacsi-Nagy, O. et al. (2026) Synthetic transcription
factors designed by domain recombination enhance CAR T cell antitumor function. *Cell* 189,
1-20, [doi:10.1016/j.cell.2026.07.054](https://doi.org/10.1016/j.cell.2026.07.054), used under
[CC BY 4.0](http://creativecommons.org/licenses/by/4.0/). They were read from that paper's
Table S1 and translated back to protein; the DNA here is this package's own.

[Highly parallel DNA synthesis and assembly](../../synthesis-and-assembly.md) explains the
method these files follow.

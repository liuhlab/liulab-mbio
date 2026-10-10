# iGGA: a barcoded combinatorial library

Iterative Golden Gate appends one list of parts a round, so a few dozen proteins become
thousands of plasmids. Each part carries an 11 bp barcode at its junction, so sequencing a
member says which protein sits at every position. One command designs the synthesis order,
draws the barcodes, simulates the rounds, sizes the colonies each round needs and writes the
bench protocol.

## When to use it

- Two or more lists of protein or DNA parts, and you want every combination as a plasmid.
- You want a member read back from sequencing alone, with nothing to look up.
- The library is too large to build member by member.

Use [Golden Gate assembly](../methods/golden-gate.md) for one defined construct instead, and
[DMX](dmx.md) to read picked wells back.

## What you need

One file. It names a FASTA of the parts, the vector the first round opens, the backbone the
library ends in, and this run's own numbers. A **project** is the chain of protocols; a
**build** is one run's choices, and the file holding a build is called `project.json`. The run
below uses files that ship with the repo. Copy them to start one of your own.

| Input | What it is |
| --- | --- |
| [project.json](../examples/ap1-library/project.json) | the build: three positions, the host, the completeness, the barcode rules, and the numbers this lab set for itself |
| [parts.fasta](../examples/ap1-library/parts.fasta) | 72 proteins, 24 for each of the N, DBD and C positions. The name says the position; no overhang, stuffer, barcode or codon is given |
| [vector.gb](../examples/ap1-library/vector.gb), [carrier.gb](../examples/ap1-library/carrier.gb) | the destination the first round opens, and the plasmid each part is seated in. Leave `carrier` out and nothing is seated, which is what a lab holding its parts wants |
| [working-vector.gb](../examples/ap1-library/working-vector.gb) | the backbone the finished library ends in. Swap in your own and name a site in it |
| [primers.tsv](../examples/ap1-library/primers.tsv), [prices.csv](../examples/ap1-library/prices.csv) | the orthogonal primer set the pool is amplified by, and what this lab pays for each thing the run buys. Replace every price row with your own quote |

## Run it

```bash
pixi run synbio igga plan docs/examples/ap1-library/project.json --out library/ \
  --working-site EGFP --prices docs/examples/ap1-library/prices.csv
```

`--working-site` says where the ccdB cassette goes in a backbone carrying none; `--prices`
points at the price list this lab holds. The run states its design in one line, then every
file:

```text
AP-1 DESynR round 3: 6435 bp, 72 parts in 3 lists, 13,824 constructs, 3 rounds, entry overhangs AGGA, AGAT, GCAT, scar TTCC, 27 amino acid changes, checks pass
library/parts.tsv
library/barcodes.tsv
library/changes.tsv
library/round-1.dna
library/round-2.dna
library/product.dna
library/block-vector-1.dna
library/block-vector-2.dna
library/block-vector-3.dna
library/working-vector-ccdb.dna
library/protocol/project.json
library/protocol/index.html
library/protocol/01-part-carrier.html
library/protocol/02-primer-plates.html
library/protocol/03-cargo-ordering-and-pool-preparation.html
library/protocol/04-cargo-creation.html
library/protocol/05-cargo-validation-barcode-ligation.html
library/protocol/05-cargo-validation-index-pcr.html
library/protocol/06-library-assembly-in-rounds.html
library/protocol/07-final-cargo-ligation.html
library/protocol/reagents.html
library/protocol/references.html
library/pool.tsv
library/pool-primers.tsv
library/oligo.dna
library/library-read-primers.tsv
```

Read that line first. 72 parts in three lists make 13,824 constructs. Each list enters on its
own overhang and every part leaves the same `TTCC` scar. The planner picked that set because it
costs the proteins fewest residues: the 27 changes in `changes.tsv` are what you pay for it.

A round opens the library with BbsI and releases one part list with BsaI, in separate tubes,
then ligates, electroporates and preps the result as the next round's destination. A 34 bp
internal stuffer is kept each round, which is what lets the next round open the product; a
capping block ends the chain instead. The 41 bp barcode block reads in reverse of the round
order.

The last round needs 195,386 colonies for a 0.99 chance that none of its 13,824 products is
missing, 14 times its products. That is **library coverage**: colonies against distinct
products, not sequencing depth. A round below its floor loses members no later round puts back.

The run ends by moving the library into a working vector, where an application gets it. The
plan reads that backbone's cargo enzyme first, then keeps every block clear of it: an enzyme a
block spells cannot be the one that admits the library.

## What it wrote

| File | What it is |
| --- | --- |
| [parts.tsv](../examples/ap1-library/parts.tsv) | the synthesis order sheet: 72 blocks, 133 to 1,149 bases, each with its barcode on the same row |
| [pool.tsv](../examples/ap1-library/pool.tsv), [pool-primers.tsv](../examples/ap1-library/pool-primers.tsv), [oligo.dna](../examples/ap1-library/oligo.dna) | the oligo pool those blocks are synthesised as, 131 oligos every one 350 nt; the 74 primers that pull batches out of it; and one oligo drawn with its three primers |
| [barcodes.tsv](../examples/ap1-library/barcodes.tsv) | which barcode names which part, and where it sits. This is what turns a read back into a list of parts |
| [changes.tsv](../examples/ap1-library/changes.tsv) | every amino acid the overhang standard moved, wild type beside synthesised |
| [round-1.dna](../examples/ap1-library/round-1.dna), [round-2.dna](../examples/ap1-library/round-2.dna) | one annotated record a round. Opens in SnapGene |
| [product.dna](../examples/ap1-library/product.dna) | one member of the library, 6,435 bases, before the move into the working vector. The product is the plasmid; an assembly product is a kit |
| [block-vector-1.dna](../examples/ap1-library/block-vector-1.dna), [-2](../examples/ap1-library/block-vector-2.dna), [-3](../examples/ap1-library/block-vector-3.dna) | the vector each position's blocks are built in, one a position |
| [working-vector-ccdb.dna](../examples/ap1-library/working-vector-ccdb.dna) | the backbone with a ccdB cassette put in at `--working-site`, 10,284 bases: what the last sitting opens |
| [library-read-primers.tsv](../examples/ap1-library/library-read-primers.tsv) | the pairs that read linkage and representation back, with each amplicon |
| [protocol/round-1-map.html](../examples/ap1-library/protocol/round-1-map.html) and one like it for each record a figure draws | that record as a map you can explore, which the figure on the protocol page opens |
| [protocol/](../examples/ap1-library/protocol/index.html) | the bench pages, and `project.json`, the whole chain as data |

Each name links to what that run wrote, published here unedited. Nothing is typed by hand, so
the same inputs always write the same bytes. Change the build and run the command again to
change the design. To change what a page says, edit the chain's own data file and render it
again:

```bash
pixi run mbio protocol render library/protocol/project.json
```

Figures, drawn by `mbio plot map`:

| Figure | What it draws |
| --- | --- |
| [vector-map.pdf](../examples/ap1-library/vector-map.pdf) | the destination the first round opens |
| [round-1-map.pdf](../examples/ap1-library/round-1-map.pdf), [round-2-map.pdf](../examples/ap1-library/round-2-map.pdf) | what each round leaves |
| [product-map.pdf](../examples/ap1-library/product-map.pdf) | one finished member, with the four enzymes this method uses labelled |
| [barcode-block.pdf](../examples/ap1-library/barcode-block.pdf) | the 41 bp barcode block and the stuffer ahead of it, every base shown |

## What the bench works through

The run is not one sitting. The `protocol/` folder holds a page for each, in order.

| Page | Handed | Leaves |
| --- | --- | --- |
| [Part carrier](../examples/ap1-library/protocol/01-part-carrier.html) | the 72 parts | one carrier plasmid a part, one part a well |
| [Primer plates](../examples/ap1-library/protocol/02-primer-plates.html) | nothing yet | a primer stock plate, and a working copy |
| [Cargo ordering and pool preparation](../examples/ap1-library/protocol/03-cargo-ordering-and-pool-preparation.html) | nothing yet | the oligo pool |
| [Cargo creation](../examples/ap1-library/protocol/04-cargo-creation.html) | the pool, the working plate, and a block vector a position | one archived well a design |
| [Cargo validation: barcode ligation](../examples/ap1-library/protocol/05-cargo-validation-barcode-ligation.html) | the archive plate, and the DMX barcode kit | clonal wells, and a call for each |
| [Cargo validation: index PCR](../examples/ap1-library/protocol/05-cargo-validation-index-pcr.html) | the archive plate, and the index primer plate | clonal wells, and a call for each |
| [Library assembly in rounds](../examples/ap1-library/protocol/06-library-assembly-in-rounds.html) | the carrier plate, and the clonal wells and their calls | the library after round 3 |
| [Final cargo ligation](../examples/ap1-library/protocol/07-final-cargo-ligation.html) | the library after round 3 | the library in a working vector |

Reading every well back has a page for each route, and the bench does one, never both;
[DMX](dmx.md) covers that job.
[The way in](../examples/ap1-library/protocol/index.html) draws the chain, lists the checks and
names every number with no source: 49 steps, 31 of which hold the bench. Each page names the
one before and after it, so you can mail one to whoever runs that sitting.

## Before you order

- Read the way in, then every page of the chain.
- Replace this build's numbers with yours: completeness, barcode rules, linkage fidelity and
  final assembly.
- Decide that each residue change in `changes.tsv` is harmless in your proteins.
- Check each round's colony floor against the cells and cuvettes you hold.
- Re-price every row against your own vendors; a **band** here is a price band, not a gel
  band.

## Tips and troubleshooting

**Your backbone carries Type IIS sites.** The
[working vector](../examples/ap1-library/working-vector.gb) here is not the file Addgene sends.
Addgene 171123 reads six BsaI and two BsmBI sites and a round would cut every one, so all eight
came out first. Two sit in TAR, the stem-loop at each end of the viral genome, which has to
fold to work: removing a site there changes a shape the virus needs. **Whether it still packages
is open** — the one published backbone carrying that change reports no titre. Measure it.

**A round comes up short of colonies.** Plate a dilution of every round beside a no-donor
control, and subtract the control. Run the round again rather than carrying on: members lost
in one round are missing from all of them after it.

**A part's protein comes back changed.** A junction spells whole codons, so the overhang
decides the residues either side of it. The planner scores every candidate set across all the
parts and keeps the cheapest, so a change in `changes.tsv` is one no other set avoided.

## Where the numbers come from

This run's settings are **one build's answers, not the method's defaults**: its completeness,
its barcode rules, its linkage fidelity of 0.9, and the 75 ng of working vector its last tube
takes. The build states no cargo mass, because cargo is never pipetted — the whole release goes
in, so the digest fixes how much meets the vector. The pages label each of those as this run's
own figure. Copy the file without replacing them and you carry someone else's settings.

Volumes, units, times and temperatures are cited step by step inside each page, and
`references.html` lists every document. Where nothing publishes a number the page prints a
hole, not a guess: 35 here, four at the bench and 31 prices, since this lab priced nine of the
bill's 40 rows. Round fidelity is measured ligation data — see
[how fidelity is predicted](../reference/ligation-fidelity.md). A **barcode** names a library
part here, where a barcode-kit member marks a well;
[barcode sets](../methods/barcodes.md) says how a set is drawn.

`pixi run synbio igga plan --help` lists the rest: what the parts FASTA holds, where a stuffer
or cassette goes, how a record name is matched to a position, and reading a ligase matrix you
hold.

## References

| Source | What it gives |
| --- | --- |
| Takacsi-Nagy, O. et al. (2026) Synthetic transcription factors designed by domain recombination enhance CAR T cell antitumor function. *Cell* 189, 1-20. [doi:10.1016/j.cell.2026.07.054](https://doi.org/10.1016/j.cell.2026.07.054) | Assembly in rounds, the internal stuffer and the 11 bp barcode length. The 72 AP-1 domains above were read from its Table S1 under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) and translated back to protein; the DNA here is this package's own |
| Lund, S., Potapov, V., Johnson, S. R., Buss, J. and Tanner, N. A. (2024) Highly parallelized construction of DNA from low-cost oligonucleotide mixtures using Data-optimized Assembly Design and Golden Gate. *ACS Synth. Biol.* 13, 745-751 | Building each part from a cheap oligo pool: the fragment split and its overhangs |
| Subramanian, S. K., Russ, W. P. and Ranganathan, R. (2018) A set of experimentally validated, mutually orthogonal primers for combinatorially specifying genetic components. *Synth. Biol.* 3, ysx008 | The orthogonal primers the pool is amplified by |
| Clarke, L. and Carbon, J. (1976) A colony bank containing synthetic ColE1 hybrid plasmids representative of the entire E. coli genome. *Cell* 9, 91-99 | The colonies a round needs to hold every member |
| Qian, Z. et al. (2026) Accelerating protein design by scaling experimental characterization. *Nat. Commun.* [doi:10.1038/s41467-026-76740-9](https://doi.org/10.1038/s41467-026-76740-9) | The vector each part is built in, and the barcode kit that marks its wells |

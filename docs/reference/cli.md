# Command line

One install gives two commands. `mbio` is the toolkit: sequence files, enzymes, primers and
cloning. `synbio` holds the methods this lab named. Every command on this site starts with one
of those two words.

Run them through pixi, from a checkout of the repo:

```bash
pixi run mbio version
pixi run synbio version
```

Both print the same number, because one install ships both.

Words in capitals are yours to fill in. `DIR` is a folder, `VECTOR` is a sequence file, and so
on. Every command below appears once as a line in a `bash` block, and a test keeps that list
and the commands in step. Each command also wraps a function you can call yourself, which
[the Python reference](python.md) covers.

## What each command does

| Command | What it does |
| --- | --- |
| `mbio version` | Print the installed version |
| `mbio codon-optimize` | Write a protein as DNA for a host, free of named sites |
| `mbio barcode-design` | Draw a set of barcodes, each far enough from the rest |
| `mbio cloning goldengate plan` | Plan a Golden Gate assembly |
| `mbio cloning gibson plan` | Plan a Gibson assembly |
| `mbio cloning restriction plan` | Plan restriction and ligation cloning |
| `mbio cloning gateway plan` | Plan the BP and LR reactions of a Gateway cloning |
| `mbio protocol render` | Turn protocol data back into a page |
| `mbio plot map` | Draw a record as a map |
| `mbio primers design` | Design a PCR primer pair |
| `mbio sequence-verify` | Check one clone's sequencing results against the plasmid it should be |
| `synbio version` | Print the installed version |
| `synbio dmx plan` | Plan a read-back of designs the lab already holds |
| `synbio igga plan` | Plan a barcoded combinatorial library |

## How positions are written

Two spellings, and they count differently.

`--site` and `--working-site` take `START-END` counted from 0, with the end one base past the
last base replaced. `--region`, on `plot map` and on `primers design`, takes `START..END`
counted from 1 with both ends included. That is how a map prints a position, so you can read a
span off a map and type it straight back in.

## Write a protein as DNA

`mbio codon-optimize` writes a protein as DNA in the codons a host uses most often. Hand it DNA
instead and it checks the coding sequence you already have. Either way it clears out any site
you forbid, by swapping a codon for one that reads as the same amino acid.

```bash
pixi run mbio codon-optimize SEQUENCE --kind protein --host e-coli-k12 --forbid BsaI --out DIR
```

| Argument or option | Default | What it is |
| --- | --- | --- |
| `SEQUENCE` | required | The protein or the coding sequence, as letters |
| `--kind` | required | What `SEQUENCE` is: `protein` or `dna` |
| `--host` | required | Codon usage table to write for: `e-coli-k12`, `human` or `mouse` |
| `--forbid` | none | An enzyme whose site must not appear. Repeat it once per enzyme |
| `--name` | `coding sequence` | What to call the sequence |
| `--out`, `-o` | none | Folder to write into |

With `--out` it writes two files: `coding-sequence.dna`, and `codon-changes.tsv` listing every
codon it moved. Without `--out` it prints the changes and the bases instead.

The enzymes and the codon usage tables this package ships are listed in
[the shipped data](data.md).

## Design a set of barcodes

`mbio barcode-design` draws barcodes that stand far enough apart that no two read as each other.
It keeps them free of the sites you name, scar included, and off a stop codon wherever the
construct reads through them.

```bash
pixi run mbio barcode-design COUNT --length 10 --out DIR
```

| Argument or option | Default | What it is |
| --- | --- | --- |
| `COUNT` | required | How many barcodes the part list needs |
| `--length` | required | How many bases name one part |
| `--out`, `-o` | required | Folder to write into |
| `--scar` | none | The cloning scar joining one barcode to the next |
| `--phase` | `0` | Bases of the barcode's first codon read before it: 0, 1 or 2 |
| `--untranslated` | off | The construct never reads the barcode, so drop the frame rules |
| `--forbid` | none | An enzyme whose site must not appear. Repeat it once per enzyme |
| `--distance` | `3` | The fewest units two barcodes stand apart |
| `--metric` | `sequence-levenshtein` | How a distance is counted. Also `hamming` |
| `--max-homopolymer` | `5` | The longest run of one base allowed. `0` lifts the cap |
| `--gc-band` | none | The share of G and C allowed, as `LOW-HIGH`, each a share of 1 |
| `--seed` | `0` | The seed the draw is made with |

It writes the set as `barcodes.tsv` and its verdicts as `checks.tsv`. The same seed and the same
rules draw the same set again.

## Plan a cloning experiment

Four methods, one command each, all under `mbio cloning`. Each reads SnapGene `.dna`, GenBank or
FASTA files, and writes five files into `--out`:

| File | What it holds |
| --- | --- |
| `product.dna` | the plasmid the design makes, annotated |
| `primers.tsv` | the oligos to order |
| `protocol.json` | the bench protocol as data you can edit |
| `product-map.html` | the product as a map you can explore, which the page's figure opens |
| `protocol.html` | that protocol as one page, which opens offline |

Gateway also writes `entry-clone.dna` and its map, `entry-clone-map.html`, where it plans a BP
reaction. Golden Gate also writes each amplicon, `NAME-amplicon.dna`, and its map, which the
PCR step's figure opens.

Every one of the four prints a summary line, then each file it wrote. Five options are shared:

| Option | Default | What it is |
| --- | --- | --- |
| `--out`, `-o` | required | Folder to write into |
| `--polymerase` | `Q5` | Polymerase for the PCRs. Also `Phusion`, `Taq` and `OneTaq` |
| `--host` | `NEB 5-alpha Competent E. coli (C2987)` | Strain the protocol names |
| `--cleanup-kit` | `T1130` (Monarch Spin PCR & DNA Cleanup Kit) | Spin-column kit the protocol names, by catalogue number or by name |
| `--name` | the parts' names, joined | What to call the product |

Gateway's `--host` starts at `TOP10 chemically competent E. coli` instead. A strain it
recognises also gets its ccdB selection checked.

`--cleanup-kit` takes a catalogue number the package knows — `T1130`, the Monarch Spin PCR &
DNA Cleanup Kit, or `D4003`, Zymo Research's DNA Clean & Concentrator-5 — or any other kit's own
name, which the page then names with no catalogue number beside it.

### Golden Gate

```bash
pixi run mbio cloning goldengate plan VECTOR INSERT --out DIR
```

Name as many insert files as you like. They go round the product in the order you give them.

| Option | Default | What it is |
| --- | --- | --- |
| `--site` | the vector's `MCS` feature | Feature name, or `START-END`, that the inserts replace |
| `--orientation` | `forward` | Which way round an insert goes: `forward` or `reverse`. Once per insert, or once for all of them |
| `--in-frame` | off | Hold every junction on a codon boundary |
| `--enzyme` | the best free one | Type IIS enzyme to use |
| `--codon-table` | `e-coli-k12` | Whose codon usage a proposed domestication picks from |
| `--ligase-matrix` | none | A ligase fidelity file you hold, to score the overhangs of an enzyme nobody has measured |
| `--prefer-ligase-matrix` | off | Score with that file even where the enzyme has a measured one |

### Gibson

```bash
pixi run mbio cloning gibson plan VECTOR INSERT --out DIR
```

| Option | Default | What it is |
| --- | --- | --- |
| `--site` | the vector's `MCS` feature | Feature name, or `START-END`, that the inserts replace |
| `--orientation` | `forward` | `forward` or `reverse`. Once per insert, or once for all of them |
| `--route` | `amplify` | How an insert is made: `amplify` or `stitch`. Once per insert |
| `--bridge` | none | A junction joined by one oligo, written `BEFORE:AFTER`. Once per junction |
| `--product` | `NEBuilder HiFi DNA Assembly Master Mix` | The assembly product on the bench |

The other two products are `Gibson Assembly Master Mix` and
`In-Fusion Snap Assembly Master Mix`.

### Restriction and ligation

```bash
pixi run mbio cloning restriction plan VECTOR INSERT --out DIR
```

The second file is the insert, or the plasmid it is cut out of.

| Option | Default | What it is |
| --- | --- | --- |
| `--enzyme` | chosen for you | An enzyme both digests use. Name it at most twice |

Name no enzyme and the plan weighs the pairs itself, then says in its summary how many it turned
down.

### Gateway

```bash
pixi run mbio cloning gateway plan CARRIER DESTINATION --out DIR
```

`CARRIER` is the entry clone. Add `--donor` and it is the insert instead, and the plan then
covers the BP reaction that makes the entry clone.

| Option | Default | What it is |
| --- | --- | --- |
| `--donor` | none | Donor vector, to plan the BP reaction |
| `--amplify` / `--no-amplify` | `--no-amplify` | Design attB primers for an insert carrying no att site, and amplify it |
| `--fusion` | `none` | Tag the insert is read into: `none`, `N-terminal`, `C-terminal` or `both` |

## Turn protocol data back into a page

Every plan writes its protocol twice: as data, and as the page rendered from it. Edit the data,
then render it again with `mbio protocol render`. The page hands the bench back the check marks
it had already ticked.

```bash
pixi run mbio protocol render SOURCE
```

`SOURCE` is one `protocol.json` file, or a folder holding a `project.json`. A single file is
written beside itself with an `.html` suffix, or wherever `--output` says. A folder is rendered
in place: the index, one page per protocol, and the reagents and references pages. Those pages
link each other, so they have to stay together, and `--output` is refused for a folder.

| Option | Default | What it is |
| --- | --- | --- |
| `--output`, `-o` | `SOURCE` with `.html` | HTML file to write. Not allowed for a folder |

## Draw a record as a map

`mbio plot map` draws a record the way SnapGene shows one: a circle or a line, with features,
primers and cut sites labelled so no label sits on another.

```bash
pixi run mbio plot map RECORD --output map.html
```

The suffix picks the format. Use `.html` for a page to look through, `.png` or `.pdf` for a
figure. Repeat `--output` to write several at once.

| Argument or option | Default | What it is |
| --- | --- | --- |
| `RECORD` | required | Sequence file: `.dna`, GenBank or FASTA |
| `--output`, `-o` | required | File to write. Repeat it to write several |
| `--region` | the whole record | A feature's name, or `START..END`. An end before the start runs across the origin |
| `--linear` | off | Open a circular record out as a line |
| `--sequence-view` | off | Show the bases with the map, at most 100,000 of them |
| `--no-features` | off | Leave the features out |
| `--no-primers` | off | Leave the primers out |
| `--no-cut-sites` | off | Leave the cut sites out |
| `--enzyme` | the shipped enzymes that cut once | An enzyme whose every cut site to draw. Once per enzyme |
| `--hide-type` | none | A feature type to leave out. Once per type |
| `--source` | off | Show the source feature |
| `--bases-per-row` | `60` | How many bases a row of the sequence view holds |
| `--one-strand` | off | Leave the bottom strand out of the sequence view |
| `--dpi` | `300` | Resolution of a PNG, in dots per inch |

It prints each file it wrote. If any label had to be hidden for want of room, it says how many.
A page keeps whatever you switched off behind its own switches; a PNG or a PDF leaves it out.

## Design a primer pair

`mbio primers design` picks a pair over a region, judges it, and writes the sheet you order
from.

```bash
pixi run mbio primers design TEMPLATE --out DIR
```

| Argument or option | Default | What it is |
| --- | --- | --- |
| `TEMPLATE` | required | A `.dna`, GenBank or FASTA record, or an indexed genome FASTA |
| `--out`, `-o` | required | Folder to write into |
| `--region` | the whole record | The amplicon, as `START..END`, or `NAME:START..END` on a genome, where it is required |
| `--assembly` | none | The genome's name, such as `hg38` |
| `--flank` | `200` | How far either side of a genome region a primer may lie |
| `--forward-tail` | none | Bases joined to the forward primer's 5' end |
| `--reverse-tail` | none | Bases joined to the reverse primer's 5' end |
| `--forward-name` | `forward` | What to order the forward primer under |
| `--reverse-name` | `reverse` | What to order the reverse primer under |
| `--target-tm` | `62` | The melting temperature the design aims for, in °C |
| `--polymerase` | `Q5` | Polymerase whose rules judge the pair |

Give `--assembly` and `TEMPLATE` is read as that genome's FASTA. Name `--region` too, as
`NAME:START..END`, because a whole genome is no amplicon. The pair is then searched against the
genome itself and moved off any off-target amplicon.

It writes `primers.tsv`. On a genome it writes `genome.tsv` too, holding every amplicon the pair
makes there, the intended one first.

## Check a clone's sequencing

`mbio sequence-verify` holds one clone's sequencing results against the plasmid it should be,
and gives each junction and each insert a verdict.

```bash
pixi run mbio sequence-verify PRODUCT RESULT --out DIR
```

| Argument or option | Default | What it is |
| --- | --- | --- |
| `PRODUCT` | required | The plasmid the clone should be, such as the `product.dna` a plan wrote |
| `RESULT` | required | A Sanger read as `.ab1`, or a whole-plasmid consensus as FASTA, GenBank or `.dna`. Give one or more |
| `--feature` | none | A feature to judge, in place of the junctions and inserts. Repeat it once per feature |
| `--out` | none | A folder to write the clone's page into, as `verification.html` |

A cloning plan marks each junction in the product it writes, so the command finds the junctions
and the inserts between them on its own. It marks the junction the vector's own bases follow,
so the vector itself is left out. On a plasmid no plan wrote, name what to judge with
`--feature`.

It prints one line for each junction and insert, with its verdict: pass, fail, warn or no
verdict. Any difference outside them is listed with the feature it falls in, and judges nothing.
Then comes a line for each result, and last the clone's line, `verified` or `not verified`.

A consensus can't show a mixed sample: the plasmid read most becomes the consensus, so each one
says so. A result where more than a tenth of the bases it trusts disagree, or fail to line up,
prints one line, that it does not read as this plasmid, and nothing else. The vendor's table of
reads at each base is not read yet, so give the consensus file rather than the folder.

With `--out`, it also writes a page to read the result by, and prints its path last. The page
gives the verdict and the same table, then a map of the plasmid with each result and each
difference on it, and a close-up of every difference. Under a Sanger read's bases it draws the
trace.

The command exits with 1 when the clone is not verified, so a script can stop on it.

## Read designs back

`synbio dmx plan` plans a read-back: clonal stock taken to one well per design, each well marked
so sequencing says which well it came from.

```bash
pixi run synbio dmx plan BUILD --out DIR
```

`BUILD` is a JSON file. It names the run, the sheet of designs, the plate of clonal stock they
are spotted from, which of the two marking routes to use, and the fragment count at or above
which a design is read back.

| Option | Default | What it is |
| --- | --- | --- |
| `--out`, `-o` | required | Folder to write into |

A read-back is one sitting, so the folder gets one page and its data: `protocol.html` and
`protocol.json`. There is no sequence file and no order sheet, because a design is read back as
a name and a count of its pieces rather than as bases.

The two routes are `barcode ligation` and `index PCR`. Barcode ligation uses the lab's own
barcode kit, whose sequences this package does not ship.
[The DMX barcode kit](dmx-barcode-kit.md) says what that file holds.

## Plan a library

`synbio igga plan` plans a barcoded combinatorial library: lists of parts joined round by round,
every part carrying a barcode that names it.

```bash
pixi run synbio igga plan BUILD --out DIR
```

`BUILD` is a JSON file naming the positions, the part and vector files, and what one run
chooses.

| Option | Default | What it is |
| --- | --- | --- |
| `--out`, `-o` | required | Folder to write into |
| `--kind` | `protein` | What the parts FASTA holds: `protein` or `dna` |
| `--site` | none | Feature name, or `START-END`, to put an internal stuffer at, where the vector carries none |
| `--working-site` | none | Feature name, or `START-END`, to put the ccdB cassette at, where the working vector carries none |
| `--pattern` | matches the position's name on a word boundary | Regular expression matching a record name, with `{position}` for a position |
| `--prices` | none | A banded price record you hold, as `.csv`, to cost the bill |
| `--ligase-matrix` | none | A ligase fidelity file you hold, to report how often the ligase joins each round's overhangs |
| `--ligase-sheet` | the only sheet | Which sheet of that file to read, for a workbook holding several |

Without `--prices` the quantities still compute, and every money cell is left as a hole naming
what went unpriced, because nothing here is ever estimated. Nothing is designed on the ligase
file either: it only reports.

It writes the sheets, a record per round, and the run's protocols as a folder:

| File | What it holds |
| --- | --- |
| `parts.tsv` | the synthesis order sheet, every block with its barcode on the same row |
| `barcodes.tsv` | which barcode names which part, and where it sits in the block |
| `changes.tsv` | every amino acid the overhang standard moved, wild type beside synthesised |
| `round-N.dna` | one record per round, up to the last |
| `product.dna` | the last round's product, the construct that stands for the library |
| `block-vector-N.dna` | one block vector per position |
| `working-vector-ccdb.dna` | the working vector the final assembly opens, cassette and all |
| `pool.tsv`, `pool-primers.tsv`, `oligo.dna` | the oligo pool to order, the primers that amplify it, and one oligo drawn as a record |
| `library-read-primers.tsv` | the pairs that read the finished library back |
| `protocol/` | the run as a folder of linked pages, one per protocol |

The block vectors, the pool and the oligo record are written only where the build names a primer
set. Without one, each block is ordered whole and needs no vector to supply its stuffers. The
working vector is written only where the build names one.

## Environment variables

| Variable | What it stands in for |
| --- | --- |
| `LIULAB_MBIO_LIGASE_MATRIX` | `--ligase-matrix`, on `cloning goldengate plan` and on `igga plan` |
| `LIULAB_MBIO_PRICES` | `--prices`, on `igga plan` |

A path typed on the command line wins over the variable. Neither file is shipped, and both are
read from a copy you hold. [Ligation fidelity](ligation-fidelity.md) says what a ligase matrix
is and where one comes from.

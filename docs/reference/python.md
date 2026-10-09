# Python reference

These are the calls a script makes. The list is kept by hand, so it stays short. Anything not
here is read from the source in `src/mbio/` and `src/synbio/`.

## What to import

The root package re-exports only `__version__`. Import a name from the module that owns it.

```python
from mbio.io import read_record
from mbio.cloning.goldengate import plan_assembly
```

A few names are spelled more than once, and each one means something different. The module
path is what tells them apart.

| Name | Resolves to |
| --- | --- |
| `Check` | `mbio.checks.Check` is the judged check, with the value it measured. `mbio.protocol.Check` is how a protocol page shows one |
| `Plan` | Four of them, one per cloning method, each with its own fields. All four are on this page |
| `Build`, `read_build` | `synbio.igga` and `synbio.dmx` each define both, with different fields. Both pairs are on this page |
| `Protocol` | `mbio.protocol.Protocol` is the bench protocol. `synbio.igga.protocols.protocol.Protocol` is one page of a library run |
| `Files` | Six of them, one per pipeline. None is on this page. Read what `write` returns |

## Open and save a record

`read_record` takes a SnapGene `.dna` file, a GenBank file or a FASTA file. `read_region`
pulls one span out of a large FASTA, through the `.fai` index beside it.

::: mbio.io.read_record

::: mbio.io.read_region

::: mbio.snapgene.write_dna

## The record

Coordinates are 0-based and half-open. The first base is 0, and a span's end is the first base
after it. A span that crosses the origin of a circular record ends past the record's length.

::: mbio.sequence.SequenceRecord

::: mbio.sequence.Feature

::: mbio.sequence.Segment

::: mbio.sequence.Primer

::: mbio.sequence.Strand

::: mbio.sequence.reverse_complement

## A check

Most of what the package judges comes back as a list of these. One check holds the name of
what was measured, the measurement, and a status. Where no sourced threshold can judge the
value, the check carries no verdict rather than a pass.

::: mbio.checks.Check

## Plan a cloning job

Four methods, four ways in. Each one returns a `Plan`. Each `Plan.write` puts the product, the
oligo sheet and the protocol into one directory.

One field name collides. `gibson.Plan.product` is the assembly kit: the reaction mix and the
rules that come with it. The plasmid the assembly makes is `gibson.Plan.plasmid`. Every other
method's `.product` is the record.

Two other fields are easy to miss. `restriction.Plan.refusals` says why every other enzyme
pair was refused. `gateway.Plan.write` adds the entry clone as a fifth file, where a BP
reaction was planned.

::: mbio.cloning.goldengate.plan_assembly

::: mbio.cloning.goldengate.Plan
    options:
      members: [product, checks, status, write]

::: mbio.cloning.gibson.plan_gibson

::: mbio.cloning.gibson.Plan
    options:
      members: [plasmid, product, checks, status, write]

::: mbio.cloning.restriction.plan_restriction

::: mbio.cloning.restriction.Plan
    options:
      members: [product, refusals, checks, status, write]

::: mbio.cloning.gateway.plan_gateway

::: mbio.cloning.gateway.Plan
    options:
      members: [product, entry, checks, status, write]

## Design primers

The `mbio primers design` command covers the usual job, and the
[command line reference](cli.md) describes it. What a script adds is control. You can name a
placement by hand, read every ranked pair instead of the winner, judge an oligo you already
hold, and read each check that fired.

`POLYMERASES` is the four polymerases the package ships. Each one carries the buffer its Tm is
computed in, and its own cycling rules.

::: mbio.primers.design_pair

::: mbio.primers.design_primer

::: mbio.primers.ranked_pairs

::: mbio.primers.evaluate_primer

::: mbio.primers.evaluate_pair

::: mbio.primers.Placement

::: mbio.primers.PrimerReport

::: mbio.primers.PairReport

::: mbio.primers.melting_temperature

::: mbio.primers.POLYMERASES

::: mbio.primers.design_pair_on_genome

::: mbio.primers.evaluate_pair_on_genome

::: mbio.primers.GenomeReport

## Design barcodes

`mbio barcode-design` draws a set from the command line. From a script you set the rules
yourself, check a set someone already holds, and ask how far apart the closest two members
stand.

::: mbio.barcodes.BarcodeRules

::: mbio.barcodes.design_barcodes

::: mbio.barcodes.check_barcodes

::: mbio.barcodes.separation

## Write a protein as DNA

Write a protein with the codons your host uses most, then take out any restriction or Type IIS
site you name by synonymous codon change. `codon_tables` says which hosts ship.

::: mbio.translate.optimize_protein

::: mbio.translate.optimize_coding_sequence

::: mbio.translate.CodingSequence

::: mbio.codons.codon_tables

## Enzymes and cut sites

`get_enzyme` finds one enzyme by its name, its trade name or an isoschizomer's name.
`free_enzymes` returns the enzymes that cut none of the records you give it, which is how you
pick an enzyme for a new site.

::: mbio.enzymes.enzymes

::: mbio.enzymes.get_enzyme

::: mbio.enzymes.Enzyme

::: mbio.sites.find_sites

::: mbio.sites.CutSite

::: mbio.sites.free_enzymes

## Draw a map

`draw_map` lays a record out as a circle or a line. `Drawing.write` then writes that layout as
an HTML page that opens offline, as a PNG or as a PDF. The layout is done once, however many
formats you ask for.

::: mbio.plot.draw_map

::: mbio.plot.Drawing

## Edit a protocol

Every pipeline writes `protocol.json` beside the HTML page rendered from it. Read the JSON,
reword a step or drop one, write it back, and render the page again.

::: mbio.protocol.read_protocol

::: mbio.protocol.Protocol

::: mbio.protocol.Step

::: mbio.protocol.write_protocol

::: mbio.protocol.write_html

## Build a library

iGGA is the lab's iterative Golden Gate method. It builds a barcoded library of every
combination of the part lists it is given, one round at a time. `plan_igga` is the way in, and
it reads what one build chooses from a `Build`. `check_library` judges a design that nothing
here planned, by the same rules.

::: synbio.igga.plan_igga

::: synbio.igga.LibraryPlan
    options:
      members: [product, coverage, checks, status, write]

::: synbio.igga.read_build

::: synbio.igga.Build

::: synbio.igga.check_library

## Read wells back

DMX reads back designs the lab already holds. One design sits in one well, and that well is
marked, sequenced and called on its own. The 96 barcode sequences do not ship, so `read_kit`
reads a copy you hold; see [the DMX barcode kit](dmx-barcode-kit.md) for the file it reads.

::: synbio.dmx.plan_dmx

::: synbio.dmx.ReadBackPlan
    options:
      members: [designs, protocol, write]

::: synbio.dmx.read_kit

::: synbio.dmx.Kit

# API reference

Built from the docstrings in `src/liulab_mbio/` and `src/liulab_synbio/`, so this page and the
code cannot drift apart. Write the docstring; this page follows.

## What to import

`liulab_mbio` itself re-exports only `__version__`. Import a name from the module that owns it:

```python
from liulab_mbio.io import read_record
from liulab_mbio.cloning.goldengate import plan_assembly
```

A few names are spelled more than once across the package on purpose. `Check`, `Junction`,
`Part` and `Files` each mean something different in every module that defines one. A
`liulab_mbio.checks.Check` is the judged check, with the value it measured. A
`liulab_mbio.protocol.Check` is how a protocol page shows one. A `Junction` or a `Files` belongs
to the cloning method that defines it. A flat re-export would have to rename one of each, and
would import every dependency the moment you imported the package. So the module path is the
name.

`liulab_mbio.cloning.goldengate` re-exports the pipeline's entry point and its result types.
`design` is **not** re-exported: reach it at `liulab_mbio.cloning.goldengate.design`.

## The examples are tests

`pixi run check` runs the `Examples` blocks in `src/liulab_mbio/`, so an example that no longer
matches its code fails the tests.

Write one where it makes the object easier to use, and leave it out where it would not.
An example nobody keeps up to date is worse than none.

Keep an example cheap, offline and deterministic. It has to give the same answer on any
machine, with no network. A line that cannot do that needs `# doctest: +SKIP` at the end of
that line.

Two things about that marker are easy to get backwards:

- It covers only the line it sits on. It does not carry to the line below. In a block that
  mixes lines that run with lines that cannot, each line that cannot run needs its own.
- A trailing comment written as plain prose looks just like a marker and is not one. Only
  the `# doctest:` form is read as one.

## The sequence model

Everything else reads and writes these. Coordinates are 0-based and half-open, and a span
across the origin of a circular record ends past the record's length — see
[the coordinates decision](adr/0001-coordinates.md).

::: liulab_mbio.sequence

## Checks

::: liulab_mbio.checks

## Files

::: liulab_mbio.io

::: liulab_mbio.snapgene

::: liulab_mbio.edits

## Enzymes, sites, codons and translation

::: liulab_mbio.enzymes

::: liulab_mbio.sites

::: liulab_mbio.codons

::: liulab_mbio.translate

## Barcodes

A barcode names one part, so that reading a product says which part it carries. This module
draws a set whose members stand far enough apart that no two read as one, and checks a set
someone already holds by the same rules. The same seed draws the same set again.

::: liulab_mbio.barcodes

## Overhangs and ligation

Whether two cut ends anneal, the rules a set of Type IIS overhangs is held to, and how well the
set should ligate. Every cloning method reads them here. `ligase` reads a fidelity matrix the
user holds on their own disk: that archive's licence forbids redistribution, so none of it ships
with the package.

::: liulab_mbio.overhangs

::: liulab_mbio.ligase

## Reactions

One tube: the molecules in it, each bound to the role it plays there, and the enzymes acting on
them. A reaction holds no verdict — what must be true inside one is the method's to judge, which
is what lets a gate say which tube a failure belongs to.

::: liulab_mbio.reaction

## Splitting a long cargo

A sequence longer than one synthesised oligo is ordered as several and joined in one pot. This
says where the cuts fall, how far each one could still move, and what overhang each leaves. It
names no method: what the fragments are dressed as to be ordered is
`liulab_mbio.bench.pools`'s, and which enzyme and which primers do the dressing is the caller's.

::: liulab_mbio.split

## Maps

`draw_map` draws a sequence record as a map. `Drawing.write` writes it as one HTML page that
opens offline, or as a PNG or a PDF. The modules under it are its steps. `layers` picks what a
record draws, with its names and colours. `circular` places those items round a circle,
`linear` along a line, and `sequence_view` base by base in rows with the map. `labels` keeps
their labels apart. `fonts` measures text in the faces the page embeds, `svg` writes the shapes,
and `page` wraps them in the page. `convert` turns the shapes into a PNG or a PDF, with every
letter drawn as its outline. `draw_plate` draws a plate's wells the same way, and `plate` lays
the wells out.

::: liulab_mbio.plot
    options:
      members: false

::: liulab_mbio.plot.drawing

::: liulab_mbio.plot.layers

::: liulab_mbio.plot.circular

::: liulab_mbio.plot.linear

::: liulab_mbio.plot.sequence_view

::: liulab_mbio.plot.labels

::: liulab_mbio.plot.fonts

::: liulab_mbio.plot.svg

::: liulab_mbio.plot.page

::: liulab_mbio.plot.convert

::: liulab_mbio.plot.plate

## Primers

Every public name in the modules below imports from `liulab_mbio.primers` too:
`from liulab_mbio.primers import design_pair` works as well as the longer path.

::: liulab_mbio.primers
    options:
      members: false

::: liulab_mbio.primers.polymerase

::: liulab_mbio.primers.thresholds

::: liulab_mbio.primers.placement

::: liulab_mbio.primers.evaluation

::: liulab_mbio.primers.design

::: liulab_mbio.primers.genome

## Protocols

The model a bench protocol is written in, and the renderer that turns one into a single
self-contained HTML page. `Check` is how a page shows a verdict, with no value; `Oligo` is one
row of the order sheet, carrying its own verdict and the checks that fired where something
judged it; `OVERVIEW_CHARS` is the character budget for a header card, and a longer value is
refused rather than truncated. `figures` is the library of named figures a step may show, each
a spec the renderer draws and never a drawing.

::: liulab_mbio.protocol
    options:
      members: false

::: liulab_mbio.protocol.model

::: liulab_mbio.protocol.figures

::: liulab_mbio.protocol.render

## Bench

The numbers any pipeline shares: DNA amounts, the reaction table filled to volume, PCR and
colony PCR, Golden Gate's own enzymes, reaction and cycling, gels, how many colonies a library
round needs and what reading one back a well at a time comes to, the checks that confirm a
clone, heat inactivation, what each material brings with it, a plate and the moves between its
wells, a price record and the bill it makes, the phenotype a clone should show, the primer
order sheet, an oligo pool as a vendor takes it, and the protocol steps any pipeline reuses.

Every public name in the modules below imports from `liulab_mbio.bench` too, `goldengate`,
`coverage` and `pools` excepted: each is imported by module, because its `REFERENCES` are one
chemistry's, one sizing rule's or one vendor order's and not this package's. A module that
cites a source keeps its own `REFERENCES`, and `liulab_mbio.bench.REFERENCES` gathers the ones
re-exported. `steps` is re-exported, and a protocol cites its `DPNI_REFERENCE` and
`PLATE_REFERENCE` only when it runs the step they belong to.

::: liulab_mbio.bench
    options:
      members: false

::: liulab_mbio.bench.amounts

::: liulab_mbio.bench.reactions

::: liulab_mbio.bench.pcr

::: liulab_mbio.bench.goldengate

::: liulab_mbio.bench.gels

::: liulab_mbio.bench.coverage

::: liulab_mbio.bench.readback

::: liulab_mbio.bench.validation

::: liulab_mbio.bench.inactivation

::: liulab_mbio.bench.materials

::: liulab_mbio.bench.plates

::: liulab_mbio.bench.prices

::: liulab_mbio.bench.phenotype

::: liulab_mbio.bench.oligos

::: liulab_mbio.bench.pools

::: liulab_mbio.bench.steps

## Cloning

A cloning method is a package under `liulab_mbio.cloning`, and `plan` is what every method's
plan shares: the files any plan writes, its status, and taking a record already read.

::: liulab_mbio.cloning
    options:
      members: false

::: liulab_mbio.cloning.plan

## Golden Gate

`plan_assembly` is the way in, and `Plan.write` puts the product, the primer sheet and the
protocol in one directory. The inserts are varargs, so `Plan.inserts` is a tuple — plural,
because one reaction joins as many inserts as the overhangs allow. The reaction and the cycling
are not here: a library build runs the same tables, so they sit in `liulab_mbio.bench.goldengate`.

::: liulab_mbio.cloning.goldengate
    options:
      members: false

::: liulab_mbio.cloning.goldengate.plan

::: liulab_mbio.cloning.goldengate.design

::: liulab_mbio.cloning.goldengate.assembly

::: liulab_mbio.cloning.goldengate.oligos

::: liulab_mbio.cloning.goldengate.steps

## Gibson assembly

`plan_gibson` is the way in, and `Plan.write` puts the plasmid, the oligo sheet and the protocol
in one directory. `design` chooses each junction's overlap and lays out the oligos that stitch a
short part or bridge two fragments; `assembly` makes the parts and joins them; `bench` is every
assembly product's own documented numbers, the reaction and the incubation; `oligos` and `steps`
carry what the sheet and the protocol need.

One name means two things in a plan here, so read it carefully. `Plan.product` is the
**assembly product** — the kit the reaction is run with, which sets the overlap rule, the
reaction, the incubation and the fragment count. The plasmid the assembly makes is
`Plan.plasmid`. Golden Gate's `Plan.product` is the record, so the two plans differ here.

::: liulab_mbio.cloning.gibson
    options:
      members: false

::: liulab_mbio.cloning.gibson.plan

::: liulab_mbio.cloning.gibson.design

::: liulab_mbio.cloning.gibson.assembly

::: liulab_mbio.cloning.gibson.bench

::: liulab_mbio.cloning.gibson.oligos

::: liulab_mbio.cloning.gibson.steps

## Restriction and ligation

`plan_restriction` is the way in, and `Plan.write` puts the same four files in one directory.
The second record it takes is the insert, or the plasmid the insert is cut out of, and what that
record carries is what picks the route: a record holding a site for each enzyme is cut, and one
holding neither is amplified with the sites on its primer tails. The modules under it are its
steps. `design` chooses the enzyme pair and says why every other pair was refused, `digest` cuts
a record and says which enzyme left each end, and `amplify` is the route through a PCR.
`ligation` joins the pieces and reads off what each junction now spells. `oligos` says what each
designed oligo is for, `bench` holds this method's own numbers with the source of each, and
`verdicts` the checks a plan carries — including the one nothing sourced can judge, which comes
back with no verdict rather than a pass. `steps` writes the protocol.

Nothing here is re-exported above the method, and only `plan_restriction` and its result types
are re-exported from `liulab_mbio.cloning.restriction` itself.

::: liulab_mbio.cloning.restriction
    options:
      members: false

::: liulab_mbio.cloning.restriction.plan

::: liulab_mbio.cloning.restriction.design

::: liulab_mbio.cloning.restriction.digest

::: liulab_mbio.cloning.restriction.amplify

::: liulab_mbio.cloning.restriction.ligation

::: liulab_mbio.cloning.restriction.oligos

::: liulab_mbio.cloning.restriction.bench

::: liulab_mbio.cloning.restriction.verdicts

::: liulab_mbio.cloning.restriction.steps

## Gateway

`plan_gateway` is the way in, and `Plan.write` puts the expression clone, the oligo sheet and
the protocol in one directory, with the entry clone beside them where a BP reaction was planned.
So a plan writes four files, or five where BP ran. The first record it takes is an entry clone,
or an insert carrying att ends when `donor` is given, or a plain insert when `amplify` is set as
well. Nothing is cut and nothing is ligated here: two att sites recombine, and the reaction
rewrites the sites themselves.

The modules under it are its steps. `att` owns the eight att sequences, which one pairs with
which, and the arithmetic a junction follows. It is not `liulab_mbio.sites`, which means enzyme
cut sites. `design` owns the attB primer tail and the PCR that puts it on an insert.
`recombination` simulates one reaction on two records. `checks` holds the verdicts that span
both reactions — including the one nothing sourced can judge, which comes back with no verdict
rather than a pass. `oligos` says what each designed oligo is for, `bench` holds this method's
own numbers with the source of each, and `steps` writes the protocol.

Two names here repay a second look. A `Junction` is the att site one reaction wrote. Its 25
bases straddle the boundary between the two records that made the product, so where the moved
DNA starts and stops is `Recombination.boundaries` rather than the junction's own span.

::: liulab_mbio.cloning.gateway
    options:
      members: false

::: liulab_mbio.cloning.gateway.plan

::: liulab_mbio.cloning.gateway.att

::: liulab_mbio.cloning.gateway.design

::: liulab_mbio.cloning.gateway.recombination

::: liulab_mbio.cloning.gateway.checks

::: liulab_mbio.cloning.gateway.oligos

::: liulab_mbio.cloning.gateway.bench

::: liulab_mbio.cloning.gateway.steps

## Libraries

`plan_igga` is the way in. It builds a barcoded library of every combination of the part
lists it is given, one round at a time. `LibraryPlan.write` puts the synthesis order sheet, the
barcode and amino-acid change tables, a record for each round, the product and a `protocol`
folder of the run's protocols in one directory. The modules under it are its steps. `method`
holds `IGGA`, the one method, checked when it is imported; `project` holds what one build
chooses, checked as it is read; `gate` judges a finished design reaction by reaction. `standard`
picks the overhang set, `parts` writes each synthesis block and `cargo` orders a block too long
to synthesise as an oligo pool. `vector` takes the destination or retrofits it, `rounds`
simulates each round, and `reads` designs the reads that judge the pool. `bench` turns the
method into amounts, `stages` says what each stage asks of the bench model, and `figures` picks
the two figures only this method needs. `protocols` is one module a protocol of the run, each
owning what its own page prints, and `chain` holds the order they run in.

How many colonies a round needs is not here: any pipeline building a library in rounds counts
them the same way, so that is `liulab_mbio.bench.coverage`.

A library is a pipeline over Golden Gate rather than a cloning method of its own — see
[the library rounds decision](adr/0004-library-rounds.md). It is `liulab_synbio`'s, because it
fixes one method's enzymes, stuffers and round order; everything it builds on is
`liulab_mbio`'s. The method is code and one build's choices are a file — see
[the method decision](adr/0010-method-in-code.md).

::: liulab_synbio.igga
    options:
      members: false

::: liulab_synbio.igga.plan

::: liulab_synbio.igga.method

::: liulab_synbio.igga.project

::: liulab_synbio.igga.gate

::: liulab_synbio.igga.standard

::: liulab_synbio.igga.parts

::: liulab_synbio.igga.cargo

::: liulab_synbio.igga.vector

::: liulab_synbio.igga.rounds

::: liulab_synbio.igga.reads

::: liulab_synbio.igga.bench

::: liulab_synbio.igga.stages

::: liulab_synbio.igga.figures

::: liulab_synbio.igga.chain

::: liulab_synbio.igga.protocols
    options:
      members: false

::: liulab_synbio.igga.protocols.protocol

::: liulab_synbio.igga.protocols.run

::: liulab_synbio.igga.protocols.primer_plates

::: liulab_synbio.igga.protocols.ordering

::: liulab_synbio.igga.protocols.validation

::: liulab_synbio.igga.protocols.creation

::: liulab_synbio.igga.protocols.assembly

::: liulab_synbio.igga.protocols.final

## DMX

DMX is its own protocol for multiplexed validation. A design sits one per well, the well is
marked, sequenced and called on its own, and identity stays with well position throughout. It
takes any cargo, and what `liulab_synbio.igga` builds is one kind of cargo among others, so
`dmx` stands beside `igga` rather than downstream of it: iGGA chains this protocol as any
caller would.

`kit` is the barcode kit a user holds, read from their own copy because the sequences are not
shipped. `method` holds the two marking routes, each route's depth floor, a well's derived
address and the pass rule, with the document each number came from. `steps` writes them up as
protocol steps. `seating` is the method's carrier step — one part a well, nothing pooled at the
end — which is why it is no **round**.

::: liulab_synbio.dmx
    options:
      members: false

::: liulab_synbio.dmx.kit

::: liulab_synbio.dmx.method

::: liulab_synbio.dmx.steps

::: liulab_synbio.dmx.seating

## The command line

The whole module, because typer makes every verb a plain function with a docstring, and
`liulab_mbio.cli:app` and `liulab_synbio.cli:app` — the two objects `[project.scripts]`
registers — are built from them. One cloning method is one sub-app under the `cloning` group,
and every plan verb is the same spine: plan, write, and report what was written.

::: liulab_mbio.cli

::: liulab_mbio.cloning.cli

::: liulab_mbio.cloning.goldengate.cli

::: liulab_mbio.cloning.gibson.cli

::: liulab_mbio.cloning.restriction.cli

::: liulab_mbio.cloning.gateway.cli

::: liulab_synbio.cli

::: liulab_synbio.igga.cli

::: liulab_mbio.protocol.cli

::: liulab_mbio.plot.cli

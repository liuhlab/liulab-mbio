# API reference

Built from the docstrings in `src/liulab_mbio/`, so this page and the code cannot drift apart.
Write the docstring; this page follows.

## What to import

`liulab_mbio` itself re-exports only `__version__`. Import a name from the module that owns it:

```python
from liulab_mbio.io import read_record
from liulab_mbio.goldengate import plan_assembly
```

Two names are spelled twice across the package on purpose — `Check` and `Junction` each mean
something different in the two modules that define them. A
`liulab_mbio.checks.Check` is the judged check, with the value it measured. A
`liulab_mbio.protocol.Check` is how a protocol page shows one. A flat re-export would have to
rename one of each pair, and would import every dependency the moment you imported the package.
So the module path is the name.

`liulab_mbio.goldengate` re-exports the pipeline's entry point and its result types.
`design` and `ligase` are **not** re-exported: reach them at
`liulab_mbio.goldengate.design` and `liulab_mbio.goldengate.ligase`.

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

## Enzymes, sites and codons

::: liulab_mbio.enzymes

::: liulab_mbio.sites

::: liulab_mbio.codons

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

## Protocols

The model a bench protocol is written in, and the renderer that turns one into a single
self-contained HTML page. `Check` is how a page shows a verdict, with no value; `Oligo` is one
row of the order sheet, carrying its own verdict and the checks that fired where something
judged it; `OVERVIEW_CHARS` is the character budget for a header card, and a longer value is
refused rather than truncated.

::: liulab_mbio.protocol
    options:
      members: false

::: liulab_mbio.protocol.model

::: liulab_mbio.protocol.render

## Bench

The numbers any cloning pipeline shares: DNA amounts, PCR and colony PCR, gels, the checks
that confirm a clone, heat inactivation, the phenotype a clone should show, and the primer
order sheet. Every public name in the modules below imports
from `liulab_mbio.bench` too. A module that cites a source keeps its own `REFERENCES`, and
`liulab_mbio.bench.REFERENCES` gathers them all.

::: liulab_mbio.bench
    options:
      members: false

::: liulab_mbio.bench.amounts

::: liulab_mbio.bench.pcr

::: liulab_mbio.bench.gels

::: liulab_mbio.bench.validation

::: liulab_mbio.bench.inactivation

::: liulab_mbio.bench.phenotype

::: liulab_mbio.bench.oligos

## Golden Gate

`plan_assembly` is the way in, and `Plan.write` puts the product, the primer sheet and the
protocol in one directory. The inserts are varargs, so `Plan.inserts` is a tuple — plural,
because one reaction joins as many inserts as the overhangs allow.

::: liulab_mbio.goldengate
    options:
      members: false

::: liulab_mbio.goldengate.plan

::: liulab_mbio.goldengate.design

::: liulab_mbio.goldengate.assembly

::: liulab_mbio.goldengate.bench

::: liulab_mbio.goldengate.ligase

::: liulab_mbio.goldengate.oligos

::: liulab_mbio.goldengate.steps

## The command line

The whole module, because typer makes every verb a plain function with a docstring, and
`liulab_mbio.cli:app` — the object `[project.scripts]` registers — is built from them.

::: liulab_mbio.cli

::: liulab_mbio.goldengate.cli

::: liulab_mbio.protocol.cli

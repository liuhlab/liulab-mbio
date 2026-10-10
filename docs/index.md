# liulab-mbio

Molecular biology design tools for sequences, enzymes, primers and cloning. You give it your
own files; it picks the enzymes, designs the oligos, builds the plasmid you should get, and
writes the bench protocol. One command plans a whole job. Everything it writes is a file you
can open — a map in SnapGene, a sheet in a spreadsheet, one page in a browser.

One install gives two commands: `mbio` for a job with a standard name, and `synbio` for a
method this lab named. [Two commands](two-packages.md) says which one a job wants.

## Pick a job

Each row is one job, with a worked run you can repeat. For a cloning job,
[Choosing a method](methods/index.md) weighs the four against each other.

| What you want | Start here |
| --- | --- |
| Join fragments in one tube, in the order you choose | [Golden Gate assembly](methods/golden-gate.md) |
| Join fragments by bases they share, with no enzyme site anywhere | [Gibson assembly](methods/gibson.md) |
| Cut and paste at a pair of sites both plasmids carry | [Restriction and ligation](methods/restriction-ligation.md) |
| Move a gene between plasmids that carry att sites | [Gateway cloning](methods/gateway.md) |
| A primer pair for a fragment, checked against a whole genome | [Primer design](methods/primers.md) |
| Barcodes far enough apart that no two can be read as one | [Barcode sets](methods/barcodes.md) |
| A protein written as DNA for a host, free of the sites you name | [Codon optimisation](methods/codon-optimisation.md) |
| A plasmid map to look at, or a figure for a paper | [Maps and figures](methods/maps.md) |

The next two are projects rather than single reactions: many plates, and one bench page for
each sitting.

| What you want | Start here |
| --- | --- |
| Every combination of several protein lists, each part barcoded | [iGGA](projects/igga.md) |
| Read back which design landed in each well of a plate | [DMX](projects/dmx.md) |

## Install it

The repo uses [pixi](https://pixi.sh) and nothing else. No pip, no conda, no uv. Clone the
repo, then:

```bash
pixi install
```

That reads `pyproject.toml` and builds the environment from the lock file, so you get the
same versions the tests ran on. Both commands are now on your path.

## Run it

Golden Gate, from a vector file and an insert file that ship with the repo. `--out` names the
folder to write into, and makes it if it is not there:

```bash
pixi run mbio cloning goldengate plan tests/data/pUC19.dna tests/data/GFP.dna --out plan/
```

It prints the design in one line, then the four files it wrote:

```text
pUC19-GFP: 3347 bp, BbsI, 2 fragments, overhangs ATGA, TGGC, fidelity 100% (measured), checks warn
plan/product.dna
plan/primers.tsv
plan/protocol.json
plan/protocol.html
```

| File | What it is |
| --- | --- |
| [product.dna](examples/pUC19-GFP/product.dna) | the finished plasmid, features carried over and both joins marked. Opens in SnapGene |
| [primers.tsv](examples/pUC19-GFP/primers.tsv) | the nine oligos to order, with length and melting temperature |
| [protocol.json](examples/pUC19-GFP/protocol.json) | the same protocol as data |
| [protocol.html](examples/pUC19-GFP/protocol.html) | the protocol as one page: no network, nothing to install |

Each name links to what that run wrote, published here unedited. Nothing in it is typed by
hand, so the same two input files always give the same four. To change the design, change the
command and run it again. To change what the page says, edit `protocol.json` and
[make the page from it again](protocols/editing.md).

[Golden Gate assembly](methods/golden-gate.md) reads that run line by line, and
[Protocols](protocols/index.md) covers the page it wrote.

## Do it from Python

```python
from mbio.cloning.goldengate import plan_assembly

plan = plan_assembly("tests/data/pUC19.dna", "tests/data/GFP.dna")
plan.write("plan/")
```

Each method has one function like that as its way in. Sequence files read into one shared
model, whatever their format:

```python
from mbio.io import read_record

record = read_record("tests/data/pUC19.dna")
```

[Python reference](reference/python.md) has the rest, and [Command line](reference/cli.md)
every verb. The code is at
[liuhlab/liulab-mbio](https://github.com/liuhlab/liulab-mbio), where `pixi run check` runs
the linters, the type checker and the tests.

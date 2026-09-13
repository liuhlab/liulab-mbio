# liulab-mbio

Molecular biology design and analysis tools for DNA sequences, enzymes, primers and cloning.

## Install it

The repo uses [pixi](https://pixi.sh) and nothing else. No pip, no conda, no uv. Clone the
repo, then:

```bash
pixi install
```

That reads `pyproject.toml` and builds the environment from the lock file, so you get the
same versions the tests ran on.

## Use it

Plan a Golden Gate cloning job from a vector file and an insert file:

```bash
pixi run liulab_mbio goldengate plan vector.dna insert.dna --out plan/
```

Four files land in `plan/`. `product.dna` is the assembled plasmid, with its features and
primers marked. `primers.tsv` is the oligos to order. `protocol.json` is the protocol written
as data. `protocol.html` is the page made from that data, one page you can follow at the bench.

If you edit `protocol.json`, turn it back into a page:

```bash
pixi run liulab_mbio protocol render plan/protocol.json
```

The same thing from Python:

```python
from liulab_mbio.goldengate import plan_assembly

plan = plan_assembly("vector.dna", "insert.dna")
plan.write("plan/")
```

Sequence files are read into one shared model, whatever their format:

```python
from liulab_mbio.io import read_record

record = read_record("vector.dna")
```

A primer pair can also be checked against a whole genome, so you learn where else it would
amplify before you order it — see [check primers on a genome](genome-check.md).

[Put GFP into pUC19](golden-gate.md) walks through one job from end to end. The
[API reference](api.md) has the full list, built from the docstrings in `src/`.

## Check your work

One command runs the linters, the type checker and the tests:

```bash
pixi run check
```

It runs every step, then prints all the failures at once. Read to the bottom before you
fix anything.

The docs site is built by a separate command, because it needs a heavier environment:

```bash
pixi run docs-build
```

## Where things live

| Path | What it holds |
| --- | --- |
| `src/liulab_mbio/` | the package |
| `tests/` | the tests |
| `docs/` | this site |
| `scripts/check.sh` | the gate every commit has to pass |
| `CONTEXT.md` | the glossary: the words this repo uses |

Some notes are written for coding agents, not for people. Conventions go under
`docs/agents/`, decision records under `docs/adr/`, and research notes under
`docs/research/`. Nothing in those three directories shows up in the menu or the search
box, and a page written there is still reachable by its own URL.

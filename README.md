# liulab-mbio

Molecular biology design and analysis tools for DNA sequences, enzymes, primers and cloning.

## Set it up

This repo uses [pixi](https://pixi.sh) and nothing else — no pip, no conda, no uv.
Clone it, then:

```bash
pixi install
```

## Use it

Plan a Golden Gate cloning job from a vector file and an insert file:

```bash
pixi run liulab_mbio goldengate plan vector.dna insert.dna --out plan/
```

Four files land in `plan/`. `product.dna` is the assembled plasmid, with its features and
primers marked. `primers.tsv` is the oligos to order. `protocol.json` is the protocol written
as data. `protocol.html` is the page made from that data, one page you can follow at the bench:
reagents, reaction tables, programs, expected bands and troubleshooting.

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
amplify before you order it:
[check primers on a genome](https://liuhlab.github.io/liulab-mbio/genome-check/).

The docs have [a worked example](https://liuhlab.github.io/liulab-mbio/golden-gate/): GFP into
pUC19, with the protocol it writes.

## Check your work

```bash
pixi run check
```

That runs the linters, the type checker and the tests. It reports every failure at once, so read to
the bottom before you fix anything.

## Read the docs

The site is at <https://liuhlab.github.io/liulab-mbio/>.
Build it yourself with `pixi run docs-build`.

## Set up your agent

Skills for coding agents live in `skills/`. Link them into each agent's own folder:

```bash
python skills/install.py --target all
```

If you work on the lab's clusters, add the shared plugin once per machine:

```text
/plugin marketplace add liuhlab/liulab-compute-skills
/plugin install lab-compute@liulab
```

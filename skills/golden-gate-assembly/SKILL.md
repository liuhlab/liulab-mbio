---
name: golden-gate-assembly
description: >-
  Plan a Golden Gate cloning experiment end to end from a vector and an insert sequence file
  (SnapGene .dna, GenBank or FASTA): pick a Type IIS enzyme with no site in the parts, design
  and check the PCR primers, simulate the assembly, design the colony PCR that validates it,
  and write an interactive HTML bench protocol with expected results. Use when the user wants
  to clone, insert or subclone a sequence into a plasmid by Golden Gate, asks for primers with
  BsaI, BsmBI, BbsI, PaqCI or SapI tails, or wants a ready-to-run cloning protocol from two
  sequence files.
---

# Golden Gate assembly

`liulab_mbio` does the design. This skill is the way in: one command turns two sequence files
into a product map, a primer order sheet and a bench protocol. Do not hand-design primers,
overhangs or band sizes beside it — the package computes them from the sequences, and a number
written by hand is one that nothing checks.

## Run it

```bash
pixi run liulab_mbio goldengate plan vector.dna insert.dna --out plan/
```

It chooses the enzyme, designs the overhangs, simulates both PCRs and the ligation, works out
the bench quantities, and designs the colony PCR and sequencing that confirm the clone. Three
files land in the directory you name:

- `product.dna` — the assembled plasmid, features carried over and each junction annotated
- `primers.tsv` — every oligo it designed, with length and Tm
- `protocol.html` — one self-contained page: reagents, reaction tables, thermocycler programs,
  expected bands, a simulated gel and troubleshooting, step by step

The command prints a summary line and the three paths. The same inputs write the same bytes, so
a protocol can be regenerated rather than edited.

## Ask only for what the files do not settle

Every option below has a default; a vector and an insert are enough on their own.

| Option | What it sets |
| --- | --- |
| `--site` | A feature name, or `START-END`, that the insert replaces. Defaults to the vector's own `MCS` feature |
| `--orientation` | `forward` or `reverse` |
| `--in-frame` | Hold the insert's junction on a codon boundary, for a fusion |
| `--enzyme` | Force one. Refused when it reads a site in either part |
| `--polymerase` | For the two PCRs |
| `--host` | The strain the protocol names |
| `--name` | What to call the product |

## From Python

When you need the numbers rather than the files:

```python
from liulab_mbio.goldengate import plan_assembly

plan = plan_assembly("vector.dna", "insert.dna")
plan.status  # "pass", "warn" or "fail" over every check and every primer
plan.write("plan/")  # the same three files
```

`plan.assembly.checks` carries the product's verdicts, `plan.colony.clones` the expected bands
per clone, and `plan.phenotype` what the product says about itself — what drives the insert,
whether anything should be translated, and how a plate reads.

## Before you hand it over

Open `protocol.html` and read it back, as `build-protocol` asks. Then tell the user what the
design chose and what it cost: the enzyme, the overhangs and their fidelity, the product
length, and the expected colony PCR bands. Say which checks warned rather than hiding them, and
say whether a fidelity number was measured or estimated — the report carries that flag, and the
two are not the same kind of number.

## When it refuses

It raises rather than guessing, and the message names the cause.

- **No enzyme is free.** Every candidate reads a site in the parts. Taking one out changes what
  a part spells, so it is the user's call; the ranking reports which sites a synonymous codon
  change could reach and which lie outside a coding sequence.
- **A junction has no overhang left.** Every candidate was refused; the message names the rule.
  A wider junction window lets the vector-side junction slide.
- **The DNA does not fit the reaction.** Concentrate the fragments, or scale the reaction up.
- **The parts do not chain.** The overhangs do not close the circle.

## The protocol side

`build-protocol` owns the page: the model, the data format and how to render one. Read it when
you need to change what a protocol says, rather than what the design computes.

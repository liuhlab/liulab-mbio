---
name: golden-gate-assembly
description: >-
  Plan a Golden Gate cloning experiment end to end from a vector and one or more insert
  sequence files (SnapGene .dna, GenBank or FASTA): pick a Type IIS enzyme with no site in the
  parts, design and score the whole overhang set, check the PCR primers, simulate the assembly,
  design the colony PCR that validates every junction, and write an interactive HTML bench
  protocol with expected results. Use when the user wants to clone, insert, subclone or join
  sequences into a plasmid by Golden Gate, asks for primers with BsaI, BsmBI, BbsI, PaqCI or
  SapI tails, or wants a ready-to-run cloning protocol from sequence files.
---

# Golden Gate assembly

`liulab_mbio` does the design. This skill is the way in: one command turns a vector and its
inserts into a product map, a primer order sheet and a bench protocol. Never invent a primer, an
overhang, a band size or any other number the package computes: it works each one out from the
sequences and checks it, and nothing checks a number you made up.

## Run it

```bash
pixi run liulab_mbio goldengate plan vector.dna insert.dna --out plan/
pixi run liulab_mbio goldengate plan vector.dna first.dna second.dna third.dna --out plan/
```

One reaction joins as many inserts as the overhangs allow. Give them in the order they go round
the product: the vector is opened across the span they replace, and each insert follows the one
before it. A single insert is the ordinary case, not a special one.

It chooses the enzyme, designs every junction's overhang together, simulates the PCRs and the
ligation, works out the bench quantities, and designs the colony PCR and sequencing that confirm
the clone. Four files land in the directory you name:

- `product.dna` — the assembled plasmid, features carried over and each junction annotated
- `primers.tsv` — every oligo it designed, with length and Tm; it stays whole when a step
  leaves the page
- `protocol.json` — the protocol as data: a draft you may edit through `build-protocol`
- `protocol.html` — the page rendered from `protocol.json`, self-contained: reagents, reaction
  tables, thermocycler programs, expected bands, a simulated gel and troubleshooting, step by step

The command prints a summary line and the four paths. The same inputs write the same bytes.

## What the fragment count changes

- **n inserts means n + 1 junctions**, and the overhangs are chosen as one set: all different,
  none the reverse complement of another, and scored together.
- **Fidelity falls as junctions are added.** The report and the protocol print the number and
  say whether it was measured or estimated.
- **The reaction and the cycling follow the count.** NEB's tables step up at three fragments, at
  seven and at fourteen, and the amounts cover every part.
- **Validation covers every junction**: a colony PCR band that reads each one, and a
  wrong-orientation lane for each insert.
- **The steps stay per experiment.** More parts add one PCR step each, not a protocol per pair.

## Ask only for what the files do not settle

Every option below has a default; a vector and one insert are enough on their own.

| Option | What it sets |
| --- | --- |
| `--site` | A feature name, or `START-END`, that the inserts replace. Defaults to the vector's own `MCS` feature |
| `--orientation` | `forward` or `reverse`. Give it once for all the inserts, or once per insert, in order |
| `--in-frame` | Hold every insert's junction on a codon boundary, for a fusion |
| `--enzyme` | Force one. Refused when it reads a site in any part |
| `--polymerase` | For the PCRs |
| `--host` | The strain the protocol names |
| `--name` | What to call the product |

## From Python

When you need the numbers rather than the files:

```python
from liulab_mbio.goldengate import plan_assembly

plan = plan_assembly("vector.dna", "first.dna", "second.dna")
plan.status  # "pass", "warn" or "fail" over every check and every primer
plan.write("plan/")  # the same four files
```

`plan.assembly.checks` carries one verdict per part beside the product's own,
`plan.overhangs.fidelity` the score of the set, `plan.colony.clones` the expected bands per
candidate clone, and `plan.phenotype` what the product says about itself — what drives the
inserts, whether anything should be translated, and how a plate reads.

A value the user gives goes in here. For another amount of DNA, replace `plan.amounts` with
`liulab_mbio.goldengate.bench.assembly_amounts` at that amount, using `dataclasses.replace`,
and write the plan again.

## Before you hand it over

Open `protocol.html` and read it back, as `build-protocol` asks. Then tell the user what the
design chose and what it cost: the enzyme, the overhang set and its fidelity, the product
length, and the expected colony PCR bands. Say which checks warned rather than hiding them, and
say whether a fidelity number was measured or estimated — the report carries that flag, and the
two are not the same kind of number.

## When it refuses

It raises rather than guessing, and the message names the cause.

- **No enzyme is free.** Every candidate reads a site in the parts. Taking one out changes what
  a part spells, so it is the user's call; the ranking reports which sites a synonymous codon
  change could reach and which lie outside a coding sequence.
- **A junction has no overhang left.** Every candidate was refused and the message names the
  rule. More junctions leave fewer overhangs, so this is likeliest on a long assembly. A wider
  junction window lets the vector-side junction slide; a junction between two inserts cannot
  move, neither part being able to spell the other's bases.
- **The DNA does not fit the reaction.** Concentrate the fragments, or scale the reaction up.
- **The parts do not chain.** The overhangs do not close the circle.

## The protocol side

`build-protocol` owns the page and the rules for changing it. When the user asks for what the
plan did not anticipate — skip the DpnI digest, add a gel purification, add a caution — edit
`protocol.json` as it says and render it again.

A change to what the plan computes — the enzyme, the polymerase, the inserts or their
orientation — goes back through `goldengate plan`, which writes `protocol.json` afresh. Within
the session, reapply the changes the user asked for to the new protocol.

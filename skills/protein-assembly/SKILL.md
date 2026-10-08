---
name: protein-assembly
description: >-
  Build a barcoded combinatorial library from lists of protein or DNA sequences by iterative
  Golden Gate with `liulab_synbio`: design every part's synthesis block, choose the overhang
  standard that costs the proteins fewest amino acids, draw a barcode for each part, accept or
  retrofit the destination vector, simulate every round, size each round's colonies for the
  completeness asked for, and write an interactive HTML bench protocol covering all the rounds. Use
  when someone wants every combination of two or more lists of protein domains as a real plasmid
  library, a combinatorial or multiplexed assembly in rounds, a synthesis order sheet for a
  domain library, or asks how many colonies a library round needs.
---

# Protein library assembly

`liulab_synbio.igga` does the design. This skill is the way in: one command turns a project
file naming lists of proteins and a vector into a synthesis order sheet, annotated records and a
bench protocol. Never invent a block, a barcode, an overhang, an amount or a colony count: the package
works each one out and checks it, and nothing checks a number you made up.

## Run it

```bash
pixi run liulab_synbio igga plan project.json --out library/
```

`project.json` is what the user writes: `name`, `positions`, `parts` and `vector` by path,
`host`, `oligo_length`, `batch_size`, `completeness`, and optionally `seed`, `reserved_extra`,
`barcode`, `primers` and `bands` for the pool, and `validate_from` with `route` for the read
back. It is checked where it is read, so a bad value fails before anything is designed.
Copy [the AP-1 project](../../docs/examples/ap1-library/project.json), a whole run with its
inputs and outputs beside it.

The method itself is code, not a file, and `docs/adr/0010-method-in-code.md` says why: its
overhangs, enzymes and stuffers are the DNA of molecules already on the shelf. Never ask the user
for them and never write them into a project.

`parts.fasta` holds every part list in one file. A record's name says which position it fills —
`N_ATF2`, `bZIP_JUN`, `C_VP64` — and a name that says no position, or two of them, is refused
naming it. Pass `--kind dna` where the sequences are already coded: those codons are checked and
kept, not written again.

`reserved_extra` names any further enzyme a step outside the rounds cuts the cargo with — seating
a part in a carrier, or a last transfer into a working vector. It **adds** to the method's own
list and never replaces it, and every block is held clear of the union, its stuffers included.
Ask the user what cuts their cargo outside the rounds, and leave it out where nothing does.

Files land in the directory you name:

- `parts.tsv` — the synthesis order sheet: every block 5' to 3', with its barcode on the same row
- `barcodes.tsv` — which barcode names which part, and its slot in the finished block
- `changes.tsv` — every amino acid the overhang standard moved, wild type beside synthesised
- `pool.tsv`, `pool-primers.tsv` — the oligo pool the blocks are built from, and the primers that
  amplify it. Both appear where the project names a primer set, and the blocks are then not
  ordered at all
- `round-1.dna` … `product.dna` — one annotated record a round, the last the whole construct
- `protocol/` — the run as a chain of protocols: `project.json` is the data, a draft you may
  edit through `build-protocol`, and `index.html` is the way into the pages rendered from it

The command prints a summary line and the paths. The same inputs write the same bytes.

## Price, and the figures

`--prices prices.csv` costs the bill from a record the user holds; the package ships none, and
every row nobody priced carries a hole rather than a guess. Ask the user for their record, never
for a number. Pass `--prices` only where there is one: the quantities compute without it.

Draw each record the run wrote with `plot-map`, which owns every drawing question. The gel and
the plate layout are already in the pages, drawn beside the step that uses them.

## What the list count and the round count change

- **One round a position.** Each round appends one part list to every member at once, so three
  lists of 24 are three rounds and 13,824 constructs, not 13,824 syntheses.
- **Constructs multiply, and so does the colony count.** `completeness` has no default: the
  chance a round may miss a member is the user's call. Every round is sized separately, because a
  round short of its floor loses members no later round can put back.
- **The overhang standard is charged to the proteins.** One standard serves the whole library, so
  a junction forces terminal residues on every member either side. `changes.tsv` is what the user
  is paying for — show it to them before they order.
- **More positions means more overhangs from one set**, so a long project is likelier to refuse.

## The vector decides one overhang

- **A vector already carrying an internal stuffer is taken as it stands**, and position one's
  entry overhang is pinned to the one that stuffer spells. DNA that exists cannot be re-chosen.
- **A vector carrying none is retrofitted** at the `--site` you name, and the stuffer put in
  carries the overhang the design chose. `plan.destination.edit` says what the insertion changed;
  report it as an edit rather than handing over a new file silently.

## From Python

```python
from liulab_synbio.igga.plan import plan_igga

plan = plan_igga("project.json")
plan.status  # "pass", "warn" or "fail" over every round's checks
plan.write("library/")
```

`plan.standard.changes` is the amino-acid bookkeeping, `plan.parts` every synthesis block,
`plan.rounds` each round's simulated product and its verdicts, `plan.coverage` the colonies each
round needs, and `plan.bench` what each round takes in the tube.

## Before you hand it over

Open `protocol/index.html` and read the pages back, as `build-protocol` asks. Then tell the user what the
design chose and what it cost: the entry overhangs, the amino acids changed, the construct count,
the colonies each round needs, and whether their vector was retrofitted. Say which checks warned
rather than hiding them.

## When it refuses

It raises rather than guessing, and the message names the cause.

- **A record's name says no position, or says two.** Rename it, or pass `--pattern`.
- **No overhang standard fits.** Every candidate was refused at some junction and the message
  names the rule. Fewer positions is the fix.
- **A part spells a site no synonymous change can remove.** The part would cut itself; the user
  decides whether to change the protein.
- **A part list is larger than its barcodes allow.** A longer barcode, or a dial turned off —
  `barcode-design` owns that call.
- **The vector reads an enzyme site outside its stuffer.** A round would cut the backbone.

## The neighbouring skills

`barcode-design` owns barcode rules and why there is no GC band. `codon-optimize` owns writing or
checking one coding sequence for a host. `primer-design` owns every primer question this plan does not
already answer; the plan designs the two library read pairs against the simulated record and
writes them with the other sheets. `build-protocol` owns the page:
when the user wants a step this plan did not anticipate, edit `protocol/project.json` and
render the folder again. A change to what the plan computes — the project's vector, host or
completeness — goes back through `igga plan`, which writes the folder afresh.

---
name: protein-assembly
description: >-
  Build a barcoded combinatorial library from lists of protein or DNA sequences by iterative
  Golden Gate with `liulab_mbio`: design every part's synthesis block, choose the overhang
  standard that costs the proteins fewest amino acids, draw a barcode for each part, accept or
  retrofit the destination vector, simulate every round, size each round's colonies for the
  coverage asked for, and write an interactive HTML bench protocol covering all the rounds. Use
  when someone wants every combination of two or more lists of protein domains as a real plasmid
  library, a combinatorial or multiplexed assembly in rounds, a synthesis order sheet for a
  domain library, or asks how many colonies a library round needs.
---

# Protein library assembly

`liulab_mbio.library` does the design. This skill is the way in: one command turns lists of
proteins, a scheme and a vector into a synthesis order sheet, annotated records and a bench
protocol. Never invent a block, a barcode, an overhang, an amount or a colony count: the package
works each one out and checks it, and nothing checks a number you made up.

## Run it

```bash
pixi run liulab_mbio library plan parts.fasta \
  --scheme scheme.json --vector vector.dna \
  --host human --coverage 10 --out library/
```

`parts.fasta` holds every part list in one file. A record's name says which position it fills —
`N_ATF2`, `bZIP_JUN`, `C_VP64` — and a name that says no position, or two of them, is refused
naming it. Pass `--kind dna` where the sequences are already coded: those codons are checked and
kept, not written again.

The scheme is data the user supplies, not something this package ships;
`docs/examples/protein-library/scheme.json` is the worked example to copy and edit, and
`docs/adr/0005-scheme-data.md` says why. Read it back to the user when they ask what the design
assumed.

Files land in the directory you name:

- `parts.tsv` — the synthesis order sheet: every block 5' to 3', with its barcode on the same row
- `barcodes.tsv` — which barcode names which part, and its slot in the finished block
- `changes.tsv` — every amino acid the overhang standard moved, wild type beside synthesised
- `round-1.dna` … `product.dna` — one annotated record a round, the last the whole construct
- `protocol.json` — the protocol as data: a draft you may edit through `build-protocol`
- `protocol.html` — the page rendered from it, covering every round as one experiment

The command prints a summary line and the paths. The same inputs write the same bytes.

## What the list count and the round count change

- **One round a position.** Each round appends one part list to every member at once, so three
  lists of 24 are three rounds and 13,824 constructs, not 13,824 syntheses.
- **Constructs multiply, and so does the colony count.** `--coverage` has no default: how much of
  a library a round may lose is the user's call. Every round is sized separately, because a round
  short of its coverage loses members no later round can put back.
- **The overhang standard is charged to the proteins.** One standard serves the whole library, so
  a junction forces terminal residues on every member either side. `changes.tsv` is what the user
  is paying for — show it to them before they order.
- **More positions means more overhangs from one set**, so a long scheme is likelier to refuse.

## The vector decides one overhang

- **A vector already carrying an internal stuffer is taken as it stands**, and position one's
  entry overhang is pinned to the one that stuffer spells. DNA that exists cannot be re-chosen.
- **A vector carrying none is retrofitted** at the `--site` you name, and the stuffer put in
  carries the overhang the design chose. `plan.destination.edit` says what the insertion changed;
  report it as an edit rather than handing over a new file silently.

## From Python

```python
from liulab_mbio.library.plan import plan_library

plan = plan_library("parts.fasta", "scheme.json", "vector.dna", host="human", coverage=10)
plan.status  # "pass", "warn" or "fail" over every round's checks
plan.write("library/")
```

`plan.standard.changes` is the amino-acid bookkeeping, `plan.parts` every synthesis block,
`plan.rounds` each round's simulated product and its verdicts, `plan.coverage` the colonies each
round needs, and `plan.bench` what each round takes in the tube.

## Before you hand it over

Open `protocol.html` and read it back, as `build-protocol` asks. Then tell the user what the
design chose and what it cost: the entry overhangs, the amino acids changed, the construct count,
the colonies each round needs, and whether their vector was retrofitted. Say which checks warned
rather than hiding them.

## When it refuses

It raises rather than guessing, and the message names the cause.

- **A record's name says no position, or says two.** Rename it, or pass `--pattern`.
- **No overhang standard fits.** Every candidate was refused at some junction and the message
  names the rule. Fewer positions, or a scheme with a different scar, is the fix.
- **A part spells a site no synonymous change can remove.** The part would cut itself; the user
  decides whether to change the protein.
- **A part list is larger than its barcodes allow.** A longer barcode, or a dial turned off —
  `barcode-design` owns that call.
- **The vector reads an enzyme site outside its stuffer.** A round would cut the backbone.

## The neighbouring skills

`barcode-design` owns barcode rules and why there is no GC band. `codon-optimize` owns writing or
checking one coding sequence for a host. `primer-design` owns every primer question, including
the primers that read the barcode block — this plan designs none. `build-protocol` owns the page:
when the user wants a step this plan did not anticipate, edit `protocol.json` and render again.
A change to what the plan computes — the scheme, the vector, the host, the coverage — goes back
through `library plan`, which writes `protocol.json` afresh.

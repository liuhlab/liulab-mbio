# Gateway cloning

Read this once `SKILL.md` has chosen Gateway.

`liulab_mbio` does the design, and one command turns an insert and the vectors into clone maps,
an oligo order sheet and a bench protocol. Never invent an att sequence, a primer tail, a DNA
amount, an incubation or a band size: the package works each one out from the sequences and the
vendor's own documented numbers, and nothing checks a number you made up.

Nothing is cut and nothing is ligated. Two att sites recombine and the reaction rewrites the
sites themselves, so there is no enzyme to choose and no overhang to score.

## Run it

Three routes, and the records you have decide which:

```bash
pixi run liulab_mbio cloning gateway plan entry.dna destination.dna --out plan/
pixi run liulab_mbio cloning gateway plan insert.dna destination.dna --donor donor.dna --out plan/
pixi run liulab_mbio cloning gateway plan gene.dna destination.dna --donor donor.dna --amplify --out plan/
```

An entry clone alone plans LR. A donor vector as well plans BP first, feeding its entry clone
into LR; that first record is then an attB-flanked fragment, which may be linear. `--amplify`
takes a plain insert carrying no att site and amplifies it onto attB ends first.

**Both vectors are the user's own files.** This package ships no vector catalogue, and it
recognises a vector by searching its bases for att sites rather than by matching a shipped
sequence — see `docs/adr/0008-att-sites-and-vectors.md`.

It designs the attB tails where an insert needs them, simulates both recombinations base for
base, works out what each reaction takes, and designs the colony PCR and the sequencing that
confirm the clone. Four files land in the directory you name, and five where BP ran:

- `entry-clone.dna` — the entry clone BP makes, where BP was planned and not otherwise
- `product.dna` — the expression clone, features carried over and each att junction annotated
- `primers.tsv` — every oligo it designed, with length and Tm
- `protocol.json` — the protocol as data: a draft you may edit through `build-protocol`
- `protocol.html` — the page rendered from it, self-contained: both reaction tables, the
  incubations, expected bands, a simulated gel and troubleshooting

The command prints a summary line and those paths. The same inputs write the same bytes.

## Ask only for what the files do not settle

| Option | What it sets |
| --- | --- |
| `--donor` | The donor vector, which plans BP first. Leave it off for an entry clone |
| `--amplify` | Design attB primers for an insert carrying no att site. Needs `--donor` |
| `--fusion` | `none`, `N-terminal`, `C-terminal` or `both`: which tag the insert reads into |
| `--polymerase` | For the attB PCR |
| `--host` | The strain the protocol names for selecting each clone |
| `--name` | What to call the expression clone |

## The junction is not scarless

A whole att site stands at each end of the insert, so the product gains 25 bases there. Where a
fusion reads through one, those bases are translated. The frame is judged only at an end
`--fusion` names, because no record says which coding sequence is a tag.

An attB tail is four G residues, the 25 bp site and the frame bases that end needs: 29 bases,
31 at the N terminus and 30 at the C. The G residues leave with the BP by-product, so the entry
clone reads the same whichever route made it.

## Both strains matter, and they are not the same one

A donor and a destination vector carry ccdB, which kills an ordinary strain, so each is grown in
a ccdB-resistant strain the protocol names. The clone is then selected in a strain ccdB kills,
which is what stops unreacted vector growing. Never one carrying F′: its `ccdA` cancels the
selection, and the plan fails that check by name.

## What carries no verdict

A check nothing sourced measures says so instead of passing. Here that is whether a reaction
plated on a ccdB-resistant strain still counter-selects: the manuals require that strain for
growing the vectors and forbid F′, and say nothing about this. The page shows it as `not judged`
with the guidance beside it. Read that as unmeasured, never as fine.

## From Python

```python
from liulab_mbio.cloning.gateway import plan_gateway

plan = plan_gateway("entry.dna", "destination.dna")
plan.status  # "pass", "warn" or "fail" over every check and every oligo
plan.write("plan/")  # the same files
```

`plan.route` is which reactions were planned, `plan.entry` and `plan.product` the two clones,
`plan.junctions` what each att junction now spells and where, `plan.checks` every verdict,
`plan.colony.clones` the expected bands per candidate colony, and `plan.lr.phenotype` how the
plate should read.

## When it refuses

It raises rather than guessing, and the message names the cause: a record carrying none of the
att sites its side of a reaction needs, naming every site looked for; a vector handed in linear;
`--amplify` with no donor vector; an insert end no annealing region fits; or DNA that does not
fit a reaction.

What it plans and fails instead of refusing: an att site inside the insert, which would
recombine where nobody meant it to; an F′ host; two clones a plate cannot tell apart; and a
fusion frame broken at either end. Read those off `plan.checks` and say so.

## Variants and quirks

`variants.md` in this directory holds the vendor routes this command does not plan — the
one-tube protocol, the two-step adapter PCR, and the entry vectors whose own sites drift. Read
it only when one of them comes up.

## Before you hand it over

Open `protocol.html` and read it back, as `build-protocol` asks. Then tell the user the route,
what each junction spells, the product length, which strain grows each vector and which selects
the clone, and the expected colony PCR bands. Say which checks warned and which nothing judged,
rather than hiding either.

## The primer side

`primer-design` owns primer questions: which placement a task calls for, how a preference such
as a maximum length goes in, and what a warning left on the sheet means. An attB primer is a
fixed tail over an annealing region, so its Tm is read from the region that anneals while its
structure is judged on the whole oligo.

## The protocol side

`build-protocol` owns the page and the rules for changing it. When the user asks for what the
plan did not anticipate — skip a step, add a cleanup, add a caution — edit `protocol.json` as it
says and render it again.

A change to what the plan computes — the vectors, the fusion, the polymerase or the route —
goes back through `cloning gateway plan`, which writes `protocol.json` afresh. Within the
session, reapply the changes the user asked for to the new protocol.

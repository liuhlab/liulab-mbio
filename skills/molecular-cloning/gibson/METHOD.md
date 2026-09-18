# Gibson assembly

Read this once `SKILL.md` has chosen Gibson.

`liulab_mbio` does the design, and one command turns a vector and its inserts into a product
map, an oligo order sheet and a bench protocol. Never invent an overlap, a primer, a volume, an
incubation or a band size: the package works each one out from the sequences and the kit's own
documented numbers, and nothing checks a number you made up.

## Run it

```bash
pixi run liulab_mbio cloning gibson plan vector.dna insert.dna --out plan/
pixi run liulab_mbio cloning gibson plan vector.dna first.dna second.dna third.dna --out plan/
```

Give the inserts in the order they go round the product. The vector is opened by PCR across the
span they replace; a vector handed in already linear is taken as it is, and then names no site.

It chooses the overlap at every junction, designs the primers that carry it as a tail, simulates
the product, sizes the one assembly reaction, and designs the colony PCR and the sequencing that
confirm the clone. Four files land in the directory you name:

- `product.dna` — the assembled plasmid, features carried over and each junction annotated
- `primers.tsv` — every oligo it designed, with length and Tm; it stays whole when a step
  leaves the page
- `protocol.json` — the protocol as data: a draft you may edit through `build-protocol`
- `protocol.html` — the page rendered from `protocol.json`, self-contained: reagents, the
  reaction table, the incubation, expected bands, a simulated gel and troubleshooting

The command prints a summary line and the four paths. The same inputs write the same bytes.

## Ask only for what the files do not settle

Every option has a default; a vector and one insert are enough on their own.

| Option | What it sets |
| --- | --- |
| `--site` | A feature name, or `START-END`, that the inserts replace. Defaults to the vector's own `MCS` feature. Refused for a vector already linear |
| `--orientation` | `forward` or `reverse`. Give it once for all the inserts, or once per insert, in order |
| `--route` | `amplify` or `stitch`, once for all or once per insert |
| `--bridge` | A junction joined by one oligo, written `BEFORE:AFTER` as the two parts are named |
| `--product` | The kit on your bench. The start of a name is enough: `nebuilder`, `gibson`, `in-fusion` |
| `--polymerase` | For the PCRs |
| `--host` | The strain the protocol names |
| `--name` | What to call the product |

## The kit decides the numbers

NEBuilder HiFi is assumed. Each kit carries its own overlap rule, reaction, incubation, molar
ratio and documented fragment count, so naming yours changes the design and not just a label.

| `--product` | What is its own |
| --- | --- |
| `nebuilder` | The default |
| `gibson` | A longer overlap at high fragment counts, and more insert at low ones |
| `in-fusion` | Takara's rules throughout: a shorter overlap, no melting temperature on it, a smaller reaction, one incubation, and neither oligo route |

The page's title names the method for any kit. The summary, the materials and the reaction step
name the kit, so say which one the plan used when you hand it over.

## Two oligo routes

- **`--route stitch`** builds a short insert out of overlapping oligos tiling both strands,
  assembled in the same reaction rather than amplified. It suits a linker, a tag or a short
  promoter, and it takes no PCR and no synthesis order. Two ceilings, and they differ: above
  **500 bases** the plan refuses, because the note allows only twelve 60-base oligos; outside
  **60 to 150 bp** it warns, because that is the window Addgene recommends the route between.
- **`--bridge`** joins two fragments sharing no homology with one oligo carrying homology to
  each end, so neither fragment needs a tailed primer and an amplicon made for something else
  goes in as it is.

Both routes' oligos are rows of the one sheet, with the purpose each serves. Neither primes
anything, so no threshold judges them and each row carries no verdict rather than a pass.

## What carries no verdict

A check nothing sourced measures says so instead of passing: an overlap's GC, and how alike two
overlaps are. Under In-Fusion, the overlap's melting temperature and the assembly DNA join them,
because Takara publishes neither. Read "not judged" as unmeasured, never as fine.

## From Python

```python
from liulab_mbio.cloning.gibson import plan_gibson

plan = plan_gibson("vector.dna", "first.dna", "second.dna")
plan.status  # "pass", "warn" or "fail" over every check and every oligo
plan.write("plan/")  # the same four files
```

`plan.product` is the kit; `plan.plasmid` is the assembled plasmid. `plan.overlaps` gives the
bases at each junction, `plan.checks` every verdict, `plan.colony.clones` the expected bands per
candidate clone, and `plan.phenotype` how the plate should read.

For an amount of vector the user holds, replace `plan.amounts` with
`liulab_mbio.cloning.gibson.bench.assembly_amounts` at their own `vector_ng`, using
`dataclasses.replace`, and write the plan again.

## Before you hand it over

Open `protocol.html` and read it back, as `build-protocol` asks. Then tell the user the kit, the
overlap at each junction, the product length and the expected colony PCR bands. Say which checks
warned and which nothing judged, rather than hiding either.

## When it refuses

It raises rather than guessing, and the message names the cause.

- **A part is too long to stitch**, or too short to take an overlap from. The message names the
  ceiling and what to do instead.
- **The kit takes no single-stranded oligo**, so it can neither stitch nor bridge. The message
  names the kits that can.
- **A bridge names a junction this assembly has not got.** Use the two part names the plan uses.
- **No insertion site.** A circular vector annotates no `MCS` feature and none was named.
- **No annealing region fits a part's end**, so no primer can carry that overlap.

## The primer side

`primer-design` owns primer questions: which placement a task calls for, how a preference such
as a maximum length goes in, and what a warning left on the sheet means. A Gibson primer is a
tail plus an annealing region, so its Tm is read from the part that anneals while its hairpins
and dimers are judged on the whole oligo.

## The protocol side

`build-protocol` owns the page and the rules for changing it. When the user asks for what the
plan did not anticipate — skip the DpnI digest, add a gel purification, add a caution — edit
`protocol.json` as it says and render it again.

A change to what the plan computes — the kit, the polymerase, the inserts, their orientation or
their route — goes back through `cloning gibson plan`, which writes `protocol.json` afresh.
Within the session, reapply the changes the user asked for to the new protocol.

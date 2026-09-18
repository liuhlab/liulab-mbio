# Restriction and ligation

Read this once `SKILL.md` has chosen restriction and ligation.

`liulab_mbio` does the design, and one command turns a vector and an insert into a product map,
an oligo order sheet and a bench protocol. Never invent an enzyme pair, a primer tail, a band
size or a ligation amount: the package works each one out from the sequences and checks it, and
nothing checks a number you made up.

## Run it

```bash
pixi run liulab_mbio cloning restriction plan vector.dna insert.dna --out plan/
pixi run liulab_mbio cloning restriction plan vector.dna source.dna --out plan/ --enzyme EcoRI --enzyme BamHI
```

The second argument is the insert, or the plasmid it is cut out of. Name no enzyme and the pair
is chosen from the shipped Type II enzymes; name one or two and that pair is planned, or refused
in one sentence. One enzyme at both ends is supported: the backbone can then close on itself, so
the plan dephosphorylates it and screens for an insert that went in the other way round.

It plans both digests and the gel that separates their fragments, works the ligation out in
picomoles, and designs the colony PCR, the diagnostic digest and the sequencing that confirm the
clone. Four files land in the directory you name:

- `product.dna` — the finished plasmid, features carried over and each junction annotated
- `primers.tsv` — every oligo it designed, with length and Tm
- `protocol.json` — the protocol as data: a draft you may edit through `build-protocol`
- `protocol.html` — the page rendered from `protocol.json`, self-contained: reagents, reaction
  tables, thermocycler programs, expected bands, a simulated gel and troubleshooting

The command prints a summary line and the four paths. The same inputs write the same bytes.

## Which route the insert takes

The command line does not show this, and the record you hand in is what decides it.

| What that record reads | What the plan does |
| --- | --- |
| a site of each enzyme | cuts the insert out of it and gel purifies it — no primer is designed |
| neither enzyme's site | amplifies it first, with a spacer and the recognition site on each primer's 5' tail, then cuts the amplicon. A circular template is taken away with DpnI |
| any other count | refuses, naming every site it found |

So an insert already sitting between usable sites in another plasmid costs no primers, and a
bare insert costs one pair.

## The junction is not scarless

The two ends came from a recognition site, so ligating them puts that site back and the product
gains its bases at each end. The plan says what each junction now spells, and where a junction
falls inside a coding sequence it reports whether the added bases hold the reading frame. Tell
the user that before they order, not after they sequence.

## Ask only for what the files do not settle

| Option | What it sets |
| --- | --- |
| `--enzyme` | An enzyme both digests use; give it once or twice. Chosen for you when you name none |
| `--polymerase` | For the insert's PCR, where one runs |
| `--host` | The strain the protocol names |
| `--name` | What to call the product |

## What carries no verdict

A check nothing sourced measures says so instead of passing. Here that is whether two enzymes
can share one tube. Their supplier sells each in a named buffer, so the same buffer for both is
an answer; where they differ, judging the pair needs figures for how much activity each keeps in
the other's buffer, and nothing this package may ship states them. The page shows the check as
`not judged` and names the supplier table to look the pair up in. Read it as unmeasured, never
as fine.

The page's title names the method. The materials and every reaction table name each enzyme by
the product its supplier sells and its catalogue number, so the user can order what the plan
assumed.

## From Python

When you need the numbers rather than the files:

```python
from liulab_mbio.cloning.restriction import plan_restriction

plan = plan_restriction("vector.dna", "insert.dna")
plan.status  # "pass", "warn" or "fail" over every check and every oligo
plan.write("plan/")  # the same four files
```

`plan.enzymes` is the pair, `plan.refusals` every pair weighed against it with the rule that
refused each, `plan.junctions` what each join spells, `plan.amounts` what goes into the ligation,
`plan.diagnostic` the bands a miniprep should give, `plan.colony.clones` the bands each candidate
colony should give, and `plan.dephosphorylates` whether the backbone needs a phosphatase.

## When it refuses

It raises rather than guessing, and the message names the cause: an enzyme reading no site in
the vector, a record neither route can get two cut ends onto, a digest that drops more of the
vector than it keeps, or ends that do not anneal. Where the offending site sits in a coding
sequence the message also names the synonymous codon change that would take it away, so the
user can decide whether to keep that enzyme.

## Before you hand it over

Open `protocol.html` and read it back, as `build-protocol` asks. Then tell the user what the
design chose and what it cost: the pair, the product length, what each junction spells, what
the frame check said, and the bands the colony PCR and the diagnostic digest should give.

Say which checks warned rather than hiding them, and say plainly which carried no verdict.

## The primer side

`primer-design` owns primer questions: which placement a task calls for, how a preference such
as a maximum length goes in, and what a warning left on the order sheet means. The plan designs
its own primers; go there for primers outside one.

## The protocol side

`build-protocol` owns the page and the rules for changing it. When the user asks for what the
plan did not anticipate — skip a step, add a gel purification, add a caution — edit
`protocol.json` as it says and render it again.

A change to what the plan computes — the enzymes, the polymerase, the insert or where it comes
from — goes back through `cloning restriction plan`, which writes `protocol.json` afresh. Within
the session, reapply the changes the user asked for to the new protocol.

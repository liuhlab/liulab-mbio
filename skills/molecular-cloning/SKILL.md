---
name: molecular-cloning
description: >-
  Choose the cloning method a job calls for and plan it with `liulab_mbio`, weighing how many
  fragments are joined, whether the junction may gain bases, what the parts already carry, and
  what the method costs in time and money. Golden Gate and Gibson assembly are planned end to
  end from a vector and its insert sequence files: the enzyme or the kit, every junction's
  overhang or overlap, the PCR primers, the simulated product, the colony PCR that validates
  every junction, and an interactive HTML bench protocol. Use when someone wants to clone,
  insert, subclone or join sequences into a plasmid, asks which cloning method fits their
  fragments, asks for primers with BsaI, BsmBI, BbsI, PaqCI or SapI tails, wants a ready-to-run
  cloning protocol from sequence files, or names a method such as Gibson, Gateway, TOPO,
  ligation-independent or restriction and ligation.
---

# Molecular cloning

`liulab_mbio` plans the cloning. This skill chooses the method; the method's own file carries
the detail, and you read it only after the choice is made.

## What this package plans

```bash
pixi run liulab_mbio cloning --help
```

That list is the supported set. Read it off the command rather than off prose, which ages. Each
method planned here has its own file, read only once you have chosen it:

| Method | Read next |
| --- | --- |
| Golden Gate | `golden-gate/METHOD.md` |
| Gibson assembly | `gibson/METHOD.md` |

## Choosing

Five things decide a cloning method. Ask the user only the ones their files do not answer.

| Method | How many fragments | Whether the junction may gain bases | What the parts carry | Speed | Cost |
| --- | --- | --- | --- | --- | --- |
| **Golden Gate** | one insert or several, in one reaction in the order you give them | the junction spells its own overhang, so the product reads as the user wrote it; `--in-frame` holds every junction on a codon boundary for a fusion | one Type IIS enzyme has to read no site in any part. Where none is free, the plan reports which sites a synonymous codon change could reach | one pot, one day: amplify the parts, cut and ligate together, transform | one tailed primer pair per part, one enzyme and a ligase — no kit |
| **Gibson assembly** | one to five inserts, in one reaction in the order you give them; above the kit's documented count it warns | the junction spells only what the two parts spell, so nothing is added and a fusion stays in whatever frame the parts are in | nothing is cut, so a part carrying a site for every Type IIS enzyme still goes in as it is | one pot, one day: amplify the parts, one 15 to 60 minute incubation, transform | one kit and one tailed primer pair per part; a short part can be stitched from oligos instead |

Where the answers rule out every row, stop rather than reaching for another mechanism.

## Say so and stop, rather than improvise

TOPO cloning, ligation-independent cloning, yeast-mediated assembly, Gateway cloning and
classical restriction and ligation are all real methods this package does not plan.
Nothing here designs their junctions, sizes their reactions or checks their products, so a
protocol written for one would be prose that nothing measured — which is the one thing a plan
from this package is not.

Tell the user which method they asked for, that it is not supported here, and what is. Where a
method in the table above would do the same job, offer it and say what changes for them.

## Every combination of two or more lists is a library

`protein-assembly` owns that job and this skill does not: it builds a barcoded combinatorial
library by iterative Golden Gate, sizes each round's colonies for the coverage asked for, and
writes one protocol covering every round. Send someone who wants every combination of lists of
domains there.

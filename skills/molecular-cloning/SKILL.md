---
name: molecular-cloning
description: >-
  Choose the cloning method a job calls for and plan it with `liulab_mbio`, weighing how many
  fragments are joined, whether the junction may gain bases, what the parts already carry, and
  what the method costs in time and money. Golden Gate is planned end to end from a vector and
  its insert sequence files: the Type IIS enzyme, the whole overhang set, the PCR primers, the
  simulated product, the colony PCR that validates every junction, and an interactive HTML bench
  protocol. Use when someone wants to clone, insert, subclone or join sequences into a plasmid,
  asks which cloning method fits their fragments, asks for primers with BsaI, BsmBI, BbsI, PaqCI
  or SapI tails, wants a ready-to-run cloning protocol from sequence files, or names a method
  such as Gibson, Gateway, TOPO, ligation-independent or restriction and ligation.
---

# Molecular cloning

`liulab_mbio` plans the cloning. This skill chooses the method; the method's own file carries
the detail, and you read it only after the choice is made.

## What this package plans

```bash
pixi run liulab_mbio cloning --help
```

That list is the supported set. Read it off the command rather than off prose, which ages.
Today it holds one method:

| Method | Read next |
| --- | --- |
| Golden Gate | `golden-gate/METHOD.md` |

## Choosing

Five things decide a cloning method. Ask the user only the ones their files do not answer.

| What decides it | Golden Gate fits when |
| --- | --- |
| **How many fragments** join | one insert or several, joined in one reaction in the order you give them |
| **Whether the junction may gain bases** | the junction can spell its own overhang, so the product reads as the user wrote it. `--in-frame` holds every junction on a codon boundary for a fusion |
| **What the parts carry** | one Type IIS enzyme reads no site in any part. When none is free, the plan reports which sites a synonymous codon change could reach |
| **Speed** | one pot, one day: amplify the parts, cut and ligate together, transform |
| **Cost** | one tailed primer pair per part, one enzyme and a ligase — no kit and no synthesised homology arms |

Where the answers rule Golden Gate out, stop rather than reaching for another mechanism.

## Say so and stop, rather than improvise

TOPO cloning, ligation-independent cloning, yeast-mediated assembly, Gibson assembly, Gateway
cloning and classical restriction and ligation are all real methods this package does not plan.
Nothing here designs their junctions, sizes their reactions or checks their products, so a
protocol written for one would be prose that nothing measured — which is the one thing a plan
from this package is not.

Tell the user which method they asked for, that it is not supported here, and what is. Where
Golden Gate would do the same job, offer it and say what changes for them.

## Every combination of two or more lists is a library

`protein-assembly` owns that job and this skill does not: it builds a barcoded combinatorial
library by iterative Golden Gate, sizes each round's colonies for the coverage asked for, and
writes one protocol covering every round. Send someone who wants every combination of lists of
domains there.

---
name: barcode-design
description: >-
  Design a set of DNA barcodes with `liulab_mbio`: every pair of one part list far enough apart
  in mismatches that no two read as one, no restriction or Type IIS site in either orientation,
  the junction with the cloning scar included, and no stop codon where the construct translates
  the barcode. The set is reproducible from a seed, and a set someone already holds is checked
  by the same rules. Use whenever someone wants barcodes designed or checked: naming the members
  of a pooled or combinatorial library, indexing a screen, a barcoded Golden Gate part, asking
  how far apart two barcodes must stand, or whether a GC or homopolymer limit is worth imposing.
---

# Barcode design

`liulab_mbio.barcodes` designs the set and judges it. Your job is to state the length, the
cloning scar, the reading frame and the enzymes. Never write a barcode from your own knowledge:
nothing checks a base you invented, and the package draws each one
(`docs/adr/0002-editable-protocols.md`).

## State the rules, then design

```python
from liulab_mbio.barcodes import BarcodeRules, check_barcodes, design_barcodes

rules = BarcodeRules(11, scar="AGCG", phase=0, forbidden=["BsaI", "BbsI"])
design_barcodes(24, rules)
```

`BarcodeRules` is checked as it is built, so an impossible request fails where it is written
rather than at the bench. Read the docstrings rather than reconstructing a call:

```bash
pixi run python -c "from liulab_mbio.barcodes import design_barcodes; help(design_barcodes)"
```

## The four things to ask the user

- **Length**, and the **cloning scar** the ligation leaves between one barcode and the next.
  Pass `scar=""` where a barcode stands alone.
- **The frame.** `phase` is how many bases of the first codon were read before the barcode, 0
  where it starts one. Pass `phase=None` where the construct never translates the barcode — a
  lineage or screen barcode outside any coding sequence — which takes away both frame rules.
- **The enzymes** the cloning uses, every one in a single `forbidden` list.
- **How many** each part list needs.

Barcode length plus scar length must be a whole number of codons whenever the barcode is
translated. That is what makes the block frame-neutral, and it is refused rather than rounded.

## Distance is held within a part list

Two barcodes need telling apart only where both can fill the same position, so the rule runs
within one part list and not across the whole library. Design each part list with its own call.
The published set this default comes from holds three mismatches within a part list and only
two across, so a library-wide rule would be tighter than the work that succeeded.

## The dials, and where they sit

| Dial | Default | Why |
| --- | --- | --- |
| `distance` | 3 | measured on the published set |
| `metric` | `sequence-levenshtein` | a Hamming rule cannot see a lost base |
| `max_homopolymer` | 5 | the lost base is likeliest inside a run |
| `gc_band` | none | convention, with one measurement against it |

`metric="hamming"` counts mismatches instead, which is what a published set was designed to and
what a short barcode may have to fall back on: the indel-aware metric costs about one base of
length at a given set size. `deletion_ambiguity(barcodes)` is the other side of that choice — the
share of one-base deletions that read as another barcode of the set — so report it with any set
designed on mismatches.

`docs/research/barcode-design.md` is the evidence for all four, and it is what to quote when a
user asks why there is no GC band: every band in the literature descends from one uncited
sentence, and the paper that imposed a band and then measured it found GC was not what
mattered. The homopolymer cap of 5 is a judgement rather than a measurement — say so, and turn
either dial where the user has a reason. Do not add a rule because it is conventional.

## Check a set someone already holds

```python
check_barcodes(["AATCATGGCCT", "ACGCTTATAAT"], rules)
```

One message per rule broken, empty where the set holds. Use it on any set a user brings, on a
published one, and on anything you were handed rather than designed — it reads each barcode on
its own and each pair on its distance.

## When it refuses

`SpaceExhaustedError`, a `ValueError`, when no further barcode can be found. The message says
how many were found and which rule rejected the most candidates. Pass that cause on: the fix is
the user's, a shorter part list, a longer barcode, or a dial turned off. `ValueError` also comes
back from `BarcodeRules` for rules nothing could be designed under, and `KeyError` for an enzyme
name nothing ships.

## Before you hand it over

Give the user the set, the rules it was designed to, and the seed — the same seed and rules
return the same set, which is what lets them regenerate it. Say which dials were on, the metric
among them. A barcode
set is what decodes their sequencing later, so what it was designed to is part of the result.

`codon-optimize` writes the coding sequence a part carries; come here for the barcode naming it.

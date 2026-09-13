---
name: primer-design
description: >-
  Design and check PCR primers for a fragment or a molecular cloning step with `liulab_mbio`:
  state how far each primer may move as a placement — anchored by a tail, near a target it must
  read across, or free inside a region — and the package searches every binding site that
  placement allows and ranks them on every check, so a primer warns only where nothing nearby
  passes. Use whenever someone wants primers designed, judged or improved: amplifying a fragment
  from a plasmid or genomic DNA, cloning primers carrying a restriction or Type IIS tail, colony
  PCR, junction or sequencing primers, a primer with a maximum length or another preference, or
  asks why a primer warns and whether a better one lies nearby.
---

# Primer design

`liulab_mbio.primers` designs the primer and judges it. Your job is to say where each primer may
lie — its **placement** — and to read back what the package found. Never write a primer, a Tm or
a band size from your own knowledge: nothing checks a number you invented, and the package
computes each one (`docs/adr/0002-editable-protocols.md`).

## Which placement the task calls for

| Case | What holds the primer | What moves | Examples |
| --- | --- | --- | --- |
| Anchored | a tail joins its 5' end | the 3' end, by length | Golden Gate, Gibson and restriction-cloning part primers |
| Near a target | a distance range from something it reads or amplifies across | both ends, inside that range | colony PCR, a junction primer, sequencing primers |
| In a region | the flanks of a region to amplify | both ends, and the pair | amplifying a fragment from a plasmid or other template |

One form covers all three. `Placement(five_prime=..., three_prime=...)` takes a `Segment` of the
positions that end may take — 0-based and half-open, crossing the origin of a circular template
where it must. Bound one end, or both; a placement bounding neither is refused.

## Design it

```python
from liulab_mbio.primers import Placement, design_pair, design_primer
from liulab_mbio.sequence import Segment, Strand

design_primer(template, position, Strand.FORWARD, placement=..., tail=..., name=...)
design_pair(template, start, end, forward_placement=..., reverse_placement=...)
```

`position`, and a pair's `start` and `end`, stay the 5' end you asked for: it breaks a tie
between equal candidates, so pass where you would have put the primer yourself. A pair's two
placements bound the amplicon between them, which is why nothing takes an amplicon size.

Each case is worked through in the docstrings, on a template you can paste. Read them rather
than reconstructing a span:

```bash
pixi run python -c "from liulab_mbio.primers import design_primer; help(design_primer)"
```

Judge what comes back, or anything the user hands you, with `evaluate_primer(primer, template)`
and `evaluate_pair(forward, reverse, template)`. `report.status` is the worst of its checks, and
`report["gc_percent"]` reads one.

## What a warning means

The design judges every binding site the placement allows, at every length, on every check —
structures, runs, binding and off-target sites included — and takes the best. So a warning on a
designed primer means **no option within its placement passed**.

Tell the user that, rather than hiding it or apologising for it: name the check and its value,
say what holds the primer there (the tail, the junction, the region they named), and that only a
wider placement or a changed constraint could clear it. A primer that warns wherever it may go
is worth ordering; one nobody can explain is not.

## A preference goes in through the thresholds

A user asking for primers no longer than 25 bases is setting a band, not choosing a primer:

```python
import dataclasses

from liulab_mbio.primers import THRESHOLDS, Band

thresholds = dataclasses.replace(THRESHOLDS, length=Band(18, 25, 15, 25))
```

Pass it as `thresholds=`, and the package still does the choosing. `THRESHOLDS_FOR` holds the
defaults per role — `"amplification"`, `"colony PCR"`, `"sequencing"` — and `Band(low, high,
warn_low, warn_high)` is what passes and what only warns. Never meet a preference by picking a
primer by hand.

## What it does not check

Binding sites and off-target sites are found on the template you pass, and nowhere else. A
primer meant for a genome needs a specificity check this package does not do — Primer-BLAST or
similar. Say so plainly; never let a clean report read as genome-wide specificity.

## When it refuses

`ValueError`, naming the cause: no annealing region fits the placement, or fits the template at
that position. Widen the span, or move the position, rather than catching it — a placement
holding nothing is a question about the task, not a failure to handle.

## A whole cloning experiment

A Golden Gate plan designs, checks and orders its own primers; `golden-gate-assembly` runs it.
Come here for primers outside a plan, or to explain a warning one left.

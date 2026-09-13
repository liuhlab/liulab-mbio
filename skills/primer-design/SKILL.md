---
name: primer-design
description: >-
  Design and check PCR primers for a fragment or a molecular cloning step with `liulab_mbio`:
  state how far each primer may move as a placement — anchored by a tail, near a target it must
  read across, or free inside a region — and the package searches every binding site that
  placement allows and ranks them on every check, so a primer warns only where nothing nearby
  passes. Primers for genomic DNA are then checked against the genome itself, offline, and
  redesigned around any off-target amplicon. Use whenever someone wants primers designed, judged
  or improved: amplifying a fragment from a plasmid or genomic DNA, cloning primers carrying a
  restriction or Type IIS tail, colony PCR, junction or sequencing primers, a primer with a
  maximum length or another preference, asks whether a pair is specific on a genome or would
  amplify somewhere else, or asks why a primer warns and whether a better one lies nearby.
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

**Where the task leaves room, spend it before you report a warning.** A colony PCR distance
range, a region's flanks, a sequencing window — widen the one the task itself left open, design
again, and keep the better result. Never change a constraint the user set, and never move a
vector cut to clear a warning: a cut moves for the overhang rules only. Then explain what is
left.

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

## A genomic primer needs a genome check

Binding sites and off-target sites are found on the template you pass and nowhere else, so a
clean report covers that template only. Check the genome as well whenever the primers amplify
genomic DNA, or must not amplify a host — colony PCR on a genomic background. A plasmid template
needs no genome check: the plasmid is the template.

The genome comes from liulab-genome. Never use `Genome(name)`, which registers and can download:

```bash
genome assembly files ce11 --json
```

It downloads nothing and prints one JSON line: `genome_files.fasta` and `genome_files.fai` are
absolute paths, and the top-level `assembly` is the name to report. The top-level `files` maps
file names to byte sizes, not paths. Exit 1 with empty stdout means the assembly is not
registered or not trusted — pass its stderr on, which names `genome assembly register`. Show the
user a `--force` it suggests rather than running it, because that can download again. Exit 2 with
`No such command 'files'` means the installed liulab-genome predates the command, which is on
main and in no release yet: say so rather than guessing a cache layout.

The search runs `ipcr`, which `pixi install` brings in — a pixi dependency, not on PyPI.

```python
from liulab_mbio.primers import Locus, design_pair_on_genome, evaluate_pair_on_genome

evaluate_pair_on_genome(forward, reverse, fasta, assembly, intended=Locus(name, start, end))
design_pair_on_genome(fasta, assembly, Locus(name, start, end), flank=200)
```

Check a pair you already hold, or design one for a region: the design places the primers in the
region's flanks, checks the best-ranked pairs in one search, and designs around an off-target
amplicon until a pair is specific. `evaluate_on_genome` checks several pairs in one search. Read
the docstrings rather than reconstructing a call.

A `Locus` names its sequence as the FASTA spells it. A combined genome keeps each component's
own chromosome names and suffixes them, so ce11 in `ce11_ecHT115` reads `I__ce11`, not
`chrI__ce11`; `details.separator` carries the separator.

## What the genome answers

- `report.off_target` is every amplicon but the intended one, each with where it lies, how long
  it is, and `made_by` — `"pair"`, or `"forward"` or `"reverse"` where one primer makes it alone.
  Give the user the size and the place: that is what a gel would show.
- `design.specific` says the pair makes only its intended amplicon, and `design.rounds` how many
  searches that took. Where it is false the rounds ran out, and the best pair found comes back
  with its off-target amplicons — a result to report, not to bury.
- `report.near_matches_checked` false means the genome was too large to search for near matches,
  so only perfect ones were found. Say so every time it is false, or a clean mouse or human
  result reads as a specificity nobody checked.

When nothing specific turns up, wider flanks give the design more room; the other option is to
accept the pair with its off-target amplicons stated. Which of the two is the user's call.

## When it refuses

`ValueError`, naming the cause: no annealing region fits the placement, or fits the template at
that position. Widen the span, or move the position, rather than catching it — a placement
holding nothing is a question about the task, not a failure to handle. A genome call also
raises on a FASTA with no index, a missing `ipcr`, or a search past its timeout, each message
naming its fix.

## A whole cloning experiment

A Golden Gate plan designs, checks and orders its own primers; `golden-gate-assembly` runs it.
Come here for primers outside a plan, or to explain a warning one left.

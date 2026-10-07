---
search:
  exclude: true
---

# Prototype: splitting a cargo for Golden Gate

Throwaway code for issue #261. It is not part of the package, nothing in `src/` imports it, and
the gate does not collect it. It exists to produce the measurements the human decides from.

```bash
pixi run python prototypes/dad-split/measure.py
```

That one command prints every table. An argument prints a subset — `old`, `m1`, `m2` or `m3`,
in any combination — which is how the long tables were run side by side.

`split.py` holds three searches over the same inputs. `measure.py` runs them over
`docs/examples/ap1-library/` and prints the tables.

## The three searches

| Name | Shape |
| --- | --- |
| `simple_split` | Cuts spaced evenly, then one greedy first-fit pass over a window of 16 bases, ranking candidates by their own on-target count. This is what `cloning.goldengate.design.design_overhangs` does today, lifted to positions it picks itself. |
| `greedy_set_split` | The same pass — same cuts, same window, first-fit, no backtracking — ranking each candidate by the fidelity of the set it would make instead of by its own on-target count. One knob changes: the key. |
| `clever_split` | An interval dynamic program over fragment length pins the fewest fragments and the exact window each cut may take, then depth-first branch and bound over the per-cut candidates maximises Pryor 2020 fidelity. |

The bound is admissible. Adding an overhang to a set can only raise the denominators of the
overhangs already in it, and every further factor is at most one, so a partial set's fidelity
is an upper bound on every set extending it.

Three things the search holds that the dynamic program alone does not. A window is a marginal:
it admits no position that no whole partition contains, but one position per window picked
independently can still break the chain. So the walk runs left to right and holds each step to
the length bound and to the backward reachability mask. And `cap` keeps the candidates nearest
an evenly spaced anchor, not a spread across the window — spread picks sit too far apart to
satisfy the chain, and the search then reports an infeasibility that belongs to the thinning.

`liulab_mbio` is read and never written. Scoring is `overhangs.fidelity`, refusal is
`overhangs.refusal` with `reserved=`, and the ligase matrix is read by `ligase.read_profile`.

## What this assumes, and does not decide

**Where the Type IIS sites sit (§6.1) is not decided here.** The prototype only needs one
number from that decision — how much of an oligo is not cargo — so it takes it as
`Budget.overhead`, set to 54 nt: two 20 nt primer sites, two 6 nt enzyme sites and a 1 nt
spacer each. Both options in §6.1 give the same arithmetic at this granularity, because the
site costs its six bases wherever it is carried. Only the sign of the trade differs: sites on
the gene are cargo the split must route around, sites on the primers are not. The prototype
assumes the second and treats the whole cargo as cuttable. Change `overhead` and every table
moves; nothing else in the code depends on the choice.

**The padding rule (§6.2) is not decided here, and nothing pads.** Fragments come out at
whatever length the dynamic program gives. Table 3 reports the widest spread between the
shortest and longest oligo of one part, which is the number the vendor's 15% floor is judged
against, so the decision has evidence rather than an assumption.

Two smaller assumptions, both parameters: the shortest fragment is 40 bases, as OMEGA uses, and
`min_distance` stays at the package's `MIN_DISTANCE` of 2 unless a table says otherwise.

## The tables

Tables 1 to 8 measure the first two searches. Three more answer what the first round could not.

| Table | Question |
| --- | --- |
| M1a, M1b | All three searches over the 72 parts, at a 300 nt and at a 350 nt oligo. |
| M1c | All three against fragment count on the 2,276 bp cargo. |
| M2 | Whether a wall-clock budget of 10, 60 or 300 seconds turns the branch and bound's `stopped early` into a proof, and what the nodes do meanwhile. |
| M3 | How the dynamic program, the greedies and the branch and bound scale to 20 kb of cargo, at a 300 nt and at a 1,000 nt oligo. |
| M4a, M4b | All four searches over the 72 parts, the fourth an uncapped branch and bound held to 15 s per part, so the proved and stopped columns mean optimality and not optimality over kept candidates. |
| M5 | Whether widening the set-fidelity greedy's window removes the refusals M1c found at 13 and 14 fragments. |

M3's cargo is AP-1 parts concatenated in frame, each trimmed to a codon boundary and the list
cycled. Random DNA would change the base composition, and the answer with it.

Every runtime is indicative: one machine, one run, three searches sharing it.

## Inputs

- Cargo: `docs/examples/ap1-library/parts.tsv` (72 blocks, 133 to 1,149 bp) and `product.dna`
  (2,276 bp).
- Ligase matrix: `reference_docs/ligation-fidelity/potapov2018/FileS03_T4_18h_25C.xlsx`, which
  the user holds and the repository does not ship. Without it the shipped BsaI matrix scores
  instead, and the run says which it used.
- Held out by name: `AGGA` and `TTCC`.

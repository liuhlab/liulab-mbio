---
search:
  exclude: true
---

# Splitting a gene for Golden Gate: what DAD is, what implements it, and what is left to build

Research note for issue #235 (parent #225). Everything below was read on **2026-10-06** unless a
line says otherwise. It establishes what the two existing implementations of the gene split
actually do, judges the one that ships code, names the search problem, and draws the line between
what `mbio` already covers and what #236 has to build.

Step 3 of cargo synthesis is the split: break each gene into oligo-sized fragments whose internal
overhangs are chosen on ligation fidelity, with `AGGA` and `TTCC` held out of the set. Nothing in
this repository does it. `docs/research/synthesis-and-assembly.md` records the hole in its Open
gaps: *"The split tool is unresolved. The method it names is not the tool its own source ran, and
neither is pointed at."* This note closes the first half of that sentence.

## 1. The method Lund cites, and the tool Lund ran

They are not the same thing, and Lund et al. never say they are.

Lund et al. 2024 (`reference_docs/synthesis_and_assembly/long_fragment_GGA/Lund2024/main.txt`,
read 2026-10-06) defines the method in its introduction:

> These overhangs can be methodically selected using ligation fidelity data in a process referred
> to as Data-optimized Assembly Design (DAD), which has been shown to support assemblies in as
> high complexity as 52 fragments. Rather than being a deterministic model in which a particular
> set of Golden Gate overhangs performs optimally, DAD selects overhangs based on compatibility
> with one another and ranks potential sets based on predicted fidelity.

and then, in the next paragraph, names the tool:

> Here, we apply DAD for gene construction from oligo pools using the fragment design NEBridge
> SplitSet Tool.

The three citations attached to the DAD sentence are refs 14–16 of that paper: Pryor et al. 2022
*ACS Synth. Biol.* (40 kb from 52 parts), Pryor et al. 2020 *PLoS One*, and Sikkema et al. 2023
*Curr. Protoc.* All three are NEB papers about overhang sets. **None of them is SplitSet**, and the
SplitSet sentence carries no citation at all.

NEB itself draws the same line, and more plainly than Lund does. From NEB's own request form for
the SplitSet Lite API (`neb.com/en-us/forms/nebridge-splitset-lite-api`, read 2026-10-06):

> The API is intended for use in identifying an optimal set of overhangs for NEBridge Golden Gate
> Assembly near specified breakpoints in a target sequence **using Data-Optimized Assembly Design
> as described in our publication** "Enabling one-pot Golden Gate assemblies of unprecedented
> complexity using data-optimized assembly design."

So the reading, now from both sides:

| | What it is |
| --- | --- |
| **DAD** | A named design principle: score a candidate overhang set against measured ligation-fidelity data, and rank candidate sets by that score. |
| **SplitSet** | One of NEB's tools that *applies* DAD. A tool, not a method; a web form, not a citation. |

The two are related as method and implementation, which answers the ticket's first question: **not
the same thing, and NEB says so.** The paper NEB anchors DAD to is *"Enabling one-pot Golden Gate
assemblies…"* — **Pryor et al. 2020, *PLoS One***, not the 2022 paper that carries DAD in its
title. Pryor 2022 agrees, attributing the naming to the earlier work: *"Recent work by our
laboratory has significantly expanded the capacity of GGA by using ligase fidelity data to select
fusion sites, a process termed data-optimized assembly design (DAD)."*

That matters directly. Pryor 2020 is the paper this repository already ships data from, as
`src/mbio/data/ligation_fidelity.json` under CC BY 4.0
(`docs/research/ligation-fidelity.md`). **The scoring function DAD is named for is already in the
package**, computed by `mbio.overhangs.fidelity`. Section 4 checks that against the one
independent implementation available.

### 1.1 Nobody publishes the search

Read Lund's definition closely: DAD *"selects overhangs based on compatibility with one another and
**ranks** potential sets based on predicted fidelity."* Ranking is scoring. **DAD names the
objective and says nothing about how the candidate sets are generated.** That is the search, and
the search is published nowhere:

- **NEB does not document it.** Across the SplitSet help page, the tool form, the Ligase Fidelity
  landing page, the NEBridge tool overview and the two code-request forms (all read 2026-10-06),
  no page says how split points or overhangs are chosen. Every statement is outcome-level — fusion
  sites are *"optimized"*, overhangs are *"selected to maximize the expected ligation fidelity"*.
- **Pryor 2022 does not document it.** The PMC full text (PMC9208013, CC BY-NC-ND, read
  2026-10-06) contains no occurrence of "algorithm", "exhaustive", "greedy", "heuristic" or
  "anneal". The closest it comes is constraining the search space, not describing the procedure:
  *"The SplitSet tool was used to select 10 fusion sites, equally spaced in the genome, predicted
  to be an overall high-fidelity set."* It points at the tools instead: *"the optimizer code used
  in the SplitSet tool is available under a non-commercial use license on request."*
- **The code exists and we cannot use it.** NEB's Overhang Optimizer, holding the `GetSet` and
  `SplitSet` modules, is released by request under terms that state: *"The tools are for internal
  research only. No commercial use, reproduction or distribution of the tools is permitted. You
  shall not cause or allow the reverse engineering, disassembly, or decompilation of a tool."*
  Reading it would contaminate anything we ship.

**There is no published algorithm to reimplement or to benchmark against.** Our splitter has to
stand on its own stated objective, which is the one DAD names and the one we already compute.

### 1.2 What SplitSet does, which is worth knowing even though we cannot call it

From NEB's tool form and help (`ligasefidelity.neb.com/splitset/`, read 2026-10-06), SplitSet takes
a sequence, a topology, an overhang length of 3 or 4, a ligation-condition dataset (15 of them,
with the enzyme folded into the dataset rather than named separately), and then three inputs that
map one-to-one onto the problem statement in section 3:

| SplitSet input | What it is in our terms |
| --- | --- |
| *"Define split regions based on: Number of fragments / Maximal fragment size / Minimal fragment size"* | the fragment-length bound, and the fragment count |
| *"Split regions — Provide non-overlapping regions in the nucleotide sequence where it can be divided"* | per-junction position windows |
| *"Excluded overhangs — Specify overhangs that must not be present in the resulting split set"* | **exactly the `AGGA`/`TTCC` hold-out our method needs** |

It returns the overhang set, the overhang coordinates, the fragments and a fidelity defined as
*"the total fraction of correct ligation events out of all ligation events"* — the same quantity as
ours. NEB separates its three tools cleanly: the **Ligase Fidelity Viewer** scores a set you
already hold, **GetSet** invents a set when the overhangs *"may be arbitrary"*, and **SplitSet**
cuts a given sequence so that its own bases form a good set. Our splitter is a SplitSet, not a
GetSet, and `overhangs.fidelity` is our Viewer.

It is still of no use to an agent. The full tool is web-only and runs as a queued job behind a
request ID; the only programmatic route is the **SplitSet Lite API**, gated by a registration form
and per-user Terms of Service, and SplitSet Lite is capped at about 20 fragments where the full
tool is not. A per-user credential is not something a package can depend on.

## 2. The OMEGA implementation, judged

### 2.1 Which repository, and why it is the right one

`https://github.com/RomeroLab/omega`, commit **`160be2f7292676f64d54dad58902d5de557d9487`**
(authored 2025-11-17, branch `main`), read 2026-10-06. Python, about 1,900 lines across seven
files in `code/`. Licence **GPL-3.0**.

Three independent links tie it to Freschlin et al. 2026:

1. The paper's own data-availability statement (`Freschlin2026/main.txt`): *"OMEGA code is
   available on GitHub and archived with Zenodo (<https://zenodo.org/records/17637683>)."* The
   supplement (`Freschlin2026/supp.txt`) is explicit: *"full implementation is available at
   github.com/RomeroLab/omega."*
2. The repository description quotes the paper's title verbatim, and the README opens *"This is
   the code to run the OMEGA program presented in Freschlin et al."*
3. The README's Zenodo badge carries concept DOI `10.5281/zenodo.17637682`; the paper cites record
   `17637683` under that concept.

No GitHub discussion thread was needed; the paper points at the code directly. The one other
candidate that surfaced, `scbirlab/ogilo`, is a different lab's oligo-pool concatenator with no
fidelity optimisation and no connection to the paper.

**The assembly method in this paper was rejected for our use in favour of the Baker lab's
three-primer scheme, and nothing here reopens that.** What is read below is the cut design only.

The README states its own caveat, which should be read before any of the judgements that follow:
*"We have not exhaustively tested this code."*

### 2.2 What it optimises

`code/predict_fidelity.py:43`. For a set of overhangs, the product over the set of

> (correct ligations for the overhang, plus correct ligations for its Watson–Crick partner)
> ÷ (every ligation either of them was seen making against the set and the set's reverse
> complements)

This is Pryor 2020's definition. It is the same quantity `mbio.overhangs.fidelity`
computes, down to counting both ends of each junction — the detail the docstring of our own
function says is what reproduces the paper's worked examples. **Two independent implementations
of the same formula agree, which is the best confirmation available that our scoring half is
right.**

Two further metrics are computed — `min_gene_fidelity` (the worst single gene in a pool) and
`min_site_fidelity` (the worst single overhang) — but they are **reported, not optimised**. The
README says so in its own words: the optimised metric *"assumes that all GG sites are being used in
a single sequential assembly — it does not fully reflect OMEGA conditions"*, while
`min_gene_fidelity` *"is the more relevant metric for OMEGA."* **OMEGA's search maximises a number
its own README calls the wrong one.** That is not a subtlety to pass over; it is the first thing to
fix in any design that borrows the shape.

### 2.3 How it searches

Simulated annealing over the joint space of cut positions and overhangs. The two are not separate
variables: every candidate overhang is a 4-mer the gene already spells at a candidate cut
(`Gene.__populate_sites`, `library_classes.py:710`), so **choosing where to cut chooses the
overhang**, and the product is scarless by construction. OMEGA never recodes to reach a better
overhang.

One move (`Gene.shuffle_site`, `library_classes.py:754`): pick a gene at random, drop one of its
cuts, split the sequence at the cuts that remain, take the **longest** resulting section, and
re-place the dropped cut anywhere in that section leaving both neighbouring fragments between a
floor and the oligo's coding space. Accept by Metropolis on a temperature ramp from 0.005 to
0.00001 across `nopt_steps` (default 1000), with a hard reset to the best state whenever a
candidate falls more than 0.01 fidelity below it. `nopt_runs` independent seeds run in parallel and
the best result is kept. A `greedy` variant (`Pool`) accepts improvements only.

### 2.4 What it gets right

Four things, and they are worth taking even though the code is not.

1. **The objective is the published one, on published data.** Nothing is invented, and the
   ligation matrices are Potapov 2018 and Pryor 2020 read from CSV, cited by PMID in the directory
   names.
2. **The move repairs feasibility instead of testing it.** Re-placing the dropped cut only inside
   the longest section, bounded by both neighbours, means every state the chain visits is a legal
   fragmentation. There is no penalty term and no infeasible region to wander into. That is the
   right way round for a hard length constraint, and it is the one piece of the search design
   worth copying outright.
3. **It simulates the product before emitting oligos.** `Gene.gene_assembles` reassembles the
   fragments overhang by overhang, translates the result, and raises if the protein differs from
   the input. A real round trip, not an assertion about the design.
4. **The padding is screened across its own boundaries.** `Gene.__add_padding` redraws random
   filler until neither the join with the primer nor the join with the payload spells the enzyme
   site or any caller-named illegal sequence. Checking a sequence that only exists where two
   designed pieces meet is easy to forget, and this does not.

### 2.5 What is unusable

These are not polish. Each is a correctness defect, verified against the commit named above.

1. **The palindrome exclusion list is wrong.** `CCGG` is a 4-base palindrome and is missing; `GGCC`
   is listed twice. The literal holds 16 entries but 15 distinct values, against 16 palindromic
   4-mers; the set difference is exactly `{CCGG}`. The same literal is pasted in two places —
   `junctions.py:162` and `library_classes.py:722` — so both paths admit 241 overhangs where 240
   are legal, and a fragment taking `CCGG` ligates to itself. Our `overhangs.refusal` derives the
   rule (`candidate == reverse_complement(candidate)`) instead of listing it, and cannot have this
   bug. A hand-written list of a thing a predicate computes is the defect; the missing member is
   only its symptom.
2. **The parameter that would hold out `AGGA` and `TTCC` is dead.** `JunctionSet.__init__` accepts
   `excluded_sites` and stores it at `junctions.py:73`; nothing reads it again anywhere in the
   repository. `__site_options` filters on `fixed_sites` and the palindrome list only. `omega.py`
   passes the value through from the command line. **A user who excludes an overhang gets no error
   and no exclusion.** This is precisely the requirement our method states for step 3, and in
   OMEGA's junction-set path it silently does not work. (The library path is better: `Gene.__populate_sites`
   does exclude `other_used_sites` and their reverse complements — so the exclusion we need exists
   in the codebase, under a different name, in the other module.)
3. **`min_size` is dead in the same shape.** `Library.__init__` stores it at
   `library_classes.py:103` and nothing reads it; all four calls to `shuffle_site` pass the literal
   `min_dist=40`. The config key and the command-line flag do nothing. Two dead parameters of the
   same shape in one small codebase is a pattern, not an accident, and it means the config file is
   not a description of what ran.
4. **The move proposal is seeded with a constant.** `Gene.shuffle_site` ends in
   `candidates.sample(n=1, random_state=42)`: the replacement cut is not drawn at random from the
   legal candidates but is a deterministic function of the candidate frame.
   `get_start_sites_range` shuffles with `random_state=42` likewise. The acceptance test is random;
   the proposal is not. A chain that cannot propose two different moves from one state is not
   sampling the neighbourhood the method describes, and the independent seeds share the proposal
   rule. Whether this costs fidelity in practice is unmeasured here — the claim is only that the
   code does not do what the method says.
5. **A bare `except: continue` hides infeasibility.** The whole candidate block in `shuffle_site`
   sits inside `try: ... except: continue` within a 50-attempt loop, and on exhaustion the function
   falls off the end and returns `None`. `permissible_pos[0]` raises `IndexError` when no legal
   position exists — the genuine "this gene cannot be split here" case — and that is
   indistinguishable from a typo. `SAPool.optimize` tests for the `None`; `Pool.optimize`
   (`library_classes.py:356`) does not.
6. **No minimum-distance rule.** OMEGA refuses exact repeats, reverse-complement repeats and
   (most) palindromes, and otherwise trusts the data. `MIN_DISTANCE = 2` in `overhangs.py` is the
   guard it has no equivalent of. An earlier draft of this note argued the guard is needed because
   the matrix is sparse — a pair whose cross-ligation was never observed contributes nothing to
   the denominator, so it would score as perfectly orthogonal. That mechanism is real but it has
   nothing to bite on here: #327 counted the shipped BsmBI matrix at **65,536 of 65,536 cells
   present**, 4,130 of them nonzero, and only 4.4% of the pairs one base apart measure zero. What
   justifies the rule is what the data says, not what it omits — a third to a half of all pairs
   one base apart cross-ligate at or above `MODEST_MISMATCH`, and two bases apart almost none do.

### 2.6 Where it breaks on a 2.3 kb cargo

Arithmetic taken from the code at this commit, for BsaI, 300 nt oligos and the 20 nt Subramanian
primers OMEGA ships.

| Quantity | Value | Where |
| --- | --- | --- |
| Coding space per fragment | 300 − 20 − 20 − (1+6)·2 − 4 = **242 nt** | `get_coding_space` |
| Coding space used to pick the fragment count | the same with `oligo_len − 24` → **218 nt** | `Library.estimate_nfrags` |
| Fragments for 2,300 nt | first `n` with `2300 // n ≤ 218` → **11** | `Library.estimate_nfrags` |
| Overhangs in the set | 10 internal junctions + 2 backbone = **12** | |
| Starting window per junction | step 210, biggest centroid gap 210, `(242−210)//2` → **±16 bases**, 33 positions | `get_start_sites_range` |

Four things follow.

- **The fragment count is chosen, not optimised.** `estimate_nfrags` takes the first value that
  fits and stops. But the fragment count is the number that predicts whether the assembly works:
  Lund's own measured table (held in `long_fragment_GGA/README.md`) runs 100% error-free clones at
  2 fragments, 84.6% at 5, 66.7% at 8, **40.0% at 12** and 0% at 16. OMEGA's splitter puts our
  cargo at 11 — one below the point where Lund's rate has already more than halved — and nothing
  in the search trades a longer oligo or a different enzyme against that count.
- **The feasibility guard is a hand-tuned constant.** If `biggest_diff` exceeded `coding_space`,
  `half_step` would go negative and `.loc[c−half:c+half]` would select nothing, failing inside the
  bare `except`. The `oligo_len − 24` in `estimate_nfrags`, commented *"subtract some bp from
  oligo_len so we do not get sequences that perfectly fit oligo_len"*, is the only thing keeping
  that from happening. It is a fudge factor, not a bound.
- **Every gene in a pool is forced to the same fragment count**, set by the longest. A design list
  of mixed lengths pays the longest gene's oligo count for every member. For a library of
  same-length variants, which is what OMEGA is for, that is deliberate — it holds assembly
  efficiency constant. For our design lists it is waste.
- **The ±16 window is narrow enough that the result depends on the starting point**, and the
  starting point is drawn with `random_state=42`.

### 2.7 What fixing it would cost, and whether to

The six defects in 2.5 are each a few lines: the palindrome list becomes a predicate, the two dead
parameters get wired through, `random_state=42` is dropped, the bare `except` becomes a typed one,
a distance rule is added. Call it a day's work, and the code would then do what it says.

That is not the reason to leave it. Two architectural facts are not patchable, and one licence
settles the question:

- **The objective is pool-wide, and ours is not.** OMEGA assembles a whole subpool in one tube
  (README, assembly protocol step 5), so all ~50 junctions in a pool must be mutually orthogonal.
  Our method retrieves each design from the pool by PCR and assembles it alone, so only one gene's
  12 overhangs need be orthogonal at once. **We have the easier problem and must not inherit the
  harder one.** Most of OMEGA's search cost — the annealing, the parallel seeds, the thousand steps
  — is paid for a constraint our method does not have.
- **The fragment count is fixed before the search rather than searched over**, and it is the number
  that predicts success.
- **GPL-3.0.** The code cannot be vendored or adapted into `liulab-mbio` whatever its quality. Read
  it, cite it, write our own. The ligation matrices it ships under `data/ligation_data/` are the
  same NEB datasets we already hold, and must be taken from the publishers' supplements rather than
  from this clone.

**Verdict.** A faithful, readable implementation of the right objective, wrapped in a search that
is over-powered for our instance and under-tested for its own, with a scattering of defects that a
rule expressed as a predicate would have made impossible. Worth reading once, in full. Not worth
depending on, and not legally available to us anyway.

## 3. The search problem, named

### 3.1 The statement

Choose `k − 1` cut positions in a sequence of length `L`. A cut's overhang is the 4-mer the
sequence already spells there, so **positions and overhangs are one variable, not two** — unless
the design may recode, which OMEGA never does and Lund does only once, before the split.

Constraints: every fragment between a floor and the oligo's coding space; no two overhangs equal or
reverse-complementary; no palindrome; the reserved overhangs held out; and, ours, a minimum Hamming
distance and a primer tail that spells no further enzyme site.

Objective: the Pryor 2020 fidelity of the resulting set — a product over the set of ratios whose
denominators depend on the whole set.

### 3.2 The size

For `L = 2300`, `k = 11`, fragments of 40–242 nt, counted by a dynamic program over positions:

| Space | Size |
| --- | --- |
| Legal cut placements | **1.23 × 10^19** |
| Subsets of 12 from the 240 legal overhangs | 5.77 × 10^19 |
| Placements pinned to OMEGA's ±16 windows (33^11) | 5.05 × 10^16 |

### 3.3 Why it is hard

- **The objective does not decompose along the sequence.** Each overhang's score divides by its
  cross-talk against every other overhang in the set. A dynamic program over cut positions would
  have to carry the set chosen so far in its state, and that state is a subset of 240. The
  interval-partition structure that makes the length constraint trivially tractable is exactly the
  structure the objective does not respect. This is the whole difficulty in one sentence.
- **The set-selection half is a clique problem.** Threshold the matrix to a 0/1 "these two do not
  cross-ligate" relation, and "is there a set of `k` mutually compatible overhangs" is a `k`-clique
  in the compatibility graph. Maximum clique is NP-hard, so no polynomial exact method should be
  expected for the general objective. *(This is an argument from the structure, not a theorem
  anyone has published about this instance.)*

### 3.4 The honest naming

**The general problem is a heuristic's problem. Our instance is a branch-and-bound's.** `k` is
10–12, not 52; the pairwise rules do almost all the pruning; and the per-junction candidate list,
once a position window is pinned, is a few dozen values.

The decomposition that follows: **pin the cut positions first**, by an interval dynamic program on
length alone — exact, cheap, and it also answers "how few fragments does this gene need?" — then
search the overhangs within a window at each pinned position. That second step is a bounded
constraint problem with about 11 variables, a few dozen values each, and only pairwise constraints.
It is not globally optimal; it is optimal given the positions. Widen the windows and re-run when no
set is found, rather than reaching for a stochastic search.

Greedy is not enough — it has no way to back out of an early overhang that blocks a later junction,
which is exactly what happens when the windows are narrow. Simulated annealing is more than this
instance needs, and brings a seed, a schedule and a convergence question we would then have to
justify. Branch and bound over per-junction candidate lists is the right class, and the repository
already runs that algorithm once (section 4).

## 4. What `mbio` already has

The ticket's claim is *"The scoring half exists; the search half does not."* Checked against the
code, the first half is right and the second is too strong.

### 4.1 The scoring half: confirmed

| Function | What it does |
| --- | --- |
| `overhangs.fidelity` | Pryor 2020's product, over the shipped matrices for BsaI, BsmBI, Esp3I, BbsI and SapI. Verified above to compute the same quantity as OMEGA's `predict_fidelity`. Also reports the weak junctions and every mismatch above `MODEST_MISMATCH`. |
| `overhangs.refusal` | Six rules, each naming itself: `length`, `palindrome`, `uniform`, `repeat`, `near-duplicate`, `site`. The palindrome rule is a predicate, not a list. |
| `overhangs.scoring` | Picks the enzyme's own matrix, a user-held ligase profile, or neither. |
| `ligase.read_profile` | Reads a ligase matrix from a copy the user holds, for an enzyme nobody has published. |
| `overhangs._by_rule` | A rule-based fallback score where there is no data at all, honestly labelled as a ranking and not a prediction. |

Everything DAD names as its objective is here, already cited, already licensed.

### 4.2 The search half: two searches exist, and neither is the missing one

| Where | Shape | Decides | Optimises |
| --- | --- | --- | --- |
| `cloning/goldengate/design.design_overhangs` | greedy, first-fit, one pass, **no backtracking**; junctions ordered least-free first | one overhang per junction; the junction may move within its own `window` | **nothing.** Free candidates are ranked by `table.count(one, reverse_complement(one))` — the overhang's own on-target count, a per-overhang proxy — and the first that `refusal` accepts is taken. `fidelity` is called once at the end to *report*, not to choose. |
| `synbio.igga.standard._settle` | **exhaustive depth-first with a lower bound** — branch and bound — over per-junction candidate lists | one overhang per junction; positions fixed | minimum total amino-acid cost. Fidelity is not in the objective. |

`_settle`'s own docstring states the structural fact the split problem shares: *"every rule is
either about one overhang or about a pair of them, so a set is allowed exactly when each of its
pairs is."* That is the property section 3.4 leans on, and it is already written down and already
exploited — in synbio, for a different objective.

### 4.3 Where the line actually falls

| Half | State |
| --- | --- |
| Score a set by Pryor 2020 fidelity | **present** |
| Refuse one candidate by rule | **present** |
| Hold a named overhang out of the set | **present in effect, absent as an argument.** `refusal` has no `exclude` parameter and `design_overhangs` has no way to say "not `AGGA`". synbio reaches it by passing the reserved overhang in `taken`, which makes it a `repeat`. It works; it is not named, and a caller has to know the trick. |
| Search overhangs at positions somebody else chose | **present twice**, in two different shapes, in two different packages, **neither on fidelity** |
| Search cut positions under a fragment-length bound | **absent.** Nothing in the repository chooses where to cut. `Junction.position` is given by the caller; `window` moves it by a handful of bases around that position, and no code relates one junction's window to another's or to any length budget. |
| Choose how many fragments a gene needs | **absent** |
| Use fidelity as an objective rather than a report | **absent** |

**The line, in one sentence: `mbio` can score any overhang set and refuse any candidate, and
can fill in overhangs at positions a caller already chose — it cannot choose the positions, and it
has never once used fidelity as an objective rather than as a report.**

That is two gaps, not one, and the second is the one the ticket does not predict. "The scoring half
exists" is true of the *function* and not of its *use*: the only existing search over overhangs
ranks candidates by a per-overhang proxy and never consults the set score it then prints.

## 5. What #236 should take, leave and build

### Take

- Pryor 2020 fidelity as the objective. Already in the package, already cited, already licensed.
- OMEGA's feasibility-repairing move, if a stochastic search is ever needed.
- OMEGA's round trip: reassemble the fragments and translate before emitting oligos.
- OMEGA's padding screen, which checks the sequence that exists only where two designed pieces
  meet.
- `_settle`'s branch-and-bound shape, which is already in this repository and already handles the
  pairwise structure.

### Leave

- The pool-wide objective. We assemble one design at a time; adopting OMEGA's constraint would make
  our problem harder for nothing.
- Simulated annealing. Over-powered at `k ≈ 11`, and it brings a seed and a schedule to justify.
- Every line of OMEGA's code: GPL-3.0.
- Hand-written exclusion lists of anything a predicate can compute.

### Build

- An interval dynamic program that places cuts under a fragment-length bound, and reports the
  fragment count it needed. Exact and cheap, and it answers "how few fragments does this gene
  need?" — the question the Lund table says predicts success.
- An `exclude` argument on `refusal`, so `AGGA`/`TTCC` is stated rather than smuggled in through
  `taken`.
- Fidelity as a ranking key inside the candidate loop, not only as a report afterwards. This is a
  change to `design_overhangs` as much as to anything new.

### Open, and not settled here

- Nothing trades a longer oligo, a different enzyme or a different reserved set against the
  fragment count, and the fragment count is the number Lund measures success against.
- Whether a junction may recode to reach a better overhang. OMEGA never recodes; Lund recodes once,
  before the split; `mbio.codons` could. Allowing it widens every candidate list and changes
  the problem.

### Settled since

Whether the splitter should enforce `MIN_DISTANCE` at all when a measured matrix is available was
open here, on the ground that the matrix is sparse and the sparsity biases the score upward.
**#327 measured it, and it should.** The premise was wrong — the matrix is complete, as §2.5.6 now
records. Two bases is the smallest separation at which every shipped matrix holds its measured
cross-ligations under `MODEST_MISMATCH`, and it is not a free rule: on the 72 AP-1 parts at a
200 nt oligo it takes the worst part from 0.7065 to 0.9829 for eight extra oligos out of 288,
while three bases refuses 32 of those parts outright. It is the lookahead the one-pass greedy
does not have, which is why dropping the branch and bound cost so little.

## 6. Provenance of files read

The branch and bound #261 measured the greedy against lives at the `archive/261-dad-split` tag,
kept there as evidence rather than as a tool.

| File | Source | Read |
| --- | --- | --- |
| `reference_docs/synthesis_and_assembly/dad-split/omega-160be2f/` | `https://github.com/RomeroLab/omega` at `160be2f7292676f64d54dad58902d5de557d9487`, GPL-3.0. Reference only: the licence forbids vendoring it here. | 2026-10-06 |
| `reference_docs/synthesis_and_assembly/Freschlin2026/` | Already held. Used for the data-availability statement that confirms the repository. | 2026-10-06 |
| `reference_docs/synthesis_and_assembly/long_fragment_GGA/Lund2024/` | Already held. Used for the DAD definition, the SplitSet sentence and the fragment-count table. | 2026-10-06 |

## Sources

- Lund, S., Potapov, V., Johnson, S. R., Buss, J. and Tanner, N. A. (2024) Highly parallelized
  construction of DNA from low-cost oligonucleotide mixtures using Data-optimized Assembly Design
  and Golden Gate. *ACS Synth. Biol.* 13(3), 745–751.
  [doi:10.1021/acssynbio.3c00694](https://doi.org/10.1021/acssynbio.3c00694). Read 2026-10-06.
- Freschlin, C. R., Yang, K. K. and Romero, P. A. (2026) Scalable and cost-efficient custom gene
  library assembly from oligopools. *Sci. Adv.* 12, eady2279.
  [doi:10.1126/sciadv.ady2279](https://doi.org/10.1126/sciadv.ady2279). Read 2026-10-06.
- RomeroLab/omega, commit `160be2f7292676f64d54dad58902d5de557d9487`.
  [github.com/RomeroLab/omega](https://github.com/RomeroLab/omega), GPL-3.0; archived at
  [doi:10.5281/zenodo.17637682](https://doi.org/10.5281/zenodo.17637682). Read 2026-10-06.
- Pryor, J. M., Potapov, V., Kucera, R. B., Bilotti, K., Cantor, E. J. and Lohman, G. J. S. (2020)
  Enabling one-pot Golden Gate assemblies of unprecedented complexity using data-optimized assembly
  design. *PLoS One* 15(9), e0238592.
  [doi:10.1371/journal.pone.0238592](https://doi.org/10.1371/journal.pone.0238592). The fidelity
  definition and the shipped matrices; see `docs/research/ligation-fidelity.md`.
- Sikkema, A. P., Tabatabaei, S. K., Lee, Y. J., Lund, S. and Lohman, G. J. S. (2023)
  High-complexity one-pot Golden Gate assembly. *Curr. Protoc.* 3(9), e882. Cited by Lund as a DAD
  source; not read in full here.
- Pryor, J. M., Potapov, V., Bilotti, K., Pokhrel, N. and Lohman, G. J. S. (2022) Rapid 40 kb
  genome construction from 52 parts through data-optimized assembly design. *ACS Synth. Biol.*
  11(6), 2036-2042. [doi:10.1021/acssynbio.1c00525](https://doi.org/10.1021/acssynbio.1c00525),
  PMC9208013, CC BY-NC-ND 4.0. Read 2026-10-06 for whether it states the search; it does not.
- NEBridge SplitSet Tool, SplitSet Lite API and Overhang Optimizer Code. NEB web pages and request
  forms at `ligasefidelity.neb.com/splitset/` and `neb.com/en-us/forms/`, read 2026-10-06. Read
  in a browser: `neb.com` returns 403 to a plain fetch, as
  `reference_docs/synthesis_and_assembly/README.md` already records.

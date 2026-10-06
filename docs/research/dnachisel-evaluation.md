---
search:
  exclude: true
---

# DnaChisel and the Edinburgh Genome Foundry stack: what to adopt, and what not to

Research note for issue #234, a sub-issue of #225. Everything below was retrieved on
**2026-10-06** from the package indexes and the source repositories themselves. It records what
DnaChisel would replace in this package function by function, what else the same group ships,
whether any of it can be a pixi dependency here, what coordinate convention it uses, and the
patterns worth taking without the dependency.

**Nothing was installed and no dependency was added.** DnaChisel 3.2.16 was read at `master`;
the version read is named wherever a claim rests on code.

## 1. The verdict, first

**Adopt nothing. Take three ideas.**

DnaChisel is adoptable on the packaging merits — bioconda `noarch`, MIT, no `python` pin, every
core dependency present for both platforms (§4). That is the part most evaluations stop at, and
it is not the part that decides this.

It fails on three counts that matter more:

1. **Its codon data is the table this repo already rejected.** DnaChisel reads
   `python_codon_tables`, whose nine bundled tables are Kazusa's. `docs/research/codon-usage.md`
   §3.2 rejected Kazusa for *E. coli* K-12 on the measurement: the entry is **14 coding
   sequences** with a zero cell, against the 4,317 this package counts from the genome.
   Adopting `CodonOptimize` means adopting that table, or passing our own and using DnaChisel
   for arithmetic we already have.
2. **It would replace one function and a half.** `translate.optimize_protein` has a
   counterpart. `translate.optimize_coding_sequence` and `sites.domesticate` do not, because
   both are defined by a restraint DnaChisel has no way to express: *a sequence already coded is
   checked, never written again*, and *a site outside a coding sequence is reported, not edited*
   (§3).
3. **It is not being released.** No release and no commit on `master` since **2025-05-10**,
   seventeen months. CI runs Python 3.12 on ubuntu only — no macOS runner, no 3.13 (§4.4).

The three ideas worth taking are in §6. None of them needs the dependency.

## 2. What was asked, and the short answers

| The ticket's question | The answer |
| --- | --- |
| What would DnaChisel replace? | `reverse_translate` and about half of `optimize_protein`. Nothing in `sites` (§3) |
| Can it be a pixi dependency at all? | Yes. bioconda `noarch`, MIT, `python` unpinned, both platforms (§4) |
| What coordinates? | **0-based and half-open, and a span across a circular origin ends past the length — the same rule as ADR 0001** (§5) |
| Does a sibling already decompose a sequence into orderable fragments? | **Yes, two do — and neither solves the DAD split, and neither is installable here.** §7, which #236 is blocked on |
| Patterns worth stealing? | Three, in §6 |

## 3. Function by function: what would be replaced

Read against `src/liulab_mbio/` at `818e4e1`.

### 3.1 `codons` — replaced in form, not in substance

| mbio | DnaChisel's counterpart | Verdict |
| --- | --- | --- |
| `CodonUsage`, `codon_usage(name)`, `codon_tables()` | `python_codon_tables.get_codons_table(name)` | **Do not replace.** Different data, and ours is the one this repo chose |
| `amino_acid(codon)` | Biopython's `CodonTable` | Already shared; both call Biopython |
| `CodonUsage.fraction`, `.per_thousand`, `.synonymous` | internal to their optimiser, not public API | Nothing to call |

The data is the whole story. `python_codon_tables` is CC0 and its README states plainly: *"All
tables are from kazusa.or.jp"*, citing Nakamura, Gojobori and Ikemura (2000). It bundles nine
organisms as CSV — *B. subtilis*, *C. elegans*, *D. melanogaster*, *E. coli*, *G. gallus*,
*H. sapiens*, *M. musculus*, *M. musculus domesticus*, *S. cerevisiae*.

`docs/research/codon-usage.md` already read Kazusa and refused it, on a measurement rather than
a licence: its *E. coli* K-12 entry is computed over 14 coding sequences. Every shipped table
here is counted from a genome instead. Taking `CodonOptimize` at its default would silently
reverse that decision.

There is a second hazard in the same package. `get_codons_table` branches on its argument:

```python
if isinstance(table_name, int) or str.isdigit(table_name):
    return download_codons_table(taxid=table_name, timeout=web_timeout)
```

A **name** reads a bundled CSV offline; an **integer or digit-string TaxID** fetches
`http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi` over plain HTTP at call time and scrapes
the HTML. DnaChisel's own `CodonOptimize` docstring repeats the warning. A design package that
reaches the network mid-design is not something to ship behind our CLI, and the test suite runs
with `filterwarnings = ["error"]` and no network.

### 3.2 `translate` — one function replaced, two not

| mbio | DnaChisel's counterpart | Verdict |
| --- | --- | --- |
| `translate(dna)` | `dnachisel.biotools.translate` | Both wrap Biopython. Nothing gained |
| `reverse_translate(protein, host=)` | `dnachisel.biotools.reverse_translate` | **Genuine overlap.** About 20 lines |
| `optimize_protein(...)` | `CodonOptimize` + `EnforceTranslation` + `AvoidPattern` | Partial. See below |
| `optimize_coding_sequence(...)` | none | **No counterpart** |
| `CodingSequence`, `Domestication` | `SpecEvaluation.locations` | Ours names the codon that moved; theirs names a span |
| `SiteNotRemovableError` | `NoSolutionError` | Theirs carries the whole failed-constraint summary. Worth reading (§6) |

`optimize_protein` is the one real candidate. DnaChisel would express it as:

```python
DnaOptimizationProblem(
    sequence=reverse_translate(protein),
    constraints=[EnforceTranslation(), AvoidPattern("BsaI_site")],
    objectives=[CodonOptimize(species="e_coli")],
)
```

and its solver is better than ours at the hard case: `EnforceTranslation.restrict_nucleotides`
builds a `MutationSpace` of per-codon synonymous variants, so every later search is
*structurally* unable to change the protein, and site removal and codon choice happen in one
space rather than one after the other. Where our `_synonymous` changes one codon and gives up if
no synonym clears the site, theirs widens the window by 5 bp and tries again
(`local_extensions = (0, 5)`). That is a real capability we do not have.

It buys that at the cost of the two functions below it.

**`optimize_coding_sequence` has no DnaChisel shape.** Its contract is the opposite of a
solver's: *a codon spelling no forbidden site is left exactly as it stands*. DnaChisel's
`CodonOptimize` is an **objective**, and an objective's job is to move every codon it can
improve. `AvoidChanges` can freeze regions, but freezing everything except the codons under a
site requires knowing the sites first — at which point the solver is doing nothing our code is
not already doing.

### 3.3 `sites` — nothing is replaced

| mbio | DnaChisel's counterpart |
| --- | --- |
| `find_sites`, `has_site`, `site_counts`, `free_enzymes` | **None.** `AvoidPattern` returns a score and spans, never a cut |
| `CutSite` with `top_cut`, `bottom_cut`, `overhang`, `certain` | **None** |
| `digest`, `Fragment` with its two overhangs | **None** |
| `insert_site`, `primer_tail` | `EnforcePatternOccurence` places a site; it designs no primer tail |
| `domesticate` | `AvoidPattern` + `EnforceTranslation` |
| `DomesticationReport.outside_cds` | **None** |

This is the decisive table. DnaChisel models a recognition site as **a pattern to avoid**. This
package models it as **an enzyme that cuts at an offset, leaving an overhang** — which is what
every cloning pipeline downstream of `sites` actually consumes. `digest` and `Fragment` have no
counterpart anywhere in the stack, and `CutSite.overhang` is what `overhangs` and `ligase` score.

Their enzyme data is Biopython's `Restriction_Dictionary`:

```python
from Bio.Restriction.Restriction_Dictionary import rest_dict

class EnzymeSitePattern(DnaNotationPattern):
    def __init__(self, enzyme_name):
        self.enzyme_site = rest_dict[enzyme_name]["site"]
```

Recognition sequence only. No cut offset in the pattern, so the overhang cannot be derived from
it. `docs/research/restriction-enzyme-data.md` records why this package carries REBASE offsets
and supplier facts instead.

`domesticate`'s refusal is also not expressible. It reports a site lying in no coding sequence
and **does not edit it**, because removing it would change what the record spells. A DnaChisel
constraint edits wherever the mutation space allows; it has no notion of "I could fix this, and
declining is the right answer."

## 4. Whether it can be a dependency here at all

Verified from the indexes on 2026-10-06, not from memory.

### 4.1 Availability: yes, on bioconda, `noarch`

| Question | Answer | Source |
| --- | --- | --- |
| conda-forge? | **No.** `api.anaconda.org/package/conda-forge/dnachisel` returns 404; no feedstock; no staged-recipes PR, open or closed | anaconda.org API, GitHub API |
| bioconda? | **Yes.** 3.2.16, uploaded 2025-05-10 | `api.anaconda.org/package/bioconda/dnachisel` |
| Which subdirs? | **`noarch` only** — so osx-arm64 and linux-64 both resolve | the same, `conda_platforms: ["noarch"]` |
| Python pin? | **None.** `depends` is `[biopython, docopt, flametree, numpy, proglog, python, python-codon-tables]`, `python` unversioned | the per-file `attrs` |
| PyPI wheel | `dnachisel-3.2.16-py3-none-any.whl`, pure python, `requires_python: null` | `pypi.org/pypi/dnachisel/json` |

Nobody ever submitted a conda-forge recipe, so its absence is not a rejection. bioconda is
already a channel here, for `ipcr`. **A pip fallback is not forced and no platform fails.**

### 4.2 Transitive weight: light, if the extras are left alone

Six core dependencies, all on conda-forge or bioconda for both platforms with Python 3.13
builds: `numpy`, `biopython` (both already ours), `proglog`, `flametree`,
`python-codon-tables`, `docopt`.

Two notes on those. `docopt` is **0.6.2, uploaded 2014-06-16** and abandoned upstream; the
conda-forge `noarch` build was re-rendered in 2024 with `python >=3.9`, so the solver accepts it
on 3.13, but it is a decade dead and DnaChisel uses it only for its own CLI. `flametree` and
`python-codon-tables` were last built in 2021 and 2022.

The `[reports]` extra is where the weight is, and the bioconda package does **not** pull it —
its `run` list is core-only. It would drag `weasyprint`, and with it `pango`, `glib` and that
whole native subtree, for PDF rendering this package does not need. Worse, `pdf-reports` pins
`jinja2 <3.1.0`, whose newest satisfying conda-forge build is **3.0.3 from 2021-11-10**. Taking
the extra would hold a shared environment at a five-year-old jinja2. If DnaChisel is ever taken,
the extra must not be.

### 4.3 Licence: MIT, verified from the file

`raw.githubusercontent.com/.../DnaChisel/master/LICENSE` reads `MIT License / Copyright (c) 2017
Edinburgh Genome Foundry, University of Edinburgh`. `pyproject.toml`, PyPI and bioconda all
agree. Compatible with this package's MIT.

`python-codon-tables` is **CC0-1.0** by its own `LICENSE` and `pyproject.toml`. bioconda's
metadata says `NIST-PD`, which is wrong; the file is authoritative. Nothing blocks
redistribution of either.

### 4.4 Maintenance: the real risk

| Signal | Value |
| --- | --- |
| Last release | **v3.2.16, 2025-05-10** — seventeen months ago |
| Last commit on `master` | **2025-05-10** (`68c0930`) |
| Last commit on any branch | 2025-05-15 (`dev`) |
| Open issues | 19, with five filed in 2026 and none closed by a release |
| CI | **Python 3.12 only, ubuntu-24.04 only.** No macOS runner, no arm runner, no 3.13 |
| Stars | 281 |

Its metadata permits Python 3.13 and it is pure Python, so it very probably works. **Nobody has
tested it there**, and a search of its issues for "3.13" returns nothing — neither a breakage
report nor a success. Taking it would make this repo the first to find out, on a platform
(osx-arm64) its CI has never run on.

That is survivable for a dependency that removes a lot of code. It is not worth it for one that
removes a reverse-translation helper.

## 5. The coordinate hazard: there isn't one

**DnaChisel is 0-based and half-open, and an origin-spanning location ends past the sequence
length — the same rule as `docs/adr/0001-coordinates.md`.** This was the question most likely to
sink adoption, and it is the question that came back cleanest.

The `Location` docstring says so outright:

> Warning: we use Python's splicing notation, so `Location(5, 10)` represents
> `sequence[5, 6, 7, 8, 9]` which corresponds to nucleotides number 6, 7, 8, 9, 10.

Confirmed three more ways in the same file: `__len__` returns `self.end - self.start`,
`extract_sequence` slices `sequence[self.start : self.end]`, and `indices` is
`range(self.start, self.end)`. Overlap is half-open too.

Their own circular example carries the origin rule:

```python
constraints = [
    dc.AvoidPattern("BsmBI_site"),
    dc.EnforceGCContent(mini=0.4, maxi=0.6, location=(1500, 2500), window=50),
]
problem = dc.CircularDnaOptimizationProblem(sequence=dna_sequence, constraints=constraints)
```

with a sequence of about 2,006 bp. `location=(1500, 2500)` crosses the origin and ends past the
length, exactly as `Segment(2683, 2689)` does on a 2,686 bp plasmid in ADR 0001. The file
boundary is the same discipline as ours, too: `to_biopython_location` and
`from_biopython_location` apply **no offset at all**, and the 1-based GenBank conversion is left
entirely to Biopython's `SeqIO` — one boundary, one layer down.

So the conversion cost the ticket asked about is **zero**. Three caveats to carry anyway:

- **`Location` is origin-unaware.** No modulo, no wrap, no length. There is no circular location
  type; the whole of circularity lives in the problem class.
- **Circularity is implemented by triplication.** `CircularDnaOptimizationProblem` builds a
  problem on `3 * sequence`, shifts every spec into all three copies, solves, and takes the
  middle third back — reconciling the copies by majority vote in a function named
  `return_the_loony`. Its own docstring calls the class an "Attempt" and carries a standing
  TODO. `location=(1500, 2500)` works only because index 2500 is real in the tripled view.
- **Plain `DnaOptimizationProblem` has no circular support.** A cross-origin site is simply
  invisible to it, which our `find_sites` is not.

ADR 0001 needs no change. DnaChisel is usable as corroboration that the convention is the one a
sibling project reached independently. Its *circular* implementation is the weaker design of the
two and should not be copied.

## 6. Patterns worth taking, with no dependency

The ticket asks for these as a suggestion. Reported, not proposed — none of them is a
restructuring, and the first two bear on #227, which is deciding what a `Check` carries.

### 6.1 A signed score with a known best, so "passes" and "optimal" are one number

`SpecEvaluation.__init__` derives both flags from the score:

```python
self.passes = score >= 0
self.is_optimal = score == specification.best_possible_score
```

Every specification declares `best_possible_score` as a class attribute — for `AvoidPattern` it
is `0`, and the score is `-len(locations)`. One number answers two questions, and it is what
lets their solver know when to stop.

Our `Check` carries `status` and `value`. It can say *pass*; it cannot say *passes, and is two
short of the best achievable here*. Our thresholds come from sources rather than from a solver,
so the fit is not automatic — but the distinction between "acceptable" and "best available" is
exactly the one `primers/thresholds` has to draw by hand today.

### 6.2 A failure carries where it failed, not just that it failed

`SpecEvaluation` carries seven fields: the specification, the problem, the score, `passes`,
`is_optimal`, **`locations`**, a message, and a free-form `data` dict. `locations` is the list of
spans that breached.

```python
def evaluate(self, problem):
    """Return score=-number_of_occurences. And patterns locations."""
    locations = self.pattern.find_matches(problem.sequence, self.location)
    score = -len(locations)
    if score == 0:
        message = "Passed. Pattern not found !"
    else:
        message = "Failed. Pattern found at positions %s" % locations
    return SpecEvaluation(self, problem, score, locations=locations, message=message)
```

That is what makes a failure *actionable* rather than merely reported: the solver localises its
next search on those spans, and the report paints them onto the GenBank output as red
`misc_feature`s. Our `Check` has `detail`, a string. A span would let a check point at a map.

Three smaller touches in the same file worth noting. The rendered failure line identifies the
specification by a reconstructable constructor-like label, not a free-form name:

```text
FAIL ┍ UniquifyAllKmers[10-1000](k:9)
     │ Score:        -2. Locations: [232-241, 233-242]
```

A long message is truncated **in the middle** — `message[:half] + " ... " + message[-half:]` —
so a breach at 400 sites still shows its first and its last. And one renderer serves three
outputs: `constraints_text_summary()` is the report, the console, and the body of the
`NoSolutionError` that is raised when resolution fails. Our `SiteNotRemovableError` builds its
message by hand in `_why`.

### 6.3 Constraint against objective is a position, not a type

The same `Specification` subclass can sit in `constraints=[...]` or in `objectives=[...]`. What
changes is whether a breach is fatal or merely costly, and the vocabulary it renders with —
constraints print `✔PASS` / `FAIL`, objectives print `✔` for optimal.

This package draws the same line in a different place, and `CONTEXT.md`'s rule that a check with
no sourced threshold carries `None` rather than a pass is a third category theirs expresses as
`enforced_by_nucleotide_restrictions`, which short-circuits to a pass with the message *"Enforced
by nucleotides restrictions"*. Worth knowing the shape exists before #227 settles.

**Not worth taking:** the `@AvoidPattern(...)` / `~CodonOptimize(...)` GenBank label syntax, by
which a whole problem is defined in a record's `misc_feature` labels. It is a neat idea, and
their own `to_biopython_feature` docstring disowns the round trip: *"this is not the intended
goal."* ADR 0002 already settled that an editable protocol JSON is where an agent intervenes.

## 7. The rest of the stack, and the answer #235 and #236 need

### 7.1 Does a sibling already decompose a sequence into orderable fragments?

**Yes. Two of them do — DnaWeaver and GoldenHinges. Neither solves the DAD split, and neither
is adoptable. #236 does not shrink to wiring.**

Issue #236 is blocked on this note for exactly this question, so it is answered in full.

#### DnaWeaver: the right shape, the wrong objective

DnaWeaver is "a route planner for DNA assembly". Its README:

> Given an arbitrary sequence, DNA Weaver will select the most adapted commercial DNA providers,
> cloning methods and parts repositories (depending on your preferences), and will design all
> necessary assembly fragments and assembly steps.

The algorithm is worth knowing, because it is a good one. `SequenceDecomposer`:

> Find the sequence cuts which optimize the sum of segment scores. […] a Networkx DiGraph whose
> nodes are indices of locations in the sequence and the 'weight' of edge `[i][j]` is given by
> `segment_score_function(segment[i][j])`.

A **shortest path over a graph of candidate cut positions**, solved by Dijkstra, or by A\* with a
heuristic its own docstring warns "yields suboptimal decompositions". A coarse pass then a local
refinement around each coarse cut, which is what makes a long sequence tractable. An assembly
station is itself a supplier, so stations nest and hierarchical plans fall out of one mechanism.

**The edge weight is a price.** In `DnaAssemblyStation.new_sequence_decomposer`, `segment_score`
returns `quote.price`; lead time is a hard constraint that prunes quotes, not a second objective.
The README's "price, duration, or assembly success probabilities" is reached by changing what the
user's pricing function returns, not by a multi-objective solver.

It is more junction-aware than that makes it sound, and the ticket deserves the correction:
`GoldenGateAssemblyMethod` registers a **`cuts_set_constraints`** — a constraint on the whole
chosen set, requiring every pair of overhangs to differ by `min_overhangs_differences` and to
differ from each other's reverse complement. So it does reason about the set, not only about one
junction.

**What it does not have is the fidelity model.** Its junction rule is Hamming distance on the
four bases. Nothing from Potapov 2018 or Pryor 2020 — those live in separate EGF packages
(`tatapov`, `Overhang`, `kappagate`) that DnaWeaver does not depend on. `src/liulab_mbio/overhangs.py`
and `ligase.py` are strictly stronger here.

That is the half #235 says is missing, restated: *"the scoring half exists; the search half does
not."* DnaWeaver is a search over the wrong cost. Converting it would mean replacing the edge
weight, which is the whole of what it does.

#### GoldenHinges: the right objective, and effectively uninstallable

Closer prior art, and the ticket should see it. GoldenHinges will

> decompose a DNA sequence into fragments with compatible overhangs for scar-less DNA assembly.
> Can also suggest minimal edits in the sequence to allow such assembly.

`OverhangsSelector.cut_sequence(sequence, equal_segments=50, max_radius=20)` cuts a long sequence
into size-balanced fragments whose 4 bp overhangs are mutually compatible, and writes out the
sequences to order with the restriction sites already flanking them. Three features have no
counterpart anywhere else read for this note:

- **`allow_edits=True`** mutates the sequence to make a decomposition possible when none exists,
  penalised so no base changes unless strictly necessary — and it honours DnaChisel annotations
  (`@AvoidChanges`, `@EnforceTranslation`) to forbid an illegal edit. #236's question 3, *what
  does it do when it cannot win*, has a published answer here.
- A `!cut` GenBank feature forces a cut in a named region.
- Its README documents feeding it **tatapov** (Potapov 2018) data to exclude weak-annealing
  overhangs and high cross-talk pairs — as a user-assembled `forbidden_overhangs` and
  `forbidden_pairs` list, not built in. The held-out set #235 needs for `AGGA` and `TTCC` is
  exactly that shape.

**And it cannot be taken.** Last commit **2022-06-08**, four years. Not on bioconda, not on
conda-forge, PyPI only. It depends on **Numberjack**, a SWIG-based constraint solver whose last
release is **1.2.0 from 2016**, which is LGPL-2.1, and whose install instructions tell you to
downgrade SWIG to v3 and build from an EGF fork. Under `CLAUDE.md`'s pixi-only rule on osx-arm64,
that is the end of it.

#### So what #236 should take

Not a dependency. Two framings, both MIT and free to read:

- **DnaWeaver's**: decomposition as a shortest path over candidate cut positions, with a cost per
  candidate fragment. Score the path on fidelity rather than on price — `overhangs` already holds
  the scorer the path would call. The caveat is in #236's own brief: fidelity is a property of
  the whole set in one pot, and a shortest path over independent edge costs cannot express a
  pairwise constraint across junctions. DnaWeaver meets that with a set-level
  `cuts_set_constraints` check; whether that is enough is what the prototype is for.
- **GoldenHinges'**: overhang-set compatibility as constraint satisfaction, and the
  edit-to-enable fallback when no clean cut set exists.

**Neither answers #235's question 1**, which is what SplitSet implements and whether it is the
method Lund cites. Nothing in this stack is SplitSet.

### 7.2 The rest, briefly

| Library | What it is | State | Bears on us? |
| --- | --- | --- | --- |
| **DnaWeaver** | decomposition over a supply network | MIT, bioconda `noarch` **v0.3.7, three releases behind PyPI**, last commit 2025-05-12 | §7.1. The framing, not the package |
| **GoldenHinges** | overhang-compatible decomposition, with edits | MIT, **PyPI only**, last commit 2022-06-08, needs Numberjack | §7.1. Read it; cannot install it |
| **DnaCauldron** | cloning simulator, Golden Gate and Gibson | MIT, bioconda `noarch` 2.0.12, last commit 2025-05-12 | Overlaps `cloning/`. See below |
| **tatapov** | the Potapov and Pryor fidelity matrices | MIT code, **CC BY-ND 4.0 data** | See below — it corroborates our rule |
| **Overhang** | describes and scores an overhang set | MIT, PyPI only, last commit 2025-05-15 | `overhangs` covers it |
| **kappagate** | predicts the good-clone fraction for 4 bp overhangs | MIT, **last commit 2020-07-24** | Adjacent to `library/coverage`. Dead |
| **genedom** | domestication to an assembly standard | MIT, bioconda 0.2.2 | Adjacent to `sites.domesticate` |
| **GeneBlocks** | common blocks and diffs between sequences | MIT, PyPI only, last commit 2025-11-05 | Adjacent to `edits`. Not a splitter |
| **DnaFeaturesViewer** | map drawing, matplotlib | MIT, 703 stars, PyPI only | `plot/` writes its own SVG. No |
| **Primavera**, **BandWitch**, **BandWagon** | verification primers, enzyme choice, gel bands | MIT | Adjacent to `bench/validation` and `bench/gels` |
| **Plateo**, **blabel**, **Caravagene** | plates, labels, construct schematics | MIT | No |

**One finding worth carrying out of this note.** `tatapov` reads the NEB fidelity data from a
separate repository, `tatapov_data`, whose README states its licence is **CC BY-ND 4.0** — no
derivatives — and `tatapov` downloads it at first use rather than vendoring it. `CLAUDE.md`'s
rule, *"ship nothing whose licence forbids it — a ligase matrix is read from a copy the user
holds"*, is the same wall, met independently, and ours is the better answer for a package that
must work offline. `docs/research/ligation-fidelity.md` §8 can cite this.

**Two gaps in the whole portfolio that this package already fills.** There is **no EGF barcode
design** of any kind — a repo search finds only a label printer. And **nothing generates a bench
protocol**: their output idiom is a PDF report of what was computed, never numbered steps,
rescaling reaction tables or a thermocycler program. `build-protocol` has no counterpart.

**DnaCauldron deserves one extra line**, as the closest thing to a duplicate of `cloning/`. It
simulates Golden Gate, Gibson, BioBrick, BASIC and LCR assemblies by finding circular paths in a
part homology graph, after Pereira et al. 2015. It does not design the primers, choose the
method, score the overhang set against a ligase profile, or write a protocol — which is most of
what a plan is here. It is the best-packaged of the stack, and the one worth remembering if we
ever want to cross-check our simulation against an independent one.

**The whole stack shares one release date, and it is not what it looks like.** DnaWeaver,
DnaCauldron, `genedom`, `easy_dna`, `BandWagon` and others were all last touched on 2025-05-12,
with the same commit message — `"Include pkg data, exclude unwanted folders"`. That is one
packaging sweep across a dozen repositories, not twelve healthy projects. Open bugs on DnaWeaver
date from 2021 and are unanswered. Read the stack as **prior art, not as dependencies that will
receive fixes.**

## 8. Recommendation, and the trade

**Adopt nothing. Take the three ideas in §6, which cost no dependency.**

The trade, stated plainly:

**What is given up.** A better codon solver. DnaChisel's `MutationSpace` searches codon choice
and site removal in one space and widens a stuck window by 5 bp, where `translate._cleaned`
changes one codon and raises `SiteNotRemovableError` if no synonym clears the site. A construct
that needs two codons moved together will fail here and succeed there. That is a real
capability, and it is the honest cost of this recommendation.

**What is kept.** The codon tables this repo measured rather than copied, and the reason #10
rejected Kazusa. A `domesticate` that declines to edit outside a coding sequence. A site model
that carries a cut and an overhang rather than a pattern. No seventeen-month-stale dependency,
no abandoned `docopt`, no package whose CI has never run on macOS or Python 3.13.

**When to reopen this.** If `SiteNotRemovableError` starts firing on real designs — a
measurement, not a hypothesis — the answer is not necessarily DnaChisel. It is a multi-codon
search in `translate`, which is perhaps 60 lines and is the part of their solver that earns its
keep. `CLAUDE.md`: *generalizable, lightweight, uncustomized — in that order.* A dependency
carrying nine Kazusa tables, a web scraper and a 2014 CLI parser, to supply 60 lines of search,
is not lightweight.

## 9. Open gaps

Not guessed above.

| Item | Why it is missing |
| --- | --- |
| Whether DnaChisel runs correctly on Python 3.13 | Nothing was installed, and upstream has never tested it. Metadata permits it; no one has confirmed it |
| Whether it runs on osx-arm64 | Same. Pure Python, so very probably, but its CI is ubuntu-only |
| How `CircularDnaOptimizationProblem`'s majority vote behaves on a real plasmid | Read, not run |
| Lund et al. 2024's actual DnaChisel invocation | Named in the ticket; their parameters were not traced to a script in this note |

## Sources

All retrieved 2026-10-06.

- Zulko, *DnaChisel*, Edinburgh Genome Foundry, version 3.2.16 at `master`:
  [github.com/Edinburgh-Genome-Foundry/DnaChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel)
  (MIT, read from `LICENSE`). `dnachisel/Location.py`,
  `dnachisel/DnaOptimizationProblem/CircularDnaOptimizationProblem.py`,
  `dnachisel/DnaOptimizationProblem/mixins/ConstraintsSolverMixin.py`,
  `dnachisel/Specification/SpecEvaluation/SpecEvaluation.py`,
  `dnachisel/SequencePattern/EnzymeSitePattern.py`,
  `dnachisel/builtin_specifications/`, `examples/common_scenarios/circular_sequence.py`
- *DnaChisel* documentation:
  [edinburgh-genome-foundry.github.io/DnaChisel](https://edinburgh-genome-foundry.github.io/DnaChisel/)
  — the built-in specifications reference and the GenBank API page
- Zulko, *codon-usage-tables* (`python_codon_tables`):
  [github.com/Edinburgh-Genome-Foundry/codon-usage-tables](https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables)
  (CC0-1.0, read from `LICENSE`), for the bundled organisms and the TaxID download path
- Nakamura, Y., Gojobori, T. and Ikemura, T. (2000) Codon usage tabulated from international DNA
  sequence databases. *Nucleic Acids Res.* 28, 292 — the source `python_codon_tables` cites
- anaconda.org API: `api.anaconda.org/package/bioconda/dnachisel` and `/files`, and the
  conda-forge 404, for availability, subdirs, `depends` and upload dates
- bioconda-recipes, `recipes/dnachisel/meta.yaml`, for `noarch: python` and the run list
- PyPI JSON API: `pypi.org/pypi/dnachisel/json`, for the wheel tag, `requires_python` and
  `requires_dist`
- GitHub API: `api.github.com/repos/Edinburgh-Genome-Foundry/DnaChisel` and its `/releases`,
  `/commits` and `/branches`, for the maintenance signals; `.github/workflows/build.yml` for
  what CI runs
- `docs/research/codon-usage.md` (issue #10), for the Kazusa decision this note holds to, and
  `docs/research/restriction-enzyme-data.md` (issue #7), for why the enzyme data here carries
  cut offsets
- Zulko, *DnaWeaver*:
  [github.com/Edinburgh-Genome-Foundry/DnaWeaver](https://github.com/Edinburgh-Genome-Foundry/DnaWeaver)
  (MIT, read from `LICENSE`), its `README.rst`,
  `dnaweaver/DnaSupplier/builtin_suppliers/DnaAssemblyStation/SequenceDecomposer.py`,
  `dnaweaver/DnaAssemblyMethod/GoldenGateAssemblyMethod.py` and `dnaweaver/SegmentSelector/`
- Zulko, *GoldenHinges*:
  [github.com/Edinburgh-Genome-Foundry/GoldenHinges](https://github.com/Edinburgh-Genome-Foundry/GoldenHinges)
  (MIT, read from `LICENCE.txt`), its `README.rst`, for `cut_sequence`, `allow_edits` and the
  tatapov-fed forbidden lists; and Numberjack's PyPI record for the 2016 release that blocks it
- Zulko, *DnaCauldron*:
  [github.com/Edinburgh-Genome-Foundry/DnaCauldron](https://github.com/Edinburgh-Genome-Foundry/DnaCauldron)
  (MIT), its `README.rst` and `dnacauldron/Assembly/builtin_assembly_classes/`
- Pereira, F. et al. (2015) *BMC Bioinformatics* 16:90, the circular-path method DnaCauldron
  names: [doi:10.1186/s12859-015-0544-x](https://doi.org/10.1186/s12859-015-0544-x)
- *tatapov* and *tatapov_data*:
  [github.com/Edinburgh-Genome-Foundry/tatapov_data](https://github.com/Edinburgh-Genome-Foundry/tatapov_data)
  — its README states **CC BY-ND 4.0**, which §7.2 reads against `CLAUDE.md`'s ligase-matrix rule
- The Edinburgh Genome Foundry organisation listing (76 public repositories) and its portfolio
  page source, for which libraries exist and which are not relevant
- `docs/adr/0001-coordinates.md`, the rule §5 compares against

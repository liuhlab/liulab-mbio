---
search:
  exclude: true
---

# Ligation fidelity data: sources, licence, and how a set of overhangs is scored

Research note for issues #11 and #18. Everything below was retrieved on **2026-09-12**. It records where
`src/mbio/data/ligation_fidelity.json` comes from, what its licence allows, how the
shipped matrices are read, and which overhang rules the code applies and on whose authority.

Issue #4's note (`docs/research/golden-gate-assembly.md`, section 4) surveyed the field and drew
the licence line. This note carries out what it decided and records what could be checked
against the paper's own published numbers.

## 1. What the data file holds

One file, `ligation_fidelity.json`, read with `importlib.resources` by
`mbio.overhangs`. It holds five count matrices, one per enzyme, plus the
attribution every copy of the data has to carry.

| Field | Meaning |
| --- | --- |
| `source` | citation, DOI, download address, retrieval date, copyright, licence, licence text |
| `format` | what the matrix is, stated in the file so a reader need not find the code first |
| `matrices[].enzyme` | the name this package ships the enzyme under |
| `matrices[].product` | the supplier's product the measurement used, such as BsaI-HFv2 |
| `matrices[].table`, `.file`, `.sheet` | which supplementary table, file and sheet it came from |
| `matrices[].overhang_length` | 4, or 3 for SapI |
| `matrices[].cycling_celsius` | the two temperatures that reaction was cycled between |
| `matrices[].observations` | every ligation event the matrix counts |
| `matrices[].counts` | the matrix itself, sparse |

The five, as built:

| Table | Product | Overhangs | Cycling | Observations | Non-zero cells |
| --- | --- | --- | --- | --- | --- |
| S1 | BsaI-HFv2 | 256 of 4 nt | 37 / 16 °C | 203,364 | 4,005 |
| S2 | BsmBI-v2 | 256 of 4 nt | 42 / 16 °C | 219,378 | 4,130 |
| S3 | Esp3I | 256 of 4 nt | 37 / 16 °C | 224,222 | 4,051 |
| S4 | BbsI-HF | 256 of 4 nt | 37 / 16 °C | 189,282 | 4,230 |
| S5 | SapI | 64 of 3 nt | 37 / 16 °C | 172,332 | 739 |

The full matrices are 256 x 256 and 64 x 64, which is 266,240 cells holding 17,155 non-zero
counts. **The shipped file is sparse**: a pair never observed is absent rather than written as
a zero, and each row is written on one line so a change to the data reads as a change to a row.
That is the whole of the format, and the file says so in its own `format` field.

## 2. Licences

### Pryor et al. 2020: CC BY 4.0, so it ships

The article's own permissions block, from the JATS XML at `journals.plos.org`:

> This is an open access article distributed under the terms of the Creative Commons
> Attribution License, which permits unrestricted use, distribution, and reproduction in any
> medium, provided the original author and source are credited.

The licence is `http://creativecommons.org/licenses/by/4.0/` and the copyright line is
`Copyright (c) 2020 Pryor et al`. Attribution is the only condition, and it is met three times
over: in the data file's `source` block, in the module docstring of
`mbio/overhangs.py`, and here.

> Pryor, J.M., Potapov, V., Kucera, R.B., Bilotti, K., Cantor, E.J. and Lohman, G.J.S. (2020)
> Enabling one-pot Golden Gate assemblies of unprecedented complexity using data-optimized
> assembly design. *PLoS One* 15(9): e0238592.
> [doi:10.1371/journal.pone.0238592](https://doi.org/10.1371/journal.pone.0238592)

**Verdict: S1-S5 Tables ship**, converted to the sparse form above. The conversion is a format
change, not a new work: every integer is the integer the spreadsheet holds.

### Potapov et al. 2018: CC BY-NC, so it does not

Both 2018 papers and the figshare Supporting Data are CC BY-NC 4.0. A non-commercial term
cannot be honoured by an MIT-licensed package, so **none of that data ships**. It is cited
where it is the authority for something:

- The measurement that the two-mismatch rule is stricter than it needs to be.
- The junction 6 (`GCCG`/`CGGC`) truncation, which is why an overhang of one base kind is
  refused here even though the fidelity data alone would pass it.
- The normalisation to 100,000 ligation events that NEB's thresholds are stated in.
- The ready-made high-fidelity overhang sets of its Table 1, which are short factual lists
  reproduced with citation in a test, not bulk data.

The analysis code at `github.com/potapovneb/ligase-fidelity` is AGPL-3.0 and was not read into
this package.

### NEB: cited, never mirrored

The Ligase Fidelity Viewer help page supplies the two thresholds and the statement of which
axis is which. Those are facts read from a cited page, in the same position issue #7 took for
the enzyme properties: no NEB file, page or table is redistributed.

## 3. Reading the workbooks without a spreadsheet library

An `.xlsx` is a zip of XML, so `scripts/build_ligation_fidelity.py` reads the three parts it
needs with `zipfile` and `xml.etree` and **no dependency is added**: `xl/workbook.xml` for the
sheet name, `xl/sharedStrings.xml` for the labels, and `xl/worksheets/sheet1.xml` for the
cells. Labels are shared strings and counts are numbers; a cell of any other type is refused
rather than guessed at.

Two guards, because a mis-mapped table would be silent and wrong:

- The sheet name must carry the product the table is supposed to measure. The five sheets are
  named `S1 Table. BsaI-HFv2`, `Table S2. BsmBI-v2`, `Table S3. Esp3I`, `Table S4. BbsI-HF`
  and `Table S5. SapI`, so this is a real check and not a formality.
- Every row label must be an overhang of the length that enzyme leaves, which is what catches
  a four-base table served where the three-base one was asked for.

A cell is placed by the letters of its reference, not by its position in the row: a row omits
the cells it has no count for, so counting along the row would shift every column after the
first gap.

## 4. Which way round the matrix is

From NEB's Ligase Fidelity Viewer help page:

> Overhangs in rows correspond to top strand. Overhangs in columns correspond to bottom strand.
> Overhangs are always written 5' to 3'.

So **the Watson-Crick entry of a row is the column spelling its reverse complement**, and the
cell where a row meets the column of the same name is a mispair. Row `TTTT` against column
`AAAA` holds 635 in the BsaI table while row `TTTT` against column `TTTT` holds 4. Getting this
backwards would read every fidelity as near zero, which is why a test pins it.

## 5. Fidelity, as the paper defines it and as the code computes it

Pryor 2020, Methods, verbatim:

> Fidelity F for a set of n overhangs {O1, O2, O3, ..., On}, was defined as a probability that
> all overhangs in the set ligate correctly to their WC pair. Fidelity estimates the fraction
> of correctly ligated products when using a given set of overhangs in Golden Gate assembly and
> was computed as follows: F = p(O1) x p(O2) x p(O3) x p(On), where p(Oi) is the probability of
> overhang Oi ligating correctly to its WC pair in a given set of overhangs.

The paper leaves one thing for a reader to work out: which cells of the matrix p(Oi) is a ratio
of. A junction is one thing with two ends. One end presents its overhang on the top strand, a
row; the other presents the reverse complement on the bottom strand, a column. **Both ends are
counted**, and the set of columns a reaction offers is the set of overhangs closed under
reverse complement, which is also what the Viewer means by "Complementary (reverse) overhangs
are automatically added to the query".

So, for a junction with top-strand overhang O and partner O' = reverse complement of O, over a
query set Q closed under reverse complement:

```text
p(O) = (count[O][O'] + count[O'][O]) / sum over c in Q of (count[O][c] + count[O'][c])
F    = product of p(O) over the junctions
```

### This was checked against the paper's own numbers

Counting only one end of each junction gives 76% where the paper says 81%; counting both gives
81%. Three of the paper's published examples come back exactly, which is the evidence that the
formula, the axis orientation and the parsing are all right together:

| Example | Paper | Computed here |
| --- | --- | --- |
| 11 plant overhangs, BsmBI-v2, 42/16 (Fig 4A) | 81% | 80.93% |
| the same set without the `GGTA`/`TACT` mispair (Fig 4A) | 92% | 92.11% |
| the same set extended to 20 by GetSet (Fig 4B) | 80% | 80.48% |

The paper also names the worst mispair of the plant set, `5'-GGTA`/`5'-TACT`, and that pair is
what this code reports first, at 22.8 ligations per 100,000 events against 8.7 for the next.

**Two of the paper's numbers do not come back exactly**, and neither is guessed at or tuned
for. Both are from the table of GetSet-generated sets used for the lac assemblies rather than
from the worked example:

| Example | Paper | Computed here |
| --- | --- | --- |
| 13 SapI three-base overhangs | 79% | 77.45% |
| 35 BsmBI-v2 overhangs | 65% | 61.06% |

Those sets were generated for an assembly that also carries a destination plasmid junction, so
the reaction they were scored in holds at least one overhang the table does not list. That is
the likeliest reading and it is **not** confirmed, so it is recorded here rather than worked
around. The three exact reproductions are the ones the tests pin.

### The thresholds are on normalised counts

The Viewer's help page calls a Watson-Crick pair **strong** above 100 normalised ligations and
a mismatch **trace** below 10, **modest** between 10 and 100, and **high-count** above 100.
Those numbers are per 100,000 ligation events, the scale Potapov 2018 normalises to, and not
counts in one matrix's own total, which is why `LigationMatrix.normalised` exists.

A measurement worth recording, pinned by a test: **no Watson-Crick pair in any of the five
shipped matrices is weak.** The lowest normalised pair sits just above the threshold, so
`FidelityReport.weak` is empty for every set scored on shipped data. A weak pair would be news.

## 6. The overhang rules, and who says so

Pryor 2020 states the conventional rules and then shows they are the wrong abstraction:

> avoiding use of palindromic overhang sequences or the same overhang pair more than once in an
> assembly reaction. In addition, most modular cloning systems also require that
> non-complementary overhang sequences have at least 2 mispaired bases. Moreover, overhangs
> that contain 100% A/T or G/C content are also often avoided.

| Rule | What it refuses | Source | Dial |
| --- | --- | --- | --- |
| `length` | anything but the enzyme's own overhang length in A, C, G, T | the enzyme record | none |
| `palindrome` | an overhang reading the same on both strands | NEB: overhangs are designed "non-palindromic (to eliminate self insert ligations)" | none |
| `repeat` | an overhang already taken, or its reverse complement | Pryor 2020, above; NEB designs overhangs "unique" | none |
| `near-duplicate` | fewer than two bases different from another overhang or its reverse complement | the modular cloning convention Pryor 2020 states | `min_distance` |
| `uniform` | an overhang that is all G/C or all A/T | Pryor 2020, above; Potapov 2018's junction 6 truncation | `allow_uniform` |
| `site` | an overhang no primer tail can carry without spelling a second site | `mbio.sites.primer_tail` | `avoid` |

**The distance rule is two bases and it is an argument.** Potapov 2018 measured the convention
as stricter than it needs to be:

> it is not necessary to ensure all overhangs have at least two bases different from all other
> overhangs, as many pairs with only a single base difference ... form very few if any mismatch
> ligation products with each other.

So it stays the default, because it is what the modular cloning standards require and what
makes a set portable, and it is the one rule a caller can lower. The data is the better
authority on any particular pair, and that is what `fidelity` is for: the rules pick a set
quickly, and the score says whether that set is any good.

Composition is not a rule beyond `uniform`. Pryor 2020 measured that efficiency is not a
function of composition at all:

> the relative efficiency of each Watson-Crick pair was not simply a function of GC content,
> and thus, difficult to predict based on the sequence composition alone.

So the 25%-GC heuristic that some systems apply is **not** implemented. `uniform` survives
because Potapov 2018 saw a real failure at an all-GC junction, a 23% drop in connections at
junction 6, and that is a truncation the fidelity data does not predict.

Where a matrix exists, free candidates are **ordered** by the strength of their own
Watson-Crick pair, so a design takes the best overhang that passes the rules rather than the
first one in alphabetical order. Ordering, not refusing: a weak pair is never forbidden, which
matters for a scarless junction, where the sequence chooses the overhang and not the designer.

## 7. The fallback for an enzyme nobody measured

PaqCI, AarI, BspQI, BpiI and BtgZI have no published matrix. **A shipped matrix of the same
overhang length stands in**, and the report names whose: `measured` is `True`,
`enzyme_specific` is `False`, and the label reads `measured with Esp3I, not specific to PaqCI`.
Decided in #320 on the measurement below, built in #366.

Which matrix stands in is read off the shipped data, not listed in a table the code carries:

| Unmeasured | Overhang | Stands in | On what ground |
| --- | --- | --- | --- |
| BspQI (GCTCTTC) | 3 | SapI | reads and cuts the same site |
| BpiI (GAAGAC) | 4 | BbsI-HF | reads and cuts the same site |
| PaqCI, AarI (CACCTGC) | 4 | Esp3I | the four-base matrix with the most ligations behind it |
| BtgZI (GCGATG) | 4 | Esp3I | the same |

The first two are barely substitutions: an isoschizomer cuts the same site the same way, so the
measurement is of both enzymes under two names. For the rest the authority is Pryor 2020's own
Discussion sentence, quoted in section 10 — the predicted fidelity "is unlikely to be
significantly impacted by the choice of Type IIS restriction enzyme".

**An earlier draft of this section rejected exactly this**, as a substitution of one enzyme's
measurement for another's, and sent the question to section 8's ligase profile instead. The
objection was never about the number but about what the package may **claim**, and
`enzyme_specific` and `label` already answer that. The measurement settles which substitution is
smaller. Reproduced independently at seed 0, 200 random sets per size:

| Overhangs | Cross-enzyme spread, mean / max | Profile deviation, mean / max |
| --- | --- | --- |
| 8 | 0.0542 / 0.1815 | 0.0763 / 0.2606 |
| 12 | 0.0785 / 0.1766 | 0.1435 / 0.3032 |
| 20 | 0.0762 / 0.1548 | 0.1930 / 0.3819 |

The gap widens with set size — 1.4x at eight overhangs and 2.5x at twenty — because the two
reactions differ systematically and the error compounds over a product of per-junction
probabilities. Section 10 reproduces the same effect from the other direction, at mean 0.052 and
worst 0.148 over the four four-base matrices. On Pryor's own eleven plant overhangs, every
shipped matrix lands within 1.9 points of the 81% the paper reports, including the three that
were not the enzyme used, while every pure-ligation profile reads 8 to 10 points high.

So the stand-in is taken ahead of a ligase profile, and only `prefer_profile` puts a profile
first. One substitution is licensed by an author; the other is licensed by nobody.

The rules remain as the last resort, and their two weights — 0.05 for each near-duplicate
partner and 0.02 for an overhang of one base kind — are **ranking choices and not
measurements**. They are module constants saying so, and a set scored that way must never be
printed beside a measured fidelity as though the two were the same kind of number. In practice
nothing shipped reaches them: every Type IIS enzyme the package holds leaves three or four
bases, and shipped matrices cover both lengths, so the rules now score only an enzyme outside
that range. That is a reachability note, not a reason to delete them.

`choose_enzyme` still prefers an enzyme with its own matrix, which stays right: an enzyme's own
measurement outranks a stand-in. Its wording may no longer say an unmeasured enzyme cannot be
scored.

## 8. A matrix the user holds

Potapov 2018 covers what Pryor does not. It profiled **T4 DNA ligase itself** across all 256
four-base overhangs, at 25 °C and at 37 °C, for 1 hour and for 18 hours, and T7 ligase besides.
Ligation is the ligase's work, so that data speaks to every Type IIS enzyme, including the ones
nobody has measured.

**It is CC BY-NC 4.0, so none of it is here.** No file from that archive is in this repository,
no test fixture is copied from it, and no number of it is transcribed into the code or the tests.
The package reads a matrix from a file **the user already holds**, so nothing is redistributed
and nothing needs licence marking. Whoever downloads the archive accepts its terms themselves.

Not to be confused with EGF's `tatapov_data` package, which redistributes a **CC BY-ND**
repackaging of this data: a different artifact under a different licence, and not what is read
here. **Section 11 flags this sentence** — `tatapov` itself is MIT and downloads its tables rather
than carrying them, and whether a separate CC BY-ND deposit exists was not re-checked.

### Where the archive is

figshare item 7267505, `sb8b00333_si_002.zip`, the Supporting Data of
[doi:10.1021/acssynbio.8b00333](https://doi.org/10.1021/acssynbio.8b00333). It holds fifteen
files. Six of them are count matrices over every overhang pair, one per condition, and those six
are the only ones this package can read:

| File | Ligase | Incubation |
| --- | --- | --- |
| `FileS01_T4_01h_25C.xlsx` | T4 | 1 h at 25 °C |
| `FileS02_T4_01h_37C.xlsx` | T4 | 1 h at 37 °C |
| `FileS03_T4_18h_25C.xlsx` | T4 | 18 h at 25 °C |
| `FileS04_T4_18h_37C.xlsx` | T4 | 18 h at 37 °C |
| `FileS06_T7_18h_25C.csv` | T7 | 18 h at 25 °C |
| `FileS08_T7_18h_37C.csv` | T7 | 18 h at 37 °C |

`FileS03` is the one to reach for, because NEB's own Viewer defaults to it: "The default
conditions are ligation at 25°C for 18 hours; these conditions have been shown to well predict
the results of Golden Gate assembly using typical cycled conditions."

### The other nine files, and why none of them is a seventh matrix

The remaining nine divide into three groups, and **not one of them is a count matrix over
overhang pairs**:

| File | What it is |
| --- | --- |
| `FileS05_HF_cycled.xlsx`, `FileS07_LF_cycled.xlsx`, `FileS09_DP_cycled.xlsx`, `FileS10_FP_cycled.xlsx` | the four ten-fragment test assemblies, cycled 5 min 37 °C / 5 min 16 °C, 30 times |
| `FileS11_HF_01h_37C.xlsx`, `FileS12_LF_18h_37C.xlsx`, `FileS13_DP_18h_37C.xlsx`, `FileS14_FP_18h_37C.xlsx` | the same four assemblies, held at 37 °C instead of cycled |
| `lac.fasta` | one 4,851-base record, the lac cassette of the twelve- and twenty-four-fragment assemblies |

**`HF`, `LF`, `DP` and `FP` are not ligases and not buffers. They name the four junction sets**
the paper designed for its ten-fragment assembly of inserts A to J: high-fidelity,
low-fidelity, deletion-prone and failure-prone. The preprint of the same work states it:

> The junctions between fragment pairs (Junctions 1 – 9) were selected to either be 9
> high-fidelity (HF) junctions, or a low-fidelity (LF) set where 9 junction pairs were chosen
> such that many mismatch ligation events were predicted (Table 1).

And, on the same page:

> Two additional sets were designed: a deletion-prone (DP) set, where junction 7 of the HF set
> was changed [...] such that deletion (and to a lesser extent, duplication) of insert G was
> predicted to result; a failure-prone (FP) set where junction 7 was replaced with the high
> fidelity but low efficiency pair.

So `FileS11` to `FileS14` are **the same kind of file as the `*_cycled` ones**, differing only in
the reaction they summarise: `_cycled` was thermocycled, the others were held at 37 °C. Opening
them confirms it. `FileS03` has one sheet, `18h @ 25C`, carrying 256 overhang labels; each of
these eight has five, `table_01` to `table_05`, labelled `Insert`, `Count`, `Fraction` and the
inserts `A` to `J`. `read_profile` refuses all eight with one message, naming the label `B` as
one that is not an overhang.

**Verdict: the four are excluded, for the reason the `*_cycled` files are excluded** — they
count assemblies, not ligation events between overhang pairs, so there is nothing in them to
score a set of overhangs against. `FileS05`, `FileS07`, `FileS09` and `FileS10` are all the
`*_cycled` files there are. `lac.fasta` is a sequence file nothing here reads, and needs no
entry above. The archive holds six matrices, not ten.

### Pointing the tool at a copy

A path under `~` is not one this project may depend on, so the example reads the copy under
`reference_docs/`, which `reference_docs/ligation-fidelity/README.md` says how to fetch. That
directory is git-ignored: the download is the reader's own, and **nothing of it ships**.

```python
from mbio.cloning.goldengate import plan_assembly
from mbio.ligase import read_profile

profile = read_profile("reference_docs/ligation-fidelity/potapov2018/FileS03_T4_18h_25C.xlsx")
plan_assembly(vector, insert, enzyme="PaqCI", profile=profile)
```

On the command line, as an option or as an environment variable:

```sh
matrix=reference_docs/ligation-fidelity/potapov2018/FileS03_T4_18h_25C.xlsx

mbio cloning goldengate plan vector.dna insert.dna --out run --ligase-matrix "$matrix"

export LIULAB_MBIO_LIGASE_MATRIX="$matrix"
mbio cloning goldengate plan vector.dna insert.dna --out run
```

Both shapes load with the standard library alone. An `.xlsx` is read by the same `zipfile` and
`xml.etree` code the build script uses, which moved into `mbio.ligase` so that
the package and the script share one reader; a `.csv` is read by `csv`. **No dependency was
added.** A file that is not a count matrix is refused with a message saying what one is, rather
than a stack trace: a header row of overhang labels, the same labels down the first column, and
a count in each cell.

The conditions are read off the file name where it is named the way the archive names one, and
`conditions=` states them for a file that has been renamed.

### What the report then says

A profile belongs to the ligase and the conditions, not to the Type IIS enzyme, and
`FidelityReport` keeps the four kinds of number apart:

| Scored by | `measured` | `enzyme_specific` | `label` |
| --- | --- | --- | --- |
| the enzyme's own shipped matrix | `True` | `True` | `measured` |
| another enzyme's, standing in | `True` | `False` | `measured with Esp3I, not specific to PaqCI` |
| a ligase profile | `True` | `False` | `measured ligase profile, not specific to PaqCI` |
| the rules | `False` | `True` | `rule-based estimate` |

`FidelityReport.stand_in` is what tells the middle two apart: it carries the product whose
matrix stood in, and is empty for a profile.

`source` names the conditions and the file, so the protocol cites the file the number came from
and the overview prints the label beside the percentage. **A ligase profile is never presented as
a measurement of the enzyme**, and a rule-based score is still never printed beside a measured
one as though the two were the same kind of number.

The enzyme's own matrix wins where there is one, because a profile stands in for a measurement
nobody has made rather than replacing one somebody has. **A profile is now the third fallback,
not the second**: section 7's stand-in comes first, since it measures the same reaction, and a
profile is reached only for an enzyme no shipped matrix shares an overhang length with.
`prefer_profile` overrides all of that, for a designer who wants one set of conditions across
several enzymes, and it is how a caller who passed `--ligase-matrix` on purpose keeps it.

**Without such a file, PaqCI is scored on Esp3I's matrix** rather than by the rules, which a
test pins.

### A second ligase's matrix, and what it may not replace

Bilotti 2022 is read the same way, and it is CC BY 4.0 rather than CC BY-NC, so a user holding it
carries no licence question at all. It is one workbook of eight sheets (section 11), so the sheet
is named:

```python
profile = read_profile(
    "reference_docs/ligation-fidelity/bilotti2022/File S1_NAR.xlsx", sheet="File S2. T7"
)
```

A workbook of several sheets, asked for none, is refused and names them, so a question about T7
cannot be answered with T4's numbers. The sheet is resolved through the workbook's own
relationships rather than by the number in `sheetN.xml`, which need not follow the order the
sheets are listed in. Its name also stands in for the conditions where the file name states none,
as `File S1_NAR.xlsx` does: `"File S7. T7 PEG"` already carries the ligase and the buffer, so
`LigaseProfile` gained no buffer field of its own.

**A T7 profile may not replace an enzyme's own one-pot matrix.** The temptation is real: an iGGA
round ligates with T7 in PEG, so a T7 PEG matrix looks like the closer match to the bench. It is
not. Scored on this branch against the shipped BbsI-HF matrix, 200 random sets per size, seed 1, a
Bilotti T7 PEG profile reads **+0.080 at six overhangs, +0.135 at eight and +0.249 at twelve** —
higher, which is to say more permissive, and growing with set size exactly as section 10's T4
profile does. Pure ligation without the Type IIS enzyme and without cycling reads optimistic
against one-pot numbers whichever ligase it measures, so swapping in the matching ligase trades a
labelling problem for an accuracy problem in the loose direction. The enzyme's own one-pot matrix
stays the score; a second ligase's matrix is read for what it compares, not for what it scores.

### What it is read for instead: one overhang's on-target rate

`mbio.overhangs.on_target` is what such a profile is read for. It reports, per overhang,
the correct Watson-Crick pair per 100,000 events on the profile's own sheet, and names the ones
below `STRONG_LIGATION`. It returns no score, nothing ranks on it, and the fidelity number above
does not move.

**The comparison is only honest at matched buffer.** The paragraphs above this sub-section reason
from T7 against T4 in standard T4 buffer, and that is the wrong pair of cells for a reaction run
in PEG. PEG raises both ligases, so a T7 PEG number read against T4 in standard buffer flatters
T7 and a T7 standard-buffer number read against T4 PEG damns it. The sheets to compare are
`File S6. T4 PEG` against `File S7. T7 PEG`, and `File S1. T4` against `File S2. T7`. Measured on
a user-held copy of `File S1_NAR.xlsx`, correct Watson-Crick pair per 100,000 events:

| Overhang | T4 | T7 | T4 PEG | T7 PEG | T7 PEG / T4 PEG |
| --- | --- | --- | --- | --- | --- |
| `AGGA` | 325.3 | 95.5 | 302.3 | 275.7 | 0.91 |
| `AGAT` | 336.7 | 66.5 | 284.6 | 175.8 | 0.62 |
| `GCAT` | 343.3 | 497.2 | 282.4 | 458.9 | 1.62 |
| `TTCC` | 400.0 | 165.0 | 358.9 | 335.7 | 0.94 |

Those four are the AP-1 example's three entry overhangs and its cloning scar. At matched buffer
T7 costs `AGAT` 38% of its rate and leaves `GCAT` better than T4 — which is why the check is an
absolute floor and not a ratio. A ratio rewards `GCAT` for no bench reason and condemns `AGAT`,
which at 175.8 is in no trouble at all.

**The A/T-rich overhangs usually cited are unreachable, so they do not justify the floor.** The
seven-fold losses quoted for `TTAA`, `TATA`, `TAAA`, `TTTA` and `AAAA` describe overhangs
`mbio.overhangs.refusal` already refuses — the first two as palindromes, the last three as
uniform. Of the 256 four-base overhangs, 216 are reachable and 40 are not (16 palindrome, 24
uniform). The whole reachable tail below the floor is nine strand pairs:

| Overhang | T4 PEG | T7 PEG | ratio |
| --- | --- | --- | --- |
| `TAGA` / `TCTA` | 173.8 | 60.8 | 0.35 |
| `TCAA` / `TTGA` | 165.2 | 63.2 | 0.38 |
| `TGAA` / `TTCA` | 212.6 | 69.3 | 0.33 |
| `CTAA` / `TTAG` | 204.2 | 71.8 | 0.35 |
| `CTTA` / `TAAG` | 238.9 | 74.2 | 0.31 |
| `ACTA` / `TAGT` | 205.9 | 85.7 | 0.42 |
| `AGAA` / `TTCT` | 276.4 | 86.1 | 0.31 |
| `AAGA` / `TCTT` | 256.8 | 86.5 | 0.34 |
| `TACA` / `TGTA` | 275.2 | 99.9 | 0.36 |

About threefold, not sevenfold, and every one of the nine is three A/T bases plus one G or C.
**The floor is `STRONG_LIGATION`**, NEB's own threshold for a strong Watson-Crick pair, already
sourced in section 5 and reused rather than added to. Over the 216 reachable overhangs on
`File S7. T7 PEG` it warns on 18 — those nine pairs, counted on both strands — leaves AP-1's
worst at 175.8 silent, and a floor of 50 would warn on nothing at all.

**What a low rate costs is colonies, not product.** Bilotti's counts are of the *correct* pair,
so a threefold loss is threefold fewer good joins. Each iGGA round's tube holds one entry overhang
and the cloning scar and nothing else — `synbio.igga.gate.LIGATION_OVERHANGS` — so its
fidelity saturates and a set comparison has nothing to discriminate at two. Yield is what is at
risk, and section 11 records
Strzelecki 2024 attributing it to duplex strength with no count matrix capturing it.

**Nothing is designed on it.** `synbio.igga.standard` keeps ranking candidates on the
enzyme's own shipped matrix. The profile is a file the user holds, so ranking on it would make the
same build yield different overhangs depending on whether that file is present, and a design
that is not reproducible from the build alone costs more than the overhangs it would save.

## 9. Regeneration

```sh
pixi run python scripts/build_ligation_fidelity.py            # downloads the five tables
pixi run python scripts/build_ligation_fidelity.py --from DIR # from files already downloaded
```

`journals.plos.org` serves each table from
`https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0238592.sNNN&type=supplementary`,
redirecting to signed storage, which `urllib` follows. Nothing else is fetched.

Tests never touch the network: `tests/scripts/test_build_ligation_fidelity.py` writes the three XML
parts of a workbook itself and reads them back, which is also what makes the guards testable.

### Which supplementary file each matrix came from

`source.url` in the data file is that template, `sNNN` and all, because one URL cannot name five
files. **The resolved id is not missing**: each matrix carries its own `file`, and
`scripts/build_ligation_fidelity.py` holds the suffix in its `TABLES` tuple. Checked on
2026-10-06 by downloading all five and reading `xl/workbook.xml`:

| Enzyme | `matrices[].file` | Suffix | The one sheet in that workbook | sha256, first 16 |
| --- | --- | --- | --- | --- |
| BsaI | `pone.0238592.s001.xlsx` | `s001` | `S1 Table. BsaI-HFv2` | `320dd058f3ca6372` |
| BsmBI | `pone.0238592.s002.xlsx` | `s002` | `Table S2. BsmBI-v2` | `7c444e99e5e4d245` |
| Esp3I | `pone.0238592.s003.xlsx` | `s003` | `Table S3. Esp3I` | `1557e62cfa4f89cd` |
| BbsI | `pone.0238592.s004.xlsx` | `s004` | `Table S4. BbsI-HF` | `56143bb445e6d84b` |
| SapI | `pone.0238592.s005.xlsx` | `s005` | `Table S5. SapI` | `dfa2df906bd306e3` |

All five addresses answered 302 to storage and then 200 with the spreadsheet content type, so a
reviewer can repeat this. **The mapping is verified rather than asserted**, and by the build
itself: section 3's first guard refuses a workbook whose sheet name does not carry the product
the table is supposed to measure, so a shipped matrix under the wrong enzyme could not have been
built. The sheet names above are what that guard matched.

The copies are under `reference_docs/ligation-fidelity/pryor2020/`, which is git-ignored, so a
fresh clone has none of them; the table above is the record. They were downloaded only to resolve
the mapping. Nothing was copied out of them — the shipped JSON is still the build script's own
output from those same addresses.

## 10. Which dataset is authoritative, and what the two numbers each mean

Read 2026-10-06. A prototype scored one design at **1.0000** on the shipped Pryor BsaI matrix and
**0.9978** on the user-held Potapov `FileS03_T4_18h_25C.xlsx`, and called the shipped matrix
"more generous and far easier to satisfy" — the implication being that it is the weaker evidence.
**That characterisation is wrong, in both directions**, and the papers and the data each say so.

### The two datasets measure different reactions

Pryor 2020, Materials and Methods, "Golden Gate assembly fidelity and bias assay". The Type IIS
enzyme is **in the tube** and the reaction is **thermocycled**:

> Reactions (20 μL final volume) with T4 DNA ligase and BsaIHF-v2 or BsmBI-v2 were carried out
> using their respective NEB Golden Gate Enzyme Mixes (2 μL) in 1X T4 DNA ligase buffer. [...]
> The reactions were cycled between 37°C and 16°C (SapI, Esp3I, BsaI-HFv2, BbsI-HF) or 42°C and
> 16°C (BsmBI-v2) for 5 minutes at each temperature for 30 cycles, and then subjected to a final
> heat-soak for 5 minutes at 60°C.

Its S1-S5 Table legends say the same in one line: "Ligation frequency for each overhang pair in
assembly reactions with BsaI-HFv2 and T4 DNA ligase", and so on for the other four.

Potapov 2018, Methods, ligation reaction. **No restriction enzyme, and a static hold**:

> In a typical ligation reaction, substrate (100 nM) was combined with 2.5 μL high concentration
> T4 DNA ligase (2000 U, 1.75 μM final concentration) in 1× T4 DNA ligase buffer in a 50 μL total
> reaction volume and incubated for 1 h or 18 h at 25°C or 37°C.

The site is cut during substrate preparation and the cut substrate purified before ligation
("confirmed to be >95% cut"), so the ligase meets ready-made ends. Cycling appears in that paper
only in its separate validation assemblies, which are real Golden Gate reactions.

The readout is the same PacBio SMRT hairpin assay in both — Pryor's substrates were "prepared as
previously described" from Potapov — so the two are comparable numbers from different chemistry.
**So the first half of the hypothesis holds: these are different reactions, and the
enzyme-specific set is the one that matches what this package actually plans.**

### The re-cutting mechanism is not in either paper — do not claim it

A tempting explanation for the gap is that in a one-pot reaction the Type IIS enzyme re-cuts a
mis-ligated junction, which is then retried, so one-pot fidelity should read higher. **Neither
paper says this.** Pryor 2020's full text was searched for `re-cut`, `recut`, `re-cleav`,
`reversib` and `proofread`: no hits. The only adjacent sentence is about yield, not fidelity:

> Suboptimal conditions, such as temperature or buffer conditions where the activity of the
> restriction enzyme is poor and cutting is inefficient relative to re-ligation, could decrease
> the assembly yield.

The caveat the papers do raise runs the other way — cutting and overhang melting are steps the
ligation-only assay does not capture, so an all-G/C overhang may under-assemble relative to what
ligation data predicts. That is the same observation `uniform` already rests on (section 6).

### And the direction is the opposite of "more generous"

Pryor 2020, Results, comparing its own one-pot data with Potapov's pure ligation:

> However, in comparison with our previous ligation fidelity study, we note **higher frequencies
> of mismatch pairs** and less bias against A/T-rich overhang sequences under Golden Gate assembly
> conditions. Presumably, this is due to differences in the reaction temperatures and buffer
> conditions between the two studies.

More mismatching means **lower** computed fidelity. The one-pot enzyme-specific data is the
stricter of the two, not the looser one. Measured on this branch, 2026-10-06, the shipped matrices
reproduce the paper's own published Golden Gate fidelities and the pure-ligation profile does not:

| Pryor 2020's own example, BsmBI-v2 | Paper | Shipped matrix | Potapov T4 18 h / 25 °C |
| --- | --- | --- | --- |
| 11 plant overhangs (Fig 4A) | 81% | 80.93% | 89.65% |
| the same set without the `GGTA` mispair (Fig 4A) | 92% | 92.11% | 96.80% |
| the same set extended to 20 by GetSet (Fig 4B) | 80% | 80.48% | 87.19% |

The profile reads 7 to 9 points **high** on all three. Over 280 random BsaI sets, seed 0, the
shipped matrix is the lower score in 220 of them, and the gap grows with the set:

| Overhangs in the set | Mean shipped minus profile | Shipped higher |
| --- | --- | --- |
| 2 | -0.0008 | 27/40 |
| 4 | -0.0184 | 16/40 |
| 6 | -0.0452 | 7/40 |
| 8 | -0.0613 | 7/40 |
| 12 | -0.1395 | 2/40 |
| 16 | -0.1633 | 1/40 |
| 20 | -0.1898 | 0/40 |

At two overhangs the two agree to a thousandth and which one is higher is a coin flip. **That is
the regime the prototype measured**: a 0.0022 difference on one small set, generalised into a
claim that reverses by six overhangs and is wrong in all forty sets at twenty.

**The one sense in which the shipped data is more forgiving** is the weak-pair floor, and it is
narrow. The lowest normalised Watson-Crick pair is 118.0 (`TTAA`) in shipped BsaI-HFv2 against
18.5 (`TTAA`) in the T4 profile, so `FidelityReport.weak` can fire on a profile and never fires on
shipped data — which is what section 5 records. That is a different statement from set fidelity,
and it is the only one the prototype's wording fits.

### What this says about the fallback

Pryor 2020's Discussion licenses carrying a matrix across Type IIS enzymes:

> Thus, the predicted fidelity of overhang sets is unlikely to be significantly impacted by the
> choice of Type IIS restriction enzyme, and this is likely broadly applicable to all Type IIS
> restriction enzymes that generate the same overhang structure including enzymes not explicitly
> tested here.

Measured here, substituting one enzyme's matrix for another's is a smaller error than substituting
the pure-ligation profile. The plant set scores 81.78% on BsaI, 80.93% on BsmBI, 81.24% on Esp3I
and 79.90% on BbsI; across 200 random eight-overhang sets the spread between those four matrices
is mean 0.052, median 0.047, worst 0.148. Against that, the profile was 7 to 9 points off on the
paper's own benchmark and up to 0.19 off at twenty overhangs.

**So section 7's rejected idea is the better-conditioned one of the two fallbacks, on this
measurement and on the paper's own sentence.** It was rejected as a substitution of one enzyme's
measurement for another's, and a substitution of a different reaction's measurement was adopted
instead — the larger error, and the one the authors do not license.

**This was decided in #320 and built in #366, and section 7 now records it**: a shipped matrix of
the same overhang length stands in ahead of a pure-ligation profile, `enzyme_specific` stays
`False`, and the label says whose measurement it is.

**PaqCI is no longer an enzyme nobody measured** — it is one whose measurement is unpublished and
unshippable (section 10). That strengthens the recommendation rather than weakening it: the gap
being filled is in what is *available*, not in what is *known*, so there is no published number
the package is approximating badly and no prospect of one arriving. It also means
`choose_enzyme`'s premise holds for a reason worth stating. Preferring an enzyme that has its own
matrix is not a claim that PaqCI ligates less predictably; it is a claim that this package can
show its work for BsaI and cannot for PaqCI.

Pryor 2020 also sets the ceiling on all of this:

> Importantly, predicted assembly fidelity should be taken as a qualitative prediction, most
> useful for comparing expected performance between alternative junction sets.

### What NEB's own tool defaults to, and what it serves

**Read in a real browser on 2026-10-06.** The tools are `.cgi`, not `.html`, which is why fetching
`run.html` failed: the Viewer is `https://ligasefidelity.neb.com/viewset/run.cgi`, and the live
site is **version 1.0**. `help.html` loads normally from a browser — **the 403 is bot-blocking,
not a missing page**, and `nebridgetools.neb.com` genuinely does not exist.

It defaults to Potapov 2018's pure-ligation profile, justified by cross-validation rather than by
matching the reaction:

> The default conditions are ligation at 25°C for 18 hours; these conditions have been shown to
> well predict the results of Golden Gate assembly using typical cycled conditions (16°C 5
> min/37°C 5 min, 30 cycles).

The help page also states the axis convention independently of the sentence section 4 quotes:
rows are the top strand, columns the bottom strand, both written 5' to 3'. **That corroborates
section 4** from a second place on the same site.

An earlier draft of this section said the tool has no Type IIS enzyme selector. **That was wrong**
— it was read off the help page rather than the tool. The `dataset` select offers **fifteen**
four-base conditions, verbatim from its option labels:

| # | Condition |
| --- | --- |
| 1-4 | T4 DNA Ligase at 25 °C / 1 h, 25 °C / 18 h, 37 °C / 1 h, 37 °C / 18 h |
| 5-6 | T7 DNA Ligase at 25 °C / 18 h, 37 °C / 18 h |
| 7 | `BsaI-HFv2 37-16 cycling` |
| 8 | `BsaI-HFv2 37 static` |
| 9 | `BsmBI-v2 42-16 cycling` |
| 10 | `Esp3I, 1x T4 DNA Ligase buffer, 37-16 cycling` |
| 11 | `BsaI-HFv2, 1x NEBridge Ligase MM, 37-16 cycling` |
| 12 | `PaqCI, 1x T4 DNA Ligase buffer, 37-16 cycling` |
| 13 | `PaqCI, 1x NEBridge Ligase MM, 37-16 cycling` |
| 14 | `BsmBI-v2, 1x NEBridge Ligase MM, 42-16 cycling` |
| 15 | `BbsI-HF, 1x NEBridge Ligase MM, 37-16 cycling` |

The `ohlen` select also offers three-base overhangs. **That list was not enumerated**, so what
three-base conditions exist is unchecked here rather than guessed.

### A PaqCI matrix exists, and it is unpublished vendor data

Rows 12 and 13 carry the internal identifiers `NBC_20250411_PaqCI_60_Cycles_T4Buffer` and
`NBC_20250411_PaqCI_60_Cycles_LigaseMM`. **So a PaqCI matrix exists and NEB holds it**, dated
2025-04-11 by its own identifier. The search-index lead section 11 recorded is now a fact.

The provenance is the decisive part. The help page cites exactly three references: both Potapov
2018 papers, [doi:10.1093/nar/gky303](https://doi.org/10.1093/nar/gky303) and
[doi:10.1021/acssynbio.8b00333](https://doi.org/10.1021/acssynbio.8b00333), and a Current
Protocols protocol, [doi:10.1002/cpz1.882](https://doi.org/10.1002/cpz1.882). **Pryor 2020 is not
cited there at all**, yet rows 7 to 10 are Pryor's enzyme-cycling sets. Rows 11 to 15 — the two
PaqCI sets and the three NEBridge Ligase Master Mix sets — **are cited to nothing**, and carry only
those identifiers.

So the help page's own description of its data is now **stale**:

> The data used comes from recent publications on ligation fidelity using T4 DNA Ligase (1, 2) and
> **represents subsets** of the data sets discussed in those papers.

Six of the fifteen conditions are in neither cited paper. That bears on trusting the tool as a
reference: it is a vendor tool serving a mix of published and unpublished data without saying
which is which.

**NEB data is all rights reserved, and section 2's standing position already covers this: cited,
never mirrored.** Nothing was downloaded from the tool, and nothing here implies the package would
take any of it. A PaqCI matrix being measured does not make it available.

The buffer axis is worth noting for its own sake. The same enzyme appears twice, under plain T4
ligase buffer and under NEBridge Ligase Master Mix — **exactly the axis Bilotti 2022 varies**
(section 11). Two independent sources now treat buffer as a condition that changes the matrix,
which is support for giving `LigaseProfile`'s conditions a buffer field.

## 11. Other fidelity data sources, and whether each may ship

Surveyed 2026-10-06. The repo holds Pryor 2020 (shipped) and reads Potapov 2018 (user-held).

| Source | What it adds | Licence | Obtainable | May ship |
| --- | --- | --- | --- | --- |
| **Bilotti et al. 2022**, *NAR* 50(8):4647-4658, [doi:10.1093/nar/gkac241](https://doi.org/10.1093/nar/gkac241) | **Five ligases** — T4, T3, T7, PBCV-1 (SplintR), human Ligase 3 — over all 256 four-base overhangs, by the same SMRT assay. PEG for three of the five | **CC BY 4.0**, confirmed on the publisher's own deposited permissions block and on Crossref (Sources) | Yes, and downloaded: `gkac241_supplemental_files.zip` at PMC, holding **one workbook of eight sheets** | **Yes** |
| Pryor et al. 2022, *ACS Synth Biol* 11(6):2036-2042, [doi:10.1021/acssynbio.1c00525](https://doi.org/10.1021/acssynbio.1c00525) | No new matrices — applies Pryor 2020 to a 40 kb, 52-part build | CC BY-NC-ND 4.0 | Supplement only | No |
| Strzelecki et al. 2024, *NAR* 52(19):e95, [doi:10.1093/nar/gkae809](https://doi.org/10.1093/nar/gkae809) | The one independent non-NEB re-measurement. Gel kinetics on 6 overhangs, BsaI-HFv2 + T4. Finds **overhang duplex strength**, not only mismatch fidelity, drives efficiency — a factor no matrix here captures | CC BY-NC (the bioRxiv preprint is no-reuse) | Paper yes, no matrix deposit | No |
| Mukundan & Madhusudhan 2025, OOGGA, [doi:10.1101/2025.06.16.659877](https://doi.org/10.1101/2025.06.16.659877) | No new measurement; scores against Potapov 2018 | Unstated | Code on GitHub | n/a |
| NEBridge Ligase Fidelity Viewer v1.0, `ligasefidelity.neb.com/viewset/run.cgi` | Fifteen four-base conditions, of which **six are in no paper** — two **PaqCI** sets and three NEBridge Ligase Master Mix sets (section 10) | Vendor tool, all rights reserved | **No** — the matrices are not downloadable | No |
| `tatapov` (Edinburgh Genome Foundry) | Nothing new. Its code is MIT and it **downloads the tables at run time** rather than vendoring them; its upstream is exactly Potapov 2018 and Pryor 2020 | MIT (code only) | Yes | The data's own licence still governs |
| GoldenHinges, DNA Chisel, kappagate | No independent dataset; annealing data via `tatapov` | MIT (code) | Yes | n/a |
| Duckworth 2023, Sikkema 2023, Lund 2024 (Springer methods chapters) | Protocols for measuring or applying this data. No dataset | All rights reserved | — | No |
| NEB patents US 12,188,011 and US 12,435,332 | PaqCI plus activator oligo beating AarI on assembly performance, and the matrix *shapes*. **No overhang-pair table found** | Patent text; claims enforceable | Text yes | No |

Four things worth carrying forward.

**Bilotti 2022 is the find, and it closes a gap this note listed as open.** Section 10's old table
said T7 ligase was covered "only [by] Potapov 2018, CC BY-NC" and therefore could not ship. That
is **no longer true**: T7, T3, PBCV-1/SplintR and human Ligase 3 are all in Bilotti 2022 under
CC BY 4.0, over all 256 four-base overhangs. It extends coverage along the axis Pryor does not —
Pryor added Type IIS enzymes under one ligase, Bilotti adds ligases.

The deposit was downloaded and opened on 2026-10-07, and it is not the shape the paper describes.
Its Data Availability says "Raw ligation product observation counts were provided as CSV formatted
data tables"; what is deposited is **one workbook, `File S1_NAR.xlsx`, of eight 256 x 256 sheets**.
That distinction was load-bearing here, because `read_profile` read the first sheet of a workbook
and nothing else, so a request for T7 returned T4 without erring. It now names the sheet, and
refuses a workbook of several sheets that names none; section 8 gives the selector.

| Sheet | Ligase | Buffer | Observations |
| --- | --- | --- | --- |
| `File S1. T4` | T4 | 1x T4 DNA ligase buffer | 317,228 |
| `File S2. T7` | T7 | 1x T4 DNA ligase buffer | 338,272 |
| `File S3. hLig3` | human Ligase 3 | 1x T4 DNA ligase buffer | 643,492 |
| `File S4. T3` | T3 | 1x T4 DNA ligase buffer | 344,420 |
| `File S5. PBCV-1` | PBCV-1 (SplintR) | 1x T4 DNA ligase buffer | 227,846 |
| `File S6. T4 PEG` | T4 | NEBNext Quick Ligation, 6% PEG 6000 | 418,184 |
| `File S7. T7 PEG` | T7 | NEBNext Quick Ligation, 6% PEG 6000 | 245,162 |
| `File S8. hLig3 PEG` | human Ligase 3 | NEBNext Quick Ligation, 6% PEG 6000 | 394,032 |

The second condition axis is the buffer: standard T4 buffer against NEBNext Quick Ligation buffer,
which has PEG, and PEG changes bias. It covers **three of the five ligases, not all five** — T4, T7
and hLig3 have a PEG sheet, T3 and PBCV-1 do not. A matrix from here needs its **buffer** recorded
next to temperature and time; the sheet name carries both, so it supplies the free-form conditions
string where the file name states none, and no buffer field was added for a value nothing computes
on.

All eight sheets are **pure ligation: 1 h at 25 °C, no Type IIS enzyme and no cycling** — the same
chemistry as Potapov 2018 and not the chemistry of a Golden Gate reaction. Section 10's measurement
therefore applies to every one of them, whichever ligase it names: a profile from this source reads
high against a one-pot matrix and must never be taken for a Golden Gate measurement. Section 8 says
what follows for the T7 sheets in particular.

**Taq ligase and E. coli ligase have no such data and probably cannot.** Taq ligase is
nick-selective rather than end-joining, so an end-joining overhang matrix for it is not a
meaningful object. Nothing post-2020 covers E. coli ligase at overhang level.

**No publication carries a PaqCI matrix, and one exists anyway.** PaqCI is absent from Potapov
2018, Pryor 2020, Pryor 2022 and Bilotti 2022; the only PaqCI evidence in the literature is the two
NEB patents, which give aggregate assembly performance against AarI and describe matrix shapes but
carry no pair table. AarI is in the same position, with no matrix anywhere.

But NEB's live Viewer serves **two PaqCI datasets**, under internal identifiers dated 2025-04-11
and cited to no publication at all (section 10). So the old wording "nobody has published one" is
true and materially incomplete: **the measurement has been made and is not public.** It is all
rights reserved, not downloadable from the tool, and so citable but never shippable — which is
where section 2 already puts anything of NEB's.

**The open tooling ecosystem rests entirely on the two datasets this repo already knows.** Every
Golden Gate overhang designer checked — `tatapov`, GoldenHinges, DNA Chisel, OOGGA — scores against
Potapov 2018 or Pryor 2020 and nothing else. No vendor other than NEB publishes overhang-pair
ligation data: nothing from IDT, Twist, Thermo, Promega or Takara.

One correction to section 8 falls out of this. That section calls `tatapov` a package that
"redistributes a **CC BY-ND** repackaging of this data". `tatapov`'s own README says the opposite
of the redistribution half — "Tatapov provides these tables (it will download them automatically
[...])" — and its licence file is MIT. **Whether a separate CC BY-ND deposit exists was not
re-checked**, so the sentence is flagged rather than rewritten.

## 12. Open gaps

| Item | Why it is missing | What is done instead |
| --- | --- | --- |
| The exact query set NEB's own tool uses | The help page states the axis convention and the default, but not the arithmetic | The reading in section 5, checked against three published numbers |
| Why two GetSet table sets score 1.5 and 4 points low | The sets as printed may not be the whole reaction | Recorded in section 5, untouched |
| A **shippable** matrix for PaqCI, and any matrix for AarI, BspQI, BtgZI | No publication carries one, re-verified 2026-10-06 against Potapov 2018, Pryor 2020, Pryor 2022 and Bilotti 2022. For PaqCI the measurement exists but is NEB's, unpublished and all rights reserved | A ligase profile where the user holds one, the rule-based fallback otherwise, each labelled. Section 10 recommends preferring another four-base Type IIS matrix over the pure-ligation profile |
| Which three-base conditions NEB's Viewer offers | The `ohlen` select was not enumerated | Nothing. Unchecked rather than guessed, and the four-base list is in section 10 |
| Whether anyone will publish the PaqCI data | Not knowable from here | Nothing. Worth re-checking if NEB publishes a successor to Pryor 2020 |

Closed since 2026-09-12:

| Item | What closed it |
| --- | --- |
| T7 DNA ligase, and ligases other than T4 | **Bilotti 2022 covers T4, T3, T7, PBCV-1/SplintR and human Ligase 3 under CC BY 4.0**, so this no longer depends on a CC BY-NC copy. Section 11 |
| Which supplementary file each shipped matrix came from | Resolved and verified by sheet name and checksum. Section 9 |
| Whether the shipped matrix is the weaker evidence | It is not, and it is the stricter of the two on every set above four overhangs. Section 10 |
| Whether NEB's tool serves a PaqCI matrix | **It serves two**, under internal identifiers citing no publication. Read in a browser; section 10 |
| The Viewer help page returning 403 | **Bot-blocking, not a missing page.** It loads in a browser, and the tools are `.cgi` rather than `.html`. Section 10 |

Still unmeasured by anything here: **overhang duplex strength** as a driver of assembly efficiency,
which Strzelecki 2024 reports and which no count matrix captures.

## Sources

Sections 1 to 9 were read on 2026-09-12. Sections 10 to 12 were read on **2026-10-06**, and each
source below says which.

Added 2026-10-06:

- Bilotti, K., Potapov, V., Pryor, J.M., Duckworth, A.T., Keck, J.L. and Lohman, G.J.S. (2022)
  Mismatch discrimination and sequence bias during end-joining by DNA ligases.
  *Nucleic Acids Research* 50(8), 4647-4658.
  [doi:10.1093/nar/gkac241](https://doi.org/10.1093/nar/gkac241). **CC BY 4.0, confirmed
  2026-10-07 on the publisher's own deposit**, not only on Europe PMC's `license: cc by` for
  PMC9071435. The permissions block in OUP's deposited JATS, served at
  `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9071435/fullTextXML`:

  > © The Author(s) 2022. Published by Oxford University Press on behalf of Nucleic Acids
  > Research. This is an Open Access article distributed under the terms of the Creative Commons
  > Attribution License (<https://creativecommons.org/licenses/by/4.0/>), which permits
  > unrestricted reuse, distribution, and reproduction in any medium, provided the original work
  > is properly cited.

  Crossref's publisher-deposited record for the same DOI carries
  `license[0].URL = https://creativecommons.org/licenses/by/4.0/`, with `content-version: vor` and
  `delay-in-days: 0`, so the version of record is CC BY from the day of publication. Attribution is
  the only condition — the same footing as Pryor 2020, which already ships. The Supplementary Data
  is **one workbook, `File S1_NAR.xlsx`, of eight sheets**, not the CSV tables the paper's own Data
  Availability describes; a copy is under `reference_docs/ligation-fidelity/bilotti2022/`, which is
  git-ignored
- Pryor, J.M., Potapov, V., Bilotti, K., Pokhrel, N. and Lohman, G.J.S. (2022) Rapid 40 kb genome
  construction from 52 parts through data-optimized assembly design. *ACS Synth. Biol.* 11(6),
  2036-2042. [doi:10.1021/acssynbio.1c00525](https://doi.org/10.1021/acssynbio.1c00525)
  (CC BY-NC-ND 4.0). No new matrices
- Strzelecki, P., Joly, N., Hébraud, P., Hoffmann, E., Cech, G.M., Kloska, A., Busi, F. and
  Grange, W. (2024) Enhanced Golden Gate Assembly: evaluating overhang strength for improved
  ligation efficiency. *Nucleic Acids Research* 52(19), e95.
  [doi:10.1093/nar/gkae809](https://doi.org/10.1093/nar/gkae809) (CC BY-NC; the bioRxiv preprint,
  [doi:10.1101/2022.09.09.507109](https://doi.org/10.1101/2022.09.09.507109), is no-reuse)
- Pryor et al. 2020's full text as JATS XML, read for the Methods and Results quoted in section 10:
  `https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0238592&type=manuscript`
- The five S1-S5 supplementary workbooks themselves, read for the mapping in section 9. Copies
  under `reference_docs/ligation-fidelity/pryor2020/`, which is git-ignored
- `tatapov`, Edinburgh Genome Foundry: its `README.rst` and `LICENSE` (MIT), read for whether it
  vendors or downloads its tables
- NEB patents US 12,188,011 and US 12,435,332, read for PaqCI; no overhang-pair table found
- NEB, *Ligase Fidelity Viewer: Help Page*,
  [tools.neb.com](https://tools.neb.com/~potapov/ligase-fidelity-viewer/help.html) — **loads now**,
  by plain request; the 403 recorded on 2026-09-12 is stale
- NEB, **NEBridge Ligase Fidelity Viewer v1.0**,
  `https://ligasefidelity.neb.com/viewset/run.cgi`, and its help page — **read in a browser**,
  which is the only way past the bot-blocking. Source of the fifteen four-base conditions, the two
  PaqCI identifiers and the axis convention in section 10. All rights reserved; nothing downloaded
- Current Protocols, [doi:10.1002/cpz1.882](https://doi.org/10.1002/cpz1.882) — the third and only
  other reference the Viewer's help page cites. Noted for completeness; not read

Could not be read on 2026-10-06, and why:

| What | Status |
| --- | --- |
| Potapov et al. 2018 as published in *ACS Synth. Biol.* | Not open access. Europe PMC gives PMID 30335370, no PMCID, "Subscription required". **The bioRxiv preprint, [doi:10.1101/322297](https://doi.org/10.1101/322297), was read instead**, and section 10's Methods quotes are from it |
| A PaqCI or AarI overhang-pair matrix, as a file | PaqCI's exists but is NEB's, unpublished and not downloadable from the tool. AarI's does not exist anywhere |
| NEB's three-base condition list | Not enumerated. See section 12 |

Two entries that were here have been **withdrawn**, because a browser reached both:

| What | What it turned out to be |
| --- | --- |
| `ligasefidelity.neb.com` and its help page returning 403 | **Bot-blocking only.** The pages load in a real browser, and the tools are `.cgi`, not `.html`: the Viewer is `/viewset/run.cgi`. Section 10 reads the condition list off it |
| `goldengate.neb.com` serving "a JavaScript shell" | The same mis-read. Nothing was missing from the markup; the fetch was being refused |

`nebridgetools.neb.com` still does not resolve, and that host does not exist.

All other sources read on 2026-09-12.

- Pryor, J.M., Potapov, V., Kucera, R.B., Bilotti, K., Cantor, E.J. and Lohman, G.J.S. (2020)
  Enabling one-pot Golden Gate assemblies of unprecedented complexity using data-optimized
  assembly design. *PLoS One* 15(9): e0238592.
  [doi:10.1371/journal.pone.0238592](https://doi.org/10.1371/journal.pone.0238592) (CC BY 4.0),
  its JATS XML, and its S1-S5 Tables
- Potapov, V. et al. (2018) Comprehensive profiling of four base overhang ligation fidelity by
  T4 DNA Ligase and application to DNA assembly. *ACS Synth. Biol.* 7, 2665-2674.
  [doi:10.1021/acssynbio.8b00333](https://doi.org/10.1021/acssynbio.8b00333) (CC BY-NC 4.0)
- Potapov, V. et al. (2018) Optimization of Golden Gate assembly through application of ligation
  sequence-dependent fidelity and bias profiling. *bioRxiv* 322297.
  [doi:10.1101/322297](https://doi.org/10.1101/322297) (CC BY-ND), the preprint of the paper
  above, read for what `HF`, `LF`, `DP` and `FP` name. Read 2026-10-06
- Potapov, V. et al. (2018) A single-molecule sequencing assay for the comprehensive profiling
  of T4 DNA ligase fidelity and bias during DNA end-joining. *Nucleic Acids Res.* 46, e79.
  [doi:10.1093/nar/gky303](https://doi.org/10.1093/nar/gky303) (CC BY-NC 4.0)
- NEB, *Ligase Fidelity Viewer: Help Page*:
  [tools.neb.com](https://tools.neb.com/~potapov/ligase-fidelity-viewer/help.html)
- NEB, *NEBridge Golden Gate Assembly Kit (BsmBI-v2)* instruction manual, NEB #E1602S/L, for
  what the NEBridge design tool guarantees about the overhangs it picks
- ECMA-376, *Office Open XML File Formats*, for the three parts of a workbook that hold a sheet
- `docs/research/golden-gate-assembly.md` (issue #4), which drew the licence line this note
  carries out, and `docs/research/restriction-enzyme-data.md` (issue #7) for the precedent:
  ship what the licence allows, cite the rest, and record what is unverified

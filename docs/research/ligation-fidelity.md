---
search:
  exclude: true
---

# Ligation fidelity data: sources, licence, and how a set of overhangs is scored

Research note for issue #11. Everything below was retrieved on **2026-09-12**. It records where
`src/liulab_mbio/data/ligation_fidelity.json` comes from, what its licence allows, how the
shipped matrices are read, and which overhang rules the code applies and on whose authority.

Issue #4's note (`docs/research/golden-gate-assembly.md`, section 4) surveyed the field and drew
the licence line. This note carries out what it decided and records what could be checked
against the paper's own published numbers.

## 1. What the data file holds

One file, `ligation_fidelity.json`, read with `importlib.resources` by
`liulab_mbio.goldengate.design`. It holds five count matrices, one per enzyme, plus the
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
`liulab_mbio/goldengate/design.py`, and here.

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
| `site` | an overhang no primer tail can carry without spelling a second site | `liulab_mbio.sites.primer_tail` | `avoid` |

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

PaqCI, AarI, BspQI, BpiI and BtgZI have no published matrix. Their overhangs are scored by the
rules instead, and `FidelityReport.measured` is `False` so that nothing mistakes the number for
a measurement.

The two weights, 0.05 for each near-duplicate partner and 0.02 for an overhang of one base
kind, are **ranking choices and not measurements**. They are module constants saying so. A set
scored this way compares with another set scored this way and with nothing else; it must never
be printed beside a measured fidelity as though the two were the same kind of number.

A tempting alternative was rejected: scoring PaqCI on BsaI's matrix, on the argument that
ligation fidelity belongs to T4 ligase and the cycling temperature rather than to the Type IIS
enzyme. It is a reasonable argument and it is still a substitution of one enzyme's measurement
for another's, so `choose_enzyme` prefers an enzyme that has its own matrix instead, and says
so in its ranking.

## 8. Regeneration

```sh
pixi run python scripts/build_ligation_fidelity.py            # downloads the five tables
pixi run python scripts/build_ligation_fidelity.py --from DIR # from files already downloaded
```

`journals.plos.org` serves each table from
`https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0238592.sNNN&type=supplementary`,
redirecting to signed storage, which `urllib` follows. Nothing else is fetched.

Tests never touch the network: `tests/test_build_ligation_fidelity.py` writes the three XML
parts of a workbook itself and reads them back, which is also what makes the guards testable.

## 9. Open gaps

| Item | Why it is missing | What is done instead |
| --- | --- | --- |
| The exact query set NEB's own tool uses | The Viewer v2 help page is 403 and the v1 page does not spell the arithmetic | The reading above, checked against three published numbers |
| Why two GetSet table sets score 1.5 and 4 points low | The sets as printed may not be the whole reaction | Recorded in section 5, untouched |
| T7 DNA ligase, and static 25 °C conditions | Only Potapov 2018 covers them, CC BY-NC | Cited, not shipped |
| A matrix for PaqCI, AarI, BspQI, BtgZI | Nobody has published one | The rule-based fallback, labelled as such |

## Sources

All read on 2026-09-12.

- Pryor, J.M., Potapov, V., Kucera, R.B., Bilotti, K., Cantor, E.J. and Lohman, G.J.S. (2020)
  Enabling one-pot Golden Gate assemblies of unprecedented complexity using data-optimized
  assembly design. *PLoS One* 15(9): e0238592.
  [doi:10.1371/journal.pone.0238592](https://doi.org/10.1371/journal.pone.0238592) (CC BY 4.0),
  its JATS XML, and its S1-S5 Tables
- Potapov, V. et al. (2018) Comprehensive profiling of four base overhang ligation fidelity by
  T4 DNA Ligase and application to DNA assembly. *ACS Synth. Biol.* 7, 2665-2674.
  [doi:10.1021/acssynbio.8b00333](https://doi.org/10.1021/acssynbio.8b00333) (CC BY-NC 4.0)
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

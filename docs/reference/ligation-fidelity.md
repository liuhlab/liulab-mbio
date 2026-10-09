# How fidelity is predicted

Every Golden Gate plan prints one fidelity figure, and so does every round of a library build.
It is a prediction. It says how cleanly the overhangs you chose should pair with each other and
with nothing else. This page says what the number counts, what stands behind it, and what it
does not tell you.

## What the number counts

The score is for a whole set of overhangs at once, never for one overhang alone.

A junction has two ends. One carries the overhang on the top strand, the other carries its
reverse complement on the bottom strand. For each junction the package reads a count matrix
twice. The correct count is those two ends ligating to each other. The total is those two ends
ligating to anything the tube holds, which is every overhang in the set and every reverse
complement of one. The junction scores correct over total. The set scores the product of its
junctions.

A longer set therefore scores lower than a short one even when every junction is clean, because
the product has more terms. One bad pair pulls down the whole set.

The report carries more than the figure. It lists any overhang whose own correct pair was seen
fewer than 100 times per 100,000 ligations, which is where NEB's Ligase Fidelity Viewer stops
calling a pair strong. It also lists every cross pair seen 10 times or more per 100,000, worst
first. Those cross pairs are the ones that would give you the wrong product.

## The measurements that ship

All five matrices come from one paper: Pryor, J.M., Potapov, V., Kucera, R.B., Bilotti, K.,
Cantor, E.J. and Lohman, G.J.S. (2020) Enabling one-pot Golden Gate assemblies of unprecedented
complexity using data-optimized assembly design. *PLoS One* 15(9): e0238592,
[doi:10.1371/journal.pone.0238592](https://doi.org/10.1371/journal.pone.0238592). It is CC BY
4.0, so the data ships with the package. [What ships](data.md) lists the file.

| Enzyme | Product measured | Overhang | Cycled between | Ligations counted |
| --- | --- | --- | --- | --- |
| BsaI | BsaI-HFv2 | 4 bases | 37 and 16 °C | 203,364 |
| BsmBI | BsmBI-v2 | 4 bases | 42 and 16 °C | 219,378 |
| Esp3I | Esp3I | 4 bases | 37 and 16 °C | 224,222 |
| BbsI | BbsI-HF | 4 bases | 37 and 16 °C | 189,282 |
| SapI | SapI | 3 bases | 37 and 16 °C | 172,332 |

Each matrix covers every overhang of its length against every other: 256 by 256 for the
four-base ones, 64 by 64 for SapI. Rows hold top-strand overhangs, columns hold bottom-strand
ones, and both are written 5' to 3'. A row's correct partner is the column spelling its reverse
complement, not the column of the same name. A pair never seen is left out of the file rather
than written as a zero.

A matrix belongs to the product named beside it, cycled between the two temperatures shown. The
same enzyme bought from somebody else does not inherit it.

## Which enzyme stands in for which

PaqCI, AarI, BspQI, BpiI and BtgZI have no published matrix. Rather than refuse a number, the
package borrows one and says whose.

It tries an enzyme that reads the same site first. An isoschizomer cuts the same site the same
way, so one measurement covers both enzymes under two names. Failing that, the enzyme takes the
shipped matrix of its own overhang length with the most ligations behind it.

| Enzyme | Site | Overhang | Scored on | Why |
| --- | --- | --- | --- | --- |
| BspQI | GCTCTTC | 3 bases | SapI | the same site |
| BpiI | GAAGAC | 4 bases | BbsI-HF | the same site |
| PaqCI | CACCTGC | 4 bases | Esp3I | most ligations at four bases |
| AarI | CACCTGC | 4 bases | Esp3I | the same |
| BtgZI | GCGATG | 4 bases | Esp3I | the same |

Pryor's Discussion is the authority for the second step: predicted fidelity "is unlikely to be
significantly impacted by the choice of Type IIS restriction enzyme". Both steps are worked out
from the shipped data, so a matrix added later changes this table without anyone editing a list.

## What a measured number is not

It is not a yield, not a pass mark, and not a measurement of your own reaction. A set at 87% can
still hand you a plate of correct colonies, and a clean set can fail for reasons ligation never
touches.

A matrix also belongs to the ligase and the conditions it ran under, not to the Type IIS enzyme.
The package keeps that straight in the label it prints beside every figure.

| Label | What scored the set |
| --- | --- |
| `measured` | the enzyme's own shipped matrix |
| `measured with Esp3I, not specific to PaqCI` | another enzyme's matrix standing in |
| `measured ligase profile, not specific to PaqCI` | a ligase matrix you supplied |
| `rule-based estimate` | no measurement covers it at all |

The last row is a ranking and not a prediction. Those rules take 0.05 off a junction for each
near-duplicate partner and 0.02 for an overhang of one base kind, and both weights are choices
rather than measurements. A set scored that way compares with another set scored that way, and
with nothing else.

## Bringing your own ligase matrix

Potapov et al. (2018) profiled T4 DNA ligase itself across all 256 four-base overhangs. Ligation
is the ligase's work, so that data speaks to every Type IIS enzyme, including the ones nobody
has measured. The archive is CC BY-NC 4.0. **None of it ships with the package and nothing from
it is redistributed.** You download your own copy, accept its terms yourself, and point the
package at the file.

One file is read, either `.xlsx` or `.csv`. It needs a header row of overhang labels, the same
labels down the first column, and a count in each cell. Anything else is refused with a message
saying what was expected. Where the file name states the conditions, as `FileS03_T4_18h_25C.xlsx`
does, they are read off the name and printed with the result.

| Command | Option | What it does |
| --- | --- | --- |
| `mbio cloning goldengate plan` | `--ligase-matrix` | scores an enzyme nobody has measured |
| `mbio cloning goldengate plan` | `--prefer-ligase-matrix` | scores on your matrix even where the enzyme has one of its own |
| `synbio igga plan` | `--ligase-matrix` | reports how often that ligase joins each round's overhangs; nothing is designed on it |
| `synbio igga plan` | `--ligase-sheet` | names the sheet to read, for a workbook holding several |

Both `--ligase-matrix` options fall back to an environment variable when you give no path:

```sh
export LIULAB_MBIO_LIGASE_MATRIX=/path/to/FileS03_T4_18h_25C.xlsx
```

Your matrix is the third choice, not the first. The enzyme's own matrix wins, then a stand-in,
then your profile. A stand-in measures the same kind of reaction and a pure ligation does not.
Pass `--prefer-ligase-matrix` to override that. The [command line reference](cli.md) lists these
options with the rest.

## Reading the score yourself

```python
from mbio.overhangs import fidelity

report = fidelity(["AGGT", "CAGC", "GCTT", "TACA"], "BsaI")
print(f"{report.value:.0%}", report.label)
```

That prints `87% measured`. The same four overhangs under PaqCI print `90% measured with Esp3I,
not specific to PaqCI`. The [Python reference](python.md) covers the rest of the report.

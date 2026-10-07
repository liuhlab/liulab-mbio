---
search:
  exclude: true
---

# Did the domestication break it? A QC panel, by failure mode

Issue #271, under the AP-1 demo spec #225, carrying #250 forward. This is the companion to
`lentiviral-tolerance.md`. That note says **what a working vector tolerates**. This one says
**how you find out whether the change broke anything** — which assay, how small a change it can
see, and, for each one, **what it would not catch**. The last of those three is the reason the
note exists.

It is written for the next working vector, whichever backbone that is. Nothing below depends on
pLVX, on BsaI, or on this lab's library.

Claims are marked, as in the companion note, and one mark is new:

- **Measured** — a primary source changed the thing and reported a number.
- **Figure-only** — the number exists only in a plot, because the paper's text states none. The
  value was read off the figure here, is approximate, and says so each time.
- **Calculated** — arithmetic done here on published values. Not a measurement.
- **Inferred** — a chain of reasoning from something measured.
- **Nothing measures this** — followed by what would close it.

Every source was read **2026-10-06**. Downloads are under
`reference_docs/synthesis_and_assembly/lentiviral-tolerance/`.

## 0. If you read nothing else

1. **The edit has a published precedent. The check does not.** One paper removed BsaI and BsmBI
   sites from lentiviral LTRs and the WPRE, and wrote in its own methods that it did not test the
   mutations for function. Section 1.
2. **"Silent" is not safe.** A two-base substitution in the HIV-1 5' UTR, in no reading frame,
   lowered genome packaging — and **neither base alone did anything**. Section 3.6.
3. **Run p24 and functional titre on the same harvest, and look at the ratio.** Across one
   published panel p24 varied 1.2-fold while titre spanned about 3,900-fold. Section 2.

## 1. The precedent is absent, and the absence is the finding

**No published account exists of anyone removing Type IIS sites from a lentiviral transfer vector
and reporting what they checked afterwards.** The companion note searched for this five ways —
for a BsaI-domesticated transfer vector, for a mutated `GGTCTC` in an LTR, through the mammalian
Golden Gate toolkits that ship lentiviral backbones, through the CRISPR-screen vector lineage,
and through the MoClo kits — and found the edit nowhere and the measurement nowhere.

Then it found the edit, and the edit makes the point sharper than the absence did.

**Haellman, V. *et al.* A versatile plasmid architecture for mammalian synthetic biology.
*Metab. Eng.* **66**, 41–50 (2021). PMID 33857582,
[doi:10.1016/j.ymben.2021.04.003](https://doi.org/10.1016/j.ymben.2021.04.003).** Its Tier-3
lentiviral vectors were, in the paper's words, "mutated to remove BsaI and BsmBI cut sites within
the LTRs' and WPRE enhancer". The same sentence continues:

> However, HindIII cut sites are too abundant to be easily removed, and **we refrained from
> testing all mutations for functionality.**

**The word "titer" does not appear anywhere in the paper.** Neither does "titre". There is no
transduction number for the domesticated backbone, against the parent or against anything else.

The sequence scan under `toolkits/VAMSyB-tier3-LTR-scan.txt`, run for #262 against the
Addgene-verified sequence of pTS1106_Tier3(Lenti) (Addgene 169661), says what was changed. Both
LTR copies carry two single-base substitutions relative to HXB2, R-relative positions **+4** and
**+80**: `GGTCTC` at R+1 becomes `GGTTTC`, which is the BsaI site gone, and `ATAAAGCTTGCCT` at
R+74 becomes `ATAAAGTTTGCCT`. **Measured**, from the sequence. The second change sits a couple of
bases past the R-region `AATAAA` hexamer and takes out a HindIII site — **inferred** from the
window the scan prints, not from any statement by the authors.

**So the state of the literature is: one group made a domestication-sized edit in the most
conserved part of a lentiviral vector, in both LTR copies, next to the poly(A) signal, and said
plainly that it did not test it.** Everything below is assembled from adjacent work — psi
mutagenesis, clinical-vector release assays, silencing studies, pooled-screen methodology — and
none of it was done on a domesticated vector.

## 2. The spine: two cheap numbers from one harvest

**Kim, S. H., Jun, H. J., Jang, S. I. & You, J. C.** The determination of importance of sequences
neighboring the psi sequence in lentiviral vector transduction and packaging efficiency.
*PLoS One* **7**, e50148 (2012). PMID 23185560, PMC3503997,
[doi:10.1371/journal.pone.0050148](https://doi.org/10.1371/journal.pone.0050148).

They inserted restriction sites around psi in a lentiviral transfer vector — the closest published
thing to a domestication edit — and measured particle output and function side by side. Table 1,
**measured**, quoted whole:

| Construct | p24 (ng/mL) | HT1080 colonies (CFU/mL) | Transduction, % of parent |
| --- | --- | --- | --- |
| EGFP (parent) | 298 ± 12 | 235,000 | 100 |
| ES2.psi | 252 ± 3 | 41,333 | 18 |
| ES3.psi | 293 ± 8 | 215,000 | 91 |
| PMM1.psi | 284 ± 6 | 215,000 | 91 |
| PBS.MAp | 282 ± 13 | 60 | 0.03 |

**p24 spans 1.18-fold across this panel. Functional titre spans about 3,900-fold.** That ratio is
**calculated** here from the table; the paper does not state it.

Three things follow, and they are the whole argument for the panel.

**p24 alone would have told them nothing.** Every construct made particles. One of them made 60
colonies per millilitre.

**The ratio localises the failure.** Normal particle output with collapsed infectivity is a
genome, packaging or reverse-transcription failure — something wrong with what is inside the
particle, or with what it does after entry. Low p24 *and* low titre is a production failure, and
you go and look at the transfection. One pair of numbers from one harvest splits the lifecycle in
half.

**Position matters at a resolution of a few bases.** ES2 put an EcoRV site at the 5' end of psi
and reached 18% of parent. ES3 left the 5' end intact, used the same 3' site, and reached 91%.
ES1 and ES2 differ from each other by **three nucleotides** at that 5' junction and both land near
20% of control (body text, MT4 transduction, **measured**). PMM1 and PMM2, at HXB2 positions 203
and 197, were both tolerated. **A few bases decide it, and nothing about the sequence told them
which few.**

## 3. The panel, by failure mode

Each entry gives the assay, what it has been measured to resolve, and what it is blind to.

### 3.1 Lost titre

Lost titre is the failure everyone tests for and the one that hides best in a pooled library. You
transduce the planned number of cells, get fewer transductants than you think, and believe you are
at 500-fold coverage when you are at 100-fold.

**The assay.** Functional titre on the same harvest as p24. Flow cytometry if the vector carries a
marker; otherwise ddPCR on genomic DNA from transduced cells. **Kandell, J., Milian, S., Snyder,
R. & Lakshmipathy, U.** *Mol. Ther. Methods Clin. Dev.* **31**, 101120 (2023). PMID 37841416,
PMC10568280, [doi:10.1016/j.omtm.2023.101120](https://doi.org/10.1016/j.omtm.2023.101120) —
targets **RRE** duplexed with GAPDH, so one assay serves every member of a library regardless of
transgene.

**Measured sensitivity.** LLOQ **1.16 copies/µL**, which the paper converts to 12 copies per µg of
genomic DNA; linear range 1.16 to 1.11 × 10³ copies/µL; %CV 5 to 17 across that range (Table 1,
three independent runs). Against the NIST single-copy reference it read **0.89 copies per cell,
3.6% CV, n = 3**, an 11% accuracy difference. For vector copy number, **Corre, G. *et al.***
*Gene Ther.* **29**, 536–543 (2022). PMID 35194185, PMC9482878,
[doi:10.1038/s41434-022-00315-8](https://doi.org/10.1038/s41434-022-00315-8) — against clonal
standards at VCN 1, 2 and 3, **ddPCR CV 4–8% against qPCR 10–15%**, LOD **0.005 VCN**, LOQ
**0.01 VCN** at 20–40 ng of genomic DNA per reaction. A 2-fold change is resolvable at n = 2.

**What it would not catch.**

- **Truncation, silencing or barcode linkage.** None of them changes titre by definition, and
  section 4 has a case where titre moved 10–30% while the RNA measurement moved 70%.
- **Which half of the lifecycle broke** — not without p24 beside it.
- **Carried-over transfer plasmid.** Corre showed directly that a psi-targeting assay amplifies
  the transfer plasmid as well as the provirus, "as confirmed by a positive PCR result". Their
  provirus-specific design puts the probe across the **U5-to-vector junction that only exists
  after reverse transcription duplicates the LTR**, and matched the psi assay's precision.
  Kandell's RRE assay is **not** described as provirus-specific and the paper states no
  carry-over caveat. **Wu, X. *et al.*** *Heliyon* **10**, e38512 (2024). PMID 39498040,
  PMC11532291, [doi:10.1016/j.heliyon.2024.e38512](https://doi.org/10.1016/j.heliyon.2024.e38512)
  handles the same problem at the bench instead: 100 U/mL Benzonase for 1.5 h drove spiked
  transfer-plasmid copies to zero.
- **Agreement with flow is not guaranteed.** Kandell reports ddPCR "generally within 30%" of flow,
  but its own Table 2 holds one sample at **181% difference** (attributed to a flow staining
  reagent). Wu reports ddPCR titres **4–5-fold higher than FACS** (p < 0.05, n = 3) and treats
  that as expected, not as error. **Two titre methods are not interchangeable; pick one and keep
  it.**

### 3.2 Truncated or mis-ended genome

A fraction of the RNA packaged into particles is not full length. This is silent in every bulk
assay, including titre.

**The assay.** Nanopore direct RNA sequencing of RNA extracted from purified virions. Every read's
3' terminus is the RNA's own 3' end, so truncation and cryptic polyadenylation sites fall out
directly rather than being inferred.

**Zeglinski, K. *et al.*** An optimized protocol for quality control of gene therapy vectors using
nanopore direct RNA sequencing. *Genome Res.* **34**, 1966–1975 (2024). PMID 39467647,
PMC11610601, [doi:10.1101/gr.279405.124](https://doi.org/10.1101/gr.279405.124).

> The file holding this paper is named `fm2-pal2024-...`. **The first author is Zeglinski.** The
> filename is wrong and the citation above is read off the article itself.

**Measured sensitivity.** On clinical WAS vectors, about **80%** of reads ended at the expected
3' end and about **10%** ended at roughly 4,600 bp, at a non-canonical `ATTACA` inside the
**WPRE**. With the artificial-polyadenylation arm, which recovers RNAs that lost their poly(A)
tail, one vector resolved: WPRE cryptic poly(A) **4.44%**, a hairpin truncation at about 500 bp
in the psi stem-loops **11.26%**, a second at about 2,300 bp at an shRNA cassette **8.61%**
(19.88% together), and spliced reads **1.31%**. **Full-length RNA was only about 60–75% even in a
working vector.** The smallest species the paper quantifies is around 1–2%, so treat roughly 2% as
the practical floor — that is **inferred** from what they report, not a stated detection limit.

The paper's repair attempt is the most useful thing in it for domestication. Mutating `ATTACA` to
`ATTTCA` and removing the shRNA cut the truncation window to **3.39%**, then **2.04%** — but the
authors write that mutation "did not eliminate premature termination", and the repaired vector
**gained a new splice event that removed most of the promoter in 19.35% of reads**, which took two
further base changes to remove. **One deliberate base change created a second defect that only
this assay could see.**

**What it would not catch.**

- **Function.** No titre was measured. The authors say so.
- **Non-polyadenylated truncations**, unless the artificial-polyadenylation arm is run as well.
  Run both.
- **Base-level changes.** SNP calling on RNA002 chemistry is unreliable; errors and modified bases
  give miscalls and indels. **Do not use this to confirm your edit — sequence the plasmid.**
- **Anything, cheaply.** Only **1.5–3%** of reads aligned to the vector, and a MinION flow cell
  yields about **10,000 vector reads**. This is the most expensive assay in the panel.

Read-through past the 3' LTR is the same failure mode seen from the other end, and for
domestication it carries a specific warning: **termination and promoter are the same bases.**
**Yang, Q., Lucas, A., Son, S. & Chang, L.-J.** *Retrovirology* **4**, 4 (2007). PMID 17241475,
PMC1802088, [doi:10.1186/1742-4690-4-4](https://doi.org/10.1186/1742-4690-4-4) measured it with a
Cre-loxP reporter line and found that **70–80% of termination activity lives in a 124-nt stretch
overlapping the NF-κB, Sp1 and TATA sites**, and that deleting NF-κB, NF-κB plus Sp1, or TATA
alone **each raised read-through**. Inserting alternative terminators — HTLV-1 poly(A), a 3'
intron, a tRNA motif — restored nothing and made read-through worse.

**Two limits on that number, and they matter here.** The paper's text gives **no absolute
read-through percentage** for any construct; every value is a bar in a figure, and the 70–80%
figure appears only as a relative statement in the abstract and discussion. And **every mutation
Yang made was a block deletion or a block restoration — no base substitution was ever tested.**
That is exactly the resolution a domestication edit works at, and it is unmeasured.

For an absolute frequency, use **Koldej, R. M. & Anson, D. S.** *BMC Biotechnol.* **9**, 86 (2009).
PMID 19811661, PMC2765960,
[doi:10.1186/1472-6750-9-86](https://doi.org/10.1186/1472-6750-9-86), whose TaqMan assay reads a
transcript spanning the 3' LTR poly(A) into the plasmid backbone against total transgene
transcript. Their Table 1, **measured**, n = 3: the parent vector sits at
**9.05 ± 0.93 × 10⁻³** — about 0.9% read-through — and the designs that replaced the HIV-1 poly(A)
with a heterologous one reached **8.7 × 10⁻²**, nearly ten times worse. That is the assay to copy,
and the result is a warning about "improving" the poly(A) while domesticating near it.

**Zaiss, A. K., Son, S. & Chang, L.-J.** *J. Virol.* **76**, 7209–7219 (2002). **PMID 12072520**,
PMC136337,
[doi:10.1128/jvi.76.14.7209-7219.2002](https://doi.org/10.1128/jvi.76.14.7209-7219.2002) is the
source of the method and of the finding that a SIN HIV-1 poly(A) is already as leaky as MLV's.
**The PMID is corrected here.** An earlier copy of the source file held the abstract of PMID
12072521 — Lecellier *et al.*, equine foamy virus, the adjacent article in the same issue. Use
12072520.

### 3.3 Silent provirus

The provirus integrates, the member is in the pool, it amplifies normally in a barcode readout,
and it contributes no phenotype. **It reads as "this variant is inactive", not "this variant is
broken".** For a screen this is the worst failure in the note.

**The assay.** Vector copy number and percent expressing, **on the same cells, at four or more
timepoints**. The gap between them is the signature.

**Herbst, F. *et al.*** *Mol. Ther.* **20**, 1014–1021 (2012). PMID 22434137, PMC3345972,
[doi:10.1038/mt.2012.46](https://doi.org/10.1038/mt.2012.46). Vector copy number per cell stayed
flat at **1.1 ± 0.004 falling to 0.9 ± 0.3** (not significant) while eGFP-positive marrow fell
from **5.4–17.6% at 4 weeks to 0.0–0.6% at 36 weeks**. Bisulfite sequencing found **18 of 18
promoter CpGs over 90% methylated** against under 3% in controls, and 5-azacytidine recovered
expression from 7.9 ± 0.7% to 14.2 ± 0.1%, about 1.8-fold — partial, not full. **Promoter choice
dominates the problem**: CMV fell 3.6–5.2-fold during differentiation, PGK and EF1α only
1.5–1.7-fold.

**Baseline, so you know what is normal.** **Jordan, A., Bisgrove, D. & Verdin, E.** *EMBO J.*
**22**, 1868–1877 (2003). PMID 12682019, PMC154479,
[doi:10.1093/emboj/cdg188](https://doi.org/10.1093/emboj/cdg188): productive infection 4% against
a latent phenotype at 0.06%, so **about 1.5% of integrations — one in 66 — land silent with a
perfectly good vector**, much of that driven by integration site rather than vector sequence.
Rescue is the confirmation step, and it is unreliable: TNF-α restored expression in all clones,
but the HDAC inhibitor trichostatin A worked **only in some**.

**Contreras, X. *et al.*** *J. Biol. Chem.* **284**, 6782–6789 (2009). PMID 19136668, PMC2652322,
[doi:10.1074/jbc.M807898200](https://doi.org/10.1074/jbc.M807898200) sets the scale of an HDAC
rescue: **5-fold** in resting primary CD4+ cells, up to **30-fold** in a reporter line, **40-fold**
in a chronically infected line.

> **Provenance, marked.** The Herbst, Jordan and Contreras numbers were read from the PMC page
> text of each article, which attributes them to figures that were not opened here. They are
> first-hand from the article text and second-hand from the figures. No OA PDF or full-text XML is
> served for any of the three.

**What it would not catch.**

- **A negative rescue does not prove the provirus is absent.** Trichostatin A fired on only some
  clones in Jordan's own hands.
- **It cannot run on a pool.** Copy number and percent expressing are population averages. A pool
  in which 5% of members are fully silenced looks like a pool at 95% expression, and so does a
  pool in which every member is 5% dimmer. Section 5 step 6 has the fix.
- **It is slow.** The Herbst signature needed 36 weeks to become obvious. Four timepoints over six
  weeks is the shortest version worth running, and that is **inferred** from the shape of their
  curve, not a published minimum.

### 3.4 Lost or skewed expression

For a Tet-on vector, two numbers that must be reported separately: induced level, and leak.

**The assay.** **Loew, R., Heinz, N., Hampf, M., Bujard, H. & Gossen, M.** *BMC Biotechnol.*
**10**, 81 (2010). PMID 21106052, PMC3002914,
[doi:10.1186/1472-6750-10-81](https://doi.org/10.1186/1472-6750-10-81) measures the same promoters
twice — transiently, and after single-copy integration from a SIN vector at MOI 0.06–0.12 — which
is the comparison a vector lab needs.

**Measured sensitivity.** **Regulation collapses on integration.** Transiently, Ptet-T6 reaches
about **50,000-fold** regulation, against Ptet-1 at 400-fold and Ptet-14 at 6,200-fold. At single
integrated copy, Figure 4B gives **4.5 × 10⁴ (± 2.8 × 10³) to 1.6 × 10⁷ (± 4.7 × 10⁵) rlu/µg for
Ptet-1 and 2.9 × 10³ (± 1.1 × 10²) to 2.3 × 10⁷ (± 3.9 × 10⁵) for Ptet-T6** — about **356-fold**
and **7,900-fold**, which is **calculated** here from the figure's own values. The body text
states only "15-fold" lower background and "about 30-fold" wider range.

**Flow cytometry cannot measure the leak, and this is the measurement that proves it.** Figure 4C
was read off the figure here. **Figure-only**, mean fluorescence of the whole population:

| Dox (ng/mL) | Ptet-1 (mfu) | Ptet-T6 (mfu) |
| --- | --- | --- |
| 0 | 2.2 | 1.9 |
| 10 | 6.9 | 4.7 |
| 30 | 94 | 179 |
| 100 | 266 | 409 |
| 1000 | 402 | 520 |

**The parental GFP-negative background is 1.83.** So the off state of Ptet-1 sits at 1.2-fold
above background and Ptet-T6 at 1.04-fold above it. The authors state outright that **no
meaningful conclusion can be drawn** from the GFP data about leakiness. Luciferase on the same
populations resolved the 15-fold difference immediately. By flow the two promoters look like
183-fold and 274-fold regulation; by luciferase they are 356-fold and 7,900-fold. **Flow
understates Ptet-T6's regulation range by roughly 29-fold** — **calculated** here from the figure.

**What it would not catch.**

- **The leak**, on any cytometer, when the off state is near background. If the leak matters, put
  a luminescent reporter on the construct.
- **Per-cell heterogeneity at partial induction.** All cells responded at 30 ng/mL Dox, but not
  homogeneously — a population mean hides that, and flow is the only thing that shows it. The two
  assays are blind to opposite things.
- **The smallest fold change a cytometer can resolve.** No source states one; it is instrument-
  and fluorophore-specific. **Nothing measures this.** *Closes it:* a 2-fold ladder on the lab's
  own machine.

### 3.5 Skewed representation

The screen runs and the answer is wrong.

**The assay.** Sequence the barcode or guide amplicon at the **plasmid pool**, after transduction,
and at the endpoint. **Normalise to the plasmid pool, not to the first cell timepoint.**

**Imkeller, K., Ambrosi, G., Boutros, M. & Huber, W.** *Genome Biol.* **21**, 53 (2020). PMID
32122365, PMC7052974,
[doi:10.1186/s13059-020-1939-1](https://doi.org/10.1186/s13059-020-1939-1) is the source of that
rule, and gives the reason: guides against essential genes are **already depleted by T0**, up to
four doublings after transduction, so a T0 reference bakes the effect you are measuring into the
denominator.

**Measured acceptance bar.** **Joung, J. *et al.*** *Nat. Protoc.* **12**, 828–863 (2017). PMID
28333914, PMC5526071, [doi:10.1038/nprot.2017.016](https://doi.org/10.1038/nprot.2017.016):
"An ideal sgRNA library should have **more than 70% perfectly matching guides, less than 0.5%
undetected guides, and a skew ratio of less than 10**", where the skew ratio is the 90th to 10th
percentile abundance ratio, judged at over 100 reads per member. Screen at **over 500 cells per
member** and **MOI under 0.3**, so most cells take one.

**Coverage against skew.** Imkeller's Table 2 pairs the two, **measured by simulation**:

| Library width (p90/p10) | Recommended coverage |
| --- | --- |
| 2.5 | 200 |
| 5 | 300 |
| 7.5 | 300 |
| 10 | 400 |
| 17 | 400 |

**Heo, S.-J. *et al.*** *Genome Biol.* **25**, 25 (2024). PMID 38243310, PMC10797759,
[doi:10.1186/s13059-023-03132-3](https://doi.org/10.1186/s13059-023-03132-3) pushes it further: a
library with a 90/10 skew ratio **under 2** screened at 100-fold coverage matched a 1,000-fold
screen (phenotype-score r² = 0.900; BAGEL2 AUC 0.949 against 0.937), and the paper concludes that
**50-fold cell coverage** suffices for such a library. **Its coverage titration at 200, 100, 50
and 10-fold is described only qualitatively in the text — the per-coverage skew numbers live in
Fig. 4B/C and Table S5, which are not in the gathered file and are therefore not quoted here.**

**What it would not catch.**

- **A silenced member.** It amplifies from genomic DNA exactly as a working one does. Barcode
  abundance is a count of integrations, not of function.
- **A quiet titre loss.** Section 7 works the arithmetic: a 2-fold titre loss is invisible in the
  skew metric. Imkeller measured the same thing from the other side — reducing **cell-splitting**
  coverage pushed up to **20%** of guides below an LFC of −1, while the same reduction in
  **transduction** coverage moved only **3%**. Transduction is not where screens lose their
  representation.
- **Barcode-to-part linkage.** Nothing in an amplicon readout says the barcode still names what it
  was designed to name. Only long-read sequencing of the plasmid pool says that.

### 3.6 Packaging and dimerisation

This is the section that justifies assaying rather than reasoning, and it is the strongest result
in the material.

**Sakuragi, S. *et al.*** Identification of a novel cis-acting regulator of HIV-1 genome packaging.
*Int. J. Mol. Sci.* **22**, 3435 (2021). PMID 33810482, PMC8036536,
[doi:10.3390/ijms22073435](https://doi.org/10.3390/ijms22073435).

**The edit.** Two consecutive bases, **`GA` to `AC` at NL4-3 positions 226 and 227**, inside psi,
downstream of the primer binding site, in the 5' UTR. Two bases. The paper places them in the 5'
UTR, so no reading frame covers them — **inferred**, one step, from the paper's own description of
the region.

**The assay.** A **competitive RT-qPCR packaging ratio**. Each mutant is co-transfected with a
control genome that differs only by silent changes in CA, and the measurement is the mutant-to-
control ratio in virions normalised to the mutant-to-control ratio in the cytoplasm. Transfection
efficiency and RNA recovery cancel. n = 3, SEM.

**Measured — but the magnitude is figure-only.** The paper's text gives **no fold number at all**,
calling the effect "low packaging ability" and "moderate rather than lethal". Figure 3 was fetched
and read here. Relative packaging efficiency, approximate values read off the plot:

| Construct | Bases changed | Relative packaging |
| --- | --- | --- |
| NLAINh (parent) | none | ~0.95 |
| #8-1 | one of the two critical bases | ~1.05 |
| #8-8 | the other critical base, in a four-base block | ~1.35 |
| **#8-2** | **both, `GA` to `AC`** | **~0.80** |
| #8 | all five discordant bases | ~0.60 |

Panel 3C, same reading: the segment carrying them drops to **~0.47** against a parent at **~0.88**.

**The finding.** The two-base change lowers packaging. **Neither single base alone does
anything** — the paper says it twice, in the results and in the discussion: "neither single-base
changes lead to phenotype conversion" and "There was no difference in packaging ability when only
one base was substituted." **A domestication edit that changes two adjacent bases is not two
single-base edits, and testing them one at a time would have found nothing.**

**Dimerisation is the second assay, and it is not the same assay.** **Clever, J. L. & Parslow,
T. G.** *J. Virol.* **71**, 3407–3414 (1997). PMID 9094610, PMC191485: SL1 loop mutants that block
dimerisation *in vitro* have **impaired infectivity at normal virion RNA content**. A packaging
ratio alone would call those vectors fine. The separating assay is a **non-denaturing-gel
Northern** on virion RNA. Clever also found that the **structures, not the sequences**, of the SL1
and SL3 stems are what packaging needs — which is what makes a compensating change thinkable at
all.

**What it would not catch.**

- **A dimerisation defect**, if only the packaging ratio is run. Clever is the counter-example.
- **A packaging defect**, if only the gel is run. Sakuragi measured no dimerisation change for the
  subtype swap and never ran the gel on the 226/227 mutants themselves.
- **Function.** Sakuragi reports **no infectivity or titre number** for these mutants; p24 was used
  only to normalise input. A packaging drop of this size has never been converted into a titre.
  **Nothing measures this.**

### 3.7 Recombination and repeat instability

Two different failures share this heading and they want different responses.

**Between the LTRs, in the plasmid prep.** **Rhode, B. W., Emerman, M. & Temin, H. M.** *J. Virol.*
**61**, 925–927 (1987). PMID 3806800, PMC254040: direct repeats of 0.85–1.3 kb deleted at high
frequency during virus replication, whether tandem or separated. **The paper states no rate**, it
is MLV-based, and it is the only primary measurement in the material. The assay is whole-plasmid
sequencing of the prep you are about to use — not of the design.

**During reverse transcription, between co-packaged genomes.** **Levy, D. N., Aldrovandi, G. M.,
Kutsch, O. & Shaw, G. M.** *PNAS* **101**, 4204–4209 (2004). PMID 15010526, PMC384719,
[doi:10.1073/pnas.0306764101](https://doi.org/10.1073/pnas.0306764101): **an average of nine
recombination events per virus** in a single round in T lymphocytes, and about **30 crossovers**
in macrophages. **But recombination requires two different genomes in one particle, so it scales
as the square of the infection rate** — with MOI and pooling, not with your edits. **The
mitigation is low MOI, not sequence design.**

**SIN repair is the case where sequence design does not help either.** Koldej & Anson measured it
by marker rescue and colony PCR, Table 2, **measured**: the parent repaired its 5' LTR in **17 of
17** colonies; a vector with U3 homology cut to the 19 bp att repaired in **13 of 22** with 18%
still SIN; a vector with **zero U3-to-U3 homology** still repaired completely in **2 of 11**, with
another 42% partial. In their words, "repair still occurred in the absence of any homology between
the U3 regions". **Reducing homology is not a cure.** That paper is also the closest thing in the
literature to the QC panel this note is describing — titre at n = 18, read-through by TaqMan across
the vector-to-flank junction, SIN repair by colony PCR — and it is worth copying wholesale. It is
**not**, as is sometimes said, a codon recoding: the changes were splice-donor mutations, a U3
deletion, a heterologous poly(A) and a 76 bp *gag* deletion.

**What it would not catch.** Whole-plasmid sequencing of the prep sees a deletion that already
happened in *E. coli*; it says nothing about what reverse transcription will do in the target
cells. Nothing in the panel measures per-library recombination rate for a modern SIN transfer
plasmid. **Nothing measures this.**

### What each assay is blind to, in one table

The point of the panel is this column, so here it is alone.

| Assay | Blind to |
| --- | --- |
| p24 ELISA | Everything that happens after the particle leaves the cell. Across Kim's panel it varied 1.2-fold while titre varied 3,900-fold |
| Functional titre | Truncation, silencing, barcode linkage, dimerisation, read-through. And a 10–30% titre change can hide a 70% RNA change (section 4) |
| ddPCR copy number | Whether the copies express. With a psi-targeting design, also whether they are provirus or carried-over plasmid |
| Nanopore direct RNA sequencing | Function; non-polyadenylated truncations without the second arm; base-level edits |
| Read-through TaqMan | Substitutions — every published dissection used block deletions |
| VCN against percent expressing | Which member of a pool went silent. Needs sorting |
| Barcode or guide NGS | Silenced members, and whether a barcode still names its part |
| Packaging ratio | A dimerisation defect at normal RNA content |
| Non-denaturing Northern | A packaging defect |
| Whole-plasmid sequencing | Anything that happens during reverse transcription |

## 4. Why that column is the point: titre and RNA come apart

**Cui, Y., Iwakuma, T. & Chang, L.-J.** Contributions of viral splice sites and cis-regulatory
elements to lentivirus vector function. *J. Virol.* **73**, 6171–6176 (1999). PMID 10364378, PMC112687,
[doi:10.1128/jvi.73.7.6171-6176.1999](https://doi.org/10.1128/jvi.73.7.6171-6176.1999).

They mutated the major splice donor from `GGTG` to `GCAG` or `GGGG`. **Measured**, and quoted:

> Mutations of SD from GGTG to GCAG (SD1) or to GGGG (SD2) resulted in **a moderate decrease (10
> to 30%) in vector titer**. […] mutations of SD **reduced expression of cytoplasmic full-length
> RNA by more than 70%** compared with wild-type vector. […] The decline in genomic RNA expression
> (70%), however, did not correlate with the reduction in vector titer (10 to 30%).

Titre averaged 4 to 7 independent experiments. **A four-base change cost almost nothing in titre
and most of the full-length cytoplasmic RNA.** If titre were the only assay, that vector passes.

Three more decouplings, all above, all pointing the same way: Zeglinski's working clinical vectors
carry only 60–75% full-length genome and titre never said so; Loew's cytometer could not see a
15-fold difference in promoter leak that a luminometer read off the same cells; Sanjana's
lentiCRISPRv2 gained roughly 10-fold in functional titre over v1 and lentiGuide-Puro roughly
100-fold (**Sanjana, N. E., Shalem, O. & Zhang, F.** *Nat. Methods* **11**, 783–784 (2014). PMID
25075903, [doi:10.1038/nmeth.3047](https://doi.org/10.1038/nmeth.3047)) — a vector-design change
whose entire published characterisation is a titre.

**One assay answers one question. The panel exists because the questions are not nested.**

## 5. The panel, cheapest first

1. **Whole-plasmid sequencing of the new prep and of the parent.** Catches every unintended edit
   and any rearrangement between the LTRs, before a cent is spent on virus. It is also the only
   step that catches an element that has reverted in the prep. **Sequence the prep, not the
   design.**
2. **One paired packaging run — parent and domesticated, same day, same mix, n ≥ 3 — giving p24
   and functional titre from the same harvest.** Two numbers; their ratio splits the lifecycle in
   half. This is the step that would have caught every effect in Kim's table.
3. **ddPCR vector copy number on transduced cells, provirus-specific.** Resolves a 2-fold change
   at n = 2 (Corre: CV 4–8%). Use a design that cannot read the transfer plasmid, or digest the
   plasmid away as Wu does.
4. **Flow for expression. For a Tet vector, induction and leak as separate numbers — and a
   luminescent reporter if the leak matters**, because Loew measured that flow cannot see it.
5. **Percent expressing divided by copy number, at four or more timepoints over six weeks**, for
   silencing. Baseline: about 1 integration in 66 goes silent with a good vector.
6. **For a library, always: barcode NGS at the plasmid pool and after transduction, normalised to
   the plasmid pool.** Plus one cheap addition that closes the pooled blind spot in section 3.3 —
   **sort expressing from non-expressing cells out of the same pool and sequence barcodes from
   both.** A member enriched in the non-expressing fraction at unchanged total abundance is
   silenced, not lost. No source does this; it follows from Herbst's signature applied to a pool.
7. **Conditional on where the edits landed.** The packaging ratio if anything sits in the 5' UTR
   or psi. Direct RNA sequencing if anything sits in WPRE, psi, or near a splice site. Read-through
   TaqMan if anything sits in U3. A non-denaturing Northern if anything sits in SL1.

**The design rule that outranks the panel: keep every comparison paired and same-day against the
undomesticated parent.** Every effect in this literature is a ratio to a parent, and batch
variance in lentiviral packaging is larger than any of these assays' CV.

## 6. If you can afford exactly one assay

**Run step 2 — p24 and functional titre from one paired harvest — and nothing else.**

Not because it is the most informative. It is not: it would have missed Cui's splice donor, every
one of Zeglinski's truncations, and all of Herbst's silencing. Run it because it is the only
single assay that **both** detects a large class of failures **and** tells you where to look next.
A titre loss with normal p24 sends you to sections 3.2 and 3.6. A loss in both sends you back to
the transfection. Any other single assay returns a number with nowhere to go.

**And if you cannot afford that**, run step 1 instead and do not make virus. A prep that is not
the sequence you designed makes every later number meaningless, and sequencing is the cheapest
thing here.

## 7. Does a quiet titre loss actually skew a pooled library?

Less than instinct suggests. **Calculated, not measured.** Titre enters a pooled library only
through coverage, and counts are Poisson, so at 500-fold coverage the sampling contribution to the
90/10 skew ratio is about **1.12**; at 250-fold — a 2-fold titre loss — about **1.18**; at 50-fold
about **1.44**. Against a cloning-derived skew of 2 to 10, **a 2-fold titre loss is invisible in
the skew metric.** The Poisson floor starts to dominate only below roughly 10 to 25-fold, which is
a 20- to 50-fold titre loss.

Two independent results agree with it. Heo's library held its skew from 200-fold down to 50-fold
coverage. Imkeller's simulations found that cell splitting, not transduction, is the dominant bias
term — 20% of guides below LFC −1 from reduced splitting coverage against 3% from the same
reduction at transduction.

**So the risk from a quiet titre loss is operational, not statistical.** You do not see it in the
skew. You see it as believing you were at 500-fold when you were at 100-fold, and then the
split-and-grow bottleneck — which *is* the dominant term — acts on a smaller population than
planned. **Measure the titre; the skew will not tell you.**

## 8. What nothing measures, and the experiment that closes it

**Nothing in the literature says what a domestication-sized edit costs a lentiviral vector.** Not
one number. The nearest approaches are Kim's restriction-site insertions around psi, Sakuragi's
two bases, and Cui's four — none of them framed as domestication, none in an LTR, and the one
group that did edit the LTRs stated that it did not test them.

**The experiment that would answer it**, and it is one packaging run:

Build the domesticated vector **and each single-site intermediate** — one plasmid per removed
site, each carrying that change alone against the parent backbone. Produce them in **one particle
pool**, same day, same mix, n ≥ 3, with the undomesticated parent in the same run. Measure the
cheap assays only: p24, functional titre, and provirus-specific copy number on transduced cells.

**The intermediates are the experiment.** The full domesticate alone gives one number and no
attribution: if titre drops 3-fold you know nothing about which of six changes did it, and
Sakuragi's result says you cannot reason your way to the answer either. The intermediates convert
one number into a per-site cost, which is reusable on the next backbone — and they are the only
design that distinguishes "this site is expensive" from "these two changes interact", which is the
specific failure mode section 3.6 documents.

**Two things to decide before running it.** Each single-site intermediate costs a synthesis and a
prep, so the panel is as large as the site count — which is the argument from the companion note's
section 1 for **choosing the enzyme before choosing the edit**. And edits that fall in the same
element should be built both separately and together, because the two-base result says the
combination is not the sum.

**Other things nothing measures here.**

- **The per-base sensitivity of 3'-LTR read-through.** Yang measured block deletions. A
  substitution has never been tested. *Closes it:* substitute, do not delete, and run Koldej's
  TaqMan ratio.
- **A 5'-to-3' RT-qPCR ratio on packaged RNA.** No primary source defines it, calibrates it, or
  states a limit. *Closes it:* spike a 3'-truncated transcript into full-length RNA at 0, 5, 10,
  25 and 50% and find the smallest distinguishable fraction.
- **The smallest fold change a cytometer resolves.** Instrument-specific; no source states one.
- **How a silenced member behaves in a pooled barcode readout.** Nobody measured it. The sort in
  step 6 would.
- **Plasmid-level recombination frequency for a modern SIN transfer plasmid.** The only primary
  measurement is from 1987, MLV-based, at kilobase scale, and states no rate.
- **What a packaging-ratio drop costs in titre.** Sakuragi never titred the mutants.

## 9. Sources, and what was left out

Twenty-eight files were gathered under
`reference_docs/synthesis_and_assembly/lentiviral-tolerance/assays/` and all twenty-eight were
read first-hand for this note. Every citation above is read off the article text, not recalled.
Two sources outside that directory are cited and were read the same way: Cui 1999 from
`elements/`, and Haellman 2021 with its sequence scan from `toolkits/`.

**Where the gathered material is weaker than it looks, stated plainly:**

- **Three sources are PMC page text with figure-referenced numbers.** Herbst 2012, Jordan 2003 and
  Contreras 2009 serve no OA PDF or full-text XML. Their numbers are marked where used.
- **One paper's numbers live only in its figures.** Sakuragi 2021 states no fold change anywhere
  in its text. Figure 3 was fetched from PMC and read here, and every value taken from it is
  marked figure-only and approximate.
- **One figure set was read directly rather than quoted.** Loew 2010's Figure 4 was rendered from
  the saved OA PDF and read, which is where the mfu table in section 3.4 comes from.
- **Four of the pooled-screen files define no skew metric at all** — Doench 2016, Hart 2017,
  Shalem 2014 and Cross 2016 — and are used here only for coverage and MOI. Their supplements,
  which hold the methods, are not in the gathered files.
- **Heo 2024's per-coverage skew numbers are in supplementary tables that were not gathered.** The
  coverage titration is reported here only as the paper's own qualitative statement.
- **Wu 2024's precision tables are empty in the extracted text**, so only the stated bound
  (CV ≤ 25%) is quoted. The same file states a titre range whose two endpoints disagree by a
  factor of 100; neither endpoint is used here.
- **Yang 2007 states no absolute read-through percentage.** Every value is a bar. Section 3.2 uses
  Koldej's TaqMan ratio for an absolute number instead.
- **One PMID is corrected.** Zaiss *et al.* 2002 is **PMID 12072520**, not 12072521 — which is
  Lecellier *et al.*, equine foamy virus, the adjacent article in the same issue. An earlier copy
  of the source file held the wrong abstract.
- **One filename is wrong.** `fm2-pal2024-nanopore-drs-vector-rna-integrity.txt` holds Zeglinski
  *et al.* 2024. The file was not renamed; the citation here is correct.
- **No vendor material enters this note.** The search-snippet files under
  `domestication-methods/verification-and-repeats/` — Thermo Stbl2 and Stbl3, NEB Stable,
  Plasmidsaurus, Eurofins, Twist and GenScript — are labelled unverified in their own text and are
  excluded by name by their own README. Nothing above cites one for a value. Where they would have
  supplied a number, section 8 records a hole instead.

## 10. What this note did not settle

- **Whether VAMSyB's domesticated LTRs work.** The paper says it did not test them. The plasmids
  are on Addgene. *Closes it:* make virus from pTS1106_Tier3(Lenti) and its parent in one paired
  run and titre both — which is step 2 of section 5, and would be the first published number of
  its kind.
- **Whether the R+80 change next to the poly(A) hexamer is free.** It sits in the element Yang and
  Zaiss both show is leaky already. Nobody has measured a substitution there.
- **What any of these assays cost.** No price, instrument time or hands-on time appears here,
  because none of it is in the sources and a measurement made on one site is not a fact about
  another. The ordering in section 5 is by rough cost class, not by a number.

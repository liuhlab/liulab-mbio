---
search:
  exclude: true
---

# Route B's per-well numbers: LevSeq's index PCR

Research note for issue #298, under #225. It fills the one deliberate omission of
`docs/research/bench-numbers.md`, which inventoried every other stage and recorded Route B's
per-well PCR as never read. Everything below was retrieved on **2026-10-06** unless a row
carries its own date.

Sections 2 to 6 answer the five questions #298 asks, in its order.

## How to read this note

Two documents carry everything: the article and its Supporting Information. The protocol lives
in the SI; the article's Materials and Methods restates parts of it in a paragraph. Five
disagreements turned up: two between the article and the SI, two inside the SI, and one between
LevSeq and the ONT protocol `bench-numbers.md` already quotes. Where sources differ, both
readings are given and neither is picked — section 7 collects them.

| Route | Works | Used for |
| --- | --- | --- |
| The article PDF and `sb4c00625_si_001.pdf`, read as `pdftotext -layout` dumps | yes | every protocol number below |
| `api.crossref.org/works/<doi>` | yes (`curl`) | the version of record's date, pages and licence policy |
| `ebi.ac.uk/europepmc/.../search?query=DOI:<doi>` | yes (`curl`) | open-access status |
| `api.figshare.com/v2/articles/28090968` | yes (`curl`) | the SI file's own licence |
| `api.biorxiv.org/details/biorxiv/<doi>` | yes (`curl`) | the preprint's licence |
| `pubs.acs.org/doi/<doi>` | **no**, HTTP 403 | nothing — the licence came from the three APIs instead |
| `neb.com/en/-/media/catalog/specifications/<letter>/<digit>/<catalog>_v1.pdf` | partly | B9004S came through; no spelling of M0267 did |

A value in a table is quoted from the cited document. A row marked **derived** is arithmetic on
quoted values, and the arithmetic is shown. Nothing here is from memory, and no number was
rounded into a cleaner one.

## 1. Sources and licences

| Document or dataset | Licence | Can we ship it? |
| --- | --- | --- |
| Long et al. 2025, *ACS Synth. Biol.* 14(1), 230–238, version of record | © ACS. Europe PMC reports `isOpenAccess: N`; Crossref lists only the ACS article-sharing policies (`doi:10.15223/policy-029`, `-037`, `-045`) | **No** — cite, never mirror |
| Supporting Information `sb4c00625_si_001.pdf`, deposited on ACS figshare as item 28090968 | **CC BY-NC 4.0**, from the figshare API on 2026-10-06 | **Not as a file or a table.** Non-commercial clashes with this repository's MIT licence, the same verdict Potapov 2018's data got. Single facts are cited and re-entered by hand |
| The bioRxiv preprint, `doi:10.1101/2024.09.04.611255` v1, posted 2024-09-04 | `cc_by_nc` from the bioRxiv API | **No**, same reason |
| The 24 nt barcode sequences (SI Tables S1–S4) | taken by the authors from the ONT native barcoding kit; ONT's own terms were not read here | **No** — and nothing in this note needs them. They stay user-held data, as #260 already decided for Route A's kit |
| `fhalab/LevSeq` software | GPL-3.0 (GitHub licence API, 2026-10-06) | **No** — a copyleft licence this repository cannot take code from |
| NEB product specification PS-B9004S v1.0, effective 10 Aug 2016 | © NEB, all rights reserved | **No** — cite, and re-enter single facts by hand, as `docs/research/restriction-ligation.md` §1 already set |

The practical verdict: **every number below may be cited and re-entered by hand; no document,
table or sequence set from these sources may be shipped as package data.** That is the same line
the repository already draws for NEB and for Addgene, and it is why the barcode plates stay a
thing the user holds rather than a file we ship.

Numbers themselves are facts, not expression. What the non-commercial terms stop is mirroring
the documents and lifting their tables wholesale, which this note does not do.

## 2. The per-well PCR

LevSeq's marking step is one colony PCR per well, straight from overnight culture — no miniprep,
no lysate step, no template purification. Article, "Colony PCR for generating barcoded
amplicons"; SI, "LevSeq Library Preparation and Sequencing" steps 1–5.

### What goes into one well

| Component | Volume | Source |
| --- | --- | --- |
| PCR master mix | 7 µL | SI step 2; article Methods |
| Barcode-linked primer mix | 2 µL | SI step 3; article Methods |
| Overnight culture | 1 µL | SI step 4; article Methods |
| **Reaction** | **10 µL** | **derived**: 7 + 2 + 1 |

The plate is a half-skirted 96-well PCR plate (USA Scientific 1402-9700, SI step 2). The article
adds that either a 96-well or a 384-well thermocycler block may be used; everything downstream
is still laid out and pooled as 96 wells.

### The master mix, as published

Stated per plate, not per well (SI step 1):

| Component | Amount per plate (µL) |
| --- | --- |
| ThermoPol Buffer (NEB B9004L) | 144 |
| 10 mM dNTPs (NEB N0447) | 28.8 |
| Taq Polymerase (NEB M0267) | 7.2 |
| Mol-Bio Grade DMSO (mp 194819) | 57.6 |
| ddH2O | 770 |

**Derived.** The five rows total 1007.6 µL, and 96 wells at 7 µL consume 672 µL, so the table is
1.5 times what one full plate needs (672 × 1.5 = 1008). The SI never says what the surplus is
for, and never says how to scale a part-plate run.

### What one 10 µL reaction holds

Every row is **derived**, by taking each component's share of the 1007.6 µL mix into 7 µL:

| Component | Per reaction | Final in 10 µL |
| --- | --- | --- |
| ThermoPol Buffer | 1.0 µL | 1x. NEB PS-B9004S v1.0 gives B9004 as a 10X concentrate: 20 mM Tris-HCl, 10 mM (NH4)2SO4, 10 mM KCl, 2 mM MgSO4, 0.1% Triton X-100, pH 8.8 at 25 °C |
| dNTPs | 0.2 µL of 10 mM each | 0.2 mM each |
| Taq | 0.05 µL | units per reaction unknown — see [Open gaps](#8-open-gaps) |
| DMSO | 0.4 µL | 4% (v/v) |
| ddH2O | 5.35 µL | — |
| Each primer | 2 µL at 1 µM | 0.2 µM, if the 1 µM plate is the one used — see §7 |
| Template | 1 µL overnight culture | — |

The arithmetic is exact: the mix is proportioned so that 7 µL of it carries round numbers.

### The cycling

Quoted from the SI's thermal-cycler table, with its two loop lines spelled out:

| Step | Temperature | Time |
| --- | --- | --- |
| 1 | 95 °C | 5 min |
| 2 | 95 °C | 20 s |
| 3 | touchdown 68 → 63.5 °C | 20 s |
| 4 | 68 °C | 1 min per kb |
| 5 | return to 2, 9 times | — |
| 6 | 95 °C | 20 s |
| 7 | 68 °C | 1 min per kb |
| 8 | return to 6, 24 times | — |
| 9 | 68 °C | 5 min |
| 10 | 4 °C | hold |

The touchdown drops 0.5 °C per cycle, 68 °C to 63.5 °C (SI note 5a). **Derived:** that is ten
values, which matches steps 2–4 running ten times in all, and the second block runs 25 times, so
the program is 35 cycles.

Extension is "minimally 1 minute per kb of gene", and the authors used 2 minutes for genes below
1 kb while developing the protocol (SI note 5b). The article's Methods restates the 300 s initial
denaturation and the 1 min/kb rule and refers the rest to the SI.

### The polymerase

Taq, NEB M0267. The article states the choice as a result: "PCR protocols are optimized for
robust amplification of the full-length gene. Best performance is obtained using Taq polymerase
and a touchdown PCR program." No proofreading enzyme is offered as an alternative, and no
comparison is reported.

### The primers, and why ours cannot be theirs

Each primer is a 24 nt barcode followed by a backbone-specific annealing region; the forward
barcode names the well and the reverse names the plate (SI, "Primer Design"). The published
layout is specific to pET-22b(+): the 3' region "binds to the cloning vector near the T7 and T7'
promotor sites", two universal sites upstream and downstream of the cloning site (SI, "Primer
Design"; article Methods). The SI says plainly that a user on a different cloning vector "must verify that the primers in the 'Primer
Design' section would still amplify the desired region", and recommends testing a new primer set
on a few wells before running it at plate scale.

For our method that means the sequences in Tables S3 and S4 are not usable as published: they
anneal to pET-22b(+), not to the working vector. What carries over is the layout — barcode at the
5' end, vector-specific annealing region at the 3' end — and every reaction number in this
section, which does not depend on which vector is amplified.

## 3. Pooling the indexed product

Two poolings, by two different rules.

**Within a plate, by equal volume.** 5 µL of each reaction into one tube, "for a total of 480 µL
of pooled PCR products per plate", one tube per plate (SI step 7). Nothing is normalised at this
stage: a well that amplified badly contributes a small amount of product in the same 5 µL.

The SI gives a 12-channel route to the same tube: 5 µL from each well of a column into one
destination well gives 40 µL per well across a row of 12, then 30 µL from each of those 12 wells
gives 360 µL in the tube (SI step 7a). **Derived:** 12 × 40 = 480 µL and 12 × 30 = 360 µL, so the
second route carries three quarters of the pooled volume forward and still holds every well. The
SI's own words for why it matters: combining 5 µL from every well "improve[s] the evenness" of
the pool. A liquid handler may do the same job (SI step 7b).

**Across plates, by equimolar amount, after the clean-up.** Each plate's gel-extracted product is
quantified and the plates are combined "in equimolar concentrations"; where every template is the
same length, equal mass is equal molarity, so equal nanograms per plate is what the SI instructs
(step 13). Its worked example: a plate at 20 ng/µL and a plate at 10 ng/µL give 5 µL + 10 µL +
5 µL water = 20 µL at 10 ng/µL.

| Number | Value | Source |
| --- | --- | --- |
| Volume per well into the plate pool | 5 µL | SI step 7 |
| Pooled volume per plate | 480 µL, or 360 µL by the 12-channel route | SI step 7, 7a |
| Ratio between plates | equimolar; equal mass where lengths match | SI step 13 |
| Enough material for any prep | 200 ng of pooled DNA per 1 kb of gene | SI step 13 |
| Library input | 200 fmol, "which translates to 120 ng for a 1-kb gene" | SI step 15 |
| Kit | ONT SQK-LSK114 on a MinION flow cell, FLO-MIN114 | SI steps 15–16; article Methods |
| What follows | ONT's own ligation protocol, linked and not restated | SI step 16 |

LevSeq stops at one tube of normalised amplicon and hands over to ONT's protocol, exactly as
Qian's route does. `docs/research/bench-numbers.md` already carries ONT's end-prep, ligation and
loading numbers from the published amplicon protocol, and they apply here unchanged with the
barcode-ligation rows dropped — Route B's barcode is in the primer, so no barcoding kit is used.

## 4. The clean-up

**There is no bead clean-up in Route B.** The clean-up is a gel: one lane per plate, cut and
extracted, which is also the quality check.

| Number | Value | Source |
| --- | --- | --- |
| Gel | 1% agarose with SYBR Gold (Thermo Fisher S11494), cast while the PCR runs | SI steps 6, 9 |
| Loaded | 100 µL of the pooled product + 20 µL 6x loading dye (NEB B7025S) | SI step 8 |
| The rest of the pool | stored at −20 °C for future use | SI step 8 |
| Ladder | 10 µL of 1 kb ladder in the flanking lanes | SI step 9 |
| Band wanted | the gene plus about 100 bp — barcodes, plus the stretch between the primer sites and the open reading frame | SI step 10 |
| Extraction kit | Zymoclean Gel DNA Recovery Kit, Zymo D4001 (SI) or D4002 (article Methods) | SI step 11; article Methods |
| Elution | 10 µL of ddH2O, explicitly **not** elution buffer | SI step 11 |
| Quantification | ng/µL on a GE NanoVue Plus, per plate | SI step 12 |

**What is measured about it: nothing.** No recovery yield, no purity figure, no comparison with a
bead clean-up, and no statement of how much of the 480 µL pool survives to the 120 ng that is
loaded. The gel's job is stated as a check — "evaluate … PCR products which should be the size of
the gene plus about 100 bp" — and the paper's only quantified consequence of skipping that check
is in the next paragraph.

**What is measured nearby, and is about the PCR rather than the clean-up.** On a ten-plate error-
prone library run "performed without optimizing PCR procedure", three plates returned under 60%
of variants sequenced and the other seven over 70% (SI Figure S3). The article attributes that
shortfall to "improper PCR amplification" and says checking the gel before sequencing mitigates
it. In the Sanger comparison, one 96-well plate read with two different barcode plates gave 94 of
94 samples on one plate and 90 on the other, with four wells below 10 reads (article; SI
Figure S2).

## 5. The limit on wells per run

**9,216 wells, and the barcode set is what sets it.** 96 forward barcodes address the wells of a
plate and 96 reverse barcodes address the plates, so one flow cell can tell apart 96 × 96 wells.
The article states it as "theoretically possible to demultiplex and sequence 9,216 variants"; the
SI states it as a limit, "the throughput is limited to 9,216 samples".

| Limit | Value | What sets it | Source |
| --- | --- | --- | --- |
| Wells per run | 9,216 | the barcode set: 96 forward x 96 reverse | article; SI, "Alternative Methods" |
| Barcode plates actually built | 8 (LevSeq01–08), so 768 wells before more are made | what the authors ordered | SI, "Preparation of LevSeq Barcode Primer Mixes"; Tables S5–S12 |
| Flongle flow cell | up to 1,600 variants recommended | flow cell output | article Methods |
| MinION flow cell | a theoretical 2,500 96-well plates at 1,000 bp | flow cell output | article Discussion |
| When to stop the run | 0.02 Gb of basecalled bases per plate; 20 plates means 0.4 Gb | enough depth per well | SI step 17; article Methods |
| Flow cell reuse | washed and reused until pores fall below 100; skip the storage buffer, which cost about 300 pores | pore count | article Methods |

**Derived:** 2,500 plates is 240,000 wells, which is 26 times the 9,216 the barcodes can name. On
a MinION the binding limit is therefore the barcode set, not sequencing output. On a Flongle the
flow cell binds first, at 1,600 wells — about 17 plates.

## 6. Where the depth floor is measured

**Per well, after demultiplexing, and not per pool.** The unit is reads assigned to a well, and
each read covers the whole amplicon, so there is no per-base depth and no pool-level figure
anywhere in either document.

| Statement | Value | Source |
| --- | --- | --- |
| Quality checklist, per well | "Alignment Count: Must be >20 for each well" | SI, "Initial Data Quality Assessment" |
| Second per-well criterion | mean error < 10% for each well | same |
| Recommended minimum depth | 20 reads, "as strand bias may impact calling accuracy" | article |
| Warning threshold in the software | fewer than 20 reads assigned to a well | article |
| Simulated accuracy | variants detected at >99% with 10 reads, up to a 20% sequencing error rate | article; SI Figure S1 |
| What "low coverage" meant in practice | fewer than 10 reads; four wells of 94 on one barcode plate | article; SI Figure S2 |
| Mixed-well detection | at 20 reads the false-negative rate is high; below 10 reads true positives are not detected at all | SI Figure S1E |

So #260's reading — twenty wanted, ten tolerable — is what the two documents support, and both
numbers are per well.

**Against what read length: none is stated, and the simulation says length does not matter.**
Varying the simulated sequence length showed "no trend … so gene length does not affect our
software's variant calling ability" (SI Figure S1C). The only place a length enters is the
sequencing budget, which assumes 1,000 bp. The amplicon itself is the gene plus about 100 bp.

**Derived**, and worth knowing because it says how much headroom the stop rule leaves: 0.02 Gb
over 96 wells is about 208 kb per well, which at a 1.1 kb amplicon is about 190 reads per well —
roughly ten times the 20-read floor. The gap is what covers uneven amplification across wells.

## 7. Where the sources disagree

Each of these is a real disagreement in the published material, not a reading difficulty.

1. **Primer concentration, by a factor of ten.** The SI's plate preparation ends by diluting the
   1 µM stock tenfold into "1X reaction-ready barcode plates … 0.1 µM for both forward and
   reverse" (step 5). But the library protocol stamps "2 µL of the 1-µM barcode-linked primer
   mix" (step 3), and the article's Methods agrees: "2 µL of 1 µM each barcoded primer mix". Two
   of three places say 1 µM, giving 0.2 µM in the reaction, which is the conventional figure and
   the one NEB's own ThermoPol specification assays at. The 0.1 µM plate would give 0.02 µM. This
   note records 0.2 µM as what two sources state and flags the third.
2. **The boundary at exactly 20 reads.** The SI's checklist requires an alignment count **above**
   20 and routes a well at exactly 20 into the suboptimal branch. The article says a warning is
   issued when **fewer than** 20 reads are assigned, which passes a well at exactly 20.
   `liulab_synbio.dmx.ROUTE_B` currently takes `inclusive=True`, matching the article. One well in
   many will land exactly on the line; the two documents judge it differently.
3. **200 fmol in nanograms.** LevSeq says 200 fmol "translates to 120 ng for a 1-kb gene" (SI
   step 15). ONT's own amplicon protocol, quoted in `docs/research/bench-numbers.md`, gives the
   same 200 fmol as 130 ng for a 1 kb amplicon. Both are conversions of the same molar figure, so
   the molar number is the one to carry.
4. **The gel kit's catalogue number.** Zymo D4001 in the SI, D4002 in the article's Methods —
   the same Zymoclean kit in two pack sizes. Which was used is not stated.
5. **The SI's own step cross-references slip by one** near the gel: its step 9 says "each tube
   made in step 9" where it means step 8, and "the agarose gel prepared in step 7" where it means
   step 6. Anyone following the printed order hits this.

Nothing here contradicts `docs/research/bench-numbers.md`, which recorded this stage as unread.
Point 3 is the only place the two notes meet, and they meet on the same 200 fmol.

## 8. Open gaps

- **Units of Taq per reaction.** The SI gives 7.2 µL per plate, which is 0.05 µL per reaction,
  and never gives the enzyme's concentration. NEB's specification for M0267 could not be read:
  `pubs`-style blocks aside, every spelling tried under
  `neb.com/en/-/media/catalog/specifications/m/0/` (`m0267s_l_v1.pdf`, `m0267s_v1.pdf`,
  `m0267l_v1.pdf`, `m0267x_v1.pdf`, `m0267s_l_x_v1.pdf`, `m0267_v1.pdf`, `m0267s_l_v2.pdf`)
  returned the same 108,707-byte generic page, and `r.jina.ai` returned an empty body for the
  product page. The route itself works — `b/9/b9004s_v1.pdf` came through it on the same day — so
  this is a filename that was not found, not a blocked route. A protocol can print the volume
  without the unit count; a unit count must not be guessed.
- **Gel-extraction recovery.** Not measured, so how much of a plate's pool reaches the flow cell
  is unknown.
- **Scaling a part-plate run.** The master mix is 1.5 times what 96 wells need and the SI does not
  say why, nor what to do for fewer wells.
- **Thermocycler detail.** No ramp rate, no lid temperature, no block requirement beyond "96-well
  or 384-well".
- **Storage and stability.** "Store at −20 °C for future use" is all that is said about the
  leftover pool; no hold time is given, for it or for the gel-extracted product.
- **Plate format other than 96.** A 384-well block is named for cycling only. Every barcode,
  pooling and plate-map number in both documents is 96-well, so Route B at 384 wells per plate is
  not published.
- **Primers for our vector.** None exist. LevSeq's are pET-22b(+)-specific and the SI hands the
  redesign to the user. What a new pair must satisfy is a primer-design question, not a number
  this note can supply.
- **Whether the 24 nt barcodes may be re-ordered as oligos.** They come from the ONT native
  barcoding kit; ONT's terms were not read for this note.
- **The mean-error criterion.** "Mean error < 10% for each well" has no definition in either
  document beyond what the software computes.
- **The preprint's SI was not compared with the published SI.** Only the version-of-record SI was
  read. If a number ever has to be redistributed rather than cited, the CC BY-NC preprint is where
  to check first.

## What this note supplies

The reaction Route B has been missing: a 10 µL colony PCR of 7 µL master mix, 2 µL barcoded
primer mix and 1 µL overnight culture, cycled 35 times with a ten-cycle touchdown from 68 °C to
63.5 °C and 1 min/kb extension; pooled 5 µL a well into one tube a plate; gel-purified per plate;
combined between plates by equal mass; and read to 20 reads a well. Every value above carries its
citation, and the three that are disagreements carry both readings.

It supplies no primer sequences, and it changes nothing in `src/`.

## Sources

- Long, Y.; Mora, A.; Li, F.-Z.; Gürsoy, E.; Johnston, K. E.; Arnold, F. H. LevSeq: Rapid
  Generation of Sequence-Function Data for Directed Evolution and Machine Learning. *ACS Synth.
  Biol.* **2025**, *14* (1), 230–238.
  [doi:10.1021/acssynbio.4c00625](https://doi.org/10.1021/acssynbio.4c00625). Published online
  2024-12-24; read 2026-10-06.
- Long, Y. et al. Supporting Information for the above, `sb4c00625_si_001.pdf`, ACS figshare item
  28090968, CC BY-NC 4.0, deposited 2024-12-24; read 2026-10-06. Sections "Oligonucleotide
  Design", "Supplementary Protocols" and "LevSeq Post Data Analysis Workflow and Interpretation";
  Figures S1–S3.
- Long, Y. et al., preprint, *bioRxiv* 2024.09.04.611255 v1, posted 2024-09-04, CC BY-NC.
  [doi:10.1101/2024.09.04.611255](https://doi.org/10.1101/2024.09.04.611255). Licence read from
  the bioRxiv API 2026-10-06; the text was not used.
- New England Biolabs, *Product Specification: ThermoPol Reaction Buffer Pack*, B9004S,
  PS-B9004S v1.0, effective 10 Aug 2016. Read 2026-10-06 from
  `neb.com/en/-/media/catalog/specifications/b/9/b9004s_v1.pdf`.
- `docs/research/bench-numbers.md` (issue #233) for every other stage, for the ONT library
  numbers this stage hands over to, and for the omission this note fills.
- `docs/research/synthesis-and-assembly-materials.md` (issue #225) for the plate maps and primer
  plates already inventoried.

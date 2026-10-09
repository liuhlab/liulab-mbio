---
search:
  exclude: true
---

# Choosing a read-back route, and what a run that reads nothing back holds

Research note for issue #473, under map #462. The AP-1 demo offers barcode ligation and index
PCR as two labelled alternatives and has to tell a reader how to pick one, or neither. Three
claims had to be sourced before that guidance could be written. This note states each claim,
what the published sources actually say, and — where they say nothing — says so plainly.

Everything below was read on **2026-10-08** from the two method papers already held under
`reference_docs/synthesis_and_assembly/`. No new document was fetched.

## How to read this note

| Marking route | Source | Note that already holds its bench numbers |
| --- | --- | --- |
| barcode ligation | Qian et al. 2026 (article, SI) | `docs/research/bench-numbers.md` Stage 3 |
| index PCR | Long et al. 2025, LevSeq (article, SI) | `docs/research/route-b-index-pcr.md` |

A quoted row is verbatim from the cited document. A row marked **derived** is arithmetic on
quoted values, with the arithmetic shown. A row marked **no source** means the documents were
searched and state nothing; it is not an invitation to supply a figure.

**One correction to #473's premise.** The issue cites `route-b-index-pcr.md:501-504` for the two
capacity figures. That file is 422 lines; the text meant is `route-b-index-pcr.md` section 5 and
`docs/research/synthesis-and-assembly.md` lines 502-503. The capacities themselves are as the
issue states, and section 1 below confirms that neither of them decides anything.

## 1. Capacity does not pick a route, and both notes already said so

| Route | Addressable wells | Source |
| --- | --- | --- |
| barcode ligation | ">330k wells to be tagged from a single 96-well plate of UMIs (24^4)" | Qian article |
| barcode ligation, as coded | 24³ wells × 24 plates = 13,824 per plate set | `synbio.dmx.method.ROUTE_LIGATION` |
| index PCR | "theoretically possible to demultiplex and sequence 9,216 variants" | LevSeq article |

AP-1 needs 1,536. Both routes clear it by a wide margin, and the smaller of the two is index
PCR's. So capacity ranks the routes the *wrong way* for the lab's claim, and cannot be its basis.

`docs/research/synthesis-and-assembly.md` reached the same conclusion under #299: "nothing picks
between them, so a build names its route … what binds either is the flow cell and one reaction
per well." That last clause turns out to be the live half, and section 2 is what it rests on.

## 2. The claim: index PCR for few samples, barcode ligation for many

**Verdict: the claim has a sourced basis, but not the one the capacity figures suggest, and the
evidence for the two halves is of very different quality.** Three separate things support it.

### 2.1 One thermocycled reaction a well — the only direct comparison between the routes

Qian's introduction names index PCR as the limitation the ligation route was built to remove:

> GOIs are barcoded. This step is typically achieved by dual-indexed PCR, which limits
> throughput and requires large numbers of thermocyclers when processing thousands of samples.
> To overcome this limitation, we developed an isothermal barcoding strategy based on GGA that
> ligates barcodes with the GOI directly in bacterial lysate.

— Qian article. The same paper states the mechanism again in its summary of what makes the
pipeline scale: "a combinatorial in-lysate isothermal barcoding framework, which eliminates the
need for thermocyclers".

**This is the whole of the direct evidence**, and it is worth being exact about three things.

1. It is an instrument-count argument, not a cost or capacity one. Index PCR cycles; barcode
   ligation sits at 37 °C for an hour and then 60 °C for five minutes, so a whole 1536-well
   plate marks in one incubator.
2. It names no number of samples. "Thousands" is as precise as the sentence gets.
3. It is the ligation route's own authors describing the alternative. LevSeq does not answer it.

The per-well volumes, both already recorded, show the same split from the other side: barcode
ligation is **3.5 µL** a well in a 1536-well plate moved by an acoustic handler (Qian SI Day
4.1, `bench-numbers.md`), index PCR is **10 µL** a well in a 96-well plate (LevSeq SI,
`route-b-index-pcr.md` section 2). **Derived:** 1,536 wells is 5.4 mL of ligation reaction
against 15.4 mL of PCR, and 1 plate against 16.

### 2.2 A flat per-run cost that only amortises over a large library

Qian's Supplementary Table 5 is the one real cost table either paper publishes. Costs are the
table's, in US dollars, for the whole library:

| Library complexity | 100 | 500 | 1,000 | 2,000 |
| --- | --- | --- | --- | --- |
| Twist oligo pool, 250-300 nt, list price | 1,030 | 2,060 | 3,090 | 4,121 |
| Cloning (BsaI + Salt T4 ligase) | 191.46 | 382.92 | 765.85 | 1,531.70 |
| Sequencing (MinION flow cell $600 + ONT library prep kit $95) | 695 | 695 | 695 | 695 |
| Plasticware | 35.35 | 101.32 | 202.64 | 405.27 |
| **Total** | **1,951.81** | **3,239.24** | **4,753.49** | **6,752.97** |
| **Per design** | **19.52** | **6.48** | **4.75** | **3.38** |
| The arrayed-fragment alternative, per design | 22.24 | 22.24 | 22.24 | 22.24 |

**The sequencing row is flat.** $695 is the same whether it reads 100 designs or 2,000, which is
why cost per design falls 5.8-fold across the table while nothing about the chemistry changes.
At 100 designs the route costs $19.52 against the alternative's $22.24 — a saving of 1.14x, which
is nothing. At 2,000 it is 6.59x.

**That is the sourced bound on the ligation route: a fixed per-run cost, paid once, that is only
worth paying when many designs share it.** It is a cost-amortisation argument, and it is the
strongest evidence either paper offers for anything in the lab's claim.

Two cautions on the table:

- **The paper's own text disagrees with it.** The article says the approach "is five and eight
  times more cost-effective … than ordering 500 and 2000 individual gene fragments"; the table
  gives 3.43x and 6.59x at those two sizes. The paper does not reconcile them. Cite the table.
- **The comparison is not against index PCR.** The arrayed-fragment column is Qian's own SAPP
  route — eBlocks, no clone isolation. No cost comparison between the two read-back routes
  exists in any source read here.

### 2.3 Each paper's own "ideally suited" sentence, and what each one is actually about

| Route | What the paper says it suits | Source |
| --- | --- | --- |
| barcode ligation | "ideally suited for large design campaigns where the time overhead is absorbed over many constructs, while the polyclonal approach from SAPP is better suited for smaller campaigns (~100-200 designs) where faster turnaround time is desired" | Qian article |
| barcode ligation | optimum "about 2000 genes, requiring the sampling of about 8000 colonies (20 × 384 well plates) to achieve >95% recovery"; "For smaller campaigns, libraries of 500 genes can be demultiplexed by manually picking 1500 colonies" | Qian article |
| index PCR | "ideally suited for small to medium-sized protein engineering laboratories" | LevSeq SI |
| index PCR | "the throughput is limited to 9,216 samples" | LevSeq SI |
| index PCR | "For industrial scale protein engineering procedures or in a sequencing facility, ParSEQ may be a more suitable pipeline" | LevSeq article |

**Read these carefully, because they are not symmetric.**

- Qian's sentence does say *large campaigns*, and gives the mechanism: a **five-day** fixed time
  overhead ("The DMX pipeline is a 5-day protocol"; "only adding five days to the pipeline")
  that, like the $695, is paid once and divided among the designs. But its contrast is with
  **SAPP's polyclonal route**, not with index PCR. Quoting it as a comparison between the two
  read-back routes would misattribute it.
- LevSeq's sentence is about **laboratory size**, not sample count — a lab "without centralized
  sequencing facilities". It does not say index PCR suits few samples. Its only sample-count
  statement is an upper limit of 9,216 and a pointer to parSEQ ("up to 36,864 variants per run",
  LevSeq SI) above it.

### 2.4 What the demo may and may not write

**May be written, with these citations:**

- Barcode ligation carries a fixed cost and a fixed delay — $695 of sequencing and five days —
  that fall per design as the library grows, from $19.52 a design at 100 to $3.38 at 2,000
  (Qian, Supplementary Table 5; Qian article).
- Index PCR needs one thermocycled reaction per well, which at thousands of samples means many
  thermocyclers; barcode ligation is isothermal and marks a whole 1536-well plate in one
  incubator (Qian article).
- Index PCR's published ceiling is 9,216 samples, and its authors point past it to parSEQ
  (LevSeq SI; LevSeq article).

**Has no source and must not be written:**

- **Any cost figure for index PCR.** LevSeq states no per-sample, per-plate, per-well,
  primer-plate or sequencing price, and no dollar comparison with Sanger. Its only cost words
  are qualitative — "a more economical option compared to Sanger", "cost-efficient alternative".
  The dollar and euro figures in its SI (parSEQ "under $3 per variant", Ramirez Rojas "£2.20 per
  sample", DuBA.flow "under €0.10 per sample", BPS "as low as $0.12 per plasmid") belong to
  **other methods** and must never be attached to index PCR.
- **A cost per sample for either route at AP-1's 1,536 wells.** Qian's table stops at 2,000
  designs and is priced per *design*, not per well; AP-1 reads 72 designs across 1,536 wells.
  Interpolating it would be a guess.
- **Hands-on hours for either route.** Qian marks hands-on versus automated steps in a figure
  only, with no hour counts in the text. LevSeq's "3-12 hours" turnaround is a citation to
  evSeq, not a LevSeq measurement, and neither paper states person-time, pipetting-step counts
  or samples per person-day.
- **A one-time price for the barcode plasmid kit or the index primer plate.** The kit is 96
  Addgene plasmids (IDs 255161-255256) and the primer set is 192 IDT oligos ordered in plate
  format at 100 µM; neither paper prices either, and `route-b-index-primers.md` section 3.3
  already recorded that no vendor page was reachable for the oligos.
- **A crossover point.** No source states a library size at which one route overtakes the other.

**A caveat the guidance should carry.** Every direct comparison above is Qian's, written by the
ligation route's authors about the route they replaced. No independent or head-to-head
measurement of the two exists in the material held here.

## 3. What a run that reads nothing back holds

`CONTEXT.md`'s **Validation floor** says that with no floor "nothing is read and the cargo stays
polyclonal". What a reader may assume about that material is narrower than it looks.

### 3.1 The one measurement, and what it was measured on

> While polyclonality due to synthesis or cloning errors is a risk with this approach,
> sequencing analysis of 929 GGA reactions with insert size below 314 bp showed that 89.3% had a
> clonal purity >= 90%.

— Qian article, Fig. 1b. Method: stabs from 929 glycerol stocks pooled, amplified with Q5, read
on an Illumina NextSeq2000, 31 million aligned reads; "clonal purity scores were calculated by
dividing the number of perfect matches by the total number of filtered reads for each reference
sequence".

| What it measures | 89.3% of wells reached >=90% clonal purity |
| --- | --- |
| The material | **SAPP**: eBlock gene fragments, Golden Gate assembled, no clone isolation |
| The insert | **below 314 bp** |
| The platform | Illumina NextSeq2000, short-read |
| What it is not | not an oligo-pool library, not a DMX product, and not a per-well guarantee |

A second measurement on the same material: expected mass within ±1 Da for 90% of expressed
products, 777 of 863 (Qian article).

### 3.2 What a reader may assume

**May assume:** an unread well is *enriched* for the intended design, and for short inserts from
arrayed fragments the enrichment is high — 89.3% of wells at >=90% purity below 314 bp.

**May not assume:**

- **That the figure transfers to a longer insert.** The paper forecloses it: "the small levels
  of polyclonality resulting from synthesis error are length-dependent, which can complicate
  using the protocol for longer proteins."
- **That it transfers to oligo-pool cargo.** It was measured on eBlocks. AP-1's cargo comes from
  an oligo pool, which carries its own synthesis error and abundance bias, and **no measured
  purity figure exists for an unread oligo-pool well in any source read here.**
- **That the material is uniform.** **No source.** Qian's only statement about pool evenness is a
  *simulation input* — Supplementary Fig. 5 assumes "the oligo abundance bias … to be normal and
  approximating the distribution reported by Twist", with a synthesis error rate of 1/3000 nt.
  Both are assumptions fed into a model, not measurements of this material, and neither may be
  quoted as one.
- **That it is mutation-free.** Nothing claims this. The 10.7% of SAPP wells below 90% purity,
  and the length dependence above, are what the sources say instead.
- **That Qian's 78% and 92% describe unread material.** They do not: "We successfully recovered
  78% of the library variants as clonal, sequence-verified constructs by screening 4608 clones
  (3x the library complexity), approaching the theoretical maximum recovery rate of 92%". Those
  are recovery rates *after* reading back, and the 92% is a simulated ceiling.

**So the demo's claim — enriched for a design, but not uniform and not mutation-free — is sound
in direction and only partly sourced in substance.** "Not mutation-free" has one anchor, on
material unlike AP-1's. "Not uniform" has none, and should be written as a caution rather than
as a measured fact.

## 4. Where each depth floor was measured

`synbio.dmx.method` states that "each floor was measured on its own library prep". That
is true, but it is not where the difference lies — **the ONT kit and flow cell are the same in
both papers**, and reading the docstring as "different platforms, so don't compare" would be
wrong for the wrong reason.

| | barcode ligation | index PCR |
| --- | --- | --- |
| Floor | consensus called where read depth **>150** | alignment count **>20** per well, 10 tolerated |
| Where stated | Qian SI, "Generation of consensus sequences" | LevSeq SI, "Initial Data Quality Assessment" |
| What the number is | a **consensus-calling parameter** in the pipeline: a consensus is generated for a well only above it, with a base included at >=51% read support | a **QC checklist criterion** a human applies per well, paired with a second one — mean error <10% |
| Sequencing kit | ONT Ligation Sequencing Kit V14, SQK-LSK114 | ONT LSK114, "used while developing this protocol" |
| Flow cell | MinION R10.4.1, FLO-MIN114 | MinION, FLO-MIN114 |
| Basecalling | Dorado 0.9.1, super-accuracy | super-accuracy model recommended |
| Read filtering before counting | chopper, >Q15 and 400-1000 bp; cutadapt UMI match of >=20 of 25 bp | not stated as a length or quality filter |
| How the amplicon was made | barcodes ligated in lysate, then the **pooled** plate amplified in three 25 µL KAPA HiFi reactions of **10 cycles**, gel-extracted | **one 10 µL Taq colony PCR per well**, 35 cycles with a touchdown, 5 µL a well pooled, gel-purified per plate |
| The sample | 4,608 wells of a 1,500-design Twist oligo-pool library, inserts 275-400 bp | a ~200-aa gene in pET-22b, amplicon the gene plus ~100 bp |
| Evidence class | stated as a pipeline parameter | 10 reads from **simulation** (2% epPCR rate, 10% nanopore error); 20 "recommended" as a hedge, "as strand bias may impact calling accuracy" |

### 4.1 Why the two numbers are not comparable, in one line each

1. **They are different kinds of number.** 150 is where software starts calling; 20 is where a
   person stops trusting. A pipeline parameter and a QC threshold are not on one scale.
2. **The amplicons have different error structures.** 10 cycles of a high-fidelity polymerase
   across a pool is not 35 cycles of Taq in each well, and the depth needed to out-vote
   amplification error follows from that.
3. **One floor filters its reads first.** Qian counts only reads passing Q15 and a 400-1000 bp
   length window with a 20-of-25 bp UMI match. LevSeq counts alignments.
4. **One is simulated and one is not.** LevSeq's 10-read figure comes from simulated data;
   LevSeq itself warns "experimental variation may reduce this accuracy".

**So the guidance must not present 150 and 20 as a strict-versus-lenient pair.** The honest
statement is that each route is run at the depth its own published pipeline was run at.

### 4.2 Two things no source states

- **Why 150.** Qian gives no derivation, no curve and no sensitivity analysis for it. It appears
  as a parameter in the consensus step and nowhere else. **No source.**
- **The depth either run actually achieved.** Neither paper reports a mean or median reads per
  well as text. LevSeq's per-well coverage appears only in Figure 3B; Qian reports none.
  `route-b-index-pcr.md` section 6 derives ~190 reads a well for index PCR from the 0.02 Gb
  stop rule, which is arithmetic on a budget, not a measurement.

## What this note supplies

The sourced basis for three sentences of demo guidance, and an explicit list of the figures that
must not appear in them. It supplies **no new number to `src/`**: every value above is either
already carried by `synbio.dmx.method` and the two bench notes, or is a cost figure the
package does not compute.

The one thing it changes about existing text is a caution: `dmx/method.py`'s "each floor was
measured on its own library prep" is right in its conclusion and loose in its reason, since both
floors came off the same ONT kit and flow cell. Section 4 is the reason.

## Sources

- Qian, Z. et al. Accelerating protein design by scaling experimental characterization.
  *Nat. Commun.* **2026**.
  [doi:10.1038/s41467-026-76740-9](https://doi.org/10.1038/s41467-026-76740-9). Article and
  Supplementary Information, including Supplementary Table 5 and Supplementary Fig. 5, held
  under `reference_docs/synthesis_and_assembly/dmx/paper/`; read 2026-10-08.
- Long, Y.; Mora, A.; Li, F.-Z.; Gürsoy, E.; Johnston, K. E.; Arnold, F. H. LevSeq: Rapid
  Generation of Sequence-Function Data for Directed Evolution and Machine Learning.
  *ACS Synth. Biol.* **2025**, *14* (1), 230-238.
  [doi:10.1021/acssynbio.4c00625](https://doi.org/10.1021/acssynbio.4c00625). Article and
  Supporting Information `sb4c00625_si_001.pdf`, held under
  `reference_docs/synthesis_and_assembly/LevSeq/papers/`; read 2026-10-08.
- `docs/research/route-b-index-pcr.md` (#298) for index PCR's per-well reaction, its 9,216-well
  limit and its depth floor.
- `docs/research/bench-numbers.md` (#233) for the barcode-ligation stage, the lysate reaction
  and the ONT read-out.
- `docs/research/route-b-index-primers.md` (#312) for the 192-oligo primer set and the vendor
  prices that could not be retrieved.
- `docs/research/synthesis-and-assembly.md` (#299) for the earlier finding that capacity picks
  no route.

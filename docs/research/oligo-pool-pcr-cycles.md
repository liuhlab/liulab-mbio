---
search:
  exclude: true
---

# Cycles for the two oligo-pool PCRs

Research note for issue #340, under #225. It is hole 11 of the eleven #268's end-to-end run
found, and the one the repository treats as worst: a number that looks sourced and is not.
Everything below was retrieved on **2026-10-07** unless a row carries its own date.

## The contradiction

`igga/steps.py`'s `_pcr1_step` and `_pcr2_step` pass no `cycles`, so `pcr_program` falls through
to the polymerase's own profile and the rendered protocol prints **30** for both. The AP-1 demo
shows it: `docs/examples/ap1-library/protocol.json` carries `"cycles": 30` in the inner stage of
PCR1 and of PCR2. The same step's troubleshooting says *"Take the fewest cycles that give a
visible band."* A reader at the bench has a printed count and a rule that contradicts it.

## How to read this note

| Route | Works | Used for |
| --- | --- | --- |
| `twistbioscience.com/content/dam/...` direct PDF fetch | yes (`curl`) | both oligo-pool amplification guides |
| The four Twist PDFs already under `reference_docs/`, as `pdftotext -layout` dumps | yes | the Multiplexed Gene Fragments rule this note separates out |
| Lund 2024, Freschlin 2026, Romanowicz 2026 text dumps already held | yes | every published cycle count below |
| `idtdna.com` product pages | partly — a product announcement, no protocol document | one contrast in [section 8](#8-open-gaps), and no number |
| An Agilent oligo-pool amplification protocol | **not fetched** | nothing |

A value in a table is quoted from the cited document. A row marked **derived** is arithmetic on
quoted values, and the arithmetic is shown. Nothing here is from memory.

## 1. The two reactions are not one reaction

| | PCR1 | PCR2 |
| --- | --- | --- |
| Template | the synthesised pool as it arrives | PCR1's product |
| Primers | the batch's forward and outer, both shared by every oligo in the batch | the batch's forward, and the block's own inner |
| What is amplified | every oligo of the batch | the pieces of one block |
| Reactions | one per batch | one per block |
| In the AP-1 demo | 1 reaction, 153 oligos at 350 nt | 72 reactions |

The demo has 72 blocks and a `batch_size` of 96, so it makes **one** batch, and PCR1 is therefore
a whole-pool amplification with a single universal primer pair. That is precisely the reaction
every oligo-pool vendor protocol is written for. A project with more blocks than `batch_size`
splits PCR1 into several reactions, each pulling its own share out of the same tube.

PCR2 is a different animal. Its template is an amplified, abundant product, and it selects one
block within it.

## 2. Twist's own rule for an Oligo Pool

Twist publishes two amplification guides for two different products. This method orders an
**Oligo Pool** — 20–350 nt, single-stranded, delivered at femtomole levels
(`docs/research/synthesis-and-assembly-materials.md`). Both of its guides give the same table:

| Oligo pool length | Cycles |
| --- | --- |
| 20–100 nt | 6–10 |
| 100–150 nt | 10–12 |
| **151–350 nt** | **12–14** |

Source: Twist FRM-001034 REV 8 p. 2, and DOC-4060 REV 1.0 "General Notes and Precautions" and
step 1. Two independently revised documents, identical numbers.

Oligo length is a project decision. The AP-1 demo is 350 nt and the method's earlier default was
300 nt, so both land in the top band: **12–14 cycles** is the sourced count for this design's
PCR1.

**How the pool arrives, and where the 20 ng/µL comes from.** DOC-4060's "Before You Begin", the
same document as the band table:

> Twist Oligo Pools are delivered as a lyophilized product pooled in a single tube. The total
> yield in ng is printed on the shipping tube label.
>
> Prepare a stock solution of your Oligo Pool by resuspending in 10 mM Tris buffer, pH 8.0 to a
> concentration of at least 20 ng/µl. Stock solution concentration (ng/µl) = Total yield (ng) /
> resuspension volume (µl).

Its component table stores the pool at "-20°C for 24 months", 4 °C for 12, and -80 °C long term.

Three things follow. The 20 ng/µL the PCR1 row above takes as template is the **stock the reader
prepares**, not a concentration the vendor ships at — DOC-4060's own amplification table asks for
"Oligo pool (20 ng/µl), 20 ng, 1 µl", and "at least 20 ng/µl" is the separate rule for making
that stock. The volume is **derived** from the label and that floor: resuspension volume (µl) =
total yield (ng) / 20, rounded down, leaves the stock at the floor or above it, so the table's
1 µl still delivers its 20 ng. And the yield is printed on the tube and written in no file, so
the protocol states the division and the reader supplies the one number in front of them — a rule
they can execute, not a hole.

The rest of the reaction, where the two revisions differ:

| | FRM-001034 REV 8 | DOC-4060 REV 1.0 |
| --- | --- | --- |
| Polymerase | KAPA HiFi HotStart (Roche) | Twist TrueAmp Polymerase Mix |
| Reaction | 25 µL | 50 µL |
| Template | 0.5 µL of 20 ng/µL — **derived:** 10 ng, matching the text's "~10 ng" | 1 µL of 20 ng/µL = 20 ng |
| Each primer | 0.3 µM | 0.3 µM |
| Initial denaturation | 95 °C 3 min | 98 °C 45 s |
| Denature / anneal / extend | 98 °C 20 s / 15 s at optimum / 72 °C 15 s | 98 °C 15 s / 60 °C 15 s / 68 °C 45 s |
| Final extension | 72 °C 1 min | 68 °C 1 min |
| Clean-up | spin column, or 1.8x SPRI above 120 nt | beads; 1.8x for 120–300 nt, 1.2x for 300–500 nt |
| Uniformity claimed | >90% of oligos within <2x of the mean | >90% within <2.5x of the mean |

**A conflict inside REV 8.** Its cycling table says "6–12 Cycles\*\*" while its own footnote
sends the reader to the chart, which gives 12–14 for our band. DOC-4060 prints no count in the
cycling table and defers to the chart alone. Take the chart: it is what both documents agree on,
and what the FAQ is written against.

DOC-4060's Appendix B answers this ticket's question in Twist's own words:

> **What is the relationship between uniformity and the PCR cycle number?** More PCR cycles would
> lead to worse uniformity. Keeping PCR cycles in the ranges outlined in the General Notes and
> Precautions section is recommended.

Three more from the same appendix bear on the design:

- 20 ng input is "the optimized number to secure good uniformity and a low dropout rate"; "using
  1 amol of each oligo would still give good uniformity and low dropout rate."
- GC should stay between 35% and 65%: "The higher the GC%, the less efficient PCR amplification
  is in the same oligo pool ... to avoid bad uniformity and high dropout."
- Another polymerase will amplify a pool, but "there may be worse uniformity and less yield with
  the same PCR cycle numbers when comparing with alternative polymerases."

That last one lands on us. The pipeline uses **Q5**. Twist states 12–14 against KAPA HiFi HotStart
or its own TrueAmp mix, and says in the same breath that the count does not carry to another
enzyme unchanged.

**The observable for over-cycling**, identical in both guides, is not a smear: "The presence of a
hump after the peak of interest indicates heteroduplexes, a result of over-amplification," remedy
"Re-try PCR with lower number of cycles." It is read on capillary electrophoresis — a Bioanalyzer
DNA 1000 chip, or a TapeStation — not on an agarose gel. The step's present troubleshooting looks
for "a smear rather than a band" on a gel, which is a coarser and different observation.

## 3. The eight-cycle rule this corrects

`docs/research/synthesis-and-assembly.md`, under "Twist — inherited", records "the amplification
rule of a hard maximum of eight cycles with a high-fidelity hot-start polymerase."

**That rule is from the wrong product.** Eight is the Multiplexed Gene Fragments figure:
DOC-001498 REV 1.0 ("the number of amplification cycles should be minimized to fewer than 8") and
DOC-4057 REV 1.0 ("Adhere to a strict maximum of 8 amplification cycles").

The two products get opposite advice because they arrive in opposite states:

| | Multiplexed Gene Fragments | Oligo Pools |
| --- | --- | --- |
| What arrives | double-stranded, 301–500 bp, ≥200 ng | single-stranded, 20–350 nt, femtomole levels |
| Twist's advice | "does not recommend performing any PCR amplification ... before use"; amplifying voids the quality guarantee | "we recommend including primer-binding sites in your oligo design and performing PCR amplification of the pool prior to usage" |
| Cycles if you do | fewer than 8, hard maximum 8 | 6–14, by length band |

So the number that reads as Twist's rule for our pool is Twist's rule for a product we do not
order, and it is wrong in the conservative direction — which is why nothing has caught it.
`reference_docs/synthesis_and_assembly/twist/README.md` already warns against carrying a figure
between the two; this is that warning firing on a tracked note.

## 4. What the published methods do

| Source | The reaction | Template | Polymerase | Cycles |
| --- | --- | --- | --- | --- |
| Lund et al. 2024 | one amplification of a Twist oligo pool, 300 nt, by a design's own primer pair | 3.33 pg/µL | Q5 Hot Start HiFi 2X (NEB M0494) | **none stated** — run on a Bio-Rad CFX96 with LAMP Fluorescent Dye (NEB B1700) at 0.5x; "Cycles were stopped before all amplifications plateaued in fluorescence" |
| Romanowicz et al. 2026 | Twist oligo subpool amplification | not stated | KAPA HiFi | **none stated** — "qPCR used to determine cycle numbers and minimize overamplification bias"; a later enrichment is "stopped in the late exponential phase" |
| Freschlin et al. 2026 | subpool PCR off a Twist oligo pool of 300 bp oligos, orthogonal primer pairs "for selective amplification" | 1 ng | KAPA HiFi HotStart ReadyMix 2x (Roche KK2601) | **35** |
| Qian et al. 2026 | one whole-library amplification off the oligo pool | 2.5 ng | KAPA HiFi HotStart Ready Mix, with 20x EvaGreen | **14** |

Two readings are worth keeping.

**The two that state no number measure instead.** Lund and Romanowicz put a dye in the reaction
and stop on the curve. Qian puts EvaGreen in the reaction *and* prints 14 — the dye is there so
the operator can watch the same thing, whether or not the printed count is kept. Qian's 14 also
sits inside Twist's 12–14 band, which is the only independent agreement this note found.

**Freschlin's 35 is the one published count above ours, and it is the one with a measured
consequence.** Their assembled libraries sit "typically within 10-fold of the median", 83% of
constructs within 10-fold of it — against Twist's incoming claim of >90% within <2x of the mean.
The widening is the whole pipeline's, not the PCR's alone, and they do not separate the two. It
is evidence that a high-cycle subpool PCR is survivable, not evidence that it is free.

The literature supplies no 30, and no single number for a nested pair. What it supplies twice
over is a **stopping rule**.

## 5. What the evenness argument actually rests on

The ticket's premise is that evenness over the pool's 153 members is what the coverage is counted
on, and over-cycling is how it is lost. That is right about over-cycling and needs one correction
about where the evenness is held.

Molar representation across a part list is **re-established after assembly, not carried through
the PCRs**. `_pool_step` instructs: "Resuspend every block and measure each concentration. Pool
the members of each part list in equal molar amounts, one tube per position," and its
troubleshooting says to pool to the lowest member rather than the mean. Every block has its own
PCR2 well and its own assembly well, so a block that amplified unevenly is corrected by weighing.

What over-cycling at PCR1 costs that no later step repairs:

- **Dropout.** A pool member lost at PCR1 has no template at PCR2 and no block at assembly. Twist
  names dropout as a consequence of cycle number and of GC. An equimolar pool cannot restore a
  member that is not there.
- **Chimeras and heteroduplexes.** Twist's documented over-amplification products. These put
  wrong sequence into a block, which quantification does not see.
- **Polymerase error**, which accumulates with cycles.

PCR2 carries a risk of its own: its forward primer is the batch's, shared by every oligo in the
batch, and only the inner primer is block-specific. One shared primer means over-cycling can pull
a neighbour's product into the well — the step already watches for a band at PCR1's length.

So the case for few cycles survives, but the mechanism is dropout, chimera and error, not the
pool's molar evenness. The rendered protocol should give the reason that is true.

## 6. Recommendation

**PCR1 — a sourced number, and a sourced rule beside it.**

- Print the count from Twist's length band, chosen from the project's own `oligo_length`:
  **12–14 cycles** at 151–350 nt, 10–12 at 100–150 nt, 6–10 at 20–100 nt. Cited to Twist
  FRM-001034 REV 8 and DOC-4060 REV 1.0. The count then follows the design rather than sitting
  fixed in the source.
- Keep the stopping rule as the instruction beside it, in its sourced form rather than the present
  wording: run the reaction on a real-time instrument with an intercalating dye and stop before
  the curve plateaus (Lund 2024; Romanowicz 2026; and Qian, who puts EvaGreen in the same
  reaction). Where no real-time instrument is at hand, the table's count is the starting point and
  the symptom is the heteroduplex hump on capillary electrophoresis, not a gel smear.
- Record the **polymerase caveat as a hole** rather than transferring it silently. Twist states
  12–14 against KAPA HiFi HotStart or TrueAmp, and its own FAQ says another polymerase may give
  worse uniformity at the same count. The pipeline runs Q5, and nobody has published a count for
  Q5 on a Twist pool — Lund used Q5 Hot Start and declined to state one.

**PCR2 — a named hole, with the stopping rule as the instruction.** No source read here covers
it. Twist's table is for amplifying the pool as delivered; PCR2's template is PCR1's product.
Qian's 10-cycle Day 4.1 is a sequencing-amplicon PCR, not this. Freschlin's 35 is a first subpool
PCR. The honest instruction is the rule alone, with no printed count.

PCR2's hole is **H30**. H27 and H28 were taken before this note reached code, and H29 records the
polymerase caveat above. H31 to H34 have been taken since, in `igga/stages.py` and
`scripts/build_working_vector.py`, so **H35** is the first free id; `bench-numbers.md` collects
the holes and says which numbers are retired rather than free:

| Field | Value |
| --- | --- |
| `missing` | no source gives a cycle count for PCR2, which pulls one block out of an already-amplified batch |
| `kind` | `unpublished` |
| `where` | PCR2, cycle count |
| `filled_by` | a pilot that titrates PCR2 against the heteroduplex hump on capillary electrophoresis, or a real-time run stopped before plateau |

`docs/research/bench-numbers.md` H12 — "No cycling for PCR1 or PCR2 of the three-primer scheme" —
is **half closed** by this note. PCR1 now has a source; PCR2 still has none.

**What must stop either way.** The program must stop falling through to the polymerase default.
Thirty is not above every published count — Freschlin's 35 is higher — so the defect is not that
it is too many. The defect is that one printed number stands where every source consulted says
the count is chosen per pool: by length, by polymerase, and by watching the reaction.

## 7. What this implies in code

Not done here; #340 is a research ticket.

1. `_pcr1_step` passes an explicit `cycles`, chosen from `project.oligo_length` against Twist's
   three-band table, with a citation on the program.
2. `_pcr2_step` passes no count and carries a `Hole` instead, as `_assembly_step` already carries
   `stages.POOL_HOLES`.
3. PCR1's troubleshooting entry changes its observable from "a smear rather than a band" to the
   heteroduplex hump on capillary electrophoresis, and gains the real-time stopping rule.
4. A hole records that Twist's count is stated against KAPA HiFi HotStart or TrueAmp while the
   pipeline uses Q5.
5. `docs/research/synthesis-and-assembly.md` drops "a hard maximum of eight cycles" from
   "Twist — inherited", which is the Multiplexed Gene Fragments figure, and points at section 2
   here.

## 8. Open gaps

- **No count published for Q5 on a Twist oligo pool.** Lund used Q5 Hot Start and stated none.
- **Nothing measures what the nested pair costs in evenness.** Freschlin measured one selective
  subpool PCR; nobody measured a PCR1 followed by a PCR2.
- **Twist's table is stated without data.** Neither guide shows the experiment behind 12–14, and
  neither says what the band's endpoints mean — whether 14 is a cap or a typical.
- **IDT was not read as a protocol.** Its public position on oPools is that they ship at picomole
  yields, so no pre-amplification is needed; that comes from a product announcement, not a
  protocol document, and no number here rests on it. It is worth knowing as a vendor choice that
  removes the problem rather than bounding it.
- **No Agilent oligo-pool amplification protocol was fetched.**

## 9. Provenance

Downloaded for this note, into `reference_docs/synthesis_and_assembly/twist/protocol/`, with
`pdftotext -layout` dumps beside them in `protocol/dump/`:

| File | Source URL | Retrieved |
| --- | --- | --- |
| `FRM-001034_AmplifyingOligoPools_REV8_singles.pdf` | `twistbioscience.com/content/dam/twistbioscience/resources/2026-01/FRM-001034-AmplifyingOligoPools-REV8%20singles.pdf` | 2026-10-07 |
| `DOC-4060_Twist_Oligo_Pools_Amplification_Protocol_REV1_singles.pdf` | `twistbioscience.com/content/dam/DOC_4060_Twist_Oligo_Pools_Amplification_Protocol_REV1_singles.pdf` | 2026-10-07 |

Already held and re-read, with provenance in `docs/research/synthesis-and-assembly.md` section 7:
Twist DOC-001498 REV 1.0 and DOC-4057 REV 1.0; Lund et al. 2024; Freschlin et al. 2026;
Romanowicz et al. 2026. Qian et al. 2026's numbers are quoted from
`docs/research/bench-numbers.md`, which read the supplement directly.

## 10. Re-checked for #476: both holes stand

Searched on **2026-10-08**, the day after everything above. Nothing a vendor publishes could
have moved in a day, so this pass is not a re-run: it reads the two vendors section 8 left open
and writes a plain verdict on each hole, so that a pilot is priced against a search rather than
against a hunch.

### H29 — a cycle count for Q5 on a Twist oligo pool

**Searched.** Both Twist oligo-pool guides held here, for the string `Q5`: no hit in either.
FRM-001034 REV 8 names KAPA HiFi HotStart and DOC-4060 REV 1.0 names Twist TrueAmp, each as the
one polymerase it recommends, and neither names a second. Then a web search for any Twist or NEB
document stating a count for Q5 on an oligo pool, and the two other pool vendors' guides.

**Found.** No Twist document mentions Q5 at all. NEB publishes no oligo-pool amplification
protocol; its Q5 routine-PCR page returns HTTP 403 to `curl` and to WebFetch, as `neb.com` HTML
does throughout this repo's research, and a routine-PCR count would not be a pool count anyway.

The nearest thing found is **GenScript's Oligo Pool User Manual**, and it is worth recording
because it names Q5 and still does not close this:

> High-fidelity polymerases such as Phusion, Q5, KAPA, are recommended for PCR amplification.
> Taq polymerase is not recommended for use in amplification of oligo pools. Please follow the
> recommendations of the polymerase supplier for optimal PCR amplification conditions.

Its typical procedure prints **20-30 cycles** from 5-100 ng of template. Three reasons it does
not transfer. It is GenScript's pool, not Twist's, and the two arrive at different amounts.
20-30 is roughly twice Twist's 12-14 for the same job, so taking it would mean contradicting the
vendor whose pool the method buys, on that vendor's own product. And the same page defers twice
over — to the polymerase supplier in the quotation above, and then to the user: "It is
recommended to carry out pilot experiments to determine the optimal reaction conditions." It is
evidence that Q5 amplifies a pool, which nobody doubted, not a count for this one.

Agilent's SurePrint guideline surfaced in the same search and was not read as a document; what
is quoted of it is stated against Herculase II, so it could not close a Q5 hole either.

**Verdict: no source. H29 stands as `unpublished`**, filled by Twist stating a count against Q5,
or by a pilot on this pool. A document search is finished with this one; only the bench is left.

### H30 — a cycle count for PCR2

**Searched.** Both Oligo Pools guides for anything about a second round off a first PCR's
product: `subpool`, `sub-pool`, `re-amplif`, `reamplif`, `nested`, `second round`. In the Oligo
Pools documents there is nothing. The only hits anywhere in the Twist material are in the
**Multiplexed Gene Fragments** guide, which mentions "amplifying out subpools" as a reason you
might amplify at all and gives no count for it, and in the MGF design guidelines, where
"split into sub-pools" is an ordering remedy for the length-spread rule and not a PCR at all.
Section 3's warning applies to both: that is the other product.

**Found.** GenScript's manual is the one document read in this whole note that speaks to a
sub-pool PCR, and what it says is the opposite of a count: raise the template and the cycles to
pull more sub-pools out of one reaction, and pilot it. Nothing anywhere states a count for
pulling one block out of a batch that has already been amplified.

**The gap is structural, not accidental.** Every vendor protocol here is written for the pool as
delivered, because that is the reaction the vendor sells a product for. PCR2's template is
PCR1's product, which no vendor has in hand. That is why section 4's published counts are all
first amplifications, and why no amount of further searching is likely to produce one.

**Verdict: no source. H30 stands as `unpublished`**, filled by a pilot titrated against the
heteroduplex hump on capillary electrophoresis, or a real-time run stopped before the curve
plateaus. The printed stopping rule remains the honest instruction.

### One lead, not followed

A Research Square preprint, *Oligo replication advantage driven by GC content and Gibbs free
energy*, studies amplification bias on a commercial oligo pool and so touches both holes. Its
body did not come back readable through the fetch, and it is a study of how bias behaves rather
than a protocol that sets a count. Recorded here as a lead; nothing in this note rests on it.

### What this pass did not change

Section 2's 12-14 for PCR1 is still the sourced count, and still stated against a polymerase the
method does not run. Section 8's open gaps stand, with one narrowed: GenScript has now been read
as a protocol and found not to carry our number, and Agilent still has not.

## Sources

- Twist Bioscience, *Amplifying Twist Oligo Pools*, FRM-001034 REV 8.
- Twist Bioscience, *Twist Oligo Pools Amplification Protocol*, DOC-4060 REV 1.0.
- Twist Bioscience, *Amplifying Twist Multiplexed Gene Fragments*, DOC-001498 REV 1.0.
- Twist Bioscience, *Multiplexed Gene Fragments Design Guidelines*, DOC-4057 REV 1.0.
- S. Lund, V. Potapov, S. R. Johnson et al., "Highly Parallelized Construction of DNA from
  Low-Cost Oligonucleotide Mixtures Using Data-Optimized Assembly Design and Golden Gate",
  *ACS Synth. Biol.* **13**, 745–751 (2024), `doi:10.1021/acssynbio.3c00694`.
- E. Freschlin et al., "Scalable and cost-efficient custom gene library assembly from oligopools",
  *Sci. Adv.* **12**, eady2279 (2026), `doi:10.1126/sciadv.ady2279`.
- Romanowicz et al., "DropSynth-Gold: Golden Gate assembly in emulsions extends multiplexed gene
  libraries to greater lengths" (2026).
- Qian et al. 2026, supplementary protocol, as read in `docs/research/bench-numbers.md`.
- GenScript, *Oligo Pool User Manual*, read 2026-10-08 from
  `genscript.com/gsfiles/techfiles/oligo-pool-user-guide.pdf`, as a `pdftotext -layout` dump. It
  carries no revision number. Not kept under `reference_docs/`: it closes nothing, and section 10
  quotes everything of it that bears on either hole.

All four Twist documents are vendor material: cite and re-enter single facts by hand, never
mirror. The same line `docs/research/restriction-ligation.md` section 1 drew for NEB.

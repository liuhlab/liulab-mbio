---
search:
  exclude: true
---

# Protein library assembly: the method, its numbers, its licence

Research note for issue #56. The paper, its supplement and Table S1 were read on
**2026-09-14**, and every sequence measurement below was computed from Table S1 on that date.
It records the method behind the library-assembly pipeline of issue #55: how the two digests
work, the overhang standard and the frame arithmetic, the bench numbers with the place each
one is written, a discrepancy between the paper's two accounts of which enzyme cuts what, and
what the licence allows.

The source is Takacsi-Nagy, O. et al. (2026) Synthetic transcription factors designed by
domain recombination enhance CAR T cell antitumor function. *Cell* 189, 1-20,
[doi:10.1016/j.cell.2026.07.054](https://doi.org/10.1016/j.cell.2026.07.054), with its
supplemental figures (`mmc1.pdf`, Figure S1A and its legend) and `Table S1.xlsx`.

Table S1 is **not** in this repository and is not package data. The measurements here were
computed from a local copy of the workbook and are reproducible from it. Where a number is
measured rather than quoted, this note says so and says over what.

## 1. Licences, and what may ship

| What | Licence | Evidence | Can we ship it? |
| --- | --- | --- | --- |
| The article, its text and figures | CC BY 4.0 | licence line, p. 1 | **Yes** — quote and reuse with attribution |
| Figure S1A and its legend | CC BY 4.0, through the article | supplemental information to a CC BY article | **Yes** — with attribution |
| Table S1 sequences | CC BY 4.0, through the article | no notice of its own; see below | **Yes** — with attribution and a link to the licence |
| The plasmids themselves | not a copyright question | "Plasmid libraries are available on Addgene", p. 16 | n/a — a material transfer, not a licence |

The article's licence line, printed on p. 1 under the copyright:

> Cell 189, 1-20, September 17, 2026 (c) 2026 The Author(s). Published by Elsevier Inc.
> This is an open access article under the CC BY license
> (`http://creativecommons.org/licenses/by/4.0/`).

Every page of the article carries an `OPEN ACCESS` rule in its running head. Under
*Materials availability*, p. 16:

> All DNA sequences are available in supplemental tables. Plasmid libraries are available on
> Addgene. Please contact the lead contact for any other reagent requests.

and under *Supplemental information*, p. 17:

> Supplemental information can be found online at `https://doi.org/10.1016/j.cell.2026.07.054`.

**Verdict on the paper: it ships.** CC BY 4.0 permits redistribution and adaptation of the
whole article, including its figures, provided the source is credited and the licence named.
That is the same grant `docs/research/ligation-fidelity.md` relied on for Pryor et al. 2020,
and it is the strongest of the licences this repository has had to read.

**Verdict on Table S1: its sequences may ship as a worked example, with attribution.** The
workbook carries no licence notice of its own — checked, and not merely assumed: no cell in
any of its ten sheets contains the words licence, copyright, rights, reserved or permission,
and its only file metadata is an author name and two timestamps. It is supplemental
information to the article, the article's licence line is unqualified, and *Materials
availability* names the supplemental tables as the place the DNA sequences are published.
Nothing scopes the CC BY grant to the typeset article and away from the files distributed with
it.

Two further points both push the same way and neither is needed to reach the verdict. A short
DNA sequence is a measurement of a molecule rather than creative expression, so the copyright
question is weak in either direction. And the authors put the sequences in a supplemental
table precisely so that others could build on them.

**What the attribution has to carry**, wherever a sequence from Table S1 appears in this
repository: the citation above, the DOI, a link to
`http://creativecommons.org/licenses/by/4.0/`, and a note of anything changed. This clears the
later ticket that plans the paper's scheme as a worked example under `docs/`: it may ship, and
the scheme may ship with real sequence in it rather than a placeholder.

## 2. The part, and how the two digests work

Every part is one synthesised block. Lengths below are measured off Table S1; the three
external-stuffer and internal-stuffer sequences are constant across all 24 rows of each sheet.

```text
5' external stuffer │ domain CDS │ internal stuffer │ barcode │ 3' external stuffer
       29 bp           27-1,046      56 or 58 bp        11 bp          29 bp
   PmeI · BbsI                    BsaI · SrfI · BsaI               BbsI · PmeI
```

The constant regions, with every enzyme site in them located by search and the overhang each
cut yields computed from the enzyme's own offsets:

| Region | Length | Sites in it | Overhang it yields |
| --- | --- | --- | --- |
| 5' external stuffer, N | 29 | PmeI at the outer edge, BbsI | `CTCC` |
| 5' external stuffer, bZIP | 29 | PmeI at the outer edge, BbsI | `GGAG` |
| 5' external stuffer, C | 29 | PmeI at the outer edge, BbsI | `CCGA` |
| internal stuffer, N | 56 | BsaI at each end, SrfI in the core | `GGAG` and `AGCG` |
| internal stuffer, bZIP | 56 | BsaI at each end, SrfI in the core | `CCGA` and `AGCG` |
| internal stuffer, C (terminal) | 58 | BsaI at each end, SrfI in the core | `CTCC` and `AGCG` |
| internal stuffer, CAR vector | 58 | BsaI at each end, SrfI in the core | `CTCC` and `AGCG` |
| 3' external stuffer, all three | 29 | BbsI, PmeI at the outer edge | `AGCG` |

The C domain's internal stuffer and the CAR vector's internal stuffer are the same 58-mer,
`GGCTCCGGAGACCTATCAAGCCCGGGCAACAATGTGCGGACGGCGTTGGTCTCTAGCG`. A C-domain part also carries a
constant 54-bp T2A between its CDS and its internal stuffer.

**Two digests, two jobs.**

- The **internal digest, BsaI + SrfI, opens the destination.** BsaI's two sites flank the
  internal stuffer and excise it, leaving the destination open on the next position's entry
  overhang at one end and on `AGCG` at the other. SrfI's single site sits in the middle of the
  stuffer, so the excised piece is cut again and cannot go back in.
- The **external digest, BbsI + PmeI, releases the donor.** BbsI's two sites, one in each
  external stuffer, cut the part out. PmeI's two sites sit outside them, at the outer edge of
  each external stuffer, so each remnant left behind on the donor backbone carries one — the
  backbone is cut twice and cannot carry through.

Measured over a whole synthesised block: BsaI 2 sites, BbsI 2, SrfI 1, PmeI 2. So the blunt
enzymes destroy the piece nobody wants, but not symmetrically: **SrfI cuts the excised internal
stuffer once** and **PmeI cuts the released donor backbone twice**.

**The mechanism was simulated on the published sequences**, cutting with each enzyme's real
offsets and ligating only where both overhangs agree. All three rounds join, and the overhangs
match at every one:

| Round | Destination, cut internally | Donor, cut externally | Match |
| --- | --- | --- | --- |
| 1 | N part, BsaI: `GGAG` / `AGCG` | bZIP part, BbsI: `GGAG` / `AGCG` | yes |
| 2 | N+bZIP, BsaI: `CCGA` / `AGCG` | C part, BbsI: `CCGA` / `AGCG` | yes |
| 3 | CAR vector, BsaI: `CTCC` / `AGCG` | assembled TF, BbsI: `CTCC` / `AGCG` | yes |

For JUN-JUN-JUN the assembled insert is 1,204 bp.

**It is not one pot.** Each round is digest, SPRI, ligate, SPRI, electroporate, grow, prep
(METHOD DETAILS, pp. e4-e5). The product deliberately keeps both enzymes' sites, and that is
what lets the next round happen.

## 3. The overhang standard, and the frame arithmetic

Three positions plus the vector junction need four overhangs. Each is written twice — once as
an internal-stuffer prefix and once at the end of a 5' external stuffer — and that redundancy
is what encodes position:

| Junction | Overhang | Written as |
| --- | --- | --- |
| CAR vector to N | `CTCC` | vector internal stuffer prefix `GGCTCC`; N 5' external stuffer ends `CTCC` |
| N to bZIP | `GGAG` | N internal stuffer prefix `GGAG`; bZIP 5' external stuffer ends `GGAG` |
| bZIP to C | `CCGA` | bZIP internal stuffer prefix `CCGA`; C 5' external stuffer ends `CCGA` |
| every part's 3' end | `AGCG` | every internal stuffer ends `GGTCTCTAGCG`; every 3' external stuffer begins `AGCG` |

A part's internal stuffer begins with the **next** position's entry overhang, so the N stuffer
reads `GGAG` and only a bZIP part can follow it. Figure S1A prints the two internal junctions
as duplexes, `GGAG` over `CCTC` and `CCGA` over `GGCT`, which agrees with the overhangs
computed above.

### Barcodes accumulate backwards

Each round inserts its barcode upstream of those already there, joined by the `AGCG` scar, so
the finished block reads C, then bZIP, then N — the reverse of the assembly order. Measured on
the assembled product: 41 bp, `C │ AGCG │ bZIP │ AGCG │ N`. In the final knock-in construct a
third scar and the vector's own 11-bp barcode follow, making 56 bp in all: **four barcodes and
three scars**, not four of each.

### Why the barcode is 11 bp

The barcode-plus-scar unit is 11 + 4 = 15 bp, which is five codons, so appending one does not
shift the frame. That is the whole of the argument, and it picks 11 out of its neighbours:
10 + 4 = 14 and 12 + 4 = 16, neither a multiple of three. Eleven is the length that is both
frame-neutral against a 4-bp scar and long enough to hold a distance rule.

### Why the terminal stuffer prefix is 6 bp rather than 4

The terminal position's internal stuffer is never excised — assembly stops there — so all
58 bp of it stay in the product, ahead of the barcode block. Checked:

- with the 6-bp prefix `GGCTCC`: 58 + 56 = **114**, a multiple of three;
- with a 4-bp prefix the stuffer would be 56 bp: 56 + 56 = **112**, which is not.

The two extra bases buy the frame back. The non-terminal stuffers are 56 bp with 4-bp prefixes
(`GGAG`, `CCGA`) because they are excised and only their 4-bp overhang survives into the
product.

### The frame runs through four neutral units

The reading frame is set upstream in the vector, at the tNGFR start codon; it is **not** the
domain's own `ATG`, and assuming it is gives the wrong answer. Measured on the assembled
JUN-JUN-JUN construct:

| Unit | Measured | Multiple of three? |
| --- | --- | --- |
| N CDS + `GGAG` scar | 749 + 4 = 753 | yes |
| bZIP CDS + `CCGA` scar | 191 + 4 = 195 | yes |
| C insert + T2A | 45 + 54 = 99 | yes |
| terminal stuffer + barcode block | 58 + 56 = 114 | yes |
| total from the N domain's first base | 1,161 = 387 codons | yes |

This is not a property of JUN. Measured across all 72 parts: every N and bZIP CDS has a length
congruent to 2 modulo 3, so CDS plus its 4-bp scar is always a multiple of three, and every C
insert has a length divisible by three, so insert plus the constant 54-bp T2A always is too.
That is the bookkeeping the fixed overhang standard charges to the proteins, and Figure S1A's
legend is explicit about paying it:

> Domains were not separated by linker sequences so minor sequence modifications affecting one
> to two amino acids were made to standardize the four-bp overhangs (see Table S1).

Measured in Table S1: of 24 N parts, **15** have a changed final pair of residues (ATF2 ED to
EE, ATF4 EK to EE, ATF6 IA to ME, and so on), 9 are unchanged. Of 24 C parts, **19** have a
changed first residue and 5 do not — every synthesised C domain begins with R, from eight
different wild-type residues. The workbook records wild type beside synthesised in both
sheets, which is the bookkeeping a pipeline has to reproduce.

The barcode block begins **one base into a codon**, so barcodes are not codon-aligned; what is
frame-neutral is the 15-bp unit, not the barcode on its own. In that frame, none of the 13,824
published blocks carries an in-frame stop, and none of the 72 barcodes introduces one in its
own slot.

## 4. The enzyme discrepancy

The paper says which enzyme cuts which molecule twice, and the two statements disagree.

**METHOD DETAILS, p. e4:**

> 1 ug of plasmid pools containing TF domains were digested in two steps: first with 2.5 uL of
> BsaI restriction enzyme (New England Biolabs, NEB) in 50 uL total volume (with CutSmart
> Buffer, NEB) for 1 hour at 37 C and then 2.5uL of SrfI was added and incubation was continued
> for a second hour. The destination backbone plasmids were digested with BbsI and PmeI (NEB)
> using the same protocol.

**Figure S1A legend:**

> To append the "bZIP" domain to the "N" domain, the "bZIP" domain was cleaved externally
> (BbsI), and the "N" domain was cleaved internally (BsaI) and then ligated. The "C" domain was
> then appended to the "N+bZIP" domain-pair analogously.

Both agree on the **pairing** — BsaI with SrfI, BbsI with PmeI. They disagree on the
**assignment**. METHOD DETAILS gives the domain pools, which are the donor, BsaI + SrfI, and
gives the destination backbone BbsI + PmeI. The figure gives the donor the external enzyme and
the destination the internal one, which is the opposite way round.

Three things corroborate the figure, and none of them is a matter of taste.

1. **The figure's own panel labels.** `BsaI BsaI` is printed under the N domain plasmid
   library's internal stuffer, and `BbsI` ... `BbsI` under the bZIP domain plasmid library's
   two external stuffers. The panel and its legend agree with each other.
2. **Where the sites actually are** (section 2, measured). BsaI's two sites flank the internal
   stuffer; BbsI's two sit in the external stuffers. SrfI lies inside the internal stuffer's
   core, so it can only shred what BsaI excises. PmeI lies at the outer edge of each external
   stuffer, so it can only shred what BbsI leaves behind. Each blunt enzyme is co-located with
   the sticky enzyme it is paired with, and therefore acts on the same molecule.
3. **Only the figure's assignment assembles.** Simulated on the published sequences, internal
   BsaI on the destination and external BbsI on the donor produce matching overhangs at every
   round. Under the METHOD DETAILS assignment nothing ligates: cutting a donor part with BsaI
   excises that part's own internal stuffer and never releases the domain, and cutting the
   destination with BbsI cuts in its external stuffers, which is not where the destination has
   to open.

**Verdict: the package follows Figure S1A.** The internal digest, BsaI + SrfI, opens the
destination; the external digest, BbsI + PmeI, releases the donor. METHOD DETAILS has the donor
and the destination swapped. This is a verdict on the evidence above rather than a preference
between two sources: the figure's reading is the only one under which the paper's own published
sequences assemble, and that is a measurement.

What this does **not** change: the reagents, the volumes and the times in METHOD DETAILS are
unaffected, because both digests use the same protocol. Only the label on each tube moves.

## 5. The bench numbers, and where each is written

Every number is quoted from the paper. None of them is a default for this package: amounts are
computed from the user's own fragment lengths, and these are recorded so that a computed number
can be checked against a real experiment.

| Step | Number | Source |
| --- | --- | --- |
| Digest input | 1 µg of pooled plasmid | METHOD DETAILS, p. e4 |
| Digest volume and buffer | 50 µL total, CutSmart | METHOD DETAILS, p. e4 |
| First enzyme | 2.5 µL, 1 hour at 37 °C | METHOD DETAILS, p. e4 |
| Second enzyme | 2.5 µL added, incubation continued a second hour | METHOD DETAILS, p. e4 |
| The other digest | "using the same protocol" | METHOD DETAILS, p. e4 |
| Clean-up after digest | SPRI at **2×** volume ratio, eluted in H2O | METHOD DETAILS, p. e4 |
| Ligase and buffer | T7 ligase with StickTogether DNA Ligase Buffer | METHOD DETAILS, p. e4; NEB M0318L, Key Resources Table, p. e2 |
| Molar ratio | 1:1 | METHOD DETAILS, p. e4 |
| Backbone per reaction | 20 ng digested backbone per 200 µL ligation | METHOD DETAILS, p. e4 |
| Ligation | 30 minutes at room temperature | METHOD DETAILS, p. e4 |
| Clean-up after ligation | SPRI at **1×** volume ratio, eluted in H2O | METHOD DETAILS, pp. e4-e5 |
| Electroporation input | up to 100 ng of purified ligation product | METHOD DETAILS, p. e5 |
| Strain | Endura electrocompetent cells (Lucigen) | METHOD DETAILS, p. e5; Cat#60242-2, Key Resources Table, p. e1 |
| Instrument | Gene Pulser XCell (BioRad) | METHOD DETAILS, p. e5 |
| Recovery | 1 hour at **30 °C** with shaking | METHOD DETAILS, p. e5 |
| Outgrowth | 12-16 hours at **30 °C** | METHOD DETAILS, p. e5 |
| Plasmid prep | ZymoPURE II | METHOD DETAILS, p. e5 |
| Enzyme products | BsaI-HF v2 R3733L, SrfI R0629L, BbsI-HF R3539L, PmeI R0560L | Key Resources Table, p. e2 |

Two of these are worth flagging because a pipeline would get them wrong by default.

- **Both growth steps are at 30 °C**, not 37 °C. That is a library precaution, and a protocol
  that defaulted to 37 °C would not be following this method.
- **The two SPRI ratios differ**: 2× after the digest, 1× after the ligation. They are not one
  number used twice.

The method's own justification for building this way, p. 3:

> Compared with synthesizing individual genes, combinatorial construction reduced synthesis
> costs by ~100-fold.

## 6. The barcode set, measured

**The metric is Hamming distance on equal-length barcodes.** Every barcode in Table S1 is
11 bp, so no alignment is involved and Hamming is well defined. Each row below says how many
barcodes were compared and, where it is a minimum over pairs, how many pairs.

| Property | Measured | Over |
| --- | --- | --- |
| Barcode length | 11 bp, every one | every barcode in the workbook |
| The three AP-1 domain pools | 24 + 24 + 24 = **72** distinct | the N, bZIP and C sheets |
| Unique 11-mers anywhere in Table S1 | **164** | all ten sheets |
| Assembled barcode block | **41 bp** | 13,824 published blocks |
| Scar at the first scar position | `AGCG`, 13,824 of 13,824 | positions 11-15 |
| Scar at the second scar position | `AGCG`, 13,824 of 13,824 | positions 26-30 |
| Distinct barcodes per slot | **24** in each of the three slots | 13,824 blocks |
| Slot order | slot 1 = C pool, slot 2 = bZIP, slot 3 = N, each matching its sheet exactly | 13,824 blocks |
| Set completeness | all 24³ = **13,824** combinations, each exactly once, no duplicates | 13,824 blocks |
| Minimum Hamming, N pool | **3** | 24 barcodes, 276 pairs |
| Minimum Hamming, bZIP pool | **4** | 24 barcodes, 276 pairs |
| Minimum Hamming, C pool | **4** | 24 barcodes, 276 pairs |
| Minimum Hamming, the three pools together | **2** | 72 barcodes, 2,556 pairs |
| Minimum Hamming, every 11-mer in Table S1 | **2** | 164 barcodes, 13,366 pairs |
| GC content | **18.2% to 81.8%** (2 of 11 to 9 of 11) | the 72, and the 164, alike |
| Longest homopolymer run | **6** | the 72, and the 164, alike |
| BsaI, BbsI, SrfI or PmeI site in a bare barcode | **none** | the 72, and the 164 |
| A new enzyme site created in context | **none** | 13,824 blocks with their real flanks |
| In-frame stop codon | **none** | 13,824 blocks, and the 72 in their own slots |

The 164 unique 11-mers break down as: 24 N, 24 bZIP, 24 C, 30 in the natural full-length AP-1
sheet, 39 in the individually synthesised DESynR sheet, 37 control barcodes in the barcode
reference sheet, 16 in the HA GD2 CAR sheet and 1 in the vector — a union of 164 after
duplicates across sheets are removed.

### The distance rule is per pool, and the paper overstates it

METHOD DETAILS, p. e7:

> Each synthesized TF domain contains a unique domain barcode (BC), with a minimum Hamming
> distance of three between any two domain BCs.

Measured, that holds **within** each pool — 3, 4 and 4 — but not **across** them. Two barcodes
in different pools are 2 apart: `AACAAATTAGT` and `AACAGATTGGT`. So the paper's "any two domain
BCs" is not met by the paper's own set, and the rule that the set actually satisfies is a
per-pool one. This matters for the package: distance has to be checked inside a position's
pool, which is also the only place it is needed, since two barcodes in different slots are
never confusable.

### No GC rule and no homopolymer rule was imposed

GC spans 18% to 82% and a homopolymer run reaches 6, so no band and no cap was applied. The
set still worked: the printed data labels in Figure S1B are 84.1% correct, 8.3% barcode
mismatches, 1.0% missing domains, 4.5% missing barcode and 2.0% other, and Figure S1C is
annotated "93% of library is >80% correctly linked". The main text, p. 3, rounds these to:

> Of reads containing three valid BCs (~95%), nearly 90% of the library was correctly linked

Note that the figure's own categories give 93.4% of reads with three barcodes, against the
text's ~95%, and 91.0% correctly linked among them, against "nearly 90%". The two accounts are
consistent to a point or two, and this note records both rather than picking one.

**That loose worked is not evidence that loose is optimal**, nor that the question was
considered. Whether a GC band or a homopolymer cap is worth imposing is for issue #58 to settle
from primary sources; nothing here should be read as a default.

## 7. Where this note contradicts issue #55

Issue #55 was written before these measurements existed. Most of it holds exactly; four things
do not.

| Issue #55 says | Measured here |
| --- | --- |
| "the 148 unique 11-mer barcodes in Table S1" | **164** unique 11-mers across the workbook, of which **72** are the three AP-1 domain pools. 148 is neither |
| "barcodes are at least 3 mismatches apart within a pool ... within a pool it is 3" | true as a floor, but tight only in the N pool. The bZIP and C pools are both **4** |
| "Each digest cuts its unwanted piece twice" | PmeI cuts the donor backbone **twice**; SrfI cuts the excised internal stuffer **once** |
| "the terminal part's stuffer prefix is six bases ... to correct the frame: 54 + 15 × 4 = 114" | the total 114 is right, the decomposition is wrong. The unit is the 58-bp terminal stuffer plus the 56-bp barcode block. 54 is the T2A length, and the block is 56 bp — **three** scars, not four |
| "Confirmed on all 13,824 rows: 41 bp, AGCG at both scar positions, 24 distinct barcodes per slot, slot order matching the C, bZIP and N sheets" | confirmed exactly, and the set is the complete 24³ product with no duplicates |
| "GC runs 18% to 82% and homopolymer runs reach 6" | confirmed |
| "free of stop codons in the construct's real frame" | confirmed, in the frame measured in section 3 |
| "the measured minimum across all 148 barcodes is 2" | the minimum across the three pools is **2**, and across all 164 it is also 2 |
| "23 N-, 27 ETS-" domains | that is the ETS and FOX library, a different one. The AP-1 library is 24 / 24 / 24 |

One number in the paper does not survive either: the main text, p. 3, gives the 72 domains as
"50-1,000 bp", while the measured CDS lengths run **44 to 1,046 bp** (N 44-1,046, bZIP 191-263,
C insert 27-942 before its constant T2A). The paper's range is a round summary, not a bound a
designer can rely on.

## Open gaps

Nothing in this section should be turned into a package default. Each is a value the paper does
not state, or a claim it makes that its own data does not support.

- **Enzyme units are not given**, only volumes: 2.5 µL of each. No stock concentration is
  stated, so units per µg cannot be derived.
- **T7 ligase amount is not given** — neither units nor volume, only the buffer and the 1:1
  molar ratio.
- **What the 1:1 molar ratio is between** is not said. Insert to backbone is the obvious
  reading; the paper does not name the two species.
- **The destination digest's temperature is not restated.** "Using the same protocol" implies
  two hours at 37 °C for BbsI and PmeI as well, but it is an inference.
- **No heat inactivation is described** at any step; the protocol goes digest straight to SPRI.
- **The SPRI bead product is not named** in the cloning paragraph, only the volume ratios.
- **No electroporation settings**: the instrument is named, the voltage, capacitance,
  resistance and cuvette gap are not.
- **No colony count or coverage target per round.** Figure S1G annotates coverage for the pilot
  library as a schematic, and reading numbers off it is not a citation.
- **Why the barcode is 11 bp is never stated.** The frame argument in section 3 is derived from
  the sequences here, not quoted from the paper.
- **Whether a GC band or homopolymer cap was considered is not stated.** Their absence from the
  set is measured; the authors' intent is not recoverable.
- **The paper's minimum Hamming distance of three is not met across pools** (measured 2), and
  the paper does not say the rule is per pool.
- **The stated domain size range 50-1,000 bp does not match** the measured 44-1,046 bp.
- **The donor backbone was not measured whole.** Table S1 does not give the pTwist Kan vector
  sequence, so "PmeI cuts the donor backbone twice" is counted on the two external-stuffer
  remnants, which is where both sites lie, and not over the entire plasmid.
- **The pilot library cannot be separated from the main one.** The barcode reference sheet is
  labelled "Pilot + Big", and the pilot's 9 × 3 × 8 design is given only in Figure S1H.
- **Table S1 carries no licence notice of its own.** The verdict in section 1 rests on the
  article's CC BY grant covering its supplemental files, which is the ordinary reading but is
  not stated inside the workbook.

## Sources

The paper, its supplement and Table S1 were read on 2026-09-14.

- Takacsi-Nagy, O., Kasinathan, S., Hartman, A., Yin, Y., Wu, L., Chen, A.Y., Moser, L.M.,
  May, A.P., Reeder, G.C., Celallos Fuentes, E., Kernick, C., Lu, J., McClellan, A.K.,
  Raposo, C.J., Terrall, B., Theberath, N.E., Yan, P.K., Xu, P., Sotillo, E., Eyquem, J.,
  Mackall, C.L., Roth, T.L. and Satpathy, A.T. (2026) Synthetic transcription factors designed
  by domain recombination enhance CAR T cell antitumor function. *Cell* 189, 1-20.
  [doi:10.1016/j.cell.2026.07.054](https://doi.org/10.1016/j.cell.2026.07.054). CC BY 4.0.
- The same paper's STAR Methods: Key Resources Table (pp. e1-e2), METHOD DETAILS, *DNA
  Constructs and HDR Templates* (pp. e4-e5), *Barcode Abundance Analysis* (p. e7) and
  *Long-Read Sequencing* (p. e9).
- The same paper's *Resource availability* and *Supplemental information* (pp. 16-17).
- Supplemental information `mmc1.pdf`: Figure S1, panel A and its legend.
- `Table S1.xlsx`, ten sheets. The sequence and barcode measurements in sections 2, 3, 6 and 7
  were computed from it; the workbook is not redistributed here.
- Creative Commons Attribution 4.0 International:
  [creativecommons.org](http://creativecommons.org/licenses/by/4.0/)
- `docs/research/ligation-fidelity.md` for the CC BY precedent this note follows, and
  `docs/research/restriction-enzyme-data.md` (issue #7) for the cut-offset convention used to
  compute every overhang above.

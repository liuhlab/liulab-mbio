---
search:
  exclude: true
---

# Materials, by the process that uses them

Answers "what do I need for this step". The method page, `docs/synthesis-and-assembly.md`,
answers the other question, "what do I buy", under `## Reagents and equipment`, and owns the
catalogue: vendors, accessions and part numbers live there and are not repeated here.

What this file adds is the per-process grouping and the sequences. Each is quoted from a source
and marked with it; anything still unknown is marked **pending** rather than guessed.

## Lab resources (build once)

| Item | What it is here |
| --- | --- |
| pCR-Blunt II-TOPO | Part carrier. KanR + Zeocin, ccdB, no BsmBI site |
| DMX0001 / DMX0002 | Parent of the DMX vector, modified before first use |
| Working vector backbone | Chosen per application; must have no BsaI or BsmBI and use AmpR/CarbR |
| RFP stuffer cassette | Fills the cassette site of the input working vector; long enough that the cut backbone separates on a gel, and fluorescent, so a colony that lost it is the one to pick |
| ccdB cassette, four versions | N and C, no N, no C, neither |
| Orthogonal primer set | 185 primers of 20 nt, 165 of them in the orthogonal set; Subramanian et al. 2018, Supplementary Table 1 |
| DMX barcode kit | 96 plasmids, four groups of 24; `docs/research/synthesis-and-assembly-barcode-kit.md` |
| DB3.1 or ccdB Survival 2 T1R | Propagating a stock whose ccdB is expressed. NEB Stable is ccdB-sensitive and will not grow one — #300 |

The DMX vector is rebuilt from its parent before first use: AmpR replaced with KanR, the four
BbsI sites removed (two in `lacI`, two in the backbone), one PmeI site added outboard of each
BsaI site. One of DMX0001's three BsaI sites sits inside AmpR, so the marker swap removes it.

Primer plates are laid out once and kept: source plates by role, plus pre-made combination
plates holding each slot's pair. A slot always means the same pair.

## Working cassette

| Item | Note |
| --- | --- |
| Input working vector | from lab resources |
| Parts | from the part collection |
| BsmBI, ligase | one reaction cuts every part |
| Retailoring primers | only for a part whose default pair is wrong: `5'-[spacer]-CGTCTC-N-[new overhang]-[annealing]-3'` |
| DB3.1 or ccdB Survival 2 T1R | while an expressed ccdB is present |

## Cargo synthesis

| Item | Note |
| --- | --- |
| Twist oligo pool | product and length chosen at design time, from Twist's own product and pricing material |
| Primer source and combination plates | the P1 / P2 / P3 roles |
| BsmBI, ligase | assembly into the DMX vector |
| ccdB-sensitive strain carrying T7 polymerase | receiving cargo; without the polymerase the cassette is never transcribed and nothing counter-selects |
| Electrocompetent BL21(DE3) | for a library |
| Glycerol | archive plates |

The DMX paper's own library primers, which carry the BsmBI sites that give `AGGA`/`TTCC`
(Supplementary Table 4):

| Primer | Sequence |
| --- | --- |
| Forward library primer | `atactacCGTCTCgaggaGGTGGATCAGGAGGTTCG` |
| Reverse library primer | `gcattacCGTCTCcggaaCCACTTCCACCGCTTCC` |

These amplify a whole library in one reaction. The three-primer scheme replaces them with the
per-batch and per-gene pairs drawn from the orthogonal set.

## Cargo validation

Two read-outs, one of which is chosen.

### DMX barcoding

| Item | Note |
| --- | --- |
| DMX barcode kit | four UMIs per well; `docs/research/synthesis-and-assembly-barcode-kit.md` |
| BsaI, ligase | barcoding runs in cell lysate |
| Acoustic dispenser | required for the 1536-well barcoding reaction |
| DMX1–DMX6 primers | three pairs, run separately and pooled |
| Nanopore sequencing | long read across design and barcodes |

| Primer | Sequence |
| --- | --- |
| DMX 1_rv | `TCATGCCGCGACGCATACTTTAGAC` |
| DMX 2_fw | `CTGGTTGGGATCAGTCGCTTAGTGC` |
| DMX 3_rv | `CTGGGTACTGTCCCAAGGACCACTC` |
| DMX 4_fw | `CGACACCGAACGTGCGACAAAACTA` |
| DMX 5_rv | `CTCTGCGACCTAGCCGGTTTATCCA` |
| DMX 6_fw | `TTTATTCAGGGCACTACCCGGAGCT` |

### Index PCR, LevSeq-style

| Item | Note |
| --- | --- |
| Barcoded primer plate | Long et al., Supporting Information Table S3 (forward) and Table S4 (reverse), the full-length barcode-linked primers; Tables S5–S12 are the plate maps and S14 the plate-to-library pairing |
| Nanopore sequencing | as above |

Their primers are backbone-specific, annealing outside the insert, so they were designed against
the paper's plasmid and not ours. Measured on 2026-10-06, both constant backbone tails match
both DMX parents at full length and each pair brackets the ccdB cassette, so the plates are
ordered as published and no redesign is needed. The DMX route does not have this problem,
because its barcodes attach by Type IIS and ignore the backbone.

## iGGA

| Item | Note |
| --- | --- |
| BbsI + SrfI | opens the destination |
| BsaI + PmeI | releases the donor and blunts its backbone |
| T7 ligase, StickTogether buffer | will not join blunt ends, which is what removes escapees |
| SPRI beads | 2X after each digest, 1X after the ligation; both eluted in water |
| Endura electrocompetent cells | recombination-deficient |
| LB agar plates, selection for the destination | one dilution plate and one no-donor control plate per round |
| Capping block | optional last round; any `AGGA`/`TTCC` block that is 1 mod 3, such as T2A |

The two digests never share a tube. The two bead ratios are not one ratio used twice: the 2X
step is what removes the 17 bp and 13 bp cut stubs, and the PmeI argument in the departures note
rests on it.

## Final assembly

| Item | Note |
| --- | --- |
| Working vector | with its cassette |
| Cargo in the DMX vector | one well, or a pool |
| BsaI, ligase | one pot |
| SPRI beads, electrocompetent cells | — |

No counter-selection reagent is needed beyond what the vectors carry: the released DMX backbone
is KanR on an Amp plate, and released ccdB kills anything that takes it.

## Library reads

| Read | Platform | Where it runs |
| --- | --- | --- |
| Linkage, design to barcodes | long read | once, on the finished iGGA library |
| Representation | short amplicon across the barcodes | after every bottleneck |

## Looking something up by kind

The method page's `## Reagents and equipment` lists everything by kind — plasmids, primers,
reagents, equipment — with its source. This file does not repeat it.

---
search:
  exclude: true
---

# Golden Gate assembly with NEBridge, overhang fidelity and other assembly methods

Research note for issue #4. Everything below was retrieved on **2026-09-12** unless a
different date is given next to the source.

## How to read this note

`neb.com` returns HTTP 403 to both `curl` and WebFetch for its HTML pages. Three routes work
and every fact below came from one of them:

| Route | Works | Used for |
| --- | --- | --- |
| `neb.com/-/media/nebus/files/manuals/*.pdf` | yes (WebFetch) | kit manuals E1601, E1602 |
| `neb.com/-/media/catalog/specifications/*.pdf` | yes (`curl`) | per-enzyme product specifications |
| Wayback snapshots of `neb.com` HTML | yes | product pages, protocol and chart pages |

Values taken from an archived snapshot carry the word **archived** and the snapshot date. That
date matters: NEB revises these tables, and two of them have already changed (see
[Two NEB revisions to watch](#two-neb-revisions-to-watch)). Anything still unverified is in
[Open gaps](#10-open-gaps) and is never guessed.

## 1. Sources and licences

| Dataset or document | Licence | Can we ship it? |
| --- | --- | --- |
| Pryor et al. 2020, *PLoS One*, S1–S5 Tables (xlsx) | CC BY 4.0 | **Yes**, with attribution |
| Potapov et al. 2018, *ACS Synth. Biol.*, Supporting Data (figshare) | CC BY-NC 4.0 | **No** — non-commercial clashes with our MIT licence |
| Potapov et al. 2018, *Nucleic Acids Research* (gky303) | CC BY-NC 4.0 | **No**, same reason |
| `potapovneb/ligase-fidelity` analysis scripts | AGPL-3.0 | **No** — do not vendor into an MIT package |
| NEB manuals, protocols, charts, specification sheets | © NEB, all rights reserved | **No** — cite and paraphrase numbers, never redistribute the files |
| REBASE enzyme records | See [Open gaps](#10-open-gaps) | Licence statement not found on the pages read |

Licence text, quoted:

- Pryor 2020: "This is an open access article distributed under the terms of the Creative
  Commons Attribution License, which permits unrestricted use, distribution, and reproduction
  in any medium, provided the original author and source are credited."
  (`journals.plos.org` article XML, `<permissions>`).
- Potapov 2018 NAR: "This is an Open Access article distributed under the terms of the
  Creative Commons Attribution Non-Commercial License
  (`http://creativecommons.org/licenses/by-nc/4.0/`) …".
- Potapov 2018 ACS Supporting Data: figshare API reports
  `{"name": "CC BY-NC 4.0"}` for both items of collection 4283045.

## 2. Reaction setup per enzyme

NEB ships two different one-pot systems, and they have different tables. Do not mix them.

- **Kits** E1601 (BsaI-HFv2) and E1602 (BsmBI-v2): one enzyme mix that already contains
  T4 DNA Ligase, plus the pGGAselect destination plasmid.
- **NEBridge Ligase Master Mix** (M1100), a 3X ligase master mix you combine with **your
  choice** of NEB Type IIS enzyme. This is the only NEB table that covers BbsI-HF, Esp3I,
  PaqCI, SapI and BspQI.

### 2.1 Kit reaction, E1601 and E1602

Both manuals give the same table (E1601 manual v5.0_6/26; E1602 manual):

| Reagent | Assembly reaction |
| --- | --- |
| Destination plasmid | 0.05 pmol (pGGAselect = 75 ng = 0.05 pmol) |
| Inserts, precloned | 0.05 pmol each |
| Inserts, amplicon | 1:1 molar ratio (insert:vector) |
| T4 DNA Ligase Buffer (10X) | 2 µl |
| NEBridge Golden Gate Enzyme Mix | 1–2 µl |
| Nuclease-free water | to 20 µl |

Footnotes that carry real constraints:

- Enzyme mix volume: "For assemblies ≤ 10 inserts, use 1 µl; for assemblies > 10 inserts, use
  2 µl."
- Volume: "Can be increased to 25 µl volume if required due to DNA component volumes; add
  additional 0.5 µl T4 DNA Ligase Buffer (10X)."
- Amplicon inserts "must possess 5´ flanking bases (6 recommended) and BsaI restriction sites
  at both ends of the amplicon and in the proper orientation."
- Molar conversions: NEB points at NEBioCalculator.

### 2.2 Kit cycling programs

| Insert number | E1601 (BsaI-HFv2) | E1602 (BsmBI-v2) |
| --- | --- | --- |
| 1 insert | 37 °C 5 min (cloning) or 37 °C 1 h (library) → 60 °C 5 min | 42 °C 5 min (cloning) or 42 °C 1 h (library) → 60 °C 5 min |
| 2–10 inserts | (37 °C 1 min → 16 °C 1 min) × 30 → 60 °C 5 min | (42 °C 1 min → 16 °C 1 min) × 30–60 → 60 °C 5 min |
| 11–20+ inserts | (37 °C 5 min → 16 °C 5 min) × 30 → 60 °C 5 min | (42 °C 5 min → 16 °C 5 min) × 30–60 → 60 °C 5 min |

The 60 °C step is not heat inactivation of the ligase; it is a background-reducing digest.
NEB, FAQ 10/11 in the two manuals: "The final incubation step at 60°C favors Type IIS
restriction enzyme cutting, in the absence of DNA ligation. Digesting any uncut or
cut/religated destination plasmid still present in the assembly reactions reduces
background."

More cycles are allowed and help: "efficiency increases dramatically from 30 cycles to 60-65
cycles, with no loss of fidelity" (E1601 manual, FAQ 14).

### 2.3 NEBridge Ligase Master Mix (M1100) — the per-enzyme table

Reaction, from the NEB M1100 protocol page (**archived, 2023-03-31**):

| Component | 2-fragment or 3–6 fragment | 7+ fragment |
| --- | --- | --- |
| NEBridge Ligase Master Mix | 5 µl | 10 µl |
| DNA fragments | 0.05 pmol each | 0.05 pmol each |
| Type IIS restriction enzyme | x µl (below) | x µl (below) |
| Nuclease-free water | to volume | to volume |
| Total reaction volume | 15 µl | 30 µl |

Enzyme amounts, from the same protocol page and from the NEBridge Ligase Master Mix Protocol
Guidelines (**archived, 2025-07-13**); the two agree except where noted:

| Enzyme | 2 fragment | 3–6 fragment | 7+ fragment |
| --- | --- | --- | --- |
| BbsI-HF (R3539) | 1 µl (20 U) | 1 µl (20 U) | 1 µl (50 U) — needs R3539M at 50 U/µl |
| BsaI-HFv2 (R3733) | 1 µl (20 U) | 1 µl (20 U) | 1 µl (20 U) |
| BsmBI-v2 (R0739) | 3 µl (30 U) | 3 µl (30 U) | 6 µl (60 U) — "Use of less enzyme will reduce performance" |
| BspQI / BspQI-HF (R0712 / R3712) | 1 µl (10 U) | 1 µl (10 U) | 2 µl (20 U) |
| Esp3I (R0734) | 2 µl (20 U) | 3 µl (30 U) | 4 µl (40 U) |
| PaqCI (R0745) | 1 µl (10 U) | 1 µl (10 U) | 2.5 µl (25 U) |
| SapI (R0569) | 1 µl (10 U) | 1 µl (10 U) | 2 µl (20 U) |

Cycling, from the M1100 protocol page (**archived, 2023-03-31**). Enzymes split into a 37 °C
group and a 42 °C group:

| Complexity | BbsI-HF, BsaI-HFv2, Esp3I, PaqCI, SapI | BsmBI-v2, BspQI |
| --- | --- | --- |
| 2 fragment, single gene cloning | 37 °C for 15 min | 15 cycles of 42 °C 1 min / 16 °C 1 min |
| 2 fragment, library construction | 37 °C for 60 min | 30 cycles of 42 °C 1 min / 16 °C 1 min |
| 3–6 fragment | 30 cycles of 37 °C 1 min / 16 °C 1 min | 30 cycles of 42 °C 1 min / 16 °C 1 min |
| 7–13 fragment | 30 cycles of 37 °C 5 min / 16 °C 5 min | 30 cycles of 42 °C 5 min / 16 °C 5 min |
| 14+ fragment | 60 cycles of 37 °C 5 min / 16 °C 5 min | 60 cycles of 42 °C 5 min / 16 °C 5 min |

Then, for every enzyme and every tier: "End Soak: Incubate at 60°C for 5 minutes, before
transformation. Chill on ice. Use 2 μl of the reaction to transform 50 μl of competent cells.
If reaction will not be used immediately for transformation, store at -20°C."

Order of addition matters enough that NEB spells it out: DNA + water first, then master mix
(mix by pipetting 3 times), then the Type IIS enzyme (mix by pipetting 5 times), on ice.

### 2.4 PaqCI needs an activator

PaqCI is a tetramer that must engage several copies of its recognition site before it cuts, so
NEB supplies a short activator duplex. The mechanism is in the open-access structure paper
(Kennedy et al. 2023, *Nucleic Acids Res.* 51, 4467–4487, CC BY): PaqCI scans for identical
7 bp targets and synapses them through its target-recognition domains, and the activator is "a
short double-stranded hairpin DNA construct spanning the PaqCI recognition site" that "does not
extend past the point of cutting", occupying spare target-recognition domains so the
endonuclease domains can re-dimerise. That paper used a 5:1 activator:enzyme molar ratio
in vitro.

For bench use, two NEB tables exist and they differ:

- Ligase Master Mix guidelines (**archived, 2025-07-13**): "Recommended PaqCI Activator :
  PaqCI ratio is 1:1 (pmol:U). Use 0.5 µl of PaqCI Activator (20 µM) for 2 and 3-6 fragments;
  1.25 µl of PaqCI Activator (20 µM) for 7+ fragments."
- Usage Guidelines for Golden Gate Assembly with PaqCI, the T4-DNA-Ligase version
  (**archived, 2021-06-15**):

  | Assembly complexity | PaqCI | T4 DNA Ligase | PaqCI Activator |
  | --- | --- | --- | --- |
  | Single insert cloning (10 min 37 °C) or library prep (60 min 37 °C) | 5 U | 200 U | 5 pmol (¼ µl of 20 µM) |
  | Simple to moderate, 2–10 fragments | 5–10 U | 200–400 U | 5 pmol (¼ µl of 20 µM) |
  | Complex, 11–20+ fragments | 10–20 U | 400–800 U | 5–10 pmol (¼–½ µl of 20 µM) |

  Same page: "The activator solution is in a Mg-free buffer for optimal long-term storage. For
  short-term working stocks, if desired, dilute an appropriate amount in 1X T4 DNA Ligase
  Buffer to achieve more easily pipettable volumes (e.g., a four-fold dilution = 5 µM, 5
  pmoles/µl activator.)"

### 2.5 Home-brew reactions (individual enzyme + T4 DNA Ligase)

NEB's own high-complexity protocol, from the technical note *Breaking through the Limitations
of Golden Gate Assembly* (NEB, 01/19, retrieved from the `neb-online.de` mirror):

| Reagent | 24-fragment assembly |
| --- | --- |
| pGGA destination plasmid, 75 ng/µl | 1 µl (75 ng) |
| 24 precloned inserts, 100 ng/µl each | 0.75 µl (75 ng) each, 18 µl total |
| T4 DNA Ligase Buffer (10X) | 2.5 µl |
| T4 DNA Ligase (M0202), 2000 U/µl | 0.5 µl (1000 units) |
| BsaI-HFv2 (R3733), 20 U/µl | 1.5 µl (30 units) |
| Water | 1.5 µl |

Program: "(5 min 37°C → 5 min 16°C) x 30 cycles followed by 5 min 60°C". For ≤10 fragments the
same note says 500 U T4 DNA Ligase and 15 U BsaI-HFv2 suffice. Caveat from the E1602 manual
(FAQ 5): home brews work but "will not match the efficiency and fidelity achieved with this
kit", because glycerol from stock enzymes has to stay at or below about 10%.

Buffer choices, from NEB's *Technical Tips for Optimizing Golden Gate Assembly Reactions*
(**archived, 2020-10-26**) and the 2024 Golden Gate brochure: T4 DNA Ligase Buffer is the
default for BsaI-HFv2, BsmBI-v2 and PaqCI; alternatives are NEBuffer r1.1 (BsaI-HFv2),
NEBuffer r2.1 (BsmBI-v2) or rCutSmart (PaqCI), each "supplemented with 1 mM ATP and 5–10 mM
DTT".

### 2.6 Per-enzyme properties

Sites and cut offsets are REBASE records cross-checked against NEB's own table in the 2024
Golden Gate brochure; they agree. Incubation, heat inactivation and methylation come from NEB
product pages (**archived, 2021-04**) and from NEB's Heat Inactivation and
Dam-Dcm-and-CpG-Methylation charts (**archived, 2021-04**). Unit definitions come from the
specification PDFs fetched directly today.

| Enzyme | Site (cut offsets) | Overhang | Golden Gate temp | Heat inactivation | Dam / Dcm / CpG |
| --- | --- | --- | --- | --- | --- |
| BsaI-HFv2 (R3733) | `GGTCTC(1/5)` | 4 nt, 5′ | 37 °C | 80 °C, 20 min | Not sensitive / impaired by some overlapping / blocked by some overlapping |
| BsmBI-v2 (R0739) | `CGTCTC(1/5)` | 4 nt, 5′ | 42 °C | 80 °C, 20 min | CpG blocked; Dam/Dcm cells unread in the archived chart |
| Esp3I (R0734) | `CGTCTC(1/5)` | 4 nt, 5′ | 37 °C | 65 °C, 20 min | Not sensitive / not sensitive / blocked |
| BbsI-HF (R3539) | `GAAGAC(2/6)` | 4 nt, 5′ | 37 °C | 65 °C, 20 min | Not sensitive / not sensitive / not sensitive |
| PaqCI (R0745) | `CACCTGC(4/8)` | 4 nt, 5′ | 37 °C | 65 °C, 20 min | Not sensitive / not sensitive / impaired by overlapping |
| SapI (R0569) | `GCTCTTC(1/4)` | **3 nt**, 5′ | 37 °C | 65 °C, 20 min | Not sensitive / not sensitive / not sensitive |
| BspQI (R0712) | `GCTCTTC(1/4)` | **3 nt**, 5′ | 42 °C per M1100 table | 80 °C, 20 min | Not sensitive / not sensitive / not sensitive |
| BtgZI (R0703) | `GCGATG(10/14)` | 4 nt, 5′ | **no NEB Golden Gate protocol** | 80 °C, 20 min | Not sensitive / not sensitive / impaired |

Notes worth carrying into code:

- REBASE prototypes and isoschizomers: BsaI's prototype is Eco31I; BsmBI's prototype is
  Esp3I (so BsmBI-v2 and Esp3I are isoschizomers with **different** optimal temperatures, 42 °C
  and 37 °C — E1602 manual FAQ 4); BbsI's prototype is BbvII; PaqCI's prototype is AarI;
  SapI is its own prototype and BspQI is the isoschizomer.
- BsmBI-v2 is only 25% active in rCutSmart (product page, archived); BsaI-HFv2, BbsI-HF, SapI
  and PaqCI are 100%.
- BtgZI's unit definition is at **60 °C** (specification PS-R0703S/L v1.0), and its terminal
  integrity is much weaker than the others: "After a 5-fold over-digestion of Lambda DNA with
  BtgZI, ~75% of the DNA fragments can be ligated … Of these ligated fragments, ~75% can be
  recut". NEB names BtgZI only as a fallback specificity (E1601 manual FAQ 4), never with a
  reaction table. Treat BtgZI as last-resort in our ranking.
- SapI, per the specification sheet, is qualified for 4-hour incubations; NEB's 3-nt overhang
  means overhang design has 64 sequences to work with, not 256.
- BsaI and BsaI-HFv2 are "blocked by overlapping dcm methylation (methylation at the C5
  position of cytosine in the sequences CCAGG or CCTGG)", which matters only when a
  `CC(A/T)` sits immediately 5′ of `GGTCTC` (E1601 manual, tip 6). This is a real constraint
  on spacer-base choice in primer tails.

### Two NEB revisions to watch

1. The E1601 manual changed its reaction amounts in revision 5.0 (6/26): earlier NEB copy
   (protocol page text, mirrored by an iGEM team, undated copy) used a **2:1** insert:vector
   ratio and grouped cycling as "2–4 inserts 37 °C 1 hr" and "5–10 inserts × 30 cycles". The
   current manual says **1:1** for amplicons and "2–10 inserts × 30 cycles". The archived
   *Insert Considerations* page (**archived, 2020-11-29**) also says 2:1 against pGGA (2,174 bp).
   Our code should follow the current manual and record which revision it followed.
2. The destination plasmid changed from pGGA (2,174 or 2,155 bp depending on the page) to
   pGGAselect (2,220 bp per the E1601 manual text, 2,155 bp per the E1602 manual insert
   section). NEB's own pages disagree on this number; do not hard-code it.

## 3. Vector prepared by PCR

- **Outward primers carrying the site.** NEB's rule is orientation, not sequence: "the
  recognition sites should always face inwards towards your DNA to be assembled" (Technical
  Tips, archived 2020-10-26). For a PCR-linearised vector the same rule applies to the
  backbone: each outward primer carries spacer bases, then the recognition site pointing back
  into the backbone, then the 4-nt (or 3-nt) overhang.
- **Bases needed 5′ of the site.** NEB, *Cleavage Close to the End of DNA Fragments*
  (**archived, 2021-04-16**), verbatim: "As a general rule and for enzymes not listed below, 6
  base pairs should be added on on either side of the recognition site to cleave efficiently.
  The extra bases should be chosen so that palindromes and primer dimers are not formed. In
  most cases there is no requirement for specific bases." The kit manuals repeat this as "5´
  flanking bases (6 recommended)". The per-enzyme table on that page (which enzymes cut with
  1, 2, 3… bp from the end) did not survive the snapshot — see [Open gaps](#10-open-gaps).
- **DpnI treatment of the template.** DpnI is the right tool because it only cuts methylated
  DNA: REBASE records DpnI as recognising `G6mA^TC`, a "Type II Methyl-directed restriction
  enzyme (recognizes the sequence only when methylated)". Plasmid template grown in a
  Dam-positive *E. coli* is cut; the PCR product, unmethylated, is not. NEB's own Golden Gate
  pages do **not** prescribe DpnI, so this step is ours to specify and to justify by the
  REBASE record, not by an NEB protocol.
- **Cleanup.** NEB is explicit that a spin column is required for multi-fragment work and
  optional-but-better for single insert: "While single insert amplicons can be used directly
  from PCR without purification under certain circumstances (see FAQ), it is always best to
  purify amplicons using spin columns". The reason to purify, from FAQ 9: carryover PCR
  polymerase fills in 5′ four-base overhangs, producing blunt ends and non-specific assembly.
  Unpurified amplicon is tolerated for single-insert cloning only "as long as insert volume is
  1 µl or less".
- **Storage.** "For long term storage at –20°C, store DNA in 10 mM Tris (pH 8.5), 1 mM EDTA
  (TE) or short-term storage in 10 mM Tris (pH 8.5), 0.1 mM EDTA" — EDTA at those levels does
  not meaningfully lower the 10 mM MgCl₂ of the ligase buffer (E1601 manual, FAQ 6).
- **Polymerase and cycle count.** "Use a high-fidelity DNA polymerase and avoid
  over-amplification… use the minimum number of cycles required… this is usually 20 cycles or
  less" (E1601 manual, FAQ 11).

## 4. Overhangs and ligation fidelity

### 4.1 The classic design rules, and their measured limits

Pryor 2020 states the conventional rules and then shows they are the wrong abstraction:
"avoiding use of palindromic overhang sequences or the same overhang pair more than once in an
assembly reaction. In addition, most modular cloning systems also require that
non-complementary overhang sequences have at least 2 mispaired bases. Moreover, overhangs that
contain 100% A/T or G/C content are also often avoided… sometimes, overhangs with 25% G/C
content are also avoided for the same reason."

Potapov 2018 measured that the two-mismatch rule is unnecessarily strict: "it is not necessary
to ensure all overhangs have at least two bases different from all other overhangs, as many
pairs with only a single base difference … form very few if any mismatch ligation products with
each other." And Pryor 2020 measured that efficiency is not a function of composition: "the
relative efficiency of each Watson-Crick pair was not simply a function of GC content, and
thus, difficult to predict based on the sequence composition alone. For example, … the
5′-ATTT/5′-AAAT overhang pair is significantly more efficient than assembly of the
5′-TTTA/5′-TAAA pair."

Two empirical findings that a scorer should not lose:

- A 100% GC junction can fail even when the fidelity data says it is fine. Potapov 2018 saw
  truncations at junction 6, the `GCCG/CGGC` pair, with "a drop of roughly 23% of connections
  at junction 6 relative to the average for all junctions, despite this sequence not being a
  low efficiency junction in the ligation fidelity data set."
- A high-fidelity but low-efficiency pair (`TAAA/ATTT`) joined with about 30% reduced incidence
  and produced truncated products. Fidelity and efficiency are two different axes.

### 4.2 Potapov et al. 2018, *ACS Synth. Biol.* 7, 2665–2674

- DOI 10.1021/acssynbio.8b00333. Method: all 256 four-base overhangs on one multiplexed
  substrate, ligated then read by PacBio SMRT sequencing.
- Conditions: substrate 100 nM, 1.75 µM T4 (or T7) DNA Ligase, 1X T4 DNA Ligase Buffer
  (50 mM Tris-HCl pH 7.5, 10 mM MgCl₂, 1 mM ATP, 10 mM DTT), 50 µl, "incubated for 1 or 18 h
  at 25 or 37 °C". A cycled Golden Gate condition was also run: "30 cycles of 1 or 5 min at
  37 °C and 1 or 5 min at 16 °C, followed by a 5 min 65 °C heat inactivation step."
- Counts are normalised to 100,000 ligation events; a typical Watson–Crick pair then sits at
  300–400 observations.
- Validation: a 10-fragment assembly with the predicted high-fidelity (HF) set gave "99.9% of
  all observed assemblies … formed by only correct Watson−Crick pairings", against the
  low-fidelity (LF) set where "74.1% of assemblies containing at least one erroneous junction".
  Predicted-versus-observed for lac cassettes: 99% (12-fragment) and 91% (24-fragment)
  predicted; about 99% blue observed for the 12-fragment set, 45 ± 5% for the deliberately
  deletion-prone set (predicted about 33%).
- **Ready-made high-fidelity overhang sets**, Table 1, "for use with cycled assembly
  (16 °C/37 °C cycles)", one member of each complementary pair shown:

  | Set | Overhangs | Estimated fidelity | Sequences |
  | --- | --- | --- | --- |
  | 1 | 15 | 98.5% | TGCC, GCAA, ACTA, TTAC, CAGA, TGTG, GAGC, AGGA, ATTC, CGAA, ATAG, AAGG, AACT, AAAA, ACCG |
  | 2 | 20 | 98.1% | AGTG, CAGG, ACTC, AAAA, AGAC, CGAA, ATAG, AACC, TACA, TAGA, ATGC, GATA, CTCC, GTAA, CTGA, ACAA, AGGA, ATTA, ACCG, GCGA |
  | 3 | 25 | 95.8% | CCTC, CTAA, GACA, GCAC, AATC, GTAA, TGAA, ATTA, CCAG, AGGA, ACAA, TAGA, CGGA, CATA, CAGC, AACG, AAGT, CTCC, AGAT, ACCA, AGTG, GGTA, GCGA, AAAA, ATGA |
  | 4 | 30 | 91.7% | TACA, CTAA, GGAA, GCCA, CACG, ACTC, CTTC, TCAA, GATA, ACTG, AACT, AAGC, CATA, GACC, AGGA, ATCG, AGAG, ATTA, CGGA, TAGA, AGCA, TGAA, ACAT, CCAG, GTGA, ACGA, ATAC, AAAA, AAGG, CAAC |

  Set 1 extends the MoClo standard set (TGCC, GCAA, ACTA, TTAC, CAGA, TGTG, GAGC). Sets 2–4
  "exclude all Watson−Crick pairs that ligate significantly below the mean" and avoid 100% GC
  overhangs. Subsets are predicted to be at least as good as the whole set. A separate table
  (S8) covers static 37 °C incubation.
- Supporting Data, figshare item 7267505, `sb8b00333_si_002.zip` (980,557 bytes, CC BY-NC 4.0),
  contains one file per condition plus the test sequence:

  `FileS01_T4_01h_25C.xlsx`, `FileS02_T4_01h_37C.xlsx`, `FileS03_T4_18h_25C.xlsx`,
  `FileS04_T4_18h_37C.xlsx`, `FileS05_HF_cycled.xlsx`, `FileS06_T7_18h_25C.csv`,
  `FileS07_LF_cycled.xlsx`, `FileS08_T7_18h_37C.csv`, `FileS09_DP_cycled.xlsx`,
  `FileS10_FP_cycled.xlsx`, `FileS11_HF_01h_37C.xlsx`, `FileS12_LF_18h_37C.xlsx`,
  `FileS13_DP_18h_37C.xlsx`, `FileS14_FP_18h_37C.xlsx`, `lac.fasta`.

  The four `T4_*` files are 257 × 257 count matrices (`Overhang` header, 256 overhangs on each
  axis). The `*_cycled` files are per-assembly summaries instead: sheets `table_01`…`table_05`
  holding insert counts and fractions, a small overhang-versus-overhang junction matrix, and an
  assembly-composition table such as `{A:B:C:D:E:F:G:H`, count, correct flag, size, fraction.

### 4.3 Potapov et al. 2018, *Nucleic Acids Research* 46, e79 (gky303)

The companion method paper, three-base overhangs, CC BY-NC 4.0. Conditions: 1.75 µM T4 DNA
Ligase, 100 nM randomised three-base 5′-overhang substrate, 1X T4 ligase buffer, 25 °C or
37 °C, 1 h or 18 h. Fidelity there is "the fraction of correct ligations divided by the total
fraction of ligations for a given overhang". Data: SRA `SRP130363`; analysis code
`github.com/potapovneb/ligase-fidelity` (**AGPL-3.0** — read it, do not vendor it).

### 4.4 Pryor et al. 2020, *PLoS One* 15(9): e0238592 — the shippable dataset

This is the one to ship. It is CC BY 4.0, it is per-enzyme, and it covers SapI's three-base
overhangs, which Potapov 2018 does not.

Conditions, verbatim from Methods: reactions of 20 µl, "with T4 DNA ligase and BsaIHF-v2 or
BsmBI-v2 … using their respective NEB Golden Gate Enzyme Mixes (2 μL) in 1X T4 DNA ligase
buffer"; "T4 DNA ligase (500 U) and SapI (15 U) or Esp3I (15 U) … in 1X T4 DNA ligase buffer";
"T4 DNA ligase (500 U) and BbsI-HF (15 U) … in 1X CutSmart Buffer supplemented with 10 mM DTT
and 1 mM ATP". Substrate at 100 nM. "The reactions were cycled between 37°C and 16°C (SapI,
Esp3I, BsaI-HFv2, BbsI-HF) or 42°C and 16°C (BsmBI-v2) for 5 minutes at each temperature for
30 cycles, and then subjected to a final heat-soak for 5 minutes at 60°C."

Fidelity, verbatim: "Fidelity F for a set of n overhangs {O₁, O₂, O₃, …, Oₙ} was defined as a
probability that all overhangs in the set ligate correctly to their WC pair… computed as
follows: F = p(O₁) × p(O₂) × p(O₃) × p(Oₙ)", where "p(Oᵢ) = N_correct / N_total". GetSet
searches for sets with an MCMC optimiser, because exhaustive search "exceeds 10¹⁴ combinations
for 10-overhang sets"; consequently "repeating a search can result in different junction sets
with similar predicted fidelities".

Files, and exactly what each is:

| Label | File | Bytes | Contents |
| --- | --- | --- | --- |
| S1 Table | `pone.0238592.s001.xlsx` | 202,549 | BsaI-HFv2 ligation frequencies, sheet `S1 Table. BsaI-HFv2`, range `A1:IW257` |
| S2 Table | `pone.0238592.s002.xlsx` | 203,006 | BsmBI-v2, same shape |
| S3 Table | `pone.0238592.s003.xlsx` | 202,302 | Esp3I, same shape |
| S4 Table | `pone.0238592.s004.xlsx` | 203,060 | BbsI-HF, same shape |
| S5 Table | `pone.0238592.s005.xlsx` | 24,232 | SapI, sheet `Table S5. SapI`, range `A1:BM65` (64 three-base overhangs) |
| S6 Table | `pone.0238592.s006.xlsx` | 13,251 | 13-fragment *lac* test system: `Fragment #`, `Sequence`, `Length (bp)` |
| S7 Table | `pone.0238592.s007.xlsx` | 13,808 | 35-fragment *lac* test system |
| S8 Table | `pone.0238592.s008.xlsx` | 10,494 | DNA substrate precursor oligonucleotides |
| S1–S3 Fig | `.s009`–`.s011` | TIFF | Mispair partners; colony PCR gels; fidelity versus set size |
| S1 Text | `pone.0238592.s012.pdf` | 79,081 | Error propagation analysis |

Download URL pattern (302-redirects to a signed Google Storage URL, `curl -L` works):

```text
https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0238592.sNNN&type=supplementary
```

Matrix format, confirmed by opening the files: row 1 is a header, cell `A1` is the literal
`Overhang`, then 256 (or 64) overhang labels across; column A holds the overhang labels down.
Cells are **integer observation counts**, not probabilities. First rows of S1:

```text
Overhang,AAAA,AAAC,AAAG,AAAT,AACA,AACC,...
TTTT,635,8,40,16,2,0,...
GTTT,3,476,4,45,0,20,...
```

Axis orientation comes from the Ligase Fidelity Viewer help page: "Overhangs in rows correspond
to top strand. Overhangs in columns correspond to bottom strand. Overhangs are always written
5' to 3'." So `ATTC` (row) against `GAAT` (column) is the Watson–Crick pair. Note the diagonal
is not the Watson–Crick entry — `TTTT`/`AAAA` is, i.e. the pairing is row against the reverse
complement of the column label.

### 4.5 What NEB's tools compute

From the Ligase Fidelity Viewer help page (`tools.neb.com/~potapov/ligase-fidelity-viewer/help.html`,
reachable today):

- Input: a comma-separated list of overhangs, 5′→3′, four canonical bases, no degenerate bases;
  "Complementary (reverse) overhangs are automatically added to the query"; duplicates and
  malformed entries are discarded.
- Condition choice selects which published dataset scores the set. "The default conditions are
  ligation at 25°C for 18 hours; these conditions have been shown to well predict the results
  of Golden Gate assembly using typical cycled conditions (16°C 5 min/37°C 5 min, 30 cycles)."
- Output: "the total fraction of correct ligation events out of all ligation events", plus a
  matrix with thresholds we can copy directly: Watson–Crick pairs above 100 normalised
  ligations are "strong", below 100 "weak" (users may wish to avoid them); mismatches below 10
  are "trace", 10–100 "modest … should be avoided", above 100 "high-count … observed on a
  frequency similar to that of proper Watson-Crick ligations".
- The three tools, per the E1601/E1602 manuals and Pryor 2020: **Ligase Fidelity Viewer**
  (score a set), **GetSet** (generate or extend a high-fidelity set, with MCMC), **SplitSet**
  (split a target sequence into high-fidelity fragments). Worked example from Pryor 2020: the
  11 standard plant-synthetic-biology overhangs (GGAG, TGAC, TCCC, TACT, CCAT, AATG, AGCC,
  TTCG, GCTT, GGTA, CGCT) score 81% with BsmBI-v2 and 42 °C/16 °C cycling; GetSet extended the
  set to 20 overhangs at 80%.
- The NEBridge Golden Gate Assembly Tool (`goldengate.neb.com`) designs primers: per the
  E1602 manual, "The primers will feature 6 bases at the 5´ end flanking the BsmBI recognition
  site, the recognition site itself, plus the 4-base overhangs… All overhangs will
  automatically be designed as non-palindromic (to eliminate self insert ligations), unique,
  and in the correct orientations."

## 5. Transformation and screening

From the E1601 and E1602 manuals (identical text):

1. Thaw 50 µl NEB 10-beta Competent *E. coli* (C3019) on ice for 10 min.
2. Add **2 µl** of the assembly reaction; flick 4–5 times.
3. Ice 30 min.
4. Heat shock 42 °C for 30 s.
5. Ice 5 min.
6. Add 950 µl room-temperature NEB 10-beta/Stable Outgrowth Medium; 37 °C for 60 min at
   250 rpm.

Plating: warm selective plates to 37 °C; "spread 50 µl of a 1:5 dilution (single inserts) or
50–100 µl (multiple inserts) of the 1 ml outgrowth onto each plate"; overnight at 37 °C, or
24–36 h at 30 °C, or 48 h at 25 °C.

Strains (E1601 manual, FAQ 15): 10-beta (C3019) or NEB Stable (C3040) for large or repeat-rich
assemblies; 5-alpha (C2987), Turbo (C2984) or T7 Express (C2566) for smaller ones;
"Subcloning efficiency cells will result in lower transformation levels and should not be used
for multi-component assemblies."

**Expected counts.** NEB's release specification for E1601 is a usable yardstick: a 5-insert
lacI/lacZ assembly transformed into T7 Express "yields greater than 250 colonies and > 80% blue
colonies when 5% of the outgrowth is spread on an IPTG/Xgal/Chloramphenicol plate and incubated
overnight at 37°C" (PS-E1601S/L v1.0; E1602's is the same shape with BsmBI-v2 at 42 °C, and
M1100's uses BsaI-HFv2 in a 15 µl reaction). The technical note gives a fuller efficiency
ladder for BsaI-HFv2 with T4 DNA Ligase:

| Fragments | Protocol (outgrowth plated) | Correct colonies/plate | Fidelity | Per full reaction |
| --- | --- | --- | --- | --- |
| 1 | 5 min 37 °C (2.5 µl) | 687 | 100% | 2,742,000 |
| 1 | 60 min 37 °C (2.5 µl) | 1,623 | 100% | 6,492,000 |
| 12 | (5 min 37 °C → 5 min 16 °C) × 30 (5 µl) | 245 | 99.5% | 489,000 |
| 24 | same × 30 (100 µl) | 78 | 90.7% | 9,792 |

NEB's 52-fragment brochure figure (Pryor et al. 2022, *ACS Synth. Biol.* 11, 2036–2042) reports
an average 727 correct and 733 incorrect colonies per 100 µl outgrowth, i.e. 49% correct.

**Blue/white.** Two opposite conventions exist and confusing them will mis-read a plate:

- NEB's own kit test systems use **reverse** blue-white: the assembly reconstructs *lacZ*, so
  **correct assemblies are blue**, and NEB deliberately uses the non-complementing T7 Express
  strain "to avoid any possibility of alpha-fragment LacZ complementation" (technical note,
  Figure 2 caption). Potapov 2018 gives plate composition for this: "Luria−Bertani broth
  supplemented with 1 mg/mL dextrose, 1 mg/mL MgCl2, 30 μg/mL Chloramphenicol, 200 μM IPTG and
  80 μg/mL X-gal", read after 18 h at 37 °C plus 8 h at 4 °C.
- Classic pUC19 screening is the other way round: an insert in the MCS disrupts *lacZα*, so
  **correct clones are white**. This is our smoke test (see below).

**Controls.** NEB's position: "Golden Gate assembly protocols do not usually call for a
negative control. However if desired, a 'no insert(s) added' reaction can be used" (FAQ 13).
The technical note describes the stronger control it used for validation: omit any single
component and no blue colonies appear.

**Screening protocol.** Colony PCR with OneTaq (M0480) or OneTaq Hot Start (M0481); Taq
(M0267) is acceptable; LongAmp (M0323/M0534) "strongly recommended" for large assemblies.
NEB's own pGGAselect screening primers, which are also its sequencing primers, are given in
the manuals: forward 65 bp upstream `5´-CTGCAGGAAGGTTTAAACGCATTTAGG-3´`, reverse 62 bp
downstream `5´-TAATACGACTCACTATAGGGAGACTC-3´`. Whatever the screen, "the correct assembly of
insert(s) should always be confirmed by sequencing of the plasmid construct across the 4 base
junctions and inserts."

## 6. Troubleshooting, from NEB

Collected from the E1601/E1602 manual FAQs and tips, the Technical Tips page
(**archived, 2020-10-26**) and the 2024 brochure's eleven tips.

| Symptom or step | NEB's guidance |
| --- | --- |
| Internal site in a part | Single-insert assembly usually survives one internal site; multi-insert does not. Switch enzyme, domesticate by site-directed mutagenesis (Q5 SDM kit E0554, NEBaseChanger), or "a junction point can be created at the internal site's recognition sequence". |
| No enzyme is free | Use a 7-base-recognition enzyme (PaqCI) — "less likely to have internal sites present in any given sequence" — or switch to NEBuilder HiFi if ≤5 inserts. |
| Background of empty vector | Keep the 60 °C end soak; it digests uncut and re-ligated destination plasmid. |
| Low transformant numbers | Increase cycles from 30 to 45–65; plate more of the outgrowth (NEB plated 2.5 µl for 1 fragment, 100 µl for 24). |
| Mis-assemblies from amplicons | Purify amplicons; primer dimers carrying Type IIS sites "would be active in the assembly reaction and result in mis-assemblies". |
| Blunt-end products from unpurified PCR | Carryover polymerase fills in 5′ overhangs; purify, or keep unpurified insert to ≤1 µl and single-insert only. |
| Over-large plasmid | "Efficiencies are highest with assembled product plasmid constructs ~ 10-12 kb"; larger needs more colonies screened. |
| Repetitive or very small/large parts | Preclone inserts <250 bp or >3 kb, or repeat-rich ones; use NEB Stable cells for unstable repeats. |
| Plasmid concentration overestimated | Make sure precloned-insert plasmid preps are RNA-free. |
| Previously working part stops working | Suspect a propagation error in *E. coli*, "usually a frameshift due to slippage in a run of a single base (e.g., AAAA)". |
| Complex assembly, >10 fragments | Insert amounts can drop from 75 to 50 ng each without hurting efficiency. |
| Reaction volume too small for the DNA | Scale the 20 µl reaction up to 25 µl (+0.5 µl 10X buffer), or scale everything down 2–3× if volumes allow. |
| PCR errors in inserts | Q5 High-Fidelity polymerase, ≤20 cycles. |

## 7. Other assembly methods, briefly

| Method | Mechanism | When it fits | What a later skill needs from us |
| --- | --- | --- | --- |
| NEBuilder HiFi / Gibson | 5′ exonuclease chews back ends, overlapping homology anneals, polymerase and ligase finish it; one 50 °C incubation | NEB's own cut-off: "if the assembly will involve 5 or less inserts" (E1601 FAQ 4). No Type IIS site needed, so no domestication | Homology-arm design (length, Tm), junction-free product simulation. Reaction table still pending — see [Open gaps](#10-open-gaps) |
| In-Fusion (Takara) | Proprietary enzyme joins 15-bp overlaps; 50 °C, 15 min | Single or few inserts into a linearised vector, no scars, no Type IIS constraint | Primer design with a 15-bp vector-homology 5′ tail, 18–25 nt template-specific 3′ part; 20 bp when joining more than two fragments |
| MoClo (Weber et al. 2011, *PLoS ONE* 6:e16765, CC BY) | Hierarchical Golden Gate with alternating enzymes: BpiI at level 0 and level 2, BsaI at level 1 | Multigene constructs built from a shared, standard part library | Standard fusion sites (AATG at the start codon, AGGT for a signal-peptide junction, GGAG/TACT/GCTT/CGCT in non-translated regions; position-specific sites such as TGCC/GCAA), and its cycling: 5 h 37 °C → 50 °C 5 min → 80 °C 10 min, or 45 × (37 °C 2 min → 16 °C 5 min) for hard constructs |
| GoldenBraid (Sarrion-Perdigones et al. 2011, *PLoS ONE* 6:e21622, CC BY) | Double loop: BsaI at level α, BsmBI at level Ω, alternating indefinitely; BbsI used for domestication | Iterative binary assembly where composite parts re-enter the other level | Level-aware enzyme choice, numbered BsaI sites versus lettered BsmBI sites, 25 digestion/ligation cycles |

Both standards matter for enzyme choice: NEB's own BbsI-HF and SapI pages list "Expanded
'assembly standards' for MoClo, GoldenBraid2.0 and other modular Golden Gate Assembly methods"
as a product note (**archived, 2021-04**).

## 8. Notes for the smoke test (pUC19 + GFP)

NEB's pUC19 GenBank record (`neb.com/-/media/.../puc19gbk.txt`, fetched today) gives the
coordinates the fixture should reproduce: `LOCUS pUC19 2686 bp DNA circular 18-OCT-2007`,
*lacZα* `complement(146..469)`, MCS `misc_feature 396..452` ("multiple cloning site
(EcoRI-HindIII)"), lac promoter signals `complement(514..519)` and `complement(538..543)`,
CAP site `complement(563..575)`, pMB1 origin `complement(867..1455)`, *bla* (AmpR)
`complement(1626..2486)`.

Consequences for the protocol the skill writes:

- **Colony colour: white.** GFP inserted into the MCS interrupts the *lacZα* reading frame, so
  correct clones are white on X-gal/IPTG and empty vector is blue. This is the opposite of
  NEB's kit test systems, which score correct assemblies as blue. Say which convention the
  plate is using.
- **Host must α-complement.** Blue/white with pUC19 needs a *lacZΔM15* host such as NEB
  5-alpha (C2987). NEB's own Golden Gate test systems deliberately use T7 Express, whose
  genotype carries `lacZ::T7 gene1` and therefore cannot complement (Potapov 2018, Methods) —
  a fine choice for their reverse screen and the wrong choice for ours.
- **Plate composition.** A citable recipe for the indicator plate is Potapov 2018's: LB with
  1 mg/mL dextrose, 1 mg/mL MgCl₂, 200 µM IPTG, 80 µg/mL X-gal, plus the relevant antibiotic —
  chloramphenicol there, ampicillin or carbenicillin for pUC19. Read after overnight 37 °C;
  their protocol adds 8 h at 4 °C before scoring colour.
- **Enzyme choice.** With BsaI blocked by a site in AmpR and BsmBI/Esp3I blocked by two sites
  in pUC19, the free enzymes in the fixture pair are BbsI and PaqCI. Both have NEB reaction
  tables via M1100 (BbsI-HF: 1 µl/20 U, 37 °C; PaqCI: 1 µl/10 U at 37 °C plus 0.5 µl activator),
  so either yields a fully cited protocol. SapI has one site in pUC19 (683) and gives 3-nt
  overhangs; prefer BbsI or PaqCI.
- **Single-insert program.** Two fragments, single-gene cloning: with the kit, 37 °C 5 min →
  60 °C 5 min; with M1100 and BbsI-HF or PaqCI, 37 °C 15 min → 60 °C 5 min. Use 15 cycles of
  42 °C/16 °C only if the chosen enzyme is BsmBI-v2 or BspQI.
- **No fluorescence expected.** The Clontech GFP CDS fixture has no ribosome binding site of
  its own, and in the MCS it sits downstream of the lac promoter on the reverse strand; the
  protocol should say the clone is not expected to fluoresce and that colour scoring is about
  *lacZα*, not GFP.
- **Colony PCR.** Per the fixture description in issue #1 (re-derive it in tests rather than
  trusting it here), M13 forward (379–395) and M13 reverse (465–481) flank the MCS, giving a
  103 bp empty-vector amplicon, so the correct-clone band is 103 bp + insert (plus tails). NEB's
  pGGAselect primers above are the wrong primers for pUC19 — do not reuse them.

## 9. Implications for this package

**For #7 (enzyme data).** The per-enzyme table in §2.6 is the field list: name, commercial
name, recognition site in IUPAC, top/bottom cut offsets in REBASE `(n/m)` form, overhang length
and end type (5′), incubation temperature, Golden Gate temperature (which differs from the
incubation temperature for BsmBI-v2 and BspQI, and from the unit-definition temperature for
BsmBI-v2 and BtgZI), heat-inactivation temperature and time, Dam/Dcm/CpG sensitivity as one of
NEB's seven states (not sensitive, impaired, blocked, and the four "by overlapping" /"by some
combinations of overlapping" variants), isoschizomers, prototype, and a
`requires_activator` flag for PaqCI. Two values need a note field rather than a number:
BsmBI-v2's Dam/Dcm cells were not readable in the archived chart, and BtgZI has no NEB Golden
Gate protocol at all.

**For #11 (overhang design).** Ship Pryor 2020 S1–S5 (CC BY 4.0) as package data under
`src/liulab_mbio/data/`, one file per enzyme, with a `LICENSE`/attribution file naming the
paper, DOI and CC BY 4.0, and a regeneration script that re-downloads from the
`journals.plos.org` URL pattern in §4.4. Convert the xlsx to a compact tabular form but keep
the integer counts, since the fidelity formula needs `N_correct / N_total`. Do **not** ship
Potapov 2018's ACS data (CC BY-NC) — it is the only source for T7 ligase and for static 25 °C
conditions, so cite it and offer the Potapov Table 1 overhang sets as literature constants
(short factual lists, reproduced with citation) rather than bulk data. Implement fidelity as
`F = Π p(Oᵢ)` with `p(Oᵢ) = N_correct / N_total` over the chosen set, matched to the enzyme's
own matrix and thermocycling temperature. Score sets with the viewer's thresholds (weak WC pair
below 100, mismatch trace/modest/high at 10 and 100 normalised counts). Keep the classic rules
as cheap pre-filters (non-palindromic, unique, no reuse of a pair) but let the data override
the "at least two mismatched bases" and GC-content heuristics, and flag 100%-GC junctions as
truncation-prone on Potapov's junction-6 evidence. Remember SapI/BspQI need the 64-overhang
three-base matrix, not the 256-overhang one, and that the matrix axes are top strand versus
bottom strand.

**For #13 (bench numbers).** Every table in §2 is quotable with a citation. Pick one system
per protocol and say which: kit E1601/E1602 (20 µl, 0.05 pmol vector, 1:1 amplicon:vector,
1–2 µl enzyme mix) or M1100 (15 µl for 2 to 6 fragments, 30 µl for 7+, 0.05 pmol each, enzyme
volume from the table, PaqCI activator when applicable). Always end with the 60 °C 5 min soak,
then 2 µl into 50 µl cells. Keep the ng↔pmol conversion in code rather than quoting NEB's
75 ng figure, which is specific to a 2,155–2,220 bp plasmid. Record the manual revision
(E1601 v5.0_6/26) next to the numbers, because NEB changed the amplicon ratio from 2:1 to 1:1
between revisions.

**For #14 (the skill).** The HTML protocol should print: the enzyme's own Golden Gate
temperature and program tier, the reaction table for whichever system the user has, the
activator line only for PaqCI, the end soak, the transformation and plating steps with the
strain requirement (α-complementing host for pUC19 blue/white), the expected colony colour
with the convention named explicitly, expected colony counts from §5 as an order-of-magnitude
check, the colony PCR band sizes from the fixture's own M13 primers, and the per-step
troubleshooting rows in §6. Cite NEB per step; the skill should never invent a number that is
not in this note or computed by the package.

## 10. Open gaps

Still unverified. None of these is guessed anywhere above.

| Item | Why it is missing | Workaround in use |
| --- | --- | --- |
| Per-enzyme rows of *Cleavage Close to the End of DNA Fragments* | Table did not render in the Wayback snapshot; live page 403 | NEB's general rule of 6 bp either side, quoted in §3 |
| NEBuilder HiFi reaction table and overlap length | `manuale2621.pdf` returns 403 to both `curl` and WebFetch; no snapshot fetched | Method described qualitatively in §7; NEB's ≤5-insert cut-off cited from the E1601 manual |
| BtgZI product page | No Wayback snapshot under either URL spelling | Specification PDF (unit definition at 60 °C), heat-inactivation and methylation charts |
| BsmBI-v2 Dam and Dcm cells | Archived chart row truncated on parse | CpG "Blocked" from the archived product page; Dam/Dcm left unstated |
| Ligase Fidelity Viewer **v2** help, `ligasefidelity.neb.com` | 403 live; only the older `tools.neb.com` v1 help page is reachable | v1 help page, quoted in §4.5; the Pryor 2020 paper describes the v2 condition list |
| `goldengate.neb.com` tool behaviour | 403 live, no snapshot fetched | The E1601/E1602 manuals describe what the tool designs, quoted in §4.5 |
| NEB screening-protocol and 52-fragment protocol pages | 403 live | Manual screening section and the brochure's 52-fragment figure |
| REBASE licence or citation terms | No statement on the enzyme pages read | Flagged for #7; check REBASE's own terms page before shipping REBASE-derived data |

## Sources

All retrieved 2026-09-12.

- NEB, *NEBridge Golden Gate Assembly Kit (BsaI-HFv2)* instruction manual, NEB #E1601S/L,
  version 5.0_6/26: [manuale1601.pdf](https://www.neb.com/-/media/nebus/files/manuals/manuale1601.pdf)
- NEB, *NEBridge Golden Gate Assembly Kit (BsmBI-v2)* instruction manual, NEB #E1602S/L:
  [manuale1602.pdf](https://www.neb.com/en/-/media/nebus/files/manuals/manuale1602.pdf)
- NEB product specifications (fetched under `/-/media/catalog/specifications/`): PS-E1601S/L
  v1.0, PS-E1602S/L v1.0, PS-M1100S/L v1.0, PS-R3733S/L v1.0, PS-R0739S/L v1.0, PS-R3539S/L
  v1.0, PS-R0745S/L v1.0, PS-R0569S/L v1.0, PS-R0703S/L v1.0, PS-R0734S/L v1.0
- NEB, *Protocol for NEBridge Ligase Master Mix (NEB #M1100)* — archived snapshot 2023-03-31:
  [web.archive.org](https://web.archive.org/web/20230331002719id_/https://www.neb.com/protocols/2021/09/14/protocol-for-nebridge-ligase-master-mix-neb-m1100)
- NEB, *NEBridge Ligase Master Mix Protocol Guidelines* — archived snapshot 2025-07-13:
  [web.archive.org](https://web.archive.org/web/20250713174145id_/https://www.neb.com/en-us/tools-and-resources/usage-guidelines/nebridge-ligase-master-mix-protocol-guidelines)
- NEB, *Usage Guidelines for Golden Gate Assembly with PaqCI* — archived snapshot 2021-06-15:
  [web.archive.org](https://web.archive.org/web/20210615031818id_/https://www.neb.com/tools-and-resources/usage-guidelines/usage-guidelines-for-golden-gate-assembly-with-paqci)
- NEB, *Cleavage Close to the End of DNA Fragments* — archived snapshot 2021-04-16;
  *Heat Inactivation* — 2021-04-20; *Dam-Dcm and CpG Methylation* — 2021-04-13;
  *Technical Tips For Optimizing Golden Gate Assembly Reactions* — 2020-10-26;
  *Insert Considerations … (NEB #E1601)* — 2020-11-29 (all via `web.archive.org`)
- NEB product pages, archived snapshots 2021-04: BsaI-HFv2 (R3733), BsmBI-v2 (R0739),
  BbsI-HF (R3539), SapI (R0569), PaqCI (R0745)
- NEB, *Golden Gate Assembly: 50+ fragment assembly now achievable* brochure, NEB142_V4.1_1024,
  and *Breaking through the Limitations of Golden Gate Assembly* technical note, 01/19, both via
  the `neb-online.de` mirror
- NEB, pUC19 GenBank record: [puc19gbk.txt](https://www.neb.com/-/media/nebus/page-images/tools-and-resources/interactive-tools/dna-sequences-and-maps/text-documents/puc19gbk.txt)
- NEB, *Ligase Fidelity Viewer: Help Page*:
  [tools.neb.com](https://tools.neb.com/~potapov/ligase-fidelity-viewer/help.html)
- Potapov, V. et al. (2018) Comprehensive profiling of four base overhang ligation fidelity by
  T4 DNA Ligase and application to DNA assembly. *ACS Synth. Biol.* 7, 2665–2674.
  [doi:10.1021/acssynbio.8b00333](https://doi.org/10.1021/acssynbio.8b00333). Supporting Data
  via figshare collection 4283045 (CC BY-NC 4.0)
- Potapov, V. et al. (2018) A single-molecule sequencing assay for the comprehensive profiling
  of T4 DNA ligase fidelity and bias during DNA end-joining. *Nucleic Acids Res.* 46, e79.
  [doi:10.1093/nar/gky303](https://doi.org/10.1093/nar/gky303) (CC BY-NC 4.0)
- Pryor, J.M. et al. (2020) Enabling one-pot Golden Gate assemblies of unprecedented complexity
  using data-optimized assembly design. *PLoS One* 15(9): e0238592.
  [doi:10.1371/journal.pone.0238592](https://doi.org/10.1371/journal.pone.0238592) (CC BY 4.0)
- Pryor, J.M. et al. (2022) *ACS Synth. Biol.* 11, 2036–2042 (52-fragment assembly; cited via
  the NEB brochure figure)
- Kennedy, M.A., Hosford, C.J., Azumaya, C.M., Luyten, Y.A., Chen, M., Morgan, R.D. and
  Stoddard, B.L. (2023) Structures, activity and mechanism of the Type IIS restriction
  endonuclease PaqCI. *Nucleic Acids Res.* 51, 4467–4487.
  [doi:10.1093/nar/gkad228](https://doi.org/10.1093/nar/gkad228) (CC BY)
- Weber, E. et al. (2011) A modular cloning system for standardized assembly of multigene
  constructs. *PLoS ONE* 6: e16765 (CC BY)
- Sarrion-Perdigones, A. et al. (2011) GoldenBraid: an iterative cloning system for
  standardized assembly of reusable genetic modules. *PLoS ONE* 6: e21622 (CC BY)
- REBASE enzyme records for BsaI, BsmBI, Esp3I, BbsI, SapI, PaqCI, BtgZI and DpnI:
  [rebase.neb.com](http://rebase.neb.com/rebase/rebase.html)
- Takara Bio, In-Fusion HD cloning user manuals and *In-Fusion Cloning* learning centre pages

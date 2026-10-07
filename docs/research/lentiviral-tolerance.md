---
search:
  exclude: true
---

# What a lentiviral vector tolerates: a domestication map

Issue #244, under the AP-1 demo spec #225, blocking #239. Written for someone domesticating a
**different** backbone later: the question it answers is not "may we change pLVX" but "which
parts of a third-generation lentiviral transfer plasmid may be changed at all, and what does it
cost to be wrong".

Two kinds of claim live here and they are marked differently. **Measured** means a primary source
that changed the thing and reported a number, or a count run here on a real sequence. **Inferred**
means a chain of reasoning from something measured; every inference is labelled, because an
unmarked one in this note propagates into every future domestication. **Nothing measures this**
means what it says, and is followed by what would close it.

## 0. Read this first, if you are in a hurry

1. **Pick the enzyme before you pick the edit.** The HIV-1 LTR has no BsmBI, PaqCI or BtgZI
   site. It has one BsaI site and one BbsI site, in each copy. Section 1 is the table; changing
   enzyme is free and changing the LTR is not.
2. **BsaI's site in the LTR is not a quirk of one plasmid.** `GGTCTC` is the first six bases of
   R in the HIV-1 reference genome, which is to say the first six bases the vector transcribes.
   Every HIV-1-derived transfer vector inherits it, twice, unless someone has removed it.
3. **One group has already made the substitution, and nobody has measured what it cost.**
   VAMSyB killed the BsaI site with a single C→T at R+4, in both LTRs, and the vector still
   packages and transduces. No titre was reported, for that vector or any other with a changed
   R. Section 2 has the edit; section 3 says what experiment would price it.

## 1. Choose the enzyme first — the cheapest domestication is the one you skip

Every count in this section was run on **2026-10-06** against GenBank **K03455.1** (HIV-1 HXB2,
9,719 bp), the reference every HIV-1-derived transfer vector descends from, and against the
Addgene-verified sequence of pLVX-TetOne-Puro-GFP (Addgene 171123) already held under
`reference_docs/synthesis_and_assembly/working-vector/`. Both strands counted. Positions are
1-based here because GenBank is; the package's own coordinates are 0-based and half-open.

### Sites in the 634 bp HIV-1 LTR

| Enzyme | Site | In one LTR | Where | In the whole HXB2 genome |
| --- | --- | --- | --- | --- |
| **BsaI** | `GGTCTC` | **1** | offset 455, the first six bases of R | 4 (two of them the two LTRs) |
| **BbsI** | `GAAGAC` | **1** | offset 25, in U3 | 7 |
| BsmBI / Esp3I | `CGTCTC` | **0** | — | 1 (offset 8856, in *env*) |
| PaqCI / AarI | `CACCTGC` | **0** | — | 1 (offset 5055, in *pol*) |
| BtgZI | `GCGATG` | **0** | — | 0 |
| SapI / BspQI | `GCTCTTC` | **0** | — | 2 |

**This is the most reusable fact in the note.** The HIV-1 LTR is already clean for BsmBI, PaqCI,
BtgZI and SapI. A method that assembles with any of those never has to touch an LTR. A method
that assembles with BsaI meets one site per LTR — two per vector — and a method that assembles
with BbsI meets one per LTR as well, in U3 rather than R.

The same count on pLVX-TetOne-Puro-GFP reproduces the eight sites #230 reported, and adds the
rest of the panel:

| Enzyme | Sites in pLVX-TetOne-Puro-GFP | Positions (1-based) |
| --- | --- | --- |
| BsaI | 6 | 455, 3725, 5731, 6473, 7113, 8721 |
| BsmBI | 2 | 3877, 5637 |
| BbsI | 3 | 25, 6399, 6683 |
| SapI | 2 | 1107, 7638 |
| PaqCI | **0** | — |
| BtgZI | 5 | 833, 3122, 3936, 4432, 6181 |

Positions 455 and 7113 are the two LTR copies of the BsaI site; 25 and 6683 the two LTR copies
of the BbsI site. **PaqCI is clean on the whole plasmid** — zero sites in 9,895 bp. Whether the
method may move to PaqCI is a method decision and belongs to #239, not here; what this note
records is that the option exists and costs no base change at all.

### The BsaI site is HIV-1, not Clontech

`GGTCTC` sits at offset 455 of the HXB2 LTR and at offset 455 of **both** pLVX LTRs, and the
71 bp window around it is identical in all three. Over the 634 bp LTR, pLVX differs from HXB2 at
**nine** positions, none of them in R. So the site is inherited, not introduced.

**Inference, marked:** any transfer vector whose LTRs derive from HXB2 or NL4-3 without
deliberate recoding carries this site twice. That follows from the reference sequence and from
R being the part of the LTR a SIN deletion does not touch; it was not checked against a panel of
backbones here. What would close it: count `GGTCTC` in the LTRs of pLKO.1, lentiCRISPRv2, pCDH
and FUGW.

### What fixes the transcription start site

The TATA box `TATAA` lies at LTR offset 427 in both copies, and R begins at 455 — 28 bp
downstream, the canonical TATA-to-+1 spacing. That is an independent, on-sequence confirmation
that **the first `G` of `GGTCTC` is +1**, the first base transcribed and the 5' end of the
packaged genome. The BsaI site is therefore not merely *in* TAR; it *is* transcription
positions +1 to +6.

### One thing the sequence says about this parent that the map does not

Both pLVX LTRs are a full 634 bp and both carry an intact U3 enhancer-promoter: two NF-κB sites
(`GGGACTTTCC` at offsets 350 and 364), an Sp1 GC-box (offset 389) and the TATA box (offset 427).
A self-inactivating vector's 3' LTR has U3 deleted; this one does not. The two copies differ at
exactly **one** base in 634, offset 24, immediately 5' of the U3 BbsI site — `A` in the 5' copy,
`C` in the 3' copy, where HXB2 has `C`.

**This matters for the map, so it is stated rather than assumed.** In a true ΔU3 SIN vector a
site in U3 exists in **one** copy and a site in R or U5 exists in **two**. In this parent both
U3s are present, so a U3 site exists in two copies as well. Whoever applies this note to another
backbone must check which case they are in before counting copies.

**Nothing read here establishes whether pLVX-TetOne is intended to be self-inactivating.** The
measurement above is of the Addgene-verified sequence only. What would close it: the
manufacturer's own map or sequence, or a passage assay for LTR-driven transcription from the
integrated provirus.

### The same count on a second HIV-1 reference

Run on **2026-10-06** against GenBank **AF324493.2** (pNL4-3, 14,825 bp, provirus plus plasmid
backbone), the other parent most vectors descend from. The 20-mer that opens R,
`GGTCTCTCTGGTTAGACCAG`, occurs **exactly twice** — at 455 and 9530, the R start of each LTR.

Two strains, same offset, same 20 bases. But the U3 `GAAGAC` is **not** conserved: HXB2 has it in
both LTRs, pNL4-3 only in the 3' copy. **The R-region BsaI site is invariant across the two
references; the U3 BbsI site is strain-variable.** Anyone applying this note to a different
backbone should count BbsI themselves and may take BsaI-in-R as given.

## 2. Has someone already published a Type IIS-clean lentiviral backbone?

This was the ticket's highest-value question. The answer has two halves and the second one is
the useful one.

**One group has published a BsaI-domesticated lentiviral transfer vector — VAMSyB — and reported
no titre for it.** The edit is below, read out of the deposited sequence base by base. **No titre
measurement after such a mutation exists to cite anywhere**, so a lab that wants BsaI assembly
directly into a lentiviral backbone now has a precedent but still no price.

**BsmBI-clean lentiviral transfer vectors are in universal use, and nobody had to domesticate
anything to get them** — because, as section 1 measures, the HIV-1 LTR never had a BsmBI site.
Two existence proofs, both read 2026-10-06:

| Source | What it shows |
| --- | --- |
| Fonseca, J. P. *et al.* A toolkit for rapid modular construction of biological circuits in mammalian cells. *ACS Synth. Biol.* **8**, 2593–2606 (2019). [doi:10.1021/acssynbio.9b00322](https://doi.org/10.1021/acssynbio.9b00322), PMID 31686495 | The Mammalian ToolKit (MTK, El-Samad lab) ships lentiviral destination vectors (Addgene 123951, 194220). Its order is "BsmBI part domestication, BsaI transcriptional unit assembly, and a final BsmBI assembly". The lentiviral backbone enters at the **final BsmBI** step; BsaI only ever acts inside part-entry vectors. The full text never mentions the LTR, TAR or a site in the viral backbone |
| Sanjana, N. E., Shalem, O. & Zhang, F. Improved vectors and genome-wide libraries for CRISPR screening. *Nat. Methods* **11**, 783–784 (2014). [doi:10.1038/nmeth.3047](https://doi.org/10.1038/nmeth.3047), PMID 25075903 | lentiCRISPRv2 (Addgene 52961) clones sgRNAs by **BsmBI digestion of the intact lentiviral transfer vector**. The whole GeCKO/Brunello library world does Golden-Gate-style BsmBI ligation into an HIV-1 backbone. The documented v1→v2 changes are titre-related, not Type IIS domestication |

**EMMA** (Martella, A. *et al.* *ACS Synth. Biol.* **6**, 1380–1392 (2017),
[doi:10.1021/acssynbio.7b00016](https://doi.org/10.1021/acssynbio.7b00016), PMID 28418644, read
2026-10-06) is a mammalian MoClo kit but delivers by PiggyBac transposon, so it ships no
lentiviral transfer vector and the question does not arise.

**The pattern, stated as the conclusion it is.** Every published mammalian Golden Gate toolkit
that ships a lentiviral vector puts that vector at a **BsmBI** step. None puts it at a BsaI
step. That is not a coincidence anyone wrote down — no source read here says "we chose BsmBI
because the LTR has a BsaI site" — but the sequence in section 1 explains it, and the absence of
any BsaI-based lentiviral toolkit is consistent with it. **Marked as inference.** What would
close it: a methods section or a correspondence from one of these groups saying why.

### The one published precedent: VAMSyB, read first-hand

**Haellman, V., Strittmatter, T., Bertschi, A., Stücheli, P. & Fussenegger, M. A versatile
plasmid architecture for mammalian synthetic biology (VAMSyB). *Metab. Eng.* **66**, 41–50
(2021). PMID 33857582,
[doi:10.1016/j.ymben.2021.04.003](https://doi.org/10.1016/j.ymben.2021.04.003).**

Issue #244 recorded this as an unverified search summary. **It is true, and it is now read
first-hand.** The published CC-BY PDF was retrieved on **2026-10-06** from the ETH Research
Collection, handle [20.500.11850/480636](https://doi.org/10.3929/ethz-b-000480636), by way of
the DSpace REST API — ScienceDirect still answers a direct fetch and an `r.jina.ai` fetch with
a captcha page.

The paper says it in one sentence. Their Tier-3 lentiviral vectors "were adapted from a
publicly available third-generation lentiviral transfer vector (Sanjana et al., 2014)" — that
is lentiCRISPR v2 — and:

> Furthermore, these vectors were mutated to remove BsaI and BsmBI cut sites within the LTRs'
> and WPRE enhancer. However, HindIII cut sites are too abundant to be easily removed, and we
> refrained from testing all mutations for functionality.

#### What actually moved, base by base

A scan run here on **2026-10-06** over the Addgene-verified full sequences of the child,
pTS1106_Tier3(Lenti) ([Addgene 169661](https://www.addgene.org/169661/), 5,259 bp), and of its
parent, lentiCRISPR v2 ([Addgene 52961](https://www.addgene.org/52961/), 14,873 bp), both
compared against GenBank K03455.1. **Measured:**

| | lentiCRISPR v2 | pTS1106_Tier3(Lenti) |
| --- | --- | --- |
| R+U5, both copies | identical to HXB2 455–634 at every base | two substitutions, the **same two in both copies** |
| R+1..+6 | `GGTCTC` — BsaI | `GGTTTC` — **C→T at R+4** |
| R+77..+82 | `AAGCTT` — HindIII | `AAGTTT` — **C→T at R+80** |
| BsaI sites, whole plasmid | 4 | 2, both in the designed A2 acceptor site |
| BsmBI sites, whole plasmid | 2 | 2, both in the designed A3 acceptor site |

So the LTR edit is **one base per site, in both copies**: `GGTCTC` → `GGTTTC` at R+4, which is
exactly the shape section 3 arrives at independently — a single change at nt 4, 5 or 6, nt 1
left alone, made in both LTRs.

Two corrections to the paper's own sentence, and both matter:

- **BsmBI was never in the LTRs** to remove, as section 1 measures. The two BsmBI sites in the
  parent are the sgRNA cloning sites, not viral sequence.
- **The WPRE did not change.** Its 589 annotated bases are byte-identical in parent and child,
  and carry no BsaI, BsmBI or HindIII site in either. Nothing was removed there.
- The **HindIII** site inside R *was* removed, in both copies, by the second C→T — which the
  paper does not say, and which sits two bases 3' of the `AATAAA` poly(A) hexamer at R+73..+78.
  The hexamer itself is untouched.

#### And nobody measured it

**This is the half that does not close.** The words "titre", "titer" and "MOI" do not appear in
the paper. The lentiviral method says only that "HEK-293T cells were transduced by adding
various amounts of virus-containing supernatant to the culture medium", then describes a
puromycin ramp to 2 µg/mL over 14 days. The whole lentiviral result is Fig. 3E–G: a polyclonal
line that answers doxycycline with SEAP and with iRFP670, N = 3 in quadruplicate.

**There is no titre, no transduction efficiency, and no side-by-side against the undomesticated
parent anywhere in the paper.** Puromycin selection hides exactly the quantity at issue: a
vector at 1% of parental titre and one at 100% both give a resistant polyclonal pool.

#### What this is worth

**Measured, and it is a real precedent:** a third-generation transfer vector carrying C→T at
R+4 in both LTRs packages into VSV-G particles with a third-generation packaging system,
transduces HEK-293T, integrates, and drives a three-cassette Tet-ON circuit that responds to
dose. The change is not lethal. That is more than this note had.

**Not measured, and the gap is unchanged:** what it costs. Section 3's gap 2 — functional titre
of any lentiviral vector with a modified R — survives VAMSyB intact. Anyone copying this edit
copies a vector that works, not a vector whose yield is known.

**Not established.** "LentiMoClo" does not exist as a published toolkit. Mobius/Loop mammalian,
MTK2, GoldenBraid mammalian and CIDAR mammalian lentiviral variants were not settled either way
— they ran past the time budget. **No Tet-on (TRE / TetOne / Tet3G) lentiviral backbone already
clean for BsaI was found.** Clean for BsmBI is expected to be common, since the LTRs contribute
none, but no specific TRE vector sequence was scanned. What would close it: scan the Addgene
sequences of the TRE lentiviral vectors directly. **Correction, 2026-10-06:** this does not need
a login. The bases and the enzyme list behind an Addgene map are served to anyone, and that is
how section 2's base-by-base comparison was run.

## 3. TAR, and the six bases

All sources in this section read **2026-10-06**.

### The lower stem is structure-critical, not sequence-critical

Three measurements say the same thing, and they are the load-bearing result of the whole note.

| Source | Measured | Result |
| --- | --- | --- |
| Klaver, B. & Berkhout, B. *EMBO J.* **13**, 2650–2659 (1994). PMID 8013464, PMC395139, [doi:10.1002/j.1460-2075.1994.tb06555.x](https://doi.org/10.1002/j.1460-2075.1994.tb06555.x) | Base substitutions in the TAR **lower stem** of replication-competent HIV-1 | Viruses dead. Revertants recovered only after >6 months; the point mutations and short deletions that rescued them "supported the notion that **base-pairing** of the lower stem region is essential" |
| Das, A. T., Klaver, B. & Berkhout, B. *J. Virol.* **72**, 9217–9223 (1998). PMID 9765469, PMC110341, [doi:10.1128/jvi.72.11.9217-9223.1998](https://doi.org/10.1128/jvi.72.11.9217-9223.1998) | Mutant "Xho+10" replaces nt **+3 to +16** with unrelated sequence — four of the six bases at issue, plus ten more | 5'-only: RNA ~50% of wild type, packaging 70%, CA-p24 similar. 3'-only: RNA 100%, packaging 70%. Both: RNA ~50%, packaging 40%. **A revertant that restored lower-stem pairing restored packaging** |
| Das, A. T., Vrolijk, M. M., Harwig, A. & Berkhout, B. *Retrovirology* **9**, 59 (2012). PMID 22828074, PMC3432602, [doi:10.1186/1742-4690-9-59](https://doi.org/10.1186/1742-4690-9-59) | Deletions on **one** side of the TAR stem versus **both** sides | "opening the 5' TAR structure through a deletion in **either side** of the stem region caused aberrant dimerization and reduced packaging… In contrast, **truncation of the TAR hairpin through deletions in both sides of the stem did not affect RNA dimer formation and packaging**" |

The 2012 result is the clean control: break the duplex and the vector suffers; shorten the
duplex while keeping it paired and it does not. **The identity of the lower-stem bases is
dispensable. The duplex is not.**

Two further constraints on any redesign:

- **A destabilised TAR damages its neighbour.** Vrolijk, M. M. *et al.* *Retrovirology* **6**, 13
  (2009). PMID 19210761, PMC2645353,
  [doi:10.1186/1742-4690-6-13](https://doi.org/10.1186/1742-4690-6-13). Measured: destabilising
  TAR extends the adjacent polyA hairpin and **inhibits polyadenylation**. Unpaired bases do not
  stay local.
- **Over-stabilising is also costly.** Da Silva Amaral, C. *et al.* *J. Virol.* (2026). PMID
  42053302, PMC13185573, [doi:10.1128/jvi.01840-25](https://doi.org/10.1128/jvi.01840-25).
  Measured: mutations that **stabilise** the lower stem of cTAR impair NC-mediated cTAR–TAR
  annealing, reduce replication, and reduce full-length ssDNA in infected cells. A redesign must
  land near wild-type free energy, not merely "paired".

### Is TAR needed at all once Tat is gone?

**The one published death of lower-stem mutants was measured with Tat driving the HIV-1 LTR** —
Klaver & Berkhout above infer the defect is "at the level of transcription from an integrated
provirus". That mechanism does not exist in a third-generation vector whose producer transcript
comes from CMV or TRE and whose cells contain no Tat. Das 2012 reports that "**complete TAR
deletion is allowed** in the context of an HIV-1 variant that does not depend on this Tat-TAR
axis for transcription".

So the transcriptional job of TAR is gone. **Three jobs survive into a Tat-independent vector,
and all three are structural:**

- **First strand transfer.** Berkhout, B., Vastenhouw, N. L., Klasens, B. I. F. & Huthoff, H.
  *RNA* **7**, 1097–1114 (2001). PMID 11497429, PMC1370158,
  [doi:10.1017/s1355838201002035](https://doi.org/10.1017/s1355838201002035). In an *in vitro*
  donor/acceptor assay, mutating the TAR **loop** of the donor RNA reduced transfer **more than
  fivefold**, and the authors conclude "it is not the RNA hairpin motif in the 5'R donor, but
  rather the **antisense motif in the ssDNA copy**… that is critical for strand transfer".
- **Dimerisation and packaging.** Das 1998 and Das 2012, numbers above.
- **Polyadenylation**, through the neighbouring hairpin. Vrolijk 2009.

### Do both LTR copies have to change?

**Which copy is inherited.** Das 1998 states that "mutations introduced in 5' TAR will be
**inherited in a dominant manner in both R regions**", and that "**3' TAR mutations will be lost
after a single replication cycle**". This is the authors' interpretation of their own experiment;
the full text was blocked this run and **no proviral sequencing panel was seen**. Treat it as
strongly supported mechanism, not as a sequencing measurement. What would close it: the PDF of
PMC110341, read for the sequencing figure.

**An asymmetric change survives — measured.** Das 1998's 5'-only and 3'-only mutants both
replicated, and on packaging the asymmetric constructs (70%) beat the symmetric double (40%).
That is a measurement against the intuition that the copies must match.

**But the modern field builds the symmetric mutant on purpose.** Da Silva Amaral 2026 state that
they put the mutations in both TARs so the constructs "**cannot form mismatches between the cTAR
sequence of the ssDNA and the TAR sequence at the 3'-end of the genomic RNA**", and that "the
lower stem of the 3' TAR hairpin was preserved in the double mutants because its opening
indirectly affects the **polyadenylation** of the viral transcripts".

That is the direct confirmation of the ticket's framing: **the two copies do different jobs in
the producer cell.** The 5' copy is the transcript's 5' end and the cTAR template; the 3' copy
is the polyadenylation context and the strand-transfer acceptor. A 5'-only change makes a
cTAR/3'-TAR mismatch at first strand transfer; a 3'-only change perturbs polyadenylation and is
then lost anyway.

### Verdict on positions 1–6

**Likely survivable, conditional on three things, and no one has measured it.**

Survivable because the lower stem is load-bearing by pairing and not by sequence (three
measurements), and because the one lethal result was a Tat-dependent transcription defect that a
third-generation vector does not have. Conditional on: preserving the duplex with a compensating
change in the partner base; keeping free energy near wild type, checked against the adjacent
polyA hairpin as well; and making the change in **both** LTRs.

**Confidence: moderate-to-high on the mechanism, low on the magnitude.** Everything above is
HIV-1 virus biology plus an *in vitro* strand-transfer assay, carried into a vector by inference.

One position deserves separate treatment. **Nucleotide 1 is +1**, the transcription start site
(section 1 fixes it from the TATA box at offset 427). Changing +1 risks initiation efficiency and
cap status, which is a different failure from a broken hairpin, and **nothing measures it**
(gap 1 below). Killing `GGTCTC` needs only one base; **inference, marked:** prefer a single
transversion at nt 4, 5 or 6 over a rewrite of all six, and leave nt 1 alone.

**One group has done this, and section 2 reads their sequence.** VAMSyB's Tier-3 lentiviral
vector carries C→T at R+4 in both LTRs, and the vector packages and transduces. It is a
transition, not a transversion, and no second change was made inside TAR — so it does **not**
satisfy the compensation condition above, and it still worked. **What it does not supply is
a number:** they reported no titre, so the precedent says the change is survivable and says
nothing about the cost. Gap 2 below is unchanged by it.

The partner base must be read off a fold of the **real** R sequence of the backbone in hand. A
complementarity scan run here on 2026-10-06 over HXB2 TAR +1..+59 finds the best 6-nt partner
for `GGUCUC` at nt 52–57 (`GGAACC`), pairing 5 of 6 with wobble allowed — close to, but not the
same as, the textbook register. **Do not assume nt 1–14 pairs with nt 44–57; fold the sequence.**

### What nothing measures, in this section

1. **Pol II initiation efficiency, +1 position or capping as a function of R nt 1–6 identity
   under a heterologous promoter.** Every transcription measurement found used the HIV-1 LTR,
   where a +1 change and a Tat-response change are confounded. *Closes it:* make the change in
   the real backbone and measure producer-cell vector RNA by RT-qPCR against a plasmid-copy
   control, plus the 5' end by 5'-RACE, plus functional titre.
2. **Functional titre of any lentiviral transfer vector with a modified or deleted TAR/R.**
   Patents claim hybrid R regions; no measured titre was retrieved. *Closes it:* a side-by-side
   titre against the parent — one experiment, and it also closes 1.
3. **A designed compensatory double mutation in the TAR lower stem of a vector**, as opposed to a
   selected revertant in virus (Das 1998) or a both-sides deletion (Das 2012).
4. **How many mismatched bases between the 5' and 3' R copies are tolerated at first strand
   transfer, and at what cost.** *Closes it:* the Berkhout 2001 *in vitro* assay, run with a
   changed donor against a wild-type acceptor.
5. **Whether the one-base asymmetry pLVX-TetOne already carries (section 1) behaves differently
   from a six-base one.** It is suggestive, not a measurement of tolerance.

## 4. The tolerance map

All sources read **2026-10-06**. A number marked **(unverified)** was taken from an automated
extraction of the journal page and was not matched against the saved text; it wants a one-line
check before anyone designs against it. Everything else was matched to a saved dump.

| Element | What it does | Class | The measurement |
| --- | --- | --- | --- |
| **U3, 5' LTR** | In a third-generation transfer plasmid it drives nothing — a heterologous enhancer replaces it | **Permissive; wholly replaceable** | pRRL fuses an RSV enhancer to HIV-1 R, pCCL a CMV one. pRRL gave 1.0 × 10⁷ TU/mL and 13,495 TU/ng p24, against 2nd-gen pHR2 at 4.1 × 10⁶ and 13,805 — higher particle titre, identical infectivity per p24 (unverified). Dull 1998 |
| **U3, 3' LTR (the ΔU3 / SIN deletion)** | Carries TATA, Sp1, NF-κB; the residual stubs carry the integrase *att* site and the element licensing 3' polyadenylation | **Mostly permissive; two sequence-critical stubs** | SIN-18 deletes −418 to −18, leaving "only 53 nucleotides… 35 nucleotides upstream… to preserve efficient recognition and processing by integrase and 18 nucleotides downstream… to govern polyadenylation". Titre essentially unchanged (1,476 ± 232 vs 1,544 ± 126 TU/ng p24, unverified); LTR transcription down **350-fold**. Zufferey 1998. Independently: a 133 bp deletion, −141 to −9, titre 5.0 × 10⁵ vs 5.7 × 10⁵ TU/mL (unverified). Miyoshi 1998 |
| **R, outside TAR — the poly(A) signal** | `AATAAA` plus cleavage at the R/U5 boundary; present at both ends, so the 5' copy must be suppressed and the 3' copy used | **Sequence-critical** (hexamer and an upstream element), **plus structure-critical** (the polyA hairpin) | Processing "require[s] sequences 76 nucleotides upstream of the AAUAAA hexamer", which bind CPSF-160 directly. Gilmartin 1995. Section 3's Vrolijk 2009 supplies the hairpin half |
| **U5** | Holds the GU/U-rich downstream element of the poly(A) site, the 3' half of the polyA hairpin, and the tRNA-primer surface | **Mixed — and the structural half is an inference** | The downstream element is sequence-defined (Gilmartin 1995). The measured U5–tRNA structural requirement is in **RSV/tRNA-Trp, not HIV-1** (Aiyar 1992). **Calling HIV-1 U5 structure-critical is an inference by analogy** and is marked as one |
| **PBS** (18 nt, tRNA-Lys3 complement) | Primes minus-strand synthesis | **Sequence-critical in effect, though not absolutely — and it reverts** | PBS swapped to five other tRNA species: viruses "were able to replicate, although with delayed replication kinetics… all mutants reverted to the wild-type PBS(3Lys) sequence". Das 1995. Even a single-base change reverts fast — 10% by 12 h, 50% within 8 days. Berkhout & Das 2015. **Reversion is the point: a PBS edit does not stay edited** |
| **Psi (SL1–SL4)** | Gag recognition and genome encapsidation | **Structure-critical**, with a sequence-critical subset of base pairs | The compensation experiment: destabilising substitutions plus "second-site compensatory mutations that would restore secondary structure" identified two hairpins that work by fold; "a subset of the hydrogen-bonded base pairs within the stems… is likely to be required for function in cis". McBride & Panganiban 1996. Three stem-loops each bind Gag independently at ~4-fold lower affinity than full psi. Clever 1995 |
| **Major splice donor** | `GT` plus consensus — and it sits **inside the SL2 stem**, so an SD point mutation is also a psi structure mutation | **Dual, and base changes are tolerated at modest cost** | "Mutations of SD from GGTG to GCAG (SD1) or to GGGG (SD2) resulted in a moderate decrease (10 to 30%) in vector titer", while abrogating splicing; cytoplasmic full-length RNA fell >70% (unverified). Cui 1999 |
| **RRE** (351 nt) | Rev-binding RNA; exports the unspliced genome | **Stem-loop IIB sequence-critical; the rest structure-critical and conformationally plastic** | Stem-loop II is "both necessary and sufficient for the binding of Rev". Malim 1990. And the fold itself carries output: wild-type RRE "exists as a mixture of both structures", and the five-stem form "promotes greater functional Rev/RRE activity". Sherpa 2015. **This is the strongest warning against silent changes outside stem IIB** — a change that shifts the conformer balance costs output without touching the Rev site. Deleting RRE cost >90% of titre (unverified). Cui 1999 |
| **cPPT/CTS** | Central plus-strand initiation and termination, making the central DNA flap | **Sequence-critical in principle; its benefit is contested** | For: "a 10-fold increase in the amount of integrated DNA… in the presence of the cPPT". Van Maele 2003; Zennou 2000, 2001; Follenzi 2000. **Against, and a tolerance map must carry it:** ten base changes within the cPPT, in two backbones — "the wild-type and cPPT-D mutant viruses replicated with similarly robust kinetics" (unverified). Dvorin 2002 |
| **3' PPT** (15 nt) | RNase-H-resistant plus-strand primer | **3' end sequence-critical; 5' end permissive** | "Changes at the 5' end of the PPT have no effect on primer function, whereas the identity of bases at the 3' end is crucial"; a primer keeping only the six 3' G residues worked like wild type. Powell & Levin 1996. In virus, a fully randomised PPT retained **12%** of wild-type infectivity and reverted quickly. Miles 2005 |
| **WPRE** | Post-transcriptional enhancer of transgene expression | **Permissive — and it carries the one published domestication precedent** | Expression up ≥5-fold, and "the WPRE did not influence titers" (unverified). Zufferey 1999. No single sub-element is essential: α, β and γ each give ~12% of full activity alone (unverified). Donello 1998. **WPRE-mut6**: mutating the woodchuck X-protein ORF start codon out — "the native form or mutated derivatives functioned equivalently". Zanta-Boussif 2009. The abstract gives no fold number; the panel is paywalled |
| **Δgag fragment** | Residual gag 5' end carrying SL4 and the gag start | **Permissive at ~30% titre cost** | A 5' gag deletion cost 30% of titre with no significant RNA effect; a combined gag/env/RRE/SD1 deletion still reached ~50% of wild type (unverified). Cui 1999 |
| **Internal promoter and env/SA region** | The heterologous cassette | **Permissive — whole-cassette swaps are free** | PGK, CMV and gene-specific internal promoters all give normal titres across Zufferey 1998, Dull 1998 and Zanta-Boussif 2009. **Nothing measures silent base changes *inside* an internal promoter** — the evidence is for swapping the whole cassette |

### Three things in this table that change how you domesticate

- **Compensatory changes are a measured technique here, not a hope.** McBride & Panganiban
  built destabilising mutations and second-site restorations in psi and read the fold back out.
  That is the same move section 3 recommends for TAR.
- **Two elements revert rather than fail.** The PBS and the 3' PPT both restore wild-type
  sequence under selection. An edit there does not produce a stable broken vector; it produces a
  vector that quietly turns back into the parent, which a sequencing check at passage 1 would
  miss.
- **The RRE punishes silent changes.** It is the one element where a change that touches no
  recognised motif can still cost output, by moving a conformer equilibrium. Treat the whole 351
  nt as off-limits, not just stem IIB.

### What nothing measures, in this section

| Open | What would close it |
| --- | --- |
| Tolerance to base change **inside** the retained 35 nt integrase stub and 18 nt poly(A) stub of a SIN 3' LTR. Zufferey 1998 kept them deliberately but ran no substitution series in them | A scanning triple-base substitution walk across both stubs in a SIN-18 backbone, read as TU/mL plus integrated copies by qPCR and 2-LTR circles |
| HIV-1 **U5** base changes in a vector. The only measured tRNA–U5 structural requirement is in RSV | Compensatory mutation pairs across the HIV-1 U5–PBS helix in a SIN vector, titred |
| Per-stem-loop packaging and titre numbers for **SL1–SL4**. The two psi papers are abstract-and-references at PMC | The PDFs of PMC190155 and PMC188876 through institutional access |
| Whether the splice donor's 10–30% titre cost is **splicing or SL2 structure** — the two are confounded because the SD sits in the stem | Second-site changes restoring SL2 pairing while keeping the SD dead. Titre recovery would settle it |
| **RRE length versus titre** in an HIV-1 vector. Only the BIV analogue was found | A nested truncation series, 351 → 234 → stem-IIB-only, in one backbone, titred, with the spliced/unspliced cytoplasmic ratio |
| Base changes in the sequence **flanking** cPPT/CTS, independently of the cPPT itself | Vary the flanks alone and measure integrated copies |
| Silent base changes **inside** an internal promoter — which is where two of pLVX's own sites sit | A reporter assay for promoter strength before and after |

## 5. Alternative and engineered LTR designs

| Design | Changed for | Measured |
| --- | --- | --- |
| **RSV-U3/HIV-R (pRRL)** and **CMV-U3/HIV-R (pCCL)** chimeric 5' LTR | Tat independence, so the third-generation system can drop *tat* | pRRL 1.0 × 10⁷ TU/mL; without Tat it kept 60% of control against 5–15% for pHR2 (unverified). Dull 1998 |
| **CMV/LTR hybrid 5' LTR plus ΔU3 3' LTR** | The same, independently | Hybrid 5' LTR *raised* titre: 9.3 × 10⁵ and 8.2 × 10⁵ against 5.7 × 10⁵ TU/mL for wild type (unverified). Miyoshi 1998 |
| **Varying ΔU3 extent (SIN-18 / -36 / -45 / -78)** | Trading residual LTR activity against titre | SIN-78, which keeps TATA and three Sp1 sites, still showed LTR transcript and a 21–30-fold GFP rise on HIV superinfection; SIN-18/-36/-45 showed none. **The deletion must reach past the TATA box, not merely past the enhancer.** Zufferey 1998 |
| **Insulated LTR (cHS4 in the ΔU3)** | Shielding the transgene from position effects | Titre falls with insert length **in the 3' LTR specifically** — the same 1.2 kb placed internally cost nothing. The block is post-entry: less reverse transcription and integration, with equal full-length mRNA in packaging cells. Urbinati 2009 |
| **Non-primate backbones (EIAV, FIV, BIV, SIV)** | Avoiding HIV-derived sequence | Only BIV was verified: its RRE maps to 312 bp of *env* and the first 104 bp of *gag* carries part of the packaging signal; a minimal BIV SIN vector was built from these. Molina 2002. **EIAV and FIV titre comparisons surfaced only in vendor PDFs and reviews and are not cited here** |
| "SIN with internal poly(A)", sequence-reduced LTRs, integrase-deficient designs | — | **Nothing primary found. Open** |

Two leads were found and **not** read, both likely to matter:

- **Wright *et al.* (2026)**, "supA-LTRs" — sequence-upgraded polyA LTRs in real lentiviral
  vectors, PMC13370168, open access. This is the closest published thing to changing R-region
  sequence in a vector and reporting titre, which is gap 2 of section 3.
- **Wright *et al.*, *Heliyon*** — major-splice-donor mutation plus a U1 snRNA enhancer,
  reporting 5.7- to 8.5-fold titre improvement. The numbers reached this note only through a
  search snippet and a vendor-hosted PDF, so they are **recorded as a lead, not cited**.

## 6. If a permissive change is wrong, how would you know?

The map above calls several things permissive. This section says how each would fail and what
catches it. For a pooled library the ordering principle is not cost but **loudness**: a change
that halves titre announces itself in no readout the bench normally runs.

All sources read **2026-10-06**.

### The two results that should change how you think about this

**A restriction site inserted near psi cost 82% of titre, and the same insertion three bases
further out cost nothing.** Kim, S. H., Jun, H. J., Jang, S. I. & You, J. C. *PLoS One* **7**,
e50148 (2012). PMID 23185560, PMC3503997,
[doi:10.1371/journal.pone.0050148](https://doi.org/10.1371/journal.pone.0050148). They put
restriction sites around psi in a lentiviral transfer vector — the closest published thing to a
domestication edit — and titred:

| Vector | p24 (ng/mL) | Titre (CFU/mL) | Transduction, % of parent |
| --- | --- | --- | --- |
| parent | 298 ± 12 | 235,000 | 100 |
| ES2.psi | 252 ± 3 | 41,333 | 18 |
| ES3.psi | 293 ± 8 | 215,000 | 91 |
| PBS.MAp | 282 ± 13 | 60 | 0.03 |

**p24 varies 1.2-fold across the panel while functional titre spans 3,900-fold.** Two lessons:
position matters at a resolution of a few bases, and **p24 alone would have told them nothing**.

**A two-base change can cost what neither base costs alone.** Sakuragi, S. *et al.*
*Int. J. Mol. Sci.* **22**, 3435 (2021). PMID 33810482,
[doi:10.3390/ijms22073435](https://doi.org/10.3390/ijms22073435). A two-base substitution
downstream of the PBS inside psi, in no reading frame, measurably lowered packaging — while
**either single base alone changed nothing**. If an edit lands in the 5' UTR or psi, assay it
rather than reasoning about it.

### The panel

| # | Failure | How it shows | The assay | Source |
| --- | --- | --- | --- | --- |
| 1 | **Lost titre** | In a library, as nothing: you transduce the planned cells, get fewer transductants, and believe you are at 500× coverage when you are at 100× | Functional titre by flow, and **ddPCR**. The universal assay targets **RRE** duplexed with GAPDH, so one assay serves a whole library; LLOQ 1.16 copies/µL, 3.6% CV against the NIST one-copy standard, agreeing with flow within 30% | Kandell, J. *et al.* *Mol. Ther. Methods Clin. Dev.* **31**, 101120 (2023). PMID 37841416, [doi:10.1016/j.omtm.2023.101120](https://doi.org/10.1016/j.omtm.2023.101120) |
| 1b | **Assay precision** | — | ddPCR CV **4–8%** against qPCR **10–15%** on cloned standards at VCN 1, 2, 3. **Use a provirus-specific design**: a psi-targeting assay also amplifies carried-over transfer plasmid and inflates titre | Corre, G. *et al.* *Gene Ther.* (2022). PMID 35194185, [doi:10.1038/s41434-022-00315-8](https://doi.org/10.1038/s41434-022-00315-8) |
| 1c | **Where it broke** | Normal particle yield, few transductions | **The p24-to-TU ratio**, as Kim 2012 above and Dull 1998 (section 4) both use it | — |
| 2 | **Truncated genome, cryptic poly(A), cryptic splicing** | A fraction of packaged RNA is not full length. Silent in every bulk assay | **Nanopore direct RNA sequencing of virion RNA** — each read's 3' terminus is the RNA's 3' end, so truncation sites fall out directly. Measured on a clinical vector: a cryptic poly(A) inside **WPRE** used by 4.44% of reads (up to 10% across vectors), two hairpin truncation sites totalling 19.88%, and **full-length RNA only ~60–75% even in a working vector**. Only ~1.5–3% of reads are on-target, so run deep | Pal, A. *et al.* *Genome Res.* **34**, 1966 (2024). PMID 39467647, [doi:10.1101/gr.279405.124](https://doi.org/10.1101/gr.279405.124) |
| 2b | **Read-through past the 3' LTR** | Transcription into flanking DNA | **A warning for domestication specifically: termination and promoter are the same bases.** 70–80% of termination activity sits in a 124-nt region overlapping the NF-κB, Sp1 and TATA sites, and deleting NF-κB, Sp1 or TATA **each raised read-through**. The HIV-1 SIN poly(A) is already as leaky as MLV's | Yang, Q. *et al.* *Retrovirology* **4**, 4 (2007). PMID 17241475, [doi:10.1186/1742-4690-4-4](https://doi.org/10.1186/1742-4690-4-4); Zaiss, A. K. *et al.* *J. Virol.* **76**, 7209 (2002). PMID 12072520 |
| 3 | **Silent provirus** | The member is in the pool, amplifies normally in a barcode readout, and contributes no phenotype. It reads as "this variant is inactive", not "this variant is broken" — **the worst failure for a screen** | **VCN and percent expressing, on the same cells, over time. The gap between them is the signature.** Measured: VCN flat at 1.1 ± 0.004 → 0.9 ± 0.3 while eGFP+ marrow fell from 5.4–17.6% at 4 weeks to 0.0–0.6% at 36 weeks, with 18 of 18 promoter CpGs over 90% methylated. **Promoter choice dominates**: CMV fell 3.6–5.2×, PGK and EF1α only 1.5–1.7× | Herbst, F. *et al.* *Mol. Ther.* **20**, 1014 (2012). PMID 22434137, [doi:10.1038/mt.2012.46](https://doi.org/10.1038/mt.2012.46) |
| 3b | **Baseline to compare against** | — | ~1.5% of integrations (1 in 66) land silent even with a good vector, and much of that is integration site rather than vector sequence. HDAC-inhibitor rescue fires on only some clones, so **a negative result does not prove the provirus is absent** | Jordan, A. *et al.* *EMBO J.* **22**, 1868 (2003), [doi:10.1093/emboj/cdg188](https://doi.org/10.1093/emboj/cdg188); Contreras, X. *et al.* *J. Biol. Chem.* **284**, 6782 (2009), [doi:10.1074/jbc.M807898200](https://doi.org/10.1074/jbc.M807898200) |
| 4 | **Lost induction or new leak** (Tet-on) | Lower induced signal, or expression without doxycycline | Fold induction and leak as **separate** numbers, judged **at single integrated copy** — regulation collapses on integration: one promoter went from ~50,000-fold transient to ~7,900-fold single-copy. **Flow cannot measure the leak**: the off state sat too close to a 1.83 mfu background for any conclusion, and luciferase on the same cells resolved it | Loew, R. *et al.* *BMC Biotechnol.* **10**, 81 (2010). PMID 21106052, [doi:10.1186/1472-6750-10-81](https://doi.org/10.1186/1472-6750-10-81) |
| 5 | **Skewed representation** | The screen runs and the answer is wrong | NGS of the barcode amplicon at **plasmid pool, then after transduction, then endpoint** — use the plasmid pool as reference, not the first timepoint. Report the 90th/10th percentile ratio; a published acceptance bar is **under 10**, with over 70% perfect-match reads. A low-skew library held its skew flat from 200× down to 50× coverage | Heo, S. J. *et al.* *Genome Biol.* **25**, 19 (2024), [doi:10.1186/s13059-023-03132-3](https://doi.org/10.1186/s13059-023-03132-3); Joung, J. *et al.* *Nat. Protoc.* (2017) |
| 6 | **Packaging or dimerisation defect** | Normal transcription, less genome per particle | **Competitive RT-qPCR packaging ratio**: virion RNA over cytoplasmic RNA, against an internal control differing only by silent changes, so transfection and recovery variation cancel. A **native-gel Northern** separates a dimerisation defect from a packaging one — SL1 loop mutants lose infectivity at **normal virion RNA content**, which the packaging ratio alone would miss | Sakuragi 2021, above; Clever, J. L. & Parslow, T. G. *J. Virol.* **71**, 3407 (1997). PMID 9094610 |
| 7 | **Recombination between repeats** | A deletion between the LTRs, or a barcode uncoupled from the part it names | Whole-plasmid nanopore sequencing. At the reverse-transcription level the rate is high — roughly 5.5 to 9 crossovers per genome per cycle in T cells and primary CD4+. **But it depends on co-packaging two different genomes, so it scales with MOI and pooling, not with your edits.** The mitigation is low MOI, not sequence design | Levy, D. N. *et al.* *PNAS* **101**, 4204 (2004), [doi:10.1073/pnas.0306764101](https://doi.org/10.1073/pnas.0306764101); Rhode, B. W. *et al.* *J. Virol.* **61**, 925 (1987) |

### Does a titre drop actually skew a library?

Less than instinct suggests, and the arithmetic is worth carrying. Titre enters a pooled library
only through coverage. Counts are Poisson, so at 500× the sampling contribution to the 90/10 skew
ratio is about 1.12; at 250× — a **2-fold** titre loss — about 1.18; at 50× about 1.44. Against a
cloning-derived skew of 2 to 10, a 2-fold titre loss is **invisible in the skew metric**. The
Poisson floor starts to dominate only below roughly 10–25×, which is a 20- to 50-fold titre loss.
This is a calculation, **not a measurement**, but it agrees with Heo's coverage titration and with
the finding that cell-splitting coverage, not transduction, is the dominant bias term
(Imkeller, K. *et al.* *Genome Biol.* **21**, 53 (2020),
[doi:10.1186/s13059-020-1939-1](https://doi.org/10.1186/s13059-020-1939-1)).

**So the risk from a quiet titre loss is operational, not statistical.** You do not see it in the
skew; you see it as believing you were at 500× when you were at 100×, and then the split-and-grow
bottleneck acts on a smaller population than planned. **Measure the titre; the skew will not tell
you.**

### The panel, cheapest first

1. **Whole-plasmid sequencing of the new prep and the parent.** Catches every unintended edit and
   any rearrangement between the LTRs, before anything is spent on virus. It is also the only
   step that catches a **reverting** element (section 4) — sequence the prep, not the design.
2. **One paired packaging run: parent and domesticated, same day, same mix, n ≥ 3. p24 and
   functional titre from the same harvest.** Two numbers, and their ratio localises any loss to
   before or after entry.
3. **ddPCR vector copy number on transduced cells**, provirus-specific. Resolves a 2-fold change
   at n = 2.
4. **Flow for expression; for a Tet vector, induction and leak — and put a luminescent reporter
   on it if the leak matters.**
5. **Percent expressing divided by copy number, at four or more timepoints over six weeks**, for
   silencing.
6. **For a library, always: barcode NGS at the plasmid stage and after transduction.** Plus one
   cheap addition that closes this lab's specific blind spot — **sort expressing from
   non-expressing cells out of the same pool and sequence barcodes from both.** A member enriched
   in the non-expressing fraction at unchanged total abundance is silenced, not lost.
7. **Conditional on where the edits landed:** the packaging ratio if anything sits in the 5' UTR
   or psi; direct RNA sequencing if anything sits in WPRE, psi or the U3 control region;
   read-through if anything sits in U3.

**The one design rule that matters more than the panel:** keep every comparison **paired and
same-day against the undomesticated parent**. Every effect in this literature is a ratio to a
parent, and batch variance in lentiviral packaging is larger than any of these assays' CV.

### The nearest thing to a worked precedent

**No paper was found that removes a restriction site from a lentiviral transfer vector and
reports what it checked afterwards.** The nearest four, in descending usefulness: Kim 2012 and
Sakuragi 2021 above; Cui 1999 (section 4), which shows titre and RNA level are decoupled — a dead
splice donor cost 10–30% of titre despite losing over 70% of cytoplasmic full-length RNA; and
**Koldej, R. M. & Anson, D. S.** *BMC Biotechnol.* **9**, 86 (2009). PMID 19811661,
[doi:10.1186/1472-6750-9-86](https://doi.org/10.1186/1472-6750-9-86) — a whole-backbone recoding
that reduced LTR-to-LTR homology, then measured titre at n = 6, 3'-LTR read-through by TaqMan
across the vector-to-flank junction, and SIN-repair rate by colony PCR. **That is the QC panel to
copy.** Its sting: SIN repair still happened with *zero* U3–U3 homology, so reducing homology is
not a cure.

### What nothing measures here

- **How much titre loss a pooled library tolerates before representation skews measurably.** The
  table above is a calculation, not an experiment. *Closes it:* titrate a barcoded pool at
  several coverages and measure recovery against input.
- **The 5'/3' RT-qPCR ratio on packaged RNA.** No primary source defines it, calibrates it
  against a known-truncated control, or states a detection limit. *Closes it:* spike a
  3'-truncated transcript into full-length RNA at 0, 5, 10, 25 and 50% and find the smallest
  distinguishable fraction.
- **Per-base sensitivity of 3'-LTR read-through.** Yang 2007 measured deletions, not
  substitutions — which is exactly the resolution a domestication edit works at.
- **The smallest fold change flow can resolve.** No source states one; it is instrument- and
  fluorophore-specific. *Closes it:* a 2-fold ladder on the lab's own cytometer.
- **How a silenced member behaves in a pooled barcode readout.** Nobody measured it. The sort
  experiment in step 6 above would.
- **Plasmid-level recombination frequency for modern SIN transfer plasmids.** The only primary
  measurement is from 1987, MLV-based, at kilobase scale, and states no number.

## 7. Using this map on the next backbone

A procedure, in the order that spends least.

1. **Count the sites for every candidate enzyme, not just the one you meant to use.** Section 1's
   table is for HIV-1; run the same count on the backbone in hand with
   `liulab_mbio.sites.find_sites`. If one enzyme is already clean, the whole question closes.
2. **Separate the sites you must clear from the sites you merely found.** Only the enzyme an
   assembly reaction actually sees matters. #230 made this point on pLVX: narrowing the
   requirement from "no BsaI and no BsmBI" to "no BsmBI" turned two research questions into two
   synonymous codon changes.
3. **Sort what is left by this map.** Coding sequence first — a synonymous change is free and
   `liulab_mbio.sites.domesticate` makes it. Then unannotated spacer. Then the elements this note
   marks permissive. Leave sequence-critical elements alone; for structure-critical ones, design
   a compensating change and check the fold.
4. **Count the copies before you count the edits.** A site in R or U5 exists in both LTRs. A site
   in U3 exists in both only if the vector is not ΔU3 SIN — check, as section 1 does, rather than
   assuming.
5. **Decide what failure you could not afford to miss, and run that assay.** Section 6 orders
   them. For a pooled library the answer is usually titre plus representation, because those are
   the two that fail quietly.

### The general lesson

Of the eight sites in pLVX, four were cheap, two were regulatory sequence nobody has measured,
and two were in the most conserved six bases of the vector. The expensive ones were not expensive
because the plasmid was badly designed — they were expensive because they sit in HIV-1 sequence
that every lentiviral vector carries. **Changing the enzyme is almost always cheaper than
changing the virus.** That is the finding to carry forward.

## 8. Sources

Every source below was read **2026-10-06**. Downloads are under
`reference_docs/synthesis_and_assembly/lentiviral-tolerance/`, whose `README.md` says which file
is which and what was left out.

### Sequence

- GenBank **K03455.1**, HIV-1 HXB2 reference genome, 9,719 bp.
  <https://www.ncbi.nlm.nih.gov/nuccore/K03455.1>
- GenBank **AF324493.2**, HIV-1 vector pNL4-3, 14,825 bp.
  <https://www.ncbi.nlm.nih.gov/nuccore/AF324493.2>
- Addgene 171123, pLVX-TetOne-Puro-GFP, Addgene-verified sequence 336792, already held for #230.
  <https://www.addgene.org/171123/>

### Toolkits

- Fonseca, J. P. *et al.* A toolkit for rapid modular construction of biological circuits in
  mammalian cells. *ACS Synth. Biol.* **8**, 2593–2606 (2019). PMID 31686495,
  [doi:10.1021/acssynbio.9b00322](https://doi.org/10.1021/acssynbio.9b00322).
- Sanjana, N. E., Shalem, O. & Zhang, F. Improved vectors and genome-wide libraries for CRISPR
  screening. *Nat. Methods* **11**, 783–784 (2014). PMID 25075903,
  [doi:10.1038/nmeth.3047](https://doi.org/10.1038/nmeth.3047).
- Martella, A. *et al.* EMMA: an extensible mammalian modular assembly toolkit. *ACS Synth. Biol.*
  **6**, 1380–1392 (2017). PMID 28418644,
  [doi:10.1021/acssynbio.7b00016](https://doi.org/10.1021/acssynbio.7b00016).
- Haellman, V. *et al.* VAMSyB. *Metab. Eng.* **66**, 41–50 (2021). PMID 33857582,
  [doi:10.1016/j.ymben.2021.04.003](https://doi.org/10.1016/j.ymben.2021.04.003). **Published
  CC-BY PDF read in full**, from the ETH Research Collection copy at
  [doi:10.3929/ethz-b-000480636](https://doi.org/10.3929/ethz-b-000480636). The supplement was
  not retrieved — it is served only from ScienceDirect, which refuses an automated fetch.
- Addgene-verified full sequences of [pTS1106_Tier3(Lenti)](https://www.addgene.org/169661/)
  (Addgene 169661) and [lentiCRISPR v2](https://www.addgene.org/52961/) (Addgene 52961), both
  read 2026-10-06 through the public sequence viewer.

### TAR

- Klaver, B. & Berkhout, B. *EMBO J.* **13**, 2650–2659 (1994). PMID 8013464, PMC395139,
  [doi:10.1002/j.1460-2075.1994.tb06555.x](https://doi.org/10.1002/j.1460-2075.1994.tb06555.x).
  **Abstract only** — the full text was blocked this run.
- Das, A. T., Klaver, B. & Berkhout, B. *J. Virol.* **72**, 9217–9223 (1998). PMID 9765469,
  PMC110341, [doi:10.1128/jvi.72.11.9217-9223.1998](https://doi.org/10.1128/jvi.72.11.9217-9223.1998).
  Read as web-page full text, **not** as the PDF; no figure panel was seen.
- Das, A. T., Vrolijk, M. M., Harwig, A. & Berkhout, B. *Retrovirology* **9**, 59 (2012).
  PMID 22828074, PMC3432602, [doi:10.1186/1742-4690-9-59](https://doi.org/10.1186/1742-4690-9-59).
  Full text saved.
- Vrolijk, M. M., Harwig, A., Berkhout, B. & Das, A. T. *Retrovirology* **6**, 13 (2009).
  PMID 19210761, PMC2645353, [doi:10.1186/1742-4690-6-13](https://doi.org/10.1186/1742-4690-6-13).
  Full text saved.
- Berkhout, B., Vastenhouw, N. L., Klasens, B. I. F. & Huthoff, H. *RNA* **7**, 1097–1114 (2001).
  PMID 11497429, PMC1370158,
  [doi:10.1017/s1355838201002035](https://doi.org/10.1017/s1355838201002035). **Abstract and
  partial page only.**
- Da Silva Amaral, C. *et al.* *J. Virol.* (2026). PMID 42053302, PMC13185573,
  [doi:10.1128/jvi.01840-25](https://doi.org/10.1128/jvi.01840-25). Full text saved.

## 9. What this note did not settle

Collected here so the next reader sees the whole hole at once. Sections 2, 3, 4 and 6 each repeat
their own entries in place.

| Open | What would close it |
| --- | --- |
| Pol II initiation, +1 and capping as a function of R nt 1–6, under a heterologous promoter | The change made in a real backbone, then producer-cell vector RNA by RT-qPCR, the 5' end by 5'-RACE, and functional titre |
| Functional titre of **any** lentiviral vector with modified or deleted TAR/R | A side-by-side titre against the parent |
| A designed compensatory lower-stem double mutant in a **vector** | Build it and titre it |
| How many 5'/3' R mismatches first strand transfer tolerates | The Berkhout 2001 *in vitro* donor/acceptor assay with a changed donor and wild-type acceptor |
| Whether a BsaI-based lentiviral toolkit avoided BsaI for the reason section 1 gives | A methods section or correspondence from the MTK or lentiCRISPR groups |
| ~~Whether VAMSyB already removed BsaI from a lentiviral LTR~~ — **closed 2026-10-06, issue #262.** They did: C→T at R+4 in both copies. They reported no titre, so this closes the precedent question and leaves the titre gap above exactly as it was | — |
| Whether any Tet-on lentiviral backbone is already BsaI-clean | Scan the Addgene sequences of TRE lentiviral vectors — no login needed, see section 2 |
| Whether pLVX-TetOne is intended to be self-inactivating | The manufacturer's map, or an LTR-transcription assay on the provirus |
| The effect of base changes in the hPGK promoter, where two of pLVX's sites sit | A reporter assay for promoter strength before and after |

Five papers central to section 3 were blocked by PMC's bot wall and read as abstracts or
web-page text rather than PDFs: PMC110341, PMC395139, PMC1370158, PMC190280 (Harrich 1996, not
read at all) and PMC191344. **A follow-up pass should retrieve them.** The one most likely to
close the titre gap was found but not read: Wright *et al.* (2026) on "supA-LTRs", sequence-
upgraded polyA LTRs in real lentiviral vectors — PMC13370168, open access.

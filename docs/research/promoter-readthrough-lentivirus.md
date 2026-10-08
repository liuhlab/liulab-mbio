---
search:
  exclude: true
---

# Read-through and promoter interference between two promoters inside the LTRs

Research note for issue #231, under #225. Everything below was read on **2026-10-06**. It answers
what `docs/research/ap1-demo-project.md` §7.3 left open: what a primary source measures when a
lentiviral transfer vector carries a constitutive promoter upstream of an inducible one, and what
that does to the inducible promoter's background without doxycycline.

The arrangement at issue is the one §7.1 committed to: between the LTRs, a constitutive promoter
driving Tet3G, puromycin resistance and mCherry, then a TRE promoter driving the cargo.

## 1. The verdict, first

**Read-through from an upstream internal promoter into a downstream transcription unit is
measured in a lentiviral vector, and it contributes to the downstream unit's expression. No
source read measures it into an inducible promoter, and none reports an uninduced background
figure for this arrangement.**

Three things are measured and three are not.

Measured:

- Transcription from an upstream internal promoter reaches the downstream transcription unit in a
  SIN lentiviral vector. Cutting it with a polyadenylation signal between the two cassettes
  *lowered* downstream expression (Tian and Andreadis 2009).
- Two internal promoters in one vector interfere, bidirectionally, and the interference is large
  (Curtin et al. 2008; Tian and Andreadis 2009).
- Which upstream constitutive promoter drives the regulator changes how much a downstream
  inducible promoter leaks, and not through the regulator's own level (Benabdellah et al. 2016).

Not measured anywhere read:

- The uninduced background of a TRE promoter sitting downstream of a constitutive promoter,
  against a matched control with the upstream promoter removed. That is the number this build
  wants and nobody has published it.
- The reverse arrangement §7.2 would need — a constitutive promoter downstream of an induced one
  in a lentivirus. Nothing read puts an inducible promoter upstream.
- Any of this with Tet3G or another rtTA transactivator in a third-generation vector. The one
  all-in-one lentiviral measurement read uses the TetR repressor, which fails in the other
  direction.

Section 5 says what experiment closes it.

## 2. What each source measures

### 2.1 Tian and Andreadis 2009 — read-through, in a lentivirus, with a number attached

*System.* A self-inactivating lentiviral vector carrying two internal transcription units in
tandem: PGK driving ZsGreen, then CMV driving DsRed (`P-G-C-R`). Both cassettes were cloned
**antisense to the viral LTRs**, deliberately, so that polyadenylation signals could be put
between them without terminating the packaged genome in the producer cell. Read out by flow
cytometry in human hair-follicle stem cells, A431 cells, bone-marrow mesenchymal stem cells and
primary keratinocytes, at matched transduction and copy number.

*What interference cost.* In the untreated `P-G-C-R` vector the PGK–ZsGreen unit was nearly
silent: mean green fluorescence 13.7 ± 0.4, against 1688.3 ± 17.3 for the same unit in the
optimised vector — over 120-fold. DsRed, the downstream unit, lost far less: 128.8 ± 1.7 against
177.8 ± 2.6, about 40%. The authors also saw cells expressing DsRed but not EGFP in every cell
type tested, so a single vector gave one gene and not the other.

*The read-through measurement.* Inserting a polyadenylation signal between the two cassettes
raised the upstream gene and **lowered the downstream one**: "while all poly(A) sequences
increased GFI significantly, no poly(A) increased red fluorescence intensity (RFI) and some of
them (rGBpA and SPA) decreased RFI compared to control lentivirus". Their explanation of their
own result: "one possible explanation may be reduction of read-through transcription by the PGK
promoter."

That is the load-bearing fact for this build. Part of what the downstream unit expressed was
not its own promoter firing; it was polymerase arriving from the upstream promoter. The
mechanism is the authors' reading of the measurement, not a separate experiment — they did not
sequence a transcript across the junction.

*Mitigation, measured.* Stacked between the two cassettes, each over the one before:

| Element, between the cassettes | Upstream gene | Downstream gene |
| --- | --- | --- |
| SPA polyadenylation signal + f1 spacer | highest of the polyA set | rGBpA and SPA lowered it |
| Tactb pause site after SPA | +20.0 ± 6.1% | — |
| cHS4 insulator after Tactb | +41.8 ± 0.8% | +9.0 ± 3.6% |
| sMAR8, downstream of the second polyA | +14.7 ± 3.2% | +30.9 ± 6.9% |

MAZ4 gave little gain upstream and lowered the downstream gene; CoTC and Tmsa were
worse than an inert spacer. The finished vector (SPA, Tactb, cHS4 between the units,
sMAR8 after the second polyA) expressed both genes as well as single-cassette vectors did, in
every cell type, and the downstream gene's level did not move when the upstream promoter was
swapped — which is what "interference eliminated" means here.

### 2.2 Curtin et al. 2008 — bidirectional, in a late-generation lentivirus

*Read as abstract only.* Gene Therapy is paywalled and no open copy was found; the full text was
not read and no number from it is quoted here.

The abstract states that interference between two internal heterologous promoters — human CMV
immediate-early and human EF-1α — in a late-generation lentiviral vector "occurred bidirectionally
with both promoters markedly impairing expression of the adjacent transcription unit". It names
no magnitude and no order effect that can be checked without the full text.

Cite it for one thing: the interference is not one-way, so putting the weaker promoter downstream
does not make the problem someone else's.

### 2.3 Eszterhas et al. 2002 — whether order matters, and the causal test

*System.* Not a lentivirus. Two reporter units, GFP and YFP, each with its own promoter, enhancer
and a strong polyadenylation signal, placed by recombinase-mediated cassette exchange into two
tagged chromosomal sites (RL5 and RL6) in mouse erythroleukaemia cells, so the comparison is one
arrangement against a single-gene construct at the *same* site. Both promoters constitutive.

*Arrangements, measured.*

| Arrangement | Result |
| --- | --- |
| Convergent | Expression "barely above the background level of fluorescence" |
| Divergent | Three- to six-fold reduction in each gene |
| Tandem | "Most strongly influenced by the integration site and the genes' orientation" |

Within tandem, the direction was not fixed. At RL5 in one orientation "the upstream gene almost
completely suppresses the downstream gene … but the upstream gene is itself also suppressed to
<50% of the single gene"; in the other orientation at the same site both fell and the upstream
gene ended up lower than the downstream one. At RL6 in one orientation the upstream gene matched
a single gene while the downstream gene was strongly suppressed.

*The causal test.* "Deletion of the upstream promoter fully restored expression of the downstream
gene to levels that are not significantly different from those with a single gene." Transcription
from the upstream unit, not its mere presence, is what suppressed the downstream one.

*What it does and does not say for us.* It answers the order question — order matters, and in
tandem the common pattern is the upstream unit suppressing the downstream one — but it adds the
qualifier that matters most for a pooled library: the size and even the sign vary with the
integration site. Every member of a lentiviral library integrates somewhere different, so leak is
a distribution across cells, not one number. It measures *suppression* of a constitutive
downstream gene. It does not measure *elevation* of a repressed one, which is our failure mode.

### 2.4 Benabdellah et al. 2016 — the closest thing to our arrangement

*System.* All-in-one Tet-On lentiviral vectors. Inside the LTRs: an internal constitutive
promoter — SFFV or human EF-1α — driving the TetR repressor, then a CMV–TetO inducible promoter
driving eGFP. So: constitutive promoter upstream, inducible promoter downstream, one packaged
transcript. Measured in 293T, human mesenchymal stromal cells and two human embryonic stem cell
lines, by flow cytometry, with leak defined as `%eGFP+(−Dox) × 100 / %eGFP+(+Dox)`, from 0 to 100.

*The upstream promoter's identity changes the downstream promoter's leak.* Swapping only the
upstream constitutive promoter, SFFV for hEF-1α, lowered leak and lowered induced output (293T
whole-cell mean fluorescence 4652 for SFFV against 1976 for hEF-1α).

The control that makes this interesting: in mesenchymal stromal cells the hEF-1α vector expressed
*less* TetR than the SFFV vector and still leaked less. The authors conclude "the lower leaking of
the CEETnl2 is not due to higher expression levels of the repressor but to the Lent-On-Plus
backbone configuration". Something other than repressor dose — the configuration of the two
promoters — sets the downstream promoter's background.

Do not stretch that further than they did. They name no mechanism, do not say "read-through", and
run no transcript assay.

*Mitigation, measured.* The Is2 insulator — the chicken β-globin HS4-650 element joined to a
synthetic 388 bp SAR, built in Benabdellah et al. 2014 and placed in the 3′ LTR, so it is
duplicated on integration and flanks the provirus.

| Cells | Leak, without Is2 | Leak, with Is2 | Fold induction, without → with |
| --- | --- | --- | --- |
| 293T, fully responsive population | 66.2 | 7.4 | 5.7 → 38.8 |
| hMSC, sorted, one copy per cell | 54 | 2.6 | 2.8 → 14 |
| hESC, sorted | 4.7 | 1.3 | 5.5 → 5.2 |

The cost they report: titre down 2–3 fold against the uninsulated vector.

*The caveats, which are large.* This is the TetR repressor, not an rtTA transactivator — TetR sits
on the operator and is displaced by doxycycline, so its leak is a failure of occupancy, while a
Tet3G vector's leak is a minimal promoter firing without its activator. And Is2 flanks the
provirus; it does not sit between the two promoters, so what it was measured to fix is chromatin
reaching in, not polymerase arriving from the cassette next door. It is evidence that this
arrangement leaks and that leak is tunable. It is not a measurement of read-through.

## 3. Promoter order: what is settled

- Interference runs both ways between two internal promoters in a lentiviral vector (Curtin 2008).
- In tandem, the upstream unit suppressing the downstream one is the common pattern, but the
  magnitude and the sign depend on the integration site (Eszterhas 2002).
- Transcription from the upstream promoter reaches the downstream unit and contributes to its
  output (Tian and Andreadis 2009).

Nothing read tests an inducible promoter upstream of a constitutive one, which is the comparison
§7.2's rejected placement would need. On the evidence above, putting the marker's constitutive
promoter downstream of the TRE does not avoid the problem — it trades a background risk for a
marker-silencing risk, and no source read measures either side of that trade.

## 4. Mitigation, ranked by what was measured

1. **Both internal cassettes antisense to the LTRs, with a polyadenylation signal between them.**
   Tian and Andreadis's answer to §7.2's objection that a terminator truncates the packaged
   genome: flip the cassettes, and the internal polyA is on the strand the producer cell does not
   package. Measured to restore the suppressed unit and to cut the read-through contribution to
   the downstream one.
2. **A pause site after the polyadenylation signal.** Tactb was best of four tested;
   MAZ4, CoTC and Tmsa were not worth their length in that vector.
3. **cHS4 between the cassettes.** Gained both genes.
4. **An insulator flanking the provirus** — Is2 in the 3′ LTR, or sMAR8 after the second polyA.
   Large effect on leak in a Tet vector (66.2 → 7.4 in 293T), at 2–3 fold titre.
5. **Choosing a weaker upstream constitutive promoter.** hEF-1α for SFFV lowered leak and lowered
   induced output with it. Measured, and not free.

None of 1 to 3 has been measured with an inducible downstream promoter.

## 5. What is still open, and the experiment that closes it

Open: how much of the committed vector's uninduced cargo expression comes from the upstream
constitutive promoter, in the cells this build uses, at one copy per cell.

The experiment, as small as it goes:

1. Build the committed vector and one matched control in which the upstream constitutive promoter
   alone is deleted, leaving the cassette's length and sequence otherwise intact. A length-matched
   promoterless control, not a shorter vector — Eszterhas's deletion control is the precedent.
2. Transduce at a multiplicity that keeps transduced cells under 30%, so most carry one
   integration, and confirm copy number.
3. Without doxycycline, measure cargo expression by flow cytometry inside the mCherry-positive
   gate. The difference between the two vectors is the read-through contribution. Report it as the
   distribution across single-copy cells, not a mean — Eszterhas's integration-site result says
   the spread is the finding.
4. Confirm the mechanism rather than inferring it: strand-specific RT-PCR or long-read sequencing
   across the TRE, showing transcripts that start in the upstream cassette and continue past the
   inducible promoter.
5. If the contribution is material, test mitigation 1 above — cassettes antisense to the LTRs with
   SPA plus Tactb between them — against the same control, and check titre.

Until step 3 has a number, the build should treat uninduced background as unquantified and
design around it: carry a no-doxycycline arm for every screen, and read each library member's
phenotype against its own uninduced well rather than against a pooled baseline.

## 6. What §7 should say now

§7.3 can stop saying "no source read measures read-through or promoter interference for either".
Three do, and §2 above is what they measure. What stays open is narrower and worth stating as
such: no source measures the uninduced background of a TRE promoter downstream of a constitutive
one. The backbone placement §7.1 chose is still the arrangement with a track record, and that
record now includes a measured cost — interference is bidirectional and read-through is real — and
a measured set of fixes nobody has yet tried on an inducible promoter.

## Sources

- Tian, J. and Andreadis, S.T. (2009) Independent and high-level dual-gene expression in adult
  stem-progenitor cells from a single lentiviral vector. *Gene Ther.* 16, 874–884.
  [doi:10.1038/gt.2009.46](https://doi.org/10.1038/gt.2009.46); full text read via Europe PMC,
  [PMC2714872](https://europepmc.org/article/MED/19440229) (author manuscript)
- Curtin, J.A., Dane, A.P., Swanson, A., Alexander, I.E. and Ginn, S.L. (2008) Bidirectional
  promoter interference between two widely used internal heterologous promoters in a
  late-generation lentiviral construct. *Gene Ther.* 15, 384–390.
  [doi:10.1038/sj.gt.3303105](https://doi.org/10.1038/sj.gt.3303105). **Abstract only** — closed
  access, full text not read
- Eszterhas, S.K., Bouhassira, E.E., Martin, D.I.K. and Fiering, S. (2002) Transcriptional
  interference by independently regulated genes occurs in any relative arrangement of the genes
  and is influenced by chromosomal integration position. *Mol. Cell. Biol.* 22, 469–479.
  [doi:10.1128/MCB.22.2.469-479.2002](https://doi.org/10.1128/MCB.22.2.469-479.2002), full text
  via [PMC139736](https://pmc.ncbi.nlm.nih.gov/articles/PMC139736/)
- Benabdellah, K., Muñoz, P., Cobo, M., Gutierrez-Guerrero, A., Sánchez-Hernández, S.,
  Garcia-Perez, A., Anderson, P., Carrillo-Gálvez, A.B., Toscano, M.G. and Martin, F. (2016)
  Lent-On-Plus lentiviral vectors for conditional expression in human stem cells. *Sci. Rep.* 6,
  37289. [doi:10.1038/srep37289](https://doi.org/10.1038/srep37289), full text via
  [PMC5112523](https://pmc.ncbi.nlm.nih.gov/articles/PMC5112523/) (CC BY 4.0)
- Benabdellah, K., Cobo, M., Muñoz, P., Toscano, M.G. and Martin, F. (2011) Development of an
  all-in-one lentiviral vector system based on the original TetR for the easy generation of Tet-ON
  cell lines. *PLoS ONE* 6, e23734.
  [doi:10.1371/journal.pone.0023734](https://doi.org/10.1371/journal.pone.0023734), full text via
  [PMC3158098](https://pmc.ncbi.nlm.nih.gov/articles/PMC3158098/) (CC BY) — the parent CEST vector
  of the 2016 study
- Benabdellah, K., Gutierrez-Guerrero, A., Cobo, M., Muñoz, P. and Martín, F. (2014) A chimeric
  HS4-SAR insulator (IS2) that prevents silencing and enhances expression of lentiviral vectors in
  pluripotent stem cells. *PLoS ONE* 9, e84268.
  [doi:10.1371/journal.pone.0084268](https://doi.org/10.1371/journal.pone.0084268) — abstract read,
  for what Is2 is and where it sits

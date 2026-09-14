---
search:
  exclude: true
---

# Barcode design: does a set need a GC band and a homopolymer cap?

Research note for issue #58. It settles one question for the barcode designer of issue #62:
beyond a minimum distance **within a part list** and freedom from a forbidden enzyme site,
should a barcode set also be constrained on GC content and on homopolymer run length?

The question is not academic. Whatever this note concludes becomes the default dials of a
designer that every later user inherits, so an unsourced default is worse than an absent one.

Everything below was read or computed on **2026-09-14**. Every count over the published barcode
set was computed from `Table S1.xlsx` on that date; every count over the 11-mer space was
computed by the arithmetic in section 7, which names the construction that produced it.

## 1. The verdicts, first

| Rule | Verdict | The number | Strength of the evidence |
| --- | --- | --- | --- |
| **GC band** | **Off by default** | none | **Convention**, and one measurement against it |
| **Homopolymer cap** | **On by default** | longest run **≤ 5** | **Measured mechanism**; the threshold is a judgement |
| Minimum distance | On, unchanged | Hamming **3** within a part list | Measured on the published set |

**The GC band is off because no measurement supports one at this length.** Every numeric GC
band in the barcode literature traces back to a single uncited sentence in Hamady et al. 2008.
The one paper that imposed such a band and then measured its consequences, Xu et al. 2009,
found GC was not the variable that mattered. The one direct measurement at a comparable tag
length, Kebschull and Zador 2015, found amplification flat out to 80% GC. Illumina and PacBio
document no GC band for a barcode at all. And the real, large GC bias that does exist, Aird
et al. 2011, is a property of the whole fragment, which an 11-mer moves by one to five points.

**The homopolymer cap is on because it covers a failure the distance rule cannot see.** A
minimum Hamming distance protects against substitutions only. The dominant residual error of
PacBio HiFi reads, which is how this barcode block is actually read, is an **indel inside a
homopolymer** — 92.0% of read discordances, one every 477 bp, measured by PacBio's own authors.
A Hamming code offers no protection against that whatever its distance, so the cap and the
distance rule guard different modes, and the cap is the only one of the two that touches this
one. It also costs nothing measurable (section 7).

**The cap is 5, not the conventional 2 or 3.** The conventional value is 454-era, and the
software that ships it says outright that it is unnecessary on Illumina. 5 is the one number
here the evidence does not pin; section 5.3 says exactly how it was chosen and what argues
against it.

**A third finding matters more than either rule.** Neither filter fixes the real failure of a
naive designer. A greedy construction that walks the 11-mer space in lexicographic order emits
`AAAAAAAAAAA` as its first barcode and a run of A-tailed, GC-poor sequences after it — and a
homopolymer cap does **not** repair this, it merely caps the run. Drawing candidates in a
seeded shuffled order does repair it, costs nothing, and keeps the result reproducible.
Section 9 shows all three sets side by side. Issue #62 should take that as the default
construction.

## 2. Licences, and what may ship

Nothing in this note is copied into package data. What ships is a **default value chosen by
this repository**, with the evidence recorded here; no vendor's or author's table is
redistributed.

| Source | Licence | Can we ship it? |
| --- | --- | --- |
| Takacsi-Nagy et al. 2026, *Cell* | CC BY 4.0, verified in `protein-library-assembly.md` | **Yes** — quote and reuse with attribution |
| Faircloth & Glenn 2012; Bystrykh 2012, *PLoS ONE* | CC BY | **Yes** — with attribution |
| Buschmann & Bystrykh 2013, *BMC Bioinformatics* | CC BY | **Yes** — with attribution |
| Aird et al. 2011, *Genome Biology* | CC BY | **Yes** — with attribution |
| Hamady et al. 2008, *Nat Methods*; Buschmann 2017, *Bioinformatics* | subscription; licence not verified here | **Quote briefly and cite**; do not vendor |
| Wenger et al. 2019, *Nat Biotechnol*; Xu et al. 2009, *PNAS*; Kebschull & Zador 2015, *NAR* | licence not verified here | **Quote briefly and cite**; do not vendor |
| Illumina and PacBio documentation | proprietary, all rights reserved | **Quote briefly with attribution**; never redistribute |
| DNABarcodes R package | GPL-2 (Bioconductor) | **No code may be copied.** Its defaults were read from source and are reported as fact, not reused |

The DNABarcodes row is the one that needs care. Reading a GPL package's source to learn what
its defaults **are** is not copying it, and this note states those defaults as a fact about the
software. No line of it may be translated into `src/`.

## 3. The one set that was measured

The strongest single piece of evidence available is a published barcode set that holds to a
distance rule and to nothing else, and that was then sequenced on two platforms. It bounds how
bad loose can be. It does not show that loose is optimal, and it is not evidence that the
authors considered the question — the paper never states that a GC or homopolymer rule was
weighed and rejected.

The source is Takacsi-Nagy, O. et al. (2026), *Cell* 189, 1-20,
[doi:10.1016/j.cell.2026.07.054](https://doi.org/10.1016/j.cell.2026.07.054), CC BY 4.0.
`docs/research/protein-library-assembly.md` records the method; this note reuses its barcode
measurements and adds the ones issue #58 needs.

### 3.1 What the set holds

| Property | The three part lists (72) | Every 11-mer in the workbook (164) |
| --- | --- | --- |
| Barcode length | 11 bp | 11 bp |
| GC content | 2 to 9 of 11, **18.2% to 81.8%** | 2 to 9 of 11, **18.2% to 81.8%** |
| Longest homopolymer run | **6** | **6** |
| Minimum Hamming, within a part list | **3** (N), **4** (bZIP), **4** (C) | n/a |
| Minimum Hamming, across part lists | **2** | **2** |
| Carries a BsaI, BbsI, SrfI or PmeI site | none | none |

So the set satisfies a distance rule **within a part list** and site-freedom, and no GC band
and no homopolymer cap was applied to it.

**The distance rule is per part list, and the paper overstates it.** METHOD DETAILS, p. e7,
claims "a minimum Hamming distance of three between any two domain BCs". Measured, that holds
within each part list but not across them: the minimum across the three lists is 2. Issue #62
should check distance within a part list, which is also the only place it is needed.

### 3.2 Which barcodes sit at the extremes

This matters more than the range does, because the identity of the extreme barcodes decides
which vendor concern the set actually tested.

| Barcode | Part list | Run | Base | GC |
| --- | --- | --- | --- | --- |
| `GCTTTTTTGGC` | bZIP | **6** | T | 5/11 |
| `TCATTGGGGGC` | N | **5** | G | 7/11 |
| `CACGGCAAAAT` | N | 4 | A | 5/11 |
| `TCTGTGCCCCG` | N | 4 | C | 8/11 |
| `ACGGCGTTTTA` | C | 4 | T | 5/11 |
| `CCTATAAGGGG` | C | 4 | G | 6/11 |
| `GAACCGAAAAG` | C | 4 | A | 5/11 |
| `TGCGGGAAAAT` | C | 4 | A | 5/11 |

**The longest run in the set is six T, not six G.** The longest G run anywhere in the 164 is
**5**, in two barcodes. That is the number any argument from two-colour chemistry has to be
measured against, and it is a weaker test than a set containing a longer G run would have been.

The seven barcodes at the GC extremes are `AACAAATTAGT`, `AATTCTTATTC`, `CTAAATTAATG` and
`TAATCACAAAT` at 2 of 11, and `CGAACCGGGCC`, `GAGGCACCGGG` and `TGGGCGGCCTG` at 9 of 11.

### 3.3 What the set achieved, verified against the paper

Issue #58 quotes roughly 95% of reads with three valid barcodes and roughly 90% correct
linkage. **Both check out, with two wrinkles.** The figures were read off the article and its
supplement rather than repeated from the ticket.

Main text, p. 3:

> Of reads containing three valid BCs (∼95%), nearly 90% of the library was correctly linked

Figure S1B prints the read categories as 84.1% correct, 8.3% barcode mismatches, 1.0% missing
domains, 4.5% missing barcode and 2.0% other. Its legend says only "Missing BC" and "Other" had
fewer than three barcodes, so reads carrying three valid barcodes are 100 − 4.5 − 2.0 =
**93.5%**, against the text's "∼95%". Among those, correct linkage is 84.1 / (84.1 + 8.3) =
**91.0%**, against the text's "nearly 90%" — and the Figure S1D legend puts it the other way
round:

> Over 90% of assembled DESynR AP-1 TFs had a correct linkage of BCs across all three domains.

So "∼95%" is really 93.5%, and the main text's "nearly 90%" is conservative against its own
supplement's "over 90%". Both accounts agree to a point or two, and the ticket's figures stand.

Two further measurements bear on dropout, which is the failure a GC rule would be imposed to
prevent:

> NGS of the plasmid library showed 99.8% of possible domain combinations were present and
> evenly distributed (∼89% within 1 log % of reads)

and, separately on the same page:

> barcoding errors did not correlate with TF size or cause dropout: >90% of the library
> maintained >80% fidelity

**The set was read on both platforms, which is why both vendors are in scope.** Barcode
abundance for the screens was counted from Illumina amplicons on "MiSeq or NextSeq instruments"
(METHOD DETAILS, p. e7) — NextSeq is two-colour chemistry. The linkage fidelity above was
measured on a PacBio **Revio** from HiFi reads (p. e9). So the 11-mers carrying a G5 and a T6
went through both a two-colour instrument and a homopolymer-sensitive long-read instrument, and
the aggregate result was 93.5% and 91.0%.

**What this does not establish.** The 8.3% barcode-mismatch rate is not nothing, and nothing
here says whether those mismatches are spread evenly or concentrated in the GC-extreme and
homopolymer-carrying barcodes. That is the measurement that would settle issue #58 outright,
and it was not available — see `## Open gaps`.

## 4. What the barcode-design literature recommends, and what it measured

Six primary sources were read in full text. **Every one of them states its GC band and its
homopolymer cap as convention.** Not one measured either.

| Source | GC band | Homopolymer cap | Justification |
| --- | --- | --- | --- |
| Hamady et al. 2008 | 40-60% | no run of 3 | **Convention**, one uncited sentence |
| Faircloth & Glenn 2012 | 40-60% | runs > 2 excluded | **Convention**, self-described "subjective" |
| Bystrykh 2012 | none set | none set | discusses; sets no cutoff |
| Buschmann & Bystrykh 2013 | 40-60% | ≤ 2 repeats | **Convention**, cites only a general claim |
| Buschmann 2017 / DNABarcodes | 40-60%, on by default | no triplets, on by default | **Convention**, "advised in the available literature" |
| Hawkins et al. 2018 (FREE) | 40-60% | no triples | **Convention**, "experimental considerations" |

### 4.1 The 40-60% band traces to one uncited sentence

Hamady et al. 2008 is where it starts. Their entire stated rationale:

> To pick our maximal set of 1544 codewords (Supplementary Data), we chose an encoding scheme
> for ATCG that resulted in the most valid 'candidate' codewords, then filtered these
> candidates to optimize PCR and sequencing performance based on GC content (40-60%), and
> eliminating consecutive triples of the same base and self-complementarity or complementarity
> to the primer.

No reference marker, no data, no citation. "To optimize PCR and sequencing performance" is the
whole of the argument. The band then reaches Bystrykh 2012, Buschmann and Bystrykh 2013,
DNABarcodes and Hawkins et al. 2018 unchanged.

Worth noting because it is often asserted otherwise: Hamady never names 454 homopolymer indels
as the reason for the no-triples rule. Later authors supply that link.

**On an 11-mer the inherited band is severe.** 40-60% admits only 5 or 6 G/C out of 11 —
4 of 11 is 36.4% and fails, 7 of 11 is 63.6% and fails. Applied to the published set it would
reject 37 of the 72 part-list barcodes.

### 4.2 The authors themselves call the filters subjective

Faircloth and Glenn 2012 state the rule:

> We did not include, in any set, sequence tags having >2 homopolymers, GC content outside the
> range 40%<GC<60%, or perfect self-complementarity.

and then describe the whole class of them:

> filtering tags based on subjective or platform-specific criteria

In their software the filters are optional, not default. Bystrykh 2012 sets no numeric cutoff
at all, treats GC as mainly a microarray concern, and says of homopolymers:

> Stringency of filtering depends on the application and sequencing chemistry. Roche 454
> pyrosequencing chemistry for instance, imposes serious restriction on this parameter. Solexa
> platforms and their recent upgrades seem to be very robust in reading homopolymers.

DNABarcodes turns all three filters on by default — `filter.triplets`, `filter.gc` and
`filter.self_complementary`, with GC kept inclusive between 0.4 and 0.6 — on the stated ground
that "in the available literature, use of all three mentioned filters is advised and therefore
enabled by default". Its own vignette then undercuts one of them:

> Some of these filters may be unnecessary on current or future platforms. For example, with
> the Illumina's Sequencing by Synthesis technology, the triplet filter is unnecessary.

### 4.3 What the papers did measure

Where these authors measured tag performance, they found no GC effect — and in two cases found
the bias lay somewhere else entirely.

- **Faircloth and Glenn 2012**, 95 primers carrying 10-nt tags: "the average fold-difference of
  read counts per sequence tagged primer did not differ from one ... suggesting that
  incorporation of edit metric sequence tags to primers did not affect amplification or
  sequencing." Their pooled **adapters** did vary, by over two orders of magnitude, which they
  attribute to "differences in ligation efficiency among individual adapters" — not to
  composition.
- **Xu et al. 2009** imposed a 40-60% band on 240,000 probes, then measured which variable
  separated good probes from bad. G+C% by performance group was 49.2 (good), 50.9 (bright),
  50.4 (medium) and 46.4 (dim), and "the differences were rather small to account for the
  disparity in their hybridization properties". What did matter was four-base runs — `CCCC`
  enriched among cross-hybridising probes, `AAAA` among dim ones. This is hybridization, not
  PCR or read counts, so it does not transfer; but it is the only paper that put its own GC
  filter to the test, and the filter failed it.
- **Kebschull and Zador 2015** is the closest measurement to this question's scale: 20-nt
  barcodes across the full GC range to 80%, with amplification efficiency essentially flat,
  slopes below 0.004.
- **Alon et al. 2011** separates arrangement from composition. Barcodes introduced by ligation
  before PCR gave "up to 100-fold differences in read counts"; barcodes carried on the PCR
  primers gave none (R² = 0.9977 ± 0.0016).
- **Berry et al. 2011** and **O'Donnell et al. 2016** both measured real barcoded-primer bias
  and attributed it to primer-template interactions and mismatch, with GC held constant. Both
  fixes are protocol changes — a two-step PCR — not barcode filters.
- **McGee et al. 2024** is the one report of a GC correlation, and it is weak and
  sample-dependent, with the extremes worst. That argues for avoiding extremes, not for a band.

## 5. What the sequencing platforms document

### 5.1 Illumina documents no GC band and no homopolymer limit

Searched across Illumina's Knowledge Base, the Index Adapters Pooling Guide v17, the Illumina
Adapter Sequences document and the Indexed Sequencing Overview Guide: **no GC percentage
requirement for an index or barcode is stated anywhere**, and **no maximum homopolymer length
is stated anywhere** — "homopolymer" and "poly-G" return no hits across the Knowledge Base.

The widely repeated "25-75% GC, no runs of four" rule is real but it is **patent claim
language**, not guidance: Illumina's US 11,028,436 B2 claims index sets where "each index
sequence ... has a guanine/cytosine (GC) content between 25% and 75%" and which "excludes any
homopolymers having four or more consecutive identical bases". A patent claim is not a
recommendation to users, and this note does not treat it as one.

### 5.2 G is dark, but the rule that follows is narrower than it is usually given

Illumina's two-colour model is confirmed in its own documentation, and G-is-dark is the
invariant across every variant of the chemistry. The MiniSeq article is the bluntest:
"Guanines (G) are permanently dark."

The premise that "A = red, C = green, T = both, G = neither" is **not** what Illumina
documents; the assignments differ between chemistries and only G-is-dark holds throughout. On
NovaSeq 6000, NextSeq 500/550 and MiniSeq it is T green, C red, A both; under XLEAP-SBS it is
T green, A blue, C both.

The hard failure condition is documented for exactly one place — **the first two cycles of the
Index Read**:

> either of the first two cycles of the Index Read must start with at least one base other
> than G ... If an Index Read starts with two G bases, signal intensity is not generated.

Beyond those two cycles, colour balance is stated as a preference, and Illumina adds that
"higher plexity pools are inherently color-balanced, so any index adapter combinations are
acceptable."

**This rule does not reach the barcodes in question.** These barcodes are not indexes. They sit
inside the amplicon and are read as template — the forward primer anneals in the internal
stuffer directly upstream of the barcode block, the reverse primer in TRAC exon 1. The i5 and
i7 indexes for these runs were "standard Illumina Nextera primers", which are Illumina's own
and already balanced. So a G-run cap on a barcode of this kind cannot be justified from the
index-read rule.

The poly-G artefact is also not what it is often taken for. Illumina's DRAGEN documentation
describes its cause:

> Poly-G artifacts appear on two-channel sequencing systems when the dark base G is called
> after synthesis has terminated. As a result, DRAGEN calls several erroneous high-confidence G
> bases on the ends of affected reads.

That is a read-end artefact of synthesis stopping, not something a G homopolymer in the
template provokes. Trimming for it is on by default. Nothing in Illumina's documentation says a
template G run causes phasing failure or base-call confidence loss.

### 5.3 PacBio: the one measured mechanism, and where the cap comes from

This is where a homopolymer cap earns its place. Wenger et al. 2019, written by PacBio authors,
measured the residual error of HiFi reads:

> The large majority of CCS read discordances are indels in homopolymer contexts: 3.4% are
> mismatches, 4.6% are indels in non-homopolymer contexts, and 92.0% are indels in
> homopolymers.

and, converting those proportions to rates:

> This equates to a mismatch every 13,048 bp in CCS reads, a non-homopolymer indel every
> 9,669 bp, and a homopolymer indel every 477 bp.

Two consequences follow for a barcode that will be read on HiFi, as this one was:

1. **The error is an indel, and a Hamming-distance rule cannot see indels at all.** Distance 3
   corrects one substitution and does nothing for a deletion, which shifts every base after it.
   So the distance rule and the homopolymer cap are not redundant with each other; the cap is
   the only one of the two that addresses the dominant long-read error.
2. **The rate is a function of run length.** Longer runs are worse, monotonically.

PacBio's own shipped barcodes were designed accordingly, though without a number:

> the PacBio barcodes were designed to have no homopolymer stretches, low pairwise sequence
> similarity ... and to be simultaneously distinguishable in both forward and
> reverse-complement orientations

That is a prohibition with no threshold and no GC figure. PacBio documents no maximum
homopolymer length and no GC band for barcodes anywhere. Oxford Nanopore documents none either.

**How the cap of 5 was chosen, and what argues against it.** No source sets a threshold. The
conventional 2 or 3 is 454-era and is explicitly disclaimed for Illumina by the software that
ships it. The only threshold measurement available, Johnson et al. 2023, puts severity much
further out: indel error rises with run length and becomes severe past runs of about 10. Against that, the published set reached 6 and worked. A cap of 5 therefore rejects
exactly one published barcode, `GCTTTTTTGGC`, which is the honest argument against it.

It is chosen anyway, for three reasons. The error mechanism is measured and monotone in run
length. The cost is nothing measurable (section 7). And "6 worked" is an aggregate claim only —
the 8.3% mismatch fraction is unattributed, so no per-barcode result shows that this particular
barcode was correct work. A cap of 6 would fire on nothing in the set and remains defensible;
the difference between them is one barcode in 72 and no measurable space.

## 6. Does an 11-mer's GC measurably bias the amplicon PCR?

Much of the GC-bias literature is about whole libraries or long fragments. An 11-mer is short
enough that the extrapolation needs checking rather than assuming, and two structural facts say
it does not carry over.

**The barcode is not in a primer.** In this scheme the forward primer anneals "directly upstream
of the BC region within the 'Internal Stuffer'" and the reverse primer in TRAC exon 1
(Figure S1G legend), so the barcode block lies **inside** the amplicon and in neither binding
site. The PacBio amplicon likewise spans the barcodes internally. The mechanism behind
barcode-driven amplification bias needs the barcode to be part of the primer, where it changes
that primer's Tm and its annealing; a barcode interior to the amplicon cannot act that way.

**An 11-mer barely moves the amplicon's GC.** Computed: a single barcode swinging across the
full measured range, 2 of 11 to 9 of 11, changes the GC content of the whole amplicon by

| Amplicon | Shift from one barcode at its extremes |
| --- | --- |
| 150 bp | 4.67 points |
| 250 bp | 2.80 points |
| 400 bp | 1.75 points |
| 600 bp | 1.17 points |

The finished barcode block carries three barcodes, so it is the larger unit. Computed over all
24³ = 13,824 published blocks, the 41-bp block runs **34.1% to 78.0%** GC — a swing of 43.9
points across the block, which is 7.20 points of a 250-bp amplicon, 4.50 of a 400-bp one and
3.00 of a 600-bp one. The longest homopolymer run inside a whole block is still 6, so
concatenation creates no longer run than its longest member already carries.

### 6.1 What has actually been measured at this scale

**The large, real GC bias is a whole-fragment effect.** Aird et al. 2011 is the standard
measurement, and it is severe:

> as few as ten PCR cycles using the enzyme formulation (Phusion HF DNA polymerase) and
> thermocycling conditions prescribed in the standard Illumina protocol depleted loci with a GC
> content > 65% to about a hundredth of the mid-GC reference loci. Amplicons < 12% GC were
> diminished to approximately one-tenth

That is a hundredfold, and it is why GC bands exist at all. But it was measured on whole
amplicons of 50-69 bp inside libraries of 180-400 bp, and the paper's own conclusion about
which composition drives it is the decisive sentence for an 11-mer:

> We found the %GC of a 250-bp window centered on the amplicons a better predictor of
> under-coverage than the %GC of the amplicons proper.

Bias tracks the composition of the whole fragment, not of a short stretch inside it. The
arithmetic above puts an 11-mer's contribution at one to five points of that figure, well
inside the flat range — which Aird et al. also found widens with a slower thermocycler ramp,
reaching "from 13% to 84% GC" on the slowest instrument they tried.

**Verdict on this section: the extrapolation does not hold.** The whole-library GC literature
does not transfer to an 11-mer interior to an amplicon, by the measurement of the very paper
that established the bias. Nothing measured supports a narrow GC band at this length, and the
weak evidence that exists points only at the extremes.

## 7. What the rules cost: the 11-mer space

**This is arithmetic run here, not a number cited.** The construction is stated because the
count depends on it.

The space is 4¹¹ = **4,194,304** 11-mers. Counts of codewords below come from a **greedy
construction in lexicographic order**: take allowed 11-mers in increasing index order, accept
one whenever it is at least distance 3 from every codeword already accepted. That is a
**lower bound** on the largest distance-3 code inside each filtered set, never the maximum. The
sphere-packing bound puts the maximum for the unfiltered space at **123,361**, and the
shortened quaternary Hamming code [11,8,3]₄ gives exactly 4⁸ = 65,536 by construction.

| Filter | 11-mers allowed | Greedy distance-3 code |
| --- | --- | --- |
| none | 4,194,304 | **65,536** |
| site-free only | 4,169,216 | 46,692 |
| site-free + run ≤ 6 | 4,165,122 | 48,356 |
| **site-free + run ≤ 5 (the default)** | **4,149,782** | **48,267** |
| site-free + run ≤ 4 | 4,079,404 | 48,155 |
| GC 3-8 of 11 only | 3,919,872 | 47,556 |
| homopolymer run ≤ 3 only | 3,791,232 | 46,164 |
| site-free + GC 3-8 + run ≤ 3 | 3,570,354 | 45,128 |
| site-free + GC 4-7 + run ≤ 3 | 3,003,954 | 41,858 |
| site-free + GC 5-6 + run ≤ 3 | 1,779,682 | 32,353 |

Site-freedom is counted against BsaI, BbsI, SrfI and PmeI in both orientations; only 25,088 of
the 4,194,304 11-mers carry any of those sites.

**Two cautions, both of which the numbers force.**

First, **do not read the unfiltered 65,536 against the rest as the cost of a filter.** That row
is a perfectly structured linear code and greedy finds it exactly; once any filter perturbs the
space the structure is gone and greedy settles near 46,000-48,000 whichever single filter is
applied. That drop is the construction losing its structure, not the filter deleting codewords.

Second, **the homopolymer cap's cost is below the resolution of this method.** Site-freedom
alone yields 46,692 while site-freedom plus a cap of 6 yields 48,356 — *more* codewords from a
*smaller* allowed set. A tighter filter cannot really raise capacity, so the spread across
46,692, 48,356, 48,267 and 48,155 is noise in the greedy construction. The honest statement is
that **the homopolymer cap is free**, not that it gains anything.

**At realistic sizes none of this binds.** A part list holds tens to a few hundred members, and
even the strictest combination in the table leaves 32,353 — more than a thousand times what a
24-member part list needs. The GC and homopolymer rules are cheap in space. The question is
whether they are justified, not whether they are affordable, and the space argument settles
nothing in either direction.

**Under the chosen defaults** — minimum Hamming distance 3, site-free, homopolymer run ≤ 5, no
GC band — the surviving count is **48,267** by lexicographic greedy, and **37,858** by the
seeded shuffled construction recommended in section 9. Both are lower bounds. The second is the
operative number, because it is the construction issue #62 should use.

## 8. What each rule would do to work that succeeded

House rule: a rule shown to fire on correct work is evidence against itself. The published set
worked, in the aggregate sense measured in section 3.3. Computed, here is what each candidate
rule would have rejected from it.

| Candidate rule | Rejects, of the 72 in part lists | Rejects, of the 164 |
| --- | --- | --- |
| GC band 3-8 of 11 (27%-73%) | 7 (9.7%) | 7 (4.3%) |
| GC band 4-7 of 11 (36%-64%) | 18 (25.0%) | 19 (11.6%) |
| **GC band 5-6 of 11 (the inherited 40-60%)** | **37 (51.4%)** | 72 (43.9%) |
| homopolymer cap ≤ 5 (the default) | **1 (1.4%)** | 1 (0.6%) |
| homopolymer cap ≤ 4 | 2 (2.8%) | 2 (1.2%) |
| homopolymer cap ≤ 3 (the inherited value) | 8 (11.1%) | 8 (4.9%) |
| homopolymer cap ≤ 2 | 29 (40.3%) | 59 (36.0%) |

The inherited pair is the damning row. A 40-60% GC band with no run of three would have thrown
out **more than half** of a barcode set that returned 93.5% and 91.0% on two platforms. The
chosen cap of 5 rejects one barcode; no GC band rejects none.

One qualification keeps section 8 from being decisive on its own: "worked" means the set
achieved those figures in aggregate, not that every individual barcode performed. If the 8.3%
mismatch fraction were concentrated in exactly these extreme barcodes, the rules would be
firing on the weakest members rather than on correct work. Nothing available here separates
those two readings.

The five other barcode sets in the same workbook are worth noting, because they are tighter
than the part lists without anyone saying so: the natural AP-1 sheet (30), the barcode reference
sheet (37), the individually synthesised DESynR sheet (39) and the HA GD2 CAR sheet (16) all
sit inside GC 3-7 of 11 with a longest run of 3, and a GC 3-8 band with a cap of 3 rejects
**none** of them. Only the three AP-1 part lists reach the extremes.

## 9. The construction matters more than either rule

This is the finding issue #62 most needs, and no source states it. A greedy designer that walks
the space in lexicographic order produces a degenerate set, and **a homopolymer cap does not
fix it**.

The first 24 codewords of an unfiltered lexicographic greedy, computed:

```text
AAAAAAAAAAA  CCCAAAAAAAA  GGGAAAAAAAA  TTTAAAAAAAA  GCACAAAAAAA  TACCAAAAAAA
ATGCAAAAAAA  CGTCAAAAAAA  TGAGAAAAAAA  GTCGAAAAAAA  CAGGAAAAAAA  ACTGAAAAAAA
```

GC runs 0 to 4 of 11 and the longest run 6 to 11. The very first barcode is `AAAAAAAAAAA`.

Adding the homopolymer cap of 5 does **not** repair this. It caps the run and leaves everything
else:

```text
AAAAACAAAAA  CCCAACAAAAA  GGGAACAAAAA  TTTAACAAAAA  GCACACAAAAA  TACCACAAAAA
```

GC now runs 1 to 5 of 11, every sequence still ends in a tail of A, and every run is exactly
the cap. A GC band would not repair it either; it would simply move the degeneracy to whatever
the band allows.

Drawing the same candidates in a **seeded shuffled order** repairs it completely, at no cost:

```text
ACGGGAGCTCT  CCATGAGGTAG  AACACAAAATG  GACAGACACCA  GCACTATTGTG  ACGGCTCACAT
AATCCCCAGGA  AAGGCTATCGT  TCAACGCCGCC  ATGGCGACCAC  GCAGGTACGGT  GAGGCATCTAG
```

GC runs 2 to 9 of 11 and the longest run is 4, with no filter on either. The whole code holds
**37,858** codewords. A fixed seed keeps it reproducible, which issue #62 requires.

So the composition of a barcode set is mostly a property of how candidates are drawn, not of
what is filtered out afterwards. Recommend the seeded shuffle as the default construction and
the filters stay nearly idle — which is the outcome restraint asks for.

## Open gaps

- **The decisive measurement was not available.** Whether the 8.3% barcode-mismatch fraction is
  enriched in the GC-extreme or homopolymer-carrying barcodes would settle this note's question
  directly. Per-construct read counts are published as Tables S2-S7; three attempts to reach
  them failed (the DOI and the Elsevier link return a redirect shell, and `cell.com` answers
  403), so the correlation was not run. Only `Table S1.xlsx` was held locally.
- **The cap of 5 is a judgement, not a measurement.** No source sets a threshold. 6 would fire
  on nothing in the published set and is equally defensible; 4 and below start to look like the
  inherited convention.
- **The set tests a G run of 5, not a longer one.** Any conclusion drawn from it about
  two-colour chemistry is bounded at that length.
- **Hamming distance does not protect against indels**, which are the dominant HiFi error. A
  Sequence-Levenshtein metric, as in Buschmann and Bystrykh 2013, would. This note keeps Hamming
  because it is what the published set holds and what the paper's own one-mismatch matching
  assumes, but the choice deserves its own ticket rather than being settled here by default.
- **Whether a GC band or homopolymer cap was considered by the paper's authors is not
  recoverable.** Their absence from the set is measured; intent is not.
- **The greedy counts in section 7 are lower bounds**, not maxima, and differences of a few
  percent between rows are construction noise rather than capacity.
- **Two vendor documents were reachable only as third-party mirrors or archives**: the 2014
  two-channel Technology Spotlight, the Illumina Adapter Sequences PDF, and the 2012 PacBio
  barcode technical note, whose live URL now 404s.

## Sources

All read on 2026-09-14.

- Takacsi-Nagy, O. et al. (2026) Synthetic transcription factors designed by domain
  recombination enhance CAR T cell antitumor function. *Cell* 189, 1-20.
  [doi:10.1016/j.cell.2026.07.054](https://doi.org/10.1016/j.cell.2026.07.054). CC BY 4.0.
  Main text p. 3; METHOD DETAILS pp. e7, e9; Figures S1B, S1C, S1D and S1G with their legends;
  `Table S1.xlsx`, from which every measurement of the published set was computed.
- Hamady, M., Walker, J.J., Harris, J.K., Gold, N.J. and Knight, R. (2008) Error-correcting
  barcoded primers for pyrosequencing hundreds of samples in multiplex. *Nat. Methods* 5,
  235-237. [doi:10.1038/nmeth.1184](https://doi.org/10.1038/nmeth.1184)
- Faircloth, B.C. and Glenn, T.C. (2012) Not all sequence tags are created equal: designing and
  validating sequence identification tags robust to indels. *PLoS ONE* 7, e42543.
  [doi:10.1371/journal.pone.0042543](https://doi.org/10.1371/journal.pone.0042543)
- Bystrykh, L.V. (2012) Generalized DNA barcode design based on Hamming codes. *PLoS ONE* 7,
  e36852. [doi:10.1371/journal.pone.0036852](https://doi.org/10.1371/journal.pone.0036852)
- Buschmann, T. and Bystrykh, L.V. (2013) Levenshtein error-correcting barcodes for multiplexed
  DNA sequencing. *BMC Bioinformatics* 14, 272.
  [doi:10.1186/1471-2105-14-272](https://doi.org/10.1186/1471-2105-14-272)
- Buschmann, T. (2017) DNABarcodes: an R package for the systematic construction of DNA sample
  tags. *Bioinformatics* 33, 920-922.
  [doi:10.1093/bioinformatics/btw759](https://doi.org/10.1093/bioinformatics/btw759), with the
  Bioconductor vignette and the package source, from which the default filter values were read.
- Hawkins, J.A. et al. (2018) Indel-correcting DNA barcodes for high-throughput sequencing.
  *PNAS* 115, E6217-E6226.
  [doi:10.1073/pnas.1802640115](https://doi.org/10.1073/pnas.1802640115)
- Xu, Q., Schlabach, M.R., Hannon, G.J. and Elledge, S.J. (2009) Design of 240,000 orthogonal
  25mer DNA barcode probes. *PNAS* 106, 2289-2294.
  [doi:10.1073/pnas.0812506106](https://doi.org/10.1073/pnas.0812506106)
- Aird, D. et al. (2011) Analyzing and minimizing PCR amplification bias in Illumina sequencing
  libraries. *Genome Biol.* 12, R18.
  [doi:10.1186/gb-2011-12-2-r18](https://doi.org/10.1186/gb-2011-12-2-r18)
- Kebschull, J.M. and Zador, A.M. (2015) Sources of PCR-induced distortions in high-throughput
  sequencing data sets. *Nucleic Acids Res.* 43, e143.
  [doi:10.1093/nar/gkv717](https://doi.org/10.1093/nar/gkv717)
- Berry, D. et al. (2011) Barcoded primers used in multiplex amplicon pyrosequencing bias
  amplification. *Appl. Environ. Microbiol.* 77, 7846-7849.
  [doi:10.1128/AEM.05220-11](https://doi.org/10.1128/AEM.05220-11)
- Alon, S. et al. (2011) Barcoding bias in high-throughput multiplex sequencing of miRNA.
  *Genome Res.* 21, 1506-1511.
  [doi:10.1101/gr.121715.111](https://doi.org/10.1101/gr.121715.111)
- O'Donnell, J.L., Kelly, R.P., Lowell, N.C. et al. (2016) Indexed PCR primers induce
  template-specific bias in large-scale DNA sequencing studies. *PLoS ONE* 11, e0148698.
  [doi:10.1371/journal.pone.0148698](https://doi.org/10.1371/journal.pone.0148698)
- Johnson, M.S., Venkataram, S. and Kryazhimskiy, S. (2023) Best practices in designing,
  sequencing, and identifying random DNA barcodes. *J. Mol. Evol.* 91, 263-280.
  [doi:10.1007/s00239-022-10083-z](https://doi.org/10.1007/s00239-022-10083-z) — the source for
  the only homopolymer threshold measurement quoted here.
- McGee, R.S., Kinsler, G., Petrov, D. et al. (2024) Improving the accuracy of bulk fitness
  assays by correcting barcode processing biases. *Mol. Biol. Evol.* 41, msae152.
  [doi:10.1093/molbev/msae152](https://doi.org/10.1093/molbev/msae152)
- Wenger, A.M. et al. (2019) Accurate circular consensus long-read sequencing improves variant
  detection and assembly of a human genome. *Nat. Biotechnol.* 37, 1155-1162.
  [doi:10.1038/s41587-019-0217-9](https://doi.org/10.1038/s41587-019-0217-9)
- Illumina, *Index Adapters Pooling Guide*, document 1000000041074 v17 (July 2026), the
  *Sequencing Chemistry* and *Color Balance* pages; *2-Channel SBS Technology*; the Knowledge
  Base articles on nucleotide diversity (#1543), PhiX spike-in (#1527), low-diversity sequencing
  (#2882), and chemistry and imaging on MiniSeq (#5628), NovaSeq X (#7970) and MiSeq i100
  (#9348); and the DRAGEN *Read Trimming* documentation for the poly-G artefact.
- Illumina, US patent 11,028,436 B2, for the 25-75% GC and four-base homopolymer claim language.
- PacBio, *Multiplexing Targeted Sequencing using Barcodes*, PN 100-114-500-01 (2012), read from
  an Internet Archive copy because the live URL returns 404; *Guidelines for Using PacBio
  Barcodes for SMRT Sequencing*; the HiFi sequencing product pages and the Revio specification
  sheet.
- Oxford Nanopore, the accuracy and basecalling platform pages.
- `docs/research/protein-library-assembly.md` (issue #56) for the method, the licence verdict
  and the measurements of the published set this note builds on.

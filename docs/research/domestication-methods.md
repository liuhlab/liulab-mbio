---
search:
  exclude: true
---

# Domesticating a working vector: which route, and when

Research note for issue #245, under spec #225. #244 asks what may be changed in a lentiviral
vector; this asks **how**, and it has to pick one. Everything below was read on **2026-10-06**
unless a source carries its own date. Coordinates are the package's: 0-based, half-open.

## The recommendation

> **Replace the span.** Take the shortest stretch of the plasmid holding every site you must
> change, order that stretch as one synthetic clonal fragment with the sites already gone, and
> put it back into the uncut parent with one NEBuilder HiFi / Gibson assembly at two junctions
> in unique sequence.
>
> **Switch to Q5 site-directed mutagenesis when fewer than three sites have to change and each
> sits in unique sequence.** Below three, a round per site is cheaper and faster than an order.
>
> **Neither applies to a site inside a long repeat.** Measure the identical context around every
> site first. A site with more than about 100 bp of identical sequence on both sides cannot be
> reached by *any* oligo — mutagenic primer, assembly overlap or Type IIS overhang — because no
> oligo placed there is unique to one copy. Such a site is not domesticated at the bench. It is
> bought as part of a whole synthesised plasmid, or designed around by assembling with the other
> enzyme.

The third clause is the one that decides real jobs, and it is checked first because it can take
the whole task away. Section 3 is the measurement behind it.

The order matters: **count repeats, then count sites, then count the span.**

| What you measured | Route |
| --- | --- |
| Any site with >100 bp identical context both sides | Not a bench job. Change the enzyme, or buy the plasmid whole (section 3.4) |
| 1–2 sites, all unique context | Q5 SDM, one round each (section 4.1) |
| ≥3 sites, all unique context, inside one span a vendor will build as one fragment | Synthesise the span, one Gibson (section 4.5) |
| ≥3 sites spread wider than one synthetic fragment | Two synthetic fragments and one Gibson, or buy the whole plasmid |

Three candidates from the ticket are **not** recommended in any case: Golden Gate reassembly
(section 5, the self-reference trap), overlap-extension PCR (section 4.4), and multi-site
mutagenesis kits (section 4.1).

## 1. What was read

| What | Where | Read |
| --- | --- | --- |
| Q5 SDM kit, protocol, size claim, success figures | NEB *Instruction Manual* #E0554S, Version 2.0_7/20 | 2026-10-06 |
| QuikChange II, Lightning, Multi, Lightning Multi — efficiency and size limits | Agilent manuals #200523 Rev E1, #210518 Rev F1, #200513 Rev F1, #210514 Rev E1 | 2026-10-06 |
| Golden Gate fragment counts, internal-site advice, precloning rule | NEB *Golden Gate Assembly Kit (BsaI-HFv2)* manual #E1601 | 2026-10-06 |
| Gibson fragment count, overlap, size, the repeat answer | NEB *NEBuilder HiFi DNA Assembly* manual #E2621/#E5520, Version 6.0_1/26 | 2026-10-06 |
| Golden Gate fidelity at 13 and 35 fragments | Pryor et al. 2020, *PLoS ONE* 15(9): e0238592, PMC7467295 | 2026-10-06 |
| 24-fragment claim (abstract only, paywalled) | Potapov et al. 2018, *ACS Synth. Biol.* 7, 2665 | 2026-10-06 |
| How a published toolkit domesticates a part | Weber et al. 2011, *PLoS ONE* 6(2): e16765, PMC3041749 | 2026-10-06 |
| Q5 and Taq error rates | Potapov & Ong 2017, *PLoS ONE* 12(1): e0169774, PMC5218489 | 2026-10-06 |
| Whole-plasmid PCR on large plasmids, narrative | Hallak et al. 2017, *PLoS ONE* 12(5): e0177788; Zhang et al. 2021, *Sci. Rep.* 11, 10577, PMC8129136 | 2026-10-06 |
| Synthesis length, complexity and turnaround | Twist gene-synthesis product page and sequence-acceptance FAQ; GenScript gene-synthesis and GenBrick pages; IDT custom-gene and gBlocks pages; GENEWIZ gene-synthesis FAQ; VectorBuilder vector-cloning page | 2026-10-06 |
| Whole-plasmid and Sanger sequencing | Plasmidsaurus plasmid page; Eurofins Genomics whole-plasmid page; Azenta/GENEWIZ whole-plasmid and Sanger pages | 2026-10-06 |
| Repeat stability in *E. coli* | Thermo Stbl3 product page C737303; Al-Allaf et al. 2012, *3 Biotech* 2, 61–70, PMC3563744 | 2026-10-06 |
| PCR template switching | Odelberg et al. 1995, *Nucleic Acids Res.* 23, 2049, PMC306983 | 2026-10-06 |
| The eight-site case, and the LTR measurement's source sequence | `docs/research/working-vector-plvx-tetone.md`; `addgene-plasmid-171123-sequence-336792.fasta` | 2026-10-06 |

Downloads are under `reference_docs/synthesis_and_assembly/domestication-methods/`. Facts already
held by this repo are cited to `docs/research/golden-gate-assembly.md` and
`docs/research/gibson-assembly.md` rather than re-sourced.

**Three fetch problems worth recording**, because the next reader will hit them. `neb.com` HTML
pages return 403 to most fetchers; the manual PDFs under `neb.com/-/media/nebus/files/manuals/`
do not, and every NEB quote here comes from a PDF. GenScript's pages are JavaScript-rendered, so
its figures below are summaries rather than verbatim text and are marked as such. Plasmidsaurus
and Eurofins returned overlapping page text (section 6).

## 2. The job, stated generally

Domestication is a list of point changes, each removing one enzyme recognition site, with the
constraint that nothing else in the plasmid moves. Three counts decide how it is done:

- **n**, the number of sites.
- **r**, the number of sites sitting inside a repeat — a stretch the plasmid carries more than
  once, identical over more than an oligo's reach.
- **s**, the length of the shortest span holding all n sites.

The worked case is pLVX-TetOne-Puro-GFP (Addgene 171123), measured in
`docs/research/working-vector-plvx-tetone.md`: 9,895 bp, **n = 8** (6 BsaI, 2 BsmBI), **r = 2**,
**s = 4,996 bp** if the two LTR sites are excluded.

## 3. Repeats decide it, and they are measured not guessed

### 3.1 The measurement

The two LTRs of pLVX-TetOne-Puro-GFP are 634 bp each, at 0–634 and 6658–7292. Measured over the
Addgene-verified sequence on 2026-10-06, they **differ at exactly one offset, 23**. The BsaI
sites sit at the same offset, 454, in both — plasmid coordinates 454 and 7112.

Walking outwards from each site until the two copies differ:

| | bp |
| --- | --- |
| Identical context upstream of the site | **430** |
| The site | 6 |
| Identical context downstream of the site | **174** |
| **Total identical window centred on the site** | **610** |

Nothing in those 610 bp tells one LTR from the other. This was computed from
`addgene-plasmid-171123-sequence-336792.fasta` and is reproducible from it.

### 3.2 Why that ends every oligo-based route at those two sites

Every route in the ticket works by putting a short oligo at or beside the site: a mutagenic
primer for SDM, a 20–30 bp homology arm for Gibson, a 4-nt overhang plus primer for Golden Gate,
an overlapping chimeric primer for SOE. All four need that oligo to bind, or to anneal, in one
place.

- **Site-directed mutagenesis.** An outward-facing primer pair at 454 binds equally at 7112.
  The intended product is the full 9,895 bp plasmid; the primer pair also makes a 6,658 bp and a
  3,237 bp product, each a deletion of one LTR and everything between. The mixture transforms,
  and the deletions are *shorter*, so they amplify and clone preferentially.
- **Gibson / NEBuilder.** A junction at 454 and a junction at 7112 have, by the measurement
  above, **identical overlap sequence** at any overlap length NEB recommends — the longest NEB
  asks for is 20–30 bp and the identical window is 610. The assembly cannot tell the two
  junctions apart, so fragment order is not determined. NEB states the general rule itself: "one
  must ensure that each DNA fragment includes a unique overlap so that the sequences may anneal
  and are properly arranged… If having repetitive sequences at the ends of each fragment is
  unavoidable, the correct DNA assembly may be produced, albeit at lower efficiency than other,
  unintended assemblies" (NEBuilder HiFi E2621 FAQ 6, quoted in
  `docs/research/gibson-assembly.md` §10). Here the repeat is not *at* the end of an overlap; it
  *is* the overlap.
- **Golden Gate.** Worse, because the discriminator is only four bases. The two junctions would
  carry the same overhang. Pryor et al. 2020 define set fidelity as the product of each
  overhang's correct-ligation probability (`docs/research/golden-gate-assembly.md` §4.4); a set
  containing one overhang twice has no defined correct partner for either, and the design tools
  reject it.
- **Overlap-extension PCR.** Same binding problem as SDM, plus the chimeric primers.

So the ticket's framing — repeats are "a specific hazard for every PCR-based route" — is
confirmed and is stronger than a hazard. At 610 bp of identical context it is not a yield
penalty, it is a well-posedness failure: there is no oligo that specifies which copy you mean.

### 3.3 The "twice and identically" requirement is a consequence, not a separate problem

If a correction is made to one LTR and not the other, the two LTRs stop being copies. Reverse
transcription's strand transfer reads one LTR and writes the other, so the lab's own note
records this as a change to the mechanism rather than a cosmetic difference
(`docs/research/working-vector-plvx-tetone.md` §1). Any route that could target one LTR
specifically would therefore have to be run twice, with a second design, and the second run
targets a plasmid whose two LTRs now differ — so its primer *is* unique, which is the one
direction in which the problem gets easier. This is possible and nobody recommends it: two
serial rounds, each needing an anchor outside the 610 bp window, each re-cloning the whole LTR.

### 3.4 What is left for a site in a repeat

Two things, and only two.

1. **Change the enzyme instead of the vector.** This is almost always right and almost always
   forgotten. The two unreachable pLVX sites are BsaI. The vector's own assembly step uses
   BsmBI, and **both BsmBI sites are outside the repeat and domesticable**
   (`docs/research/working-vector-plvx-tetone.md` §1). Clearing BsmBI and recording the BsaI
   sites, rather than clearing both, costs nothing and removes the two hard sites from the job.
   Section 5 is why the two-enzyme split exists in the first place.
2. **Buy the whole plasmid.** The repeat moves from being your problem to being the vendor's —
   and section 4.5 is whether they accept it.

And the repeat charges you on that route too. Twist's published complexity rule flags **"Direct
repeats >200 bp"** as High complexity (*Gene synthesis FAQ, sequence acceptance criteria*, read
2026-10-06). A 634 bp perfect repeat is three times that threshold. The consequence is not a
refusal — Twist's "only hard rules" are homopolymers ≥14 bp and CcdB — but the High-complexity
price tier, 15 business days instead of 10, and no Express option. Separately, Twist's clonal
genes top out at **7.0 kb**, so a 9,895 bp plasmid is off its catalogue entirely.

**This is the point of the whole note.** The 634 bp repeat is the one feature that penalises
every route at once: it blocks the oligo-based routes outright, and it pushes the buy-it-whole
route off the commodity tier onto quote-only long-construct services (section 4.5, section 7).
The route that dodges it is the one that never includes the repeat in anything you order or
amplify.

### 3.5 The general test

This generalises past lentivirus. Before choosing a route, for every site, walk outwards until
the sequence stops being unique in the plasmid. If that window is longer than the oligo the
route needs — roughly 25 bp for Gibson, 45 bp for a mutagenic primer, 4 bp plus a primer for
Golden Gate — the route cannot address that site. Common offenders besides retroviral LTRs:
AAV ITRs, tandem promoter or terminator copies, repeated selection cassettes, and two copies of
the same tag.

## 4. The routes, measured

### 4.1 Site-directed mutagenesis

**One site at a time, with NEB's kit, is the only mutagenesis worth running on a 10 kb plasmid.**

NEB Q5 Site-Directed Mutagenesis Kit, *Instruction Manual* #E0554S, Version 2.0_7/20. The
mutagenic primers face outwards from the site, back to back; the PCR copies the whole plasmid;
a single 5-minute room-temperature KLD step kinases, ligates and DpnI-digests, and 5 µl of it
goes straight into cells. Template is "1–25 ng of plasmid DNA". The manual's size claim, stated
once in the Introduction: the kit "ensures robust results with plasmids up to at least 14 kb in
length".

Its success figures, verbatim, and what they cover:

> For substitutions and deletions, thousands of colonies were observed with > 90% of the
> colonies having incorporated the desired mutation.

And the caption of its Figure 4:

> In all three cases, over 90% of the resultant colonies contained the desired mutation(s).

Those three cases are **one change each**, on 5.8, 6.7 and 7.0 kb plasmids. **No multi-site
correction is measured in that manual**, and nothing published that was found here measures
success against plasmid size for any kit. The kits cannot be ranked on measured data; they are
compared below on what their own manuals permit.

**The multi-site kits are out of spec for a plasmid this size, by their own manuals.** Agilent
QuikChange Multi (#200513, Rev F1): "virtually any plasmid of up to 8 kb is a suitable
template"; it "introduces mutations at three different sites simultaneously in the 4-kb
QuikChange Multi control plasmid with greater than 50% efficiency", with at most 5 mutagenic
primers. QuikChange Lightning Multi (#210514, Rev E1) gives "greater than 55% efficiency" on the
same 4 kb, 3-site control and the same 8 kb cap. The single-site kits claim more — QuikChange II
(#200523, Rev E1) "greater than 80% efficiency in a single reaction", QuikChange Lightning
(#210518, Rev F1) "greater than 85%" with plasmids of 4–14 kb tested — but they change one site.
Agilent also notes that its SURE 2 cells are "not recommended for use with mutagenized plasmids
greater than 10 kb".

So for eight sites on a 9,895 bp plasmid: the multi kits are over their stated size limit and
over their stated primer count, and the single-site kits mean **eight serial rounds**.

**Two papers say the size problem is real, and neither is a benchmark.** Hallak et al. 2017,
*PLoS ONE* 12(5): e0177788, working on plasmids from 2.6 to 17 kb, states that whole-plasmid
copying methods "generally work best for plasmids under 3.1 kb" and that "Introducing a mutation
in plasmids larger than 8 kb usually requires subcloning". Zhang et al. 2021, *Sci. Rep.* 11,
10577 (doi:10.1038/s41598-021-89884-z), targeting a 13.3 kb lentiviral plasmid, reports that
three conventional whole-plasmid methods — including one run with Q5 — gave nothing: "the other
three conventional methods could not give any positive plasmid… we couldn't gain any mutant
through these methods", while their own two-PCR method reached "the rate of the plasmids with a
correct size reached 93%". **Both are narrative, on one target, with the authors' own
implementations as the controls.** Neither is a rate for the commercial kits, and neither is
cited here as one.

**No mutagenesis manual read here mentions repeats at all.** All five Agilent and NEB manuals
were searched for "repeat", "homologous" and "recombin"; the only hit is a generic troubleshooting
entry on "unwanted deletion or recombination of plasmid DNA", not tied to repeated sequence and
not tied to whole-plasmid PCR across one. The hazard in section 3 is therefore **undocumented by
the vendors**, not merely unquantified.

### 4.2 PCR the plasmid in pieces, reassemble by Golden Gate

**Ruled out for a vector being cleared of a Type IIS enzyme, and the reason is section 5.** Two
further facts rule it out even where the self-reference problem does not apply.

The method itself is well documented and works at far higher fragment counts than eight. Pryor
et al. 2020, *PLoS ONE* 15(9): e0238592 (PMC7467295), reports a 13-fragment SapI assembly where
"on average 91% of the observed transformants were blue, indicating uptake of a correct assembly
product" against a predicted fidelity of 79%, and a 35-fragment BsmBI-v2 assembly where "on
average 71% of the observed transformants harbored accurately assembled constructs, compared to
a theoretical prediction of 65%", with ">700 transformants harboring correct assembly products
per μL of the assembly reaction". Potapov et al. 2018, *ACS Synth. Biol.* 7, 2665–2674
(doi:10.1021/acssynbio.8b00333), states in its abstract that the ligation profile "enabled
accurate and efficient assembly of a lac cassette from up to 24-fragments in a single reaction"
— the paper is paywalled and its per-fragment percentages were not read. NEB's own framing
(E1601 manual, Q8): Golden Gate "is normally used for insert assemblies of 5–10 or more
fragments".

So fragment count is not the objection. These are the objections:

- **Identical overhangs at the repeat.** Section 3.2.
- **NEB tells you not to PCR a repeat.** E1601 manual: precloning "is superior for inserts
  < 250 bp or > 3 kb, or those containing repetitive elements that might accumulate errors
  during PCR amplification". A lentiviral backbone is exactly that.
- **Two of the fragments would be 152 bp and 94 bp** (section 8), and NEB's precloning advice
  covers inserts under 250 bp as well.

### 4.3 PCR the plasmid in pieces, reassemble by Gibson

Better than Golden Gate here — no Type IIS enzyme is involved, so the self-reference trap does
not apply, and no scar is left. It is still not the recommendation, because it needlessly PCRs
the whole plasmid.

NEBuilder HiFi (#E2621 manual, Version 6.0_1/26) on fragment count, Q5:

> has been used to efficiently assemble up to eleven, 0.4 kb inserts into a vector at one time.
> However, we recommend the assembly of five or fewer inserts into a vector in one reaction, in
> order to produce a clone with the correct insert. A strategy involving sequential assembly can
> be used if all of the fragments cannot be assembled in a single reaction.

The FAQ page adds: "Alternatively, consider Golden Gate Assembly for assemblies of >6 fragments."
Eight corrections at one junction each is **eight inserts**, above NEB's recommended five.
Overlap is "15–30 nt… with a Tm equal to or greater than 48°C", and "Longer overlaps (20-30 base
pairs) increase assembly efficiency". Size is not a problem: Q4 records the mix cloning "a 12 kb
DNA fragment into a 7.4 kb plasmid… totaling up to 19 kb".

And NEB's repeat answer, Q6, verbatim — read the middle sentence, because it is the
recommendation in NEB's own words:

> Is this method applicable to the assembly of repetitive sequences? A6: Yes. However, one must
> ensure that each DNA fragment includes a unique overlap so that the sequences may anneal and
> are properly arranged. **The repetitive sequence can also be internalized in the first stage of
> a two-stage assembly strategy.** If having repetitive sequences at the ends of each fragment is
> unavoidable, the correct DNA assembly may be produced, albeit at lower efficiency than other,
> unintended assemblies.

"Internalized" is the whole idea of section 4.5: keep the repeat **inside** an uncut piece, so it
is never at a junction and never in an oligo. The span route does this by leaving the repeat in
the parent backbone, which is not even cut.

**PCR exposure is the other cost.** Q5 is NEB's recommendation for assembly fragments, and its
error rate is published: 5.3 × 10⁻⁷ (± 0.9 × 10⁻⁷) substitutions per base per doubling, against
Taq at 1.5 × 10⁻⁴ (Potapov & Ong 2017, *PLoS ONE* 12(1): e0169774, PMC5218489; NEB's own pages
give the ratio as "~280 times higher than Taq"). The same paper notes that for polymerases this
accurate, "DNA damage from cycling, not polymerase error, is a major contributor to mutations" —
so the exposure that matters is **how many bases you copy, how many times**, not the rate alone.
Reassembling by Gibson copies all 9,895 bp. The span route copies about 7.1 kb once, and the
2,748 bp of the span not at all.

### 4.4 Overlap-extension PCR

Not recommended, and the reason is arithmetic rather than reputation. SOE joins fragments by
chimeric primers and a second PCR, then one ligation. For n corrections it needs n+1 pieces and
n chimeric primer pairs, and the whole product is amplified again after joining — so every base
of the 10 kb passes through polymerase twice rather than once. It has the same repeat problem as
every other route in section 3, plus a second copy of it in the joining PCR. **No published
success rate for SOE on a plasmid of this size was found**, so it is not ranked against the
others on evidence; it is set aside because it is strictly more PCR for the same result as a
Gibson, and Gibson is the better-documented reaction.

### 4.5 Order the corrected span, assemble it back

This is the recommendation, and it is a hybrid rather than any single row of the ticket's table.
Do not re-synthesise the plasmid and do not PCR it into n+1 pieces. Instead:

1. Take the shortest span holding every site you must change.
2. Order that span as one synthetic clonal fragment, with every site already removed.
3. Linearise the parent by PCR or by a single unique-cutter digest at two points just outside the
   span, in sequence that is unique in the plasmid.
4. Join them with NEBuilder HiFi / Gibson. Two junctions, two overlaps, one reaction.

What it buys: **one assembly, one transformation, one sequencing run, and zero mutagenic
primers**. The corrections are not made at the bench at all — they arrive already made, and the
only bench risk is the two junctions and the one PCR of the backbone. Every one of the n
corrections is verified by the vendor before it ships.

Gibson is the assembler, not Golden Gate, for the reason in section 5.

## 5. The self-reference trap

**A vector being cleared of BsaI cannot be assembled by BsaI**, because the enzyme would cut the
sites you have not removed yet, and because the sites you have removed are the junctions. The
rule is obvious once stated and it is worth stating, because it quietly removes a whole row of
the ticket's table.

**No source read here states it in those words.** What exists is NEB's answer to the question it
is the answer to, from the Golden Gate Assembly Kit (BsaI-HFv2) manual #E1601, Q3 and Q4:

> internal sites need to be eliminated by site-directed mutagenesis

and, asked the same thing again:

> Either use site-directed mutagenesis to eliminate the internal sites, screen your sequences for
> the absence of other Type IIS restriction sites… such as BbsI, SapI/BspQI or BtgZI…, or
> consider another assembly approach such as NEBuilder HiFi DNA Assembly if the assembly will
> involve 5 or less inserts.

Three options, and NEB lists them in that order: fix the sites, **use a different enzyme**, or
**stop using Type IIS**. The second and third are the live ones for domestication work, and the
third is why section 4.5 assembles with Gibson.

**The published toolkits take the second.** Weber et al. 2011, the MoClo paper (*PLoS ONE* 6(2):
e16765, PMC3041749), domesticates a part by PCR with mismatch primers and then assembles it with
a *different* enzyme — Figure 2B: "Removal of a BsaI site in a fragment of interest is done by
amplifying two fragments… Primers pr2 and pr3 span the BsaI recognition site and introduce a
single nucleotide mismatch… As all primers have BpiI recognition sites in their 5′ extensions,
the PCR fragments are cloned with a BpiI-based Golden Gate cloning reaction." The paper reports
doing this across 11 ORFs — "we cloned a number of level 0 modules, removing at the same time all
internal interfering type IIS recognition sites" — and **gives no success rate**. Pollak et al.
2019 (Loop assembly, *New Phytol.*, doi:10.1111/nph.15625) defines domestication the same way,
as removal "by the introduction of synonymous mutations", and likewise publishes no domestication
success rate. Lee et al. 2015 (Yeast Toolkit) and Martella et al. 2017 (EMMA) are paywalled and
their domestication methods were not read.

So the precedent is real but it is a *recipe*, not a *rate*: **the established toolkits
domesticate by mismatch-primer PCR plus an orthogonal-enzyme assembly, and none of them reports
how often it works.** That is why this note does not rank the routes on published success rates
for domestication — there are none to rank.

Two practical consequences:

- Clearing BsaI, assemble with BsmBI or BbsI, and vice versa. NEB E1601 Q3 confirms the pair is
  interchangeable when the sequence is clean of both: "If your sequences have neither BsaI nor
  BsmBI sites, either can be used as both support even 24 fragment assemblies."
- Clearing **both** BsaI and BsmBI, as a vector for a two-enzyme scheme must be, neither is
  available. You are down to a third enzyme or to Gibson. **This is the common case in this lab's
  work and it is the simplest argument for Gibson as the default assembler in section 4.5.**

## 6. What each route owes in verification

Removing a site is easy. Proving nothing else moved is the work, and the routes owe different
amounts of it.

**Whole-plasmid nanopore sequencing has made this cheap enough that there is no excuse.**
Plasmidsaurus (read 2026-10-06) states a 1-day target turnaround, plasmids to 300 kb, and
"Consensus accuracy is often above Q60, which corresponds to 99.9999%, or one error per
1,000,000 bases", at a **list price of $15 per sample under 25 kb**. Eurofins Genomics quotes
the same tier structure and states "Raw read accuracy is 98.3% and the consensus accuracy for
SNPs is 99.6% with 50X coverage", 1 business day, 2.5–300 kb. Azenta/GENEWIZ offers next-morning
whole-plasmid sequencing under 25 kb. **Caution: the two page fetches that produced the
Plasmidsaurus and Eurofins numbers returned overlapping wording, so which figure belongs to
which service should be confirmed on the page before either is printed in a protocol.**

**Sanger does not cover a 10 kb plasmid cheaply.** GENEWIZ states Sanger "read lengths up to
~1000 bases. A typical read will provide 800 bases Phred20" (read 2026-10-06). At 800 bases that
is **13 reactions for one strand of 9,895 bp**, before any primer walking that fails. Against a
$15 list whole-plasmid read, Sanger is only worth running on a single short span.

**A diagnostic digest checks size and arrangement, not bases.** Addgene's own protocol page
describes the fragment pattern as something that "can indicate if the plasmid contains the
expected size insert", covering total size, insert and backbone size, orientation and fingerprint
identity, and points onward to "more expensive forms of plasmid verification, such as DNA
sequencing". **No source read here states outright that a digest misses point mutations and
small indels** — that follows from what a digest measures, but it is an inference, not a quote.

What each route owes:

| Route | What can go wrong that is invisible without sequencing | Minimum verification |
| --- | --- | --- |
| Q5 SDM, n rounds | A polymerase error anywhere in the 9,895 bp copied, in **every** round | Whole-plasmid read **after each round**, not only at the end |
| PCR in pieces + assembly | Polymerase error in any piece; junction scars; a repeat-driven rearrangement | One whole-plasmid read, plus a digest to catch a rearrangement the assembler mis-ordered |
| SOE | All of the above, with the product amplified twice | As above |
| Synthetic span + one Gibson | Error in the two junctions and in the one PCR of the backbone. The span itself is vendor-verified | One whole-plasmid read |

The asymmetry is the argument. **n rounds of whole-plasmid PCR means n chances to put a silent
error anywhere in 10 kb, and n whole-plasmid reads to catch them.** A synthetic span means one.

**A lentiviral plasmid also owes a stability check that has nothing to do with the route.**
Thermo describes Stbl3 as "designed for cloning direct repeats found in lentiviral expression
vectors as these cells reduce the frequency of homologous recombination of direct repeats found
in lentiviral vectors", by the recA13 allele (product C737303, read 2026-10-06). NEB's NEBuilder
manual (Version 6.0_1/26) puts it generally: "Some DNA structures, including inverted and tandem
repeats, are selected against by E. coli." Al-Allaf et al. 2012, *3 Biotech* 2, 61–70
(PMC3563744), measured an HIV-based lentiviral plasmid being structurally changed or lost in
500 mL Stbl2 cultures but not in 5 mL, and reports for Stbl3 that "the percentage of
plasmid-carrying Stbl3 cells was estimated to be 100 % throughout the entire 12 h interval".
**Grow every intermediate and every final clone in a recA-deficient strain, in a small culture,
whatever route produced it.**

## 7. Cost and elapsed time, for the eight-site case

Price is an input here, not a fact — a real price is negotiated and not public (#241). Every
figure below is a **list price read 2026-10-06** and is given only to show the structure. What
travels between labs is the **quantities**: how many oligos, how many reactions, how many
sequencing runs, how many serial waits.

### 7.1 Cost structure, not totals

| Route | What you buy | What scales with | What scales with plasmid size |
| --- | --- | --- | --- |
| Q5 SDM | 2 oligos, 1 PCR, 1 KLD, 1 transformation, 1 miniprep, 1 whole-plasmid read — **per site** | **n**, linearly, and serially | extension time; error probability per round |
| PCR + assembly | n+1 primer pairs, n+1 PCRs, 1 assembly, 1 read | **n** | fragment count and PCR length |
| SOE | n+1 PCRs plus n joining PCRs | **n** | twice over |
| Synthetic span | 1 synthetic fragment, 1 backbone PCR, 1 assembly, 1 read | **span length**, not n | not at all — the parent is never copied whole |
| Whole clean plasmid | 1 synthetic plasmid | nothing | **plasmid length, and repeat content** |

The fourth row is why the recommendation is what it is: **its cost does not grow with the number
of sites.** Eight corrections inside one 2.7 kb span cost the same as one.

### 7.2 The list prices, labelled

Twist Bioscience clonal genes, product page read 2026-10-06. **List price, per construct, plus a
vector fee of $35 (cloning vector), $55 (expression) or $75 (custom).**

| Length | Low/Medium complexity | High complexity |
| --- | --- | --- |
| 0.3–1.8 kb | 9¢/bp | 12¢/bp |
| 1.8–3.2 kb | 13¢/bp | 18¢/bp |
| 3.2–7.0 kb | 18¢/bp | 23¢/bp |

Turnaround, same page: Express 4–7 business days, Standard 10, **High complexity 15**, and
Express is not offered for High complexity. Clonal range 0.3–7.0 kb. (Twist's FAQ page says
300–5,000 bp, which contradicts the product page; the discrepancy is unresolved.)

Beyond 7 kb the commodity tier ends. GenScript's GenBrick covers 8–200 kb at a list of **$0.39/bp
for 8–15 kb**, 10–30 business days by size. IDT routes long and difficult constructs to its
Ansa-manufactured tier, 7,500–50,000 bp, "as soon as 25 business days", which claims support for
"difficult sequences, including extreme GC content, repeats, secondary structure, promoters, ITRs
and poly(A) tails" but states "Custom Vectors not supported at this time". Azenta/GENEWIZ claims
construction "up to 150 kb" with a ">99.8% project success rate, including high GC content,
repeats, and complex elements", noting such sequences "may require the FLEX service, adjusted
timelines and charges" and publishing **no price at all** — "General pricing requires contacting
your sales representative". VectorBuilder builds lentiviral transfer vectors to order and states
"De novo synthesis is required in less than 20% of our projects", with pricing that "correlates
with the fragment length, sequence complexity (e.g., high GC content, simple repeats, or
segmental repeats)" and is quote-only.

**No vendor read here states a refusal rule for a long perfect direct repeat, and none mentions
retroviral LTRs at all.** Twist's ">200 bp direct repeat" flag is the only number anyone
publishes. Everything else is "contact us".

### 7.3 Elapsed time, as the quantities imply

| Route, 8 sites | Serial waits | Where the time goes |
| --- | --- | --- |
| Q5 SDM, 1 site/round | **8** | Per round: under 2 h bench (NEB E0554 manual, for its 6.7 kb control), then overnight plating, pick, overnight culture, miniprep, sequencing return. A round is not finished until its read comes back |
| PCR + Gibson | 1 | n+1 PCRs in parallel, 1 assembly, 1 read |
| Synthetic span + Gibson | 1 order + 1 | 10 business days (Twist Standard, low/medium) in parallel with nothing, then one assembly week |
| Whole clean plasmid | 1 order | 15+ business days at Twist High complexity if it fit, which at 9.9 kb it does not; 10–30 at GenBrick; quote-only elsewhere |

The eight-round column is the one people underestimate. **Eight serial rounds, each gated on a
sequencing result, is the slowest route on this table by a wide margin**, and it is the one most
often chosen because each round looks cheap.

## 8. The worked case, end to end

pLVX-TetOne-Puro-GFP, 9,895 bp, eight sites. Applying the rule in order.

**Step 1, count repeats.** Two sites — BsaI at 454 and 7112 — sit inside 610 bp of identical
sequence (section 3.1). **They come off the list.** Per section 3.4 the fix is to change the
enzyme, not the vector: the vector's own assembly step uses BsmBI, both BsmBI sites are outside
the repeat, and the lab's requirement should be narrowed to "no BsmBI outside the cassette, BsaI
counted and recorded" — a decision the working-vector note already flags as open.

**Step 2, count the remaining sites and their span.** Six left — BsaI at 3724, 5730, 6472 and
8720, BsmBI at 3876 and 5636.

Measured on the sequence, the spacing is uneven, and that itself rules out the one-junction-per-
site layouts:

| From | To | bp |
| --- | --- | --- |
| 454 | 3724 | 3270 |
| 3724 | 3876 | **152** |
| 3876 | 5636 | 1760 |
| 5636 | 5730 | **94** |
| 5730 | 6472 | 742 |
| 6472 | 7112 | 640 |
| 7112 | 8720 | 1608 |
| 8720 | 454 | 1629 |

A Golden Gate or Gibson reassembly with one junction per site makes **a 152 bp fragment and a
94 bp fragment**. With NEB's recommended 20–30 bp overlaps on both ends, the 94 bp fragment is
more than half overlap. This is a second, independent reason not to reassemble this plasmid from
PCR pieces.

**Step 3, pick the span.** Five of the six lie between 3724 and 6472 — **2,748 bp, repeat-free,
and squarely inside Twist's 1.8–3.2 kb Low/Medium tier**. Both ends are in unique sequence: 3724
is in the hPGK promoter, 6472 is in the 206 bp between the WPRE and the 3' LTR that no annotation
claims. Order those 2,748 bp with all five sites removed; NEBuilder it into the parent
backbone, PCR-linearised outward from 3724 and 6472.

**Step 4, mop up.** The sixth site, 8720, is in AmpR, outside the LTRs, in unique sequence, and a
synonymous codon change clears it. **One Q5 SDM round**, by the ≤2-sites clause of the rule.

**Why not extend the span to cover 8720 too?** Because 3724→8720 is 4,996 bp and crosses the
3' LTR, so the order would carry a 634 bp perfect repeat. That moves it from the 1.8–3.2 kb
Low/Medium tier to the 3.2–7.0 kb **High complexity** tier and from 10 to 15 business days, to
save one SDM round. On list prices that is roughly $357 against roughly $1,149 — and the longer
order would also put a second, differing copy of the LTR into the plasmid, which section 3.3 says
not to do. **Take the shortest repeat-free span, not the span that covers the most sites.**

**Net:** one synthetic fragment, one assembly, one SDM round, two whole-plasmid reads, two sites
deliberately left alone. Against eight serial SDM rounds and eight reads.

## 9. Open gaps

- **Nobody has quantified PCR across a long perfect direct repeat.** Odelberg et al. 1995,
  *Nucleic Acids Res.* 23, 2049 (PMC306983) establishes the mechanism — Taq generates recombinants
  within a single extension by switching template, and "Recombination is reduced several fold
  when the complementary template strands are physically separated" — but **no source found here
  measures mispriming, template switching or deletion frequency for a template carrying two
  identical ~600 bp repeats.** Section 3's conclusion rests on the uniqueness argument, which is
  arithmetic, not on a measured failure rate.
- **No mutagenesis kit vendor documents the repeat hazard at all** (section 4.1). The absence is
  itself the finding.
- **No systematic benchmark of SDM success rate against plasmid size exists.** The kits are
  compared here on their own stated limits only, never ranked on reputation.
- **No published success rate for overlap-extension PCR** on a plasmid of this size was found.
- **Whether a vendor will build a 9.9 kb lentiviral plasmid with intact 634 bp LTRs, and at what
  price and turnaround, is quote-only** at every vendor read. It was not quoted. If the lab ever
  wants the buy-it-whole number properly, that is one email to GenScript, Azenta and
  VectorBuilder with the actual sequence, and it would settle the one route this note prices
  worst.
- **Which whole-plasmid-sequencing figures belong to Plasmidsaurus and which to Eurofins** needs
  confirming on the pages (section 6).
- **Addgene does not state that a digest misses point mutations.** The table in section 6 marks
  that as inference.
- The two LTR BsaI sites are recorded as unreachable **by the bench routes**. Whether they are
  safe to change at all is a separate open question, already held in
  `docs/research/working-vector-plvx-tetone.md` §1: no source read there measures what a base
  change at TAR positions 1–6 costs.

## Sources

Licences follow the rule `docs/research/restriction-enzyme-data.md` §2 set: NEB, Agilent, Thermo
and vendor pages are **cited, never mirrored**, and no number here is a build input. The
open-access papers are CC BY and may be quoted with attribution.

### Manuals and kit inserts

- NEB, *Q5 Site-Directed Mutagenesis Kit Instruction Manual*, #E0554S, Version 2.0_7/20,
  <https://www.neb.com/-/media/nebus/files/manuals/manuale0554.pdf>.
- NEB, *NEBuilder HiFi DNA Assembly* manual, #E2621/#E5520, Version 6.0_1/26,
  <https://www.neb.com/en/-/media/nebus/files/manuals/manuale2621_e5520.pdf>.
- NEB, *Golden Gate Assembly Kit (BsaI-HFv2)* manual, #E1601,
  <https://www.neb.com/-/media/nebus/files/manuals/manuale1601.pdf>.
- Agilent, *QuikChange II* #200523 Rev E1; *QuikChange Lightning* #210518 Rev F1; *QuikChange
  Multi* #200513 Rev F1; *QuikChange Lightning Multi* #210514 Rev E1.
- Takara, *PrimeSTAR GXL DNA Polymerase* manual R050A.

### Papers

- Pryor, J. M. et al. Enabling one-pot Golden Gate assemblies of unprecedented complexity using
  data-optimized assembly design. *PLoS ONE* **15**(9): e0238592 (2020). PMC7467295. CC BY 4.0.
- Potapov, V. et al. Comprehensive profiling of four base overhang ligation fidelity by T4 DNA
  ligase and application to DNA assembly. *ACS Synth. Biol.* **7**, 2665–2674 (2018).
  [doi:10.1021/acssynbio.8b00333](https://doi.org/10.1021/acssynbio.8b00333). **Abstract only** —
  paywalled, no PMC copy.
- Potapov, V. & Ong, J. L. Examining sources of error in PCR by single-molecule sequencing.
  *PLoS ONE* **12**(1): e0169774 (2017). PMC5218489. CC BY 4.0.
- Weber, E. et al. A modular cloning system for standardized assembly of multigene constructs.
  *PLoS ONE* **6**(2): e16765 (2011). PMC3041749. CC BY.
- Pollak, B. et al. Loop assembly: a simple and open system for recursive fabrication of DNA
  constructs. *New Phytol.* (2019).
  [doi:10.1111/nph.15625](https://doi.org/10.1111/nph.15625). Read via the bioRxiv preprint text.
- Hallak, L. K. et al. A modified universal restriction-free method for rapid and efficient
  site-directed mutagenesis. *PLoS ONE* **12**(5): e0177788 (2017).
- Zhang, L. et al. *Sci. Rep.* **11**, 10577 (2021).
  [doi:10.1038/s41598-021-89884-z](https://doi.org/10.1038/s41598-021-89884-z). PMC8129136.
- Al-Allaf, F. A. et al. Remarkable stability of an instability-prone lentiviral vector plasmid in
  *Escherichia coli* Stbl3. *3 Biotech* **2**, 61–70 (2012). PMC3563744.
- Odelberg, S. J. et al. Template-switching during DNA synthesis by *Thermus aquaticus* DNA
  polymerase I. *Nucleic Acids Res.* **23**, 2049–2057 (1995). PMC306983.
- Horton, R. M. et al. Engineering hybrid genes without the use of restriction enzymes: gene
  splicing by overlap extension. *Gene* **77**, 61–68 (1989).
  [doi:10.1016/0378-1119(89)90359-4](https://doi.org/10.1016/0378-1119(89)90359-4). **Cited for
  the method only** — the full text was not read, and no success rate was taken from it.

### Vendor pages

All read 2026-10-06. Every price is a list price.

- Twist Bioscience, *Gene Synthesis* product page and the gene-synthesis sequence-acceptance FAQ.
- GenScript, *Gene Synthesis* and *GenBrick* pages. **JavaScript-rendered — summarised, not
  verbatim.**
- IDT, *Custom Gene Synthesis*, *gBlocks* FAQ, and the gene-synthesis turnaround FAQ.
- Azenta/GENEWIZ, *Gene Synthesis FAQs*, *Whole Plasmid Sequencing*, *Sanger Sequencing*.
- VectorBuilder, *Vector Cloning* service page.
- Plasmidsaurus, *Plasmid* page. Eurofins Genomics, *Whole Plasmid Sequencing*.
- Thermo Fisher, *One Shot Stbl3 Chemically Competent E. coli*, C737303.
- Addgene, *Diagnostic digest* protocol page. Addgene's terms are non-commercial; the page is
  cited, and nothing from it is shipped.

### Within this repo

- `docs/research/working-vector-plvx-tetone.md` — the eight sites, the LTR annotation, the TAR
  question. No number in it is recomputed here except the LTR identity window, which is new.
- `docs/research/golden-gate-assembly.md` §3, §4.4 — PCR-linearised vectors, Pryor's fidelity
  definition.
- `docs/research/gibson-assembly.md` §10 — the NEB and Takara repeat positions.
- `docs/research/restriction-ligation.md` §1 — the licence rule this note follows.

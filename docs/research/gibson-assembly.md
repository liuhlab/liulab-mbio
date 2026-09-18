---
search:
  exclude: true
---

# Gibson assembly, its assembly products and the two oligo routes

Research note for issue #119, the first sub-issue of spec #116 under #74. Everything below was
retrieved on **2026-09-18** unless a different date is given next to the source. It answers the
fourteen questions the ticket asks, and it is where every number in `liulab_mbio.cloning.gibson`
must come from.

## How to read this note

`neb.com` returns HTTP 403 to `curl` and to WebFetch for its HTML pages, and
`docs/research/golden-gate-assembly.md` records three routes that work. Today two of those three
worked and one did not, and a fourth trap cost the previous note an open gap:

| Route | Worked today | Used for |
| --- | --- | --- |
| `neb.com/-/media/nebus/files/manuals/*.pdf` | yes (plain `curl`) | E2621/E5520 and E2611/E5510 manuals |
| `neb.com/-/media/catalog/specifications/<a>/<b>/*.pdf` | yes (plain `curl`) | E2611 and E2623 product specifications |
| `neb.com/-/media/nebus/files/application-notes/*.pdf` | yes (plain `curl`) | the bridging-oligo application note |
| `web.archive.org` | **no** | nothing — the Internet Archive answered "Internet Archive services are temporarily offline" all day |
| `protocols.io` PDF export | yes | NEB's own E2611 and E2621 reaction protocols, under CC BY |
| `takarabio.com` | yes (plain `curl`) | every In-Fusion manual and page |
| `neb-online.de` mirror | yes | the NEBuilder HiFi brochure and product page |

**The trap that produced the previous note's open gap is the filename, not the host.** The
NEBuilder HiFi manual is `manuale2621_e5520.pdf`, not `manuale2621.pdf`. A wrong filename and a
blocked request return the same thing — an identically sized HTML page with status 403 — so a
missing file reads as a block. Search for the exact filename before concluding a manual is
unreachable, and compare response sizes: two 403s of exactly equal length are the generic page. The
same
applies to `manuale2611_e5510.pdf`, which is revision 3.0_1/26 — one revision newer than the
`manuale2611.pdf` the same directory still serves. Both were downloaded, and both are quoted below
with their revision.

Everything from the Internet Archive is therefore absent from this note; where a claim would have
needed it, the claim is in [Open gaps](#19-open-gaps) instead. Nothing here is guessed. Where a
sentence states something no source states, it says that it is inference.

The three assembly products the spec names are abbreviated below as **NEBuilder HiFi**
(NEB #E2621/#E5520), the **Gibson Assembly Master Mix** (NEB #E2611/#E5510) and **In-Fusion**
(Takara In-Fusion Snap Assembly, which replaced In-Fusion HD).

## 1. Sources and licences — what may ship

This is the ticket's question 14, answered once for everything below.

| Dataset or document | Licence | Can we ship it? |
| --- | --- | --- |
| NEB instruction manuals E2621/E5520 v6.0_1/26 and E2611/E5510 v3.0_1/26 | © NEB, all rights reserved | **No** — cite and paraphrase the numbers, never redistribute the file |
| NEB product specifications PS-E2611S/L v2.0, PS-E2623S v1.0 | © NEB | **No** — same |
| NEB application note *Bridging dsDNA with a ssDNA Oligo…*, 06/17 | © 2018 NEB, all rights reserved | **No** — same |
| NEB *NEBuilder HiFi DNA Assembly* brochure | © 2024 NEB, all rights reserved | **No** — same |
| NEB's own protocols on protocols.io (E2611 v2, E2621 v2) | **CC BY** | **Yes**, with attribution — see below |
| Takara In-Fusion manuals, protocols-at-a-glance and web pages | © Takara Bio Inc., all rights reserved; site Terms of Use forbid republishing | **No** — cite and paraphrase, short quotations only, never mirror the document |
| Gibson, D.G. et al. (2009) *Nat. Methods* 6, 343–345 | © 2009 Nature America, closed access | **No** — cite and paraphrase; the copy read is a third-party mirror |
| Gibson, D.G. (2011) *Methods Enzymol.* 498, 349–361, PMC7149801 | © Elsevier, free to read and re-use with acknowledgement under a revocable grant | **Quote with attribution**; do not ship as data |
| De Saeger, J. et al. (2022) bioRxiv 2022.02.10.479870 | **CC BY-NC 4.0** | **No** for an MIT package — non-commercial clashes, as it does for Potapov 2018; quote with attribution |
| Addgene's Gibson assembly protocol page | Terms of Use: informational, non-commercial use, attribution required | **No** — paraphrase with citation; do not ship the text |
| Thermo Fisher GeneArt user guide MAN0019062 Rev. D and the oligo-stitching white paper | © Thermo Fisher Scientific, all rights reserved | **No** — cite and paraphrase |

Licence text, quoted:

- protocols.io, on both NEB protocol pages: "License: This is an open access protocol distributed
  under the terms of the Creative Commons Attribution License, which permits unrestricted use,
  distribution, and reproduction in any medium, provided the original author and source are
  credited." The reaction table, the optimal-quantity paragraph and the 650-dalton formula all
  appear there, so **the NEBuilder and Gibson reaction tables are quotable under CC BY** even
  though the manual that also carries them is not. Attribution: "New England Biolabs (2022)
  NEBuilder HiFi DNA Assembly Reaction (E2621), protocols.io, doi:10.17504/protocols.io.bfhrjj56".
- NEB manuals, last page: "Products and content are covered by one or more patents, trademarks
  and/or copyrights owned or controlled by New England Biolabs, Inc (NEB)." and
  "© Copyright 2026, New England Biolabs, Inc.; all rights reserved."
- Gibson 2011, the PubMed Central copy: "Copyright © 2011 Elsevier Inc. All rights reserved."
  followed by Elsevier's COVID-19 resource-centre grant — "Elsevier hereby grants permission to
  make all its COVID-19-related research that is available on the COVID-19 resource centre —
  including this research content — immediately available in PubMed Central and other publicly
  funded repositories… with rights for unrestricted research re-use and analyses in any form or by
  any means with acknowledgement of the original source. These permissions are granted for free by
  Elsevier for as long as the COVID-19 resource centre remains active." The grant is conditional and
  revocable, so this note quotes the recipe with attribution and the package must not ship it as
  data.
- De Saeger 2022, the bioRxiv preprint: "CC-BY-NC 4.0 International license". The published version
  (ACS Synth. Biol. 11, 2214–2220, doi:10.1021/acssynbio.2c00072) returned 403 and its licence line
  is unread.
- Addgene, Terms of Use: "Content may not be reproduced, duplicated, copied, sold, resold, or
  otherwise exploited for any commercial purpose without the express written consent of Addgene…
  You are free to use Content for your informational, noncommercial purposes, so long as you retain
  all copyright, Marks…, and other proprietary notices contained on such Content."
- Takara, printed in every manual read: "Our products are to be used for Research Use Only… Our
  products may not be transferred to third parties, resold, modified for resale, or used to
  manufacture commercial products or to provide a service to third parties without our prior
  written approval." Site Terms of Use: "You must not reproduce, distribute, modify, create
  derivative works of, publicly display, publicly perform, republish, download, store or transmit
  any of the material on our Website, except that you may print or download one copy of each page
  of the Website for your own personal, non-commercial use".

**The verdict is the same one the earlier notes reached, with one addition.** Vendor numbers are
facts and may be stated in our code and protocols with a citation; vendor documents may not be
mirrored. The addition is that NEB's protocols.io postings carry the two reaction tables under
CC BY, so those two tables — and nothing else from NEB — may be reproduced verbatim with
attribution. "In-Fusion", "Snap Assembly", "Gibson Assembly" and "NEBuilder" are trademarks
(GIBSON ASSEMBLY® is Codex DNA's, per the NEB manuals' trademark block) and must not be used as a
name for anything this package produces.

The downloaded documents are under `reference_docs/gibson/`, which is git-ignored.

## 2. What the three products are, and how they differ

**NEBuilder HiFi and the Gibson Assembly Master Mix** share a mechanism. Both manuals describe
three activities in one buffer, in the same words (E2621 manual v6.0, p. 2; E2611 manual v3.0,
p. 2):

> The exonuclease creates single-stranded 3´ overhangs that facilitate the annealing of fragments
> that share complementarity at one end (the overlap region) … The polymerase fills in gaps within
> each annealed fragment … The DNA ligase seals nicks in the assembled DNA.

The difference NEB states is the polymerase: "The NEBuilder HiFi DNA Assembly Master Mix utilizes a
high-fidelity polymerase… the assembled products from NEBuilder HiFi DNA Assembly Master Mix and
NEBuilder HiFi Cloning Kit will typically result in more colonies with higher accuracy" (E2621
manual, FAQ 3). The same FAQ adds "There are also no licensing fee requirements from NEB with the
NEBuilder products", and the E2611 manual's introduction records why: "Gibson Assembly was
developed by Dr. Daniel Gibson and his colleagues at the J. Craig Venter Institute and licensed to
NEB by Synthetic Genomics, Inc."

**In-Fusion is a different mechanism and must not inherit NEB's numbers.** Takara's own
description (In-Fusion Cloning FAQs, retrieved 2026-09-18): "The In‑Fusion enzyme mix generates
single-stranded 5' overhangs at the termini of the cloning insert and linearized cloning vector.
These overhangs are annealed at the sites of complementarity, and the recombinant circular
construct is rescued in E. coli." Takara states one limit outright — "In‑Fusion Cloning does not
allow for the covalent assembly of linear DNA molecules" — so a linear product is not an option and
a circular one is finished by the host rather than by the tube. Whether a circular In-Fusion product
leaves the tube covalently sealed is asserted only by NEB, whose brochure scores In-Fusion "No" in a
"Covalently sealed?" column; that is a vendor's claim about a competitor and is treated as such
here.

One more In-Fusion constraint with no NEB counterpart, from the FAQ: "Non-phosphorylated
oligonucleotides are compatible with In‑Fusion Cloning. However, 3' exonuclease activity in the
In‑Fusion enzyme mix requires terminal 3' OH groups."

## 3. Overlap length

The ticket's question 1. Each product states a rule, and the three rules disagree; no product
states a rule that depends on fragment **length**, only on fragment **count**.

| Product | 2–3 fragments | 4–6 fragments | Documented minimum | Documented maximum |
| --- | --- | --- | --- | --- |
| NEBuilder HiFi | 15–20 bp | 20–30 bp | 12 bp ("productive"), 15 bp recommended | 30 bp in the current manual |
| Gibson Assembly Master Mix | 15–25 nt | 20–80 bp with the 1 h protocol | 12 bp ("productive"), 15 bp recommended | 25 bp at 15 min, 40 bp at 60 min |
| In-Fusion | 15 bp | 20 bp (more than two fragments) | 12 bp, "not recommended" below | 21 bp, "not recommended" above |

Quoted, so the bands can be defended:

- NEBuilder HiFi, reaction-table footnotes (E2621 manual v6.0, p. 8): "To achieve optimal assembly
  efficiency, design 15–20 bp overlap regions between each fragment." (2–3 fragments) and "To
  achieve optimal assembly efficiency, design 20–30 bp overlap regions between each fragment with
  equimolarity of all fragments (suggested: 0.05 pmol each)." (4–6 fragments). FAQ 24 repeats it:
  "The 15 minute assembly reaction protocol is recommended for assembly of 2–3 fragments that are
  flanked by 15–20 nt overlaps. The 1 hour assembly protocol is recommended for the assembly of 4+
  fragments, flanked by 20–30 nt overlaps."
- NEBuilder HiFi, FAQ 7: "Productive assembly has been achieved for DNA fragments with as little as
  a 12 bp overlap, however, it depends on the GC content of the overlap. We recommend using at
  least 15 bp overlaps, or more, for dsDNA assembly with a Tm ≥ 48°C (AT pair = 2°C and GC pair =
  4°C). Increasing the length of overlap between fragments also reduces the amount of DNA needed
  for assembly."
- NEBuilder HiFi, FAQ 8: "Both the quantity of 5´ exonuclease in the NEBuilder HiFi DNA Assembly
  Master Mix and a 15 minute (recommended for 2-3 fragments) assembly reaction time have been
  optimized for the assembly of DNA molecules with 15-20 bp overlaps. If assembly reaction time is
  increased to 60 minutes (recommended for 4-6 fragments), overlaps of 20-30 bp may be used".
- NEBuilder HiFi at the extremes, FAQ 3: "When large DNA (> 10 kb) or multiple fragments (4+) need
  to be assembled, increasing the overlap region to 30 bp improves the efficiency of assembly and
  transformation."
- Gibson Assembly Master Mix (E2611 manual v3.0, p. 3): "we suggest using a 15–25 nt overlap with a
  Tm equal to or greater than 48°C". Its FAQ 6: "The quantity of 5´ exonuclease in the Gibson
  Assembly Master Mix and a 15 minute assembly reaction time have been optimized for the assembly
  of DNA molecules with ≤ 25-bp overlaps. If assembly reaction time is increased to 60 minutes,
  overlaps up to 40-bp may be used". Its FAQ 23: "The 1 hour assembly protocol is recommended for
  the assembly of up to 6 fragments, flanked by 20–80 nt overlaps."
  **Those last two disagree inside one manual**: FAQ 6 caps a 60-minute Gibson reaction at 40 bp
  and FAQ 23 allows 80 nt for the same protocol. Version 3.0_1/26 carries both. Take 40 bp as the
  figure with a mechanism behind it — it is the one tied to the exonuclease quantity — and record
  that NEB's own manual says 80 elsewhere.
- In-Fusion (In-Fusion Snap Assembly User Manual, 060822, p. 3): the enzyme "fuses DNA fragments …
  efficiently and precisely by recognizing 15-bp overlaps at their ends". p. 12: "In-Fusion allows
  you to join two or more fragments … as long as they share 15 bases of homology at each end."
- In-Fusion at higher fragment counts, same manual, p. 12: "When joining more than two fragments
  (including the linearized vector), we strongly recommend increasing the homologous region to
  20 bp. We have found that this modification yields 5- to 7-fold more transformant colonies, while
  still maintaining high cloning accuracy."
- In-Fusion's hard bounds, In-Fusion Cloning FAQs: "Current In‑Fusion Cloning reaction conditions
  favor a 15-bp homologous overlap for single-insert cloning, and a 20-bp homologous overlap for
  multiple-insert cloning. We do not recommend using overlaps shorter than 12 bp or longer than
  21 bp."
- In-Fusion requires the overlap to be terminal, and says so: "15-bp complementary regions must be
  located at the termini of adjacent DNA fragments or they will not be joined by In‑Fusion Cloning."

**Where the overlap is taken from.** Both NEB manuals allow either fragment to carry it, or a
split: "the overlap sequence can be composed of nucleotides which belong to only one fragment … or
it can be split between the two adjacent fragments in any combination" (E2621 manual, p. 3). Both
add the rule that matters for a reusable backbone: "if the same PCR-generated vector will be used
for assembly of various inserts, then the entire overlap sequence must originate from the vector
sequence and must be added to primers that will be used to amplify the insert", and the same is
required whenever the vector is cut with a restriction enzyme rather than amplified (E2621 manual,
pp. 5–6). Takara allows the same choice — "between two adjacent fragments, only one homologous
overlap is required for the In‑Fusion reaction. This overlap can be located on either of the
fragments" — with one restriction: "splitting the overlap between an insert and vector can only be
done if the vector is linearized via inverse PCR."

**One NEB revision to watch.** The E2621 manual's own revision history records that revision 3.0
(6/20) "Updated sizes of DNA fragments with varied overlaps to (15–30 bp)". NEB's protocols.io
posting, last modified 2022-02-23, still says the method "has been used to assemble either
single-stranded oligonucleotides or different sizes of DNA fragments with varied overlaps (15–80
bp)". Follow the manual, and record which revision was followed.

## 4. Overlap Tm, and its relation to the incubation

The ticket's question 2.

- **Both NEB products specify a Tm for the overlap itself, and specify how to compute it.** E2611
  manual v3.0, p. 3: "we suggest using a 15–25 nt overlap with a Tm equal to or greater than 48°C
  (assuming A-T pair = 2°C and G-C pair = 4°C)". The E2621 manual says the same with its own length
  band. That parenthesis is the Wallace rule, not a nearest-neighbour Tm, and it is the only method
  NEB names for an overlap. By that rule a 15 nt overlap clears 48 °C only when at least nine of its
  bases are G or C, and a 40%-GC sequence needs 18 nt. (That arithmetic is derived here from NEB's
  own formula, not quoted, and it is why the length band and the Tm floor have to be applied
  together.)
- **The length rule and the Tm rule are stated together and are not independent.** NEB gives both
  in one sentence, and the worked example in both manuals lengthens the overlap until the Tm rule
  is met: "The length of the overlap sequence is determined by the number of nucleotides needed to
  reach a Tm ≥ 48°C" (E2621 manual, p. 5).
- **Nothing states the relation between the overlap Tm and the 50 °C incubation.** A 48 °C floor
  sits just under the 50 °C incubation, and it is tempting to read one as derived from the other.
  No source read says so. Treat the connection as unstated.
- **In-Fusion states no Tm for the overlap at all.** Its only Tm rule is for the gene-specific part
  of the primer (§5). Takara is explicit about what it does and does not know about overlap
  composition; see §10.

## 5. Primer shape

The ticket's question 3, which also asks whether `docs/research/primer-design-and-pcr.md` is right.
**It is right, and it can be sharpened.**

- Both NEB manuals describe the same two-part primer: "an overlap sequence, required for the
  assembly of adjacent fragments" added "at the 5´-end of the primer", plus "a gene-specific
  sequence, required for template priming during PCR" added "at the 3´-end of the primer after the
  overlap sequence" (E2621 manual, p. 3).
- The Tm is read from the 3′ part, which is what the primer note already quotes: "The Tm of the 3´
  gene-specific sequence of the primer can be calculated using the Tm calculator found on the NEB
  website" (E2621 manual, p. 3, pointing at `tmcalculator.neb.com`; the E2611 manual points at the
  older `neb.com/TmCalculator` URL). So the Ta for the PCR comes from the annealing region, as
  `liulab_mbio.primers` already models it.
- **The length band NEB gives is on a figure, and the figure is in the older manual.** The E2611
  manual **version 2.0_10/21**, Figure 4A, labels the two parts "Overlap Region (≥15 nt)" and
  "Gene-specific Sequence (18-25 nt)". The same figure in versions 3.0_1/26 and in the E2621 manual
  is an image whose text does not extract, so 18–25 nt is cited to version 2.0 and no further.
- **Neither NEB product states a total primer length.** Not found in any source read.
- **In-Fusion states the most, and states it as hard bands** (In-Fusion Snap Assembly User Manual,
  060822, p. 12). The 3′ portion should "be between 18–25 bases in length and have a GC-content
  between 40–60%", should "have a melting temperature (Tm) between 58–65°C. The Tm difference
  between the forward and reverse primers should be ≤4°C, or you will not get good amplification",
  and should "not contain identical runs of nucleotides. The last five nucleotides at the 3' end of
  each primer should contain no more than two guanines (G) or cytosines (C)." The same page says
  where the Tm is read from, in the same place NEB does: "The Tm should be calculated based upon
  the 3' (gene-specific) end of the primer, and NOT the entire primer. If the calculated Tm is too
  low, increase the length of the gene-specific portion of the primer until you reach a Tm of
  between 58–65°C."
- Takara also gives the only total-length figure any vendor gives, and it is a purification
  threshold rather than a design limit: primers "longer than 45 nucleotides … may need to be PAGE
  purified; however, we usually find this is unnecessary."
- Hairpins and dimers: Takara's rule is on the primer, not the overlap — "Avoid complementarity
  within each primer to prevent hairpin structures, and between primer pairs to avoid primer
  dimers." NEB says nothing on this for assembly primers; the primer note's position (judge
  structure on the whole oligo, Tm on the annealing region) stands as inference, now with Takara's
  half of it quoted.
- Extra bases between the two parts are allowed by both vendors, for frame or for a restriction
  site: "additional nucleotides may be added between the overlap region and gene-specific sequence
  region of inserted fragment to restore the pre-existing restriction site, to introduce a new
  restriction site, or keep the translation of the fusion protein in frame" (E2621 manual, p. 6);
  Takara says the same and adds that those bases "are not part of the 15 bases of sequence
  homology."

## 6. The reaction, per product

The ticket's question 4. All three reactions are one tube, no thermocycling.

**NEBuilder HiFi** (E2621 manual v6.0, p. 8, and the CC BY protocols.io copy):

| Component | 2–3 fragments | 4–6 fragments | Positive control |
| --- | --- | --- | --- |
| Recommended DNA molar ratio | vector:insert = 1:2 | vector:insert = 1:1 | — |
| Total amount of fragments | 0.03–0.2 pmol, X µl | 0.2–0.5 pmol, X µl | 10 µl |
| NEBuilder HiFi DNA Assembly Master Mix | 10 µl | 10 µl | 10 µl |
| Nuclease-free water | 10−X µl | 10−X µl | 0 |
| Total volume | 20 µl | 20 µl | 20 µl |

**Gibson Assembly Master Mix** (E2611 manual v3.0, p. 8):

| Component | 2–3 fragments | 4–6 fragments | Positive control |
| --- | --- | --- | --- |
| Total amount of fragments | 0.02–0.5 pmol, X µl | 0.2–1 pmol, X µl | 10 µl |
| Gibson Assembly Master Mix (2X) | 10 µl | 10 µl | 10 µl |
| Nuclease-free water | 10−X µl | 10−X µl | 0 |
| Total volume | 20 µl | 20 µl | 20 µl |

Both tables carry the same footnote about scale: "If greater numbers of fragments are assembled,
increase the volume of the reaction, and use additional NEBuilder HiFi DNA Assembly Master Mix"
(E2621); the E2611 wording is "additional Gibson Assembly Master Mix may be required". Both cap
unpurified PCR product at 20% of the reaction volume — 4 µl in a 20 µl reaction.

**In-Fusion** (In-Fusion Snap Assembly User Manual, 060822, p. 7, Table 3):

| Component | Cloning reaction | Negative control | Positive control |
| --- | --- | --- | --- |
| Purified PCR fragment | 10–200 ng | — | 2 µl of 2 kb control insert |
| Linearized vector | 50–200 ng | 1 µl | 1 µl of pUC19 control vector |
| 5X In-Fusion Snap Assembly Master Mix | 2 µl | 2 µl | 2 µl |
| Deionized water | to 10 µl | to 10 µl | to 10 µl |

In-Fusion's reaction is **10 µl, not 20**, and its master mix is **2 µl of a 5X**, not half the
reaction. Two volume constraints come with it: "For reactions with larger combined volumes of
vector and PCR insert (>7 μl of vector + insert), double the amount of enzyme premix, and add dH2O
for a total volume of 20 μl", and, for gel-purified input, "The total volume of purified vector and
insert should not exceed 5 μl."

The In-Fusion Snap Assembly EcoDry format is the same table with the master-mix row removed (the
enzyme is a lyophilised pellet in the tube), per the EcoDry user manual (062022), Table 2.

## 7. Incubation

The ticket's question 5.

| Product | Temperature | Time | Changes with fragment count? |
| --- | --- | --- | --- |
| NEBuilder HiFi | 50 °C | 15 min (2–3 fragments) or 60 min (4–6) | yes |
| Gibson Assembly Master Mix | 50 °C | 15 min (2–3 fragments) or 60 min (4–6) | yes |
| In-Fusion | 50 °C | 15 min | **no** |

- NEB, both manuals: "Incubate samples in a thermocycler at 50°C for 15 minutes when 2 or 3
  fragments are being assembled or 60 minutes when 4-6 fragments are being assembled."
- NEB's range, E2621 FAQ 12: "The reaction has been optimized at 50°C, but it has been shown to
  work at temperatures between 40°C and 50°C." The E2611 wording is the same.
- NEB's longer-incubation rule, E2621 FAQ 11: "Reaction times less than 15 minutes are generally
  not recommended. Extended incubation times (up to 4 hours) have been shown to improve assembly
  efficiencies in some cases. Do not incubate the assembly reaction overnight."
- Takara is the opposite on both points. In-Fusion Snap Assembly User Manual, p. 8: "Incubate the
  reaction for 15 min at 50 °C, then place on ice", with the note "The In-Fusion reaction is
  completed within the required 15-min incubation. Longer incubation times do NOT increase cloning
  efficiency, even with multiple-insert cloning reactions." Its FAQ gives a mechanism: an increase
  in reaction time "may generate uneven single-stranded regions at the ends of the cloning insert
  and vector, resulting in inefficient annealing of the homologous overlaps, thus reducing cloning
  efficiency."
- Nothing in any source makes the incubation depend on fragment **length**.

## 8. How much DNA

The ticket's question 6.

**NEB, both products, quantities** (E2621 manual p. 8; E2611 manual p. 7):

- NEBuilder HiFi: "a total of 0.03–0.2 pmol of DNA fragments when 1 or 2 fragments are being
  assembled into a vector, and 0.2–0.5 pmol of DNA fragments when 4–6 fragments are being
  assembled."
- Gibson Assembly Master Mix: "a total of 0.02–0.5 pmols of DNA fragments when 1 or 2 fragments are
  being assembled into a vector and 0.2–1.0 pmoles of DNA fragments when 4–6 fragments are being
  assembled."
- Both add: "Efficiency of assembly decreases as the number or length of fragments increases."

**NEB, molar ratios.** NEBuilder HiFi's table states "Vector:insert = 1:2" for 2–3 fragments and
"Vector:insert = 1:1" for 4–6, with the footnote "Optimized cloning efficiency is 50–100 ng of
vector with 2-fold molar excess of each insert. Use 5-fold molar excess of any insert(s) less than
200 bp." The Gibson Assembly Master Mix asks for more insert at low fragment counts: "Optimized
cloning efficiency is 50–100 ng of vector with 2-3 fold molar excess of each insert. Use 5-fold
molar excess of any insert(s) less than 200 bp. To achieve optimal assembly efficiency using in 4-6
fragment assemblies, use a 1:1 molar ratio of each insert:vector." The E2621 usage notes restate
the floor as a concentration rule: "the concentration of assembly fragments should be at least 2
times higher than the concentration of vector."

**NEB, ng to pmol.** Both manuals give the same formula and the same two worked values:

> pmol = (weight in ng) x 1,000 / (base pairs x 650 daltons)
> 50 ng of 5,000 bp dsDNA is about 0.015 pmol
> 50 ng of 500 bp dsDNA is about 0.15 pmol

This is the 650-dalton rule that `docs/research/primer-design-and-pcr.md` §6.3 already records, and
it disagrees with NEBioCalculator's own constant (`36.04 + 615.94 × bp`) by about 5.5%. The
disagreement is NEB's, not ours, and both manuals point at NEBioCalculator in the same paragraph as
the 650-dalton formula. **Follow the primer note's decision** — compute with NEBioCalculator's
constant and print which constant was used — and record here that a protocol quoting NEB's own
worked examples will be about 5.5% off whatever the package computes.

The E2621 manual adds a second NEB tool for this job: "We recommend using our web tool, NEBuilder®
Protocol Calculator, available at nebuildercalculator.neb.com, to calculate optimal amounts of
input DNA sequences based on the length and concentration of each input fragment." Its behaviour is
an open gap (§19).

**In-Fusion, quantities.** Takara asks for mass, not moles, and its two manuals word the total
differently — a discrepancy worth carrying:

- In-Fusion Snap Assembly User Manual (060822), p. 7: "good cloning efficiency is achieved when
  using 200-ng combined amount of vector and inserts in a 10-μl reaction, regardless of their
  lengths. More is not better. If the size of the PCR fragment is shorter than 0.5 kb, maximum
  cloning efficiency may be achieved by using less than 50 ng of fragment."
- In-Fusion HD Cloning Kit User Manual (102518), p. 10, and the Snap multiple-insert
  protocol-at-a-glance (071320), p. 3: "50–200 ng of vector and inserts respectively" — i.e. per
  fragment, not combined.
- Per-length bands, from the multiple-insert protocol-at-a-glance and the HD manual: insert
  "<0.5 kb: 10–50 ng, 0.5 to 10 kb: 50–100 ng, >10 kb: 50–200 ng"; vector "<10 kb: 50–100 ng,
  >10 kb: 50–200 ng". These footnotes are absent from the Snap manual's Table 3.

**In-Fusion, molar ratios** (Snap manual, p. 7): "for optimal results under standard conditions use
an insert to vector molar ratio of 2:1. When performing a cloning reaction with two or more
inserts, the molar ratio of each of the multiple inserts should still be 2:1 with regards to the
linearized vector … The molar ratio of two inserts with one vector should be 2:2:1." The
multiple-insert protocol-at-a-glance gives three exceptions: 1:1 "if an insert is large with
respect to your linearized vector"; 3–5:1 "for cloning small DNA fragments (150–350 bp)"; and
5–15:1 "for cloning of short synthetic oligos (50–150 bp)", where "the optimal molar ratio must be
determined empirically."

**In-Fusion gives no ng-to-pmol formula.** Its molar-ratio calculator page states the totals —
"The total DNA amount (insert + vector) provided by the calculator is 200 ng, which is optimal for
a 10-µl In-Fusion Cloning reaction", "We recommend using at least 50 ng of a DNA insert", and
"Seamless cloning is not recommended for inserts less than 50 nucleotides" — and its inline script
splits 200 ng between vector and inserts in proportion to length at a fixed 2:1 molar ratio. The
page also caps the tool: "The tool currently supports cloning reactions with up to five inserts."

## 9. Fragment count

The ticket's question 7.

| Product | Demonstrated | Recommended ceiling | What changes above 2–3 fragments |
| --- | --- | --- | --- |
| NEBuilder HiFi | eleven 0.4 kb inserts | **5 inserts** | overlap 15–20 → 20–30 bp, 15 → 60 min, 0.03–0.2 → 0.2–0.5 pmol, ratio 1:2 → 1:1 |
| Gibson Assembly Master Mix | twelve 0.4 kb inserts | **5 inserts** | overlap 15–25 → 20–80 nt, 15 → 60 min, 0.02–0.5 → 0.2–1 pmol, ratio 2–3:1 → 1:1 |
| In-Fusion | five 1 kb inserts | **5 inserts** (the calculator's cap) | overlap 15 → 20 bp, plate 1/5–1/3 of the transformation instead of the usual 1/100–1/5; incubation and reaction unchanged |

- NEBuilder HiFi, FAQ 5: "NEBuilder HiFi DNA Assembly Master Mix has been used to efficiently
  assemble up to eleven, 0.4 kb inserts into a vector at one time. However, we recommend the
  assembly of five or fewer inserts into a vector in one reaction, in order to produce a clone with
  the correct insert. A strategy involving sequential assembly can be used if all of the fragments
  cannot be assembled in a single reaction."
- Gibson Assembly Master Mix, FAQ 3: the same sentence with "up to twelve 0.4 kb inserts".
- The **five-insert cut-off the Golden Gate note cites from the E1601 manual is the same number**,
  and it is stated in NEB's own assembly manuals as well. It is a count of *inserts*, so six pieces
  including the vector.
- In-Fusion, FAQs: "we have successfully tested multiple-fragment cloning with up to five inserts"
  and "cloning up to four 1-kb fragments simultaneously is as easy as cloning a single fragment."
  What changes above two fragments is the overlap (§3) and the plating volume: "For cloning
  reactions with more than two fragments, we recommend plating a larger volume (1/5–1/3 of each
  transformation reaction)" (Snap manual, p. 9).
- Assembled size, which is a separate ceiling: NEBuilder HiFi "has been used to clone a 12 kb DNA
  fragment into a 7.4 kb plasmid in E. coli, totaling up to 19 kb"; the Gibson Assembly Cloning Kit
  "has been used to clone a 15 kb DNA fragment into a 5.4 kb plasmid … totaling up to 20.4 kb".
  Both switch host above 10 kb (NEBuilder) or 15 kb (Gibson) to NEB 10-beta. In-Fusion: "DNA
  inserts up to 15 kb have been successfully cloned into pUC19".

## 10. Mismatch, repeat and GC tolerance

The ticket's question 8.

**Mismatches.** Only NEBuilder HiFi documents tolerating them, and it is the product's selling
point. The manual (p. 6): "One of the unique features of the NEBuilder HiFi DNA Assembly Master Mix
is the ability to remove both 3´ and 5´ end flap sequences upon fragment assembly… This allows
fragments generated by restriction enzyme digestion to assemble while eliminating the remaining
restriction site sequences on both the 5´ and 3´ ends in the fragment junction." FAQ 25 covers the
common case: "The additional A base at the 3´ end of PCR product will be removed during DNA
assembly if it becomes a mismatched residue once fragments anneal." The brochure's comparison table
scores "3´- and 5´-end mismatch" assembly as `+++` for NEBuilder HiFi, `++` for GeneArt Gibson
Assembly and "does not perform" for In-Fusion Snap — vendor comparison data, from NEB.

No source read documents a mismatch **inside** an overlap for any of the three products. Takara's
manuals and FAQ say nothing about mismatch tolerance at all.

**Repeats.** Both NEB products give the same answer, and it is a design rule rather than a
tolerance: "one must ensure that each DNA fragment includes a unique overlap so that the sequences
may anneal and are properly arranged… If having repetitive sequences at the ends of each fragment
is unavoidable, the correct DNA assembly may be produced, albeit at lower efficiency than other,
unintended assemblies" (E2621 FAQ 6; E2611 FAQ 4 is the same). Both also give a concrete rule for
the commonest repeat, FAQ 17/15: a 15-nt overlap made entirely of His-tag repeats is refused — "you
must flank the His-tag sequence on both sides with at least 2 nucleotides that are not part of the
His-tag repeating sequence… You should avoid repeating sequences at the end of an overlap"; the
E2621 version adds "Alternatively, intersperse CAC and CAT his codons to interrupt this repetitive
sequence."

Takara's position is about the vector rather than the junction: "Internal recombination events at
sites other than those adjacent to the vector linearization site are extremely rare. Therefore,
even if your desired region of homology is present more than once in the vector sequence, unwanted
recombination events are unlikely to occur", and In-Fusion is "routinely" used with LTR- and
ITR-carrying viral vectors.

**Palindromes.** The one quantified junction penalty either NEB manual gives, from both
troubleshooting sections: "Avoid overlaps with highly palindromic sequences, as they may cause up
to a 10-fold reduction in recombinant colonies. When assembling fragments into a multiple cloning
site (MCS) of a cloning vector, it is strongly recommended that restriction endonuclease sites be
located at the edges of the MCS to avoid overlap regions with highly-palindromic sequences."

**GC extremes.** NEB ties GC only to the minimum overlap length, through the Tm rule and through
FAQ 7's "it depends on the GC content of the overlap". Takara is the only vendor with numbers, and
disowns them in the same breath: "We have no specific data showing variability of current In‑Fusion
Cloning master mix performance depending on the GC content of the 15-bp overlap. However, the
following results were obtained using a previous version of the kit—In‑Fusion Advantage: 15-bp
homologous overlaps with GC content of 20–40% had little or no effect on the In‑Fusion Advantage
cloning efficiency. 15-bp homologous overlaps with GC content of 60–80% showed a reduced In‑Fusion
Advantage cloning efficiency in certain cases."

**Junction fidelity.** Takara, FAQ: "We have not seen any base slippage, base addition, or base
deletion with the In‑Fusion Cloning enzyme. We have cloned and sequenced over 4,000 separate clones
… and have rarely seen any evidence of errors at the cloning junctions (<2%)."

## 11. Template removal and cleanup

The ticket's question 9.

**DpnI is NEB's answer, and both manuals dose it identically** (E2621 manual p. 7; E2611 manual
p. 7). It is conditional, not routine: "When using circular plasmid DNA as a template, it is
important to use a minimal amount of DNA (e.g., 0.1–0.5 ng of plasmid template per 50 µl PCR
reaction) in order to reduce the template background after transformation. If higher amounts of
plasmid template must be used… it is recommended to digest the PCR product with DpnI".

> 1. In a total 10 µl reaction, mix 5–8 μl of PCR product with 1 μl of 10X CutSmart® Buffer and
>    1 μl (20 units) of DpnI.
> 2. Incubate at 37°C for 30 minutes.
> 3. Heat-inactivate DpnI by incubating at 80°C for 20 minutes.

The reason DpnI works is the same one the Golden Gate note records: DpnI cuts only Dam-methylated
DNA, so plasmid template is cut and the PCR product is not — stated here by NEB itself: "DpnI
cleaves only E. coli Dam methylase-methylated plasmid DNA, but does not cleave the PCR product,
since it is not methylated." Unlike Golden Gate, **this is NEB's own prescription**, not ours.

For a vector that is a PCR product, DpnI is also NEB's answer to empty-vector background: "To
significantly reduce the background of unwanted vector-only colonies, the vector should be a PCR
product, rather than a restriction fragment. If background continues to be a problem, the PCR
amplified vector can be treated with DpnI to remove the template carry-over" (E2621 FAQ 21).

**Cleanup is not required by NEB, within limits.** "PCR product purification is not necessary as
long as the product is > 90% pure. You can add unpurified PCR product directly from the PCR
reaction into the assembly reaction, for up to 20% of the total reaction volume (e.g., PCR products
should account for 4 µl, or less, in a 20 µl NEBuilder HiFi DNA assembly reaction)." Column
purification is recommended, not required: "Column purification of PCR products may increase the
efficiency of both high-fidelity DNA assembly and transformation by 2–10 fold and is highly
recommended when performing assemblies of three or more PCR fragments or assembling longer than
5 kb fragments." Gel purification is required only for a contaminated product: "If non-specific DNA
fragments are obtained, you will need to purify the target fragment from the agarose gel".

**Takara requires purification and offers an enzyme instead of a column.** FAQ: "Yes, the
PCR-amplified DNA must be purified prior to In‑Fusion Cloning." The Snap manual, p. 7: "If a single
band of the desired size is obtained, you can EITHER spin-column purify … OR treat your PCR product
with Cloning Enhancer… However, if non-specific background or multiple bands are visible on your
gel, isolate your target fragment by gel extraction, then spin-column purify." Cloning Enhancer is
dosed where DpnI is not (Cloning Enhancer Protocol-At-A-Glance, 072823): "Add 2 μl of Cloning
Enhancer to 5 μl of the PCR product", then "37°C 15 min / 80°C 15 min / 4°C Hold", extended to
20 min "If you used more than 100 ng of DNA as a template in the PCR reaction". Takara mentions
DpnI twice and doses it never: "treating the PCR product with DpnI before purification will help to
remove contaminating template DNA."

**Preparing a linear vector.** Both NEB manuals prefer a double digest to a single one — "Double
digestion of vector DNA with two restriction endonucleases is the best approach to reduce the uncut
vector background" — and both warn that a single enzyme leaving a GC-rich four-base overhang can
re-circularise: "avoid restriction enzymes that leave four-base single-stranded overhangs rich in
C/G (i.e., CCGG overhang). These overhangs may self-anneal to form the transformable form of the
vector molecule." For inverse PCR of a vector both NEB manuals say "Generally, 10–100 pg of a vector
is recommended as a template in the inverse PCR reaction." Takara adds a spacing rule for a double
digest: the sites should "have at least 5 bases between them", and warns that "overnight restriction
digests are not advisable."

## 12. The original method, and the home-brew reaction

Gibson, D.G. et al. (2009) is the source for the enzymology and for the recipe a lab can mix
itself. Gibson, D.G. (2011), *Methods in Enzymology* 498, 349–361, restates the same recipe in
cleaner text and is the copy quoted here where the two agree, because the 2009 PDF read was a
third-party mirror whose font mangles micro signs and degree signs.

**The three enzymes, and why a 5′ exonuclease** (Gibson 2009, Figure 1 legend and p. 343):

> T5 exonuclease removed nucleotides from the 5′ ends of double-stranded DNA molecules,
> complementary single-stranded DNA overhangs annealed, Phusion DNA polymerase filled the gaps and
> Taq DNA ligase sealed the nicks. T5 exonuclease is heat-labile and is inactivated during the
> 50 °C incubation.

The same page says why that enzyme choice makes a one-pot reaction possible:

> Exonucleases that recess double-stranded DNA from 5′ ends will not compete with polymerase
> activity. Thus, all enzymes required for DNA assembly can be simultaneously active in a single
> isothermal reaction. Furthermore, circular products can be enriched as they are not processed by
> any of the three enzymes in the reaction.

Taq DNA polymerase works in place of Phusion, "but the latter is preferable as it has inherent
proofreading activity for removing noncomplementary sequences (for example, partial restriction
sites) from assembled molecules" — which is the same property NEBuilder HiFi sells as flap removal.

**5X isothermal (ISO) buffer**, Gibson 2009 Online Methods: "25% PEG-8000, 500 mM Tris-HCl pH 7.5,
50 mM MgCl₂, 50 mM DTT, 1 mM each of the four dNTPs and 5 mM NAD". To make 6 ml, per Gibson 2011:
3 ml of 1 M Tris-HCl pH 7.5, 150 µl of 2 M MgCl₂, 60 µl each of 100 mM dGTP, dATP, dTTP and dCTP,
300 µl of 1 M DTT, 1.5 g PEG-8000 and 300 µl of 100 mM NAD, aliquoted and stored at −20 °C "for up
to 1 year". (The 2009 print of that last figure reads "300 ml"; 2011 reads 300 µl, and 300 µl is
the value consistent with 5 mM NAD in 6 ml.)

**Assembly master mixture**, Gibson 2011 §5 step 2:

> An enzyme–reagent master mixture is prepared by combining 320 μl of 5× ISO reaction buffer,
> 0.64 μl of 10 U/μl T5 exo (Epicentre), 20 μl of 2 U/μl Phusion® pol, 160 μl of 40 U/μl Taq lig
> (NEB), and water up to a final volume of 1.2 ml. Fifteen microliters of this enzyme-reagent mix
> can be aliquoted and stored at − 20 °C for up to 2 years.

**The reaction**: "a 20-μl reaction consisting of 5 μl DNA and 15 μl enzyme–reagent master mixture",
"incubated at 50 °C for 1 h", then "Samples are diluted 1:5 with sterile water" before
transformation. DNA input: "Approximately 10–100 ng of each DNA segment is used in equimolar
amounts. For 5–8 kb DNA fragments, 25 ng substrate DNA is ideal… It is best to use ≤ 1 kb DNA
fragments in 5- to 10-fold excess."

**Exonuclease dose is what changes with overlap length**, not time or temperature — and the two
papers give different thresholds, which a later ticket must not average:

- Gibson 2009 Online Methods: "For overlaps shorter than 150 bp, 0.2 U µl⁻¹ T5 exonuclease is used.
  For overlaps larger than 150 bp, 1.0 U µl⁻¹ T5 exonuclease was used." Of the premix: "The
  exonuclease amount is ideal for the assembly of DNA molecules with 20–150 bp overlaps."
- Gibson 2011 §5 step 2: "This exonuclease amount is ideal for overlaps that are ≤ 80 bp. For
  overlaps that are ≥ 80 bp, 3.2 μl exonuclease is used in the mixture."

**Overlap length in the papers.** The first demonstration used "~450 base pairs"; 40 bp was then
shown to work "when we reduced the concentration of T5 exonuclease". The design recommendation is
Gibson 2011 §2: "DNA molecules are designed such that neighboring fragments contain at least 40 bp
of overlapping sequence. However, as much as 500 bp can be used if the procedures are slightly
modified". **Neither paper states a minimum below 40 bp**; the 15-nt figure everyone quotes is a
vendor and third-party number (§3).

**Incubation**: "Incubations were performed at 50 °C for 15 to 60 min (60 min was optimal)."

**Scale.** Three 5 kb fragments into an 8 kb BAC; 144 kb + 166 kb + an 8 kb BAC into a 318 kb
product; four quarter-genomes into "the complete synthetic 583-kb M. genitalium genome". The stated
limit: "The size limit for in vitro DNA assembly is not known, but products as large as 900 kb have
been observed", with cloning in *E. coli* the practical ceiling at "several hundred kilobases".
**No maximum fragment count appears in either paper.**

## 13. A fourth product: GeneArt Gibson Assembly HiFi

Thermo Fisher's GeneArt Gibson Assembly HiFi is the product #74 pointed at, and it is a clean test
of the spec's claim that a fourth product is a row of data. Its numbers differ from NEB's in every
column, so it belongs in the data and not in the code (user guide MAN0019062 Rev. D, 6 June 2025,
and the article "DNA Cloning Tips – Build Clones with DNA Fragments using Gibson Assembly").

- **Reaction**: "0.08 pmol vector, 0.08 pmol each fragment up to a total of 10 µL combined" plus
  10 µl of GeneArt Gibson Assembly HiFi Master Mix, water to 20 µl. From unpurified PCR: up to 1 µl
  of each, total PCR product ≤ 4 µl.
- **Incubation**: 50 °C, "15 minutes [1-3 inserts] / 60 minutes [4-5 inserts]".
- **Ratio**: "For maximum cloning efficiency, use a 1:1 molar ratio of vector:insert(s)" —
  equimolar at every fragment count, unlike NEB.
- **ng to pmol**: `pmols = (weight in ng) x 1000 / ((fragment length in bp) x 660)`. That is a
  **third constant**: 660 Da/bp against NEB's 650 and NEBioCalculator's `36.04 + 615.94 × bp`. Three
  vendors, three constants, about 7% apart end to end.
- **Overlap**, from the article's table: 1–2 fragments, ≤ 8 kb → 20–40 bp; 8–32 kb → 25–40 bp. 3–5
  fragments, ≤ 8 kb → 40 bp; 8–32 kb → 40–100 bp. 6+ fragments, 100 bp to 100 kb → 50–100 bp. The
  2025 user guide collapses this to 1–5 fragments, ≤ 8 kb → 20–40 bp; > 8 kb → 40–100 bp. **Quote
  whichever is cited and date it**; the two disagree.
- **Primer shape**: "Use PCR primers ~65 nt in length (20-40 nt for the requisite homology at the 5'
  end, and 18-25 nt specific to the DNA element)" — the only source read that states a total primer
  length. Also "Selecting overlapping regions of homology with Tm >50°C can improve efficiency",
  which is a second Tm floor and a different one from NEB's 48 °C.
- **Fragment count**: the HiFi kit is "Capable of assembling up to five DNA fragments plus a vector
  in one reaction"; the white paper adds an EX kit "for assembly of up to 15 fragments".

## 14. Oligo stitching

The ticket's question 10: a part short enough is ordered as overlapping oligos and assembled in the
reaction rather than amplified. Both NEB products take single-stranded oligos directly; In-Fusion
takes only a double-stranded insert. Every number below is different from every other, and the
published ones are the strongest.

**Gibson 2011 §6 is the primary protocol**, citing Gibson 2010 for the method:

| Value | Gibson 2011 |
| --- | --- |
| Oligo length | 60 bases |
| Overlap between neighbours | 20 bp |
| Oligos per reaction | "only eight to twelve 60-base oligos are assembled at one time" |
| End oligos | carry a 20-bp overlap to the PCR-amplified pUC19 termini |
| Worked example | eight 60-mers → 340 bp total, of which "only 284 bp of unique sequence… is synthesized" |
| Oligo stock | 50 µM in TE, "pooled in groups of 8 or 12 and diluted in TE buffer to a per-oligo concentration of 180 or 75 nM, respectively" |
| Reaction | 5 µl DNA + 15 µl master mix, 50 °C for 1 h, then diluted 1:5 |

Two things to carry into code. First, **there is no separate annealing step or ramp** — the oligos
go into the same isothermal reaction as everything else, and nothing in Gibson 2011 anneals them
first. Second, the reason for the 8–12 ceiling is error rate, not chemistry: "To ensure that
error-free molecules are obtained at a reasonable efficiency, only eight to twelve 60-base oligos
are assembled at one time", and "the errors originating from the chemical synthesis of the oligos
are weeded out by DNA sequencing".

For longer builds the papers stage the work rather than adding oligos: Gibson 2010's abstract
describes "a one-step, isothermal assembly method for synthesizing DNA molecules from overlapping
oligonucleotides. The method cycles between in vitro recombination and amplification until the
desired length is reached", and used it to synthesise "the entire 16.3-kilobase mouse mitochondrial
genome from 600 overlapping 60-mers". That staged method is not this one, and its full text is
closed (§19).

**The vendors' numbers**, all different:

- NEB, both products, FAQ 10/8: "Can ssDNA oligonucleotides be assembled with dsDNA fragments? Yes.
  However, the optimal concentration of each oligonucleotide should be determined. As a starting
  point, we recommend using 45 nM of each oligonucleotide that is less than or equal to twelve
  60-base oligonucleotides containing 30-base overlaps." Desalted oligos are fine: "Standard,
  desalted primers may be used."
- De Saeger et al. 2022, peer-reviewed, using NEBuilder HiFi: de novo assembly from oligos alone,
  "assembled with 1 μL of each 100 nM oligo in a 10 μL NEBuilder assembly reaction", oligos
  "resuspended to 100 µM and diluted to 0.3 pmol/µl". Its own ceiling, measured: "The cloning yield
  dropped with five oligos, but the cloning efficiency was 60-100% for all four designs…, making
  this a viable approach to building DNA fragments at least ~270 bp long."
- Thermo Fisher's white paper (`COL24373 0920`) uses shorter oligos in an otherwise normal assembly:
  30-nt stitching oligonucleotides, "Each stitching oligonucleotide contained evenly split homology
  between adjacent fragments" — 15 bp each side — with fragments and vector at 0.04 pmol each in
  20 µl and 15 min at 50 °C. Its recommendation is a concentration, and it is measured: "For optimal
  cloning efficiency with 30 nt stitching oligonucleotides, we recommend using a concentration
  between 30 nM and 45 nM", from a series at 20, 30, 45 and 120 nM scoring 64%, 94%, 92% and 74%.
- In-Fusion does not do this. Takara's oligo route is a double-stranded insert: "The smallest insert
  successfully cloned with In‑Fusion Cloning was a 50-bp synthetic oligonucleotide (including two
  15-nt homologous overlaps with the vector termini)", at an "oligo-to-vector molar ratio… 5–15:1"
  for 50–150 bp, and "Seamless cloning is not recommended for inserts less than 50 nucleotides."

**The size ceiling above which a PCR is recommended instead** is stated only by Addgene, and as a
window rather than a number: stitching "is used when the part to be inserted is too long to include
on overlapping PCR primers (>60 bp) but too short to make its own part (<150 bp)". Gibson's own
ceiling is the oligo count, and the arithmetic that follows from it is: *n* oligos of length *L*
overlapping by *v* tile *n*(*L*−*v*) + *v* bases, which reproduces Gibson's worked 340 bp from
8 × 60-mers at 20 bp, and gives 390 bases for NEB's twelve 60-mers at 30-base overlaps. That
arithmetic is derived here, not quoted.

## 15. Bridging oligo

The ticket's question 11: one oligo carrying homology to two fragments that share none.
**NEBuilder HiFi has NEB's own protocol for it, peer-reviewed evidence and a bench protocol; Thermo
documents the same trick for GeneArt; the Gibson Assembly Master Mix documents ssDNA oligos in the
reaction but not this use of them; In-Fusion documents nothing.**

**NEB's own protocol** is the application note *Bridging dsDNA with a ssDNA Oligo and NEBuilder HiFi
DNA Assembly to create an sgRNA-Cas9 Expression Vector* (Hsieh, P., NEB, dated 06/17):

- Design: "Design an ssDNA oligo containing the target sequence (19-21 bases) of sgRNA flanked by
  25 bases of sequence at both ends." The worked oligo is 70 bases, and counts 25 + 20 + 25.
- Reaction: "Prepare the ssDNA oligo in 1X NEBuffer 2 to a final concentration of 0.2 µM. Assemble a
  10 µl reaction mix with 5 µl of ssDNA oligo (0.2 µM), 30 ng of restriction enzyme-linearized
  vector and ddH2O. Add 10 µl of NEBuilder HiFi DNA Assembly Master Mix to the reaction mix, and
  incubate the assembly reaction for 1 hour at 50°C." That is 1 pmol of oligo against 30 ng of a
  9,424 bp vector — a large molar excess, and the only dose NEB publishes for this route.
- Result: "Of the 10 clones sequenced, 9 contained the target sequence in the correct orientation."

NEB's product literature states the capability in general terms: the NEBuilder brochure lists "the
ability to bridge two double-stranded DNA fragments with a single-stranded DNA oligo (data not
shown)" among its advantages over GeneArt Gibson Assembly and In-Fusion Snap, and the NEB GmbH
product page words it as "Bridge two double-stranded fragments with a synthetic single-stranded DNA
oligo for simple and fast construction (e.g. linker insertion or gRNA libraries)".

**The peer-reviewed evidence is De Saeger et al. 2022**, the paper the Jacobs lab points at — its
protocols page carries no Gibson protocol, only the sentence "Oligo stitching relies on the use of
ssDNA oligos in a Gibson assembly reaction. Hence, the method is PCR-free. We reported this method
in ACS Synthetic Biology in De Saeger et al., 2022." From the paper:

> We used NEBuilder HiFi master mix with two 44-nt oligos, one for the left and one for the right
> flank… The oligos match 20 bp of the backbone on one end and 20 bp of the part on the other end,
> with the 4 bases in the middle forming the new Golden Gate overhang… We obtained hundreds to
> thousands of clones per assembly using homemade chemically-competent DH5α cells. The resulting
> cloning efficiency – i.e., the ratio of correct clones - ranged from 80 to 100%.

Its reaction: 10 µl, half master mix, "0.04 pmol was used for each of the insert(s) and the
backbone. Oligos were designed to have homology with 20 bp at each side of the junction", 50 °C for
1 h. Two design facts a checker can use: the rule of thumb is "n+1 oligos (where n equals the number
of parts)", and orientation does not matter — "all assemblies are equally efficient and accurate…,
indicating that the orientation of the oligos does not affect cloning efficiency."

**A third bench protocol**, by Thomas Jacobs (April 2015, hosted by Addgene), uses a 60-mer with
20-bp overlaps on each side at 0.2 pmol, equimolar with two PCR amplicons and ~100 ng of a
14,349 bp vector, in 10–20 µl for 1 hour. It carries a caveat worth keeping: "I have frequently
observed SNPs immediately adjacent to the overlaps with Gibson; NEBuilder is highly recommended",
and "Pools of oligos can be used in a single reaction… I've pooled 11 oligos with success."

**Thermo Fisher describes the same trick in its own words**, as the general case rather than a
special one: overlaps "can be introduced using stitching oligonucleotides that share half of the
sequence with one fragment and the other half with the adjacent fragment. In this way, the
oligonucleotides work like a bridge between two DNA fragments", which "enables virtually any
possible combination between DNA fragments that do not originally share any homology". Its count
rule: "One-fragment and two-fragment cloning require two and three stitching oligonucleotides,
respectively" — the same *n*+1 as De Saeger.

**Homology per side, across the four sources**: 25 bases (NEB), 20 bp (De Saeger, Jacobs), 15 bp
(Thermo's 30-nt oligo). Nothing below 15.

**In-Fusion: not supported, and record that as the answer.** No Takara source read mentions a
bridging or linker oligo. Two documented facts explain why it should not be attempted blind: "3'
exonuclease activity in the In‑Fusion enzyme mix requires terminal 3' OH groups", and "In‑Fusion
Cloning does not allow for the covalent assembly of linear DNA molecules." NEB's comparison table
scores "ssOligo & vector" as not performing for In-Fusion Snap; that is NEB's data about a
competitor and is weaker evidence than Takara's own silence.

## 16. Expected result

The ticket's question 12.

**Vendor release criteria**, the firmest numbers because they are what the lot is tested against:

- Gibson Assembly Master Mix, PS-E2611S/L v2.0, effective 19 Feb 2024: six 0.05 pmol fragments of
  pUC19 (five of 400 bp and one of 2,780 bp, each with a 40 bp overlap), 50 °C for 60 min, 2 µl
  transformed into NEB 5-alpha, "yields greater than 100 white colonies on an ampicillin plate with
  IPTG/X-Gal after overnight incubation at 37°C".
- NEBuilder HiFi, PS-E2623S v1.0 and the E2621 manual's Specification section: six fragments at
  0.05 pmol each — four of 1,000 bp, one of 1,152 bp with an 80 bp overlap, and a 3,373 bp vector
  with a 20 bp overlap — 50 °C for 60 min, and "Greater than 100 blue colonies were observed when
  1/10 of the outgrowth (500 µl) was spread on a plate". (NEB's own transformation protocol makes
  the outgrowth 1 ml, so "1/10" and "500 µl" cannot both be right; quote the criterion, not the
  parenthesis.)
- The screen is **reverse blue/white in the NEBuilder test system** — "Successfully assembled
  fragments produce an intact lacZ gene in the pACYC184 vector, and yield blue colonies" — and
  **classic white-is-correct in the Gibson test system**, which uses pUC19. The Golden Gate note
  records the same trap; say which convention a plate is using.

**Fractions of correct clones**, by fragment count and source:

| Source | Assembly | Correct clones |
| --- | --- | --- |
| Gibson 2009 | 3 × 5 kb fragments into an 8 kb BAC | "nine out of ten colonies tested had the predicted 15-kb insert", from 4,500 colonies |
| Gibson 2009 | 144 kb + 166 kb + 8 kb BAC | "several hundred clones, and 5 out of 10 colonies screened had the correct insert size" |
| Gibson 2009 | junction fidelity, two-step thermocycled | "30 cloned DNA molecules (210 repaired junctions)… revealed only 4 errors… about 1 error per 50 DNA molecules joined" |
| In-Fusion, Takara FAQ | 2 × 1 kb | 2,128 colonies at 1/5 plated, 10/10 correct |
| In-Fusion | 3 × 1 kb | 83 colonies, 7/10 |
| In-Fusion | 4 × 1 kb | 31 colonies, 8/10 |
| In-Fusion | 5 × 1 kb | 14 colonies, 4/10 |
| In-Fusion, general | single insert | "Cloning efficiency is at least 95% for a single insert into a vector" |
| GeneArt HiFi white paper | 4 × 2 kb into a vector, 11 kb product, 30 bp overlaps, 0.04 pmol each | "average cloning efficiency of three independent assembly experiments was 96%… with an average of 1,617 colony forming units" |
| GeneArt EX article | 25 kb and 50 kb inserts, 50 bp homology | 3/16 and 1/16 full-length; "cloning efficiency was low (between 6 and 20%)" |
| De Saeger 2022, NEBuilder HiFi + bridging oligos | part conversion | "hundreds to thousands of clones per assembly"; 80–100% correct |
| NEB app note, NEBuilder HiFi + bridging oligo | sgRNA into a 9.4 kb vector | 9/10 correct |

Two conclusions a plan can state honestly. Correct-clone fraction falls with fragment count — the
In-Fusion table is the only source that measures the fall across a series, from 10/10 at two
fragments to 4/10 at five — and colony count falls much faster than accuracy does, by two orders of
magnitude across the same series.

**What goes into the transformation.** NEB: "Add 2 µl of the chilled assembled product to the
competent cells", 30 min on ice, "Heat shock at 42°C for 30 seconds", 2 min on ice, 950 µl SOC,
37 °C 60 min at 250 rpm, "Spread 100 µl of the cells onto the selection plates", overnight at 37 °C.
Takara: "Add 2.5 µl of the In-Fusion reaction mixture to the competent cells" with a hard cap — "DO
NOT add more than 5 μl of the reaction to 50 μl of competent cells. MORE IS NOT BETTER. Using too
much of the reaction mixture inhibits the transformation" — and a competence floor, "competent cells
with a transformation efficiency ≥1 x 10^8 cfu/ug". NEB's own cells are specified at "1–3 x 10^9
colonies formed/μg" of pUC19.

## 17. Troubleshooting

The ticket's question 13. Both NEB manuals carry the same five-heading table and Takara carries
three; the rows below are theirs, condensed to the symptom each is filed under.

| Symptom | NEB (E2621 and E2611 manuals) | Takara (In-Fusion Snap manual, Table 4) |
| --- | --- | --- |
| No colonies | Run the positive control; "Analyze the reaction on an agarose gel. An efficient assembly reaction will show assembled products of the correct size and the disappearance of fragments"; check the overlap is long enough; avoid highly palindromic overlaps (up to a 10-fold loss); raise fragment and vector concentration; keep unpurified PCR at ≤20% of the reaction; consider insert toxicity and a low-copy vector; "Test the success of the DNA assembly by performing PCR with primers that flank the assembled product" | Transformed with too much reaction — never more than 5 µl into 50 µl; cells sensitive to the enzyme, so "dilute the In-Fusion reaction with TE buffer 5–10 times"; cells not competent (≥1 × 10⁸ cfu/µg); "Regions of homology were not long enough for efficient cloning of >2 fragments at once" → 15 bp to 20 bp; DNA concentration too low → 200 ng combined in 10 µl; gel purification carried contaminants → purified vector + insert ≤ 5 µl; primer sequences wrong |
| Clones without the insert (empty vector) | "PCR products may carry over large quantities of uncut plasmid template. To remove plasmid template, treat PCR products with DpnI"; a restriction-cut vector carries uncut plasmid, so double-digest, and "avoid restriction enzymes that leave four-base single-stranded overhangs rich in C/G (i.e., CCGG overhang)"; raise units, time, or gel-purify the linear vector | "Incomplete linearization of your vector… If necessary, recut your vector and gel purify"; template plasmid carried through — linearise the template before PCR, or DpnI-treat the product; "Be sure that your antibiotic plates are fresh (<1 month old)" |
| Wrong product | "Make sure that your PCR product is a single band of the correct size. If the PCR product is contaminated with non-specific bands, it is necessary to gel purify"; consider toxicity and a low-copy vector; "Consider using NEB Stable Competent E. coli (NEB #C3040) for inserts that contain repetitive sequences" | "If your PCR product is not a single distinct band, then it may be necessary to gel purify the PCR product to ensure cloning of the correct insert" |
| Small colonies | "Some recombinant proteins are not well-tolerated by E. coli… Use a low copy number vector (i.e., pACYC184) or a vector with tight control of protein expression. When assembling into the pUC19 vector, make sure that your gene is not in frame with lacZ alpha fragment" | — |
| Deletions between repeats | Not a row in either NEB manual. The nearest statements are FAQ 6's "each DNA fragment includes a unique overlap", the His-tag rule (§10), and "Some DNA structures, including inverted and tandem repeats, are selected against by E. coli" | **No such row exists** in any Takara manual read. Its nearest statement is the junction-fidelity FAQ: "We have not seen any base slippage, base addition, or base deletion with the In‑Fusion Cloning enzyme… rarely seen any evidence of errors at the cloning junctions (<2%)" |

Two rows from elsewhere that a protocol should carry, because they are the fixes a bench scientist
reaches for and no vendor troubleshooting table names:

- Secondary structure in the overlap, from Addgene: "Avoid strong secondary structures in the
  homology region. Hairpins in this region can significantly reduce the efficiency of two homologous
  ends annealing."
- Single-stranded ends being chewed past their target, from Addgene: "Add Extreme Thermostable
  Single-Stranded DNA-Binding protein (ET SSB) to the isothermal reaction mix. ET SSB protects 3'
  ssDNA ends from the ssDNA-specific endonuclease activity of T5 Exonuclease", citing Rabe & Cepko
  2020, which is unread (§19).

Thermo's user guide adds more of the same kind: use two enzymes rather than one to linearise,
extend the digest "(2-3 hours to overnight)", and limit UV exposure of a gel-purified fragment — "leaving the gel on the gel tray when exposing to UV light, using low UV power, and
minimizing the time the gel is exposed" — because damaged ends do not assemble.

## 18. Implications for this package

**For the design module.** The overlap rule is per product and per fragment count, and it is two
rules at once: a length band and a Tm floor computed by the Wallace rule (§3, §4). Design should
take the product's band, extend the overlap until NEB's 48 °C floor is met, and refuse below the
product's documented minimum — 12 bp for the two NEB products, 12 bp for In-Fusion, and nothing
below 15 bp recommended by anyone. The overlap may be taken from either side or split, except where
the vector is reused or restriction-cut, in which case it must come entirely from the vector (§3).
In-Fusion's rule is its own — 15 bp for one insert, 20 bp for more than two fragments, never above
21 — and must not inherit NEB's numbers.

**For the checks.** Sourced bands exist for: overlap length per product and fragment count;
overlap Tm ≥ 48 °C by the AT=2/GC=4 rule; a palindromic overlap, which NEB quantifies at up to a
10-fold loss of colonies; a repeat at the end of an overlap, which both NEB manuals refuse with the
His-tag rule; fragment count against the product's documented maximum of five inserts. Two checks
the spec names have **no sourced threshold**: "two overlaps alike enough to mis-assemble" is
supported only by NEB's qualitative "each DNA fragment includes a unique overlap", and GC extremes
are supported only by Takara's disowned In-Fusion Advantage data. Per `CONTEXT.md`, those two carry
no verdict until something measures them.

**For the bench layer.** Three reaction tables, three incubations, three ratio rules, three
ng-to-pmol constants (§6, §7, §8, §13). Keep the conversion in code with NEBioCalculator's constant
as the primer note decided, and print which constant was used, because a protocol that quotes NEB's
worked examples alongside a computed pmol will differ by about 5.5%, and Thermo's by about 7%.
DpnI has NEB's own dose and time (§11) — unlike Golden Gate, this step is not ours to justify.

**For the two oligo routes.** Stitching has published numbers: 60-base oligos, 20 bp overlaps, 8 to
12 per reaction, no separate annealing step, and a window from Addgene of roughly 60 to 150 bp for
when it beats a primer tail or a synthesis order (§14). Bridging has published numbers too: 20–25 bp
of homology per side, one oligo more than the number of parts, orientation irrelevant, and a large
molar excess of oligo (§15). In-Fusion supports neither as described: its short-insert route is a
double-stranded synthetic oligo of at least 50 bp carrying the usual 15-bp overlaps, and it
documents no bridging oligo at all. A plan that names In-Fusion should refuse both routes rather
than substitute NEB's numbers.

**For the order sheet.** A stitching or bridging oligo primes nothing, so none of the three oligo
roles in `liulab_mbio.primers.thresholds` judges it. What the sources hold such an oligo to is:
standard desalted purity is enough — "Standard, desalted primers may be used" (NEB), "Gel or HPLC
purification of oligonucleotides is not required" (Takara) — a length near 60–70 bases, and a
concentration. No source gives it a Tm, a hairpin or a dimer threshold. A fourth role would have
almost nothing to put in it; a row carrying no verdict, with its concentration and its purity, is
what the sources support.

**For the plan's status.** In-Fusion's measured series (§16) is the only source that quantifies how
a correct-clone fraction falls with fragment count. It is the honest basis for warning a user who
asks for five inserts, and it is In-Fusion's number, not NEB's.

## 19. Open gaps

Still unverified. None of these is guessed above.

| Item | Why it is missing | Workaround in use |
| --- | --- | --- |
| Everything that needed an archived snapshot | `web.archive.org` returned "Internet Archive services are temporarily offline" all day | Live PDFs and mirrors only; no claim here rests on a snapshot |
| NEBuilder Assembly Tool (`nebuilder.neb.com`) behaviour | 403 live, no snapshot | The manuals describe what it designs; the overlap rule is taken from the manual instead |
| NEBuilder Protocol Calculator (`nebuildercalculator.neb.com`) behaviour | 403 live, no snapshot | NEB's own pmol formula and reaction table, quoted in §6 and §8 |
| NEB's HTML protocol and FAQ pages, including the bridging-oligo protocol page | 403 live to `curl`, including with `?pdf=true` | The application note PDF, which carries the same protocol with amounts |
| NEB product specification for E2621/E5520 | No filename spelling found under `/-/media/catalog/specifications/e/2/` that returns 200 | PS-E2623S v1.0, whose functional test is the same assembly, and the E2621 manual's Specification section |
| Gibson 2009 Supplementary Information | Springer ESM endpoint returns 403 | The Online Methods carry the ISO recipe; Supplementary Table 1 (the error count), Table 2 (primers) and Fig. 6 (900 kb products) are unread |
| Gibson 2010 (`nmeth.1515`) full text | Closed access, no repository copy | Gibson 2011 restates the oligo protocol and cites 2010 for it; only the 2010 abstract is quoted |
| Rabe & Cepko 2020, the ET SSB result | bioRxiv returned HTTP 429 twice | Addgene's one-sentence summary, quoted in §17 and marked as second-hand |
| The published licence line of De Saeger 2022 | `pubs.acs.org` returns 403 | The bioRxiv preprint's CC BY-NC 4.0 |
| A mismatch tolerance **inside** an overlap, for any product | No source read documents one | NEBuilder's end-flap removal is about ends, not the overlap interior; no check is proposed |
| A threshold for two overlaps being too alike | No vendor or paper quantifies cross-annealing between overlaps | NEB's qualitative "unique overlap" rule; the check carries no verdict |
| Whether the 48 °C overlap Tm floor is derived from the 50 °C incubation | Stated by nobody | Both numbers are used as given; the connection is not asserted |
| Addgene's protocol page date | No "last updated" or copyright line in the page HTML | Cited by URL with today's retrieval date |

## Sources

All retrieved 2026-09-18.

- NEB, *NEBuilder® HiFi DNA Assembly Master Mix / NEBuilder HiFi DNA Assembly Cloning Kit*
  instruction manual, NEB #E2621S/L/X and #E5520S, version 6.0_1/26:
  [manuale2621_e5520.pdf](https://www.neb.com/-/media/nebus/files/manuals/manuale2621_e5520.pdf)
- NEB, *Gibson Assembly® Master Mix / Gibson Assembly® Cloning Kit* instruction manual,
  NEB #E2611S/L and #E5510S, version 3.0_1/26:
  [manuale2611_e5510.pdf](https://www.neb.com/-/media/nebus/files/manuals/manuale2611_e5510.pdf);
  version 2.0_10/21 at `manuale2611.pdf`, read for the primer figure
- NEB product specifications, fetched under `/-/media/catalog/specifications/<a>/<b>/`:
  PS-E2611S/L v2.0 (19 Feb 2024), PS-E2623S v1.0 (01 Apr 2019)
- New England Biolabs (2022) *NEBuilder HiFi DNA Assembly Reaction (E2621)*, protocols.io,
  [doi:10.17504/protocols.io.bfhrjj56](https://dx.doi.org/10.17504/protocols.io.bfhrjj56) (CC BY)
- New England Biolabs (2022) *Gibson Assembly® Master Mix – Assembly (E2611)*, protocols.io,
  [protocols.io](https://www.protocols.io/view/gibson-assembly-master-mix-assembly-e2611-bdd8i29w) (CC BY)
- Hsieh, P. (NEB) *Bridging dsDNA with a ssDNA Oligo and NEBuilder® HiFi DNA Assembly to create an
  sgRNA-Cas9 Expression Vector*, application note, 06/17:
  [appnote…pdf](https://www.neb.com/en/-/media/nebus/files/application-notes/appnote_bridging_dsdna_with_ssdna_oligo_and_nebuilder_hifi_dna_assembly_to_create_sgrna-cas9_expression_vector.pdf)
- NEB, *NEBuilder® HiFi DNA Assembly: the next generation of DNA assembly and cloning* brochure,
  © 2024, and the NEB GmbH product page, both via the
  [neb-online.de](https://www.neb-online.de/en/cloning-synthetic-biology/dna-assembly/nebuilder-dna-assembly/) mirror
- Gibson, D.G., Young, L., Chuang, R.-Y., Venter, J.C., Hutchison, C.A. and Smith, H.O. (2009)
  Enzymatic assembly of DNA molecules up to several hundred kilobases. *Nat. Methods* 6, 343–345.
  [doi:10.1038/nmeth.1318](https://doi.org/10.1038/nmeth.1318) (closed access; a third-party mirror
  was read)
- Gibson, D.G. (2011) Enzymatic assembly of overlapping DNA fragments. *Methods Enzymol.* 498,
  349–361. [doi:10.1016/B978-0-12-385120-8.00015-2](https://doi.org/10.1016/B978-0-12-385120-8.00015-2),
  full text via Europe PMC, PMC7149801
- Gibson, D.G. et al. (2010) Chemical synthesis of the mouse mitochondrial genome. *Nat. Methods* 7,
  901–903. [doi:10.1038/nmeth.1515](https://doi.org/10.1038/nmeth.1515) (closed access; abstract only)
- De Saeger, J., Vermeersch, M., Gaillochet, C. and Jacobs, T.B. (2022) Simple and efficient
  modification of Golden Gate design standards and parts using oligo stitching,
  bioRxiv [2022.02.10.479870](https://doi.org/10.1101/2022.02.10.479870) (CC BY-NC 4.0); published
  as *ACS Synth. Biol.* 11, 2214–2220,
  [doi:10.1021/acssynbio.2c00072](https://doi.org/10.1021/acssynbio.2c00072)
- Jacobs, T. (2015) *Using NEBuilder with ssDNA oligos*, hosted by Addgene
- VIB Jacobs lab, [protocols page](https://jacobslab.sites.vib.be/en/protocols)
- Addgene, [Gibson assembly protocol](https://www.addgene.org/protocols/gibson-assembly/) and
  [Terms of Use](https://www.addgene.org/terms-of-use/)
- Takara Bio, *In-Fusion® Snap Assembly User Manual* (060822), *In-Fusion® Snap Assembly
  Multiple-Insert Cloning Protocol-At-A-Glance* (071320), *In-Fusion® Snap Assembly EcoDry™ User
  Manual* (062022), *In-Fusion® HD Cloning Kit User Manual* (102518), *In-Fusion® HD Multiple-Insert
  Cloning Protocol-At-A-Glance* (121416), *Cloning Enhancer Protocol-At-A-Glance* (072823), all
  under [takarabio.com/documents](https://www.takarabio.com/documents/)
- Takara Bio, [In-Fusion Cloning FAQs](https://www.takarabio.com/learning-centers/cloning/in-fusion-cloning-faqs),
  the In-Fusion molar ratio calculator, the primer-design tools page and the website
  [Terms of Use](https://www.takarabio.com/terms-of-use)
- Thermo Fisher Scientific, *GeneArt™ Gibson Assembly® HiFi and EX Cloning Kits* user guide,
  MAN0019062 Rev. D, 6 June 2025
- Thermo Fisher Scientific, oligo-stitching white paper `COL24373 0920`, September 2020, and the
  article *DNA Cloning Tips – Build Clones with DNA Fragments using Gibson Assembly*

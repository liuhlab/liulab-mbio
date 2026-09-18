---
search:
  exclude: true
---

# Gateway cloning: att sites, the BP and LR reactions, and what may ship

Research note for issue #121. Everything below was retrieved on **2026-09-18** unless a
different date is given next to the source.

## How to read this note

Thermo Fisher serves its Invitrogen user guides from three hosts and only one of them answers
a script. Every vendor fact below came from a PDF fetched on that host.

| Route | Works | Used for |
| --- | --- | --- |
| `documents.thermofisher.com/TFS-Assets/LSG/manuals/*.pdf` | yes (`curl`) | every Invitrogen user guide cited here |
| `assets.thermofisher.com/TFS-Assets/LSG/manuals/gatewayman.pdf` | yes (`curl`) | the same Gateway Technology guide, second copy |
| `tools.thermofisher.com/content/sfs/manuals/*.pdf` | no — HTTP 403 | — |
| `www.thermofisher.com` HTML pages | yes (200) | catalogue and technical-resource pages |
| `patents.google.com/patent/<number>/en` | yes (`curl`) | the att sequences and the primer experiments |
| `eutils.ncbi.nlm.nih.gov/entrez/eutils/*` | yes (`curl`) | the phage lambda genome, J02459, and deposited plasmids |
| `genome.cshlp.org/content/*/*.full.pdf` | no — HTTP 403 | — |
| `pmc.ncbi.nlm.nih.gov/articles/*/pdf/*` | no — returns an HTML interstitial | — |

Guessing a manual's filename is how the set above was found: the host answers a wrong name
with HTTP 403 and an XML body, and a right one with `application/pdf`. A web search for
`documents.thermofisher.com/TFS-Assets/LSG/manuals` plus the product name returns real
filenames.

**The vendor's sequences are in pictures, not in text.** The recombination-region diagrams
that carry the att sequences are raster images inside the PDFs, so `pdftotext` returns the
captions and none of the bases. They were read by rendering the page with
`pdftoppm -r 600 -png -f <page> -l <page>`, cropping the diagram band with ImageMagick and
reading the crop.

**Two warnings about reading a sequence off Google Patents.** The page renders the description
twice, once entity-annotated and once plain, and the two copies do not always agree — the
annotated copy of US 7,670,823 B1 prints attR1 with eight A's where the plain copy prints six.
And the specification wraps a long sequence with hyphens, which are line breaks and not bases.
Both traps were resolved here by arithmetic rather than by picking a copy: the same run of A's
appears inside the patent's own attP2, and the whole 233 bp attP lands on λ J02459
27586–27818 with four differences. Anything that did not close that way is in
[Open gaps](#19-open-gaps).

Anything still unverified is in [Open gaps](#19-open-gaps) and is never guessed.

## 1. Sources and licences

| Document or dataset | Licence | Can we ship it? |
| --- | --- | --- |
| US 7,670,823 B1, FIG. 9 and SEQ ID NO:1–8 (the att sequences) | no 37 CFR 1.71(e) copyright notice anywhere in the document; USPTO marks patent full text CC PDM 1.0 | **Yes** — the sequence strings, with the patent cited; not blocks of its prose or its figures |
| NCBI GenBank records, including J02459 and the deposited Gateway plasmids | NCBI places no restrictions on the data | **Yes**, with the accession cited |
| Salim, Feger & Busso 2016, *Data in Brief* 9:946–955 (attB1 and attB2 only) | CC BY 4.0 | **Yes**, with attribution |
| Invitrogen / Thermo Fisher user guides (every PDF cited here) | © Life Technologies / Thermo Fisher, all rights reserved | **No** — cite and paraphrase, never redistribute the files |
| Thermo Fisher catalogue pages (the prices in [17](#17-what-it-costs-and-how-long-it-takes)) | same | **No** — cite the figure and the date, never mirror the page |
| Hartley, Temple & Brasch 2000, *Genome Research* 10:1788–1795 | © 2000 Cold Spring Harbor Laboratory Press | **No** — cite and quote, do not redistribute |
| Addgene, "Plasmids 101: Gateway Cloning", and Addgene's sequence records | all rights reserved, non-commercial informational use only | **No** — cite and paraphrase; no figures, no sequences |
| SnapGene, "Gateway Cloning Technique", and its plasmid set | © 2026 SnapGene; commercial entities must ask permission | **No** — cite by URL and paraphrase |

Licence text, quoted:

- USPTO, *Terms of use for USPTO websites*: "Patents are published as part of the terms of
  granting the patent to the inventor. Subject to limited exceptions reflected in 37 CFR
  1.71(d) & (e) and 1.84(s), the text and drawings of a patent are typically not subject to
  copyright restrictions." The same page cautions that part of a patent's text or drawings may
  be under copyright, which is why the exception was checked: the 37 CFR 1.71(e)
  authorization paragraph is **absent** from US 7,670,823 B1 and from the six sibling patents
  read for this note, and no `©` notice appears in any of them. USPTO's own open-data catalogue
  marks "Patent Grant Full Text (1976 – Present)" with
  `"license": "https://creativecommons.org/publicdomain/mark/1.0"`.
- NCBI Data Usage Policy, *Molecular Data Usage*: "Therefore, NCBI itself places no
  restrictions on the use or distribution of the data contained therein. Nor do we accept data
  when the submitter has requested restrictions on reuse or redistribution." The same
  paragraph adds that a submitter may still claim rights and that NCBI cannot assess such a
  claim — a disclaimer of warranty rather than a restriction.
- Salim, Feger & Busso 2016, permissions block: "© 2016 The Authors […] This is an open access
  article under the CC BY license (`http://creativecommons.org/licenses/by/4.0/`)."
- Hartley 2000, page footer: "10:1788–1795 ©2000 by Cold Spring Harbor Laboratory Press
  ISSN 1088-9051/00 $5.00".
- Addgene Terms of Use, last updated 24 January 2023, *IP Restrictions*: "Addgene and its
  contributing third parties retain all copyrights and/or applicable rights to all text, data,
  graphics, sounds, and any other media provided on the Site. Content may not be reproduced,
  duplicated, copied, sold, resold, or otherwise exploited for any commercial purpose without
  the express written consent of Addgene". Its sequence pages now also require a login.
- SnapGene plasmid-set notice: "This material may be used without restriction by academic,
  nonprofit, and governmental entities […] Commercial entities must contact GraphPad Software,
  LLC for permission and terms of use." A licence that discriminates by who the licensee is
  cannot be sublicensed under MIT, which grants everyone the same rights.
- Thermo Fisher, *Website and Mobile Application Terms of Use*, effective 17 June 2016, §4:
  "All intellectual property rights are reserved unless granted in an express written
  license." The page grants only a licence to access the site; there is no grant to copy or
  redistribute, so the manuals stay in git-ignored `reference_docs/`.

**Two limits on the "yes" for the att sequences.** First, its footing is that a DNA sequence
is a fact, and no source read here addresses DNA sequences and copyright directly — the step
that *is* citable is the general one, *Feist v. Rural Telephone*, 499 U.S. 340 (1991): "That
there can be no valid copyright in facts is universally understood." Second, copyright is not
the only right a patent carries; the Limited Use Label License printed on every Clonase sheet
restricts the **product**, and says "The buyer cannot modify the recombination sequence(s)
contained in this product for any purpose." Printing a sequence and practising an invention
are different questions, and only the first is settled here.

## 2. The att sites, base for base

Each att site has a **25 bp recombination region** — the part that both partners share and
that the reaction rewrites. `attB` is that region and nothing else. `attP`, `attL` and `attR`
are the same 25 bp with arms attached; the lengths are the vendor's:

> | Site | Length | Found in… |
> | --- | --- | --- |
> | attB | 25 bp | Expression vector, Expression clone |
> | attP | 200 bp | Donor vector |
> | attL | 100 bp | Entry vector, Entry clone |
> | attR | 125 bp | Destination vector |

(*Gateway Technology with Clonase II* user guide, MAN0000470, part number 25-0749, revision
date 2 April 2012, page 5, "Characteristics of the Modified att Sites".)

**The orientation convention this note uses**, because primary sources print site 2 both ways
and label both "attB2": every site below is written 5'→3' **in its own orientation**, so that
site 1 and site 2 line up base for base. Under that convention the two sites of a molecule
point at each other across the insert, which is how the vendor's diagrams draw them. US
7,670,823 B1 prints attB2 as `ACCCAGCTTTCTTGTACAAAGTGGT` in FIG. 9 and as
`GGGGACCACTTTGTACAAGAAAGCTGGGT` in its Example 20 primer; the two are reverse complements of
each other, and the second is the one that matches this table.

The 25 bp recombination region of all eight sites:

| Site | 25 bp recombination region | Total length | Arms |
| --- | --- | --- | --- |
| attB1 | `ACAAGTTTGTACAAAAAAGCAGGCT` | 25 | none |
| attB2 | `ACCACTTTGTACAAGAAAGCTGGGT` | 25 | none |
| attP1 | `CCAACTTTGTACAAAAAAGCTGAAC` | 233 | 75 bp before, 133 bp after |
| attP2 | `CCAACTTTGTACAAGAAAGCTGAAC` | 233 | 75 bp before, 133 bp after |
| attL1 | `CCAACTTTGTACAAAAAAGCAGGCT` | 100 | 75 bp before |
| attL2 | `CCAACTTTGTACAAGAAAGCTGGGT` | 100 | 75 bp before |
| attR1 | `ACAAGTTTGTACAAAAAAGCTGAAC` | 125 | 100 bp after |
| attR2 | `ACCACTTTGTACAAGAAAGCTGAAC` | 125 | 100 bp after |

And the four long sites in full, from US 7,670,823 B1, FIG. 9, "Recombination Site Nucleotide
Sequences", restated in the specification as SEQ ID NO:3–8 and re-oriented to the convention
above:

```text
attP1 (233 bp, SEQ ID NO:3, reverse-complemented from FIG. 9)
CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTTTTATAATG
CCAACTTTGTACAAAAAAGCTGAAC
GAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTA
AAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA

attP2 (233 bp, SEQ ID NO:4)
CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTCTTATAATG
CCAACTTTGTACAAGAAAGCTGAAC
GAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTA
AAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA

attL1 (100 bp, SEQ ID NO:7)      attL2 (100 bp, SEQ ID NO:8)
= attP1[1..75] + CCAACTTT        = attP2[1..75] + CCAACTTT
  + GTACAAAAAAGCAGGCT              + GTACAAGAAAGCTGGGT

attR1 (125 bp, SEQ ID NO:5)      attR2 (125 bp, from SEQ ID NO:6)
= ACAAGTTT + GTACAAAAAAGCTGAAC   = ACCACTTT + GTACAAGAAAGCTGAAC
  + attP1[101..200]                + attP2[101..200]
```

**Every one of those was checked twice before it was written down.**

1. *Against the vendor's own diagrams.* The BP and LR recombination-region figures in the
   2003 *Gateway Technology* user guide, MAN0000282 pages 6 and 7, spell the 25 bp region of
   all eight sites (and abbreviate the arms as `N75` and `N100`). All eight agree with the
   patent, base for base.
2. *Against the phage.* Hartley 2000's Methods says where attP comes from: "The attP sites in
   pDONR203 correspond to lambda coordinates 27586 through 27818." Fetching J02459 and
   aligning puts the patent's 233 bp attP exactly on **27586–27818** and leaves only four
   differences for attP1 and five for attP2 — listed in
   [3](#3-what-an-att-site-is-made-of). One of them is λ 27630 C→G, which is precisely what
   the same paragraph of Hartley says was done: "base 27630 has been changed to a G to remove
   an NdeI site". That a base the patent prints lands on a coordinate a different paper names
   is the strongest single check in this note.

For **attB1 and attB2** there is also a third, openly licensed source: Salim, Feger & Busso
(2016), *Data in Brief* 9:946–955, CC BY 4.0, which prints them in the same orientation as
this table — "the 25 nucleotide att B1 sequence (ACAAGTTTGTACAAAAAAGCAGGCT)" and "the 25
nucleotide att B2 sequence (ACCACTTTGTACAAGAAAGCTGGGT)". The vendor prints the same bases
inside its primer diagrams (MAN0000470 pages 13 and 14), and Hartley 2000's Methods gives the
adapters in full: "[A] = ggggacaagtttgtacaaaaaagcaggctcattta-actttaagaaggagatatatacc, in
which sequences in bold type correspond to attB1"; "[B] = ggggaccactttgtacaagaaagctgggt".

**Four discrepancies to carry, none of them guesses.**

- **attP is 200 bp or 233 bp depending on which document you read.** The vendor's length table
  says 200; the patent, Hartley's λ coordinates and the arithmetic all say 233. They are
  reconcilable — 200 is the part of attP that takes part in the reaction and comes back in the
  LR by-product, and the extra 33 is the piece attR drops
  ([6](#6-what-lr-writes)) — but a plan should measure the record it is handed rather than
  assume either. attL at 100 and attR at 125 are not in dispute.
- **FIG. 9's attR2 is 135 bp, not 125**, because it carries ten extra bases of the vector
  around it (`GTCGACCTGC`, a SalI site) past the end of the site. The table above trims them;
  the 125 bp form is what the vendor's length table and the attR1 entry both agree with.
- **attL1 and attP1 disagree at position 55** inside the patent itself: attP1 has `G` there
  and attL1 has `A`, although attL1's first 83 bases are attP1's. λ has `A`, and deposited
  GenBank attL1 features are mostly `A` with a minority `G`. The vendor's own note is probably
  the explanation: "site-specific point mutations have been made to some att sites to increase
  recombination efficiency. As a result, sequence variations may exist among the att sites.
  For example, the pDONR201 attP1 sequence varies slightly from the pDONR221 attP1 sequence."
  (MAN0000470 page 5.) **This position is not invariant, and a plan must not match on it.**
- **An older patent prints a different attB2.** US 6,143,557 A gives
  `AGCCTGCTTTCTTGTACAAACTTGT`, which shares the 15 bp core with the site above but carries
  wild-type flanks. US 7,670,823's is the one a real expression clone carries.

**Which donor vector.** The diagrams are drawn for pDONR221 and pDONR/Zeo, and the manual
warns that another donor vector may differ: "If you are performing a BP recombination
reaction using a donor vector other than pDONR221 or pDONR/Zeo, note that the sequences of
the recombination regions may vary slightly but the mechanism of recombination remains the
same." (MAN0000470 page 6.)

**Verdict on shipping: yes, and it changes almost nothing.** [1](#1-sources-and-licences)
gives the footing — the patent carries no copyright notice, and every one of the eight strings
also occurs verbatim in deposited GenBank plasmids, which NCBI places no restrictions on. So
a test fixture may be built from these bases and a check may look for them. But a plan should
still **find** its att sites in the user's own records rather than assert them, because the
two discrepancies above are the vendor telling us the sequences drift.

## 3. What an att site is made of

The 25 bp region is **5 bp + 15 bp + 5 bp**, and the vendor names those parts:

> Mutations have been introduced into the short (5 bp) regions flanking the 15-bp core
> regions of the attB sites to minimize secondary structure formation in single-stranded
> forms of attB plasmids (e.g. phagemid ssDNA or mRNA).
>
> Mutations have been made to the core regions of the att sites to eliminate stop codons and
> to ensure specificity of the recombination reactions to maintain orientation and reading
> frame.
>
> A 43 bp portion of the attR site has been removed to make the in vitro attL x attR reaction
> irreversible and more efficient (Bushman et al., 1985).

(All three, MAN0000470 page 5, "Modifications to the att Sites".)

The patent names the same parts and says what each one does. US 7,670,823 B1, Detailed
Description:

> …the 15 bp core region (GCTTTTTTATACTAA) (SEQ ID NO:9) which is identical in all four
> wildtype lambda att sites, attB, attP, attL and attR […] Analogously, the core regions in
> attB1, attP1, attL1 and attR1 are identical to one another, as are the core regions in
> attB2, attP2, attL2 and attR2.
>
> …the seven bp overlap region (TTTATAC, which is defined by the cut sites for the integrase
> protein and is the region where strand exchange takes place) that occurs within this 15 bp
> core region (GCTTT TTTATAC TAA).
>
> …mutants […] in which substitutions have been made within the first three positions of the
> seven bp overlap (TTTATAC) have been found […] to strongly affect the specificity of
> recombination, mutant nucleic acid molecules in which substitutions have been made in the
> last four positions (TTTATAC) only partially alter recombination specificity, and mutant
> nucleic acid molecules comprising nucleotide substitutions outside of the seven bp overlap,
> but elsewhere within the 15 bp core region, do not affect specificity of recombination but
> do influence the efficiency of recombination.

Lining the real sites up against J02459 shows exactly which bases were changed. The phage's
orientation, with the Gateway rows reverse-complemented from the table in
[2](#2-the-att-sites-base-for-base):

```text
             1    5    10   15   20   25
lambda wt    GTTCAGCTTTTTTATACTAAGTTGG
attP1        GTTCAGCTTTTTTGTACAAAGTTGG
attB1        AGCCTGCTTTTTTGTACAAACTTGT
attP2        GTTCAGCTTTCTTGTACAAAGTTGG
attB2        ACCCAGCTTTCTTGTACAAAGTGGT
             [flnk][     core 15    ][flnk]
                        [overlap]
```

Three things fall out, and each is checked rather than assumed:

1. **attB1 and attP1 have identical 15 bp cores** (`GCTTTTTTGTACAAA`), and so do attB2 and
   attP2 (`GCTTTCTTGTACAAA`) — which is the patent's sentence above, confirmed on the
   sequences. What tells attB from attP is the 5 bp flanks, which is what the manual says was
   mutated and why an attB site is 25 bp and an attP site is 233.
2. **The 7 bp overlap** sits at positions 11–17 of the 25 bp region in the phage's
   orientation. Wild type is `TTTATAC`; att1 is `TTTGTAC`; att2 is `CTTGTAC`. They differ at
   the **first** position of the overlap — the patent's "strongly affect the specificity"
   window. Written in this note's orientation the overlap is positions 9–15 of the 25 bp
   region: `GTACAAA` for att1 and `GTACAAG` for att2.
3. **Strand exchange happens across that overlap.** It is why the BP and LR diagrams shade the
   two strands seven bases out of step: on the top strand the product switches source between
   position 8 and position 9, and on the bottom strand between 15 and 16 — the two edges of
   the overlap.

Across the whole 233 bp, the differences from λ 27586–27818 are these and no others:

| λ coordinate | λ base | attP1 | attP2 | What it is |
| --- | --- | --- | --- | --- |
| 27630 | C | G | G | removes an NdeI site (Hartley 2000, Methods) |
| 27729 | T | — | C | the specificity base: first position of the 7 bp overlap |
| 27732 | A | G | G | inside the 7 bp overlap |
| 27736 | T | A | A | inside the 15 bp core, outside the overlap |
| 27752 | A | — | G | in the arm, outside the core |
| 27764 | T | C | — | in the arm, outside the core — position 55 of attP1 |

The last row is the one the patent contradicts itself on, discussed in
[2](#2-the-att-sites-base-for-base).

**The attR deletion.** The manual says "A 43 bp portion of the attR site has been removed to
make the in vitro attL x attR reaction irreversible and more efficient". The sequences say
**33**: attR is attP's arm truncated after 100 of its 133 bases, and what is dropped is
`AATCAACTACTTAGATGGTATTAGTGACCTGTA`, which is 33 bp. Hartley agrees with the sequences, not
the manual — "The attR sites in the Destination Vectors lack the bases between 27586 and
27618", a span of 33 bases. Use 33.

## 4. Which site pairs with which

The vendor states the pairing and nothing more:

> - attB1 sites react only with attP1 sites
> - attB2 sites react only with attP2 sites
> - attL1 sites react only with attR1 sites
> - attL2 sites react only with attR2 sites

(MAN0000470 page 5.)

So there are exactly two reactions, and each is a pair of recombinations:

| Reaction | Substrates | Products | Proteins |
| --- | --- | --- | --- |
| BP | attB1 × attP1 and attB2 × attP2 | attL1 … attL2 (entry clone) and attR1 … attR2 (by-product) | Int, IHF |
| LR | attL1 × attR1 and attL2 × attR2 | attB1 … attB2 (expression clone) and attP1 … attP2 (by-product) | Int, IHF, Xis |

The proteins are from Hartley 2000: "This system carries out two reactions: (1) attB ×
attP → attL + attR mediated by the integrase (Int) and integration host factor (IHF)
proteins and (2) attL × attR → attB + attP mediated by Int, IHF, and excisionase (Xis). Thus,
the direction of the reactions is controlled by providing different combinations of proteins
and sites."

**Orientation.** Site 1 and site 2 are different sites, and that is what fixes the insert's
direction. Hartley 2000: "Note that attB1 will recombine with attP1 but not attP2, thereby
maintaining orientation of the DNA segment during recombination. Aberrant recombination
events have not been identified in hundreds of sequenced RC clones." The paper's abstract
puts the consequence plainly: "The resulting subclones maintain orientation and reading frame
register, allowing amino- and carboxy-terminal translation fusions to be generated."

**A site the wrong way round.** No source read for this note reports what happens when an att
site sits in the reverse orientation in a Gateway substrate — see
[Open gaps](#19-open-gaps). What the sources do support is the shape of the check a plan
should make: a site pairs only with its own number, and the two sites of a molecule must
point the same way relative to the segment between them, because the segment is what moves.

## 5. What BP writes

Read off the BP diagram, MAN0000282 page 6, drawn for an attB-PCR product and pDONR221 or
pDONR/Zeo. The manual's own key: "Shaded regions correspond to those sequences transferred
from the attB-PCR product into the entry clone following recombination. Note that the attL
sites are composed of sequences from attB and attP."

Substrates, top strand as the page prints them:

```text
attB-PCR product   GGGGACAAGTTTGTACAAAAAAGCAGGCT--[insert]--ACCCAGCTTTCTTGTACAAAGTGGTCCCC
pDONR221 / pDONR/Zeo
                   ---N75-CCAACTTTGTACAAAAAAGCTGAAC-N100--[ccdB-CmR]--N100-GTTCAGCTTTCTTGTACAAAGTTGG-N75---
```

Products:

```text
entry clone        ---N75-CCAACTTTGTACAAAAAAGCAGGCT--[insert]--ACCCAGCTTTCTTGTACAAAGTTGG-N75---
by-product         GGGGACAAGTTTGTACAAAAAAGCTGAAC-N100--[ccdB-CmR]--N100-GTTCAGCTTTCTTGTACAAAGTGGTCCCC
```

Base for base, in each site's own orientation, writing the 25 bp region as
`[1..8] + [9..25]`:

| Junction | Spells | 1..8 from | 9..25 from |
| --- | --- | --- | --- |
| attL1 | `CCAACTTT` + `GTACAAAAAAGCAGGCT` | attP1 | attB1 |
| attL2 | `CCAACTTT` + `GTACAAGAAAGCTGGGT` | attP2 | attB2 |
| attR1 | `ACAAGTTT` + `GTACAAAAAAGCTGAAC` | attB1 | attP1 |
| attR2 | `ACCACTTT` + `GTACAAGAAAGCTGAAC` | attB2 | attP2 |

That is the whole arithmetic: the boundary on the top strand falls between position 8 and
position 9 of the 25 bp region, which is the 5' edge of the 7 bp overlap
([3](#3-what-an-att-site-is-made-of)). attL keeps the attP arm that was upstream of the
crossover — 75 bp — giving 75 + 8 + 17 = 100 bp. attR keeps the first 100 bp of the 133 bp arm
that was downstream, giving 8 + 17 + 100 = 125 bp, and drops the last 33
([3](#3-what-an-att-site-is-made-of)). Both totals match the vendor's length table. The
diagram's `N100` is that 100 bp on the attR side and the full 133 on the attP side; only the
attP figure is drawn short.

**What the insert keeps.** Positions 9–25 of each attB site, so an entry clone carries the
insert with `GTACAAAAAAGCAGGCT` on one side and `GTACAAGAAAGCTGGGT` on the other, each
preceded by attP's `CCAACTTT`. The four G residues of the PCR primer do **not** enter the
entry clone: the diagram puts them in the by-product.

## 6. What LR writes

Read off the LR diagram, MAN0000282 page 7, drawn for a pENTR/D-TOPO entry clone and the
pcDNA6.2/V5-DEST destination vector. The manual's key: "Shaded regions correspond to those
sequences transferred from the pENTR/D-TOPO entry clone into the expression clone following
recombination. Note that the attB sites are composed of sequences from attL and attR sites."

Substrates and products, top strand:

```text
entry clone        ---N75-CCAACTTTGTACAAAAAAGCAGGCT--[gene]--ACCCAGCTTTCTTGTACAAAGTTGG-N75---
destination vector ---ACAAGTTTGTACAAAAAAGCTGAAC-N100--[ccdB-CmR]--N100-GTTCAGCTTTCTTGTACAAAGTGGT---

expression clone   ---ACAAGTTTGTACAAAAAAGCAGGCT--[gene]--ACCCAGCTTTCTTGTACAAAGTGGT---
by-product         ---N75-CCAACTTTGTACAAAAAAGCTGAAC-N100--[ccdB-CmR]--N100-GTTCAGCTTTCTTGTACAAAGTTGG-N75---
```

Base for base:

| Junction | Spells | 1..8 from | 9..25 from |
| --- | --- | --- | --- |
| attB1 | `ACAAGTTT` + `GTACAAAAAAGCAGGCT` | attR1 | attL1 |
| attB2 | `ACCACTTT` + `GTACAAGAAAGCTGGGT` | attR2 | attL2 |
| attP1 | `CCAACTTT` + `GTACAAAAAAGCTGAAC` | attL1 | attR1 |
| attP2 | `CCAACTTT` + `GTACAAGAAAGCTGAAC` | attL2 | attR2 |

LR gives back exactly the attB1 and attB2 that the PCR primers carried in. Both reactions use
one rule and it runs in both directions: **positions 1–8 of the 25 bp region come from one
partner and positions 9–25 from the other**. The expression clone therefore carries a full
25 bp attB site at each end of the insert — this junction is a cloning scar, not a scarless
one, and the next section is what it costs a fusion.

**What LR does not give back is a whole attP.** The by-product's attP is attL's 83 bases plus
attR's 117, so it is 200 bp where the donor vector's was 233 — the 33 bp attR deletion
([3](#3-what-an-att-site-is-made-of)) is missing from it. That is the manual's stated purpose
for the deletion, "to make the in vitro attL x attR reaction irreversible and more efficient",
and it is also the likeliest reason its length table says attP is 200 bp: 200 is the attP that
takes part and comes back, and 233 is the attP the patent prints and λ carries. This note does
not decide which number a given donor vector holds; a plan should measure the record it is
given.

## 7. The attB primer tail

The forward primer, quoted from MAN0000470 page 13:

> 1. Four guanine (G) residues at the 5′ end, followed by:
> 2. The 25-bp attB1 site, followed by:
> 3. At least 18–25 bp of template- or gene-specific sequences.

and the reverse primer, page 14, the same three items with "The 25-bp attB2 site" in place of
attB1 and "18–25 bp" in place of "at least 18–25 bp".

So a tail is **29 bases** — `GGGG` plus the 25 bp site — before the annealing region, and 30
or 31 where a fusion needs a frame base ([8](#8-reading-frame-across-a-junction)). Written
out:

| Tail | Bases | Length |
| --- | --- | --- |
| attB1, no fusion | `GGGGACAAGTTTGTACAAAAAAGCAGGCT` | 29 |
| attB1, N-terminal fusion | the same, then two more bases | 31 |
| attB2, no fusion | `GGGGACCACTTTGTACAAGAAAGCTGGGT` | 29 |
| attB2, C-terminal fusion | the same, then one more base | 30 |

Hartley 2000 gives the same figure in one phrase: "PCR products flanked by attB sites can be
generated by incorporating attB sites (25 base + 4 G residues) at the 5′ end of PCR primers."

**Why the four G residues.** Neither vendor guide gives a reason; the patent gives the
measurement behind the rule. US 7,670,823 B1, Example 9, Results:

> In initial experiments, primers for amplifying tetR and ampR from pBR322 were constructed
> containing only the tetR- or ampR-specific targeting sequences, the targeting sequences plus
> attB1 (for forward primers) or attB2 (for reverse primers) sequences shown in FIG. 9, or
> the attB1 or attB2 sequences with a 5′ tail of four guanines. […] These results demonstrated
> that primers containing attB sequences provided for a somewhat higher number of colonies on
> the tetracycline and ampicillin plates. However, inclusion of the 5′ extensions of four or
> five guanines on the primers in addition to the attB sequences provided significantly better
> cloning results […] These results indicate that the optimal primers for cloning of PCR
> products using recombinational cloning will contain the recombination site sequences with a
> 5′ extension of four or five guanine bases.

Two things to carry from that. The claim is **four or five**, not four; and the numbers behind
"significantly better" are in FIGS. 66 and 67, which are figure images this note did not
read — so the direction is sourced and the size of the effect is not. The extension itself
does not reach the clone: the BP diagram puts the four G residues in the by-product.

**The shortest tail reported to work.** A whole attB tail is not the only way to build the
product. Where the full primer would run past 70 bp, the vendor splits it in two:

> We recommend using this protocol to produce attB-PCR products if your PCR primers are
> greater than 70 bp. To use this protocol, you will need to have 2 sets of PCR primers, one
> set for the gene-specific amplification and a second set to install the complete attB
> sequences (adapter-primers attB1 and attB2).
>
> Include 12 bases of the attB1 or attB2 site on the 5′ end of each primer, as appropriate.
>
> - attB1 forward: 5′-AA AAA GCA GGC TNN - template-specific sequences-3′
> - attB2 reverse: 5′-A GAA AGC TGG GTN - template-specific sequences-3′
>
> - attB1 adapter: 5′-G GGG ACA AGT TTG TAC AAA AAA GCA GGC T -3′
> - attB2 adapter: 5′-GGG GAC CAC TTT GTA CAA GAA AGC TGG GT -3′

(MAN0000470 pages 47–48, "Preparing attB-PCR Products Using attB Adapter PCR".) So the
shortest **vendor-published** gene-specific tail is 12 bases of the att site plus the frame
base, with the rest installed by a universal adapter in a second round of PCR. The product
going into BP still carries the full 29.

The patent tested the whole range and the answer depends on how the primers are balanced. US
7,670,823 B1, Example 20, "PCR Cloning Using Universal Adapter-Primers", tried "overlaps of
various lengths from 6 bp to 18 bp" against a 256 bp hemoglobin target and concluded that
"gene-specific primers with overlaps of 10 bp to 18 bp can be used successfully". The stated
fall-off:

| Overlap | Gene-specific : adapter primer | Result, quoted |
| --- | --- | --- |
| 18 bp and 15 bp | 10 : 10 pmol | "successfully amplify predominately full-length PCR product" |
| 12 bp or less | 10 : 10 pmol | "smaller intermediate products containing one or no universal attB adapter predominated the reactions" |
| 12 bp | 3 : 30 pmol | "the overlap necessary to obtain predominately full-length PCR product was reduced to 12 bp" |
| 11 bp | 2 : 40 pmol | "generated predominately full-length PCR products with gene-specific primers containing an 11 bp overlap" |

So **11 bp is the shortest overlap reported to work, and only at 2 pmol against 40** — the
figure to quote is not a flat "10 bp", it is an overlap paired with a stoichiometry. Note also
that the overlap is a property of the *gene-specific* primer; every product that reaches BP
still carries a complete 29-base attB tail installed by the adapter.

## 8. Reading frame across a junction

Hartley 2000, Figure 1D, prints the junction as codons, which settles the arithmetic in one
picture:

```text
        |<------------ attB1 ------------>|          |<------------ attB2 ------------>|
   ...  aca agt ttg tac aaa aaa gca ggc tnn  [gene]  nac cca gct ttc ttg tac aaa gtg gtn  ...
         T   S   L   Y   K   K   A   G   -            -   P   A   F   L   Y   K   V   V
```

Reading it, with the vendor's rules quoted beside it:

**At the N terminus, attB1 contributes 25 bases and a fusion adds two more**, making 27 — nine
codons, of which the first eight are fixed: Thr-Ser-Leu-Tyr-Lys-Lys-Ala-Gly. The ninth is
attB1's final `T` plus the two added bases, and the vendor constrains them:

> The attB1 site ends with a thymidine (T). If you wish to fuse your PCR product in frame with
> an N-terminal tag, the primer must include two additional nucleotides to maintain the proper
> reading frame with the attB1 region. These two nucleotides cannot be AA, AG, or GA, because
> these additions will create a translation termination codon.

(MAN0000470 page 13.) So an N-terminal fusion is in frame when the destination vector's ATG
places attB1's first base at a codon boundary, and the two added bases are anything but `AA`,
`AG` or `GA`. The vendor's cross-check is the same frame stated from the vector's side: "Keep
the -AAA-AAA- triplets in the attR1 site in frame with the translation reading frame of the
fusion protein." (MAN0000470 page 14.)

**At the C terminus, attB2 contributes 25 bases on the coding strand, a fusion adds one
before them and the vector supplies one after**, making 27 again — nine codons, of which
seven are fixed: Pro-Ala-Phe-Leu-Tyr-Lys-Val. The first codon is the added base plus attB2's
`ac`; the last is attB2's `gt` plus the vector's next base, which is Val whatever that base
is. The rules:

> The primer must include one additional nucleotide to maintain the proper reading frame with
> the attB2 region.
>
> Any in-frame stop codons between the attB2 site and your gene of interest must be removed.
>
> If you do not intend to fuse your PCR product in frame with a C-terminal tag, your gene of
> interest or the primer must include a stop codon.

(MAN0000470 page 14.) And the vector-side cross-check: "Keep the -TTT-GTA (TAC-AAA on the
complementary strand) triplets in the attR2 site in frame with the translation reading frame
of the fusion protein." (MAN0000470 page 15.)

**So the stop codon is the whole of the C-terminal decision.** A C-terminal fusion must have
no stop between the gene and attB2; anything else must have one, and the gene or the primer
supplies it. The manual is explicit that this is what the core mutations were for: they were
made "to eliminate stop codons" ([3](#3-what-an-att-site-is-made-of)), which is why the eight
and seven fixed codons above read as amino acids at all.

## 9. ccdB counter-selection

**What it does.** The vendor states the target and cites the paper:

> The CcdB protein interferes with E. coli DNA gyrase (Bernard and Couturier, 1992), thereby
> inhibiting growth of most E. coli strains (e.g. DH5α, TOP10). When recombination occurs […]
> the ccdB gene is replaced by the gene of interest. Cells that take up unreacted vectors
> carrying the ccdB gene or by-product molecules retaining the ccdB gene will fail to grow.
> This allows high-efficiency recovery of the desired clones.

(MAN0000282 page 8, "ccdB Gene"; the same paragraph, shortened, in MAN0000470 page 8.) The
mechanism itself is in the paper, not the manual — Bernard & Couturier (1992), *J. Mol. Biol.*
226:735–745, abstract: "These results strongly suggest that the CcdB protein, like quinolone
antibiotics and a variety of antitumoral drugs, is a DNA topoisomerase II poison." The same
abstract is where `gyrA462` comes from: "This mutation is located in the gene encoding the A
subunit of topoisomerase II and produces an Arg462----Cys substitution in the amino acid
sequence of the GyrA polypeptide. Hence, the mutation was called gyrA462."

**Which strain grows a ccdB vector.** Three have been sold, and they do not all work the same
way:

| Strain | Genotype, as printed | How it resists |
| --- | --- | --- |
| DB3.1 (discontinued) | `F- gyrA462 endA1 ∆(sr1-recA) mcrB mrr hsdS20(rB-, mB-) supE44 ara14 galK2 lacY1 proA2 rpsL20(Smr) xyl5 ∆leu mtl1` | `gyrA462` — the gyrase no longer binds CcdB |
| ccdB Survival T1R, C7510-03 (discontinued) | `F- mcrA ∆(mrr-hsdRMS-mcrBC) Φ80lacZ∆M15 ∆lacX74 recA1 ara∆139 ∆(ara-leu)7697 galU galK rpsL (StrR) endA1 nupG tonA::Ptrc-ccdA` | the `ccdA` antidote, expressed from the chromosome |
| ccdB Survival 2 T1R, A10460 (current) | `F- mcrA Δ(mrr-hsdRMS-mcrBC) Φ80lacZΔM15 ΔlacX74 recA1 araΔ139 Δ(ara-leu)7697 galU galK rpsL (StrR) endA1 nupG fhuA::IS2` | **not stated** — see [Open gaps](#19-open-gaps) |

(DB3.1 from MAN0000282 page 18; C7510-03 from its product sheet, part C751003.pps, revision
11 June 2004; A10460 from MAN0000761 revision 2.0.) The current strain's printed genotype
carries neither `gyrA462` nor a `ccdA` cassette, and no manual says what replaced them.

**The rule, and what it is for.** MAN0000470 page 36:

> To propagate and maintain your destination vector, you must use ccdB Survival T1R E. coli.
> […] Note: Do not use general E. coli cloning strains including OmniMAX 2-T1R, TOP10 or DH5α
> for propagation and maintenance as these strains are sensitive to CcdB effects.

**Which strain must never be used to select the clone.** The prohibition the manuals state is
about F′, not about ccdB resistance:

> Do not use E. coli strains that contain the F′ episome (e.g. TOP10F′) for transformation.
> These strains contain the ccdA gene and will prevent negative selection with the ccdB gene.

(MAN0000470 page 24, and again on page 31 for LR.) What happens if you do is a troubleshooting
row: "LR Reaction: High background in the absence of the entry clone", cause "LR reaction
transformed into an E. coli strain containing the F′ episome and the ccdA gene" (MAN0000470
page 42).

**No manual forbids plating a BP or LR reaction on a ccdB-resistant strain.** Thermo's support
page frames it the other way round: "growing non-recombined vector requires special cells (One
Shot ccdB Survival 2 T1R Competent Cells) […] On the other hand, general E. coli cloning
strains including TOP10 or DH5a may be used for plating the BP or LR reaction". That the
counter-selection would be lost on a resistant host follows from the mechanism, but it is an
inference rather than a quote, and a check that makes it should say so.

**The check that proves ccdB still works.** MAN0000470 page 38: transform 10–50 pg of the
destination vector into each of two strains, and "The destination vector should give 10,000
times more colonies in ccdB Survival T1R cells than in OmniMAX 2-T1R cells. Any ratio less
than 10,000:1 indicates either an inactive ccdB gene or contamination of the plasmid prep with
another antibiotic-resistant plasmid."

## 10. The BP reaction

Two generations, and they are not interchangeable. Clonase II pre-mixes the buffer and halves
the volume; the vendor says so: "The Gateway BP Clonase II enzyme mix combines the proprietary
enzyme formulation and 5X BP Clonase Reaction Buffer previously supplied as separate
components […] into an optimized single tube format" and "Do not use the protocol for Gateway
BP Clonase II enzyme mix" with the older product (MAN0000470 page 22).

| Component | BP Clonase II, 10 μL | BP Clonase, 20 μL |
| --- | --- | --- |
| attB-PCR product or linearised attB expression clone | 1–7 μL, 20–50 fmol (~15–150 ng) | 1–10 μL, 40–100 fmol (~30–300 ng) |
| Donor vector, 150 ng/μL | 1 μL | 2 μL |
| 5X BP Clonase Reaction Buffer | not supplied — it is in the mix | 4 μL |
| TE buffer, pH 8.0 | to 8 μL | to 16 μL |
| Enzyme mix | 2 μL | 4 μL |
| Incubation | 25 °C, 1 hour | 25 °C, 60 minutes |
| Stop | 1 μL proteinase K (2 μg/μL), 37 °C, 10 min | 2 μL proteinase K, 37 °C, 10 min |

(Clonase II from product sheet 11789.II.pps, revision 31 October 2010, page 3, and MAN0000470
page 23; Clonase I from 11789.pps, revision 10/31/10, page 3. The controls are on the same
MAN0000470 table: the positive control replaces the substrate with 2 μL of pEXP7-tet at
50 ng/μL, and the negative control is the sample reaction with the enzyme mix left out.)

**How much DNA, and why that much.** MAN0000470 page 21:

> - An equimolar amount of attB-PCR product (or linearized attB expression clone) and the
>   donor vector
> - 50 femtomoles (fmol) each […] is preferred, but the amount of attB-PCR product used may
>   range from 20–50 fmol
> - Note: 50 fmol of donor vector (pDONR221 or pDONR/Zeo) is approximately 150 ng
> - For large PCR products (>4 kb), use at least 50 fmol of attB-PCR product, but no more than
>   250 ng

with two caps: no more than 250 ng of donor vector, and no more than 0.5 μg of DNA in total,
"as excess DNA will inhibit the reaction". The conversion the manual prints is
`ng = (fmol)(N)(660 fg/fmol)(1 ng/10⁶ fg)`, worked as "(50 fmol)(2500 bp)(660 fg/fmol)(1
ng/10⁶ fg) = 82.5 ng of PCR product required" (page 22). A plan should compute that rather
than quote 150 ng, which is specific to a 4.3–4.8 kb donor.

**Longer, for a large insert.** "the length of the recombination reaction can be extended up
to 18 hours. An overnight incubation typically yields 5–10 times more colonies than a 1 hour
incubation. For large PCR products (≥5 kb), longer incubations […] are recommended."
(MAN0000470 page 23.) The product sheet puts the intermediate step in: 4–6 hours gives
"typically 2-3 fold more colonies" (11789.II.pps page 2).

**Topology matters.** "For BP reactions, the most efficient substrates are linear attB
products (PCR products or expression clones) and supercoiled attP-containing donor vectors."
(11789.II.pps page 2.)

**One number the manual contradicts itself on.** MAN0000470 page 21 opens "in a 20 μL BP
recombination reaction" and then caps everything per 10 μL, which is what its own protocol on
page 23 uses. MAN0000291 page 5 has the same paragraph with "10 μL". Read it as 10 μL for
Clonase II.

**No half-reaction exists.** The only scaling statement in any manual is "BP Clonase II enzyme
mix is supplied as a 5X solution. If you wish to scale the reaction volume, make sure the BP
Clonase II enzyme mix is at a final concentration of 1X" (11789.II.pps page 3).

**What the first paper used**, which is the only fully disclosed BP reaction anywhere: 300 ng
of the attP plasmid plus 2 μL of purified PCR product in 20 μL "that contained 4 μL BP Clonase
in 25 mM Tris HCl pH 7.5, 22 mM NaCl, 5 mM EDTA, 5 mM spermidine HCl, 1 mg/mL BSA. After
incubation for 60 min at 25°C, proteinase K (4 μg in 2 μL) was added, and each reaction was
incubated at 37°C for 20 min." (Hartley 2000, Methods.)

## 11. The LR reaction

The same shape, different amounts:

| Component | LR Clonase II, 10 μL | LR Clonase, 20 μL |
| --- | --- | --- |
| Entry clone | 1–7 μL, 50–150 ng | 1–10 μL, 100–300 ng |
| Destination vector, 150 ng/μL | 1 μL | 2 μL |
| 5X LR Clonase Reaction Buffer | in the mix | 4 μL |
| TE buffer, pH 8.0 | to 8 μL | to 16 μL |
| Enzyme mix | 2 μL | 4 μL |
| Incubation | 25 °C, 1 hour | 25 °C, 60 minutes |
| Stop | 1 μL proteinase K, 37 °C, 10 min | 2 μL proteinase K, 37 °C, 10 min |

(Clonase II from MAN0001032 revision A.0, page 3, and MAN0000470 page 32; Clonase I from
11791.pps, revision 02/26/03, page 3. The positive control is 2 μL of pENTR-gus at 50 ng/μL,
an entry clone carrying β-glucuronidase; the negative control omits the enzyme mix.)

**The amount has a ceiling with a named symptom.** MAN0001032 page 2:

> We recommend using 50–150 ng entry clone per 10 μL reaction. Highest colony yields are
> typically obtained using 150 ng entry clone and 150 ng destination vector. Do not use >150 ng
> entry clone as you may obtain colonies containing multiple DNA molecules (often with an
> associated "small colony" phenotype). Using <50 ng entry clone will generate fewer colonies.

**Topology, the other way round from BP.** "For LR recombination reactions, the most efficient
substrates are supercoiled attL-containing entry vectors and supercoiled attR-containing
destination vectors. For large (>10 kb) entry clones or destination vectors, linearizing […]
may increase the efficiency by up to 2-fold." (MAN0001032 page 2.) And the reaction can run to
18 hours, which the manual recommends for plasmids ≥ 10 kb (MAN0000470 page 32).

**LR Clonase II Plus is a different protocol, not a faster one.** MAN0001087 revision A.0
prescribes **16 hours** at 25 °C — not one hour extendable — and still needs proteinase K, 1 μL
at 37 °C for 10 minutes. Its table is in fmol: 10 fmol of each entry clone and 20 fmol of
destination vector per 10 μL. It ships no positive control, and it states no single-fragment
protocol at all; the support page allows the use but the guide does not table it.

**What the first paper used**: about 200 ng of entry clone with 300–400 ng of linearised
destination vector in 20 μL "containing 4 μL of LR Clonase, 50 mM of Tris HCl pH 7.5, 50 mM of
NaCl, 0.25 mM of EDTA, 2.5 mM of spermidine HCl, and 0.2 mg/mL of BSA", stopped the same way
(Hartley 2000, Methods). Note that this buffer is **not** the BP one — the two recipes differ,
and the paper adds that its Table 1 reactions "used earlier versions of DNAs, proteins, and
reaction conditions".

## 12. What is published about the enzyme mixes

**The proteins are stated; the buffer is not.** Every product sheet names the contents:

| Mix | Proteins, as the vendor lists them |
| --- | --- |
| BP Clonase and BP Clonase II | "the bacteriophage lambda recombination protein Integrase (Int), the E. coli-encoded protein Integration Host Factor (IHF)" |
| LR Clonase, LR Clonase II and LR Clonase II Plus | "the bacteriophage lambda recombination proteins Integrase (Int) and Excisionase (Xis), the E. coli-encoded protein Integration Host Factor (IHF)" |

(11789.pps, 11789.II.pps, 11791.pps, MAN0001032, MAN0001087, each in its "Description".) The
kit-contents tables then print one word in the composition column for both the enzyme mix and
the reaction buffer: **"Proprietary"** (MAN0000282 page viii; MAN0000470 page vii). By
contrast the same tables give the proteinase K stop solution in full — "2 μg/μL in: 10 mM
Tris-HCl, pH 7.5 / 20 mM CaCl2 / 50% glycerol" — and the purification solution as "30% PEG
8000/30 mM MgCl2". No manual states a unit definition, a protein concentration, or a salt
composition for either Clonase, and the buffers are not sold separately.

**So the composition is known only from outside the manuals**, at three different strengths:

1. **Hartley 2000's Methods** gives both reaction buffers at working strength, quoted in
   [10](#10-the-bp-reaction) and [11](#11-the-lr-reaction). This is peer-reviewed and it is
   the strongest source for the buffer.
2. **Landy (2015)**, *Microbiology Spectrum* 3(2):MDNA3-0051-2014, gives the requirement from
   the biology rather than the product: "Integrative recombination between supercoiled attP
   and linear attB requires the virally encoded integrase (Int) […] and the host-encoded
   accessory DNA bending protein integration host factor (IHF) […]. Excisive recombination
   between attL and attR to regenerate attP and attB additionally requires the phage-encoded
   Xis protein (which inhibits integrative recombination) […] and is stimulated by the
   host-encoded Fis protein." Two things there are worth carrying: the supercoiling
   requirement, which Hartley never mentions, and **Fis**, a fourth protein that no Clonase
   description lists.
3. **Invitrogen's own patent, US 8,241,896 B2** (Brasch, Cheo, Hartley, Temple; Life
   Technologies; granted 2012-08-14), Example 18, "Optimization of GATEWAY Clonase Enzyme
   Compositions", goes furthest: "Clonase: 50 ng/μl IntH6 and 20 ng/μl IHF, admixed in 25 mM
   Tris-HCl (pH 7.5), 22 mM NaCl, 5 mM EDTA, 1 mg/ml BSA, 5 mM Spermidine, and 50% glycerol",
   with a 40 μL reaction taking "8 μl Clonase (400 ng IntH6, 160 ng IHF)", and the LR version
   adding "40 ng of His 6 -carboxy-tagged Xis (Xis H6) in addition to the IntH6 and IHF". At
   1X the buffer matches Hartley's BP recipe component for component.

   **Three caveats travel with that third one** and belong beside any number taken from it: it
   was read from the Google Patents HTML rendering and has not been checked against the USPTO
   grant PDF; a patent is a legal disclosure, not a measurement anyone reproduced; and it
   describes the formulation at the 1998–99 priority date, not today's Clonase II, which the
   vendor calls a reformulation.

**Nothing here licenses a home-made Clonase.** No peer-reviewed reconstitution of either mix
was found, and the Limited Use Label License on every Clonase sheet is explicit about the
product: "The buyer cannot modify the recombination sequence(s) contained in this product for
any purpose."

## 13. Markers, and what each step selects

A Gateway vector carries two kinds of marker and they do different jobs:

> - ccdB gene located between the two attP sites for negative selection
> - Chloramphenicol resistance gene located between the two attP sites for counterselection
> - Kanamycin or Zeocin resistance gene for selection in E. coli
> - pUC origin for replication and maintenance of the plasmid in E. coli

(*Gateway pDONR Vectors* user guide, MAN0000291, part number 25-0531, revision date 29 March
2012, page 1, "Features".) A destination vector reads the same with attR in place of attP:
"Chloramphenicol resistance gene (CmR) located between the two attR sites for
counterselection; The ccdB gene located between the attR sites for negative selection;
Ampicillin resistance gene for selection in E. coli" (*E. coli Expression System with Gateway
Technology*, MAN0000278, revision date 8 August 2012, page 1).

So a vector's **backbone marker** — the one outside the att cassette — is what the plate
selects for, and the **cassette** carries ccdB plus CmR, which the reaction throws away. The
markers the manuals state for the common vectors:

| Vector | Size | Backbone marker | Cassette |
| --- | --- | --- | --- |
| pDONR221 | 4761 bp | kanamycin | ccdB 1197–1502, CmR 1825–2505, between attP1 570–801 and attP2 2753–2984 |
| pDONR/Zeo | 4291 bp | Zeocin | ccdB 1197–1502, CmR 1847–2506, between attP1 570–801 and attP2 2754–2985 |
| pDONR201 | not stated | kanamycin | not stated |
| pDEST14 / 15 / 17 / 24 | 6422 / 7013 / 6354 / 6961 bp | ampicillin | ccdB + CmR between attR1 and attR2 |
| pcDNA-DEST40 | 7143 bp | ampicillin (plus neomycin for mammalian cells) | ccdB + CmR, cassette 918–2601 |
| pcDNA3.2/V5-DEST, pDEST26, pDEST27 | 7711 / 7481 / 8123 bp | ampicillin (plus neomycin) | ccdB + CmR between attR1 and attR2 |
| Vector Conversion reading-frame cassette A / B / C.1 | 1711 / 1713 / 1714 bp | none — the user's vector supplies it | ccdB + CmR between attR1 and attR2 |

(pDONR coordinates from MAN0000291 page 14; pDEST sizes from MAN0000278 page 27; pcDNA-DEST40
from MAN0000223 revision 3.0; the mammalian series from MAN0000279 revision date 2 March 2012;
the cassettes from MAN0000469 page 16.)

**What each step plates on.** After BP, the donor vector's backbone marker: "pDONR221 | LB +
50 μg/mL Kanamycin; pDONR/Zeo | Low Salt LB + 50 μg/mL Zeocin selective antibiotic"
(MAN0000291 page 8). After LR, the destination vector's: "LB agar plates containing the
appropriate antibiotic to select for expression clones (e.g. ampicillin)" (MAN0000470 page
45). Zeocin carries a medium condition with it: "for Zeocin to be active, the salt
concentration of the bacterial medium must remain low (<90 mM) and the pH must be 7.5"
(MAN0000291 page 8).

**Why the two markers have to differ.** The vendor states it as a rule:

> Most entry vectors contain the kanamycin resistance gene for selection. For maximal
> compatibility within the Gateway Technology, we recommend that your vector not contain a
> kanamycin resistance marker. If this is unavoidable, you will need to perform the LR
> recombination reaction with an entry clone that carries a selection marker other than the
> kanamycin resistance gene.

(MAN0000469 page 1, and the same text in MAN0000470 page 33, which names the two escapes:
pDONR/Zeo for a Zeocin entry clone, or pCR8/GW/TOPO for a spectinomycin one.) The failure it
prevents is in the troubleshooting table: unreacted entry clone co-transforms and grows, and
the way to tell is that "small colonies often only grow on the selective plates used to select
for unreacted entry clones" (MAN0000470 page 41).

**What ccdB does about the rest.** MAN0000470 page 8:

> The presence of the ccdB gene allows negative selection of the donor and destination (and
> some entry) vectors in E. coli following recombination and transformation. […] When
> recombination occurs […] the ccdB gene is replaced by the gene of interest. Cells that take
> up unreacted vectors carrying the ccdB gene or by-product molecules retaining the ccdB gene
> will fail to grow. This allows high-efficiency recovery of the desired clones.

And a second, confirmatory screen rides on the cassette's other gene, because the product lost
it: "True expression clones will be ampicillin-resistant and chloramphenicol-sensitive.
Transformants containing a plasmid with a mutated ccdB gene will be both ampicillin- and
chloramphenicol-resistant. […] A true expression clone will not grow in the presence of
chloramphenicol." (MAN0000223 page 5, at 30 μg/mL; the same paragraph appears in MAN0000279
page 12 and MAN0000302 pages 14–15.)

**A plan can check all of this from the records themselves.** Nothing above needs a vector's
name. Find the att sites and the record classifies itself — attP means donor, attL entry, attR
destination, attB an expression clone or an attB substrate. Require a `ccdB` feature and a
chloramphenicol-resistance feature between the two att sites of a donor or destination record.
Read the plate antibiotic off the resistance feature that lies outside the cassette. Refuse a
pairing whose two backbone markers are the same. That is the whole selection check, and it is
a feature query on the user's own file.

## 14. How many colonies, and how many are right

**How many.** Two numbers, repeated in every document that states one:

> An efficient BP recombination reaction will produce >1500 colonies if the entire BP reaction
> is transformed and plated.
>
> An efficient LR recombination reaction will produce >5000 colonies if the entire LR reaction
> is transformed and plated.

(11789.II.pps page 3 and MAN0001032 page 3; MAN0000470 pages 26 and 32 repeat them with the
condition attached — "If you use E. coli cells with a transformation efficiency of
≥1 × 10⁸ cfu/μg".) That efficiency is the floor the reagents assume: "Any competent cells with
a transformation efficiency of >1.0 × 10⁸ transformants/µg may be used."

The transformation itself is 1 μL of the reaction into 50 μL of cells, 250 μL of S.O.C., and
20 μL plus 100 μL plated, with a 1:10 dilution before plating (MAN0000470 page 26). For LR
Clonase II Plus the counts are per fragment count: 2,000–15,000 colonies for 2 fragments,
1,000–5,000 for 3, and 50–500 for 4 (MAN0001087 page 2).

**How many are right.** No manual answers this, which is the largest gap in this note. What a
manual gives instead is a phenotype: the BP positive control makes entry clones that express
tetracycline resistance, and "True entry clones should be tetracycline-resistant"
(MAN0000470 page 24). The percentages that circulate are marketing — ">90% of the colonies
contain the Entry clone with the gene of interest in the correct orientation" in the 2022
brochure, "95% cloning efficiency" on the technology page — and a support answer gives 90% for
2-fragment MultiSite, "usually 70-90%" for 3-fragment, and "as high as 80% and as low as 30%"
for 4-fragment.

**The countable figures are in the paper**, and they are what a protocol should quote:

| What was counted | Result | Where |
| --- | --- | --- |
| BP: transformants from a tetR PCR product, scored on tetracycline | 195 of 197 | Hartley 2000, Results |
| LR: transformants from a tetR entry clone, same score | 96 of 102 | Hartley 2000, Results |
| LR into 12 different destination vectors: plasmid of expected size | 48 of 48 | Hartley 2000, Table 1 |
| the same, by restriction pattern | 24 of 24 | Hartley 2000, Table 1 |
| Negative controls without entry clone | "200- to 15,000-fold fewer colonies" | Hartley 2000, Results |

The paper also names the background's cause, which the vendor's troubleshooting agrees with:
"Background colonies contain inactive or deleted ccdB genes."

**And the one size effect anybody published.** "The efficiency of the in vitro RC reactions
decreases with increasing size of the DNAs involved, as judged by the number of colonies
produced. This effect can be minimized by using equal moles of DNA and by incubating for
longer times." (Hartley 2000, Discussion.) That is the evidence behind every "incubate
overnight for a large insert" line in [10](#10-the-bp-reaction) and
[11](#11-the-lr-reaction).

## 15. Validating a clone

**What the vendor publishes as primer sequences.** Three pairs, and only three:

| Primer | Sequence | Reads | Source |
| --- | --- | --- | --- |
| M13 Forward (−20) | `GTAAAACGACGGCCAG` | into an entry clone from pDONR221 or pDONR/Zeo | MAN0000291 page 11 |
| M13 Reverse | `CAGGAAACAGCTATGAC` | the other end of the same | MAN0000291 page 11 |
| Primer 1 | `CACATTATACGAGCCGGAAGCAT` | through the attR1 junction of a home-made destination vector | MAN0000469 page 12 |
| Primer 2 | `CAGTGTGCCGGTCTCCGTTATCG` | through the attR2 junction of the same | MAN0000469 page 12 |
| GW1 | `GTTGCAACAAATTGATGAGCAATGC` | out of attL1, in pCR8/GW/TOPO only | MAN0000437 page vii |
| GW2 | `GTTGCAACAAATTGATGAGCAATTA` | out of attL2, in pCR8/GW/TOPO only | MAN0000437 page vii |

The M13 sites sit at 537–552 and 3026–3042 in pDONR221, against attP1 at 570 and attP2 ending
at 2984 — so an M13 read reaches the insert only after crossing the att site. The manual puts
a number on that: the M13 sites are "located upstream and downstream of the attL1 and attL2
sites, respectively, requiring that at least 130 base pairs of vector-encoded DNA be read
before reaching the insert DNA" (MAN0000437 page 3).

**GW1 and GW2 do not generalise, and the manual says so.** They read from inside the attL
sites, but only because pCR8/GW/TOPO's attL2 was mutated for it: "Although other Gateway entry
vectors containing attL1 and attL2 sites are available, the GW1 and GW2 primers are only
suitable for use in sequencing inserts cloned into pCR8/GW/TOPO. This is because three
nucleotides within the attL2 site in pCR8/GW/TOPO have been mutated" (MAN0000437 page 13). A
plan that designs its own junction primer must therefore design it against the user's record,
not against an att site it assumes.

**There is no vendor attB1 or attB2 sequencing primer.** What is published under those names —
the adapter and forward or reverse primers in [7](#7-the-attb-primer-tail) — are PCR primers
for building the substrate.

**Restriction digest first, PCR second.** MAN0000291 page 10:

> Analyze the entry clones by restriction analysis to confirm the presence and correct
> orientation of the insert. Use a restriction enzyme or a combination of enzymes that cut once
> in the vector and once in the insert.
>
> You may also analyze positive transformants using PCR. Use a primer that hybridizes within
> the vector […] and one that hybridizes within your insert. You will have to determine the
> amplification conditions. If you are using this technique for the first time, you may want
> to perform restriction analysis in parallel. Artifacts may be obtained because of mispriming
> or contaminating template.

The colony PCR it describes is a lysis step and a normal amplification: 5 colonies each
resuspended in 50 μL of PCR mix, "Incubate the reaction for 10 minutes at 94°C to lyse cells
and inactivate nucleases", "Amplify for 20–30 cycles", 10 minutes at 72°C, then a gel
(MAN0000291 page 10). Per-cycle times are left to the user.

**Sequencing conditions, where the vendor gives them** (MAN0000470 page 28): at least 500 ng
of DNA, 5–50 pmol of primer, dye-terminator chemistry, and 95°C 5 min, then 30 cycles of 96°C
10–30 s / 50°C 5–15 s / 60°C 4 min.

**When sequencing stops being optional.** The one-tube BP-and-LR protocol skips the entry
clone, so "expression clones obtained using this protocol will be derived from entry clones
that are not unique. You will need to sequence your expression clone to confirm its identity."
(MAN0000470 page 45.) For the ordinary two-step route the manuals call it optional: "Optional:
To confirm that your gene of interest is in frame with the appropriate tag (if any), you may
sequence your expression construct." (MAN0000279 page 12.)

**A diagnostic band that means a specific failure.** "Entry clones migrate as 2.2 kb
supercoiled plasmids" is the signature of a BP reaction that cloned attB primer-dimers
(MAN0000470 page 44).

## 16. Troubleshooting

Only one document in the set carries troubleshooting tables for the reactions themselves: the
*Gateway Technology with Clonase II* user guide, MAN0000470, pages 40–44. The BP and LR
Clonase II product sheets have none; their equivalent advice is in "Important guidelines"
bullets, quoted where it adds a number. Columns are the vendor's: observation, reason,
solution.

**Few or no colonies, but the transformation control worked** (page 40):

| Reason | Solution |
| --- | --- |
| Incorrect antibiotic used to select for transformants | "Check the antibiotic resistance marker and use the correct antibiotic to select for entry clones or expression clones." |
| Reactions were not treated with proteinase K | "Treat reactions with proteinase K before transformation." |
| Used incorrect att sites for the reaction | "Use an entry clone (attL) and a destination vector (attR) for the LR reaction. Use an expression clone (or attB-PCR product) and a donor vector (attP) for the BP reaction." |
| Clonase II is inactive, or too little was used | Test another aliquot; store at −20°C or −80°C; "Do not freeze/thaw the Gateway Clonase II enzyme mix more than 10 times." |
| Used the wrong Clonase II | "Use the Gateway LR Clonase II enzyme mix for the LR reaction and the Gateway BP Clonase II enzyme mix for the BP reaction." |
| Too much attB-PCR product in a BP reaction | "Reduce the amount […] Remember to use an equimolar ratio of attB-PCR product and donor vector (i.e. ~50 fmol each)." |
| attB substrate ≥ 5 kb | "Incubate the BP reaction overnight." |
| Too much entry clone in an LR reaction | "Use equal fmol of destination vector and entry clone." |
| Destination vector or entry clone > 10 kb | "Incubate the LR reaction overnight. Linearize the destination vector or the entry clone. Relax the destination vector with topoisomerase I." |

**Two sizes of colony** (page 41): after LR, the small ones are usually unreacted entry clone
co-transforming, and restreaking tells them apart because "small colonies often only grow on
the selective plates used to select for unreacted entry clones". The fixes are 50 ng of entry
clone per 10 μL reaction, 1 μL into the transformation, and ampicillin raised to 300 μg/mL.
After BP, the same observation means the donor vector's ccdB gene has mutated or deleted — and
the tell is that "the negative control will give a similar number of colonies".

**High background in the absence of the entry clone** (page 42): the strain carries F′ and
therefore `ccdA`, so use one that does not; or the destination vector has lost part of ccdB,
so propagate it on the backbone antibiotic plus 15–30 μg/mL chloramphenicol and verify it
before use; or a solution is contaminated with another plasmid of the same resistance.

**High background of Zeocin-resistant transformants** (page 41): selection was not on Low Salt
LB.

**attB-PCR cloning** (pages 43–44):

| Observation | Reason | Solution |
| --- | --- | --- |
| Few or no colonies, though both controls worked | attB primers designed wrong | "Make sure that the attB PCR primers include four 5′ terminal Gs and the 25 bp attB1 or attB2 site." |
| the same | primers contaminated with incomplete sequences | HPLC- or PAGE-purified oligos, or the attB adapter PCR protocol |
| the same | product not purified enough | "Gel purify your attB-PCR product to remove attB primers and attB primer-dimers." |
| the same, product > 5 kb | too few attB molecules | "Increase the amount of attB-PCR product to 40–100 fmol per 20 μL reaction. Note: Do not exceed 500 ng DNA per 20 μL reaction. Incubate the BP reaction overnight." |
| the same | insufficient incubation | "Increase the incubation time of the BP reaction up to 18 hours." |
| Entry clones migrate as 2.2 kb supercoiled plasmids | the BP reaction cloned primer-dimers | PEG/MgCl2 or gel purification; a hot-start polymerase; redesign the primers to avoid mutual priming |

**Two numbers the 2003 guide gives differently.** MAN0000282 says "~100 fmol each" where
MAN0000470 says ~50 fmol, and "100 ng per 20 µl reaction" where MAN0000470 says 50 ng per
10 μL. Both are the same concentration, but quote the revision beside the number.

## 17. What it costs and how long it takes

**Time.** No manual in the set states an end-to-end figure, so this is the sum of the steps
each of which is quoted, and the overnights are left as overnights because the manuals never
put an hour on one.

| Step | Stated time | Source |
| --- | --- | --- |
| BP incubation | 1 hour at 25°C, extendable to 18 hours | MAN0000470 page 23 |
| LR incubation | 1 hour at 25°C, extendable to 18 hours | MAN0000470 page 32 |
| Proteinase K stop, either reaction | 10 minutes at 37°C | MAN0000470 pages 23, 32 |
| Transformation | 30 min on ice, 30 s at 42°C, 2 min on ice, 1 hour outgrowth at 37°C | MAN0000470 page 26 |
| Plates | overnight at 37°C | MAN0000470 page 26 |
| Picked colonies | overnight culture | MAN0000291 page 10 |
| PEG/MgCl2 purification of the PCR product | 15 min at 10,000 × g | MAN0000470 page 17 |
| Sanger cycle sequencing | 95°C 5 min, then 30 × (96°C 10–30 s, 50°C 5–15 s, 60°C 4 min) | MAN0000470 page 28 |

One reaction from set-up to plates is therefore about **2 h 45 min at the bench** on the
1-hour incubation, or about 19 h 45 min on the 18-hour one. A full BP then LR is two of those
blocks and four overnights — plates and culture after each — which is **five working days by
step count**, and that is a count of overnights rather than a time anyone published.

The vendor's own shortcut halves it. The one-tube protocol runs BP for 4 hours and LR for 2
in the same tube, about **7 h 45 min** to plating, and costs both yield and certainty: "fewer
expression clones will be obtained (at least 10–20% of the total number of expression clones)"
and "You will need to sequence your expression clone to confirm its identity" (MAN0000470
pages 45–46).

**Cost.** Thermo Fisher's US catalogue pages, read on 2026-09-18 without a login, so these are
public list prices and not what an institution pays:

| Catalogue number | Product | Size | List price (USD) |
| --- | --- | --- | --- |
| 11789020 | BP Clonase II Enzyme mix | 20 reactions | 407.00 |
| 11789100 | BP Clonase II Enzyme mix | 100 reactions | 1,724.00 |
| 11791020 | LR Clonase II Enzyme mix | 20 reactions | 450.00 |
| 11791100 | LR Clonase II Enzyme mix | 100 reactions | 1,792.00 |
| 12538120 | LR Clonase II Plus enzyme | 20 reactions | 1,062.00 |
| 12536017 | pDONR221 vector | 6 μg | 377.00 |
| A10460 | One Shot ccdB Survival 2 T1R competent cells | 11 × 50 μL | 349.00 |

Dividing gives about **USD 43 of enzyme per gene** for one BP plus one LR at the 20-reaction
list price, before vector, cells, media and sequencing. That division is this note's
arithmetic, not a vendor figure.

**The cost nobody bills for**, which only SnapGene's guide states plainly: the junction is a
scar, and the method has a set-up tax. "Large scars are left in your molecular constructs. For
simple applications of single fragment Gateway cloning, the scars have little impact." And:
"Gateway cloning has high barriers to entry. […] you may need to make several Gateway
compatible entry clones until you can fully take advantage of the system." That is the real
shape of the trade — Gateway is expensive for one construct and cheap for the fifth one into
the same entry clone.

## 18. Implications for this package

**For the att module.** What it owns is the arithmetic, not a catalogue:
positions 1–8 of a 25 bp region come from one partner and positions 9–25 from the other, both
directions ([5](#5-what-bp-writes), [6](#6-what-lr-writes)). Written that way, BP and LR are
one function and a direction, and the test is the round trip: BP then LR on the same insert
must give back the attB sites it started with. The eight sequences **may** ship
([1](#1-sources-and-licences)), which is what makes the built test fixtures of the spec
possible — but the plan should still recognise a vector by finding its sites rather than by
matching a constant, because the vendor says the sequences drift between products and
[2](#2-the-att-sites-base-for-base) shows two places where they do. A search that tolerates a
mismatch outside the 15 bp core, and none inside the 7 bp overlap, is the shape the sources
support.

**For the frame check.** Nine codons at each end, seven or eight of them fixed
([8](#8-reading-frame-across-a-junction)). At the N terminus, refuse `AA`, `AG` and `GA` as
the two added bases. At the C terminus, refuse a stop between the gene and attB2 for a fusion
and require one otherwise. Both are rules the manual states, so both carry a verdict.

**For the selection check.** Read the markers off the records
([13](#13-markers-and-what-each-step-selects)): a ccdB and a CmR feature inside the cassette,
one resistance feature outside it, and the two backbone markers of a pairing must differ. Add
the host rules, which are the part a bench cannot guess.

**For the protocol page.** Every number it prints is in this note or computed:
the reaction tables, the incubations and the proteinase K stop; the plate antibiotic named
from the vector's own marker; the expected colony counts as an order-of-magnitude check; the
chloramphenicol counter-screen at 30 μg/mL; the troubleshooting rows in
[16](#16-troubleshooting). Where a value is an Open gap, the page says so rather than
printing one.

## 19. Open gaps

Still unverified. None of these is guessed anywhere above.

| Item | Why it is missing | Workaround in use |
| --- | --- | --- |
| **What fraction of colonies is correct after BP or LR.** No manual states one | The figure exists only as marketing (">90%" in the 2022 brochure; "95%" and "Up to 95%" on the technology page) and in a support answer (90%, 70-90%, and 30-80% for MultiSite) | Hartley 2000's counted fractions: 195 of 197 tetracycline-resistant after BP, 96 of 102 after LR, 48/48 and 24/24 by size and digest. Cite those, not the percentages |
| **How ccdB Survival 2 T1R resists CcdB.** Its printed genotype has neither `gyrA462` nor a `ccdA` cassette | Every manual asserts only "is resistant to CcdB effects" | Name the strain the vendor names and do not explain the mechanism ([9](#9-ccdb-counter-selection)) |
| **What happens when an att site sits the wrong way round.** No source read here reports it | Hartley says only that aberrant events "have not been identified in hundreds of sequenced RC clones" | Refuse a substrate whose sites do not pair by number and point the same way, and say what was looked for ([4](#4-which-site-pairs-with-which)) |
| **No manual forbids plating a BP or LR reaction on a ccdB-resistant strain** | The only stated prohibition is against F' and `ccdA` strains | State the F' rule as a quote and the ccdB-resistant-host rule as an inference, marked as one |
| **The size of the four-G effect** | The supporting data is FIGS. 66 and 67 of US 7,670,823 B1, which are figure images | The direction is quoted in [7](#7-the-attb-primer-tail); no number is printed |
| **Clonase composition at product strength.** The manuals say "Proprietary" | US 8,241,896 B2's Materials block gives concentrations, but it was read from the Google Patents rendering rather than the USPTO grant PDF, and it describes the 1998-99 formulation, not Clonase II | Quote Hartley 2000's buffers, which are peer-reviewed, and carry the patent's figures only with those caveats attached ([12](#12-what-is-published-about-the-enzyme-mixes)) |
| **No half-reaction or scale-down protocol** | Not published | "Scale the volume, keep the enzyme mix at 1X" is the only statement ([10](#10-the-bp-reaction)) |
| **pDONR201's size, ccdB coordinates and map** | No manual for it exists on the route that works | Only its kanamycin marker and its first-generation attP sites are stated |
| **Whether DNA sequences specifically are uncopyrightable** | No source addresses it; the citable step is the general one about facts | [1](#1-sources-and-licences) states the argument and where it stops being cited |

## Sources

All retrieved 2026-09-18.

- Invitrogen / Thermo Fisher user guides, all fetched under
  `documents.thermofisher.com/TFS-Assets/LSG/manuals/`: *Gateway Technology*
  (MAN0000282, part 250522, revision 1.0, back cover 23 September 2003); *Gateway Technology
  with Clonase II* (MAN0000470, part 25-0749, revision date 2 April 2012); *Gateway pDONR
  Vectors* (MAN0000291, part 25-0531, 29 March 2012); *Gateway Vector Conversion System with
  One Shot ccdB Survival 2 T1R* (MAN0000469, part 25-0748, 23 March 2012); *E. coli Expression
  System with Gateway Technology* (MAN0000278, part 25-0517, 8 August 2012); *Mammalian
  Expression System with Gateway Technology* (MAN0000279, part 25-0518, 2 March 2012);
  *Gateway pcDNA-DEST40* (MAN0000223, revision 3.0); *pAd/CMV/V5-DEST and pAd/PL-DEST*
  (MAN0000302, part 25-0544, 4 January 2012); *pCR8/GW/TOPO TA Cloning Kit* (MAN0000437, part
  25-0706, 23 March 2012); *One Shot ccdB Survival 2 T1R* (MAN0000761 revision 2.0); *One Shot
  ccdB Survival T1R* (part C751003.pps, 11 June 2004)
- Product information sheets, same host: *BP Clonase II* (11789.II.pps, revision 31 October
  2010), *BP Clonase* (11789.pps, 10/31/10), *LR Clonase II* (MAN0001032 revision A.0),
  *LR Clonase* (11791.pps, 02/26/03), *LR Clonase II Plus* (MAN0001087 revision A.0, 26
  January 2022, and 12538.pps, 11 July 2006)
- Thermo Fisher, *Gateway Recombination and Seamless Cloning Support - Getting Started*;
  the Gateway technology page; the catalogue pages for 11789020, 11789100, 11791020, 11791100,
  12538120, 12536017 and A10460; *Website and Mobile Application Terms of Use*, effective
  17 June 2016
- Brasch, M., Cheo, D., Hartley, J. and Temple, G. *Compositions for use in recombinational
  cloning of nucleic acids*, US 7,670,823 B1, issued 2 March 2010, and US 8,241,896 B2,
  granted 14 August 2012; *Recombinational cloning using engineered recombination sites*,
  US 6,143,557 A. Read from `patents.google.com`
- USPTO, *Terms of use for USPTO websites*, and `uspto.gov/data.json`
- Hartley, J.L., Temple, G.F. and Brasch, M.A. (2000) DNA cloning using in vitro site-specific
  recombination. *Genome Res.* 10, 1788-1795.
  [doi:10.1101/gr.143000](https://doi.org/10.1101/gr.143000) (Cold Spring Harbor Laboratory
  Press)
- Landy, A. (2015) The lambda integrase site-specific recombination pathway.
  *Microbiol. Spectr.* 3(2):MDNA3-0051-2014.
  [doi:10.1128/microbiolspec.MDNA3-0051-2014](https://doi.org/10.1128/microbiolspec.MDNA3-0051-2014)
- Bernard, P. and Couturier, M. (1992) Cell killing by the F plasmid CcdB protein involves
  poisoning of DNA-topoisomerase II complexes. *J. Mol. Biol.* 226, 735-745 (abstract only);
  Bernard, P. et al. (1993) *J. Mol. Biol.* 234, 534-541 (abstract only); Dao-Thi, M.H. et al.
  (2005) *J. Mol. Biol.* 348, 1091-1102 (abstract only); Miki, T. et al. (1992) *J. Mol. Biol.*
  225, 39-52 (abstract only - a different `gyrA` allele, `Gly214Glu`, and not the source of
  `gyrA462`)
- Salim, L., Feger, C. and Busso, D. (2016) Data set for describing the elaboration of a
  compatible Gateway-based co-expression vector set. *Data in Brief* 9, 946-955.
  [doi:10.1016/j.dib.2016.11.013](https://doi.org/10.1016/j.dib.2016.11.013) (CC BY 4.0)
- NCBI: phage lambda genome, GenBank [J02459](https://www.ncbi.nlm.nih.gov/nuccore/J02459),
  fetched through E-utilities; deposited Gateway plasmids OQ700931, PP373655, KM880130,
  GU574775, OK571349 and OR569804; *NCBI Data Usage Policy*
- Addgene, *Plasmids 101: Gateway Cloning* (12 January 2017) and Terms of Use (24 January
  2023); SnapGene, *Gateway Cloning Technique* and legal disclaimers

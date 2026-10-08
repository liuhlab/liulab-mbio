---
search:
  exclude: true
---

# The working vector: pLVX-TetOne-Puro-GFP measured

Issue #230, under the AP-1 demo spec #225. `docs/research/ap1-demo-project.md` section 2.2 lists
four requirements a lentiviral parent has to meet and names none. The lab has now named
**pLVX-TetOne-Puro-GFP**, Addgene 171123. This note measures the real map against all four.

Every number below was counted on the map itself, with `liulab_mbio.sites` and
`liulab_mbio.translate`, on **2026-10-06**. Coordinates are the package's: 0-based, half-open.
Addgene's own `features.json` is 1-based inclusive and was converted once, at the boundary.

**Two requirements fail as written and one cannot be settled from this map.** The failures are
reported below with what meeting them costs. Sections 0 to 6 looked at no other parent, as #230
asks.

**Section 7 was added for #485**, under map #462, and answers a narrower question: of the eight
Type IIS sites section 1 counts, which can be domesticated and at what cost. It also reports the
Addgene search #485 asked for. Nothing in sections 0 to 6 was recomputed.

## 0. What was read, and the one thing that could not be

| What | Where | Read |
| --- | --- | --- |
| Catalogue record: backbone, size fields, marker, insert | <https://www.addgene.org/171123/> | 2026-10-06 |
| Addgene-verified full sequence 336792 | <https://www.addgene.org/browse/sequence/336792/> | 2026-10-06 |
| The bases | `…/336792/f8e81b51-…/addgene-plasmid-171123-sequence-336792-sequence-lines.json` | 2026-10-06 |
| The 32 annotations | the same path, `-features.json` | 2026-10-06 |
| The map image | the same path, `-map.png` | 2026-10-06 |
| Packaging limit | Kumar, M. et al. *Hum. Gene Ther.* **12**, 1893–1905 (2001), [doi:10.1089/104303401753153947](https://doi.org/10.1089/104303401753153947) — abstract only, PubMed 11589831 | 2026-10-06 |

The media path's stem is
`https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/336792/f8e81b51-dcce-4b08-86f1-1871c4158f37/`.
The files are under `reference_docs/synthesis_and_assembly/working-vector/`.

**No GenBank file was downloaded.** Addgene puts every sequence download behind a login, and
this run has no account. The bases were read out of the Sequence Analyzer's own render, where
each base is one SVG glyph, by `rebuild-sequence.py` beside the files. Two independent passes
over that render agree base for base.

The read was checked rather than trusted: at the coordinates Addgene's annotation gives,
`EGFP` reverse-complements to `ATGGTGAGC…TACAAGTAG` and `AmpR` to `ATGAGT…TGGTAA`, both clean
reading frames with one stop at the end and none inside. A base lost or gained anywhere before
position 8,580 would have broken the second of those.

**One measurement is open.** The analyzer renders **9,895 bp**; Addgene's catalogue field
*Total vector size* says **9,894 bp**, and its *Backbone size w/o insert* 9,173 plus *Insert
Size* 720 gives 9,893. The three disagree. Every length below is from the render, which is the
only one of the three that is a sequence rather than a typed-in field. What would close it: the
GenBank file behind the login.

## 1. Enzyme sites — **fails**

Requirement: no BsaI and no BsmBI outside the designed cassette.

**Eight sites: six BsaI, two BsmBI.** Seven of the eight lie between the LTRs, so they travel
into the cell. Counted with `liulab_mbio.sites.find_sites`; every one is certain.

| Enzyme | Site | Strand | Where it sits | Between the LTRs |
| --- | --- | --- | --- | --- |
| BsaI | 454–460 | + | 5' LTR, offset 454, first bases of TAR | yes |
| BsaI | 3724–3730 | + | hPGK promoter | yes |
| BsaI | 5730–5736 | − | PuroR, coding | yes |
| BsaI | 6472–6478 | − | between the WPRE and the 3' LTR, in no annotated element | yes |
| BsaI | 7112–7118 | + | 3' LTR, offset 454, first bases of TAR | yes |
| BsaI | 8720–8726 | − | AmpR, coding | no |
| BsmBI | 3876–3882 | + | hPGK promoter | yes |
| BsmBI | 5636–5642 | + | PuroR, coding | yes |

For reference, the same count gives **3 BbsI** (24, 6398, 6682), **2 SapI** (1106, 7637) and
**0 PaqCI**.

### What domestication costs, site by site

Four of the eight are cheap and four are not.

- **Coding, so a synonymous codon change does it**: BsaI 5730 and BsmBI 5636 in PuroR, BsaI 8720
  in AmpR. `liulab_mbio.sites.domesticate` is the tool, and PuroR and AmpR both translate
  cleanly, so it has a frame to work in.
- **No annotated element**: BsaI 6472, in the 206 bp between the end of the WPRE (6452) and the
  start of the 3' LTR (6658). Nothing on the map claims that stretch. Cheapest site on the list,
  and the only one where nothing argues against a change.
- **hPGK promoter**, BsaI 3724 and BsmBI 3876: regulatory sequence, not coding. A base change
  there is not synonymous in any sense the package can check, and the hPGK promoter is what
  drives the transactivator. Nothing read here measures how either position affects its
  strength.
- **The two LTR sites, BsaI 454 and 7112 — domestication may not touch these.** Both sit at the
  same offset, 454, in a 634 bp LTR, and their 36 bp neighbourhoods are identical:
  `GCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGA`. `GGTCTC` here is the first six bases of TAR, which
  begins the R region, which is the start of the packaged transcript. Three separate things make
  this expensive: the two LTRs are copies (the whole 634 bp differ at **one** base, offset 23),
  so a change has to be made twice and identically or reverse transcription's strand transfer is
  changing; R is transcribed, not merely replicated; and TAR is a structured hairpin where a
  base change is a change to a fold. **No source read here measures what a base change at TAR
  position 1–6 costs.** That is the hole, and it is the requirement's whole weight: until it is
  closed, the parent carries two BsaI sites that cannot be called removable.

### What the failure means for the method

The cargo arrives through BsaI and the final transfer uses BsmBI, so both enzymes see this
backbone. The honest reading of this table is that **requirement 2.2 as written is not met and
is not cheaply met**: three sites are a day's work, one is free, two are a research question and
two are a measurement nobody here has.

A narrower requirement would be met today. The sites that matter are the ones an assembly
reaction sees — BsmBI for the transfer into this vector, BsaI only if cargo is ever built in
place. **Both BsmBI sites are domesticable** (PuroR coding, hPGK promoter). Whether the
requirement should be narrowed to "no BsmBI outside the cassette, BsaI counted and recorded" is
a method decision, not this note's.

## 2. Bacterial marker — **passes**

Requirement: AmpR or CarbR, so the last transfer does not share a marker with the KanR DMX
backbone (departure D11).

**AmpR, and it is intact.** The map annotates `AmpR` at 8580–9441 on the reverse strand with its
own promoter at 9441–9546. Translated, that is 861 bp, 287 codons, `MSIQHFRVAL…ASLIKHW*`, one
stop and it is the last codon. Addgene's catalogue record says *Bacterial Resistance(s):
Ampicillin, 100 μg/mL*.

Ampicillin against the DMX backbone's kanamycin is a clean split. Nothing further is needed.

## 3. Polyadenylation signals — **fails as written**

Requirement: no polyadenylation signal anywhere between the 5' LTR and the 3' LTR, the stuffer
included; screen `AATAAA` and `ATTAAA` on both strands.

Run literally, the screen finds **16 hexamers on the plasmid** and **six strictly between the
LTRs** (634–6658), so the requirement fails. Run against what the requirement is protecting —
the packaged transcript, which reads the plus strand from the 5' LTR — **two of the six are in
scope and four are not**.

### The map's orientation is why

This is the measurement the rest of the section rests on. The parent holds its two cassettes
back to back, facing away from each other:

```text
5'LTR ─ Ψ ─ RRE ─ cPPT ─ [SV40 pA ← EGFP ← TRE3GS] ─ [hPGK → Tet-On 3G] ─ [SV40 → PuroR] ─ WPRE ─ 3'LTR
                              2186    2513   3242      3626       4155      4912    5250     5863   6658
```

`TRE3GS promoter` (3242–3610) and `EGFP` (2513–3233) are both annotated reverse, and the SV40
poly(A) signal at 2186–2321 sits at the far end of them. **The induced cassette reads on the
minus strand.** The transactivator and puromycin cassettes read forward, on the genome's own
strand, and end at the WPRE with no terminator of their own.

So the SV40 poly(A) signal, the one piece of this plasmid whose job is to terminate, is pointed
away from the packaged genome. That is a property of the parent the lab chose, not an accident
the method arranged, and it is what makes the parent usable at all.

### The six hexamers

| Position | Motif | Strand | Where | In the genome's path |
| --- | --- | --- | --- | --- |
| 880 | ATTAAA | + | gag fragment, no annotated element | **yes** |
| 1603 | AATAAA | + | between the RRE and gp41, no annotated element | **yes** |
| 2148 | AATAAA | − | just past the cPPT/CTS | no |
| 2192 | AATAAA | − | inside the annotated SV40 poly(A) signal | no |
| 5152 | AATAAA | − | SV40 promoter / SV40 ori | no |
| 5949 | ATTAAA | − | WPRE | no |

The four minus-strand hits do not terminate a plus-strand transcript. Position 2192 is the SV40
signal itself, doing the job the EGFP cassette needs it for.

The two plus-strand hits are both in the vector's native HIV-1 segment, the part every
third-generation transfer plasmid carries, and neither is annotated as a signal. They are
hexamers, not signals: a functional poly(A) site also needs a downstream GU/U-rich element, and
**nothing read in this run measures either position**. What would close it: a source that
measures 3'-end processing in a lentiviral transfer vector's gag/RRE region, or a direct
measurement of full-length genome on this backbone.

Two things are worth saying plainly about the requirement itself. First, **a hexamer screen can
never pass a lentiviral vector**: the LTR's own R region carries `AATAAA` by design — this map
has it at 526 and again at 7184, the same offset in both LTRs — and that signal is the one the
3' LTR must use. Second, the stuffer is in scope and clears: the span the method replaces,
`EGFP` at 2513–3233, holds no `AATAAA` or `ATTAAA` on either strand, so whatever replaces it
starts from zero and only the cargo has to be screened.

## 4. Size — **passes, with the limit stated as the source states it**

| | bp | From |
| --- | --- | --- |
| The parent, whole plasmid | **9,895** | the render; the catalogue field says 9,894 (section 0) |
| `EGFP`, the span the method replaces | **720** | the annotation, 2513–3233 |
| LTR outer edge to outer edge, the parent | **7,292** | 0–7292, the two 634 bp LTR annotations |
| The insert that replaces EGFP | **about 2,416** | section 2.3 and 3 of the AP-1 note: cargo about 2,330 at its largest, plus the 78 bp C-terminal part, plus `TATG` and `CTAA` |
| The plasmid afterwards | **about 11,591** | 9,895 − 720 + 2,416 |
| LTR to LTR afterwards | **about 8,988** | 7,292 − 720 + 2,416 |

Only the last row is the one that matters, and it is the largest library member. For the floor:
the parent with EGFP taken out and nothing put back is **6,572 bp** between the LTRs, so every
member lands between that and about 9.0 kb.

### The limit, as its source states it

Kumar et al. cloned bacterial DNA fragments of several lengths into three positions of a
lentiviral vector, pseudotyped with VSV-G, and titred. Their finding, quoted from the abstract:
"viral titers decreased semi-logarithmically with increasing vector length" and "**there appears
to be no absolute packaging limit** because measurable titers were obtained even when the
proviral length was in excess of 18 kb", with the caveat that "the very low titers of the larger
vectors will be of limited utility".

So there is no cliff to clear, and the requirement's "room for about 2.4 kb" is not a pass/fail
test but a cost. At about 9.0 kb between the LTRs the design is well inside the range Kumar
measured, and the price is some titre, paid against a parent that already sits at 7.3 kb empty.

**What is open: how much titre.** The abstract gives the shape of the curve and not its slope.
The slope is in the paper's figures, which are paywalled and were not read. What would close it:
the full text, or a titre measurement on this backbone at two cargo lengths. A library is
pooled, so a titre loss costs representation — section 9.9 of the AP-1 note is where that bound
is set, and it should be revisited once the slope is known.

## 5. Can the parent's own GFP be the stuffer? — **yes, and it should be**

Section 2.2 specifies an RFP stuffer. Section 7.1 chose mCherry as the transduction marker. Red
stuffer against a red marker means colour cannot separate a vector whose cassette was never
replaced from a cell that was simply transduced. Keeping the EGFP already in the parent fixes
that, and it is being cut out anyway.

Four things were measured, and all four hold:

- **EGFP carries no site of any enzyme in play.** Across its 720 bp: 0 BsaI, 0 BsmBI, 0 BbsI, 0
  PaqCI, 0 SapI, counting both strands. It can sit between the cassette's two BsmBI sites
  without being cut in the middle, which is the whole requirement on a stuffer.
- **It translates.** `MVSKGEELFTGV…GMDELYK*`, one stop, at the end.
- **It carries no poly(A) hexamer** on either strand (section 3), so retaining it does not
  reintroduce the risk requirement 3 is about.
- **It is already in the right place and the right orientation**, under the TRE3GS promoter, on
  the minus strand with the SV40 poly(A) signal beyond it.

The build is cheaper too. An RFP stuffer has to be sourced and inserted; the EGFP needs only its
two flanks replaced — `TATG` and a BsmBI site at one end, a BsmBI site and `CTAA` at the other —
which is one PCR off the parent with tailed primers.

**One thing this does not buy, and the RFP would not have bought it either.** The stuffer sits
under a TRE promoter, which is a minimal CMV promoter and is not read in *E. coli*. So neither a
green stuffer nor a red one colours a colony, and colony colour is not a screen for an
unreplaced vector at the cloning bench. Where the colour reads is in the transduced cell after
induction: green means the cassette was never replaced, red-only means it was. If a colony-level
screen is wanted, that is a bacterial promoter on the stuffer or a counter-selection marker, and
it is a separate decision. **Nothing read here measures TRE3GS activity in *E. coli*.**

**This is a departure from the method page** — section 2.2 says RFP — and wants an entry in
`docs/research/synthesis-and-assembly-departures.md` saying: the stuffer is the parent's own
EGFP, not an RFP, because the transduction marker section 7.1 chose is red and two reds cannot
be told apart. That file was not edited here; #269 recorded it there as **D17**, with the
colony-screen caveat below.

## 6. The verdict, and what it would take

| Requirement | Verdict | What it would take |
| --- | --- | --- |
| 1. No BsaI, no BsmBI outside the cassette | **fails** — 6 BsaI, 2 BsmBI | Eight base changes clear all eight sites, and **seven of the eight have a warrant**: four are free, three carry a published precedent and none carries a price. One, the BsmBI site in the hPGK promoter, has no source at all. Section 7 is the site-by-site table |
| 2. AmpR or CarbR | **passes** | nothing |
| 3. No poly(A) signal between the LTRs | **fails the literal screen; the genome's own strand carries two unannotated hexamers** | A source on 3'-end processing in the gag/RRE segment, or a direct full-length-genome measurement. The screen itself should be narrowed to the plus strand, since the LTR's own `AATAAA` makes a two-strand screen unpassable |
| 4. Room for about 2.4 kb | **passes** | The titre cost at about 9.0 kb between the LTRs. Kumar et al.'s full text, or a titre on this backbone |

**No alternative parent was considered here.** #230 says not to pick one, and nothing in this
section argues for it: requirement 1's expensive half is the LTR, which every lentiviral transfer
vector carries, and requirement 3's is the HIV-1 segment, which every one of them carries too. A
different parent changes the four cheap sites, not the two hard ones. **Section 7 went and
measured that**, against twenty-two deposits, and it holds.

## 7. Site by site: can each one be domesticated?

Issue #485, under map #462. Section 1 counted the eight sites; this section says of each one
whether a substitution is available and what it costs. Every count was reproduced on
**2026-10-08** against `tests/data/pLVX-TetOne-Puro-GFP.gb`, the record this repo ships, with
`liulab_mbio.sites`. The eight positions, strands and elements come back identical to section 1.

Four sequences were read the same day and the placements they give are this section's new
measurements: GenBank **K03455.1** (HIV-1 HXB2) and **AF324493.2** (pNL4-3), to say which bases
are inherited viral sequence; GenBank **NG_008862.1** (human *PGK1* RefSeqGene), to place the two
hPGK sites against that promoter's own transcription start; and the Addgene-verified sequence of
**Addgene 239691**, an hPGK promoter someone has already domesticated.

The verdict uses the ticket's three words. **Free** means a substitution costs nothing anyone can
measure. **Risky, with a source** means a published vector carries the same change at the same
position, and nobody priced it. **No source** means neither holds — ADR 0017's hole, not a claim
that the base cannot be changed.

### The table

| Enzyme | Position | Strand | Element | The substitution | Verdict |
| --- | --- | --- | --- | --- | --- |
| BsaI | 454–460 | + | LTR at 0–634, R +1..+6, the first six bases of TAR | R+4 `C`→`T` at 457 | **risky, with a source** |
| BsaI | 3724–3730 | + | hPGK promoter, *PGK1* −333 or −254 | base 2 `G`→`A` at 3725 | **risky, with a source** |
| BsmBI | 3876–3882 | + | hPGK promoter, *PGK1* −181 or −102 | one of seventeen | **no source** |
| BsmBI | 5636–5642 | + | PuroR, codon 129 | `GTC`→`GTG`, V | **free** |
| BsaI | 5730–5736 | − | PuroR, codon 160 | `GAG`→`GAA`, E | **free** |
| BsaI | 6472–6478 | − | the linker joining the WPRE to the HIV-1 remnant | one of seventeen | **free** |
| BsaI | 7112–7118 | + | LTR at 6658–7292, R +1..+6, the first six bases of TAR | R+4 `C`→`T` at 7115 | **risky, with a source** |
| BsaI | 8720–8726 | − | AmpR, codon 238 | `GGG`→`GGC`, G | **free** |

**Seven of the eight have a warrant. One does not.** That is the answer to the ticket's part (a),
and it is a better answer than section 6 carried: the two LTR sites now have a published
precedent, the WPRE-side site is better than unannotated, and of the two hPGK sites the BsaI one
has a precedent too. What is left standing is the single BsmBI site at 3876.

Addgene's own annotation names **both** LTR copies `3' LTR`, and this record reproduces it. The
copy at 0–634 is the 5' one.

### What the package does unaided

`sites.domesticate(record, ["BsaI", "BsmBI"])` makes **three** changes and refuses **five**. The
three are the coding sites in the table; the five come back as `outside_cds`, because taking a
site out of a non-coding element changes what the record spells and the package leaves that
choice to its caller. `docs/adr/0013-coding-by-bases.md` is why a feature's own bases, not its
GenBank key, decide whether it is coding — all 32 annotations on this map are `misc_feature`, and
PuroR and AmpR are still read in frame.

| Site | Coding sequence | Codon | Change | Amino acid |
| --- | --- | --- | --- | --- |
| BsmBI 5636 | PuroR | 129 | `GTC` → `GTG` | V |
| BsaI 5730 | PuroR | 160 | `GAG` → `GAA` | E |
| BsaI 8720 | AmpR | 238 | `GGG` → `GGC` | G |

Each is the codon *E. coli* K-12 uses most often among those that take the site away and spell no
new BsaI or BsmBI site. The protein is unchanged at every position. These three are **free**.

### One rule that holds at every site

Changing the first base of `GGTCTC` to `C` spells `CGTCTC`, which is BsmBI; changing the first
base of `CGTCTC` to `G` spells BsaI. Measured by putting all three alternative bases at all six
positions of each of the five non-coding sites and recounting: **17 of the 18 substitutions at a
site clear it, and the eighteenth swaps the enzyme and leaves the count at eight.** A sweep
against both enzymes therefore has seventeen choices per site, not eighteen. Nothing here is
forced onto that eighteenth base, so the rule costs this plasmid nothing — but it is one way a
hand-written domestication of a Type IIS site quietly fails.

### BsaI 6472 sits in a cloning linker, not in the virus — free

Section 1 called this one "in no annotated element". It is better than that. Walking the plasmid
against HXB2 in 20-mers puts the **first base of inherited HIV-1 sequence at 6480** (HXB2 8908),
and the WPRE ends at 6452. The 28 bases between them — `CTGGAATTAATTCTGCAGTCGAGACCTA` — carry a
PacI site and a PstI site and match neither the woodchuck WPRE nor HIV-1. They are the linker the
two were joined with, and **the BsaI site ends two bases before HIV-1 begins**.

The nearest element downstream is the 3' polypurine tract, the plus-strand primer, at 6641–6656,
**163 bases** away. Powell and Levin showed the 3' end of that tract is sequence-critical, so the
distance is worth measuring rather than assuming. Any of the seventeen substitutions does the job.

One observation fell out of the same walk, recorded because it is measured and not because it
bears on this ticket: this plasmid's 3' polypurine tract reads `AAAAGAAAAGAGGGG`, where **both**
HXB2 and pNL4-3 read `AAAAGAAAAGGGGGG` — one `G`→`A` at tract position 11, plasmid position 6651.
It is in Addgene's own render of the deposit, so it is the plasmid and not this repo's copy of it.
Powell and Levin put a change at that tract's 5' end among the tolerated ones, and the vector is
in use. Nobody has measured what it costs.

### The two LTR sites: one base, in both copies, with one published precedent

Both sit at LTR offset 454 in a 634 bp LTR, which `docs/research/lentiviral-tolerance.md` section
1 fixes as transcription position **+1**, from the TATA box 28 bases upstream. So `GGTCTC` is R +1
to +6: the first six bases the vector transcribes, and the first six bases of the TAR hairpin's
lower stem. The two LTR copies differ at exactly one base in 634, and not at this one, so a change
has to be made twice and identically.

What the literature allows at these positions is settled in that note's section 3 and is not
re-derived here. Its result in one line: **the TAR lower stem is load-bearing by base-pairing and
not by base identity**, so a substitution there changes a fold rather than a motif.

| Source | What it measured | Bearing on R +1..+6 |
| --- | --- | --- |
| Klaver & Berkhout, *EMBO J.* **13**, 2650–2659 (1994) | Substitutions in the TAR lower stem of replication-competent HIV-1 | Dead; revertants restored pairing, not sequence. The defect was transcription from the integrated provirus, which a Tat-free third-generation vector does not have |
| Das, Klaver & Berkhout, *J. Virol.* **72**, 9217–9223 (1998) | Mutant Xho+10 replaces nt +3 to +16 — four of these six bases and ten more | Vector RNA about 50% of wild type, packaging 70%; a revertant that restored lower-stem pairing restored packaging |
| Das, Vrolijk, Harwig & Berkhout, *Retrovirology* **9**, 59 (2012) | Deletions on one side of the TAR stem against both sides | Opening the duplex hurt packaging; shortening it while keeping it paired did not |
| Haellman *et al.*, *Metab. Eng.* **66**, 41–50 (2021) — VAMSyB | A third-generation transfer vector carrying `GGTCTC`→`GGTTTC`, **C→T at R+4, in both LTRs** | It packages, transduces, integrates and drives a dose-responsive circuit. **No titre appears anywhere in the paper** |

So a substitution is available, it is one base, and one group has already made it and shipped the
vector. What nobody has is the price: no titre has been published for any lentiviral vector with
a changed R. That is why these two are **risky, with a source** rather than **free**.

Three conditions ride with them, all from `lentiviral-tolerance.md` section 3: change **both**
copies identically; prefer nt 4, 5 or 6 and leave nt 1 alone, because nt 1 is the transcription
start site; and check the resulting fold against the neighbouring poly(A) hairpin rather than
assuming a register. The published edit meets the first two. It meets the third by not being
tested.

### The two hPGK sites: native promoter, and one of them already domesticated

**The fragment is the real promoter, not a construct.** The 511 bases Addgene annotates as `hPGK
promoter`, 3626–4137, are a contiguous byte-identical copy of **NG_008862.1 positions 4649–5159**,
with no break anywhere across the walk. Both sites are therefore native human *PGK1* sequence, and
neither is a cloning scar inside the annotation.

**Where the sites sit depends on which RefSeq you number from, and the gap is 79 bases.** The 5'
end of **NM_000291.3** lands at fragment position 353 and the 5' end of **NM_000291.4** at
position 432; both anchor exactly. Singer-Sam *et al.* reported three transcription start points
for this gene in 1984, so the ambiguity is in the biology as well as in the annotation.

| Site | On NM_000291.3 | On NM_000291.4 |
| --- | --- | --- |
| BsaI 3724 | −254 | −333 |
| BsmBI 3876 | −102 | −181 |

That 79 bases is not a rounding difference, because two published landmarks fall inside it.
McBurney *et al.* mapped the **mouse** *Pgk-1* core promoter to the 120 bases upstream of the
start site, with a 320 bp orientation-independent activator above it, and Adra *et al.* report the
mouse and human promoters are homologous. Yang *et al.* place two DNase I footprints at centres
roughly **−360** and **−130**. So on one numbering the BsmBI site is clear of the core promoter
and on the other it is inside it, and on one numbering the BsaI site is a footprint's half-width
away and on the other it is not. **Nothing read here settles which numbering the footprints used.**

**Nothing measures a single-base substitution anywhere in this promoter.** No scanning or
saturation substitution series for human *PGK1* with a reporter readout was found. What exists is
footprinting, which locates protein contacts and measures no expression, and mouse deletion
mapping at 120 and 320 base resolution. Pfeifer *et al.* found eight in vivo protein–DNA contact
regions across 450 bases of this promoter, four of them Sp1 consensus; **their coordinates are in
figures behind a paywall and were not retrieved**, and they are what would settle both sites at
once.

For a base rate rather than a position, the nearest measurement is saturation mutagenesis of other
human promoters. Kircher *et al.* tested 31,243 single substitutions across twenty disease-linked
regulatory elements, ten of them promoters. **About 15% changed expression significantly**; of
those, the median change was 20–24%, and **only about 2% of all substitutions reached two-fold**.
Significant positions clustered in transcription-factor binding sites rather than spreading
evenly. *PGK1* was not among the elements tested, so carrying that rate here is an inference and
is marked as one.

#### The BsaI site has a precedent, read first-hand

**Addgene 239691, `pKG1367_pPV1-hPGK`**, deposited by the Galloway lab for Peterman *et al.*,
*Nucleic Acids Research* **53**, gkaf528 (2025), carries the Mutation field "Removal of BsaI cut
sites" on an insert named `hPGK promoter`. Its Addgene-verified sequence was read on **2026-10-08**
and compared base by base with this plasmid's fragment. **Measured:**

- Across the 66 bases around the BsaI site, the two differ at **one** position: `GGTCTC` →
  `GATCTC`, a `G`→`A` at the **second base of the site**, this plasmid's position 3725.
- Across the 66 bases around the BsmBI site, the two are **identical**. `CGTCTC` is still there.
  They removed BsaI and left BsmBI, which is what the Mutation field says.

So an hPGK promoter carrying exactly the change this site needs is deposited, published and in
use, and the paper that deposited it profiles transcription and translation from hPGK among other
promoters. **What it does not supply is a number:** nothing read here reports that promoter's
output against the undomesticated hPGK, so the precedent says the change is survivable and says
nothing about its cost. That is the same shape as the LTR precedent, and it earns the same
verdict: **risky, with a source**.

#### The BsmBI site has none

Nobody has published an hPGK promoter with `CGTCTC` at this position removed — the one group who
domesticated this promoter left it standing. Its verdict is **no source**. That is not a claim the
base cannot be changed: it is one base, seventeen substitutions clear it, and the general promoter
base rate above says most single substitutions in a human promoter cost little. It is a claim that
nothing licenses a number, and ADR 0017 is why that is an acceptable place for a plan to stop.

Two things sharpen the choice of base, and both are inferences rather than measurements. The
promoter is a CpG island — 62 CpGs in 511 bases, 69.3% GC — and its active state goes with
complete demethylation, so **prefer a substitution that neither destroys nor creates a CpG**. At
the BsaI site, bases 1, 2, 3 and 5 qualify, because base 6 pairs with the `G` that follows it; at
the BsmBI site, bases 3 to 6 qualify, because bases 1 and 2 are themselves a CpG. The precedent's
`G`→`A` at base 2 satisfies this. Second, this promoter drives the Tet-On 3G transactivator, so a
loss here moves the whole doxycycline dose–response of the induced cassette and not just its
ceiling.

### The sweep, measured end to end

**Eight bases clear all eight sites.** The three codon changes the package makes, plus `C`→`T` at
457 and 7115, `G`→`A` at 3725, `C`→`G` at 3881 and `G`→`A` at 6474, give a plasmid with **0 BsaI
and 0 BsmBI**. The two LTR copies still differ at the one base they differed at before, and the
hPGK promoter still holds 62 CpGs. Counted on the edited record on 2026-10-08.

Two things the sweep does not clear, and neither is requirement 1's business as section 1 writes
it: **BbsI stays at three**, two of them the LTR copies at offset 25 in U3, and SapI stays at two.
PaqCI, SrfI and PmeI are zero before and after.

### Part (b): is there a backbone with fewer conflicts?

**No.** Twenty-two candidate deposits were counted the same way on **2026-10-08**, whole plasmid,
both strands, circular, straight from each one's Addgene-verified sequence. The pipeline reproduces
this plasmid's 6 BsaI and 2 BsmBI exactly, which is what makes the other counts worth reading.

| Plasmid | Addgene | bp | BsaI | BsmBI | Equivalent function? |
| --- | --- | --- | --- | --- | --- |
| pLVX-TetOne-Puro-GFP, the parent | 171123 | 9,895 | 6 | 2 | — |
| pCW57-MCS1-2A-MCS2 | 71782 | 7,709 | 5 | 3 | yes |
| pCW57.1 | 41393 | 9,354 | 6 | 3 | yes |
| pLIX_402 | 41394 | 9,394 | 6 | 3 | yes |
| pLIX_403 | 41395 | 9,396 | 7 | 3 | yes |
| pLV-TRE-DEST-UbC-tTSrtTA-IRES-Puro | 253781 | 12,399 | 6 | 3 | yes |
| pInducer21 | 46948 | 12,145 | 4 | 2 | no — no antibiotic marker |
| pInducer20 | 44012 | 11,809 | 5 | 2 | no — neomycin, not puromycin |
| Tet-pLKO-puro | 21915 | 10,633 | 5 | 2 | no — shRNA only, no site for an ORF |
| LT3GEPIR | 111177 | 10,272 | 5 | 2 | no — miR-E cassette |
| pJT050, pJT126 | 161925, 161926 | 8,921, 8,990 | 4 | 2 | no — blasticidin, and recruitment rather than TRE expression |
| pLV-TRE3G-hYAP1-CMV-mCherry-T2A-Puro | 176850 | 10,089 | 5 | 1 | no — no transactivator on the plasmid |
| pLENTI_TetOn_Hu_BIRC5 | 136348 | 8,916 | 3 | 2 | no — hygromycin, and needs a separate rtTA vector |

**Nothing that is genuinely equivalent comes in lower than eight.** The closest functional
equivalent ties at eight with a different split, and every candidate scoring five to seven fails
on function: no marker, the wrong marker, no room for an open reading frame, or no transactivator
on the plasmid.

**The two-site floor is now measured rather than argued.** Section 6 reasoned that a different
parent changes the cheap sites and not the hard ones, because every lentiviral transfer vector
carries the LTR. Searching all twenty-two full sequences for the TAR 20-mer
`GGTCTCTCTGGTTAGACCAG` finds it **exactly twice in every one of them**, one copy per LTR, including
in the three-BsaI vector above. **Not one deposit has a recoded LTR.** The best any of them reaches
is two LTR sites plus one more, and that vector is not all-in-one.

**No Type IIS-domesticated lentiviral Tet-On vector exists.** The Mammalian ToolKit ships
lentiviral destination vectors whose BsaI and BsmBI sites are their cloning site by design, and
the Bintu lab's BsmBI backbones are the same case. Nobody ships an assembled Tet-On lentiviral
vector that is site-free.

So section 6's expectation holds, and the search that measured it also closed the open question
`lentiviral-tolerance.md` section 9 left standing — "whether any Tet-on lentiviral backbone is
already BsaI-clean". None is.

### What is still open after this section

| Open | What would close it |
| --- | --- |
| **The titre cost of the R+4 change.** VAMSyB is a precedent, not a number | A side-by-side titre against the undomesticated parent. One experiment, and it prices both LTR sites at once |
| **The output cost of either hPGK change.** The Galloway deposit is a precedent for the BsaI one and a number for neither | A reporter assay, or Peterman *et al.*'s own measurement read against an undomesticated hPGK if they ran one |
| **Which numbering the hPGK footprints use**, and so whether either site lands in one | Pfeifer *et al.* 1990's figure coordinates. The paper is paywalled at `genesdev.cshlp.org`, is not in PMC, and was not retrieved |
| **Whether the one-base difference already standing between this plasmid's two LTR copies** behaves like a deliberate asymmetry | It is suggestive, not a measurement |

## Sources

- <https://www.addgene.org/171123/> and <https://www.addgene.org/browse/sequence/336792/>, both
  read 2026-10-06. The sequence is listed there as an Addgene-verified full sequence. Backbone
  pLVX-TetOne-Puro, manufacturer Clontech.
- Kumar, M., Keller, B., Makalou, N. & Sutton, R. E. Systematic determination of the packaging
  limit of lentiviral vectors. *Hum. Gene Ther.* **12**, 1893–1905 (2001).
  [doi:10.1089/104303401753153947](https://doi.org/10.1089/104303401753153947). **Abstract only**
  — the full text is paywalled and was not read.
- `docs/research/ap1-demo-project.md` sections 2.2, 2.3, 3 and 7, for the four requirements, the
  cassette layout and the marker choice. No number in it is recomputed here.
- GenBank **K03455.1** (HIV-1 HXB2, 9,719 bp), **AF324493.2** (pNL4-3, 14,825 bp),
  **NG_008862.1** (human *PGK1* RefSeqGene), **NM_000291.3** and **NM_000291.4**, all read
  2026-10-08 from NCBI, for section 7's placements.
- Addgene 239691, `pKG1367_pPV1-hPGK`, Galloway lab, Mutation field "Removal of BsaI cut sites".
  <https://www.addgene.org/239691/>, Addgene-verified sequence 477086, read 2026-10-08. Its
  publication: Peterman, E. L. *et al.* High-resolution profiling reveals coupled transcriptional
  and translational regulation of transgenes. *Nucleic Acids Res.* **53**, gkaf528 (2025).
  [doi:10.1093/nar/gkaf528](https://doi.org/10.1093/nar/gkaf528), PMID 40530694. **Abstract only**
  — the comparison in section 7 is of the deposited sequence, not of the paper's figures.
- `docs/research/lentiviral-tolerance.md` sections 1, 2 and 3, for the LTR, TAR and VAMSyB
  measurements, and its section 8 for those sources in full. No number in it is recomputed here.
- Singer-Sam, J. *et al.* Sequence of the promoter region of the gene for human X-linked
  3-phosphoglycerate kinase. *Gene* **32**, 409–417 (1984).
  [doi:10.1016/0378-1119(84)90016-7](https://doi.org/10.1016/0378-1119(84)90016-7), PMID 6099325.
  **Abstract only.**
- Pfeifer, G. P., Tanguay, R. L., Steigerwald, S. D. & Riggs, A. D. In vivo footprint and
  methylation analysis by PCR-aided genomic sequencing. *Genes Dev.* **4**, 1277–1287 (1990).
  [doi:10.1101/gad.4.8.1277](https://doi.org/10.1101/gad.4.8.1277), PMID 2227409. **Abstract only**
  — the full text is paywalled, is not in PMC, and the footprint coordinates were not retrieved.
- Yang, T. P., Singer-Sam, J., Flores, J. C. & Riggs, A. D. DNA binding factors for the CpG-rich
  island containing the promoter of the human X-linked PGK gene. *Somat. Cell Mol. Genet.* **14**,
  461–472 (1988). [doi:10.1007/BF01534712](https://doi.org/10.1007/BF01534712), PMID 3175764.
  **Abstract only** — it gives the two footprint centres, not their extents.
- McBurney, M. W. *et al.* The mouse Pgk-1 gene promoter contains an upstream activator sequence.
  *Nucleic Acids Res.* **19**, 5755–5761 (1991).
  [doi:10.1093/nar/19.20.5755](https://doi.org/10.1093/nar/19.20.5755), PMID 1945853.
  **Abstract only.** And Adra, C. N., Boer, P. H. & McBurney, M. W. *Gene* **60**, 65–74 (1987).
  [doi:10.1016/0378-1119(87)90214-9](https://doi.org/10.1016/0378-1119(87)90214-9), PMID 3440520,
  **abstract only**, for the mouse–human homology that carries it across.
- Kircher, M. *et al.* Saturation mutagenesis of twenty disease-associated regulatory elements at
  single base-pair resolution. *Nat. Commun.* **10**, 3583 (2019).
  [doi:10.1038/s41467-019-11526-w](https://doi.org/10.1038/s41467-019-11526-w), PMID 31395865,
  PMC6687891. **Full text read**, for section 7's base rate. *PGK1* is not among its elements.
- Powell, M. D. & Levin, J. G., for the 3' polypurine tract, by way of
  `docs/research/lentiviral-tolerance.md` section 4, which cites it in place.
- The twenty-two Addgene deposits of section 7's part (b), each read 2026-10-08 through the public
  sequence viewer. Their catalogue numbers are in the table.
- `reference_docs/synthesis_and_assembly/working-vector/README.md`, for which file is which.

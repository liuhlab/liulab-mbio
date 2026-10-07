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
reported below with what meeting them costs. No other parent was looked at, as #230 asks.

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
| 1. No BsaI, no BsmBI outside the cassette | **fails** — 6 BsaI, 2 BsmBI | Three synonymous changes and one free one clear four. Two in the hPGK promoter are untested changes to regulatory sequence. **Two in TAR, one per LTR, cannot be called removable without a source** |
| 2. AmpR or CarbR | **passes** | nothing |
| 3. No poly(A) signal between the LTRs | **fails the literal screen; the genome's own strand carries two unannotated hexamers** | A source on 3'-end processing in the gag/RRE segment, or a direct full-length-genome measurement. The screen itself should be narrowed to the plus strand, since the LTR's own `AATAAA` makes a two-strand screen unpassable |
| 4. Room for about 2.4 kb | **passes** | The titre cost at about 9.0 kb between the LTRs. Kumar et al.'s full text, or a titre on this backbone |

**No alternative parent was considered.** #230 says not to pick one, and nothing here argues for
it: requirement 1's expensive half is the LTR, which every lentiviral transfer vector carries,
and requirement 3's is the HIV-1 segment, which every one of them carries too. A different
parent changes the four cheap sites, not the two hard ones.

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
- `reference_docs/synthesis_and_assembly/working-vector/README.md`, for which file is which.

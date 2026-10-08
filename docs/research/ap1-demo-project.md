---
search:
  exclude: true
---

# The AP-1 demo project

Specification for issue #217, under the spec at #210. One real project against
`docs/synthesis-and-assembly.md`, end to end, so development has something concrete to build
toward.

**Everything decided here is decided for this project only.** The method page keeps carrying
each question as open, and nothing below is a package default. Where this note answers a
decision the method defers, the heading says so.

The evidence behind every option is in `docs/research/synthesis-and-assembly.md`; the source the
cargo comes from is in `docs/research/protein-library-assembly.md`; the barcode rules are in
`docs/research/barcode-design.md`. Those notes are cited rather than restated, and no
measurement in them is contradicted here. Numbers this note derives are arithmetic over the
layout in section 3 and the two inputs in section 4.1, and are marked **derived** where they
appear. A price is a list price with its date, which is a reference and never a fact about what
this project pays.

## 1. The cargo input, and what is held back

The project builds a combinatorial AP-1 transcription-factor library over three positions:
**N**, **DBD** and **C**. The source calls the middle position bZIP; this project says DBD, and
they are the same position.

**The input is amino acid sequences and a position per sequence. Nothing else.** For each AP-1
family member: three protein sequences, each labelled N, DBD or C.

| Measured against the source | Value |
| --- | --- |
| Parts per position, counted on each sheet | **24** |
| Parts in all | 24 + 24 + 24 = **72** |
| Library members, three positions | 24³ = **13,824** |
| Published domain sequences, summed | 23,103 bp |
| Domain coding lengths | N 44-1,046 bp, DBD 191-263 bp, C 27-942 bp |

Twenty-four parts per position is what is measured; that they are 24 AP-1 family members each
split three ways is the source's own framing and the issue's, and nothing in this project turns
on it. What the design needs is 24 sequences per position and 72 in all.

Those counts are measured in `docs/research/protein-library-assembly.md` (section 6 and section
7), not taken from the issue. The same note records that the source's main text gives the
domains as "50-1,000 bp" while its own sequences measure 44-1,046 bp; the measured range is the
one this project designs against.

**Held back, and never read by the design:** the source's codon choices, its overhangs, its
internal and external stuffers, its barcode set, its padding, its primers, and the residue
changes it made at domain boundaries. Those are what this package is for, and feeding them in
would defeat the test. They are opened only in section 8, after the design is frozen and
written out.

The amino acid sequences themselves are read from the source's Table S1, which
`docs/research/protein-library-assembly.md` section 1 clears for reuse under CC BY 4.0 with
attribution. The attribution travels with the design output.

## 2. The vectors

### 2.1 The DMX vector

The lab-resources build: a rebuilt DMX parent, BbsI-free, one PmeI site outboard of each BsaI
site. Which parent is a method decision and not this project's, and this project uses whatever
the lab-resources build settles on.

`scripts/build_dmx_vector.py` now builds it from DMX0001 and the demo is planned against the
result. The worry recorded here — that inserting PmeI outboard of the right-hand cassette BsaI
lands in or abuts the `6xHis` coding sequence — is measured and wrong: the window runs 416..503
and `6xHis` ends at 432, so the site goes at 437 and sacrifices nothing. The marker is KanR:
departure D11 forces the swap, because the final transfer moves cargo into an AmpR or CarbR
working vector. `docs/research/dmx-destination.md` carries both.

### 2.2 The working vector: a Tet-on lentiviral backbone

A third-generation lentiviral transfer plasmid carrying the assembly cassette
`[TATG.BsmBI]─[stuffer]─[BsmBI.CTAA]` under a TRE promoter, with the transactivator and the
transduction marker expressed constitutively from a second promoter between the LTRs.

**The parent is named: pLVX-TetOne-Puro-GFP, Addgene 171123.** #230 measured the real map
against all four requirements, on 2026-10-06, in the package's own coordinates.
`docs/research/working-vector-plvx-tetone.md` holds every count and where it came from; the
requirements are restated below in the form that measurement left them in, and no number is
recomputed here.

1. **No BsmBI and no cargo-enzyme site outside the designed cassette.** Measured: six BsaI, two
   BsmBI. #256 settled the requirement by narrowing it — a site only has to go where its enzyme
   shares a reaction with this vector — and chose **PaqCI as the cargo enzyme**, which this map
   carries **zero** of across all 9,895 bp. **Two sites remain to remove, both BsmBI:** 5636 in
   PuroR, synonymous and free; 3876 in the hPGK promoter, edited in place, the only untested
   change on the vector and readable as the stuffer's doxycycline response. All six BsaI sites
   stay, the two in the LTRs included, because BsaI never shares a reaction with this vector.
   **Passes, at two edits.**
2. **AmpR or CarbR**, because the DMX backbone the cargo leaves is KanR and the last transfer
   must not share a marker (departure D11). Measured: AmpR, intact, translating with one stop
   and it the last codon. **Passes.**
3. **No polyadenylation signal between the LTRs on the plus strand**, the stuffer included. The
   genome is packaged as one transcript from the 5' LTR, so an internal terminator on that
   strand truncates it. **The plus strand is a narrowing, and the reason is below.**
4. **Room for a cargo of up to about 2.4 kb** (**derived**, section 3) plus the working
   cassette. **This is a cost, not a gate** — also below.

#### Why requirement 3 reads one strand and not two

**A two-strand hexamer screen cannot pass any lentiviral vector.** The LTR's own R region
carries `AATAAA` by design — on this map at 526 and again at 7184, the same offset in both
copies — and the 3' LTR's copy is the signal the vector has to use. A screen that fails on the
thing making the vector work is the screen being wrong, not the vector. Narrowed to the plus
strand it tests what it was written for: termination of the packaged transcript.

Run literally, `AATAAA`/`ATTAAA` on both strands gives 16 hexamers on the plasmid and six
strictly between the LTRs. **Two of the six are in scope**, both plus-strand, both unannotated,
both inside the native HIV-1 gag/RRE segment every third-generation transfer plasmid carries:

| Position | Motif | Where it sits |
| --- | --- | --- |
| 880 | `ATTAAA` | gag fragment, no annotated element |
| 1603 | `AATAAA` | between the RRE and gp41, no annotated element |

**Whether either is functional is open.** A hexamer is not a signal — a working poly(A) site
also needs a downstream GU/U-rich element — and nothing read in either run measures either
position. The other four hits are minus-strand and do not terminate a plus-strand transcript;
one of them is the SV40 signal doing the job the stuffer cassette needs it for (section 2.3).
Neither hexamer is replaceable in any case, because both sit in the segment the vector needs, so
what the open item decides is whether this parent carries a known risk or an inert six-mer — not
what gets built.

#### Why requirement 4 is a cost and not a gate

| | bp | From |
| --- | --- | --- |
| The parent | **9,895** | the sequence render; see the length disagreement below |
| The stuffer span replaced, `EGFP` | **720** | the annotation at 2513-3233 |
| LTR outer edge to outer edge, before | **7,292** | the two 634 bp LTR annotations |
| The insert that replaces the stuffer | about **2,416** | sections 2.3 and 3, at the largest member |
| LTR to LTR, after | about **8,988** | 7,292 − 720 + 2,416 |

Kumar et al. 2001, *Hum. Gene Ther.* **12**, 1893–1905 — **abstract only; the full text is
paywalled and was not read** — reports titres falling "semi-logarithmically with increasing
vector length" and finds **no absolute packaging limit**, measurable past 18 kb. There is no
cliff to clear at about 9.0 kb between the LTRs, so the requirement buys a price rather than a
pass: some titre, paid against a parent already 7.3 kb empty.

**Why the price lands harder here than it would elsewhere: a pooled library pays titre loss as
representation.** Every bottleneck resamples 13,824 members, and a lower titre is a smaller
sample of them. Section 9.9 sets that bound and should be revisited once the slope is known.
**How much titre is open** — the abstract gives the curve's shape and not its slope, which is in
the paywalled figures.

#### Two findings recorded, neither of them settled

- **pLVX's LTRs appear not to be ΔU3 SIN.** #244 read an intact U3 enhancer-promoter in both
  copies — two NF-κB sites, an Sp1 GC-box and a TATA box — where a self-inactivating vector has
  U3 deleted from its 3' LTR. This section assumed a SIN parent, and the reading is against that
  assumption. It bears on copy-counting during domestication, on read-through, and on the
  poly(A) screen above. **Measured and open, not settled:** the read is of an on-screen
  reconstruction of the Addgene-verified sequence and may be an artifact of it. What would close
  it: the manufacturer's own map or sequence, or an assay for LTR-driven transcription from the
  integrated provirus.
- **The parent's length disagrees with itself by up to 2 bp.** The render gives **9,895**;
  Addgene's *Total vector size* field says **9,894**; its *Backbone size w/o insert* 9,173 plus
  *Insert Size* 720 gives **9,893**. Every length above is the render's, the only one of the
  three that is a sequence rather than a typed-in field. In a repo where every span is 0-based
  and half-open, a 2 bp disagreement is not cosmetic. **Closed by the login-gated GenBank file,
  which only the lab can pull.**

### 2.3 The working cassette

`[TATG]` — no N-terminal part — ccdB cassette — C-terminal part — `[CTAA]`.

The start codon is the `ATG` of `TATG`, supplied by the backbone's promoter. The "no N part"
version of the ccdB cassette is used, which adds one Gly.

The C-terminal part is `[T|TCC]─[6xHis, 18 bp]─[T2A, 54 bp]─[GG]`, 78 bp counting its 5'
overhang. The trailing `GG` plus the `C` of `CTAA` spells the final residue, and `TAA` stops
translation — which is how a part obeys the no-stop-codon rule and still ends an open reading
frame. 78 is a multiple of three, as the part rule requires (**derived**).

**The T2A does no work in this project.** Translation stops immediately after it. It is carried
because the part is a lab resource whose later uses put something after the tag, and because the
same part serves the marker-as-a-part placement section 7 rejects. Whether the project would
rather have an induction reporter there is in section 10.

#### Which way the cassette points

**The induced cassette reads on the minus strand**, measured on the parent: `TRE3GS` (3242-3610)
and `EGFP` (2513-3233) are both annotated reverse, with the SV40 poly(A) signal at 2186-2321
beyond them. The transactivator and puromycin cassettes read forward and end at the WPRE with no
terminator of their own.

So the one terminator on the plasmid points away from the packaged genome. That is what makes
requirement 3 pass in substance, and not merely on a narrowed screen: the only element whose job
is to stop transcription is pointed the other way. It is a property of the parent the lab chose,
not something the method arranged.

#### The stuffer is the parent's own EGFP

The method page asks for RFP. Section 7.1 chose a red transduction marker, and red against red
means colour cannot separate a vector whose cassette was never replaced from a cell that was
only transduced. The EGFP already in the parent is being cut out anyway, so it serves instead.

Measured across its 720 bp, both strands: **0 BsaI, 0 BsmBI, 0 BbsI, 0 PaqCI, 0 SapI**. It
translates `MVSKGEELFTGV…GMDELYK*`, one stop and it the last codon. It carries **no poly(A)
hexamer on either strand**, so retaining it reintroduces nothing requirement 3 is about. It is
already under the TRE promoter in the right orientation, so only its two flanks are replaced —
one PCR off the parent with tailed primers, against sourcing and inserting an RFP.

This is departure **D17** in `docs/research/synthesis-and-assembly-departures.md`.

**The caveat, because it changes what the stuffer is for.** The stuffer sits under a TRE
promoter, which is a minimal CMV promoter and is not read in *E. coli*. **No stuffer colour
screens a bacterial colony** — an RFP stuffer would not have screened one either. The colour
reads in the transduced cell after induction: green means the cassette was never replaced,
red-only means it was. A colony-level screen would need a bacterial promoter on the stuffer or a
counter-selection marker, and that is a separate decision. Nothing read measures TRE3GS activity
in *E. coli*.

## 3. The cargo layout, and what it costs the protein

Every part is one iGGA cargo block, synthesised from the oligo pool:

```text
[AGGA] [fragment] [internal stuffer, 34 bp] [barcode, 11 bp] [TTCC]
```

with the 34-mer `AGGAAAGTCTTCAGCCCGGGCAGAAGACAATTCC` from the method page. Part length is
97 bp + fragment, so **141 bp to 1,143 bp** across the 72 parts (**derived**).

After the two rounds of section 6 a library member reads:

```text
AGGA | Nfrag | AGGA | DBDfrag | AGGA | Cfrag | stuffer | BC_C | TTCC | BC_DBD | TTCC | BC_N | TTCC
```

so the barcode block reads C, DBD, N — the reverse of the order the rounds added them — with
`TTCC` as the cloning scar. Largest member: 2,251 bp of fragment plus 79 bp of tail plus the
leading overhang, about **2.33 kb** (**derived**).

### 3.1 The frame, and the one residue every junction charges

The frame is set upstream at `T|ATG`. The entry overhang reads `A|GGA`: its first base closes
the preceding codon and `GGA` is a whole codon. Each fragment therefore begins at a codon
boundary and its length is ≡ 2 mod 3, so a fragment holding *m* residues is **3m − 1 bp**
(**derived**). The source's own N and bZIP sequences measure the same congruence against their
own 4-bp scars, which is the cross-check; its C position differs only because T2A follows there.

Two consequences, and they are this project's answer to a gap the research note records as
"no rule fits a fragment's boundary codons to the entry overhang":

1. **One Gly is inserted at every junction.** `GGA` is a whole codon and is added, not
   substituted. Native residues are kept.
2. **The last residue of every fragment must have a codon ending in `A`.** Eight do not:
   **C, D, F, H, M, N, W, Y** (**derived** over the standard table). A fragment ending in one of
   those has that one residue replaced.

The replacement table this project uses, one residue each, chosen by chemical class:

| Last residue | Replaced by |
| --- | --- |
| C | S |
| D | E |
| F | L |
| H | Q |
| M | L |
| N | Q |
| W | L |
| Y | L |

The matrix that would justify each pick is **open** — none was read in this run, and the package
holds no substitution matrix. How many of the 72 parts need a replacement is computed at design
time and is a prediction section 8 tests.

**This is the equivalent of the source's residue changes.** The source changed one to two amino
acids per domain to standardise four overhangs; ours changes at most one per fragment and adds
one Gly per junction, because one overhang pair serves every position (departure D3). The
alternative — overwriting each fragment's first residue with Gly instead of inserting it,
which is what the source did at its `C|CGA` junction — buys back one residue of length and costs
a native residue, and length is not scarce here.

### 3.2 What the tail costs

Between the C domain and the His tag the construct translates the retained stuffer, three
barcodes and two scars: 34 + 11 + 4 + 11 + 4 + 11 = **75 bp, 25 residues** (**derived**). The
research note's figure of 45 bp counts the three barcode-and-scar units only and not the
retained stuffer; the two agree on the barcodes.

The stuffer translates `RKVFSPGRRQF` at this frame and spells no stop (**derived**). The
barcodes are drawn stop-free at the same frame, which is section 9.5.

**No capping block.** The last round keeps its stuffer, because the C position needs its barcode
and a capping block carries none. The working cassette supplies the tag and the stop instead.

### 3.3 Domestication

Every fragment is reverse-translated from its amino acid sequence with the `human` codon table —
the expression host is a human cell line — and then cleared of:

- the five reserved enzymes, as eight motifs counting both orientations: BsaI, BsmBI, BbsI,
  SrfI, PmeI;
- any of the three orthogonal primer sites assigned to its slot, or else the part moves slot,
  as the method page's step 5 allows;
- **any polyadenylation signal, `AATAAA` and `ATTAAA`, on the sense strand.** The cargo rides
  inside the packaged lentiviral transcript, so a terminator inside it truncates the genome.
  This is a motif class the package's domestication does not currently take, and the citation
  for the screen is **open** — the mechanism is standard lentiviral practice and no source for
  it was read in this run.

## 4. Cargo synthesis

One oligo pool holds all 72 parts. The oligo layout is section 9.1, the padding section 9.2, the
primer split section 9.3 and the batch size section 9.4.

### 4.1 The two inputs

Everything in this section is arithmetic over two numbers, and both are the project's rather
than the method's. Neither the oligo length nor the number of oligos is fixed by the method.

| Input | This project | Where it is set |
| --- | --- | --- |
| Oligo length, counting the whole oligo | **350 nt** | section 9.2 — a project decision, changed from 300 because the lab changed it |
| Blocks to be synthesised | **72**, 133 to 1,149 bp, 30,519 bp in all | the design the pipeline wrote, measured 2026-10-06 and tabled in `docs/research/synthesis-and-assembly.md` |

Change either and every number below changes with it. Nothing below is edited by hand.

### 4.2 Span per oligo — **derived**

The oligo is `[P1][BsmBI][span][BsmBI][pad][P2][P3]` (section 9.1). Its overhead:

| Part of the oligo | nt | Why |
| --- | --- | --- |
| P1, P2, P3 | 3 × 20 = **60** | three mutually orthogonal 20-mers, one per role (section 9.3) |
| Two recognition sites and spacers | 2 × (6 + 1) = **14** | 6 nt recognition and 1 nt spacer at each site |
| Overhead in all | **74** | |

**The quantity to budget is the templated span, which carries its own overhang at each end**, and
an overhang is 4 nt of the part rather than overhead. The two are easy to confuse and the
difference is a fragment:

```text
span  = oligo length − 60 nt of primer sites − 14 nt of recognition site and spacer
    at 300 nt:   300 − 74 = 226 nt, holding 218 nt between its overhangs
    at 350 nt:   350 − 74 = 276 nt, holding 268 nt between its overhangs
```

**None of the 74 scales with oligo length**, so the layout holds unchanged at any length.

**Each internal overhang belongs to two fragments at once**, so an oligo contributes
`span − 4` = **272 nt** of the part at 350, not 268. Summed over *n* pieces,
`sum(span) = L + 4(n − 1)`, which gives section 4.3's piece count. Counting the overhang as
overhead instead — the error an earlier revision of this section made — overstates the piece
count for a part whose length lands just past a multiple of the span.

276 nt is a ceiling, not a length every oligo spends. The padding of section 9.2 is the slack
between a real span and 276, and it sits outboard of the 3' cut.

The research note quotes 246 bp per 300-mer from Freschlin and `synthesis-and-assembly-departures.md`
D10 quotes about 235. Both describe a two-primer layout and neither is this project's budget.

### 4.3 Fragments per part — **derived**, and the reason to prefer 350

A block of *L* bp needs ⌈(*L* − 4) / (span − 4)⌉ fragments, by the sharing rule of section 4.2.
Counted over the 72 blocks of section 4.1:

| Fragments a part needs | At 226 nt span (300) | At 276 nt span (350) | Lund's measured success |
| --- | --- | --- | --- |
| 1 | 17 parts | 21 parts | not an assembly |
| 2 | 34 | 34 | **100%** |
| 3 | 6 | 6 | **93.8%** |
| 4 | 6 | 10 | not measured |
| 5 | 8 | 1 | **84.6%** |
| 6 | 1 | **none** | not measured |
| Worst part | **6 fragments** | **5 fragments** | |
| Oligos at the floor | **173** | **152** | |
| Oligos the design spends | **174** | **153** | |

**The two counts differ by one block, and the gap is real.** The row above it is the arithmetic
floor: the fewest oligos the length budget could ever need. A cut also has to spell a legal
overhang, and N_ATF7 is 1,092 bp, which is exactly 4 x 276 less the three shared overhangs. At
four fragments every piece is forced to 276 nt, so all three cut positions are forced, and the
first spells `GCCG` — one base kind, which is refused because such a junction truncates. No
legal set exists at four, so that block spends five.

Lund's is the only measured curve: of designs assembled from that many fragments, the share with
a perfect clone among four colonies is 100% at 2, 93.8% at 3, 84.6% at 5, 66.7% at 8, 40.0% at
12 and 0% at 16. It has no point at 4 or at 6, and nothing here interpolates one.

**The fragment count is the reason to prefer 350.** At 300 one part needs six fragments, past
every anchor Lund measured above 84.6%. At 350 no part exceeds five, so **all 72 sit at or above
Lund's 84.6% point**, and 55 of them need two fragments or fewer — the band Lund measured at
100%. Seven parts drop a fragment, and 21 oligos leave the pool.

Where the fragment count turns over, since the oligo length is an input (**derived**, over the
same 72 blocks): any oligo of **307 nt or more** holds every part to five fragments. 350 clears
that by 43 nt. Holding every part to four would need an oligo of **365 nt**, which is past the
top of the price band (section 4.5).

### 4.4 Oligo count — **derived**

Summing the fragments per part: **152 oligos at 350 nt** at the floor and **153 as designed**,
against 173 and 174 at 300. One part takes one to five oligos. These are counted over the blocks the design actually wrote rather than
estimated, which is what #220 was asked for and has now done.

All four counts in sections 4.3 and 4.4 moved when #257 corrected the span arithmetic. They are
computed by the cargo designer, not typed here; #265 is where that lands.

### 4.5 Cost — the two bands are close, and the pilot's own cost is still **open**

Every figure here is a **list price with its date**. A list price is a reference, not a fact
about what this project would pay. The two recorded, from the vendor's oligo-pool price table
captured **2026-09-17**:

| Pool | Length band | List price | Per oligo |
| --- | --- | --- | --- |
| 18,000 oligos | 251-300 nt | $10,004 | $0.556 |
| 18,000 oligos | 301-350 nt | $12,505 | $0.695 |

Per usable base, at this project's own span rather than Freschlin's (**derived**). An oligo
delivers `span − 4` bases of the part, by section 4.2:

```text
at 300 nt:   $10,004 / 18,000 / 222 nt = $0.00250 per usable base   =  $2.50 per usable kb
at 350 nt:   $12,505 / 18,000 / 272 nt = $0.00255 per usable base   =  $2.55 per usable kb
```

**The two bands are within 2.0% of each other per usable base**, so 350 is not the dearer
choice. The 301-350 band lists at exactly 1.25× the 251-300 band for the same 18,000 oligos, and
272 / 222 = 1.225× buys almost all of that back. "350 costs 25% more" is true of the pool and
false of the sequence.

**The pilot's own cost is not known, and no tier figure stands in for it.** 152 oligos sit far
below 18,000, the smallest pool size any price recorded here covers, so **this project's
synthesis cost is open** — the figures above price a pool it is not ordering. The project
exercises the design path, not the pool economics, and should share a pool with other designs if
one is going out.

**Two things to flag rather than bury:**

- **350 nt is the top of its price band.** One base past it, 351 nt, is a different band whose
  price is not recorded here. The design has **no slack above 350 at all**, so anything that
  lengthens the oligo — a fourth primer role, a longer recognition site, a wider overhang — is a
  tier change and not a rounding error. Section 9.2's padding rule is what holds the line.
- **The length is the lab's, not the method's.** Section 9.2 reads 350 because the lab moved it
  from 300. The method page fixes neither number.

## 5. Validation of the parts

Route: **index PCR**, the method page's Route B.

- Four colonies picked per part, which is Lund's anchor and the only one measured: 343 of 458
  genes error-free at four colonies.
- 72 parts × 4 = **288 wells**, which is one 384-well pick plate and three 96-well index plates
  (**derived**). Picking fills a quarter of the pick plate at a time, which is what gives three
  full index plates rather than four part-filled ones.
- Nothing picks the route (#299): the project names it, and this demo ships a project file for
  each so both are exercised.
- Well-to-barcode mapping is the published LevSeq plate map, ordered as given, recorded per
  plate.
- The pass criterion is section 9.6.

The finished library is **not** read per well. It is read for linkage once and for
representation at each bottleneck (section 6).

## 6. The rounds

**Three positions, two iGGA rounds.** The issue says three rounds; the third is the step that
seats position one, and it is not an iGGA round. Spelled out, because #220 has to build the
right number of reactions:

| Step | Destination | Donor | Enzymes | Distinct products |
| --- | --- | --- | --- | --- |
| Seat every part | DMX vector | all 72 parts, one design per well | BsmBI | 72 designs |
| Round 1 | N part list | DBD part list | BbsI + SrfI / BsaI + PmeI | 576 |
| Round 2 | N+DBD library | C part list | BbsI + SrfI / BsaI + PmeI | 13,824 |
| Final assembly | working vector | the finished library | PaqCI, one pot | 13,824 |

All 72 parts are cloned into the DMX vector at the seating step, because our donor carries no
release sites of its own and cannot be a donor until it is cloned (the constraint recorded under
open decision 6.7).

**The pipeline counts differently, and it is right to.** `liulab_synbio.igga` models one round
a position, so it plans three for this project. Measured on 2026-10-06: its first round opens the
destination with the internal enzyme, ligates the N part list released by the external enzyme,
transforms it and sizes it for 24 products — the same reaction as the two after it. That is not
the seating step, which runs BsmBI one well a part and makes no library. The pipeline has no
carrier to seat a part in, so this table is not reachable through it; the gap is the fifth item
of #223.

After round 2: the linkage read, then the representation read. Representation is read again
after the final assembly, after packaging and after transduction — each is a bottleneck that
resamples the library.

Donor pools, clean-up and the per-round bound are sections 9.7, 9.8 and 9.9.

## 7. The transduction marker

A constitutive mCherry, reporting that a cell was transduced whether or not the cargo is induced.

### 7.1 The choice: a backbone cassette, beside the transactivator

**Chosen.** mCherry is expressed from the same constitutive promoter as the transactivator,
inside the backbone, upstream of the TRE promoter and its assembly cassette.

Why:

- **The source now exists, and it argues against the alternative.** Putting the marker in a part
  means a second internal promoter inside the packaged transcript. When this was first written
  the cost of that was uncited, so the placement was rejected as not committable; #231 has since
  measured it, and the measurement is a cost **against** a second internal promoter — Curtin et
  al. 2008 find interference between two internal promoters running both ways and large. The
  decision stands, on better grounds than it was made on.
- **It needs no new part class.** The marker-as-a-part placement needs the terminal part exempted
  from the no-stop-codon and in-frame rules, which is a change to the method, not a choice
  inside it. The first project to exercise the method should not also be the first to bend it.
- **Swappability buys nothing here.** The argument for the part placement is that the marker
  swaps like any other part. This project has one marker and never swaps it. The cost of the
  backbone placement — rebuilding the backbone to change the marker — is paid once by a project
  that does not change it.

What the choice costs, stated plainly: the marker is fixed in the vector. A later project that
wants puromycin resistance or a blue marker rebuilds the backbone, or revisits this decision
with a source in hand.

### 7.2 What the rejected placement would have required

Recorded so the decision can be reversed cheaply. Were the marker a C-terminal part carrying its
own cassette, the specification would have to add: the terminal part is exempt from the
no-stop-codon and in-frame rules — a new part class or a stated exception for the last position
— and the cassette must carry no internal polyadenylation signal, since an internal promoter in
a lentivirus works by relying on the 3' LTR and a terminator would truncate the packaged genome.

**That last premise is defeated, and recorded here so it is not reused.** Tian and Andreadis
2009 clone both internal cassettes **antisense to the LTRs**, which puts an internal
polyadenylation signal on the strand the producer cell does not package. A terminator is
therefore possible in a lentiviral vector. It is not free — see section 7.3 — but "impossible"
was wrong.

### 7.3 Read-through is measured; the number this project wants is not

Most of this is closed by #231; `docs/research/promoter-readthrough-lentivirus.md` has the
detail.

Measured: transcription from an upstream internal promoter reaches the downstream unit in a SIN
lentiviral vector and contributes to its output (Tian and Andreadis 2009); interference between
two internal promoters runs **both ways** and is large (Curtin et al. 2008); and which upstream
constitutive promoter drives the regulator changes how much a downstream inducible promoter
leaks, and not through the regulator's own level (Benabdellah et al. 2016).

Still unmeasured by anyone: the uninduced background of a TRE promoter downstream of a
constitutive one, against a matched control with the upstream promoter removed.

**#256 decided not to wait for it.** Eszterhas et al. 2002 is why: the magnitude *and the sign*
of interference move with integration site, so in a pooled library with random integration leak
is a distribution, not a number, and no option makes it a constant. The read-out has to survive
the spread either way. So this project carries **a no-doxycycline arm for every screen**, and
reads each member against its own uninduced well rather than a pooled baseline. That is a
requirement here, not an option.

Three mitigations stay reachable, with what they cost, should the first screen show leak
dominating:

| Mitigation | Measured | Cost |
| --- | --- | --- |
| Is2 insulator in the 3' LTR | 293T leak 66.2 → 7.4, induction 5.7 → 38.8 | **titre down 2–3 fold**, and measured with a TetR repressor, not an rtTA activator |
| Both cassettes antisense to the LTRs, terminator between | Restores the suppressed unit, cuts read-through | Never measured with an inducible downstream promoter; rebuilds the vector |
| Weaken the upstream constitutive promoter | hEF-1α for SFFV lowered leak | Lowers induced output with it |

Section 5 of #231's note gives the experiment that would close the number, if it is ever worth
running.

## 8. The comparison against the source

The source's own designed parts are the **check, not the input**. They are opened only after
this project's design is frozen and written to disk, and the comparison is run against that
frozen output. **A design re-run after the source has been read is a new design and says so.**

The comparison runs at four levels, and agreement means something different at each.

| Level | Agreement means | Identity expected? |
| --- | --- | --- |
| Protein | our per-fragment protein equals the input, except where section 3.1 forces one replacement | yes, except at boundaries |
| Boundary residues | every residue we changed is forced by our rule, and every residue they changed is forced by theirs | no — the two rules differ |
| DNA | both designs are frame-correct, stop-free and clear of their own reserved enzymes | **no** — two independent codon choices differ at most wobble positions |
| Barcodes | both sets meet their own stated rules | **no** — six of the source's 72 carry a stop at our frame offset, so the sets must differ |

Overhangs, stuffers, the enzyme assignment and the absence of external stuffers are already
accounted for by departures D1 to D15 in
`docs/research/synthesis-and-assembly-departures.md`. A difference on that list needs nothing
but its citation.

**What a difference obliges, by kind:**

1. **Already a recorded departure.** Cite it. Nothing else.
2. **A property both designs claim and ours fails** — a reserved enzyme site we left in, a frame
   that does not close, a stop inside a barcode. Our rule is wrong. File an issue, fix the rule,
   re-run the design, and say in the issue that the source found it.
3. **They did something our rules do not require.** Establish which rule forces it. If their own
   overhang standard forces it and ours does not, record it as a departure and move on. If it is
   forced by something our rules also face and we missed, it is case 2.
4. **We did something they did not.** Same test, the other way round.
5. **Unexplained either way.** It is open. It goes into a research note, not into code, and no
   package default changes until it is closed.

**The one difference already waiting:** the source changed one to two amino acids per domain to
standardise four overhangs — measured, 15 of 24 N parts have a changed final residue pair and 19
of 24 C parts a changed first residue, every synthesised C domain starting with R. Our equivalent
is section 3.1. Agreement here means our changed-residue set is explained by our one overhang
pair and theirs by their four, with the junction count matching the positions. A difference where
we changed nothing and no rule of ours covers the boundary is case 2, and the strongest single
result this project can produce.

## 9. Every deferred decision, answered for this project

Each heading names the open decision in `docs/research/synthesis-and-assembly.md` section 6.
**Each answer is this project's, not the method's**, except where a later ticket promoted one.
Ticket #258 took 9.3 and 9.4 into the method, so the method page now carries those two.

### 9.1 Where the Type IIS sites sit on the oligo — on the oligo, inboard of P1 and P2

**#257 promoted this to the method, so this section now reports rather than decides.** Every
oligo reads `[P1][BsmBI][span][BsmBI][pad][P2][P3]`, with both sites inward-facing and
templated. The two outermost cuts across a part yield `AGGA` and `TTCC`; every internal cut
yields the overhang the fragment split chose. One rule serves the vector interface and the
internal junctions alike.

Why not on the primers: all of a gene's fragments share one P1 and P2, so a primer-borne site
would give every fragment of that gene the same overhang and no internal junction could be
expressed. Oligo length was never the deciding term — #257 measured the saving at two oligos
across this project's 72 blocks. This is the Baker diagram's arrangement, which draws the cuts
on the oligo, inboard of P1 and P2.

### 9.2 The padding rule — pad every oligo to a uniform 350 nt

**#257 promoted the rule to the method; the length stays this project's.** Filler sits between
the 3' BsmBI recognition site and P2, outboard of the cut, so it never enters the product. It is screened for the eight enzyme motifs, for the orthogonal primer sites,
and at its two junctions for a motif the join creates.

Why uniform: Freschlin — which pads to a uniform length — lost whole replicates to an enzyme
site inside a payload, which is what the screen is for, and one length makes a pool-uniformity
check unfailable. Not the vendor's 15% uniformity rule: that belongs to Multiplexed Gene
Fragments, a different product, and no vendor publishes a length-spread rule for an ssDNA oligo
pool. The rule was argued against that floor in #261 — the widest within-part spread, 16%,
measured against a 15% gate that never applied — and the conclusion is unchanged, since every
oligo is padded to one length anyway. Only its justification was borrowed. Binning by length
buys nothing at 152 oligos.

**Why 350 and not 300: the lab changed it.** This is a project decision and not a method
constant, and it is the only input section 4 takes besides the block list. What it buys is
measured in section 4.3 — no part needs more than five fragments at 350, where one needs six at
300 — and what it costs is measured in section 4.5, which is 2.0% per usable base. 350 is the
top of its price band, so the length cannot drift upward without crossing a tier.

### 9.3 The orthogonal primer split — 96 inner, 35 and 34 outer

| Role | Primers | What it indexes |
| --- | --- | --- |
| P2, inner reverse | **96** | one gene inside a batch |
| P1, shared forward | **35** | the batch, with P3 |
| P3, outer reverse | **34** | the batch, with P1 |

96 inner primers makes one batch exactly one 96-well PCR2 plate, which is the Baker shape and
the reason the number is 96 rather than any other. Capacity is 35 × 34 = 1,190 batches of 96
genes (**derived**), which stops being the limit long before the reaction count does.

The split totals the **165** rows the primer supplement keeps. The abstract says 166; that
discrepancy is unresolved in the research note and this project uses the 165 sequences held.

### 9.4 Batch size — one batch, all 72 parts

One PCR1, then 72 PCR2 reactions using 72 of the 96 inner primers. **152 pieces in the PCR1
tube** — section 4.4 derives the count, and this section does not restate it — which is at the
low end of the Baker anchor of hundreds of oligos to a well.

The method's rule is to hold pieces per PCR1 roughly constant; **this project sets that constant
at whatever section 4.4 derives**, which at 350 nt is 152 pieces. Freschlin's evenness
measurement — subpools under 16 genes overabundant,
20 or more underrepresented — compares subpools against each other inside one pool, and a single
batch has none to compare against. If evenness turns out to limit anything, the fallback is four
equal batches of 18, since equal sizes are what that measurement actually argues for.

### 9.5 The barcode set's enzyme-site freedom — on, with the rest of the rules from the barcode note

Eleven bases, and the set is drawn fresh:

| Rule | This project |
| --- | --- |
| Distance metric | **Sequence-Levenshtein** |
| Minimum distance | **3**, within each part list |
| Enzyme-site freedom | **on** — all eight motifs, both orientations, checked in context with the flanking scar and stuffer |
| In-frame stop | none, at this layout's frame offset |
| GC band | **none** |
| Homopolymer cap | longest run ≤ **5** |
| Draw order | seeded shuffle, reproducible from the seed |

Site freedom is on because it is measured: an 11-mer drawn to the distance rules alone carries
one of the eight motifs 0.90% of the time, which is one library member in 112 cut apart in its
own round, and Freschlin lost whole replicates to exactly that.

The metric, the absent GC band, the cap of 5 and the shuffled draw order are the verdicts of
`docs/research/barcode-design.md`, and the method page now carries them. Enzyme-site freedom was
the method page's one open barcode decision; #260 closed it **on**, so every row above is now the
method's and this project restates rather than chooses them.

The source's 72 barcodes are not reused. The reason is the distance rule and the site screen
above, and **not** a stop codon: an earlier revision of this section said six of the 72 carry one
at this design's frame offset, and that is wrong.

**Measured 2026-10-06**, twice. The barcodes were read off each sheet's own barcode column of
Table S1 — column 6 on the C sheet, column 5 on the N and bZIP sheets — and the frame was read
off `docs/examples/ap1-library/product.dna`, the record the pipeline wrote. Every barcode sits at
**codon position 1**: the retained 34-base stuffer starts at a codon boundary, translates
`RKVFSPGRRQF` and leaves one base over. So the codons wholly inside a barcode are its bases 2-4,
5-7 and 8-10, and the codon spanning into it opens with the scar's last base — `C` for ours,
`G` for theirs — which no stop codon does. **None of the 72 carries a stop**, with our scar or
theirs, in a bare barcode-plus-scar unit or in a whole 41-base block. Translating the product's
whole reading frame through the tail finds no stop either.

Six is what the same 72 give at **codon position 0**, which is where this design does not put
them. #220 measured none independently on the same day, and #224 re-measured it.

### 9.6 What counts as a passing well — LevSeq's criterion, plus an exact match

**#260 made both halves of this the method's rule, so this section now reports rather than
decides.** The depth floor travels with the route: the read-out is Route B, so the floor is
LevSeq's — twenty reads wanted and ten tolerable, a per-base binomial test at a 5% false
discovery rate, and a well with more than one significant mutation called mixed.

**On top of it, a well passes only on an exact match to the designed part** across the whole
designed region — both entry overhangs, the fragment, the stuffer and the barcode. A silent
mismatch inside a barcode mislabels a library member for the rest of the project, and the
linkage read cannot recover what the barcode no longer names. Mixed wells fail. A well below the
floor is **low**, which is no verdict rather than a failure: it is read again or picked again,
and reformatting does not compact it out.

Four colonies per part are picked (section 5); a part with no passing well is re-picked from the
same archive spot before it is re-synthesised.

### 9.7 How donor pools are built — once, up front, from validated wells

All 72 parts are validated in one batch, so each position's passing wells are combined into one
equimolar pool before round 1, and the same pool is re-digested whenever that position is the
donor. Only the destination changes between rounds.

This is the source's structure with our per-well validation in front of it. Building each
round's pool from that round's wells would buy nothing, because there is only one validation
batch and each position is a donor exactly once.

### 9.8 Clean-up — the source's two ratios

SPRI at **2×** after each digest and **1×** after each ligation. Not one ratio used twice: the
2× step is what removes the 17 bp and 13 bp cut stubs, and the whole PmeI argument in the
departures note rests on it. One ratio throughout is supported by nothing read.

### 9.9 The per-round bound — mass in, colonies counted out

Both, because they bound different things.

- **In:** the source's numbers — 20 ng of digested backbone per 200 µL ligation, up to 100 ng of
  purified ligation product electroporated.
- **Out:** a plated dilution from each round's recovery, counted, and the round carried forward
  only when the colony count reaches the floor for that round's own complexity: 183 colonies for
  the 24 products of round 1, 6,306 for 576 after round 2, and 195,386 for 13,824 after round 3
  and again at the final assembly.

Each round also plates a no-donor control, carried through the ligation from the same digest,
and the gate is net colonies. **Both plates are a departure**, settled in #259: the source plates
nothing at any round. The reason is there, not here.

**The coverage multiple is closed: this project states a completeness of 0.99, not a multiple
— decided in #297.** 300x was Qian's, measured at one transformation whose purpose was pickable
clones for an arrayed collection. This library is never picked: it goes to a pooled screen as a
pool, and no member is re-identified. At 300x the three rounds would have asked for 7,200, then
172,800, then 4,147,200 colonies, and nothing read here justified any of the three.

A multiple was the wrong parameter and not only the wrong value: one ratio against each round's
own products buys a different risk every round. The completeness is what holds across them, so
`liulab_mbio.bench.coverage` computes each round's floor from it — 183, 6,306 and 195,386
colonies, which no whole multiple reproduces. It is a **floor, not a sufficiency claim**: it
assumes every member equally represented, and synthesis skew breaks that. What the screen
downstream needs is still open, and stays in section 10.

### 9.10 Also open, answered here

- **Read-out route:** index PCR for the parts (section 5); representation reads, not per-well
  reads, for the library.
- **DMX parent:** a method decision, not this project's (section 2.1).
- **Colony picker model:** a purchase decision, untouched.
- **Library-read primers:** decided in #258, and not by this project. Two pairs against
  sequence already constant: representation reads from the retained stuffer to the vector past
  the final `TTCC`, linkage reads the vector either side of the whole cargo. The package designs
  both against the simulated record. The two universal flanking primers turned out to have
  sequences after all — they are the barcode kit's two unmatched constant regions.

## 10. Open, and what would close it

Nothing here becomes a package default, and nothing here is guessed at.

- **The parent's own length, by up to 2 bp.** 9,895 from the render, 9,894 from Addgene's size
  field, 9,893 from its backbone-plus-insert fields (section 2.2). Every span here is 0-based
  and half-open, so this is not cosmetic. **Closed only by the login-gated GenBank file, which
  the lab can pull and this run cannot.**
- **Whether the two in-scope poly(A) hexamers are functional**, `ATTAAA` at 880 and `AATAAA` at
  1603 (section 2.2). Closed by a source measuring 3'-end processing in a transfer vector's
  gag/RRE segment, or by a full-length-genome measurement on this backbone.
- **Whether the hPGK promoter tolerates a change at 3724 or 3876**, the two sites domestication
  would touch there. hPGK drives the transactivator, a base change in a promoter is synonymous
  in no sense the package can check, and nothing read measures either position. The decision
  built on it is #256's.
- **Whether pLVX is self-inactivating.** #244 read an intact U3 in both LTRs (section 2.2).
  Measured and open: closed by the manufacturer's own sequence, or by an assay for LTR-driven
  transcription from the integrated provirus.
- **How much titre the cargo costs.** Kumar et al. give the curve's shape and not its slope, and
  the figures are paywalled. A pooled library pays it as representation (section 2.2).
- **Read-through and promoter interference.** No source read. Closed by one primary source,
  before any project uses the marker-as-a-part placement.
- **The polyadenylation screen.** The motif list `AATAAA`/`ATTAAA` is standard practice and
  carries no citation here. Closed by a source, and by domestication learning to take a motif
  list rather than only an enzyme list.
- **The substitution table in section 3.1.** Chosen by chemical class, justified by no matrix.
  Closed by naming the matrix and re-deriving the eight picks from it.
- **What the screen downstream needs.** The colony floor is settled (#297): this project states
  a 0.99 completeness and the package sizes every round for it. What representation a screen of
  this library would ask for is not settled, and a floor is not an answer to it. Closed by naming
  the screen, or by a source measuring diversity loss across an iGGA round.
- **The synthesis cost.** 152 oligos sit far below 18,000, the smallest pool any recorded price
  covers, so no price here is this project's. The list prices and the per-usable-base comparison
  are section 4.5; neither stands in for the pilot's own cost.
- **The orthogonal set's size**, 165 or 166. The 165 held are what the split uses.
- **The MOI and the cell number at transduction.** The completeness floor gives the number of
  integrants wanted; the MOI that delivers them at mostly one integrant per cell is not fixed
  here.
- **Whether the T2A should do work.** As specified it is followed immediately by a stop, so it is
  carried and unused (section 2.3). If an induction reporter downstream of the cargo is wanted,
  that is a different C-terminal part and a different decision.
- **Size skew.** Library members span about 0.5 to 2.3 kb, and every bottleneck can favour the
  short ones. The representation read measures it; no magnitude is asserted.

## Sources

Every measurement cited above is recorded in one of these, and none is recomputed here.

- `docs/research/synthesis-and-assembly.md` — the method's sources, its departures, the review's
  measurements, and section 6, the open decisions this note answers.
- `docs/research/protein-library-assembly.md` — Takacsi-Nagy et al. 2026, its method, its
  measured sequences and its CC BY 4.0 licence, which is what lets the AP-1 sequences be used.
- `docs/research/synthesis-and-assembly-departures.md` — departures D1 to D15, worked base by
  base.
- `docs/research/working-vector-plvx-tetone.md` — pLVX-TetOne-Puro-GFP measured against the four
  requirements of section 2.2, with Kumar et al. 2001 on packaging length.
- `docs/research/lentiviral-tolerance.md` — what a lentiviral vector tolerates, and the reading
  that pLVX's LTRs are not ΔU3 SIN.
- `docs/research/barcode-design.md` — the distance metric, the absent GC band, the homopolymer
  cap and the draw order.
- `docs/research/codon-usage.md` — the `human` table, counted from hg38 through liulab-genome.
- `docs/synthesis-and-assembly.md` — the method itself, which keeps every question above open.

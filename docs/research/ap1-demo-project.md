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
layout in section 3 and are marked **derived** where they appear.

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

The lab-resources build, unchanged: a rebuilt DMX parent, KanR, BbsI-free, one PmeI site
outboard of each BsaI site. Which parent is a method decision and not this project's; the
measurement that bears on it — inserting PmeI outboard of the right-hand cassette BsaI lands in
or abuts the `6xHis` coding sequence of the first parent — is in the research note, and this
project uses whatever the lab-resources build settles on.

### 2.2 The working vector: a Tet-on lentiviral backbone

A third-generation lentiviral transfer plasmid carrying the assembly cassette
`[TATG.BsmBI]─[RFP stuffer]─[BsmBI.CTAA]` under a TRE promoter, with the transactivator and the
transduction marker expressed constitutively from a second promoter between the LTRs.

Requirements the backbone has to meet, all checkable on its own map:

- No BsaI and no BsmBI site outside the designed cassette; domesticate any that are there.
- AmpR or CarbR, because the DMX backbone the cargo leaves is KanR and the last transfer must
  not share a marker (departure D11).
- No polyadenylation signal anywhere between the 5' LTR and the 3' LTR, the stuffer included.
  The genome is packaged as one transcript from the 5' LTR, so an internal terminator truncates
  it.
- Room for a cargo of up to about 2.4 kb (**derived**, section 3) plus the working cassette.

**No parent is fixed here.** No lentiviral map was read in this run, so naming one would be an
uncited choice. #220 picks the parent from a real map against the four requirements
above and records which.

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

Payload per 300 nt oligo: 300 − 60 nt of primer sites − 22 nt of Type IIS site and overhang =
**218 nt** (**derived**). This is lower than the 246 bp per 300-mer the research note quotes
from Freschlin, because that figure is measured on a design carrying two primer sites and ours
carries three. The research note's per-kb cost keeps its own denominator.

Oligos needed: 30,087 bp of designed sequence over 72 parts at 218 nt of payload each, with
per-part rounding, gives **roughly 150 to 210 oligos** (**derived**, an estimate; #220 computes
the exact count from the real sequences). Each part takes one to six oligos.

**The pool is far below the smallest tier priced.** The only prices recorded are for 18,000
oligos — $10,004 at 251-300 nt — so this project's synthesis cost is **open**, and ~200 oligos
in an 18,000-oligo pool is mostly empty space. The project exercises the design path, not the
pool economics, and should share a pool with other designs if one is going out.

## 5. Validation of the parts

Route: **index PCR**, the method page's Route B.

- Four colonies picked per part, which is Lund's anchor and the only one measured: 343 of 458
  genes error-free at four colonies.
- 72 parts × 4 = **288 wells**, which is one 384-well pick plate and three 96-well index plates.
  288 ≤ 384, so the route rule — at or below one plate of samples, index PCR — puts this project
  on index PCR (**derived**).
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
| Final assembly | working vector | the finished library | BsaI, one pot | 13,824 |

All 72 parts are cloned into the DMX vector at the seating step, because our donor carries no
release sites of its own and cannot be a donor until it is cloned (the constraint recorded under
open decision 6.7).

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

- **The alternative is not committable in this run.** Putting the marker in a part means a
  second internal promoter inside the packaged transcript, and the read-through and
  promoter-interference cost of that is a known trade that needs a source before it is
  committed. No source for it was read here and this run downloads none, so choosing it would
  put an uncited claim in the first project built against the method.
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

### 7.3 The read-through question is open

Both placements put two promoters between the LTRs, so read-through from the upstream promoter
into the downstream one, and the interference that can follow, is a question either way. The
backbone placement is the arrangement with a track record; the part placement is not. **No
source read in this run measures read-through or promoter interference for either.** It is
recorded as open, not answered. What would close it: one primary source measuring expression
from an internal constitutive promoter downstream of an induced one in a lentiviral vector,
read into `docs/research/`, before any project commits to the part placement.

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
**Each answer is this project's, not the method's.** The method page keeps carrying the question.

### 9.1 Where the Type IIS sites sit on the oligo — on the oligo, inboard of P1 and P2

Every oligo reads `[P1][BsmBI][payload][BsmBI][pad][P2][P3]`, with both sites inward-facing and
templated. The two outermost cuts across a part yield `AGGA` and `TTCC`; every internal cut
yields the overhang the fragment split chose. One rule serves the vector interface and the
internal junctions alike.

Why not on the primers: the inner primers come from an orthogonal set carrying no site, and
retailoring that set or exempting one role costs more than 22 nt of oligo. This is the Baker
diagram's arrangement, which puts the cuts inboard of P1 and P2.

### 9.2 The padding rule — pad every oligo to a uniform 300 nt

Filler sits between the 3' BsmBI recognition site and P2, outboard of the cut, so it never
enters the product. It is screened for the eight enzyme motifs, for the orthogonal primer sites,
and at its two junctions for a motif the join creates.

Why: uniform length removes the vendor's 15% uniformity question rather than managing it, and
Freschlin — which pads to a uniform 300 nt — lost whole replicates to an enzyme site inside a
payload, which is what the screen is for. Binning by length buys nothing at ~200 oligos.

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

One PCR1, then 72 PCR2 reactions using 72 of the 96 inner primers. About 200 pieces in the PCR1
tube, which is inside the Baker anchor of hundreds of oligos to a well.

The method's rule is to hold pieces per PCR1 roughly constant; **this project sets that constant
at about 200 pieces**. Freschlin's evenness measurement — subpools under 16 genes overabundant,
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
`docs/research/barcode-design.md`, and the method page now carries them. Enzyme-site freedom is
the method page's one open barcode decision, which is why this section decides it for this
project.

The source's 72 barcodes are not reused: six of them carry a stop at this design's frame offset,
measured.

### 9.6 What counts as a passing well — LevSeq's criterion, plus an exact match

The read-out is LevSeq's, so the criterion is LevSeq's: twenty reads wanted and ten tolerable,
wells below that marked low, a per-base binomial test at a 5% false discovery rate, and a well
with more than one significant mutation called mixed.

**On top of it, this project requires an exact match to the designed part** across the whole
designed region — both entry overhangs, the fragment, the stuffer and the barcode. A silent
mismatch inside a barcode mislabels a library member for the rest of the project, and the
linkage read cannot recover what the barcode no longer names. Low and mixed wells fail.

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
  only when the colony count reaches the coverage target for that round's own complexity: 576
  distinct products after round 1, 13,824 after round 2 and at the final assembly.

**This reverses the method page's "rounds run unchecked" for this project.** The research note
already records that sentence as reading a silence as a statement and as contradicting the
representation read two steps later, and that the mechanism it rests on is cited to no source
read here. A dilution plate is cheap; losing a round's diversity is not.

The coverage **multiple** is where this answer is weakest. The only gate read is Qian's, at its
one transformation, stated in the research note as a colony count against library complexity
above 300. This project inherits that multiple and marks it as inherited: it was measured for a
single transformation of a different library, not for an iGGA round. #220 should read the exact
form from the note before coding it, and the multiple belongs in section 10 as much as here.

### 9.10 Also open, answered here

- **Read-out route:** index PCR for the parts (section 5); representation reads, not per-well
  reads, for the library.
- **DMX parent:** a method decision, not this project's (section 2.1).
- **Colony picker model:** a purchase decision, untouched.
- **Library-read primers:** open. The method's primer inventory has no pair that reads an
  assembled library's barcode block, and the two universal flanking primers have **no published
  sequence**. A pair must be designed against the cargo's own constant flanks — the entry
  overhangs and the retained stuffer are constant across all 13,824 members, so one exists. This
  project does not design it; #220 does, or it is designed as its own ticket.

## 10. Open, and what would close it

Nothing here becomes a package default, and nothing here is guessed at.

- **The lentiviral parent.** No map read. Closed by reading one against the four requirements in
  section 2.2.
- **Read-through and promoter interference.** No source read. Closed by one primary source,
  before any project uses the marker-as-a-part placement.
- **The polyadenylation screen.** The motif list `AATAAA`/`ATTAAA` is standard practice and
  carries no citation here. Closed by a source, and by domestication learning to take a motif
  list rather than only an enzyme list.
- **The substitution table in section 3.1.** Chosen by chemical class, justified by no matrix.
  Closed by naming the matrix and re-deriving the eight picks from it.
- **The coverage multiple in section 9.9.** Inherited from a single-transformation gate in
  another method.
- **The synthesis cost.** The pool is far below the smallest tier priced; no price is recorded
  for a ~200-oligo pool.
- **The orthogonal set's size**, 165 or 166. The 165 held are what the split uses.
- **The MOI and the cell number at transduction.** The coverage rule gives the number of
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
- `docs/research/barcode-design.md` — the distance metric, the absent GC band, the homopolymer
  cap and the draw order.
- `docs/research/codon-usage.md` — the `human` table, counted from hg38 through liulab-genome.
- `docs/synthesis-and-assembly.md` — the method itself, which keeps every question above open.

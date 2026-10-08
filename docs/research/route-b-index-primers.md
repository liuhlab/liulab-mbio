---
search:
  exclude: true
---

# Route B's index primers: where they bind, and what marks a well

Research note for issue #312, under #225. It answers the question
`docs/research/route-b-index-pcr.md` left open — "**Primers for our vector.** None exist." —
and the answer is that the premise was wrong. Everything below was retrieved or measured on
**2026-10-06**.

## The short answer

Route B does not read the working vector. It reads the **DMX vector**, which is a pET-derived
backbone, and **LevSeq's own published annealing regions bind it verbatim and uniquely**, at
`[245, 276)` forward and `[883, 912)` reverse on DMX0001. Nothing has to be designed.

What does have to be decided is which 192 marks go on the 5' ends of those two regions. The lab
head's steer — Lund et al., NEB index sequences, dual-index oligos held as lab stock — is
**correct and has direct precedent in a paper this repository already cites**: Lund et al. 2024
appended NEBNext-for-Illumina barcode sequences to vector-annealing primers and read the
amplicons on Oxford Nanopore. There is no platform mismatch, because an index sequence is read
out of the sequence by software, not by the flow cell.

## How to read this note

Every coordinate is the package's own: 0-based, half-open, measured by running
`liulab_mbio.io.read_record`, `liulab_mbio.sites.find_sites` and
`liulab_mbio.primers.placement.find_binding_sites` over the Addgene record this repository
already holds. A row marked **derived** is arithmetic on measured or quoted values, and the
arithmetic is shown. Nothing here is from memory.

| Route | Works | Used for |
| --- | --- | --- |
| `reference_docs/.../dmx/addgene/addgene-plasmid-247434-sequence-494299.dna`, read by the package | yes | every coordinate and every length below |
| `reference_docs/.../LevSeq/papers/sb4c00625_si_001.txt`, the `pdftotext -layout` dump | yes | the two annealing regions, as printed in the SI's "Primer Design" |
| `reference_docs/.../long_fragment_GGA/Lund2024/main.txt` | yes | Lund's §4.3 "Isolate Sequencing Confirmation" |
| `reference_docs/.../long_fragment_GGA/Lund2024/sb3c00694_suppl-data/sb3c00694_si_005.dat` | yes | Lund's own cPCR annealing regions, as annotated features |

## 1. Sources and licences

| Document or dataset | Licence | Can we ship it? |
| --- | --- | --- |
| Long et al. 2025, *ACS Synth. Biol.* 14(1), 230–238 | © ACS; Europe PMC reports `isOpenAccess: N` | **No** — cite, never mirror. Verdict carried from `route-b-index-pcr.md` §1 |
| LevSeq Supporting Information `sb4c00625_si_001.pdf` | CC BY-NC 4.0 | **Not as a file or a table.** Non-commercial clashes with this repository's MIT licence. Single facts are cited and re-entered by hand |
| The two LevSeq annealing regions quoted in §2 | **a DNA sequence carries no licence** | Yes, as sequence. They are 31 and 29 bases printed in the SI's running text, re-entered by hand; Tables S3 and S4 are not reproduced and are not needed |
| Lund et al. 2024, *ACS Synth. Biol.* 13, 745–751, and its SI | © ACS; all authors were NEB employees and NEB funded the work | **No** — cite, never mirror. Already the verdict in `docs/research/dad-split.md` |
| Lund's vector records `sb3c00694_si_004.dat`, `si_005.dat` | the records are ACS supplementary data; **the bases are not licensable** | The two 18 and 20 nt annealing regions may be quoted. The records are not shipped |
| Addgene plasmid 247434 (DMX0001) | Addgene's UBMTA governs the **plasmid**, not its bases | The sequence is used freely. See `AGENTS.md`, "A DNA sequence carries no licence" |
| The ONT native barcoding kit's 24 nt barcodes | taken by LevSeq from the kit; ONT's terms were not read here | **No** — they stay user-held data, as #260 decided for Route A's kit |
| NEBNext oligo documentation | © NEB, all rights reserved | **No** — cite, and re-enter single facts by hand |

The line is the one the repository already draws, and §2's result depends on it: **a sequence may
be used, a document may not be mirrored.** Two short annealing regions printed in running text
are sequence.

## 2. Where the primers bind

### 2.1 The reframing holds: Route B reads the DMX vector

`docs/synthesis-and-assembly.md`, "Cargo validation (optional)", opens: *"The cargo sits in the
DMX vector either way; this section only reads each well."* Three other places in the same
document agree:

- "**Runs in the DMX vector**; the working vector takes only the finished cargo" (iGGA cargo).
- "Assemble into the DMX vector" is step 9 of cargo construction; the working vector appears only
  at the final transfer.
- `docs/research/ap1-demo-project.md` §2.1 puts all 72 parts in the DMX vector at the seating
  step.

So issue #312's framing — "our working vector is not pET-22b(+)" — names the wrong vector.
pLVX-TetOne-Puro-GFP never holds an unvalidated cargo. **The vector Route B amplifies is the DMX
vector, rebuilt from DMX0001 (Addgene 247434).** That is the whole of the collapse.

### 2.2 DMX0001 is a pET-style expression vector

Measured from the Addgene record, 5,839 bp, circular:

| Feature | Span | Strand |
| --- | --- | --- |
| T7 promoter | `[263, 282)` | + |
| lac operator | `[282, 307)` | — |
| RBS | `[321, 344)` | — |
| ccdB | `[477, 783)` | + |
| 6xHis | `[802, 820)` | + |
| T7 terminator | `[892, 940)` | + |
| f1 ori | `[976, 1432)` | + |
| AmpR | `[1524, 2385)` | − |
| ori | `[2611, 3200)` | + |
| rop | `[3629, 3821)` | − |
| lacI | `[4629, 5712)` | − |

That is pET architecture — T7 promoter, *lac* operator, RBS, cargo site, 6xHis, T7 terminator,
f1 ori, rop, lacI — and it is the same architecture LevSeq's primers were written against.
LevSeq's host is pET-22b(+), also AmpR; DMX0001 is AmpR as supplied, which the method swaps for
KanR.

### 2.3 The cassette geometry, and what stays constant

Site positions, measured:

| Enzyme | Recognition site | Cut | Overhang |
| --- | --- | --- | --- |
| BsaI | `[355, 361)` + | top 362 | `AGGA` at `[362, 366)` |
| BsmBI | `[367, 373)` − | top 362 | `AGGA`, the same cut |
| BsmBI | `[785, 791)` + | top 792 | `TTCC` at `[792, 796)` |
| BsaI | `[797, 803)` − | top 792 | `TTCC`, the same cut |

So the cargo occupies `[366, 792)` and the ccdB stuffer that sits there in the stock is **426 bp**
(derived: 792 − 366). Everything outside that span is the same in every well, which is the
property Route B needs. The T7 promoter, *lac* operator and RBS all sit **upstream** of 362 and
survive cargo insertion; the 6xHis and T7 terminator sit **downstream** of 796 and survive too.
Only ccdB is replaced.

### 2.4 LevSeq's own primers bind DMX0001, verbatim

The LevSeq SI prints the two annealing regions in running text, not only in Tables S3 and S4:

> F: 5' - XXXXXXXXXXXXXXXXXXXXXXXX**ATCTCGATCCCGCGAAATTAATACGACTCAC** - 3'
>
> R: 5' - XXXXXXXXXXXXXXXXXXXXXXXX**GCCCCAAGGGGTTATGCTAGTTATTGCTC** - 3'

The 24 `X` are the barcode; the bold region is "a backbone-specific primer that binds to the
cloning vector near the T7 and T7' promotor sites" (SI, "Primer Design").

Run against DMX0001 with `find_binding_sites`, both match **exactly and once**:

| Annealing region | Length | Binds DMX0001 at | Strand | Other sites | Tm | GC |
| --- | --- | --- | --- | --- | --- | --- |
| `ATCTCGATCCCGCGAAATTAATACGACTCAC` | 31 nt | `[245, 276)` | forward | none | 71.0 °C | 45.2% |
| `GCCCCAAGGGGTTATGCTAGTTATTGCTC` | 29 nt | `[883, 912)` | reverse | none | 72.8 °C | 51.7% |

Tm is `liulab_mbio.primers.polymerase.melting_temperature` at the package's Taq conditions. Both
sit above LevSeq's touchdown ceiling of 68 °C, which is what a touchdown from 68 °C down to
63.5 °C is for — so **#298's cycling applies unchanged**.

**This is the answer to #312's first question.** The sites are 121 bp upstream and 120 bp
downstream of the cargo (derived: 366 − 245 = 121; 912 − 792 = 120), they are vector sequence, and
they are therefore identical in every well by construction.

### 2.5 Two other pairs also bind, and are worth knowing about

| Pair | Forward site | Reverse site | Amplicon on empty DMX0001 | Constant flank |
| --- | --- | --- | --- | --- |
| **LevSeq** (SI "Primer Design") | `[245, 276)` | `[883, 912)` | 667 bp | 241 bp |
| **Lund** cPCR, from `si_005.dat` | `[226, 244)` | `[953, 973)` | 747 bp | 321 bp |
| Canonical T7 / T7-terminator sequencing primers | `[263, 283)` | `[878, 897)` | 634 bp | 208 bp |

Lund's own annealing regions are `GCGTAGAGGATCGAGATC` (18 nt, Tm 61.3 °C) and
`ATTCGCCAATCCGGATATAG` (20 nt, Tm 61.3 °C), annotated in their deposited pET acceptor as
"cPCR FWD annealing region" and "cPCR REV annealing region". **Both bind DMX0001 exactly and
once as well** — Lund's acceptor is a domesticated pET28a, so the shared backbone is why. Their
lower Tm matches their published 53 °C anneal, so the two published designs are internally
consistent and not interchangeable in their cycling.

The canonical pair is the shortest amplicon of the three, but it has no published precedent for
this use and its forward primer *is* the T7 promoter, so it would prime any T7 construct in the
lab. **LevSeq's pair is the recommendation**: it is published for exactly this purpose, it binds
our vector unchanged, and #298 already holds its reaction and cycling numbers.

### 2.6 Amplicon size, and how it differs from LevSeq's

**Derived.** With LevSeq's pair and LevSeq's 24 nt marks, the band is **cargo + 289 bp**
(241 bp of constant flank + 2 × 24 nt). With 8 nt or 10 nt marks it is cargo + 257 bp or
cargo + 261 bp.

LevSeq's own expected band is "the size of the gene plus about 100 bp" on pET-22b(+). **Ours is
about 190 bp larger**, because the DMX cassette's entry overhangs sit further from the T7 sites
than pET-22b(+)'s cloning site does. Nothing is wrong; the gel expectation in a Route B protocol
must be computed from the vector, not copied from LevSeq.

### 2.7 What the DMX rebuild must not disturb

The method rebuilds DMX0001 before first use: AmpR → KanR, the four BbsI sites removed, one PmeI
site added outboard of each cassette BsaI site. Measured, none of the first two touches the
amplicon — AmpR is at `[1524, 2385)` and the BbsI sites at 4130, 4490, 4878 and 5217, all outside
`[245, 912)`. **The PmeI insertions do land inside it.**

The left BsaI starts at 355 and the RBS ends at 344, so there are only **11 bp** of room between
them for an 8 bp `GTTTAAAC`; `docs/research/ap1-demo-project.md` §2.1 already records that the
right-hand insertion "lands in or abuts the `6xHis` coding sequence", which is `[802, 820)`.

**The constraint this note adds:** each PmeI insertion must fall strictly between its BsaI site
and the nearer annealing region — inside `(276, 355)` on the left and inside `(803, 883)` on the
right. An insertion outside those windows breaks a primer site or changes which strand it is on.
Both insertions lengthen the constant flank, by 8 bp each if nothing else moves, so the band in
§2.6 becomes cargo + 305 bp for the as-built vector. That is a check the package can run and
should: after rebuilding, confirm each annealing region is still present exactly once.

A second check belongs beside it: **no cargo may contain either annealing region.** The cargo is
codon-optimised coding sequence plus an 11 bp barcode and a 34 bp stuffer, so a 29–31 nt match is
most unlikely, but "most unlikely" is not a guarantee and the check is one call to
`find_binding_sites`.

## 3. What marks a well

Three leads, in the order the lab head gave them.

### 3.1 Lund et al.: NEB index sequences, read on Nanopore

This is the finding that settles the question, and it is in a paper the repository already holds.
Lund et al. 2024 §4.3, "Isolate Sequencing Confirmation", quoted in full where it matters:

> Ninety six (96) unique barcodes were each appended to primers annealing upstream of the
> constructed operon. Sixteen (16) unique barcodes were each appended to primers annealing
> downstream of the constructed operon. **Barcodes were selected from barcodes used in NEBNext
> Multiplex Oligos for Illumina (Catalog No. E73325, NEB).** The resulting set was 96 forward ×
> 16 reverse barcodes, enabling the multiplexed sequencing of 1536 isolates in a single pool.

And the read-out, from the same paragraph:

> An amount of 1.2 µg of DNA was then used as input for a Ligation Sequencing Kit (SQK-LSK109,
> Oxford Nanopore …) with the NEBNext Companion Module for Oxford Nanopore Technologies Ligation
> Sequencing (Catalog No. E7180, NEB).

So Lund took **Illumina index sequences** and **read them on Nanopore**. Demultiplexing was done
in software — `minibar`, "based on the specified combination of sequences of barcodes and primers
at the 5' and 3' ends for each isolate" — with consensus by `Amplicon_sorter`.

Their per-well reaction, for comparison with #298's LevSeq numbers:

| | Lund 2024 §4.3 | LevSeq (#298) |
| --- | --- | --- |
| Template | 0.8 µL of a colony picked into 10 µL of 25% glycerol / 50% LB | 1 µL overnight culture |
| Reaction | 8 µL | 10 µL |
| Polymerase | LongAmp Hot Start Taq 2X Master Mix, NEB M0533 | Taq, NEB M0267 |
| Primers | cherry-picked from a 384PP plate at 100 µM, to 0.5 µM final | 2 µL of a 1 µM mix, to 0.2 µM |
| Cycling | 94 °C 1 min; 35 × (94 °C 15 s, 53 °C 15 s, 65 °C 3 min); 65 °C 10 min | 95 °C 5 min; touchdown 68 → 63.5 °C; 35 cycles; 68 °C 5 min |
| Clean-up | pooled equal volumes, 1.8× SPRIselect, NanoDrop | pooled 5 µL a well, gel-extracted per plate |
| Library | SQK-LSK109 + NEBNext E7180 | SQK-LSK114 |

**One error in the source, recorded rather than repaired.** Catalogue number "E73325" is not an
NEB catalogue number. The likely intent is E7335, NEBNext Multiplex Oligos for Illumina (Index
Primers Set 1), but Lund never lists the 112 barcode sequences they used — the SI's only barcode
table, Table S1, holds the 20 nt *sequence-retrieval* barcodes drawn from Subramanian 2018, which
are a different set for a different job. **Lund's sequencing barcodes are not published.** They
cannot be copied; only the approach can.

### 3.2 NEB's index kits, and what NEB publishes

Read from NEB's own instruction manuals, which served over
`neb.com/-/media/nebus/files/manuals/`. Every `neb.com` HTML product and FAQ page returned HTTP
403, so the manuals are the source throughout — which is the stronger one anyway.

| Family | Catalogue | What is supplied | Indices |
| --- | --- | --- | --- |
| Multiplex Oligos for Illumina, **Index Primers Sets 1–4** (manual v8.0_6/24) | E7335, E7500, E7710, E7730 | adaptor, USER enzyme, universal primer, 12 single-index primers per set, in tubes | 12 per set, 48 in all; **6 nt**, single index |
| Multiplex Oligos for Illumina, **Dual Index Primers Set 1** (manual v6.0_6/24) | E7600 | 8 i5 primers (i501–i508) and 12 i7 primers (i701–i712), in tubes | **20 distinct**, 8 nt, giving 96 combinations |
| Multiplex Oligos for Illumina, **96 Unique Dual Index Primer Pairs Sets 1–5** (manual v10.1_1/26) | E6440, E6442, E6444, E6446, E6448 | a single-use 96-well plate, each well a pre-mixed i5 + i7 pair | **96 unique pairs per set**, 8 nt; 480 across five sets |
| **Unique Dual Index UMI Adaptors**, DNA Sets 1–4 (manual v4.0_6/26) | E7395, E7874, E7876, E7878 | a plate of pre-annealed adaptors, added by **ligation, not PCR** | 96 UDI per set, 8 nt, each with a 12 nt UMI |

**What NEB publishes, and what it withholds.** This divides exactly along the format line.

- **Tube kits: the whole molecule.** The E7335/E7500/E7710/E7730 and E7600 manuals print each
  primer in full, phosphorothioate bonds and all, with its expected index read. NEB's i501, for
  example, is printed as `AATGATACGGCGACCACCGAGATCTACAC` + `TATAGCCT` +
  `ACACTCTTTCCCTACACGACGCTCTTCCGATC*T`.
- **Plate kits: the index only.** The E6440 and E7395 manuals tabulate, per well, the i7 index ID
  and its 8 nt sequence and the i5 index ID with its 8 nt sequence in both orientations — and
  nothing else. The E7395 manual's own revision history, entry 3.0 (9/25), says why, in NEB's
  words: *"Deleted adaptor sequence in NEBNext Adaptor for Illumina Overview because it was
  inaccurate. **The correct adaptor sequence is proprietary.**"*

**That withholding costs us nothing**, and this is the point. Route B needs the *index*, not the
adaptor. The 8 nt index sequences are published in the same manual that calls the adaptor
proprietary.

**Illumina publishes the equivalent, free and without a login.**
`support-docs.illumina.com/SHARE/AdapterSequences/` is the live form of document 1000000002694.
Its "IDT for Illumina UD Indexes" page tabulates UDP0001 onward as **10 nt** i7 and i5 pairs, in
both forward and reverse-complement orientation — UDP0001 is i7 `CGCTCAGTTC`, i5 `TCGTGGAGCG`.
Sixty pairs were read firsthand on 2026-10-06; the table may run longer.

### 3.3 Ordering a dual-index set as lab stock

**Structurally, yes, and it is ordinary catalogue work.** A Route B primer is an index of 8 to
10 nt plus an annealing region of 29 to 31 nt: **37 to 41 nt**, well inside a standard synthesis
scale, and 192 of them is two plates.

**But no vendor page was retrieved firsthand, and no price.** `idtdna.com` returned an unbroken
HTTP 302 redirect loop on every path tried, with and without a browser user-agent, and
`www2.idtdna.com` failed to connect; `twistbioscience.com` served navigation only without
JavaScript. This looks like bot protection rather than a login wall, but that cannot be told from
outside. What came only through a search tool's summary of IDT's pages — a 96-well V-bottom plate
format, a minimum of 48 oligos per plate, and a 25 nmol scale covering 15 to 60 bases at a
guaranteed 10 nmol yield — is **a lead to verify, not a citable fact**, and it is not used in any
arithmetic below.

The route to a number is a written quote, which is a primary source and dates itself. See
[Open gaps](#8-open-gaps).

## 4. Platform coherence

**There is no platform mismatch, and the reason is worth stating plainly.**

An Illumina dual-index primer is `P5 — i5 — Read1 site` and `P7 — i7 — Read2 site`. NEB names the
structure itself in the E7600 manual — *"i7 primers contain indices that are adjacent to the P7
sequence; i5 primers contain indices that are adjacent to the P5 sequence"* — and its published
sequences bear it out: every i5 begins `AATGATACGGCGACCACCGAGATCTACAC` and every i7 begins
`CAAGCAGAAGACGGCATACGAGAT`.

What makes such a primer Illumina-specific is that P5/P7 tail, which grafts the molecule to the
flow cell for bridge amplification, and the Read1/Read2 sites a sequencing primer anneals to.
**It is not the 8 or 10 nt index in the middle.** An index is a short stretch of arbitrary DNA
occupying a defined slot. It carries no platform.

An ONT library is made the other way round. From ONT's own chemistry technical document: the DNA
is *"repaired and the ends dA-tailed in preparation for the dT-tailed native barcode ligation"*,
then *"the samples are pooled together and sequencing adapters are ligated to the barcode ends"*.
The native barcode sits inside fixed flanks — `AAGGTTAA`–barcode–`CAGCACCT` forward — and what
makes a molecule readable is the ONT adapter's motor protein, not any sequence a primer could
carry. No P5, P7, Read1 or Read2 appears anywhere in ONT's barcoding chemistry. **An NEBNext
index primer put on an ONT flow cell as a library would be inert**: nothing to bind, no motor.

NEB sells no barcodes or adapters for Nanopore at all. Its two ONT products, the Companion
Modules E7180 and E7672 — E7180 being the one Lund used — are **enzyme modules only**: repair
mix, end-prep mix, Quick T4 ligase. The barcodes always come from somewhere else.

Route B's primer layout, from LevSeq, is `mark — annealing region`. There is no P5, no P7 and no
Read site, because an ONT library is made by ligating ONT adapters onto a blunt-ended amplicon
afterwards. Dropping an Illumina index sequence into the `mark` slot therefore gives an
ONT-compatible primer that happens to use a well-separated index set someone else has already
validated. **That is exactly what Lund did**, and they read it on a MinION.

So the lab head's suggestion does **not** imply an Illumina read-out. Both readings in #312's
question 3 were offered; the second is the right one, and it has precedent.

Two things do have to be checked, and neither is a platform question:

1. **Separation under the right error model.** Illumina index sets are separated for substitution
   errors, and are additionally constrained by colour balance — NEB's manual explains that a
   four-channel instrument reads A and C on one laser and G and T on the other, and that
   *"if this color balance is not maintained, sequencing the index read could fail"*. None of that
   constrains a Nanopore use, but it means the set's distance structure was chosen for a different
   problem. Nanopore's dominant error is the indel, which is why LevSeq's 24 nt barcodes come from
   the ONT native barcoding kit and why Route A's DMX barcodes are chosen on Sequence-Levenshtein
   distance — see `docs/research/barcode-design.md`. An 8 nt set at Hamming distance 3 is not the
   same guarantee as a 24 nt set at Levenshtein distance 3. `liulab_mbio.barcodes` already measures
   this, so the question is answerable rather than open: **any candidate set is run through the
   package's own distance check before it is ordered.** Length is the thing to watch — 8 nt is a
   third of what LevSeq uses.
2. **No Type IIS site, in either orientation, across the junction.** The mark sits 5' of the
   annealing region inside an amplicon that is never cut, so this matters less here than for a
   cargo barcode — but the method's enzymes are many and the check is free.
3. **An index ID is not a sequence.** `UDP0001`, `S762` and `i501` are vendor nomenclature.
   Carrying an ID across is a traceability choice, and it invites the misreading that the kit was
   used. The repository should carry the sequence and say where it came from.

## 5. How many are needed, and what the stock costs

`src/liulab_synbio/dmx.py` builds Route B as `INDEX_WELLS = 96` forward marks addressing the well
and 96 reverse marks addressing the plate: 9,216 wells on 192 primers.

| Scheme | Forward | Reverse | Wells | Primers to hold |
| --- | --- | --- | --- | --- |
| `ROUTE_B` as coded | 96 | 96 | 9,216 | 192 |
| LevSeq as published | 96 | 96 | 9,216 | 192, built as 8 plates (LevSeq01–08), so 768 wells before more are made |
| Lund as published | 96 | 16 | 1,536 | 112 |

**Lund's set does not reach `ROUTE_B`'s shape.** 96 × 16 is 1,536, one sixth of 9,216. Lund never
needed more. So the lab head's steer gives the *method* but not the *set*: reaching 96 × 96 means
drawing 192 marks, not 112, and the source they are drawn from has to hold at least 96 distinct
sequences on each side.

**Which published sets hold 96 on each side.** Measured against §3.2's inventory:

| Source of index sequences | Distinct forward | Distinct reverse | Reaches 96 × 96? |
| --- | --- | --- | --- |
| NEB E7335 + E7500 + E7710 + E7730 | 48 single indices in all | — | **No** |
| NEB E7600 | 12 i7 | 8 i5 | **No** — 96 combinations, not 9,216 |
| NEB E6440 family, one set | 96 i7 | 96 i5 | **Yes**, see the caveat |
| Illumina UD indexes, as read | 60 i7 | 60 i5 | **No**, from the 60 pairs read firsthand |
| ONT native barcoding kit, as LevSeq used it | 96 | 96 | **Yes** — this is LevSeq's own answer |

**The caveat on the E6440 row is load-bearing and is derived, not quoted.** A "unique dual index"
set means each well carries an i7 and an i5 used nowhere else in the set, so a 96-well set holds
96 distinct i7 and 96 distinct i5 sequences. That inference was not checked against the printed
tables well by well, and it should be before anything is ordered — it is one pass over a table
already published in the manual.

Note what the kit does and does not give. **As sold, an E6440 plate gives 96 combinations**,
because the pairs come pre-mixed in a well. Reaching 9,216 means re-entering the 192 index
sequences by hand into our own primers and ordering them as a free factorial. That is a different
product from the kit, built from the kit's published sequences — which §1 says is allowed, and
§3.2 says is possible, and is exactly what Lund did at smaller scale.

**What a stock costs: not answered.** No vendor price was retrieved from a primary source
(§3.3). What can be said without one: the order is **192 oligos of 37 to 41 nt**, one synthesis,
and it is reusable across every project because the annealing regions are vector constants rather
than cargo-specific. LevSeq's own practice is the comparison — they built 8 plates, not 96, and
covered 768 wells; a lab that will never run more than a few plates at once need not buy 192
primers to start.

## 6. Validation before a full run

**Yes, and both published sources say so in their own way.**

LevSeq's SI is explicit, and it is the sentence #312 was built on:

> If using a different cloning vector, the user must verify that the primers in the "Primer
> Design" section would still amplify the desired region.

`route-b-index-pcr.md` §2 records that it "recommends testing a new primer set on a few wells
before running it at plate scale".

**What this note changes about that requirement.** The verification LevSeq asks for is "would
still amplify the desired region", and §2.4 has now done it in silico: both regions bind the
target vector exactly once, in the right orientation, spanning the cargo. That is the whole of
the stated requirement, met before any bench work. What it does **not** cover, and what a few-well
test still buys:

- the rebuilt DMX vector, not DMX0001 as supplied, is what the bench will hold (§2.7);
- the cycling, which is LevSeq's and has never been run on this flank length;
- the band size, which is about 190 bp larger than LevSeq's and has to be recognised on a gel
  (§2.6).

LevSeq's own measurement is the argument for not skipping it: on a ten-plate run "performed
without optimizing PCR procedure", three plates returned under 60% of variants sequenced, and the
article attributes that to "improper PCR amplification" (`route-b-index-pcr.md` §4). A few-well
test costs one row; a bad plate costs a flow cell lane.

**Recommendation: the protocol prints it as a step**, conditioned on the primer set being new to
this vector, not on every run.

## 7. Whether the package designs them

**Recommendation: the annealing regions are package-held constants, the marks are user-held data,
and the package checks both against the simulated record rather than designing either.**

Three reasons, in order.

**The design problem does not exist.** `primers/` places a primer where a placement allows and
ranks candidates on every check. There is nothing to rank here: two published annealing regions
already bind the target vector exactly once at a Tm that matches published cycling. Designing a
third pair would be work whose output is strictly worse — unpublished, unvalidated, and no
shorter. The package's own restraint rule applies: a narrower thing will do.

**The marks are not ours to choose.** Route A's barcodes are already user-held, read from a copy
the user holds, because the kit is a physical thing someone bought (`dmx.py`, and #260). A plate
of 192 index primers is the same kind of object. Whether they come from an NEB kit, from an ONT
kit or from a custom IDT plate, the lab holds them and the package must not pretend to have
invented them. Shipping any vendor's index table is also the licence verdict §1 reaches.

**What the package should do instead is check.** Three checks, all of which `primers/` and
`barcodes` already support:

| Check | Call | Catches |
| --- | --- | --- |
| Each annealing region binds the simulated DMX vector exactly once, right strand, spanning the cargo | `primers.placement.find_binding_sites` | a rebuild that broke a site (§2.7) |
| No cargo contains either annealing region | the same call, per design | a well that would amplify twice |
| The mark set is far enough apart under Sequence-Levenshtein, and carries no enzyme site | `liulab_mbio.barcodes` | an Illumina set borrowed without re-checking for indels (§4) |

So `primers/` is used, but for evaluation, not placement. That keeps Route B's shape — "a marking
shape, a plate and a depth floor" — and adds the one thing it was missing, which is a pair of
constants and a check that they still hold.

**What this does not decide.** Which mark set the lab actually buys is a bench decision, and #312
is not closed by this note.

## 8. Open gaps

- **Whether an E6440 set really holds 96 distinct i7 and 96 distinct i5.** §5 derives it from what
  "unique dual index" means; it was not checked well by well against the printed tables. One pass
  over a published table settles it, and nothing should be ordered before it is.
- **The price of 192 plate oligos.** No vendor page was reachable: `idtdna.com` returned a 302
  redirect loop on every path, `www2.idtdna.com` would not connect, and Twist's page served
  navigation only. The plate format, minimum and scale figures in §3.3 reached this note only
  through a search summary and are **not** primary. A written quote is the fix.
- **Whether the 24 nt ONT native barcodes may be re-ordered as oligos.** Carried unchanged from
  `route-b-index-pcr.md`. ONT's terms were not read for either note, and this matters because the
  ONT set is the only verified source of 96 marks on each side that was *designed for the error
  model Route B reads under*.
- **The as-built DMX vector has not been measured.** Every coordinate here is DMX0001 as supplied.
  The rebuild's two PmeI insertions land inside the amplicon, and §2.7 states the windows they must
  fall in, but no record of the rebuilt vector exists yet to check against.
- **LevSeq's full primers were not read.** Tables S3 and S4 hold the 192 assembled sequences for
  pET-22b(+). They were not opened, are CC BY-NC, and are not needed — the SI prints the annealing
  regions in running text, which is what §2.4 used.
- **Lund's 112 sequencing barcodes are not published** (§3.1), and their catalogue number
  "E73325" is not an NEB number. Their approach transfers; their set does not.
- **The amplicon has never been run.** Both annealing regions bind in silico at a Tm that matches
  published cycling. Nothing here is a bench measurement, which is the whole of §6's argument for
  a few-well test.
- **Illumina's permission wording was not read firsthand.** The often-quoted line about the
  sequences being provided "for the sole purpose of understanding and publishing the results of
  your sequencing experiments" came only through a search summary. §1 does not rely on it: it
  rests on the repository's own rule that a sequence carries no licence.
- **Route B at 384 wells a plate** stays unpublished, as #298 recorded.

## Sources

- Long, Y.; Mora, A.; Li, F.-Z.; Gürsoy, E.; Johnston, K. E.; Arnold, F. H. LevSeq: Rapid
  Generation of Sequence-Function Data for Directed Evolution and Machine Learning. *ACS Synth.
  Biol.* **2025**, *14* (1), 230–238.
  [doi:10.1021/acssynbio.4c00625](https://doi.org/10.1021/acssynbio.4c00625). Supporting
  Information `sb4c00625_si_001.pdf`, sections "Oligonucleotide Design" and "Ordering Barcode
  Linked Primers". Read 2026-10-06.
- Lund, S.; Potapov, V.; Johnson, S. R.; Buss, J.; Tanner, N. A. Highly Parallelized Construction
  of DNA from Low-Cost Oligonucleotide Mixtures Using Data-Optimized Assembly Design and Golden
  Gate. *ACS Synth. Biol.* **2024**, *13*, 745–751.
  [doi:10.1021/acssynbio.3c00694](https://doi.org/10.1021/acssynbio.3c00694). Section 4.3 and
  supplementary data file `sb3c00694_si_005.dat`. Read 2026-10-06.
- Addgene plasmid 247434 (DMX0001), SnapGene record
  `addgene-plasmid-247434-sequence-494299.dna`, as held under
  `reference_docs/synthesis_and_assembly/dmx/addgene/`. Measured 2026-10-06.
- Qian, Z. et al. Accelerating protein design by scaling experimental characterization. *Nat.
  Commun.* (2026). [doi:10.1038/s41467-026-76740-9](https://doi.org/10.1038/s41467-026-76740-9).
  The DMX vector and its barcode kit.
- New England Biolabs, *NEBNext Multiplex Oligos for Illumina (Index Primers Sets 1–4)*, manual
  v8.0_6/24 (E7335, E7500, E7710, E7730); *Dual Index Primers Set 1*, manual v6.0_6/24 (E7600);
  *96 Unique Dual Index Primer Pairs Sets 1–5*, manual v10.1_1/26 (E6440 family); *Unique Dual
  Index UMI Adaptors DNA Sets 1–4*, manual v4.0_6/26 (E7395 family). Read 2026-10-06 from
  `neb.com/-/media/nebus/files/manuals/`. Every `neb.com` HTML page returned HTTP 403.
- Illumina, *Illumina Adapter Sequences*, document 1000000002694, HTML form at
  `support-docs.illumina.com/SHARE/AdapterSequences/`, pages "IDT for Illumina UD Indexes" and
  "Illumina Unique Dual Indexes". Read 2026-10-06, no login. The versioned PDF 404s at every
  version probed.
- Oxford Nanopore Technologies, *Chemistry Technical Document*,
  `nanoporetech.com/document/chemistry-technical-document`. Read 2026-10-06, for the native
  barcode flanks and the ligation workflow.
- `docs/research/route-b-index-pcr.md` (#298) for every reaction number this note does not
  restate, and for the gap it fills.
- `docs/research/barcode-design.md` for the distance rules a mark set is held to.
- `docs/research/dad-split.md` for Lund's licence verdict, already settled.

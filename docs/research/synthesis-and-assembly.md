---
search:
  exclude: true
---

# Highly parallel DNA synthesis and assembly: the evidence behind the method

Research note for issue #212, under the spec at #210. It holds the evidence behind
`docs/synthesis-and-assembly.md`: every source read and when, what the design took from each,
where it departs and why, what the review measured, and every design decision still open with
its live options.

Four companion documents sit beside it and are cited rather than restated. The first three were
written while the method was designed and are tracked here for the first time; the fourth is
the existing note for the paper the assembly rounds come from.

| Document | What it holds |
| --- | --- |
| `docs/research/synthesis-and-assembly-materials.md` | Every plasmid, primer, enzyme and strain, grouped by the process that uses it |
| `docs/research/synthesis-and-assembly-barcode-kit.md` | The 96-plasmid plate-barcode kit: how it chains, what each group primes, what to order |
| `docs/research/synthesis-and-assembly-departures.md` | Each departure from Takacsi-Nagy et al., worked through base by base |
| `docs/research/protein-library-assembly.md` | Takacsi-Nagy et al. itself: its method, its numbers, its licence |

The last of those is the account of the iterative-assembly paper. Nothing here repeats it;
section 3 says only what this design takes from it and section 4 where it departs.

## How to read this note

**Dates.** A source's read date is the date the local copy was retrieved, taken from that
copy's timestamp. Where a directory note records a date of its own, that date wins. Timestamps
on three files predate the project because the publisher's date survived the download; those
rows carry the directory's retrieval window instead.

**Measurements.** Section 5 carries every number the review computed, with what it was computed
over. A number with no source is not in this note and does not belong in the method page.

**The demo project.** `docs/research/ap1-demo-project.md` specifies one real project against the
method page and answers section 6's open decisions for that project alone. Nothing it decides is
a default here.

**Paths.** Downloaded files are named in section 7 so the provenance survives the file. They
are git-ignored and a fresh clone does not have them, so nothing outside section 7 cites one:
claims cite the source itself.

## 1. What the method is, in one paragraph

Cargo fragments are synthesised in one oligo pool, amplified out of it gene by gene with a
nested three-primer scheme, assembled by Golden Gate into a counter-selection vector, split one
clone per well, and read by nanopore. Validated cargo then goes two ways: through a working
vector that gives it an application, or through rounds of iterative Golden Gate that chain
fragments with a barcode at every junction. Four published systems meet here and none of them
covers the whole of it.

## 2. Sources, and the date each was read

| Source | What it is | Read | Drawn from |
| --- | --- | --- | --- |
| Qian et al. 2026, *Nat. Commun.* | The counter-selection vector, the plate-barcode kit, the barcoding and sequencing protocol | 2026-10-01 to 2026-10-05 | yes |
| Takacsi-Nagy et al. 2026, *Cell* | Iterative Golden Gate in rounds, with a barcode per part | 2026-09-14 | yes |
| Lund et al. 2024, *ACS Synth. Biol.* | Genes built from a cheap oligo pool by Golden Gate, one gene per reaction | 2026-10-01 to 2026-10-05 | yes |
| Subramanian et al. 2018, *Synth. Biol.* | The mutually orthogonal 20-mer primer set | 2026-10-01 | yes |
| Twist product and protocol material | Oligo pool and gene fragment products, lengths, prices, design rules | 2026-09-14 to 2026-09-17 | yes |
| Baker lab three-primer scheme | The nested-PCR variant, as a diagram sent to us | 2026-10-01 | yes |
| Correspondence with the Baker lab | Which primer set, which read-out, what scale | 2026-10-05 | yes |
| Long et al. 2025, *ACS Synth. Biol.* (LevSeq) | Index PCR per well, pooled per plate, read by nanopore | 2026-10-05 | yes |
| Molecular Devices QPix material | Colony picker models | 2026-10-03 | background only |
| Freschlin et al. 2026, *Sci. Adv.* (OMEGA) | Oligo-pool gene libraries, one reaction per subpool | 2026-10-05 | **no** |
| Romanowicz et al. 2026, *bioRxiv* (DropSynth-Gold) | Pooled synthesis in emulsion droplets, one gene per droplet | 2026-10-05 | **no** |

Full citations are in [Sources](#sources).

### Read but not drawn from

**Freschlin et al. 2026** and **Romanowicz et al. 2026** were read while comparing assembly
strategies and neither shaped the design. They are listed here rather than in the method page's
reference table, which means "what the design took".

Both measure what the method page's goal asserts, so both are worth citing the moment that goal
is restated. Freschlin measures oligo-pool gene assembly to 2,598 bp, gives $1.50 to $14 per
gene for constructs up to 2.6 kb, and gives about 246 bp of usable coding sequence per 300 nt
oligo — the denominator the page's per-kb figure is missing. Romanowicz measures architectures
from five 300-mers to twelve 350-mers, about 3 kb, with per-design success rates.

Two of their findings bear on open decisions in section 6 and are cited there: Freschlin pads
every oligo to a uniform 300 nt and screens the padding for enzyme sites, and Freschlin lost
whole replicates to enzyme sites inside a barcode payload.

**Molecular Devices QPix material** supports one equipment row and nothing else. The two
candidate models are real; which to buy is undecided and is not a design question.

## 3. What the design takes from each source, and what is ours

### Qian et al. — inherited

- The counter-selection vector family and its nested cassette: BsaI and BsmBI cutting the same
  four bases, `AGGA` on one side of the counter-selection cassette and `TTCC` on the other.
- The `AGGA`/`TTCC` entry overhang standard itself.
- The 96-plasmid plate-barcode kit, used as supplied: four groups of 24, one barcode per group
  per well, their chaining overhangs, and the three primer pairs that read them.
- Barcoding in heat-lysed lysate on a throwaway aliquot, so a barcoded molecule never grows.
- The acoustic dispense into a 1536-well plate, and the compression of four 384-well plates
  into one.
- The strain rule for the counter-selection cassette, inherited as written. See section 8.

### Qian et al. — ours

- Rebuilding the vector before first use: the marker swap to KanR, removing its four BbsI
  sites, adding one PmeI site outboard of each BsaI site.
- The accession formula for the whole kit, inferred from the four plasmids we hold.
- The mapping of barcode group to primer, derived from the published sequences.
- Using the kit's BsaI for three jobs rather than one.

### Takacsi-Nagy et al. — inherited

The paper's own account is `docs/research/protein-library-assembly.md`. What this design takes
from it:

- The round structure: a destination opened internally, a donor released externally, the
  product reopening for the next round.
- The idea of an internal stuffer carrying the enzyme that opens the destination, with a blunt
  cutter inside the piece that leaves.
- The 11 bp barcode length, and a barcode per part accumulating at the junction.
- Endura electrocompetent cells, recovery and growth at 30 °C each round.
- T7 ligase with its branded buffer.

### Takacsi-Nagy et al. — ours

- The enzyme assignment, the stuffer sequence, the entry overhangs, the barcode set, the read
  platforms and the acronym `iGGA`. None of these is the paper's. Section 4 and the departures
  note give each one's reason.

### Lund et al. — inherited

- The shape of the pipeline: codon-optimise and domesticate, split the gene into oligo-sized
  pieces with overhangs chosen on ligation fidelity, pad, append primer sites, order one pool,
  amplify per gene, assemble by Golden Gate.
- Fixed terminal overhangs so one acceptor vector serves every design, with the split tool
  choosing only the internal ones.
- Adding the Type IIS site and the terminal overhang to the gene **before** the split.
- The error-free ceiling and the fragment-count curve, as evidence rather than as a rule.

### Lund et al. — ours

- The counter-selection vector as the acceptor, in place of their two acceptors.
- The nested three-primer scheme in place of their one amplification per pool, which is the
  Baker lab's variant rather than Lund's.

### Subramanian et al. — inherited

- The orthogonal 20-mer primer set itself, as the inventory the three-primer scheme draws on.

### Subramanian et al. — ours

- The split of that set across the three roles. Undecided; see section 6.3.

### Baker lab — inherited

- The three-primer architecture: an oligo reads `[P1][fragment][P2][P3]`, the outer pair pulls
  a subpool, the inner pair nested inside it pulls one gene, and the Type IIS cuts sit inboard
  of P1 and P2.
- The scale anchors: hundreds of oligos to a well for a subpool, one to eight for one gene.

### Baker lab — ours

- Everything downstream of the second PCR. The diagram is the only description of the scheme we
  hold and it stops at the Golden Gate.

### Twist — inherited

- Which products exist, their length bands, their pool-size tiers and their prices.
- The uniformity gate, and the product it belongs to. A pool whose shortest member is more than
  15% below the longest is refused by the ordering interface — but that is a **Multiplexed Gene
  Fragments** rule, dsDNA at 301–500 bp, not an oligo-pool one. Verbatim, from
  `twist/protocol/dump/DOC4057_Multiplexed_Gene_Fragments_Design_Guidelines_REV1.txt` lines
  43–45: "Length Variation | ≤ 15%\* | The shortest fragment in your pool must be no more than
  15% shorter than the longest fragment." **No vendor publishes a length-spread rule for an
  ssDNA oligo pool at all.** Both Twist Oligo Pools documents held here give only "20–350
  nucleotides" and state no spread, so nothing a pool is ordered against forces a padding rule.
- The flank specification, and the amplification rule: a high-fidelity hot-start polymerase, and
  a cycle count banded by the pool's length — 6–10 at 20–100 nt, 10–12 at 100–150 nt, 12–14 at
  151–350 nt. `oligo-pool-pcr-cycles.md` section 2 gives the two guides it is read from, and its
  section 3 separates it from the Multiplexed Gene Fragments maximum, which is another product's
  and does not apply to a pool.

### Twist — ours

- Any per-kb figure. Twist prices per pool and states none. See section 5.

### Long et al. (LevSeq) — inherited

- The index-PCR read-out: one barcoded primer pair per well, the reverse barcode naming the
  plate, pooled per plate and read by nanopore.
- The published barcode-linked primer plates, which are ordered as given.
- The read-depth and per-base criteria as one of the two candidate pass rules.

### Long et al. — ours

- Using them against this vector family rather than the paper's. The match was measured; see
  section 5.

## 4. Departures from a source, and the reason for each

### From Takacsi-Nagy et al. 2026

`docs/research/synthesis-and-assembly-departures.md` works each one through base by base. The
summary, one line each:

| # | Departure | Reason |
| --- | --- | --- |
| D1 | The two enzyme roles move up one enzyme | The plate-barcode kit is used as supplied and fixes BsaI at the cargo interface |
| D2 | BbsI in the internal stuffer, where the paper has BsaI | D1 pushes the stuffer onto a third enzyme; PaqCI needs an activator oligo on the one digest where a miss is invisible |
| D3 | One entry overhang pair every round, not one per position | `AGGA`/`TTCC` is the published standard; position-specific overhangs would need a different vector per position |
| D4 | No external stuffers on the cargo | The vector supplies BsaI in the right places already |
| D5 | PmeI sits in the vector, not in the cargo | D4 left no external stuffer to carry it; paid once in the backbone rather than on every block |
| D6 | Cargo must be clear of five enzymes, not zero | Our cargo passes through more hands than theirs |
| D7 | No position-specific terminal molecule | One fragment rule serves every position, because no fragment of ours ends the protein |
| D8 | No T2A, so the barcode block is translated inside the protein | A T2A inside the cargo would end translation before the C-terminal part |
| D10 | Fragments come from an oligo pool, not clonal plasmids | The whole point of the method; costs a sequence-and-reformat step they do not pay |
| D11 | Both markers change | The last transfer moves cargo between two vectors, which must not share a marker |
| D12 | The last transfer is one pot with counter-selection | The cargo already presents the working vector's interface |
| D13 | Bacterial expression host, and a working-cassette layer they have no counterpart for | Different application |
| D14 | A second barcode system, and BsaI doing three jobs | The plate-barcode kit is added on top of the per-fragment barcode |
| D15 | Extra strains | Forced by the counter-selection cassette |

Three further departures the review established and the departures note does not yet carry:

- **The barcode set is drawn to an indel-aware metric, with a homopolymer cap the paper has
  none of.** The paper states a minimum Hamming distance of 3 between domain barcodes and
  imposes no GC band and no cap; measured, its own set spans GC 2 to 9 of 11 and reaches a run
  of 6. Ours is a minimum Sequence-Levenshtein distance of 3 within a part list, a homopolymer
  run of 5, and still no GC band. The reason is the read-out: 92.0% of the residual
  discordances of a HiFi read are indels in homopolymers, one every 477 bp, and no Hamming
  distance sees a deleted base at any value. It costs nothing measured — each of the paper's
  three part lists already stands 3 apart under the indel-aware metric, and the cap rejects one
  of its 72 barcodes. `docs/research/barcode-design.md` holds every measurement, and the
  method page carries the rules.
- **Both library reads move to nanopore.** The paper ran linkage on PacBio HiFi and
  representation on Illumina. Our equipment list has one nanopore sequencer. The reason is
  equipment, not biology, and the move is undocumented in the departures note.
- **The source contradicts itself on which molecule each digest opens.** Its Methods put
  BsaI with SrfI on the donor; its own supplementary figure legend puts BsaI on the
  destination. The design follows the figure, which is the assignment that simulates. The
  disagreement is recorded here because neither the paper nor the method page records it.

Departures D9 and D16 were open decisions when the departures note was written. D9, the barcode
rules, is decided and in the method page. D16, whether an assembled clone is sequenced whole,
is **not** decided and the departures note's claim that it is should not be relied on; the open
form is section 6.9.

### From Qian et al. 2026

- **The vector is rebuilt rather than used as supplied.** The marker swap is forced by D11. The
  BbsI removal is forced by D2. The PmeI additions are forced by D5.
- **The route from a transformation to a clonal well is ours.** Qian goes transformation,
  dilution plates for a colony count, bulk plate, picker. The method page spots a frozen
  polyclonal archive and regrows it, which has no counterpart in any source and no stated
  dilution, titre or plating density. Its handling is ours as well: duplicate plates, drawn one
  way and never re-frozen. This is invention and should be marked as such or replaced. See
  section 6.6.
- **A well's barcode combination is derived, not drawn.** Qian generates the combinations
  randomly in a script and records them. Ours factorises the well's address across the four
  groups, group 4 carrying the plate, so the combination is arithmetic and two plates stay apart
  on one flow cell without a file to keep. See section 6.6.
- **The library primers are replaced.** Qian amplifies a whole pool in one reaction with a
  primer pair carrying the BsmBI sites that give `AGGA`/`TTCC`. The three-primer scheme replaces
  that pair with orthogonal 20-mers, which carry no site. Nothing then puts a Type IIS site on
  the oligo. This is the method's one genuine gap; see section 6.1.

### From Lund et al. 2024

- **One acceptor becomes the counter-selection vector**, so the cargo arrives already able to
  be barcoded and released.
- **One amplification per pool becomes two nested ones**, following the Baker lab rather than
  Lund. It costs more reactions, N plus one per subpool, and buys a primer inventory that stops
  growing with the library.
- **The split tool is named but not fixed.** Lund calls NEB's SplitSet tool; they do not
  implement the design method they cite. The method page names the method, not the tool. See
  section 6.1.

### From Long et al. (LevSeq)

- **The primers are used against a different backbone.** The published set is specific to the
  paper's expression vector. Measured against ours they match at full length, so they are
  ordered unchanged. See section 5.
- **An exact match is required on top of the criterion.** LevSeq calls mutations it does not know
  in advance, so its per-base test is what it has. We know the design, so a well passes only on
  an exact match across the whole designed region, and a silent mismatch inside a barcode
  mislabels a member for the rest of the project. See section 6.6.

### From Twist

- **No departure.** The vendor's limits are constraints, not choices. The method page's price
  and per-kb figures depart from the published tiers, which is an error rather than a design
  departure; see section 5.

## 5. The review's measurements

Five independent reviewers read the method page and its supporting documents on one lens each.
Every number below was computed on **2026-10-06** unless another date is given, and the column
says what each was computed over. Numbers computed earlier and carried in from the drafting
session are dated separately.

### Computed against the plasmid maps

| Measurement | Computed over | Date |
| --- | --- | --- |
| PmeI sites: 0 in each counter-selection vector parent | The two depositor maps, 5,839 bp and 5,842 bp | 2026-10-06 |
| SrfI sites: 1 in each parent, inside the counter-selection coding sequence | The same two maps | 2026-10-06 |
| BsaI 3, BsmBI 2, BbsI 4, PaqCI 0 in each parent; one BsaI inside the AmpR marker | The same two maps | 2026-10-05 |
| Four barcode plasmids: BsaI 2, BsmBI 0, BbsI 0, PaqCI 0, AmpR, no counter-selection cassette | The four maps held | 2026-10-06 |
| Part carrier: BsmBI 0, **BsaI 1**, and its own selection is the counter-selection cassette | The carrier's map | 2026-10-06 |
| The barcode chain runs `TTCC`-BC4-`TCAG`-BC3-`CCTT`-BC2-`GTTC`-BC1-`AGGA` on the map strand | BsaI sites in the four barcode maps, against all 96 rows of the barcode table | 2026-10-06 |
| The published index-PCR primers' two constant backbone tails match both parents at full length, each pair bracketing the counter-selection cassette | The two parents against the published primer tables | 2026-10-06 |
| Inserting PmeI outboard of the right-hand cassette BsaI lands in or abuts the `6xHis` coding sequence | Annotated features of the first parent | 2026-10-06 |

### Computed against the published sequences

| Measurement | Computed over | Date |
| --- | --- | --- |
| The internal stuffer is 34 bp, `AGGA` at 0-3, SrfI at 13-20, BbsI cutting to `TTCC` at 30-33, and 34 + 11 is 0 mod 3 | The 34-mer | 2026-10-05 |
| SrfI's site lies inside the piece BbsI releases, so after BbsI it cuts only a discarded fragment | The same 34-mer | 2026-10-06 |
| After the two digests the only complementary end pairs are the two intended ones; one circular product | The destination and donor ends | 2026-10-05 |
| Six of the paper's 72 barcodes carry a stop at our frame offset; fifteen at offset 1; none at the paper's own offset | All 72 published 11 bp barcodes | 2026-10-05 |
| The paper's 72 domains, 23,103 bp, carry zero BbsI, BsaI, SrfI and PmeI, against 11 BsmBI and 15 PaqCI | The published domain sequences | 2026-10-05 |
| An 11-mer drawn to the stated barcode rules alone carries one of the eight enzyme motifs **0.90%** of the time, 1,798 of 200,000 draws; the union bound gives 0.0089 | 200,000 random 11-mers against eight motifs | 2026-10-06 |
| All 96 barcode rows are 97 bp with the same layout and 96 distinct 25 nt payloads; one overhang pair per group across all 24 members | The barcode table | 2026-10-06 |
| Two of the eight constant regions flanking the payload match no published primer | All six primers and their reverse complements against all 96 rows | 2026-10-06 |
| The paper's internal stuffer carries BsaI and SrfI, not PaqCI, which appears nowhere in it | The paper's text, supplement and constants | 2026-10-06 |
| The orthogonal primer supplement has 165 rows kept of 185; the paper's abstract says 166 | The primer supplement | 2026-10-06 |

### Computed against vendor material

| Measurement | Computed over | Date |
| --- | --- | --- |
| 18,000 oligos at 251-300 nt cost **$10,004**; $12,505 is the 301-350 nt band; 18,000 at 20-120 nt cost $4,056 | The oligo pool price table, captured 2026-09-17 | 2026-10-06 |
| Per-kb cost runs $1.85 (251-300 nt) to $2.82 (301-350 nt), synthesis only | The same table, with Freschlin's 246 usable bp per 300-mer as the denominator | 2026-10-06 |

Three of these settle a claim in the method page rather than inform one, and are the
corrections ticket #214 carries: PmeI is added and not removed, SrfI needs no edit, and the
index-PCR primers need no redesign.

## 6. Open design decisions

Each lists what is on the table and the evidence for each option, so a decision is made rather
than rediscovered. A section opening **Decided in #N** is closed and that ticket is where the
argument was had; the rest are open.

Each heading names the experiment the decision belongs to, because the two are judged
differently and a rule written for one does not transfer. **DAD-GGA-DMX** builds cargo and
reads each well back, so a member is picked and its identity read. **iGGA** mixes those
parts in rounds and yields a pool nobody picks from, judged by representation and linkage
reads over the whole of it. `CONTEXT.md` defines both. A decision belonging to neither
says so.

### 6.1 Where the Type IIS sites sit on the oligo — DAD-GGA-DMX

**Decided in #257: templated on the oligo, inboard of P1 and P2.** Every oligo reads
`[P1][BsmBI][span][BsmBI][pad][P2][P3]`, both sites inward-facing. The two outermost cuts across
a part yield `AGGA` and `TTCC`; every internal cut yields the overhang the split chose. One rule
serves the vector interface and the internal junctions alike.

| Option | Evidence |
| --- | --- |
| On the gene, before the split | What Lund does: the site and the terminal overhang go on the gene, then the split, then padding, then the primer sites. It makes the site part of the designed sequence, which the split tool then has to respect. |
| On the amplification primers | What Qian does: the library primer pair carries the BsmBI sites that give `AGGA`/`TTCC`. It costs nothing in oligo length, but the three-primer scheme's inner primers are drawn from an orthogonal set with no site, so the set would have to be retailored or one role exempted. |

**The second option cannot express an internal junction at all.** Every fragment of one gene
carries the same P1 and P2 — that is what lets one PCR2 pull all of a gene's pieces — so a
primer-borne site gives every fragment of that gene the same overhang. Qian can do it because
Qian's members are not split. The cost of retailoring the primer set was never the deciding
term.

**And the 14 nt it would save buys almost nothing.** Counted over the 72 blocks of
`ap1-demo-project.md` §4.1 at a 350 nt oligo: 152 oligos templated against 150 with the sites on
the primers, the worst part at five fragments either way.

An earlier revision of this section read the Baker lab diagram as the second option. That is
wrong, and `ap1-demo-project.md` §9.1 read the same diagram correctly: the cuts are drawn on the
oligo, between P1 and the fragment and between the fragment and P2, which is the first option.

### 6.2 The padding rule — DAD-GGA-DMX

**Decided in #257: pad every oligo to the project's one oligo length.** Filler sits between the
3' BsmBI recognition site and P2, outboard of the cut, so it never enters the product. It is
screened for the eight enzyme motifs, for the orthogonal primer sites, and at both of its
junctions for a motif the join creates.

| Option | Evidence |
| --- | --- |
| Pad every oligo to one length | Freschlin pads to a uniform 300 nt and screens the padding for enzyme sites, having lost whole replicates to a site inside a payload. |
| Pad to a target, then append constant flanks | Lund pads to about 260 nt before appending the 20 nt primer sites, giving 300 nt oligos. |
| Bin by length into sub-pools | Romanowicz builds fixed-architecture libraries by length instead. The vendor will discount split sub-pools. |

Uniform padding makes a pool-uniformity check unfailable, so none is shipped. That reason is the
package's own: the vendor's 15% spread rule is a Multiplexed Gene Fragments rule and does not
reach an ssDNA oligo pool, which is ordered against no spread rule at all (section 3). The second
option reaches the same uniform oligo by a different route and differs only in where the filler
sits. Binning buys nothing at 153 oligos and costs a second pool and a second PCR1.

**The oligo length itself stays a project input**, so this is a method rule with one project
number in it, not a method constant.

**Rejected here, and recorded because the argument looked good:** moving the filler outboard of
P2, between P2 and P3, so PCR2 drops it with the P3 region and every digest stub falls to a
constant 27 nt instead of up to 170 nt. What competes in a one-pot assembly is the number of cut
ends, which is one stub per oligo either way; only the stub's length changes. No source measures
a long stub behaving worse, and the method lists no clean-up between PCR2 and assembly for a
stub length to matter to.

### 6.3 The split of the orthogonal primer set across its three roles — DAD-GGA-DMX

**Decided in #258: 96 inner, 35 P1, 34 P3.** Many outer pairs, few inner.

| Role | Primers | What it indexes |
| --- | --- | --- |
| P2, inner reverse | 96 | one gene inside a batch |
| P1, shared forward | 35 | the batch, with P3 |
| P3, outer reverse | 34 | the batch, with P1 |

96 inner primers makes one batch exactly one 96-well PCR2 plate, which is the Baker shape and
the reason the number is 96 rather than any other. Capacity is 35 × 34 = 1,190 batches, far
past any library this method reaches. An even three-way split was supported by nothing.

A method constant, not a project choice: the plates are laid out once and a slot always means
the same pair. The set is **not yet ordered**, so the split is an order specification — and P3
is the cheapest role to under-order, since 12 and 12 still give 144 batches.

The set's own size stays 165 or 166: the supplement keeps 165 rows, the abstract says 166. The
split is over the rows held. A 166th, if it surfaces, is a spare and forces no re-split.

### 6.4 Batch size — DAD-GGA-DMX

**Decided in #258: 96 genes a batch, divided equally.** 96 is the method's cap and it is
physical — one inner-primer plate, one PCR2 plate. Above 96 genes a library divides into equal
batches rather than full batches plus a remainder, which is what Freschlin's evenness
measurement argues for.

The pieces-per-PCR1 budget is a **project input**. No source gives a number, and gene length is
what a project rationally chooses it from. The method ships the band Baker's own two anchors
derive — 96 genes at one to eight oligos a gene is 96 to 768 pieces a PCR1 — and a project
chooses inside it.

The cap and the piece rule do not compete: 96 binds for short genes, the budget for long ones.
The band's floor sits under Baker's "hundreds of oligos to a well", which under-loads PCR1
rather than overloading it, and no source names that as a failure.

### 6.5 The barcode set's enzyme-site freedom — iGGA

**Decided in #260: on. The barcode is cargo, so it obeys the cargo's enzyme rule.** One list, not
two: the barcode is read against D6's five enzymes, in both orientations, with its flanking scar.
A shorter list of its own would be a second rule to keep in step with the first.

The measurement decides it. An 11-mer drawn to the distance rules alone carries one of the eight
motifs 0.90% of the time, 1,798 of 200,000 draws — one library member in 112 cut apart in its own
round — and Freschlin lost whole replicates to exactly that. The alternative was to screen
afterwards, which is cheaper to state and needs a screen nobody has described.

**The code already does this; the method page is what changes.** `BarcodeRules` carries
`forbidden`, `scar` and `phase`, `find_sites` searches both strands, and
`library.parts.barcode_rules` passes the scheme's enzymes and the method's scar, naming a
junction hit as one. The page listed three rules and not this one. The known limit stays: a site
needing bases from two different barcodes is not found, because those are neighbours only once
the block is built.

Related and separable, and no longer open: the distance metric, the GC band and the homopolymer
cap. A Hamming floor does not imply a Sequence-Levenshtein floor, so the metric is a decision of
its own, and the page now states the three verdicts of `docs/research/barcode-design.md` —
Sequence-Levenshtein 3 within a part list, no GC band, and a homopolymer run of 5. Neither a
band of GC 3 to 7 of 11 nor a cap of 3 came from a source: both are the convention that note
traces to one uncited sentence, and the published set meets neither. Section 4 carries the
departure from the source that is left once they go.

### 6.6 The pass criterion for a well — DAD-GGA-DMX

**Decided in #260, and it is two rules rather than one.** A read is either good enough to call or
it is not, and a call either matches the design or it does not. Neither route changes the second,
so the judgement is shared and only the floor travels with the route.

**The floor travels with the route, because each number was measured on its own.** Qian's
consensus depth above 150 was measured on a pooled amplicon carrying four UMIs, where 1,536 wells
demultiplex in-read. LevSeq's twenty wanted and ten tolerable were measured on one amplicon per
well, with the index on the primer. Picking one of the two for both routes would apply a figure
to a library prep that never produced it. Each is a default a project may raise. The wanted depth
is exceeded; the tolerable one is reached, so a well at exactly ten reads warns rather than going
unjudged.

**Below the floor a well gets no verdict, not a fail.** That is already `CONTEXT.md`'s rule for a
`Check` — one no sourced threshold judges carries no verdict and says so — and here it reaches
the bench, because reformatting compacts out the wells that failed. Compacting out a well nobody
read throws away a design that may be clean. LevSeq's "low" is this state, and a low well is
re-read or re-picked rather than discarded.

**A pass is an exact match across the whole designed region** — both entry overhangs, the
fragment, the stuffer and the barcode. A well with more than one consensus, or more than one
significant mutation, is mixed and fails. A silent mismatch is not a pass: it leaves the protein
right and the barcode wrong, and a barcode that no longer names its member cannot be recovered by
the linkage read, which is the only thing that reads the finished pool. This is Qian's exactness
over LevSeq's statistics, and section 4 records it as a departure from LevSeq.

**A barcode combination reaches a well by arithmetic, not by a recorded draw.** The well's
address factorises across the barcode axes, and one axis is the plate. Route B's is published:
LevSeq uses 96 forward barcodes for the well and 96 reverse for the plate, 9,216 wells from 192
primers. Route A's follows the same shape — three DMX groups address the well, 13,824
combinations against the 1,536 needed, and group 4 is the plate, one barcode across a whole plate
from a reservoir. That is what tells two plates apart on one flow cell, by construction rather
than by a file, and it lets a demultiplexer check an address instead of trusting one. Qian
generates the combinations randomly in a script and records them; section 4 carries the
departure.

**Colonies per design: four by default, and a project input.** Lund's four gave 343 of 458 genes
error-free and it is the only measured anchor, but it is not fixed anywhere. The same table runs
100% clean at 2 fragments, 84.6% at 5, 40% at 12 and 0% at 16, so four is right for a short gene
and thin for a long one. The protocol states each design's predicted chance of a clean colony
from that curve and that design's own fragment count, so a project raising the number reads the
reason rather than guessing it. A design with no passing well is re-picked from the same archive
spot before it is re-synthesised.

Underneath all of it is the route to clonal wells itself, which section 4 marks as invention.
That stands: this decision judges wells, it does not source the spot-and-regrow step that fills
them.

### 6.7 How donor pools are built — iGGA

**Decided in #259: one pool per position, pooled equimolar from individually prepped plasmids,
once.** This is the paper's own structure, read from METHOD DETAILS p. e4: "For each TF family,
individual plasmids were combined into equimolar pools of 'N', 'DBD' (bZIP, Forkhead, or ETS)
and 'C' domain sequences or of full-length gene sequences."

The up-front-or-per-round fork dissolves. A pool belongs to a position, not to a round: N is the
round-1 destination, DBD donates once, C donates once. Both options name the same three pools.

Pooling is by DNA mass from one prep per part, not by culture volume. Every part sits in the
same backbone and differs only by insert, so pooling equal volumes of culture is equimolar in
cells rather than in plasmid — and a prep per part leaves a reusable part collection rather than
a cost.

**The constraint either option inherited stands:** our donor carries no release sites of its
own, so unlike the paper's it cannot be a PCR product or a gene fragment. It must be cloned
first.

**Whether a pool's members were read first is not an iGGA question.** An iGGA round consumes a
pool of plasmids and asks nothing about their provenance. How a part reaches that pool, and
whether it is validated on the way, belongs to cargo synthesis and the DMX read-out — 6.6 and
6.12.

### 6.8 Clean-up ratios — iGGA

**Decided in #259: the paper's two ratios — 2X after each digest, 1X after the ligation.** Read
from METHOD DETAILS pp. e4-e5: "Digestion products were SPRI-purified at a 2X volume ratio and
eluted in H2O", then "Ligation products were SPRI-purified at a 1X volume ratio and eluted in
H2O."

The two do different jobs, and the departures note rests the whole PmeI argument on the 2X step
removing the 17 bp and 13 bp stubs. One ratio throughout is supported by nothing read.

Both are constants in `liulab_synbio.igga.bench`, landed by #280.
`synthesis-and-assembly-materials.md` now carries them beside the beads, which is the half of
this decision that was missing.

### 6.9 The per-round bound — iGGA

**Decided in #259: mass in, a titre plate out, and no inherited coverage multiple.**

**In — the paper's masses, and nothing to choose.** 1 µg of plasmid pool per digest, 20 ng of
digested destination per 200 µL ligation, up to 100 ng of purified ligation product
electroporated. All three are constants in `liulab_synbio.igga.bench`.

**Out — a departure, recorded as one.** The paper plates nothing. Digest, SPRI, ligate, SPRI,
electroporate, recover, grow, prep, next round: no colony count appears anywhere in it, and no
round is gated on one. Ours plates a dilution of each round's recovery, and a no-donor control
carried through the ligation from the same digest, and gates the round on net colonies. Both
plates grow during the 12-16 hour outgrowth, so the cost is two plates and no extra day.

Why depart. The design argues parental background away — BbsI cuts twice, SrfI cuts whatever
BbsI missed, T7 refuses the resulting blunt ends, and the 2X SPRI takes the stubs — and the
control plate is what measures that argument rather than trusting it. Residual parental is not
noise: it is the previous round's library, one position short, and it rides into the linkage
read as truncated members. A round that collapsed cannot be repaired by a later one.

**The 300x multiple is dropped, and a completeness replaces it — decided in #297.** 300x is
Qian's, measured at a single transformation whose purpose was pickable clones for an arrayed
collection. iGGA yields a pool that goes to a pooled screen; no member is ever picked or
re-identified, so the gate does not transfer. A multiple does not transfer between rounds
either: one ratio against each round's own products is a different risk every round. What holds
across rounds is the completeness — the chance no member of that round is missing — so that is
what a project states, and `liulab_synbio.igga.coverage` computes each round's floor from it by
Clarke & Carbon's rule. It refuses a default, and reports the multiple the floor works out at
beside the chance a named product is missing. It is a floor rather than a sufficiency claim: it
assumes every member equally represented, and synthesis skew breaks that.

**What judges the library is the read, not the plate.** The paper's own quality claim is barcode
sequencing of the plasmid library for representation, and a long-read amplicon spanning the gene
and its barcodes for linkage: about 95% of reads carried three valid barcodes, nearly 90% of the
library was correctly linked, and over 90% of it held above 80% fidelity. Those are the method
page's steps 7 and 8. The titre plate gates a round; the reads judge the library.

**How many electroporations a round takes is not predicted.** One takes up to 100 ng, and the
paper names no count. No source read gives a ligation-product efficiency for Endura, and the 30 °C
recovery the method uses has no vendor figure at all — `bench-numbers.md` already records that
gap. So the protocol states the round's colony target and says to pool electroporations until the
measured titre reaches it, and the package predicts no number. A lab that measures its own
efficiency holds it the way it holds a ligase matrix or a price record.

**D16 is answered by the paper.** "The sequences of a sample of individual clones from both
libraries were confirmed by whole plasmid next-generation sequencing (Plasmidsaurus)" — a sample
of a finished library, not a per-round product check and not every member. The departures note
and the method page were describing different objects.

### 6.10 The primers for the library reads — iGGA

**Decided in #258: two pairs, each against sequence that is already constant.** Nothing is
added to the oligo — 350 nt is the top of its price band, so a fourth primer role is a tier
change, not a rounding error.

| Read | Forward | Reverse |
| --- | --- | --- |
| Representation, across the barcode block | the retained 34 bp stuffer — a method constant | the vector past the final `TTCC` — per project |
| Linkage, design to barcodes | the vector, before the first `AGGA` — per project | the vector, past the final `TTCC` — per project |

The stuffer is constant across every member because the last round keeps it: the C position
needs its barcode and a capping block carries none. The vector anchors are per project because
the working vector is the user's own — mbio designs a pair against a record, synbio asks for it
with the method's constraints.

Both pairs are designed by the package against the simulated finished record, and are not
written here. Reuse of `dmx0`/`dmx7` was the other option and is not needed; their sequences
are now identified, and they read a barcoded DMX amplicon, not an assembled library.

### 6.11 What the part carrier brings into the assembly — neither; the working vector

Measured: the part carrier has no BsmBI site and **one BsaI site**, and its own selection is the
counter-selection cassette. A part that fits its default overhang pair enters the
working-cassette reaction as the whole carrier plasmid, so both of those enter with it. Only the
marker separates the two outcomes — the carrier is KanR, the working vector AmpR or CarbR.

The seated plasmid itself is no longer unsimulated. The carrier is supplied opened at one blunt
point: `dmx/addgene/pCR-Blunt II-TOPO.dna` carries a single head-to-head `GCCCTT`/`AAGGGC` run,
and the Zero Blunt TOPO user guide (`bench/thermo/zeroblunttopo_man.pdf`) says topoisomerase I
cleaves after 5'-CCCTT on each strand, which puts the blunt point at **offset 336** of the
3,519 bp circle, between the EcoRI sites at 324 and 342. `liulab_synbio.seating.seat` inserts a
part there and hands back the carrier plasmid: a 52 bp part gives a circular record of 3,571 bp
keeping all 12 features, with two BsmBI sites and a digest that releases the part on its own
overhang pair, backbone whole. The carrier's one BsaI site is at **797**, inside `ccdB`
**585-888**, so an edit removing it would fall in the counter-selection cassette's own coding
sequence. What is still unmeasured is the reaction, not the record.

| Option | Evidence |
| --- | --- |
| Leave it, and let the marker do the work | The method page already adds a default-pair part as the carrier plasmid, with no digest and no PCR of its own, which is the whole saving of holding parts as plasmids. Nothing yet shows what the carrier's BsaI site and counter-selection cassette do once they are in that reaction. |
| Amplify every part out of its carrier | Only the part then enters, and the retailoring primer form already exists for exactly this shape. It costs a reaction per part, which is what entering as the plasmid was written to avoid. |
| Domesticate the carrier once | One edit to a lab resource built once removes the BsaI site for every part afterwards. The counter-selection cassette is the carrier's own selection and cannot be removed the same way. |

### 6.12 When a design needs validation — DAD-GGA-DMX

**Decided in #295: polyclonal stands, and a project that wants a read names a fragment-count
floor.** The three options this section carried were never three choices. Whether a design is
read back is already a project choice — `CONTEXT.md` says so — so what was open was the default
and how the project states it.

**Polyclonal by default, because the common case never reads a member.** A library headed for a
pooled screen takes its identity from that screen's own sequencing, so a per-design read buys a
cleaner input rather than a verdict. The Baker lab skips validation for most of its libraries for
exactly that reason, by correspondence with the group relayed to this repo. That is **practice,
not a rule**, which is why it sets a default and constrains nothing: Lund read all 458 of their
genes, and the AP-1 demo reads all 72 of its designs.

**One fragment-count floor replaces the fork.** A project states the fragment count at or above
which a design is read. Omitted, nothing is read; `0` reads every design. The three options
become its three settings, so none of them has to win.

**The floor is the project's number, and none ships.** Lund's curve gives a design's chance of a
clean colony, not the chance worth paying to check. That second number needs what an error costs
downstream, which is the screen, which the method cannot see — the same missing input that
made #267 refuse a coverage multiple, and #299 a route threshold. A shipped floor would be an
invented constant. The curve is what a project reads to pick the number, and the protocol prints
it per design.

**It is a per-design quantity, which is why the switch is not a boolean.** The chance runs from
certain to hopeless across one project's own designs, so a project-wide on-off would read every
two-fragment design that never needs it, or skip every long one that does.

The fragment-count table is Lund's, held in `long_fragment_GGA/README.md`, and
`liulab_synbio.dmx.CLEAN_COLONY_CURVE` carries it. #306 restores the two anchors it dropped,
and #305 builds the floor.

### 6.13 The pass mark for the library reads — iGGA

**Decided in #346. The representation read takes the pooled-screen bar; linkage fidelity stays
open.** Joung et al. 2017 set it for a plasmid library counted by a barcode amplicon before a
screen — under 0.5% of members undetected, a 90th/10th percentile skew ratio under 10, judged at
over 100 reads a member — and section 3.5 of `docs/research/vector-qc-panel.md` quotes it with
its citation. Each is a default a project may tighten and may not loosen, as a read depth is in
`liulab_synbio.dmx`.

**It transfers because the counted amplicon is length-matched.** The usual objection is that this
library's members span about 0.5 to 2.3 kb where an sgRNA is 20 nt, so amplification bias differs.
It does not apply: the representation read amplifies the barcode block alone, anchored on the
retained internal stuffer, so every member's amplicon is the same length. The size skew lives in
the library, and this read is what measures it.

**The read depth is computed, not stated.** At 100 reads a member, AP-1's 13,824 combinations
take about 1.4 M reads — a fraction of one run. Imkeller's Table 2 then prices the skew in screen
coverage downstream: p90/p10 2.5 wants 200-fold, 5 wants 300-fold, 10 wants 400-fold.

**Linkage fidelity keeps no mark.** Joung's ">70% perfectly matching guides" measures a guide's
synthesis fidelity, not whether a barcode still names its part, so it does not transfer.
Takacsi-Nagy's Figures 1D and 1E stay attribution: what one library reached, not a bar. H28 keeps
that one quantity and loses the other two.

### Also open, and smaller

- **The design order** — DAD-GGA-DMX. For iGGA cargo the barcode and the 34 bp stuffer sit
  inside the synthesised block, so both have to exist before the gene is split; the page
  splits first.
  Which step moves depends on 6.1, since the Type IIS sites carry the same constraint.
- **The read-out route** — DAD-GGA-DMX. **Decided in #260: the toolkit models both**, as two
  marking steps behind one judgement, and in #299: nothing picks between them, so a project
  names its route. The capacity ceiling once attributed to the
  primer set is wrong: forward and reverse barcodes combine freely and no dual index is needed,
  so Route B reaches 9,216 wells on 192 owned primers and Route A 24⁴ on the kit, and what binds
  either is the flow cell and one reaction per well.
- **Which vector parent** to rebuild — both, since both run in the DMX vector. The two
  differ only in which strand carries the cassette and release the same cargo, so the rule
  as written does not discriminate.
- **The colony picker model** — DAD-GGA-DMX. The QPix 420 or the QPix FLEX, a purchase
  decision rather than a design one. The method page names neither, only that the choice is
  open.
- **Whether acronyms are expanded on first use** — neither; the page itself, now published
  in the site navigation rather than kept as a working file.

**Where the added PmeI sites go is closed.** It recorded that there is no free space on one
parent without sacrificing an annotated tag. Measured on the rebuilt DMX0001, the windows a
blunt site may sit in — each bounded by a releasing cut and the primer that reads a well — are
**277..363 and 416..503**. `6xHis` ends at 432, so 433..503 is free of every annotation and a
site at 437 sacrifices nothing. The worry came from inserting immediately outboard of the BsaI
site; the rule only asks for a site between that cut and the primer's footprint.
`docs/research/dmx-destination.md` carries the rebuild and what it measures.

## 7. Provenance of downloaded files

Every file below is under `reference_docs/synthesis_and_assembly/` and is git-ignored, except
the two marked **tracked**, which the destination's rebuild reads and which therefore ship
converted to GenBank. For the rest, this table is what survives the file. Retrieved dates are
the local copies' timestamps; a window is given where a publisher's own date survived the
download.

| File | Source | Retrieved |
| --- | --- | --- |
| `dmx/paper/` article and supplementary PDFs, with their text dumps | Qian et al. 2026, doi:10.1038/s41467-026-76740-5 | 2026-10-01 to 2026-10-03 |
| `dmx/dmx-barcodes.tsv` | Extracted from the same paper's Supplementary Table 3 | 2026-10-05 |
| `dmx/addgene/addgene-plasmid-247434-*.dna`, **tracked** as `tests/data/dmx0001.gb` | Addgene 247434, depositor's map, first vector parent | 2026-10-04 |
| `dmx/addgene/addgene-plasmid-247435-*.dna` | Addgene 247435, depositor's map, second vector parent | 2026-10-04 |
| `dmx/addgene/addgene-plasmid-255161-*.dna` and `.dna.seq` | Addgene 255161, barcode kit group 1 index 1 | 2026-10-02 |
| `dmx/addgene/addgene-plasmid-255185-*.dna` | Addgene 255185, barcode kit group 2 index 1 | 2026-10-05 |
| `dmx/addgene/addgene-plasmid-255209-*.dna` | Addgene 255209, barcode kit group 3 index 1 | 2026-10-05 |
| `dmx/addgene/addgene-plasmid-255256-*.dna` and `.dna.seq` | Addgene 255256, barcode kit group 4 index 24 | 2026-10-02 |
| `dmx/addgene/pCR-Blunt II-TOPO.dna`, **tracked** as `tests/data/pcr-blunt-ii-topo.gb` | The part carrier's published map: Thermo Fisher / Invitrogen catalogue map, Zero Blunt TOPO kit | 2026-10-05 |
| `dmx/jason_email.md` | Correspondence with the Baker lab, quoted with permission; not redistributable | 2026-10-05 |
| `prot-assembly/` article, supplemental figures, Table S1 and its sheet dumps | Takacsi-Nagy et al. 2026, doi:10.1016/j.cell.2026.07.054, CC BY 4.0 | 2026-09-14 |
| `long_fragment_GGA/Lund2024/` article, supporting information and four supplementary workbooks | Lund et al. 2024, doi:10.1021/acssynbio.3c00694 | 2026-10-01 |
| `long_fragment_GGA/Lund2024/sb3c00694_suppl-data/*.dat` | The same paper's supplementary vector records; the files carry the publisher's 2025-03-24 date | 2026-10-01 |
| `long_fragment_GGA/Lund2024/primer set paper*.pdf` and `*.xlsx` | Subramanian et al. 2018, doi:10.1093/synbio/ysx008 | 2026-10-01 |
| `long_fragment_GGA/Baker Lab Strategy.png` | Sent by the Baker lab; the only description of the three-primer scheme we hold | 2026-10-01 |
| `LevSeq/papers/` article and supporting information | Long et al. 2025, doi:10.1021/acssynbio.4c00625 | 2026-10-05 |
| `Freschlin2026/` article, supplement and Data S1 | Freschlin et al. 2026, doi:10.1126/sciadv.ady2279; the workbook carries the publisher's 2026-04-20 date | 2026-10-05 |
| `Romanowicz2026/` preprint and supplementary figures | Romanowicz et al. 2026, doi:10.64898/2026.05.29.728538, version posted 2026-06-01 | 2026-10-05 |
| `twist/oligo_pool/*.png`, `twist/GeneFragments/*.png` | Screenshots of Twist's own product and pricing pages | 2026-09-14 to 2026-09-17 |
| `twist/GeneFragments/MultiplexedGeneFragments/` price sheet, product sheet and order form | Twist product material | 2026-09-14 to 2026-09-17 |
| `twist/protocol/` four guides, with their text dumps | Twist documents DOC-001498, DOC-001499, DOC4057 and DOC_4045 | 2026-09-14 to 2026-09-17 |
| `moleculardevices/` two QPix brochures | Molecular Devices product material | 2026-10-03 |

Every `*.txt` beside a PDF is a `pdftotext -layout` dump of it, made so the text can be
searched; every `dump/` and `form_dump/` is derived the same way. None of these files may be
redistributed except the two CC BY papers and their supplements.

The two tracked files were converted once with `build_dmx_vector.write_genbank`, which writes
each feature as `/label=` and nothing else. They carry the bases, the topology and the features;
they drop SnapGene's enzyme set, its auto-matched primer library and the vendor notes, and each
record was renamed for the LOCUS line, which takes neither a space nor `®`. Redo one by reading
the `.dna` with `liulab_mbio.io.read_record` and writing it back through that function; the
bases, the topology and the one coding sequence the rebuild reads are what must survive, and
`tests/scripts/test_build_dmx_vector.py` holds the rest to it.

One file in that tree is ours rather than downloaded: a short script that rebuilds the sheet
dumps of one supplementary workbook. It stays beside the workbook, because the workbook is what
it reads and the workbook is not redistributable. Every other document the lab wrote for this
method is tracked, as the three documents named at the top of this note.

Twist prices move. The screenshots are dated and any figure taken from them carries that date.

## 8. Two things the review got wrong, and one inherited claim

Recorded so they are not acted on.

- **A missing equipment row is not a defect.** One lens compared the equipment table against a
  five-row version that had been superseded; the instruction afterwards was that an
  electroporator is basic equipment and should not be listed. The page is right.
- **The strain claim is inherited, not a transcription error.** The source states that target
  plasmids must be propagated in a strain resistant to the counter-selection cassette and names
  one. That strain's published genotype carries no resistance allele. The resolution is not a
  resistance question at all: the cassette is transcribed from a T7 promoter under a *lac*
  operator with the repressor on the plasmid, so a strain without T7 polymerase never
  transcribes it. The consequence the method page does not state is that the strain
  **receiving** cargo must carry T7 polymerase, or there is no counter-selection.
- **Two lenses disagreed about which supplementary table holds the index-PCR primers.** They do
  not conflict: the raw barcodes are in the first two tables, the full-length barcode-linked
  primers in the next two, and the plate maps in the eight after that. The method page pointed
  at the plate maps, which is wrong either way, and the source's own main text makes the same
  mistake.

## Open gaps

Nothing here should become a package default.

- **No source documents the three-primer scheme in text.** One diagram and one email are all we
  hold; everything downstream of the second PCR is ours.
- **The split tool is unresolved.** The method it names is not the tool its own source ran, and
  neither is pointed at.
- **No source states a per-kb cost.** The figure in the method page is ours and is synthesis
  only.
- **The accession formula for the barcode kit is inferred from four points** and is unconfirmed
  against the depositor's own kit listing. The sequences do not rest on it.
- **The two universal flanking primers' sequences are ours by derivation**, not the source's.
  The names are the paper's; it publishes no sequence, and the two unmatched constant regions of
  the barcode kit are the only candidates. A wrong call shows as a failed amplification.
- **No vendor or catalogue number is recorded for any reagent**, and one strain is named only
  by a property. Identity matters for the branded ligase buffer, the cloning kit and the
  competent cells.
- **No polymerase is named anywhere**, for either amplification, the amplicon reactions or the
  colony PCR, although every source treats the choice as load-bearing.
- **The protein cost of the assembly junction is never stated** on the page. Three rounds is 45
  bp, fifteen residues, translated inside every library member, and no rule fits a fragment's
  boundary codons to the entry overhang.
- **Nothing lyses a culture and nothing is purified before the amplicon reactions.** Qian
  heat-lyses at 98 °C for 30 minutes as its own step, then minipreps the pooled plate and puts
  20 ng of purified DNA into each of the three reactions. The page goes from growth to an
  acoustic transfer to a crude-lysate reaction. Bench detail, owned by the detail protocol.
- **Route B pools per plate and then consumes one pool.** Nothing combines or normalises the
  per-plate pools, and nothing tells two plates apart on one flow cell. LevSeq carries plate
  identity in the reverse barcode and combines each plate's pool in equimolar amounts.
- **The compacted plate's status is unstated** — whether it is frozen, becomes the archive, or
  supersedes the one cargo synthesis makes.
- **Two length conventions are in use for the same kind of block**, one counting the 5'
  overhang only and one counting through both. They agree in effect; applying the wrong one
  builds a capping block one base short.

## Sources

Read dates are in section 2 and provenance in section 7.

- Qian, Z. et al. (2026) Accelerating protein design by scaling experimental characterization.
  *Nat. Commun.* [doi:10.1038/s41467-026-76740-5](https://doi.org/10.1038/s41467-026-76740-5),
  with its supplementary information, Supplementary Table 3 (the barcode payloads) and
  Supplementary Table 4 (the primers).
- Takacsi-Nagy, O. et al. (2026) Synthetic transcription factors designed by domain
  recombination enhance CAR T cell antitumor function. *Cell* 189, 1-20.
  [doi:10.1016/j.cell.2026.07.054](https://doi.org/10.1016/j.cell.2026.07.054). CC BY 4.0.
  Read in full for `docs/research/protein-library-assembly.md`, which is the account of it.
- Lund, S., Potapov, V., Johnson, S. R., Buss, J. and Tanner, N. A. (2024) Highly parallelized
  construction of DNA from low-cost oligonucleotide mixtures using Data-optimized Assembly
  Design and Golden Gate. *ACS Synth. Biol.* 13, 745-751.
  [doi:10.1021/acssynbio.3c00694](https://doi.org/10.1021/acssynbio.3c00694). All authors were
  employed by the enzyme vendor, which funded the work.
- Subramanian, S. K., Russ, W. P. and Ranganathan, R. (2018) A set of experimentally validated,
  mutually orthogonal primers for combinatorially specifying genetic components. *Synth. Biol.*
  3, ysx008. [doi:10.1093/synbio/ysx008](https://doi.org/10.1093/synbio/ysx008), with
  Supplementary Table 1 (the primers and their keep or drop calls) and Supplementary Table 2
  (the cross-talk matrix).
- Long, Y., Mora, A., Li, F.-Z., Gürsoy, E., Johnston, K. E. and Arnold, F. H. (2025) LevSeq:
  rapid generation of sequence-function data for directed evolution and machine learning.
  *ACS Synth. Biol.* [doi:10.1021/acssynbio.4c00625](https://doi.org/10.1021/acssynbio.4c00625),
  with its supporting information. Preprint:
  [doi:10.1101/2024.09.04.611255](https://doi.org/10.1101/2024.09.04.611255).
- Freschlin, C. R., Yang, K. K. and Romero, P. A. (2026) Scalable and cost-efficient custom gene
  library assembly from oligopools. *Sci. Adv.* 12, eady2279.
  [doi:10.1126/sciadv.ady2279](https://doi.org/10.1126/sciadv.ady2279). Read, not drawn from.
- Romanowicz, K. J., Hinton, S. R., Villegas, N. and Plesa, C. (2026) DropSynth-Gold: Golden
  Gate assembly in emulsions extends multiplexed gene libraries to greater lengths. *bioRxiv*
  [doi:10.64898/2026.05.29.728538](https://doi.org/10.64898/2026.05.29.728538). Read, not drawn
  from.
- Twist Bioscience product pages, price sheets, the multiplexed gene fragment order form, and
  the four protocol guides DOC-001498, DOC-001499, DOC4057 and DOC_4045. Vendor material, not
  redistributable.
- Molecular Devices QPix colony picker brochures. Vendor material, not redistributable.
- The Baker lab's three-primer diagram and the correspondence that came with it. Shared
  directly, not published.
- Addgene depositor maps for the two vector parents, four barcode plasmids and the part
  carrier, by plasmid number.

## The AP-1 demo, planned and compared against the source

Run on **2026-10-06** for issue #220, over the 72 AP-1 proteins of
`docs/research/ap1-demo-project.md` — 24 N, 24 DBD, 24 C — with the inputs in
`docs/examples/ap1-library`. The design was frozen and written to disk before the source's own
parts were opened, so the comparison below is a check and not an input. The source is
Takacsi-Nagy et al. 2026, Table S1, CC BY 4.0.

### What the input was, and what it was not

Amino acid sequences and a position each, and nothing else. The wild-type protein of every
domain was recovered from the full-length natural TF sheet of Table S1 — the sheet that carries
no design choice — and cross-checked against the two wild-type residue columns on the N and C
sheets. They agree on all 48. No overhang, stuffer, barcode or codon was given to the planner.

### What the pipeline chose

| | Planned |
| --- | --- |
| Entry overhangs | `AGGA` into N, `AGAT` into DBD, `GCAT` into C |
| Cloning scar | `TTCC` |
| Internal stuffer | 34 bp at every position |
| Barcodes | 72 distinct 11-mers, minimum Hamming 4 in each pool, longest run 4, GC 18-82% |
| Residues moved | 27, over 26 of 48 junction termini |
| Synonymous codons moved by domestication | 55 |
| Coding bases over the 72 parts | 23,103 |
| Blocks | 133 to 1,149 bp, 30,519 bp in all |
| Product | 2,276 bp, 13,824 distinct constructs, every check passing |
| Colonies for 0.99 completeness | 183, then 6,306, then 195,386 — the floor 6.9 settled on, in place of 300x |

23,103 bp is the source's own published total, to the base. The retained tail reads 75 bases,
25 codons, no stop, opening `RKVFSPGRRQF` — the specification's section 3.2 derived both.

### The comparison, difference by difference

| What | Theirs | Ours | Classification |
| --- | --- | --- | --- |
| Overhang count | four: `CTCC`, `GGAG`, `CCGA`, scar `AGCG` | four: `AGGA`, `AGAT`, `GCAT`, scar `TTCC` | **our defect** — the method asks for one pair at every position (D3) and the planner cannot express it |
| Which overhangs | fixed by their own standard | `AGGA` pinned by the destination, the rest chosen for fewest residues moved | two valid designs |
| Internal stuffer | 56 bp, 58 bp terminal | 34 bp throughout | recorded departure D2 |
| External stuffers | 29 bp each end, on the block | 29 bp each end, on the block | **our defect** — D4 puts them in the carrier, and the planner has no way to say so |
| Barcode length | 11 | 11 | agreement |
| Barcode distance | minimum Hamming 3, 4, 4 per pool | 4, 4, 4 | two valid designs |
| Barcode set | their 72 | 72 others, none shared | two valid designs; both sets are stop-free at our frame |
| Residues at the N junction | 15 of 24 changed | 15 of 24 changed, 6 of them the same part | two valid designs |
| Residues at the middle junction | not recorded on their sheet | 11 of 24 changed | two valid designs |
| Residues at the C junction | 19 of 24 first residues changed, all to R | none — our terminal position is charged nothing | two valid designs |
| Reserved enzymes | their blocks carry 2 BsaI, 2 BbsI, 1 SrfI, 2 PmeI by design | the same, and no BsmBI | **our defect, fixed** — see below |

### The residue rule, and its prediction

The specification predicted the rule and its exclusion set. Both hold, measured.

A 4-base overhang donates its first base to close the upstream part's last codon and spells one
whole codon with the other three. So the last residue of every fragment must have a codon ending
in the donated base, and the residue the overhang spells overwrites what follows. At `AGAT` the
donated base is `A`, and the residues no codon of which ends in `A` are **C, D, F, H, M, N, W,
Y** — the specification's eight, derived independently. The prediction that a part is changed
exactly when its last residue is one of those is **right on 24 of 24 N parts**, and the same
test at the C junction, where `GCAT` donates `G` and the excluded set is **C, D, F, H, I, N, Y**,
is right on 24 of 24 DBD parts.

This is the equivalent of the source's own residue changes, and it is reached the same way:
both designs overwrite a boundary residue rather than insert one. The specification's section
3.1 prefers inserting a Gly and keeping every native residue; the planner has no such option,
which is recorded below.

### Two things the specification says that this run does not reproduce

Both were taken up by #224 on the same day. One correction landed; one turned out to be larger
than it read.

- **Section 9.5 says six of the source's 72 barcodes carry a stop at our frame offset.**
  Measured here: **none** do, at the retained stuffer's phase of 1, with our scar or theirs, in
  a bare barcode-plus-scar unit and in a whole 41-base block alike. The reason for drawing a
  fresh set stands on the distance rule and the site screen; this one measurement does not.
  #224 measured it a third time, reading the frame off `product.dna` rather than off the
  specification, and found none again: six is what the same 72 give one base out of frame.
  Section 9.5 now says so.
- **Section 6 says three positions is two iGGA rounds.** The planner reports three. This run read
  the difference as wording; #224 measured it and it is not. The planner's first round opens the
  destination with the internal enzyme, ligates the N part list released by the external enzyme,
  transforms it and sizes it for 24 products — the same reaction as rounds 2 and 3, and a
  24-member library at the end of it. The method's seating step is a different reaction: BsmBI,
  one well a part, 72 parts into a carrier, and no library made. **One round a position is what
  the pipeline models**, and the method's accounting needs a carrier the pipeline has none of.
  That is a design question, filed as the fifth item of #223.

### What the planner could not do

Five gaps, each measured on this run. The first is fixed; the rest are open.

1. **An enzyme reserved by a step outside the rounds had nowhere to be named.** The scheme took
   one internal enzyme, one external and the blunt choppers, so BsmBI — which seats a part in
   the DMX carrier — reached neither domestication nor the barcode rules. Measured: **four BsmBI
   sites survived in the 72 designed coding regions**, one of which a later cut would have
   destroyed. Fixed in #220: a scheme now takes a `reserved` list, held clear everywhere in a
   block and refused in the scheme's own stuffers.
2. **One overhang pair cannot serve every position.** The standard holds every junction's
   overhang to a pairwise distance rule, which is right for one pot and not for rounds: each
   round's tube holds one entry overhang and the scar, and no other pair ever meets.
3. **The external stuffers cannot move to the carrier.** D4 leaves our cargo with none, and the
   scheme requires both. The order sheet therefore carries 58 bp a part that this design would
   not buy.
4. **A junction cannot add a residue.** Its codons are always charged to a neighbouring part, so
   the designed protein never grows.
5. **A destination may not carry the external enzyme's sites.** The DMX backbone does, because
   the same plasmid is a donor in another tube. The demo therefore runs on a minimal destination
   rather than a DMX build.

Domestication still takes an enzyme list and not a motif list, so the polyadenylation screen of
section 3.3 was not run. That gap was already open.

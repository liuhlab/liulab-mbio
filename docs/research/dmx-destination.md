---
search:
  exclude: true
---

# The rebuilt DMX destination

What `scripts/build_dmx_vector.py` does to DMX0001, why each step is there, and what the result
measures. Written for #358, which is #331 absorbing #334 and #353; #374 added the marker swap
of §2.3. Every number below was produced by running the script; none was typed in from
somewhere else.

## 1. The parent

`reference_docs/synthesis_and_assembly/dmx/addgene/addgene-plasmid-247434-sequence-494299.dna`
— DMX0001, Addgene 247434, the depositor's own map, **5,839 bp circular**. The provenance table
in `docs/research/synthesis-and-assembly.md` §7 records the download. DMX0002 (Addgene 247435,
5,842 bp) sits beside it and measures the same way; the rebuild uses the first.

Measured with `mbio.sites`, not read off the map:

| | sites |
| --- | --- |
| BsaI | 3 — 355 F, 797 R, 1664 R (inside AmpR) |
| BsmBI | 2 — 367 R, 785 F |
| BbsI | 4 — 4130 F, 4490 F, 4878 R, 5217 R |
| PaqCI | 0 |
| PmeI | 0 |
| SrfI | 1 — 557 F, inside ccdB |

Both enzymes give up the same piece: `BsaI 362-792 (430 bp) AGGA/TTCC`, and
`BsmBI 362-792 (430 bp) AGGA/TTCC`. Those are the method's own cargo overhangs, which is what
makes the parent a destination already — and what lets the script find the cassette by digesting
rather than by a coordinate.

## 2. The five steps

### 2.1 The parent's cassette becomes the method's internal stuffer

`vector.released_cargo` reads the 430 bp piece off a BsaI digest, and `edits.replace` puts the
method's 34-base internal stuffer there. The parent's own BsaI sites at 355 and 797 are **kept**:
they are what releases the cargo at the end, and under the tube model of
`docs/adr/0014-one-record-one-tube.md` a destination carries them outboard on purpose.

The parent's BsmBI sites lay inside the replaced piece, so both go with it, and so does the ccdB
cassette and its SrfI site. The SrfI the result reads is the internal stuffer's own.

*The method page's step 2 says the rebuild keeps the ccdB cassette. It cannot: the piece a round
gives up is 34 bases, and ccdB is 306. Counter-selection on a round's destination is a separate
question from this record and is left open.*

### 2.2 A blunt site goes outboard of each releasing cut

The method asks for one PmeI site between each BsaI cut and the primer that reads a well, and no
more. The script inserts `GTTTAAAC` 12 bases before the cassette and 33 bases after it, locating
both off the cassette it has just read rather than off a coordinate.

`gate.check_dmx_vector` is what holds those two offsets honest. It computes the clear windows off
the record — each primer's footprint and the releasing cut it faces bound them — and on the
result they are **277..363 and 416..503**, with the two sites at 350 and 437 inside them.

This closes the "Where the added PmeI sites go" decision, which recorded that there is no free
space on one parent without sacrificing an annotated tag. There is. `6xHis` ends at 432 on the
rebuilt record and the right-hand window runs to 503, so a site at 437 sacrifices nothing. The
worry came from inserting immediately outboard of the BsaI site; the rule only asks for a site
between that cut and the primer's footprint, which leaves 87 bases to choose from.

### 2.3 AmpR leaves, KanR arrives, and the third BsaI site goes with them

Departure D11 forces this step. The final transfer moves cargo out of this vector and into a
working vector that is AmpR or CarbR, so the two must not share a marker. The marker is read as
well as plated: `igga.stages.selection_for` names the drug off the record's own marker, so an
AmpR destination would plate every round on the carbenicillin D11 says does not carry over.

**What goes in.** The `NeoR/KanR` of `pCR-Blunt II-TOPO`, in the same reference directory —
the Zero Blunt TOPO carrier, whose guide `bench-numbers.md` already sources this method's
50 µg/mL kanamycin from. Measured off that map and the bases under it: `[1236, 2031)` forward,
**795 bp**, `ATG` to `TGA`, no internal stop, and **no site of any of the six method enzymes**
in either orientation.

**Where it starts and ends.** Exactly where the parent's `AmpR` coding sequence did: `[1144,
2005)` on the record as the three steps above leave it, which is DMX0001's `[1524, 2385)`
carried through them. Both ends are a boundary one of the two deposited maps already draws, so
this step chooses no coordinate of its own.

**Which promoter drives it.** DMX0001's own. The parent draws `AmpR promoter` at `[2385, 2490)`
butted against the coding sequence it drives, and the end of it is the native *bla* ribosome
binding site: `AATATTGAAAAAGGAAGAGT` reads straight into `ATGAGTATTCAA`. Only the reading frame
behind that changes, so the new start codon keeps the old one's spacing from the binding site
and nothing is claimed about expression that the parent was not already carrying. The feature
is renamed `KanR promoter`, because a map that still said `AmpR` on a kanamycin vector is the
confusion this step exists to remove.

**Which orientation.** The parent's. `AmpR` reads on the reverse strand, so the coding sequence
goes in reverse-complemented and reads on the reverse strand under the promoter already pointing
at it. The script takes the strand off the feature it replaces rather than naming one, so the
step holds on DMX0002 too, whose `AmpR` reads forward.

**What it takes with it.** The parent's third BsaI site, at 1664, is the only site of any method
enzyme in the whole of `[1524, 2490)`, so the swap removes it. The synonymous codon that used to
stand in for the swap is no longer changed, and `domesticate` is left with nothing to do for
BsaI. The four parental BbsI sites — 4130, 4490, 4878 and 5217 — all lie outside the swapped
span, so the swap does not touch them and the next step still has all four to clear.

### 2.4 BbsI leaves the backbone

Four parental sites. Two (4878, 5217 on the parent; 4432, 4771 after the edits above) lie inside
`lacI`, and `mbio.sites.domesticate` takes both away with one synonymous codon each,
`GAC` to `GAT`.

The other two (4130, 4490 on the parent; 3684, 4044 after) lie in **no annotated feature at
all** — backbone between `rop` and `lacI`. `domesticate` reports such a site and leaves it,
because changing it changes what the record spells and no method decides that on its own. Here
the script decides it, and says so: one base each, `G` to `A` at 3684 and at 4044. It reads the
site left to right and A, C, G, T in turn and takes the first change that removes the site
without spelling a new site of any method enzyme, so the same parent gives the same record.

## 3. What the result measures

**5,393 bp circular**, named `DMX-iGGA`.

| | sites |
| --- | --- |
| BsaI | 2 — 363 F, 409 R |
| BsmBI | 0 |
| BbsI | 2 — 376 R, 392 F |
| PaqCI | 0 |
| PmeI | 2 — 350 F, 437 F |
| SrfI | 1 — 383 F |

- `BbsI digest` → 370-400 (30 bp) AGGA/TTCC: what a round opens it on.
- `BsaI digest` → 370-404 (34 bp) AGGA/TTCC, so `vector.released_cargo` is `Segment(370, 404)`.
- `KanR` is at `[1144, 1939)` reverse, with `KanR promoter` at `[1939, 2044)` reverse behind it.
- `bench.phenotype.selection_marker` reads that `KanR`, `stages.selection_for` names
  **50 µg/mL kanamycin** for it, and `stages.holes_for` raises no `H22`.
- `vector.destination_vector` takes it as it stands, with its outboard BsaI and PmeI.
- `vector.donor_cassette` takes it as a donor backbone, the same record in the other tube.
- `gate.check_dmx_vector` passes with two blunt sites, "in 277..363 or 416..503".
- `gate._read_by_a_well_primer` is `True`: both LevSeq primers bind, at 245 F and 883 R on the
  parent, flanking the cassette.

Every one of these reproduces the numbers #331's resolution comment predicted.

## 4. Why the record ships

The rebuilt vector is a plasmid the lab builds, not one that is deposited anywhere. The parents
are deposited and are already held. A DNA sequence carries no licence, so none was weighed —
an MTA or a UBMTA governs a plasmid, never its bases. #287 settled the same include-or-not for
the ccdB cassette and decided that synbio ships the method's own DNA, with the input vector
staying a build input. This record is the method's own molecule by that test, and it ships as
the worked example's destination, written by the script rather than by hand.

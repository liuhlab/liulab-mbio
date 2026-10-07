---
search:
  exclude: true
---

# The rebuilt DMX destination

What `scripts/build_dmx_vector.py` does to DMX0001, why each step is there, and what the result
measures. Written for #358, which is #331 absorbing #334 and #353. Every number below was
produced by running the script; none was typed in from somewhere else.

## 1. The parent

`reference_docs/synthesis_and_assembly/dmx/addgene/addgene-plasmid-247434-sequence-494299.dna`
— DMX0001, Addgene 247434, the depositor's own map, **5,839 bp circular**. The provenance table
in `docs/research/synthesis-and-assembly.md` §7 records the download. DMX0002 (Addgene 247435,
5,842 bp) sits beside it and measures the same way; the rebuild uses the first.

Measured with `liulab_mbio.sites`, not read off the map:

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

## 2. The four steps

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

### 2.3 BbsI leaves the backbone

Four parental sites. Two (4878, 5217 on the parent; 4498, 4837 after the edits above) lie inside
`lacI`, and `liulab_mbio.sites.domesticate` takes both away with one synonymous codon each,
`GAC` to `GAT`.

The other two (4130, 4490 on the parent; 3750, 4110 after) lie in **no annotated feature at
all** — backbone between `rop` and `lacI`. `domesticate` reports such a site and leaves it,
because changing it changes what the record spells and no method decides that on its own. Here
the script decides it, and says so: one base each, `G` to `A` at 3750 and at 4110. It reads the
site left to right and A, C, G, T in turn and takes the first change that removes the site
without spelling a new site of any method enzyme, so the same parent gives the same record.

### 2.4 The third BsaI site leaves AmpR

The remaining BsaI site is at 1284, inside `AmpR`, and `domesticate` takes it away with one
synonymous codon, `GGG` to `GGC`.

**The marker is not swapped, and that is a defect rather than a choice.** The method page's step 2
swaps AmpR for KanR, and names the third BsaI site as what that swap takes with it. One
synonymous base takes the same site, so the site is not the reason to swap. The bench is:
departure D11 is forced, because the final transfer moves cargo from this vector into a working
vector that is AmpR or CarbR, and the two must not share a marker.

The marker is also read now. `igga.stages.selection_for` names the drug a record's marker
selects and `stages.holes_for` answers `H22` wherever it can, so this record plates every round
on carbenicillin — the drug D11 says does not carry over. A KanR cassette is available offline:
`pCR-Blunt II-TOPO` in the same reference directory carries one, 795 bp and free of every method
enzyme. #374 carries the swap and the numbers it moves.

## 3. What the result measures

**5,459 bp circular**, named `DMX-iGGA`.

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
staying a project input. This record is the method's own molecule by that test, and it ships as
the worked example's destination, written by the script rather than by hand.

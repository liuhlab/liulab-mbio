---
search:
  exclude: true
---

# How our design departs from the paper

Every way the iGGA design in `docs/synthesis-and-assembly.md` differs from Takacsi-Nagy et al.,
what forced each change, and whether the reason survives checking. The paper itself, and what
was measured from it, are in `docs/research/protein-library-assembly.md`. The design's other
sources, and every decision still open, are in `docs/research/synthesis-and-assembly.md`.
Drafted by one agent, then checked by two more working apart: one re-derived every fact, one
attacked every reason. Corrections from both are folded in.

Eighteen differences were found. The sixteen below are checked and sound, with nothing left to
decide. D9, the barcode rules, is decided and lives in the method document. D16 is **withdrawn**:
the paper answers it and our design matches, so it was never a difference. D18 was added by the
decision in #259 on what bounds a round.

Pointers into the method document name its sections, so they survive an edit to it. Everything
here was simulated with `liulab_mbio`, not read off the page.

## D1. The two enzyme roles moved up one enzyme

Paper: the internal stuffer carries BsaI and opens the destination; the external stuffers carry
BsaI and release the donor. Ours: the internal stuffer carries BbsI (`### iGGA cargo`) and
BsaI releases the donor (`### iGGA pipeline`).

Why: the DMX plate-barcode kit is used as supplied (`### Lab resources (build once)`), and its
96 plasmids present
`AGGA` and `TTCC` cut by BsaI. That fixes BsaI at the cargo interface, which pushes the internal
stuffer onto a third enzyme. The DMX paper's own use of BsaI is not what forces it, since we
rebuild that backbone anyway.

## D2. BbsI in the internal stuffer, where the paper has BsaI

PaqCI needed NEB's activator oligo. It is a multi-site enzyme, so it has to engage two copies of
its site to cut one, and a substrate down to its last site stalls. Our destination started with
two and had one left after the first cut. That is 5 pmol of a quarter-microlitre reagent on the
one digest where a miss is invisible, since an uncut destination is KanR, carries a real cargo
and is valid input next round.

Simulated, BbsI does the same job: the opened destination is identical, `TTCC` and `AGGA` ends
with 16 and 12 bp stuffer stubs where PaqCI left 19 and 15, and SPRI removes either. The stuffer
drops from 38 bp to 34. The donor carries its two BbsI sites through its own BsaI + PmeI digest
untouched; the two digests never share a tube.

Domestication, measured: DMX0001 and DMX0002 each carry the same four BbsI sites, two in `lacI`
and two in unannotated backbone, so four changes once on a vector already being rebuilt.
pCR-Blunt II-TOPO and all four barcode plasmids carry none, so the KanR swap brings none in.

Cargo now has to be clear of a 6 bp site instead of a 7 bp one. The paper already pays this:
all 72 of its domains, 23,103 bp including one of 1,046 bp, carry zero BbsI, zero BsaI, zero
SrfI and zero PmeI, against 11 BsmBI and 15 PaqCI it never had to remove. Recoded for E. coli
the cost is about 1.5 extra synonymous changes per 2.5 kb. Non-coding fragments are unaffected
by the choice; `### Cargo synthesis` step 2 owns them.

A destination that escapes BbsI is still linearised by SrfI into blunt ends, so the parent
could only return if the ligase joined blunt. The T7 ligase in `## Reagents and equipment` does
not.

## D3. One overhang pair every round instead of one per position

Paper: four overhangs, `CTCC` into the vector, `GGAG` into position 2, `CCGA` into position 3,
`AGCG` at every far end. A part's internal stuffer presents the next position's overhang, so an
N part can only be followed by a bZIP part. Order is fixed by the DNA, and it is a cycle: the C
part presents `CTCC` again, so a fourth round with an N pool would work. Three rounds is
protocol, not sequence.

Ours: `AGGA` and `TTCC` at every round (`### iGGA cargo`). Order is set by which pool goes in
the tube (`### iGGA pipeline` step 1).

Why: `AGGA` and `TTCC` are the published SAPP standard (`### Parts`). Position-specific
overhangs would need a different DMX vector per position.

The pot is safe. The four ends present are `TCCT` and `TTCC` from the opened destination and
`AGGA` and `GGAA` from the released donor. The only complementary pairs among them are the two
intended ones, so a donor cannot go in backwards, two donors cannot chain, and neither molecule
closes on itself. Simulated, one circular product.

Any pool fits any round once every fragment obeys the length rule (`### iGGA cargo`). That constraint is
uniform across pools rather than per position, which is why ours needs no terminal molecule (D7).

## D4. Our cargo has no external stuffers

Paper: every synthesized block is 29 bp external stuffer, CDS, internal stuffer, 11 bp barcode,
29 bp external stuffer, with BbsI and PmeI in the external stuffers. The C block also carries a
54 bp T2A. Ours is `[AGGA]─[fragment]─[internal stuffer]─[11 bp barcode]─[TTCC]` and nothing
else (`### iGGA cargo`); the BsaI sites belong to the DMX cassette (`### DMX vector`).

Why: the DMX vector already supplies BsaI in the right places, so external stuffers would
duplicate it. The oligo-pool budget is the weaker reason. 58 bp is about a quarter of one
300-mer, roughly 14% of a pool over 10,000 fragments. Real but not decisive.

A consequence worth naming: the paper's donor carries its own release sites, so a donor can be a
PCR product or a gene fragment. Ours cannot be a donor until it is cloned.

## D5. PmeI sits in the vector, not in the cargo

The paper carries PmeI in each synthesized block's external stuffers, measured at offset 0 of
every 5' stuffer and offset 21 of every 3' stuffer, across all 72 blocks, with none anywhere
else in a block. Ours sits in the DMX vector, outboard of each BsaI site (`### DMX vector`).

Same order of sites either way: PmeI outermost, then the enzyme that releases the donor, then
the cargo's overhang. Only the molecule differs, because D4 left us with no external stuffers to
put it in.

Why it matters. PmeI's job is to strip the donor backbone of the overhangs that would let it
recapture its own cargo. Cut either layout and the backbone comes out blunt at both ends, and
the only pieces that can bridge it back to the cargo are stubs small enough for a 2X SPRI:

```text
paper,  BbsI + PmeI          ours,  BsaI + PmeI
 2408 bp  BLUNT / BLUNT       1408 bp  BLUNT / BLUNT    backbone, cannot touch cargo
  262 bp  GGAG  / AGCG         199 bp  AGGA  / TTCC     cargo
   25 bp  AGCG  / BLUNT         17 bp  TTCC  / BLUNT    stub
   21 bp  BLUNT / GGAG          13 bp  BLUNT / AGGA     stub
```

Ours is cheaper. They pay 58 bp of external stuffer on every one of 72 blocks to carry PmeI and
BbsI; we pay once, in the backbone, and every cargo is covered.

The paper never states the reason. PmeI appears in its key resources table and in one methods
sentence pairing it with BbsI, and nowhere else; Figure S1A's
legend names only BbsI and BsaI. The placement is theirs and measured from Table S1. The reason
is ours, derived by digesting both layouts.

## D6. Our cargo must be clear of five enzymes

Paper: every one of the 72 published blocks contains, by design, 2 BsaI, 2 BbsI, 1 SrfI and
2 PmeI. Ours must be clear of BsaI, BsmBI, BbsI, SrfI and PmeI outside its own stuffer
(`### Vectors and parts summary`).

Why: our cargo passes through more hands. BsmBI puts it into the DMX vector, BsaI takes it out,
BbsI and SrfI run the rounds, PmeI kills the donor backbone. Cost: a non-coding fragment cannot
be recoded at all, which `### Cargo synthesis` step 2 already hands back to the user.

## D7. The paper needs a position-specific terminal molecule, we do not

Paper: every block must present its graft point, the first base of the BsaI overhang buried in
its internal stuffer, two bases into a codon. N and bZIP get there by sizing the insert.
C cannot, because it carries the 54 bp T2A and ends the protein, so its insert and its T2A are
both whole codons. It makes up the two bases with a `GG` pad in front of the stuffer's overhang,
and that pad is the whole of the 58 bp against 56: measured, the three stuffers are identical
from the overhang onward. Delete the pad and a stop lands at codon 366.

So the C pool is a physically distinct molecule, sized by where it sits in the chain. It is not
terminal, though. Its graft overhang is `CTCC`, which is what an N block enters on, so a fourth
round would work (D3).

Ours: one stuffer, identical in every position, its overhang at its first base, no pad, retained
in every round. There is no terminal form.

Measured on a finished three-round plasmid with the stuffer retained:

```text
BbsI + SrfI   2392 bp  TTCC / AGGA     reopens, a fourth round fits
BsaI          1438 bp  TTCC / AGGA     backbone
               984 bp  AGGA / TTCC     cargo, into the working vector
```

So our chain is a true cycle, as theirs is, and the finished product is still releasable.

Why ours can be uniform: one fragment rule (`### iGGA cargo`) serves every position, because no
fragment
of ours ends the protein the way the paper's C block does, and the stuffer sits at the same
phase every round. That is what makes D3's any-pool-in-any-round real rather than nominal, and
it is the one place our design is simpler than the paper's rather than merely different.

## D8. No T2A, so the barcodes are read inside the protein

Paper, measured as offsets from the N-domain ATG in the assembled product:

```text
N_cds 0   bZIP 753   C_insert 948   T2A 993   C_stuffer 1047
BC_C 1105   BC_bZIP 1120   BC_N 1135   first stop 1170
```

The 54 bp T2A sits between the C insert and the stuffer, so the stuffer and all three barcodes
come off as a separate 41-residue peptide and the transcription factor ends at the 2A.

Ours has no T2A anywhere, and the cargo sits in one ORF (`### Working vector`), so the whole
barcode
block is translated inside the fusion protein. Three rounds is 45 bp, fifteen residues.

Why: forced. A T2A inside the cargo would end translation before the C-terminal part, so every
C-terminal tag, degron or localisation signal would be lost. The barcodes have to accumulate at
the assembly junction, which is inside the cargo, so there is nowhere outside the ORF to put
them.

Cost: fifteen residues of per-member sequence in the middle of every protein, and any stop codon
in a barcode truncates it. Their barcodes are translated too, and all 72 are stop-free at their
own phase, which is offset 2 into the barcode. At ours, offset 0, six of the 72 carry a stop, so
their set is not reusable. Ours truncates the protein; theirs truncates a junk peptide.

A last round can still put a T2A in, by swapping the stuffer for an `AGGA`/`TTCC` block that is
1 mod 3 counted from its `AGGA` (`### iGGA cargo`). `AGGA` plus the paper's 54 bp T2A is 58,
which is
1 mod 3 with no padding. That ends the chain, since the product then carries no stuffer to
reopen.

## D10. The DNA comes from an oligo pool

Paper: each domain codon-optimized and individually synthesized as a clonal plasmid in pTwist
Kan, then pooled equimolar per family. Measured CDS lengths 27 to 1,046 bp.

Ours: fragments come from a Twist pool, assembled by DAD, split into wells by DMX and read by
ONT (`## Goal`, `### Cargo synthesis`, `### Cargo validation (optional)`). The oligo length is a
project input, not a method constant.

The budget fits, and the figure to use is the method page's, not this one. An earlier revision
read "about 235 bp usable per 300-mer" off a two-primer layout; the three-primer scheme spends
60 nt on primer sites and 14 nt on recognition sites and spacers, so a 300-mer carries a 226 nt
span and a 350-mer a 276 nt one. `ap1-demo-project.md` §4.2 derives it and #257 settled where it
lives. What ours carries and the paper's did not is a per-clone sequence-and-reformat step
before assembly can start.

## D11. Both markers had to change

Paper: domain plasmids are pTwist Kan. The destination in round 1 is itself a domain plasmid, so
donor and destination share a marker there too. The paper never states this; it follows.

Ours: the DMX vector becomes KanR by replacing the AmpR of DMX0001 or DMX0002
(`### Lab resources (build once)`), and the working vector it feeds is AmpR or CarbR
(`### Vectors and parts summary`).

Why: forced. Measured, DMX0001, DMX0002 and all four barcode plasmids are AmpR as supplied. The
final step moves cargo from the DMX vector into the working vector, so the two must not share a
marker. One of DMX0001's three BsaI sites sits inside AmpR, so the swap removes it for free.

## D12. The last transfer is a different reaction

Paper: the assembled library goes into the CAR vector by the same two-digest chemistry, using
that vector's own internal stuffer. No counter-selection; `ccdB` appears nowhere in the paper.

Ours: staged in one tube. BsaI and PmeI release the cargo, heat kills them, then PaqCI and ligase
place it in a working vector whose ccdB cassette the cargo replaces
(`### Final assembly`).

Why: the cargo already presents `AGGA`/`TTCC`, which is the working vector's interface. Both
byproducts are covered. The released DMX backbone is KanR on an Amp plate, and released ccdB
kills anything that takes it. The document never spells that out.

## D13. Different host, and a layer the paper has no counterpart for

Paper: the product is a repair template for knock-in at the human TRAC locus in primary T cells.
Ours goes into a bacterial expression vector, with promoter and RBS either in the backbone or in
the N-terminal part (`#### Input working vector`), and libraries in BL21(DE3)
(`## Reagents and equipment`).

Consequence: ours needs the whole working-cassette layer, reusable N- and C-terminal parts with
their own overhangs and their own carrier vector (`### Parts`). The iGGA cargo has to
present `AGGA`/`TTCC` cleanly to drop into it.

## D14. A second barcode system, and BsaI with three jobs

Paper: one 11 bp barcode per domain plus one in the knock-in vector. BsaI does one thing.

Ours: the same per-fragment barcode plus the DMX plate-barcode kit, 96 plasmids in four groups
of 24 (`### Vectors and parts summary`). Measured, they chain exactly as that row says: group 1
gives `AGGA`/`GTTC`,
group 2 `GTTC`/`CCTT`, group 3 `CCTT`/`TCAG`, group 4 `TCAG`/`TTCC`. So BsaI on the DMX vector
does three jobs: plate barcoding, releasing cargo as the iGGA donor, releasing cargo into the
working vector.

These do not conflict. Plate barcoding runs in heat-lysed cell lysate on a throwaway aliquot and
the product is PCR-amplified for nanopore, so it never has to grow. The culture in the well is
the stock and is untouched. Our document does not say this, and it belongs in the DMX pipeline
steps.

We have no equivalent of the paper's vector barcode, so our final read has fewer barcode
positions than theirs.

## D15. Extra strains, forced by ccdB

The iGGA rounds themselves are unchanged: Endura, recover and grow at 30 °C every round
(`### iGGA pipeline`; the temperature itself belongs to the detail protocol, not the method
page), exactly as the paper. The extra strains belong to the DMX entry step and the
final transfer, which are layers the paper has no counterpart for. A ccdB-sensitive strain
and NEB Stable for anything carrying ccdB, both in `## Reagents and equipment`.

## D16. Withdrawn: the paper sequences whole plasmids too

Carried first as a difference and then as an open question, because the method document and an
earlier draft of this note were describing different objects. The paper settles it: "The
sequences of a sample of individual clones from both libraries were confirmed by whole plasmid
next-generation sequencing (Plasmidsaurus)" — a sample of a finished library, not a check on
each round's product, and not every member.

Ours does the same, so nothing departs. A round's product is a pool and holds no clone to read.
The objects that get a whole-plasmid read are the validated part wells and any single construct
archived as a stock.

## D17. The working vector's stuffer is the parent's own EGFP, not an RFP

Method page: `Stuffer cassette: RFP` (`#### Input working vector`). Ours: the EGFP already
carried by pLVX-TetOne-Puro-GFP, the parent #230 measured. The paper has no counterpart — it
builds no working vector and its destination carries no stuffer cassette — so this departs from
the method page rather than from Takacsi-Nagy et al.

Why: the transduction marker is mCherry (`docs/research/ap1-demo-project.md` section 7.1). A red
stuffer against a red marker means colour cannot separate a vector whose cassette was never
replaced from a cell that was simply transduced. The EGFP is being cut out of the parent anyway.

Measured on the parent, 2026-10-06, across the 720 bp the method replaces, both strands: **0
BsaI, 0 BsmBI, 0 BbsI, 0 PaqCI, 0 SapI**, so it sits between the cassette's two BsmBI sites
without being cut in the middle. It translates `MVSKGEELFTGV…GMDELYK*`, one stop and it the last
codon. It carries **no `AATAAA` or `ATTAAA` on either strand**, so retaining it reintroduces
nothing the poly(A) requirement is about. It is already under the TRE3GS promoter in the right
orientation, with the SV40 poly(A) signal beyond it.

Cheaper to build, too: an RFP stuffer is sourced and inserted, while the EGFP needs only its two
flanks replaced — `TATG` and a BsmBI site at one end, a BsmBI site and `CTAA` at the other —
which is one PCR off the parent with tailed primers.

**The caveat, and it applies to the RFP equally.** The stuffer sits under a TRE promoter, which
is a minimal CMV promoter and is not read in *E. coli*. **No stuffer colour screens a bacterial
colony**, green or red. The colour reads in the transduced cell after induction: green means the
cassette was never replaced, red-only means it was. A colony-level screen needs a bacterial
promoter on the stuffer or a counter-selection marker, which is a separate decision. Nothing
read measures TRE3GS activity in *E. coli*.

Every count above is from `docs/research/working-vector-plvx-tetone.md` section 5, and none is
recomputed here.

## D18. Every round is plated; the paper plates none

Paper: digest, SPRI, ligate, SPRI, electroporate, recover one hour at 30 °C, grow 12-16 hours,
prep, next round. No colony count appears anywhere in it. Ours plates a dilution of each round's
recovery, and a no-donor control carried through the ligation from the same digest, and gates
the round on net colonies (`### iGGA pipeline`, step 5).

Why. The design's correctness rests on the parental destination not surviving: BbsI cuts twice,
SrfI cuts whatever BbsI missed, T7 refuses the resulting blunt ends, and the 2X SPRI takes the
stubs. The control plate measures that chain rather than assuming it. What survives is not
noise — it is the previous round's library, one position short, and it would reach the linkage
read as truncated members. A round that collapsed cannot be repaired by a later one.

Cost: two plates a round, both growing during the outgrowth, so no extra day.

What it does not change. The library is still never read per member, and no colony is ever
picked from these plates. It goes to a pooled screen as a pool. The plate bounds a round; the
representation and linkage reads judge the library.

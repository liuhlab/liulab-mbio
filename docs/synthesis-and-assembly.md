# Highly Parallel DNA Synthesis and Assembly

This is the strategy the pipeline is built against. Design specifics are settled when the code
is written and the first project runs. An item marked undecided here is work the package will
own, not work forgotten. The evidence behind each choice, and every open decision, is in
`docs/research/synthesis-and-assembly.md`. One project that answers each open decision for
itself is `docs/research/ap1-demo-project.md`.

**This page covers two experiments.** **DAD-GGA-DMX** builds cargo from an oligo pool and reads
each well back, so one design can be named at a time. **iGGA** mixes those parts in rounds and
gives a pool. Nobody picks a member of that pool, and nothing in it is read for identity again.
Reads across the whole pool judge it. Keep the two apart: a rule that suits picked wells does not
suit a pool. `CONTEXT.md` defines both.

## Goal

1. Synthesise fragments up to 2.5 kb, at a scale of 10³ to 10⁴ designs, from the cheapest oligo
   pool a vendor sells.
2. A working vector that turns one cargo into any application, in a single Golden Gate reaction.
3. Iterative Golden Gate (iGGA) that chains cargo with a barcode at every junction, in a planned
   or a random order.

## Vectors and parts

What each molecule is. What is done with it is in `## Experiment sections`.

### Vectors and parts summary

Two tables, one row per molecule. `Sites to remove` means sites outside that molecule's own
designed cassettes. `Cut by` names the enzyme that acts on the molecule and what the cut does,
wherever the site itself sits.

#### What each molecule is

| Item | Source and purpose | Selection | Common overhangs |
| --- | --- | --- | --- |
| Part carrier vector | pCR-Blunt II-TOPO; stores each part | KanR | — |
| Parts | reusable in-frame elements: tag, linker, signal peptide, localization signal, degron | — | N-term pair TATG/AGGA; C-term pair TTCC/CTAA |
| DMX vector | holds cargo | KanR | AGGA, TTCC |
| DMX barcode kit | 96 barcodes in four groups of 24; one from each group marks a well for ONT sequencing | AmpR | chain AGGA→GTTC→CCTT→TCAG→TTCC, group 1 to 4 |
| Cargo | Twist oligo pool fragments, via DAD and DMX | — | AGGA, TTCC |
| iGGA cargo | cargo built for iterative random assembly | — | AGGA, TTCC |
| Working vector | chosen backbone plus working cassette; gives cargo its application | AmpR/CarbR | TATG, CTAA (cassette); AGGA, TTCC (cargo) |
| Final vector | working vector with cargo installed | AmpR/CarbR | AGGA, TTCC (cargo junction) |

#### What cuts each molecule

| Item | Sites to remove | Cut by |
| --- | --- | --- |
| Part carrier vector | BsmBI | — |
| Parts | BsaI, BsmBI, the cargo enzyme | BsmBI (release the part) |
| DMX vector | BsaI, BsmBI, BbsI | BsmBI (insert cargo), BsaI (add a plate barcode, release cargo, release the iGGA donor), PmeI (blunt the donor backbone) |
| DMX barcode kit | — | BsaI (release a barcode) |
| Cargo | BsaI, BsmBI, the cargo enzyme | BsmBI (insert into the DMX vector) |
| iGGA cargo | BsaI, BsmBI, BbsI, SrfI, PmeI, the cargo enzyme | BsmBI (insert into the DMX vector), BsaI (release as the donor), BbsI (open the destination), SrfI (cut a destination BbsI missed) |
| Working vector | BsmBI, and the cargo enzyme | BsmBI (insert the working cassette), the cargo enzyme (insert cargo) |
| Final vector | — | — |

`docs/research/synthesis-and-assembly-barcode-kit.md` has all 96 barcodes.

#### The cargo enzyme is chosen per working vector

Every other enzyme above is fixed by something physical: parts carry BsmBI flanks, the DMX
vector carries BsaI outboard of BsmBI, the barcode kit releases with BsaI. The working vector
is the user's own, so its cargo enzyme is not fixed and is searched for per vector. The search
reads the working vector alone and runs before any block is designed, because every block has to
be free of what it returns.

What binds the choice:

- Distinct from BsmBI, so inserting cargo does not re-open the cassette
- A 4 nt 5' overhang, so `TATG`/`AGGA`/`TTCC`/`CTAA` are unchanged
- **No site in any molecule sharing its reaction, outside what it is meant to cut** — the
  reaction-scoped rule, which the gate reads
- A reaction and cycling protocol the package ships

Three other molecules share the final-assembly pot and so inherit the choice: cargo, iGGA cargo,
and any part in the working cassette. A 7 bp cutter costs them little — PaqCI's site is about
four times rarer than the 6 bp sites already on their lists.

The ccdB cassette's inner sites are whichever enzyme this returns, by synthesis or by the same
re-tailoring PCR a part uses. Where the search returns nothing, the sites that blocked each
candidate are the finding; see `docs/research/lentiviral-tolerance.md` for one vector class
where the answer is already known.

### Parts

#### Form and carrier

- Synthesis or proofreading PCR, flanked by inward-facing BsmBI sites and a default overhang: `[BsmBI.<5' overhang>]─part─[<3' overhang>.BsmBI]`
- Three cassette junctions: `TATG`/`AGGA` N-term part, `AGGA`/`TTCC` ccdB cassette and cargo, `TTCC`/`CTAA` C-term part
- Each part holds one pair, either the N-term or the C-term
- No part carries a BsaI or BsmBI site outside its own designed cassette
- Carrier: pCR-Blunt II-TOPO, whose backbone carries no BsmBI site
- A part has two usable forms: the carrier plasmid, and a PCR product carrying a different BsmBI flank, `5'-[spacer]-CGTCTC-N-[new overhang]-[annealing]-3'`

#### Rules a part obeys

- In frame, a multiple of 3 bp long, counted including its 4 bp overhang at 5'
- No stop codon
- Every overhang reads as one base closing a codon plus one codon: `T|ATG`, `A|GGA`, `T|TCC`, `C|TAA`
- The junction adds two amino acids

### Working vector

#### Input working vector

- A backbone with no BsmBI and no cargo-enzyme site, AmpR/CarbR, carrying `[TATG.BsmBI]─[stuffer cassette]─[BsmBI.CTAA]` at the cassette site
- Any other enzyme's sites are counted and recorded, not removed. A site only has to go where its enzyme shares a reaction with this vector
- Stuffer cassette: RFP
- The `ATG` of `TATG` is the start codon where the backbone carries a promoter and RBS before it. Otherwise the start sits in an N-terminal part

#### Working cassette

- `[TATG]─[N-terminal parts]─[ccdB cassette]─[C-terminal parts]─[CTAA]`
- Either side is optional, and the ccdB cassette comes in a version to match each case
- ccdB cassette (for later cargo insertion), supplied as a part in four versions, where `CE` is
  the cargo enzyme chosen for this working vector:
  - N and C parts: `[BsmBI.AGGA.CE]─[promoter+ccdB]─[CE.TTCC.BsmBI]`
  - no N part: `[BsmBI.TATG]─GG─[AGGA.CE]─[promoter+ccdB]─[CE.TTCC.BsmBI]`
  - no C part: `[BsmBI.AGGA.CE]─[promoter+ccdB]─[CE.TTCC]─GG─[CTAA.BsmBI]`
  - neither: `[BsmBI.TATG]─GG─[AGGA.CE]─[promoter+ccdB]─[CE.TTCC]─GG─[CTAA.BsmBI]`
  - All four expose `AGGA`/`TTCC` for the cargo; a skipped side adds one Gly
- The cassette's promoter is `J23119` with a `B0034` RBS, so it counter-selects in any
  ccdB-sensitive strain. The DMX cassette's T7 arrangement below does not

### DMX vector

- The modified DMX vector has no BsaI, BsmBI, BbsI or SrfI outside the ccdB cassette, carries one PmeI site outboard of each BsaI site, uses KanR
- The PmeI site blunts the donor backbone in iGGA; the DMX pipeline does not use it
- DMX ccdB cassette: `[PmeI]─[BsaI.AGGA.BsmBI]─[promoter+ccdB]─[BsmBI.TTCC.BsaI]─[PmeI]`
- The ccdB cassette reads from a T7 promoter under a *lac* operator, with the repressor on the
  plasmid. Only a strain carrying T7 polymerase counter-selects
- So a DMX stock grows in an ordinary strain because its toxin is silent there, not because the
  strain resists ccdB. Those are two different things, and only the first one holds here

### Final vector

- A working vector whose ccdB cassette has been replaced by cargo, joined at `AGGA`/`TTCC`
- AmpR/CarbR, where the DMX backbone it was released from is KanR

### iGGA cargo

Runs in the DMX vector; the working vector takes only the finished cargo.

- iGGA cargo: `[AGGA]─[fragment]─[internal stuffer]─[11 bp barcode]─[TTCC]`
- Each fragment, counted from the first base of its `AGGA`, is a multiple of 3
- Internal stuffer `[AGGA.BbsI]─[SrfI]─[BbsI.TTCC]`, 34 bp: `AGGAAAGTCTTCAGCCCGGGCAGAAGACAATTCC`
- Barcode: 11 bp, minimum Sequence-Levenshtein distance 3 within a part list, no in-frame stop codon, no homopolymer run over 5, and no site of the five cargo enzymes on either strand, read with the scar either side of it
- The stuffer is retained every round
- A capping block is any `AGGA`/`TTCC` block that is 1 mod 3 counted from its `AGGA`, such as T2A. It ends the chain

`docs/research/barcode-design.md` is the evidence for each barcode rule.

## Experiment sections

Major steps, the decisions attached to them, and anything designed that must be followed
exactly. Volumes, temperatures, times, catalogue numbers and anything a person follows
without deciding belong to the detail protocol.

Every plasmid is confirmed by whole-plasmid sequencing before it becomes a stock. Only a
pooled screen skips it.

### Lab resources (build once)

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 1 | Clone each part into the carrier | BsmBI-flanked part | part plasmid | which overhang pair the part gets |
| 2 | Rebuild the DMX vector backbone | DMX0001 or DMX0002 | KanR DMX vector, BbsI-free, PmeI added | which parent |
| 3 | Build the input working vector | chosen backbone | backbone carrying the stuffer cassette | the backbone; where the cassette goes; where promoter and RBS sit |
| 4 | Sequence each candidate whole | candidates from steps 1 to 3 | verified stocks | which clone becomes the stock |
| 5 | Lay out the primer plates | orthogonal primer set | source plates, then pre-made combination plates | which primer fills which slot |
| 6 | Stock the DMX barcode kit | DMX barcode kit | 96 plasmids in four groups | — |

**Step 5.** A slot always means the same primer pair. A gene whose own sequence contains its
slot's primer site moves slot.

### Working vector (per application)

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 1 | Choose the cassette's parts | part collection | N-terminal and C-terminal parts, or neither | what each side carries, and whether either is skipped |
| 2 | Take the matching ccdB cassette | four versions | the one fitting which sides are used | the cargo enzyme, and so which cassette set |
| 3 | Retailor an overhang | part plasmid | BsmBI-tailed PCR product | whether a part's default pair fits this cassette |
| 4 | Assemble | input working vector, parts | working vector | — |
| 5 | Screen and pick | transformants | candidate clones | — |
| 6 | Sequence each candidate whole | candidate clones | working vector stock | which clone becomes the stock |

Check total length and part compatibility before assembling, and check that the two amino
acids each junction adds do no harm in this construct. A part that fits its default pair is
added as the carrier plasmid, with no digest and no PCR of its own.

**Undecided:** what the carrier's own BsaI site and ccdB do to a part added that way — one BsaI
site, inside the counter-selection cassette. See `docs/research/synthesis-and-assembly.md`.

### Cargo synthesis

Builds full-length cargo from an oligo pool and puts it in the DMX vector. From PCR2 onward,
identity is held by well position.

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 1 | Collect input sequences | protein or DNA sequences | the design list | which genes, and how long |
| 2 | Check compatibility | design list | cleared list | whether to recode or drop a gene carrying a restriction or primer site |
| 3 | Split each gene with DAD | cleared list | fragments, internal overhangs | oligo product and length |
| 4 | Draw the barcode set — iGGA cargo only | part lists | one barcode per fragment | — |
| 5 | Assign primer slots | fragments | oligo layout | batch size |
| 6 | Order the pool | oligo layout | oligo pool | pool size |
| 7 | Run PCR1, per batch | pool | one batch's pieces | — |
| 8 | Run PCR2, per gene | batch | one gene's pieces | — |
| 9 | Assemble into the DMX vector | pieces, vector | cargo in vector | — |
| 10 | Transform and archive as cells | assembly product | glycerol archive, one design per well | polyclonal, or go on to validation |

**Step 3.** `AGGA` and `TTCC` are held out of the overhang set DAD may choose.

**Steps 3 and 5, the three-primer scheme.** Every oligo is
`[P1][BsmBI][span][BsmBI][pad][P2][P3]`. PCR1 uses the pair (P1, P3) shared by a whole batch of
genes and pulls that batch out of the pool. PCR2 uses (P1, P2) to pull one gene out of the batch
and drops the P3 region.

**Steps 3 and 5, the two BsmBI sites.** Both face inward and are part of the ordered sequence.
The two outermost cuts across a part give `AGGA` and `TTCC`; every internal cut gives the
overhang the split chose. They cannot move onto the primers: all of a gene's fragments share one
P1 and P2, so a primer would give every fragment the same overhang.

**Steps 3 and 5, the budget.** The span between the two cuts carries its own overhang at each
end, and each internal overhang belongs to two fragments at once. A gene of `L` bases needs:

```text
span  ≤ oligo length − 60 (three primer sites) − 14 (two recognition sites and spacers)
pieces = ceil((L − 4) / (span − 4))
```

At 350 nt the span is 276 bases and holds 268 between its overhangs. **Both numbers are
computed, never typed:** count the overhangs once and a 1,149 base part needs five pieces, not
six.

**Step 5, padding.** Every oligo is padded to the same length, which the build sets. Filler
goes between the second BsmBI site and P2, outside the cut, so it never reaches the product. It
is screened for the reserved enzyme sites, for the primer sites, and at both of its joins.

**Step 4.** The set is drawn fresh, to the rules in `### iGGA cargo`.

**Step 5, batch size.** At most 96 genes a batch — one inner-primer plate, one PCR2 plate — and
a larger library divides into equal batches, never full batches plus a remainder. Hold pieces
per PCR1 tube roughly constant and let batch size fall as genes lengthen; the budget is a
build input, inside the 96 to 768 pieces the Baker anchors derive.

**Step 10.** Polyclonal by default.

**The archive.** Cells, not DNA, in duplicate plates. DNA is prepared when a set is drawn, one
prep for the set. No per-well normalisation.

### Cargo validation (optional)

The cargo sits in the DMX vector either way; this section only reads each well. Identity stays
with well position throughout. Choose the route per build. Each route states its own capacity
below, and a build that reads designs back names one.

**Which designs are read.** None, unless the build asks. A library headed for a pooled screen
takes its identity from that screen, so most libraries are never read back one design at a time.
A build that wants the read names a fragment count. Every design in that many pieces or more
is read, and the rest stay polyclonal. Set the count to zero to read every design. The chance of
a clean colony falls as a gene is cut into more pieces, so this is the line the build is really
drawing. The protocol prints each design's own chance beside it.

#### Getting to clonal wells

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 1 | Spot the archive as an array | archive plate | one spot per design | spots per tray |
| 2 | Grow | tray | colonies per spot | — |
| 3 | Pick into 384-well plates | tray | clonal cultures, position kept | colonies per design |
| 4 | Grow | 384-well plates | cultures ready to sample | — |

The picked plate is where the shared section ends, on both routes. Each route's own first step is
the move out of it: barcode ligation takes four of these plates into one 1536-well plate, index
PCR takes one into 96-well plates. Neither route amplifies or lyses in the plate the colonies were picked into.

**Step 3, colonies per design.** Four by default, which is Lund's anchor and the only measured
one, and a build may pick more. Four colonies gave a clean copy of 343 of 458 genes. The same
work shows why four is not always right: a gene in two pieces came out clean every time, one in
five pieces 84.6% of the time, one in twelve 40%, and one in sixteen never. The protocol prints
each design's own chance from that table, so a build that needs more picks knows it before the
plates are poured.

#### Barcode ligation

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 5 | Compress four 384-well plates into one 1536-well plate | cultures | cell lysate per well | — |
| 6 | Add one barcode from each of the four groups | DMX barcode kit | a combination per well | — |
| 7 | Barcode in lysate | lysate, barcodes | barcoded construct per well | — |
| 8 | Pool, then amplify the three primer pairs separately | barcoded pool | three amplicon pools | — |
| 9 | Sequence | amplicons | reads | — |
| 10 | Basecall and demultiplex | reads | a sequence per well | the depth floor |
| 11 | Reformat | per-well calls | compacted plate | — |

#### Index PCR

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 5 | Sample one picked plate into 96-well plates | picked plate | a plate of reactions per 96 wells | — |
| 6 | Amplify each well with one barcoded primer pair | reactions, prepared primer plate | barcoded amplicon per well | — |
| 7 | Pool per plate and clean up | amplicons | one pool per plate | — |
| 8 | Sequence | pool | reads | — |
| 9 | Demultiplex | reads | a sequence per well | the depth floor |
| 10 | Reformat | per-well calls | compacted plate | — |

**Step 5.** One quarter of the picked plate goes to one 96-well plate. A 384-well plate's wells
sit at half a 96-well plate's spacing, so one well in four lines up under a standard multichannel
head, and the plate is covered in four passes. Picking fills one quarter at a time for the same
reason: a part-filled picked plate then gives full 96-well plates, no reverse barcode is spent on a
plate that is mostly empty, and a design's colonies stay together on one plate. An acoustic handler
does the same move where one is already booked, which is how barcode ligation moves out of the
same plate.

**Step 6, the primer plate.** The barcoded pairs are built once and kept as lab stock, by their
own preparation protocol. A run calls for a prepared plate and does not build one.

#### Which barcodes mark a well

The well's address is split across the barcode sets, and one set carries the plate. Index PCR
uses 96 forward barcodes for the well and 96 reverse for the plate, which reaches 9,216 wells on
the 192 primers already held. Barcode ligation uses three DMX groups for the well and the fourth
for the plate, so one barcode goes across a whole plate from a reservoir.

Two plates on one flow cell are told apart this way. The combination is worked out from the
well, not looked up in a file, so a demultiplexer can check an address instead of trusting one.

#### What counts as a pass

Two questions, in order.

Is the read deep enough to call? If not, the well has **no verdict**. It is read again or picked
again, and it is not a failure. The floor comes with the route, because each one was measured on
its own. Barcode ligation wants more than 150 reads. Index PCR wants more than twenty, and a well
at ten reads or more is still called, with a warning. The wanted depth is one to pass; the tolerable one is a
depth to reach. A build may raise the wanted depth. It cannot move the tolerable one, which
comes with the route.

Does the call match the design exactly — both entry overhangs, the fragment, the stuffer and the
barcode? Anything less fails, a silent change included, because a barcode that no longer names
its member cannot be put right later. A well with more than one consensus is mixed, and mixed
fails.

#### After either route

Reformatting compacts out the wells that failed. A well with no verdict is not compacted out.
A design with no passing well is picked again from the same archive spot before anyone
re-synthesises it. The polyclonal route skips this section.

### iGGA pipeline

Runs in the DMX vector. iGGA cargo is one kind of cargo and may be polyclonal or validated like
any other.

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 1 | Choose the pools and their order | one part list per round | a round plan | which list goes first |
| 2 | Digest the destination and the donor apart | destination pool, donor pool | opened destination, released donor | — |
| 3 | Ligate | cut ends | circular product | — |
| 4 | Electroporate, recover and grow | ligation | the round's pool | — |
| 5 | Plate a dilution, and a no-donor control | the recovery | the round's colony count | run the round again, or carry it on |
| 6 | Prep the pool | the round's pool | next round's destination | — |
| 7 | Repeat steps 2 to 6 | — | finished library | whether a last round swaps the stuffer for a capping block |
| 8 | Read linkage | finished library | barcode combination → cargo table | carry it forward, or rebuild a round |
| 9 | Read representation | finished library | combinations seen, and how evenly | at least 99.5% seen, skew under 10 |

**Step 1.** A pool is one position's parts. Each part is prepped on its own, then mixed in equal
amounts of DNA. A pool is built once, and cut again each round it supplies.

**Step 2.** Destination BbsI + SrfI, donor BsaI + PmeI, in separate tubes.

**Step 3.** 20 ng of cut destination per 200 µL, with donor at a 1:1 molar ratio.

**Step 4.** Up to 100 ng of the cleaned ligation goes into the cells.

**Step 5.** The control is the cut destination taken through the ligation with no donor added.
Subtract its colonies from the count. Both plates grow while the pool does, so they cost no
extra day. How many colonies a round needs is set by the build, not here: it follows from how
well the work downstream has to see each member. The source plates nothing at all, and
`docs/research/synthesis-and-assembly.md` says why this page does.

**Step 8.** The amplicon spans the whole cargo, so linkage is a long read: a barcode and the
coding bases it names have to arrive in one molecule. The plan designs the pair.

**Step 9.** Representation is read again after every later bottleneck. The amplicon spans the
barcode block alone, so it is a short read and cheap to repeat; the plan designs that pair too.

### Final assembly

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 1 | Pick the working vector | working vector stock | — | which application |
| 2 | Pick the cargo | DMX library or chosen wells | — | which cargo, and how it is pooled |
| 3 | Assemble with the cargo enzyme | both | final vector | whether the cargo is released in its own tube first |
| 4 | Clean up and electroporate | reaction | the library | how many colonies |
| 5 | Read representation | the library | what survived the transfer | at least 99.5% seen, skew under 10: screen it, or transfer again |

**Step 3.** One pot where the cargo enzyme is also what releases cargo from the DMX backbone.
Otherwise two stages in one tube: release with BsaI and PmeI, heat-kill, then add the working
vector, the cargo enzyme and ligase. A Gibson link or a re-tailoring PCR does the same job at a
different price. This is the one reaction where the working vector meets DMX-derived material,
because the iGGA rounds all finish first.

Neither by-product needs a check.

## Reagents and equipment

Detail is in `docs/research/synthesis-and-assembly-materials.md`, by the process that uses it,
and `docs/research/synthesis-and-assembly-barcode-kit.md`.

### Equipment list

| Equipment | Needed for | Model | Where |
| --- | --- | --- | --- |
| Acoustic liquid handler | DMX barcoding; primer combination plates | Echo 525 | HTS core |
| Automated liquid handler | combination plates, reformatting | Biomek i7 | planned purchase |
| Nanopore sequencer | validation read-out, library reads | MinION Mk1D | to buy |
| Colony picker | picking into 384-well plates, validation route only | — | undecided |

### Plasmid list

| Plasmid | Role | Source |
| --- | --- | --- |
| pCR-Blunt II-TOPO | part carrier | Zero Blunt TOPO kit |
| DMX0001 | DMX vector parent | Addgene 247434 |
| DMX0002 | alternative parent, cassette on the opposite strand | Addgene 247435 |
| DMX barcode kit | 96 plasmids, four groups of 24; one barcode per group marks a well | Addgene; accession range in `docs/research/synthesis-and-assembly-barcode-kit.md` |
| ccdB cassette, four versions | N and C, no N, no C, neither | synthesised, held as parts |
| Working vector backbone | per application | the user's own |

### Primer list

| Class | What it is | Source |
| --- | --- | --- |
| Part retailoring primers | `5'-[spacer]-CGTCTC-N-[new overhang]-[annealing]-3'` | designed per part |
| Orthogonal set, roles P1 / P2 / P3 | 165 twenty-mers of a tested 185, split 96 inner P2 / 35 P1 / 34 P3 | Subramanian 2018, Supplementary Table 1 |
| DMX library pair | carries the BsmBI sites giving `AGGA`/`TTCC` | Qian 2026, Supplementary Table 4 |
| DMX1–DMX6 | three pairs for the barcoded amplicon | Qian 2026, Supplementary Table 4 |
| dmx0 / dmx7 | universal primers flanking the design | named by Qian 2026, sequences derived in `docs/research/synthesis-and-assembly-barcode-kit.md` |
| Library read pairs | representation and linkage, across the barcode block | designed per library |
| Index PCR plate | one pair per well, ordered as published | Long 2025, Tables S3 and S4 |

### Reagent list

| Kind | Items |
| --- | --- |
| Type IIS enzymes | BsaI-HFv2, BsmBI, BbsI, and the cargo enzyme chosen for the working vector (PaqCI needs its activator oligo) |
| Blunt cutters | SrfI, PmeI |
| Ligases | T4 DNA ligase; T7 DNA ligase with StickTogether buffer |
| Clean-up | SPRI beads; a column clean-up kit; a plasmid prep kit |
| Cloning kit | Zero Blunt TOPO |
| Strains | DB3.1 or ccdB Survival 2 T1R (growing a stock whose ccdB is expressed); NEB Stable (a DMX stock, whose ccdB is T7-silent); a ccdB-sensitive strain carrying T7 polymerase (receiving cargo); electrocompetent BL21(DE3) (libraries); Endura (iGGA rounds) |
| Media and selection | low-salt LB; carbenicillin, kanamycin; glycerol for the archive |

**NEB Stable and Endura are ccdB-sensitive.** Neither published genotype carries a resistance
allele. That sensitivity is the point: it is why the ccdB cassette kills every vector that keeps
it after assembly. So a stock whose ccdB is expressed has to be grown in DB3.1 or ccdB Survival
2 T1R instead. Neither vendor names the allele those two are sold on.

## References

| Source | What it gives |
| --- | --- |
| Qian, Z. et al. Accelerating protein design by scaling experimental characterization. *Nat. Commun.* (2026). [doi:10.1038/s41467-026-76740-9](https://doi.org/10.1038/s41467-026-76740-9) | The DMX vector, the barcode kit and the barcoded ONT read-out |
| Lund, S., Potapov, V., Johnson, S. R., Buss, J. & Tanner, N. A. Highly parallelized construction of DNA from low-cost oligonucleotide mixtures using Data-optimized Assembly Design and Golden Gate. *ACS Synth. Biol.* 13, 745–751 (2024). [doi:10.1021/acssynbio.3c00694](https://doi.org/10.1021/acssynbio.3c00694) | Building a gene from a cheap oligo pool: the fragment split, its overhangs chosen on ligation fidelity, and fixed terminal overhangs |
| Subramanian, S. K., Russ, W. P. & Ranganathan, R. A set of experimentally validated, mutually orthogonal primers for combinatorially specifying genetic components. *Synth. Biol.* 3, ysx008 (2018). [doi:10.1093/synbio/ysx008](https://doi.org/10.1093/synbio/ysx008) | The orthogonal primer set the P1 / P2 / P3 roles draw from |
| Baker lab three-primer scheme | The nested PCR that demultiplexes the pool; unpublished, a diagram sent to us |
| Takacsi-Nagy, O. et al. Synthetic transcription factors designed by domain recombination enhance CAR T cell antitumor function. *Cell* 189, 1–20 (2026). [doi:10.1016/j.cell.2026.07.054](https://doi.org/10.1016/j.cell.2026.07.054) | Assembly in rounds: the round structure, the idea of an internal stuffer, the 11 bp barcode length, and the cells and temperature each round |
| Long, Y., Mora, A., Li, F.-Z., Gürsoy, E., Johnston, K. E. & Arnold, F. H. LevSeq: rapid generation of sequence-function data for directed evolution and machine learning. *ACS Synth. Biol.* (2025). [doi:10.1021/acssynbio.4c00625](https://doi.org/10.1021/acssynbio.4c00625) | The index PCR read-out, as the alternative to DMX barcoding |
| Twist Bioscience oligo pool and gene products | Pool lengths, sizes and prices |

How our design departs from Takacsi-Nagy is in
`docs/research/synthesis-and-assembly-departures.md`.

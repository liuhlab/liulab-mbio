# Highly Parallel DNA Synthesis and Assembly

## Goal

1. Large scale (10³–10⁴), long fragment (up to 2.5 kb) synthesis using the cheapest oligo pool (300 bp, 18,000 oligos for $12,500, Twist), ~$2.5 per kb (or lower) plus reagents and consumables.
2. Working vector system turns cargo fragments into different applications via a single Golden Gate Assembly (GGA).
3. Iterative GGA (iGGA) further enables deterministic or random assembly of fragment cargos with barcodes.

## Vectors and parts

What each molecule is. What is done with it is in `## Experiment Sections`.

### Vectors and parts summary

| Item                | Source and purpose                                                                    | Selection  | Common overhangs                             | Sites to remove outside designed cassettes | Cut by                                                                                                                         |
| ------------------- | ------------------------------------------------------------------------------------- | ---------- | -------------------------------------------- | ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------ |
| Part carrier vector | pCR-Blunt II-TOPO; stores each part                                                   | KanR       | -                                            | BsmBI                                      | -                                                                                                                              |
| Parts               | reusable in-frame elements: tag, linker, signal peptide, localization signal, degron  | -          | N-term pair TATG/AGGA; C-term pair TTCC/CTAA | BsaI, BsmBI                                | BsmBI (releases the part in the assembly reaction)                                                                             |
| DMX vector          | holds cargo                                                                           | KanR       | AGGA, TTCC                                   | BsaI, BsmBI, BbsI                          | BsmBI (insert cargo), BsaI (barcode, release cargo, iGGA donor), PmeI (blunt the donor backbone)                               |
| DMX barcode vector  | 96 barcodes in four groups of 24; one from each group marks a well for ONT sequencing | AmpR       | chain AGGA→GTTC→CCTT→TCAG→TTCC, group 1 to 4 | - (used as supplied)                       | BsaI (plate barcoding)                                                                                                         |
| Cargo               | Twist oligo pool fragments, via DAD and DMX                                           | -          | AGGA, TTCC                                   | BsaI, BsmBI                                | BsmBI (insert cargo to vector)                                                                                                 |
| iGGA cargo          | cargo built for iterative random assembly                                             | -          | AGGA, TTCC                                   | BsaI, BsmBI, BbsI, SrfI, PmeI              | BsmBI (into DMX vector), BsaI (release as donor), BbsI (open internal stuffer), SrfI (cut the released stuffer)                |
| Working vector      | chosen backbone plus working cassette; gives cargo its application                    | AmpR/CarbR | TATG, CTAA (cassette); AGGA, TTCC (cargo)    | BsaI, BsmBI                                | BsmBI (insert working cassette), BsaI (insert cargo)                                                                           |
| Final vector        | working vector with cargo installed                                                   | AmpR/CarbR | AGGA, TTCC (cargo junction)                  | -                                          | -                                                                                                                              |

The barcode chain is measured from the four plasmids we hold;
`docs/research/synthesis-and-assembly-barcode-kit.md` has all 96.

### Parts

#### Form and carrier

- Synthesis or proofreading PCR, flanked by inward-facing BsmBI sites and a default overhang: `[BsmBI.<5' overhang>]─part─[<3' overhang>.BsmBI]`
- Three cassette junctions: `TATG`/`AGGA` N-term part, `AGGA`/`TTCC` ccdB cassette and cargo, `TTCC`/`CTAA` C-term part
- Each part holds one pair, either the N-term or the C-term
- No part carries a BsaI or BsmBI site outside its own designed cassette; the carrier backbone carries no BsmBI site (e.g., pCR-Blunt II-TOPO)
- Carrier: pCR-Blunt II-TOPO (Zero Blunt TOPO kit, KanR + Zeocin, ccdB)
- A part has two usable forms: the carrier plasmid, and a PCR product carrying a different BsmBI flank, `5'-[spacer]-CGTCTC-N-[new overhang]-[annealing]-3'`

#### Rules a part obeys

- In frame, a multiple of 3 bp long, counted including its 4 bp overhang at 5'
- No stop codon
- Every overhang reads as one base closing a codon plus one codon: `T|ATG`, `A|GGA`, `T|TCC`, `C|TAA`
- The junction adds two amino acids

### Working vectors

#### Input working vector

- A backbone with no BsaI or BsmBI site, AmpR/CarbR, carrying `[TATG.BsmBI]─[stuffer]─[BsmBI.CTAA]` at the cassette site
- Stuffer: RFP cassette, long enough that the cut backbone separates on a gel
- The `ATG` of `TATG` serves as the start codon where the backbone carries a promoter and RBS before it; otherwise the construct's start sits in an N-terminal part

#### Working cassette

- `[TATG]─[N-terminal parts]─[ccdB cassette]─[C-terminal parts]─[CTAA]`
- Either side is optional, and the ccdB cassette comes in a version to match each case
- ccdB cassette (for later cargo insertion), supplied as a part in four versions:
  - N and C parts: `[BsmBI.AGGA.BsaI]─[promoter+ccdB]─[BsaI.TTCC.BsmBI]`
  - no N part: `[BsmBI.TATG]─GG─[AGGA.BsaI]─[promoter+ccdB]─[BsaI.TTCC.BsmBI]`
  - no C part: `[BsmBI.AGGA.BsaI]─[promoter+ccdB]─[BsaI.TTCC]─GG─[CTAA.BsmBI]`
  - neither: `[BsmBI.TATG]─GG─[AGGA.BsaI]─[promoter+ccdB]─[BsaI.TTCC]─GG─[CTAA.BsmBI]`
  - All four expose `AGGA`/`TTCC` for the cargo; a skipped side adds one Gly

### DMX vectors

- The modified DMX vector has no BsaI, BsmBI, BbsI or SrfI outside the ccdB cassette, carries one PmeI site outboard of each BsaI site, uses KanR
- The PmeI site blunts the donor backbone in iGGA; the DMX pipeline does not use it
- DMX ccdB cassette: `[PmeI]─[BsaI.AGGA.BsmBI]─[promoter+ccdB]─[BsmBI.TTCC.BsaI]─[PmeI]`
- BsmBI inserts cargo; BsaI adds plate barcodes or releases cargo
- The ccdB cassette reads from a T7 promoter under a *lac* operator, with the repressor on the
  plasmid. Only a strain carrying T7 polymerase counter-selects

### Final vectors

- A working vector whose ccdB cassette has been replaced by cargo, joined at `AGGA`/`TTCC`
- AmpR/CarbR, where the DMX backbone it was released from is KanR

### iGGA cargo

Runs in the DMX vector; the working vector takes only the finished cargo.

- iGGA cargo: `[AGGA]─[fragment]─[internal stuffer]─[11 bp barcode]─[TTCC]`
- Each fragment, counted from the first base of its `AGGA`, is a multiple of 3
- Internal stuffer `[AGGA.BbsI]─[SrfI]─[BbsI.TTCC]`, 34 bp: `AGGAAAGTCTTCAGCCCGGGCAGAAGACAATTCC`
- Barcode: 11 bp, minimum edit distance 3 over the whole set, no in-frame stop codon, GC 3 to 7 of 11, no homopolymer over 3
- The paper's 72 barcodes do not satisfy those rules here: six carry a stop codon at our frame offset
- The stuffer is retained every round, so the chain is a cycle and the product stays releasable
- A capping block is any `AGGA`/`TTCC` block that is 1 mod 3 counted from its `AGGA`, such as T2A. It ends the chain, because the product then carries no stuffer to reopen

## Experiment Sections

Major steps, the decisions attached to them, and anything designed that must be followed
exactly. Volumes, temperatures, times, catalogue numbers and anything a person follows
without deciding belong to the detail protocol.

Every plasmid is confirmed by whole-plasmid sequencing before it becomes a stock. Only a
pooled screen skips it.

### Lab resources (build once)

Built once and reused by every project afterwards.

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 1 | Clone each part into the carrier | BsmBI-flanked part | part plasmid | which overhang pair the part gets, by the side it is more likely to sit on |
| 2 | Rebuild the DMX vector backbone | DMX0001 or DMX0002 | KanR DMX vector, BbsI-free, PmeI added | which parent |
| 3 | Build the input working vector | chosen backbone | backbone carrying the stuffer cassette | the backbone; where the cassette goes; whether promoter and RBS sit in the backbone or in an N-terminal part |
| 4 | Whole-plasmid sequencing | candidates from 1–3 | verified stocks | which clone becomes the stock |
| 5 | Lay out the primer plates | orthogonal primer set | source plates, then pre-made combination plates | which primer fills which slot |
| 6 | Stock the barcode kit | DMX barcode plasmids | 96 plasmids in four groups | — used as supplied |

**Step 2** replaces AmpR with KanR, removes the four BbsI sites (two in `lacI`, two in the
backbone), adds one PmeI site outboard of each BsaI site, and confirms nothing else cuts.

**Step 5.** A slot always means the same primer pair, so plate maps are constants rather than
per-project files. A gene whose own sequence contains its slot's primer site moves slot; the
assignment is not redrawn. Combination plates are kept pre-made and a plate is drawn when
needed.

Anything carrying ccdB is propagated in NEB Stable, which carries no T7 polymerase and so never
transcribes the cassette.

### Working cassette (per application)

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 1 | Choose the cassette's parts | part collection | N-terminal and C-terminal parts, or neither | what each side carries, and whether either is skipped |
| 2 | Take the matching ccdB cassette | four versions | the one fitting which sides are used | follows from step 1 |
| 3 | Retailor an overhang | part plasmid | BsmBI-tailed PCR product | only when a part's default pair is wrong for this cassette |
| 4 | Assemble | input working vector, parts | working vector | — one BsmBI reaction cuts every part |
| 5 | Screen and pick | transformants | candidate clones | — colonies that lost the stuffer's fluorescence |
| 6 | Whole-plasmid sequencing | candidate clones | working vector stock | which clone becomes the stock |

Check total length and part compatibility before assembling, and check that the two amino
acids each junction adds do no harm in this construct. A part that fits its default pair is
added as the carrier plasmid, with no digest and no PCR of its own.

### Cargo synthesis

Builds full-length cargo from an oligo pool and puts it in the DMX vector. From PCR2 onward,
identity is held by well position.

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 1 | Collect input sequences | protein or DNA sequences | the design list | which genes, and how long |
| 2 | Check compatibility | design list | cleared list | recode or drop any gene carrying a restriction site or a primer site |
| 3 | DAD fragment split | cleared list | fragments, internal overhangs | oligo product and length; `AGGA`/`TTCC` are held out of the set DAD may choose |
| 4 | Design the barcode set — iGGA cargo only | part lists | one barcode per fragment | — drawn fresh to the rules in `### iGGA cargo` |
| 5 | Assign primer slots | fragments | oligo layout | batch size |
| 6 | Order the pool | oligo layout | oligo pool | pool size, a budget choice rather than a cap |
| 7 | PCR1, per batch | pool | one batch's pieces | — |
| 8 | PCR2, per gene | batch | one gene's pieces | — |
| 9 | GGA into the DMX vector | pieces, vector | cargo in vector | — |
| 10 | Transform and archive as cells | GGA product | glycerol archive, one design per well | polyclonal, or go on to validation |

**Steps 3 and 5, the three-primer scheme.** Every oligo is `[P1][fragment][P2][P3]`. PCR1 uses
the pair (P1, P3) shared by a whole batch of genes and pulls that batch out of the pool. PCR2
uses (P1, P2) to pull one gene out of the batch and drops the P3 region. Because groups are
selected by a *pair*, capacity is combinatorial and is not the limiting factor.

**Step 5, batch size.** PCR2 is one reaction per gene whatever we choose, and that is the real
ceiling on scale. The only free choice is how many genes PCR1 pulls at once. What crowds a
PCR1 tube is pieces, not genes — every piece in the batch competes for the same primer pair,
and a gene fails if any one of its pieces drops out. So hold pieces per PCR1 tube roughly
constant and let batch size fall as genes lengthen, rather than fixing a gene count.

**Step 10.** Polyclonal by default, for a pooled screen. Validate as fragment count rises.

**The archive.** Cells, not DNA: duplicate plates, drawn one way, never re-frozen. DNA is
prepared only when a set is drawn, one prep for the set. No per-well normalisation.

### Cargo validation (optional)

The cargo sits in the DMX vector either way; this section only reads each well. Identity stays
with well position throughout.

#### Getting to clonal wells

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 1 | Spot the archive as an array | archive plate | one spot per design | wells per tray |
| 2 | Grow | tray | colonies per spot | — |
| 3 | Pick into 384-well, on the QPix | tray | clonal cultures, position kept | colonies per design |
| 4 | Grow to confluence | 384-well plates | cultures ready to lyse | — |

#### Route A — DMX barcoding

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 5 | Compress four 384-well plates into one 1536-well, on the Echo | cultures | cell lysate per well | — |
| 6 | Add one barcode from each of the four groups, on the Echo | barcode kit | a combination per well | which combination marks which well |
| 7 | Barcoding GGA, in lysate | lysate, barcodes | barcoded construct per well | — |
| 8 | Pool, then PCR the three primer pairs separately | barcoded pool | three amplicon pools | — |
| 9 | Sequence | amplicons | reads | — |
| 10 | Basecall and demultiplex | reads | a sequence per well | what counts as a pass |
| 11 | Reformat | per-well calls | compacted plate | which wells carry forward |

#### Route B — index PCR

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 5 | Colony PCR, one barcoded primer pair per well | cultures | barcoded amplicon per well | which pair marks which well |
| 6 | Pool per plate and clean up | amplicons | one pool per plate | — |
| 7 | Sequence | pool | reads | — |
| 8 | Demultiplex | reads | a sequence per well | what counts as a pass |
| 9 | Reformat | per-well calls | compacted plate | which wells carry forward |

#### After either route

Choose the route per project: at or below one plate of samples, index PCR; above it, DMX
barcoding.

Reformatting only compacts out the wells that failed; every well already has its identity. The
polyclonal route skips this section, and its plate is already a clean array.

### iGGA pipeline

Runs in the DMX vector. iGGA cargo is one kind of cargo and may be polyclonal or validated like
any other.

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 1 | Choose the pools and their order | one part list per round | a round plan | which list goes first; order sets position in the construct |
| 2 | Split digest | destination pool, donor pool | opened destination, released donor | — |
| 3 | Ligate | cut ends | circular product | — |
| 4 | Electroporate, recover and grow at 30 °C | ligation | the round's pool | — |
| 5 | Prep the pool | the round's pool | next round's destination | — |
| 6 | Repeat 2–5 | | finished library | whether a last round swaps the stuffer for a capping block |
| 7 | Linkage read | finished library | barcode combination → cargo table | carry it forward, or rebuild a round |
| 8 | Representation read | finished library | combinations seen, and how evenly | — |

**Step 2.** Destination and donor are cut in separate tubes: destination BbsI + SrfI, donor
BsaI + PmeI. The two digests never share a tube.

**Rounds run unchecked**, as the paper runs them. Nothing is plated and nothing is read between
rounds. The chemistry removes its own failures — a destination that escapes BbsI is still cut
blunt by SrfI, and the ligase does not join blunt — and the linkage read names which position
lacked its barcode, so a failed round is localised afterwards rather than missed. Rounds are
pooled and cheap, so paying a few of them to learn that is a better trade than a check at each.

**Step 7, the linkage read.** Each fragment's barcode is attached at synthesis, so the mapping
is designed rather than discovered; the long read verifies it survived synthesis and assembly.
Without it the barcode is a label with nothing behind it, because the screen returns barcodes
only. It runs here because linkage cannot change afterwards — cargo moves into the working
vector as one intact piece.

**Step 8, representation.** Measured again after every later bottleneck, since each
electroporation resamples the pool.

### Final assembly

| # | Step | In | Out | Decision |
| --- | --- | --- | --- | --- |
| 1 | Pick the working vector | working vector stock | — | which application |
| 2 | Pick the cargo | DMX library or chosen wells | — | which cargo, and how it is pooled |
| 3 | One-pot BsaI | both | final vector | — |
| 4 | Clean up and electroporate | reaction | the library | how many colonies, from the coverage wanted |
| 5 | Representation read | the library | what survived the transfer | screen it, or transfer again |

Neither by-product needs a check: the released DMX backbone is KanR on an Amp plate, and
released ccdB kills anything that takes it.

## Reagents and equipment

Detail is in `docs/research/synthesis-and-assembly-materials.md`, by the process that uses it,
and `docs/research/synthesis-and-assembly-barcode-kit.md`.

### Equipment list

| Equipment | Needed for | Where |
| --- | --- | --- |
| Acoustic liquid handler | DMX barcoding; primer combination plates | **Echo 525**, HTS core |
| Automated liquid handler | combination plates, reformatting | **Biomek i7**, planned purchase |
| Nanopore sequencer | validation read-out, library reads | **MinION Mk1D**, to buy |
| Colony picker | picking into 384-well, validation route only | **QPix 420 or FLEX**, undecided |

### Plasmid list

| Plasmid | Role | Source |
| --- | --- | --- |
| pCR-Blunt II-TOPO | part carrier | Zero Blunt TOPO kit |
| DMX0001 | DMX vector parent | Addgene 247434 |
| DMX0002 | alternative parent, cassette on the opposite strand | Addgene 247435 |
| DMX barcode kit, 96 plasmids | one barcode per group marks a well | Addgene; accession range in `docs/research/synthesis-and-assembly-barcode-kit.md` |
| ccdB cassette, four versions | N and C, no N, no C, neither | synthesised, held as parts |
| Working vector backbone | per application | the user's own |

### Primer list

| Class | What it is | Source |
| --- | --- | --- |
| Part retailoring primers | `5'-[spacer]-CGTCTC-N-[new overhang]-[annealing]-3'` | designed per part, as needed |
| Orthogonal set, roles P1 / P2 / P3 | 185 twenty-mers, 165 of them in the orthogonal set | Subramanian 2018 |
| DMX library pair | carries the BsmBI sites giving `AGGA`/`TTCC` | Qian, Supplementary Table 4 |
| DMX1–DMX6 | three pairs for the barcoded amplicon | Qian, Supplementary Table 4 |
| dmx0 / dmx7 | universal primers flanking the design | Qian |
| Index-PCR barcoded plate | one pair per well, if that read-out is chosen | Long et al., Tables S3 and S4; ordered as published |

### Reagent list

| Kind | Items |
| --- | --- |
| Type IIS enzymes | BsaI-HFv2, BsmBI, BbsI |
| Blunt cutters | SrfI, PmeI |
| Ligases | T4 DNA ligase; T7 DNA ligase with StickTogether buffer |
| Clean-up | SPRI beads; a column clean-up kit; a plasmid prep kit |
| Cloning kit | Zero Blunt TOPO |
| Strains | NEB Stable (anything with ccdB); a ccdB-sensitive strain carrying T7 polymerase (receiving cargo); electrocompetent BL21(DE3) (libraries); Endura (iGGA rounds) |
| Media and selection | low-salt LB; carbenicillin, kanamycin; glycerol for the archive |

## References

| Source | What it gives |
| --- | --- |
| Qian, Z. et al. Accelerating protein design by scaling experimental characterization. *Nat. Commun.* (2026). [doi:10.1038/s41467-026-76740-5](https://doi.org/10.1038/s41467-026-76740-5) | The DMX vector, the barcode kit and the barcoded ONT read-out |
| Lund, S., Potapov, V., Johnson, S. R., Buss, J. & Tanner, N. A. Highly parallelized construction of DNA from low-cost oligonucleotide mixtures using Data-optimized Assembly Design and Golden Gate. *ACS Synth. Biol.* 13, 745–751 (2024). [doi:10.1021/acssynbio.3c00694](https://doi.org/10.1021/acssynbio.3c00694) | Building a gene from a cheap oligo pool: the fragment split, its overhangs chosen on ligation fidelity, and fixed terminal overhangs |
| Subramanian, S. K., Russ, W. P. & Ranganathan, R. A set of experimentally validated, mutually orthogonal primers for combinatorially specifying genetic components. *Synth. Biol.* 3, ysx008 (2018). [doi:10.1093/synbio/ysx008](https://doi.org/10.1093/synbio/ysx008) | The orthogonal primer set the P1 / P2 / P3 roles draw from |
| Baker lab three-primer scheme | The nested PCR that demultiplexes the pool; unpublished, a diagram sent to us |
| Takacsi-Nagy, O. et al. Synthetic transcription factors designed by domain recombination enhance CAR T cell antitumor function. *Cell* 189, 1–20 (2026). [doi:10.1016/j.cell.2026.07.054](https://doi.org/10.1016/j.cell.2026.07.054) | Assembly in rounds: the round structure, the idea of an internal stuffer, the 11 bp barcode length, and the cells and temperature each round |
| Long, Y., Mora, A., Li, F.-Z., Gürsoy, E., Johnston, K. E. & Arnold, F. H. LevSeq: rapid generation of sequence-function data for directed evolution and machine learning. *ACS Synth. Biol.* (2025). [doi:10.1021/acssynbio.4c00625](https://doi.org/10.1021/acssynbio.4c00625) | The index-PCR read-out, as the alternative to DMX barcoding |
| Twist Bioscience oligo pool and gene products | Pool lengths, sizes and prices |

How our design departs from Takacsi-Nagy is in
`docs/research/synthesis-and-assembly-departures.md`.

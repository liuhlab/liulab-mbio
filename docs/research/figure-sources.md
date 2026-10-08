---
search:
  exclude: true
---

# Figure sources for four designs

Four designs are hard to read as prose: the Baker lab's three-primer PCR, the Golden Gate that
puts a cargo gene into its vector, DMX, and iGGA. Each already has a published picture. We draw
our own, because a drawing of ours takes parameters, crops to one step and lights one reaction;
a published image shows what it showed for its own paper and nothing in it turns on or off.

This note records what each published figure holds, closely enough to draw an equivalent
without it, what a drawing of ours would have to vary, and whether the package already holds
the data to draw it. It decides nothing.

Everything cited is a file under `reference_docs/`, which a fresh clone does not have. Each
source's own citation is in the Sources section at the end.

## Summary

| Design | The figure | Can the package draw it today? |
| --- | --- | --- |
| Three-primer PCR | `long_fragment_GGA/Baker Lab Strategy.png`, an unnumbered diagram | Yes |
| Cargo GGA ligation | Qian 2026 Supplementary Fig. 4, with Fig. 4b and Lund 2024 Figure 1 beside it | Yes |
| DMX | Qian 2026 Fig. 4a–d | In part: panels b and the plate steps, not the workflow track |
| iGGA | Takacsi-Nagy 2026 Figure S1A | Yes, bar the between-row layout |

All four are covered by what is already downloaded. Nothing needs fetching. One caveat is in
**Where a source is thin**.

## 1. Baker's three-primer PCR design

**File.** `reference_docs/synthesis_and_assembly/long_fragment_GGA/Baker Lab Strategy.png`. It is
a standalone diagram, not a numbered figure: that directory's README records that it came from
the Baker lab directly and that no paper documents the scheme.

### What it shows

Three rows, each a caption on the left and a drawing on the right.

**Row 1, PCR1.** Caption: "All fragments for 96 genes per well / hundreds oligos per well". One
oligo is drawn as a horizontal bar of four abutting boxes, left to right: `Primer 1` (blue),
`Gene 1 fragment 1` (pale yellow), `Primer 2` (mauve), `Primer 3` (grey). Two red arrows mark the
pair in use: a forward arrow above the bar, starting over `Primer 1` and pointing right, and a
reverse arrow below the bar, starting at the right of `Primer 3` and pointing left. PCR1 is
therefore the outer pair, `Primer 1` with `Primer 3`, and its product still carries the `Primer 2`
site.

**Row 2, PCR2.** Caption: "All fragments for single gene per well / ~1-8 oligos per well". The
same four-box bar, and the same forward arrow over `Primer 1`, but the reverse arrow now sits
under `Primer 2`. The pair is nested inside row 1's product. `Primer 3` is drawn but unused.

**Row 3, Golden Gate Assembly.** Caption: "All fragments for single gene + backbone per well".
Three stacked bars, each `Primer 1` | `Gene 1 fragment N` | `Primer 2` for N = 1, 2, 3 — the
`Primer 3` box is gone, consumed by PCR2. A thick red vertical tick sits at each boundary between
a primer box and the fragment: the Type IIS cuts, inboard of both primer sites. Below the three
bars, a grey ring labelled `Backbone DNA` carries two red ticks on its upper arc. A long black
arrow runs right, labelled "Type IIS restriction enzyme + T4 DNA ligase", to a second grey ring
labelled `Assembled DNA` whose top arc holds three pale-yellow boxes seated side by side — the
three fragments joined in order.

### What a reader takes away

One oligo is read by two nested primer pairs. The outer pair pulls a whole 96-gene subpool out of
the chip; the inner pair, sharing the same forward primer, pulls one gene out of that subpool.
Because both Type IIS cuts lie inboard of the primer sites, the primer sites leave with the cut
and the fragments join without them. The well's complexity falls from hundreds of oligos to one
to eight between the two PCRs, which is the whole point of the second reaction.

The `long_fragment_GGA` README states the layout in text as `[P1][fragment][P2][P3]` and the cost:
N + 1 reactions per 96 genes, bought in exchange for a primer inventory that stops growing with
the library. `dmx/jason_email.md` confirms P1, P2 and P3 are drawn from the same orthogonal primer
set Lund 2024 uses, and that the lab runs 300 nt pools with one to three fragments per gene.

### What ours would have to vary

- **Fragment count.** The drawing fixes three. A real design has as many as the split chose.
- **Oligos per well at each stage.** "hundreds" and "~1-8" are this library's numbers; ours are
  the batch size and the per-block split, both computed.
- **Which primer fills which role, by name.** The three roles are constant; the sequences are
  allotted from the orthogonal set per build.
- **The enzyme and the overhangs.** The drawing says only "Type IIS restriction enzyme" and marks
  the cuts as ticks. Ours knows the enzyme and every junction's four bases.
- **The filler.** The diagram omits the stuffer that pads an oligo to the pool's one length.
  Lund 2024 Figure 2A (below) is where that appears.
- **Which PCR a step is showing.** A protocol step for PCR1 wants PCR2 dimmed, and the reverse.

Lund 2024 Figure 2A, on page 2 of the article PDF, is the companion worth knowing about: four
oligo layouts drawn as coloured bars, each reading `Barcode 1` | `BsaI` | `Fragment sequence` |
`BsaI` | and then `Barcode 2` and `Stuffer` in the two possible orders, labelled "200 nt outside",
"200 nt full", "300 nt outside" and "300 nt full". `Barcode 2` is printed upside down to show it
lies on the other strand. It is the layout question our oligo drawing answers by construction.

### Can the package draw it

Yes. `liulab_synbio.igga.cargo` allots the three roles — it names them `inner`, `forward` and
`outer` — groups blocks into batches each holding a forward and an outer primer, and exposes
`inner_pairs()` as "the batch forward, the block inner", which is exactly row 2's pair.
`liulab_mbio.bench.pools` builds each oligo as a record carrying its primer sites, its recognition
sites and its filler. `liulab_mbio.plot.linear` draws a record with its features as named arrows
and each primer as a thin arrow at its binding site. Both rows and both arrow pairs are package
data already.

## 2. Cargo GGA ligation

Three figures cover this at three zoom levels. The base-level one is the most useful.

**File.** `reference_docs/synthesis_and_assembly/dmx/paper/41467_2026_76740_MOESM1_ESM.pdf`,
**Supplementary Fig. 4**, "Multi-use cloning strategy for DMX vector", on page 20 of that PDF.

### What it shows

Two base-level panels, one above the other, with a black down arrow between them.

**Top, the empty vector.** About sixty bases of the DMX vector drawn as two rows of monospace
letters, top strand over bottom strand. Five call-out boxes label the stretch, left to right:
`BsaI` (pale yellow), `BsmBI` (blue), `ccdB gene` (black box, white letters), `BsmBI` (blue),
`BsaI` (pale yellow). The top strand reads `ATGG` `GGTCTC` `A` `AGGAG` `GAGACG`, then the ccdB
body, then `CGTCTC` `G` `TTCCG` `GAGACC`, with the complement beneath it; the recognition
hexamers are shaded in their box's colour. The ccdB body is not spelled out — it is a black ladder
of rungs between the two strands, the conventional "and this continues" glyph. Thin vertical bars
cut through both strands at the BsaI cut positions, outboard of the BsmBI pair. Under the left end
of the duplex a row of pentagon chevrons reads `M` `G` `S` `Q` `G`: the vector's own N-terminal
codons, in frame.

**Bottom, after the gene goes in.** The same flanks, with the black ccdB ladder replaced by a
green one inside a green box labelled `GOI`. The two `BsmBI` boxes are gone, consumed by the cut;
only the two `BsaI` boxes remain. The left chevrons still read `M G S Q G`; a new chevron row at
the right reads `G` `S` `D` `H` `W`, then a `6x His` box and a yellow star for the stop.

### What a reader takes away

The vector carries a nested cassette: BsaI outside, BsmBI inside, ccdB between. BsmBI opens it for
pool entry, on `AGGA` and `TTCC`. The gene lands in frame between a fixed N-terminal `MGSQG` and a
C-terminal `GSDHW`, a His tag and a stop. BsaI survives into the product, which is what later
barcodes it or swaps it into another vector. ccdB kills anything that keeps the cassette, so an
empty vector never reaches a colony. The caption says this in one line: BsmBI clones the library
in; BsaI barcodes, or swaps the gene out on `AGGA`/`TTCC`.

### The same step, zoomed out

**Qian 2026 Fig. 4b**, page 6 of the article PDF, is the plasmid-level version. `Gene library` is a
stack of offset coloured fragments; `DMX vector` is a ring carrying a skull for ccdB and two short
coloured blocks at the cassette. A square bracket joins the two down to a stack of rings, each ring
now carrying the fragment, with the coloured end blocks still on it. From that stack, one arrow
goes down to a bar reading `Design` | `His` labelled "Expression", and a second goes right to the
barcoding branch and to a bar pair `Design`|`GFP` and `His`|`Design` labelled "1-pot transfer to
other vectors".

**Lund 2024 Figure 1, step 6**, page 2 of that article's PDF, is the same reaction without the
vector detail. Three columns, one per gene, each a different colour and a different fragment
count. In each, the amplified barcoded fragments (coloured bars with a magenta block at one end
and a green or amber one at the other) collapse by a down arrow into one continuous multi-segment
bar flanked by two black bars standing for the vector. Steps 1 to 5 above it are "Codon optimize",
"Use SplitSet", "Append barcodes", "Order all as pool" — drawn as a tangle of coloured squiggles,
with arrows converging into it and diverging out — and "Amplify assemblies".

### What ours would have to vary

- **Which enzyme does which job.** Our method searches a cargo enzyme per vector; the figure's
  BsmBI and BsaI are DMX's choices.
- **The two entry overhangs.** `AGGA`/`TTCC` is DMX0001's. The `dmx` README records that DMX0002
  carries the same cassette on the opposite strand, so its cut ends read `GGAA` and `TCCT` — the
  same two overhangs seen from the other side.
- **The flanking codons and the tag.** `MGSQG` and the 6×His are this vector's; ours are whatever
  the user's record spells there.
- **Fragment count, and so the internal junctions and their four bases.** The figure shows a
  single `GOI` ladder and no internal junction at all.
- **Which junction is lit.** A step that is checking one junction wants the others dimmed.
- **Whether a counter-selection cassette is there.** A working vector's cassette is the user's.
- **The bases themselves.** This is the whole advantage: the published panel can only ever show
  Qian's sequence, while ours is drawn over the record in hand.

### Can the package draw it

Yes, and more directly than any of the other three. `liulab_mbio.plot.sequence_view` already draws
this exact view: a position ruler, the top strand, a rail, the bottom strand, each cut site drawn
through the top strand and through the bottom strand where each enzyme cuts it so the overhang
shows, each feature as a bar under the bases, and a CDS's translation above its bar with each
amino acid centred on its codon and a stop in red. That is Supplementary Fig. 4's whole content.
The DMX vector record itself is in the repo: `scripts/build_dmx_vector.py` rebuilds the destination
from `tests/data/dmx0001.gb`, and the `dmx` README's site table was counted with
`liulab_mbio.sites` rather than read off a map.

## 3. DMX

**File.** `reference_docs/synthesis_and_assembly/dmx/paper/Qian et al., Accelerating protein design
by scaling experimental characterization.pdf`, **Fig. 4**, "DMX pipeline", page 6. Four panels.

### Panel a, the pipeline

Ten numbered stops on a serpentine track: five along the top running left to right, five along the
bottom running right to left, with a timeline ribbon between the two rows.

1. **Gene library amplification.** An oligo chip drawn as a grid of coloured spots, an arrow to a
   PCR tube.
2. **GGA into dual-use DMX vector.** A tube.
3. **Colony picking.** A colony-picker instrument, an arrow to an agar plate with colonies, then
   to a stack labelled `Culture plates`.
4. **Plate duplication and compression.** `Culture plates` marked `4x` compressing to `1x`
   `Reaction plates`, by ECHO.
5. **Cell lysis.** A plate with a thermometer, `95 °C`, `20 min`.
6. **Barcode transfer.** A plate whose wells carry four coloured barcode labels, transferring into
   the reaction plate.
7. **Isothermal GGA in cell-lysate.** A plate with `37 °C`, `1 h`.
8. **Pool barcoded genes.** Several rings drawn into one tube.
9. **Sequence library.** A nanopore sequencer with a trace.
10. **Sequence-verified ready-to-use genes.** A plate grid of blue squares with a few orange ones,
    under a legend reading `Sequence:` `correct` / `incorrect`, and an arrow back out to new plates.

The ribbon between the rows carries one arrow box per step, pink for hands-on time and black for
automated, each holding its duration: 2 h, 2 h, O/N then 1–3 h, O/N then 1–2 h, 20 min, 1–2 h, 1 h,
1 h, 1–3 days.

**One discrepancy worth knowing.** The caption says the lysis is "98 °C for 20 minutes"; the icon
in the panel reads 95 °C. Both are in the same figure.

### Panel b

The plasmid-level cargo GGA ligation, described in section 2 above.

### Panel c, the Sankey

Four stages left to right, the width of each ribbon its share. `All wells` (one solid blue bar) →
93% `Complete barcode`, 7% `Incomplete barcode` → 70% `Correct gene`, 30% `Incorrect gene` → 92%
`1 gene / well`, 8% `>1 gene / well`. N = 4608 wells.

### Panel d, the cost

Grouped bars of total cost in dollars against library size in number of designs, DMX on the left
of a vertical divider and eBlock on the right, each bar stacked into DNA cost, sequencing cost and
cloning cost. The eBlock bars run to roughly $50,000 at the largest library; the DMX bars stay low
across the range. It is a crossover argument drawn as two families of bars.

### What ours would have to vary

- **The plate format and the compression ratio.** 4-to-1 into 384 is this lab's; `liulab_mbio.bench`
  holds the format parameter and the moves between wells.
- **Which marking route the build chose.** The published panel only knows the DMX barcode route;
  a build may read its wells by index PCR instead, and then steps 6 and 7 are different steps.
- **The barcode groups.** Four groups of 24 is the kit's; the figure hard-codes `BC1`–`BC4`.
- **Every incubation and duration.** The ribbon's numbers are this protocol's; ours are computed,
  and a step's figure should agree with the step's own table.
- **Panel c's percentages.** They are a measurement of one run of 4608 wells. Ours would be this
  build's own counts, or nothing — a Sankey with no numbers behind it is decoration.
- **Panel d entirely.** Cost is a price record the user holds, so the bars are per-user, and the
  comparison arm is whatever the user is comparing against.
- **Which step is lit**, when the figure sits inside one step of a protocol rather than above all
  of them.

### Can the package draw it

In part.

- **Panel b**: yes, as section 2 says.
- **The plate steps of panel a** — 3, 4, 6 and 10: yes. `liulab_mbio.plot.plate` draws a plate's
  wells, and the protocol renderer already inlines one. Step 10's grid of correct and incorrect
  wells is a plate coloured by a per-well verdict, and `liulab_synbio.dmx` holds a well's derived
  address and the pass rule that decides the colour.
- **Panel c**: only once a run has been counted. The shape is a generic Sankey; nothing in
  `liulab_mbio.plot` draws one today.
- **Panel d**: the numbers come from `liulab_mbio.bench.prices`, which prices a bill from a tariff
  the user holds. The chart itself is not drawn anywhere today.
- **The serpentine track, the instrument icons and the timeline ribbon**: not package data at all.
  That is a drawing of a workflow, and the protocol model has no figure block to put one in.

## 4. iGGA

**File.** `reference_docs/synthesis_and_assembly/prot-assembly/Takacsi-Nagy et al. 2026 - mmc1.pdf`,
**Figure S1A**. The figure is on page 2 of that PDF; its caption, "Figure S1. Related to Figure 1,
Molecular and assay development for DESynR TF", runs on page 3.

### What it shows

One descending column of linear part diagrams. The left margin names the round; the right margin
names the library each row belongs to. A library of many members is drawn as three identical bars
offset behind one another inside a dotted ellipse; a single molecule is one bar.

**Row 1, `N Domain Plasmid Library`.** `5' External Stuffer` | `N Domain` (pink) | `Internal
Stuffer` (blue) | `BC` (pink) | `3' External Stuffer`. Two `BsaI` labels with hooked arrows sit
under the two edges of the internal stuffer — the internal digest that opens this plasmid as the
destination.

**Row 2, `bZIP Domain Plasmid Library`.** The same layout in green, with two `BbsI` hooked arrows
under its two external stuffers — the external digest that releases it as the donor. Two red
dashed lines run from the edges of row 1's internal stuffer down to the ends of row 2's donor: the
matching overhangs. The margin reads `Cloning Round 1`.

**Row 3, `C Domain Plasmid Library`.** The same again in orange, with an extra `2A` box between
`C Domain` and `Internal Stuffer`, and `BbsI` arrows under its external stuffers. Red dashed lines
again run down from row 2. Margin: `Cloning Round 2`. A thick black down arrow leaves it.

**Row 4, `Fully Assembled DESynR AP-1 TF Library`.** `5' External Stuffer` | `N Domain` | `bZIP
Domain` | `C Domain` | `2A` | `Internal Stuffer` | `BCs` | `3' External Stuffer`, with `BbsI` arrows
under both external stuffers. Two small call-outs under the bar spell two junctions as duplexes:
`—GGAG—` over `—CCTC—` at the N/bZIP junction, and `—CCGA—` over `—GGCT—` at the bZIP/C junction.
Margin: `Cloning Round 3`.

**Row 5, `TRAC Exon Knock-In Vector`.** The last round's destination: `5' Homology Arm` | `2A` |
`tNGFR` | `2A` | `CAR (CD19-28ζ)` | `2A` | `Internal Stuffer` | `BC` | `3' Homology Arm`, with two
`BsaI` arrows under its internal stuffer. Red dashed lines run from row 4 into it; a thick black
down arrow leaves it.

**Row 6, `Knock-In Cassette`.** The product: `5' Homology Arm` | `2A` | `tNGFR` | `2A` |
`CAR (CD19-28ζ)` | `2A` | `N Domain` | `bZIP Domain` | `C Domain` | `2A` | `Internal Stuffer` |
`BCs` | `3' Homology Arm`.

**The inset.** A shaded wedge blows row 4's `BCs` out to the right: a small bar reading `Internal
Stuffer` | `BC C` | `BC bZIP` | `BC N` | `BC Vec`, with `Fwd Seq Primer` arrowed in from above
between two `BsaI` hooks, and `Cloning Scar` labelling the last join. Captioned "Within Continuous
Reading Frame". The barcodes stack in reverse order of the rounds — the last part's barcode is
nearest the stuffer — which is what makes one read name the whole combination.

### What a reader takes away

One picture holds the method. Two enzymes with fixed jobs: one opens the destination at its
internal stuffer, the other releases the donor from its external stuffers. The internal stuffer
carries the next position's entry overhang, which is how a part knows its place. One barcode is
appended per round, and they collect into a single block that stays in frame and is read by one
primer. The same operation repeats until the library enters its final vector.

**Follow the figure, not the text.** The `prot-assembly` README records that METHOD DETAILS on
page e4 has the donor and the destination swapped, and that only the figure's assignment produces
matching overhangs when simulated on the published sequences.

### What ours would have to vary

- **The number of positions**, and so the number of rows and rounds. Three is this build's.
- **The enzyme names** — a scheme names its internal, external and blunt enzymes.
- **The entry overhangs.** `GGAG` and `CCGA` are spelled out in the figure for one build; a design
  chooses its own and has to show them.
- **The barcode length**, 11 bp here, and **how many barcodes have collected** by a given round.
  The inset is a function of the round number, not a fixed drawing.
- **The constants between a domain and its stuffer.** The C position carries a 2A here; another
  scheme need not.
- **Which round is lit.** A protocol step for round 2 wants rounds 1 and 3 dimmed — this is the
  clearest case in all four designs for a parameter the published image cannot offer.
- **The final destination.** A TRAC knock-in vector is this paper's; ours is the user's working
  vector.

### Can the package draw it

Yes, bar one piece of layout. `liulab_synbio.igga.rounds` simulates each round and annotates the
records; `write_records` writes one record per round, and `plan_igga` also writes a block vector
per position and the product. Each part block is built with its CDS, its stuffers and its barcode
as features. The scheme's enzymes and entry overhangs are the method's own data. So every row of
Figure S1A is a record the package already writes, `liulab_mbio.plot.linear` draws a record as a
bar of named feature arrows, and `liulab_mbio.plot.sequence_view` draws the junction call-outs at
base level.

What is missing is only the figure's arrangement: several records stacked as rows with dashed
correspondence lines drawn between one row's cut and the next row's end. Nothing in
`liulab_mbio.plot` lays out more than one record today.

## Where a source is thin

All four designs are covered. One deserves a flag rather than a fetch.

**The three-primer scheme has a picture but no paper.** `Baker Lab Strategy.png` is the only
description of it we hold. The `synthesis_and_assembly` README records this as an open gap: no
published sequence for the two universal flanking primers, and no source text beyond that image
and `dmx/jason_email.md`. Qian 2026 does not document the scheme. So the drawing can be copied
faithfully, but there is nothing to check its details against, and no citation for a reader who
wants more. The gap is in the research record, not in this note's ability to describe the picture
— and nothing in `reference_docs/` would close it, because the scheme appears to be unpublished.

## Sources

Qian, Z. et al. Accelerating protein design by scaling experimental characterization.
*Nat. Commun.* 17, 9945 (2026).
[doi:10.1038/s41467-026-76740-9](https://doi.org/10.1038/s41467-026-76740-9). Fig. 4 and
Supplementary Fig. 4.

Takacsi-Nagy, O. et al. Synthetic transcription factors designed by domain recombination enhance
CAR T cell antitumor function. *Cell* 189, 1–20 (2026).
[doi:10.1016/j.cell.2026.07.054](https://doi.org/10.1016/j.cell.2026.07.054). Figure S1A, in the
supplemental figures.

Lund, S., Potapov, V., Johnson, S. R., Buss, J. & Tanner, N. A. Highly parallelized construction
of DNA from low-cost oligonucleotide mixtures using Data-optimized Assembly Design and Golden
Gate. *ACS Synth. Biol.* 13, 745–751 (2024).
[doi:10.1021/acssynbio.3c00694](https://doi.org/10.1021/acssynbio.3c00694). Figures 1 and 2A.

Subramanian, S. K., Russ, W. P. & Ranganathan, R. A set of experimentally validated, mutually
orthogonal primers for combinatorially specifying genetic components. *Synth. Biol.* 3, ysx008
(2018). [doi:10.1093/synbio/ysx008](https://doi.org/10.1093/synbio/ysx008). The primer set the
three roles are allotted from.

The three-primer diagram and the correspondence behind it came from the Baker lab directly and are
not published. They are held as `long_fragment_GGA/Baker Lab Strategy.png` and
`dmx/jason_email.md`, and the correspondence is quoted with permission.

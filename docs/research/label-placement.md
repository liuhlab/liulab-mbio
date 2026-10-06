---
search:
  exclude: true
---

# Label placement: keeping labels apart on a map and a sequence view

Research note for issue #83. Everything below was read or measured on **2026-09-16**.

It asks what one rule keeps every label clear of every other label, and of the drawing, on the
circular map, the linear map and the sequence view of the planned `liulab_mbio.plot` module,
however crowded the record. The tool that draws is not chosen yet (#77, #78, #80), so the rule has
to hold whichever tool draws.

## 1. The verdict, first

| View | The rule | Guaranteed by construction | Measured |
| --- | --- | --- | --- |
| **Circular map** | Labels in two columns, right and left of the circle, in position order. Each column is spread to the nearest non-overlapping positions (least squares). Each label touches an ellipse round the drawing, and a straight leader runs from its site to that point. The ellipse grows until every leader is straight, and taller when a column does not fit. | No label overlaps a label or the drawing. No leader crosses another label or enters the drawing. | 0 violations on every real map tried and 400 random crowded ones. Leaders never crossed each other once all were straight; that is measured, not proven. |
| **Linear map**, whole or region | Labels in rows above the line, each with a vertical leader to its site. A label sits above every label whose site its own box reaches: a staircase. Names that do not fit on a feature can go below the line by the same rule, mirrored. | No label overlaps a label or the drawing, and no leader crosses another label. A layout always exists. | 0 violations. pUC19 with 108 labels took 59 rows; lambda's 188 feature names took 162. |
| **Sequence view** | The linear map's rows, per wrapped line, with names above their cut. Features and primers are bars, packed into rows by earliest start. | As for the linear map. Bars use the fewest rows possible. | 0 violations. The worst line of pUC19 took 4 rows for the shipped enzymes and 24 for its 99 unique 6+ cutters. |
| **Text width**, all three | Measure from one pinned font file's advance widths. Every output draws with that file, and SVG `textLength` pins each label as a backstop. | Every box is as wide as its text. | Chrome with the file embedded: within 0.09%. Without it: 16% narrower to 9.5% wider. |
| **Too many labels** | Merge the names at one cut, then grow the canvas up to a cap, then hide the lowest priority with a visible notice. | Nothing overlaps at any size. Nothing is dropped without saying so. | Merging alone cut pUC19's 105 labels to 44. |

**No tool examined guarantees this.** DNA Features Viewer and pyGenomeViz keep labels off each
other on a linear map but run leaders through them; pyCirclize, SeqViz and Open Vector Editor use
heuristics that overlapped or grouped; CGView.js, Circos and Vega hide what they cannot place, and
say nothing (section 5). So the rule is our own, in Python, and it runs before anything draws: it
hands the drawing code boxes and leaders in points. That fits the route #78 recommends, SVG written
by hand from Python and converted to PNG and PDF, and it fits matplotlib just as well.

**It is SnapGene's rule, written down.** SnapGene draws its circular map as two columns of labels
hugging the circle, its linear map and sequence view as staircases of rows, merges names at one cut
on a map, and hides what does not fit with a notice (section 3). Each matches a model from the
label-placement literature in which a clear layout can always be built (section 4). Placing labels
freely comes with no such promise.

**Width is the weak point, not geometry.** Every guarantee assumes a label is as wide as it was
measured. That holds only when every output draws with the font file that was measured (section 6).

**A test asserts on the layout's numbers, not on pixels.** Section 9 gives four assertions that run
in milliseconds.

## 2. How this was read and measured

- **Records**, read by `liulab_mbio.io.read_record`: pUC19 from SnapGene's plasmid collection
  (2,686 bp, 9 features), pUC19-GFP from `docs/examples/` (3,347 bp, 12 features, 4 primer binding
  sites), and phage lambda from the same collection (48,514 bp, 188 features) as an extreme.
- **Enzyme sets.** The 28 shipped enzymes, found with `sites.find_sites`, cut pUC19 once at 16
  sites, which crowds nothing. To crowd the map as a user naming enzymes would, three sets came from
  Biopython 1.88's `CommOnly` list of commercial enzymes: those with a site of 6 or more bases that
  cut once (99 on pUC19, 84 on pUC19-GFP; SnapGene calls this set "Unique 6+ Cutters"), and, as a
  stress case, every site of every commercial enzyme.
- **Throwaway environment** from conda-forge and bioconda, outside the repository: Python 3.13,
  matplotlib 3.11.1, uharfbuzz 0.56.1, DNA Features Viewer 3.1.5, pyCirclize 1.10.1, pyGenomeViz
  1.7.0. Browser widths came from Chrome 153 on macOS.
- **The prototype** laid labels out in points from measured widths, drew them with matplotlib at 200
  and 300 dpi, and checked every guarantee on the result (sections 7 and 9).
- **Measured or read.** Every count, width, size and time is measured and marked so. Algorithms and
  limits quoted from source code, documentation or papers are read, with a link. SnapGene was
  observed in images from its own help centre, Addgene's SnapGene-drawn maps and a screenshot the
  human shared. Times are one run on one laptop.

## 3. What SnapGene does

### 3.1 The circular map

Observed in SnapGene's help images
([Highlight an Enzyme Site](https://support.snapgene.com/hc/en-us/articles/10383876076948-Highlight-an-Enzyme-Site),
[Changes in 8.0](https://support.snapgene.com/hc/en-us/articles/34768892944916-Changes-to-Map-and-Sequence-Views-in-SnapGene-8-0))
and in Addgene's pUC19 map ([Addgene 50005](https://www.addgene.org/50005/sequences/)), which
Addgene says is "powered by GSL Biotech's SnapGene Server software"
([Addgene help](https://help.addgene.org/hc/en-us/articles/115005662726-How-does-Addgene-create-plasmid-maps)):

- **Two columns.** Sites on the right half are labelled right of the circle, left-aligned; sites on
  the left half, left of it, right-aligned. Addgene's map writes `AlwNI (236)` on the right and
  `(2443) AhdI` on the left.
- **In position order, hugging the circle.** A crowded run, such as the 16 labels round pUC19's
  multiple cloning site, becomes a contiguous stack whose inner edge follows a curve outside the
  circle. Straight grey leaders fan from the sites to the stack, and none visibly crosses another.
- **One label per cut position.** Enzymes cutting at one position share a label: `SacI - SstI
  (406)`, `AvaI - BsoBI - TspMI - XmaI (1039)`, and `KpnI - NcoI`, whose sites overlap so that
  both cut after the same base. The position is "the upstream 3' nucleotide of the top strand after
  digestion"
  ([coordinates](https://support.snapgene.com/hc/en-us/articles/10383760625172-Show-or-Hide-Restriction-Site-or-Primer-Coordinates-on-a-Map)).
- **Features on the arrow when the name fits**, otherwise a boxed name that joins the column, as
  `CAP binding site` and `lac operator` do. Overlapping features are tiled "in order of longest (top)
  to shortest (bottom)"
  ([Prioritize Display](https://support.snapgene.com/hc/en-us/articles/10383910479124-Prioritize-Display-of-Features-in-Map-View)),
  and font size changes "also affect feature widths (to accommodate larger feature labels)"
  ([Font Size](https://support.snapgene.com/hc/en-us/articles/10383794482836-Change-the-Font-Size)).
- **Primers** in purple, `name (start .. end)`, join the same columns.

### 3.2 What SnapGene hides, and how it says so

Read: "Feature labels that do not fit in the display area will be hidden. To see more feature
labels, either reduce the font size or enlarge the SnapGene window"
([Display Feature Labels](https://support.snapgene.com/hc/en-us/articles/10383722725524-Display-Feature-Labels-Below-or-Inside-a-Map)).
Since 8.0, "If SnapGene is unable to display all annotation items due to space constraints,
warnings (which can now be hidden) are now displayed in the bottom right corner of the window. You
can now click the 'Details' link"
([Changes in 8.0](https://support.snapgene.com/hc/en-us/articles/34768892944916-Changes-to-Map-and-Sequence-Views-in-SnapGene-8-0)).
The screenshot there reads `1 feature is hidden · Details`; the older window read
`1 feature is not displayed`. SnapGene does not abbreviate a name in any image seen.

### 3.3 The linear map

Observed in a SnapGene Viewer screenshot of a 9,951 bp plasmid the human shared, and in a help image
of NC_001709
([Prioritize Display](https://support.snapgene.com/hc/en-us/articles/10383910479124-Prioritize-Display-of-Features-in-Map-View)):
enzyme names, primer names and the boxed names of small features sit above the line in about 17
rows. Each has a vertical leader with a short elbow down to its position, and close neighbours rise
in a diagonal staircase, so no two overlap even where about 10 sites crowd 200 bp. Below the line,
features are arrows or boxes, named inside when the name fits and underneath otherwise. A circular
record drawn opened has `• • •` at each end.

This is the panorama labelling model of section 4.2: rows above a horizon, vertical leaders that
cross no other label, and labels that slide sideways along their leader, which is what the elbow
allows.

### 3.4 The sequence view

Observed in the sequence view image of
[Changes in 8.0](https://support.snapgene.com/hc/en-us/articles/34768892944916-Changes-to-Map-and-Sequence-Views-in-SnapGene-8-0),
pUC19's cloning site: names sit above their cut, each on a short tick. **Names at one cut are
stacked into a block** (`SstI` over `SacI`, `KpnI` over `XmaI`), not merged as on the map. A name
hangs right of its tick, or left of it, as both blocks do. A name that would cover a neighbour's
tick rises a row (`SmaI` over `BamHI`, `SalI` over `AccI`), and `HincII`
reaches its tick by a short elbow. The same panorama model, per line.

## 4. Published methods and what they guarantee

### 4.1 Free placement promises nothing; restricted models do

Placing labels anywhere next to their features "is commonly related to NP-hard geometric independent
set problems", while external labelling, with labels outside the drawing and leaders to them,
"focuses primarily on optimizing crossing-free, minimal-length leaders"
([Wallinger et al. 2026](https://arxiv.org/abs/2603.08657)). The survey of external labelling names
the usual hard constraints: "crossing-free leaders, non-overlapping labels, and avoidance of
obstacles in the illustration" ([Bekos, Niedermann, Nöllenburg 2019](https://arxiv.org/abs/1902.01454)).
Those three, plus width, are what section 7 guarantees.

### 4.2 Rows above a line: panorama labelling

[Gemsa, Haunert and Nöllenburg 2011](https://www1.pub.informatik.uni-wuerzburg.de/pub/haunert/pdf/GemsaEtAl2011.pdf)
(journal version [2015](https://dl.acm.org/doi/10.1145/2794299)) place "disjoint unit-height
rectangular labels" above a panorama "in k different rows", each "connected to its label by a
vertical leader that does not intersect any other label". A label may slide sideways as long as its
leader still meets it. Read:

- **A layout always exists.** Lemma 1: every instance has a feasible labelling with ⌈n/2⌉ rows,
  and some need that many. The construction is two staircases meeting in the middle.
- **Fewest rows** (MinRow) is solved exactly in O(k\*n³) time, and **most labels in K rows**
  (MaxLabels) in O(kn³), both by dynamic programming. Weighting labels by importance makes the
  second weakly NP-hard, with a pseudo-polynomial algorithm. "Solutions for realistically-sized
  instances are computed instantaneously."

SnapGene's linear map and sequence view fit this model (sections 3.3 and 3.4).

### 4.3 Labels round a circle

- **Focus regions.** [Fink et al. 2012](https://www1.pub.informatik.uni-wuerzburg.de/pub/fink/paper/fhssw-alfr-InfoVis12.pdf)
  place labels at the boundary of a circle with straight or Bézier leaders, with algorithms that
  "rule out crossings between leaders". Where there are "more sites than space for labels", the
  sites are prioritised, or clustered under one stacked label that expands on a click.
- **Orbital labelling.** Labels in a ring round a circular figure, each with "a short,
  crossing-free leader line". Minimising leader length is solved in polynomial time for labels of
  one size, but is "weakly NP-hard" for labels of different sizes
  ([Wallinger et al. 2026](https://arxiv.org/abs/2603.08657), building on
  [Bonerath et al. 2024](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.GD.2024.22)).
  Map labels are text of different widths, so no exact optimum is on offer there.
- **Order first.** For labels on one side, Bekos et al. 2007 work in two phases, as the
  [survey](https://arxiv.org/abs/1902.01454) describes: "First, sites and labels are matched so that
  they have the same vertical order. Second, leader crossings are iteratively resolved". SnapGene's
  columns are this model with the side bent round the circle.

None of these needs an optimum. What the map needs is a layout in which the hard constraints
cannot fail, and a good-looking one among those.

### 4.4 Spreading a column

Placing labels in a fixed order, each at least one line apart, as close as possible in least squares
to where their sites point, is isotonic regression once each target has the cumulative spacing
subtracted from it. The pool-adjacent-violators algorithm solves isotonic regression exactly in
linear time ([Best and Chakravarti 1990](https://link.springer.com/article/10.1007/bf01580873)).
Clipping the result to the column's top and bottom keeps both the order and the spacing: about 15
lines of Python in the prototype.

### 4.5 Bars in rows

Features and primers drawn as bars under the bases need no leaders. Taking them in order of start
and putting each in the first row it fits is optimal: the rows used equal the most bars covering
one base ("Earliest-start-time-first algorithm is optimal",
[Kleinberg and Tardos, section 4.1](https://www.cs.princeton.edu/~wayne/kleinberg-tardos/pdf/04GreedyAlgorithmsI.pdf)).

## 5. What the tools do

The measured columns are one drawing each of pUC19 with its 9 features and 99 unique 6+ cutters,
108 labels, and label boxes taken from matplotlib's text extents. A repeat run differed by under 3%.

| Tool | How it places labels (read) | Guarantee | Measured |
| --- | --- | --- | --- |
| DNA Features Viewer 3.1.5 | Measures each label with matplotlib, treats it as an interval, and gives it the lowest level no overlapping label holds, longest first ([source](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/blob/master/dna_features_viewer/compute_features_levels.py)). Long labels are cut "with an ellipsis" ([source](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/blob/master/dna_features_viewer/GraphicRecord/MatplotlibPlottableMixin.py)). | Labels apart on a linear map; leaders not | Linear: 0 overlapping pairs, 1,631 leader–label crossings, 3 labels over a feature. Circular: 55 pairs, 27 labels cut off the figure |
| pyCirclize 1.10.1 | Tries angle shifts, then radius steps, until a label clears those placed before it; stops after 1,000 tries, and skips everything above 200 labels with a warning ([source](https://github.com/moshi4/pyCirclize/blob/main/src/pycirclize/annotation.py), [config](https://github.com/moshi4/pyCirclize/blob/main/src/pycirclize/config.py)) | None | 371 overlapping pairs, 440 leader–label crossings, 3.2 s |
| pyGenomeViz 1.7.0 | Sorted by position, pushes each label up until it clears those before it, with no cap; above 200 labels per track it silently does nothing ([source](https://github.com/moshi4/pyGenomeViz/blob/main/src/pygenomeviz/genomeviz.py)) | Labels apart up to 200; leaders not | 0 pairs, 1,405 leader–label crossings |
| CGView.js 1.8.2 | Widths from the browser's `measureText`. The first 50 labels, favourites then largest features, move outward until clear; the rest are drawn only if clear ([source](https://github.com/stothard-group/cgview-js/blob/main/src/Annotation.js)) | Labels apart; the rest hidden without a notice | read only |
| SeqViz 3.10.24 | Labels on one side within 15 px grouped as `name,+N`, expanded on hover; groups of up to 4 spread if neighbours allow. Width estimated as 7.2 px per character ([source](https://github.com/Lattice-Automation/seqviz/blob/develop/src/Circular/Labels.tsx)) | None: horizontal overlap is not checked | read only |
| Open Vector Editor | Four quadrants, stacked at a fixed spacing along a circle; labels past reach merged into `+N,` groups, long labels cut with `..`. Width from character count, pinned in the page by `textLength` ([source](https://github.com/TeselaGen/tg-oss/blob/master/packages/ove/src/CircularView/Labels/relaxLabelAngles.js), [labels](https://github.com/TeselaGen/tg-oss/blob/master/packages/ove/src/CircularView/Labels/index.js)) | None stated | read only |
| Circos | "If a label's position results in overlap with another label, the label is drawn at the same angular position but is radially shifted out." A label that does not fit its ring "is not drawn", reported only in a debug log. Snuggling "is heuristic", and "Circos does not check overlap between elements in the label track and scale ticks" ([text 1](https://circos.ca/documentation/tutorials/2d_tracks/text_1/), [text 2](https://circos.ca/documentation/tutorials/2d_tracks/text_2/)) | Labels apart in a track; the rest hidden without a notice | read only |
| Vega `label` transform | A label that fits nowhere gets opacity 0 ([#78's note](https://github.com/liuhlab/liulab-mbio/blob/research/offline-html-maps/docs/research/offline-html-maps.md), [docs](https://vega.github.io/vega/docs/transforms/label/)) | Labels apart; the rest hidden without a notice | read only |

plotly, Bokeh and Altair place no labels: they drew #78's enzyme names over each other from the same
coordinates ([#78's note](https://github.com/liuhlab/liulab-mbio/blob/research/offline-html-maps/docs/research/offline-html-maps.md)).

**One measuring trap.** A matplotlib `Annotation`'s window extent includes its arrow. Counting
pyCirclize's overlaps with that extent gave 1,877 pairs; with `Text.get_window_extent`, the text
alone, 363. A test on a matplotlib figure must measure the text, not the annotation.

## 6. Measuring text without a browser

A layout is only as good as its widths. Measured on 425 labels, every feature, primer and site name
of the five records, in DejaVu Sans at 8 pt:

| Compared | Difference (measured) |
| --- | --- |
| Chrome, font file embedded in the page, against HarfBuzz shaping of the same file | at most 0.09% (0.004 px) |
| matplotlib's PDF output metrics against HarfBuzz | at most 0.15% |
| matplotlib's PNG renderer (Agg) at 300 dpi against HarfBuzz | at most 3.1% (2.8 pt), mean 1%; at 72 dpi up to 15% on short labels |
| Chrome with matplotlib's SVG font list and no embedded font, drawing at Lucida Grande's widths on macOS, against DejaVu Sans | 16% narrower to 9.5% wider; 24 labels wider by over 0.5%, by up to 5.2 px |
| Arial or Helvetica against DejaVu Sans | median 11% narrower, up to 70% wider (`J`) |
| Verdana against DejaVu Sans | median 4% wider, up to 54% wider |

HarfBuzz "is used in Android, Chrome, ChromeOS, Firefox" and much else
([HarfBuzz](https://github.com/harfbuzz/harfbuzz)), so shaping the font file in Python predicts what
a browser draws with that file. #78 measured the same gap from the other side: Chrome's default sans
drew its 42 labels "between 15% narrower and 5% wider than matplotlib measured them", while
monospaced fonts agree at 0.602 of the font size per character
([#78's note](https://github.com/liuhlab/liulab-mbio/blob/research/offline-html-maps/docs/research/offline-html-maps.md)).

**The converters do not read a font from the page.** CairoSVG picks a font by family name through
cairo's `select_font_face`
([source](https://github.com/Kozea/CairoSVG/blob/main/cairosvg/text.py)), so it draws with
whatever the machine has installed. vl-convert "uses installed system fonts and accepts additional
font directories" ([README](https://github.com/vega/vl-convert/blob/main/vl-convert-python/README.md)).
Both honoured `textLength` in #78's measurement.

**So the widths hold only under three conditions.** Measured and read together, they are:

1. **Pin one font file** and measure every label from its advance widths: fontTools or HarfBuzz for
   hand-written SVG, or matplotlib's own metrics if matplotlib draws.
2. **Make every output draw with that file.** The HTML page embeds it: DejaVu Sans cut down to 101
   characters, printable ASCII and the few others these labels use, was 26.8 KB as TTF and 14.2 KB
   as WOFF2 (measured). The converter is given the same file. DejaVu's licence allows redistribution and embedding; it forbids selling
   the fonts alone and requires renaming a modified font and keeping the notice
   ([licence](https://dejavu-fonts.github.io/License.html)). A shipped font is package data, so the
   repository's rule applies: rebuilt by a script, sourced in a note.
3. **Pin each label with `textLength`.** "The user agent will ensure that the text does not extend
   farther than that distance"
   ([MDN](https://developer.mozilla.org/en-US/docs/Web/SVG/Reference/Attribute/textLength)). If a
   font is ever substituted, glyphs squeeze or spread, but the box stays true.

Padding each label by 1.5 pt a side at 7 pt absorbed the PNG renderer's hinting: no label drawn
at 200 or 300 dpi ran past its box in any prototype layout (measured).

## 7. The rule, prototyped and checked

The prototype is plain geometry in points. matplotlib only measured and drew. Section 9 lists
the checks behind the "violations" columns.

### 7.1 Circular map

1. Features overlapping on the ring go to separate rings, longest outermost, as SnapGene tiles them.
   A name goes on its arrow when its measured box fits inside the arrow; otherwise it becomes a boxed
   outside label.
2. Names at one cut merge into `A - B (position)`.
3. A label goes to the right column if its site is on the right half, else the left.
4. Each column is spread in position order: targets where the sites point, at least one line apart,
   least squares by pool-adjacent-violators, clipped to the column's height (section 4.4).
5. Each label's inner edge touches an ellipse whose semi-axes start just outside the drawing; the
   touching point is the label's port. A straight leader runs from the site to the port.
6. The ellipse widens until no straight leader dips into the drawing, and heightens until each
   column fits. Past a cap, labels are hidden lowest priority first, enzyme sites before primers
   before features, and a notice names the count. A leader that still dips at the cap goes round
   by short chords outside the drawing.

**Why it cannot fail.** Right labels lie at x ≥ +3 pt and left labels at x ≤ −3 pt, so the columns
never meet, and the spread keeps each pair within a column a line apart. Every label lies outside
the ellipse and the whole drawing inside it, so no label touches the drawing. Every leader joins a
point inside the ellipse to a point on it, and an ellipse is convex, so the leader stays inside and
meets no other label. The growth in step 6 keeps leaders out of the drawing. Leaders crossing each
other is not ruled out by this argument.

Measured. The real maps kept the ellipse at its first width, just outside the drawing, and routed
dipping leaders by chords; the merged pUC19 map was also widened, and the random maps always were.

| Map | Labels | Violations | Leader crossings | Canvas |
| --- | --- | --- | --- | --- |
| pUC19, 99 unique 6+ cutters, names apart | 105 | 0 | 704, among the chord leaders | 5.9 × 11.0 in; the ellipse 2.5 times taller than wide |
| the same, names merged | 44 | 0 | 152; 0 once widened until all leaders were straight | 10.2 × 4.8 in; widened, 11.4 × 6.2 in |
| the same, names apart, height capped | 81 shown, 24 hidden, notice drawn | 0 | 464 | 5.9 × 7.8 in |
| pUC19-GFP, 84 cutters and 4 primers, names apart | 96 | 0 | 0 | 6.6 × 10.6 in |
| the same, names merged | 39 | 0 | 0 | 9.8 × 4.6 in |
| lambda's 188 features drawn circular | 172 outside, 16 on arrows | 0 | 327 | 6.6 × 20.8 in |
| 400 random maps, 10 to 100 sites, half clustered, merged, ellipse widened until straight | up to 100 | 0 | 0 | — |

The merged pUC19 map, once all leaders were straight, looks like Addgene's SnapGene map of the same
plasmid. Layout took 7 to 31 ms per real map; the slowest random one, with the search for the
ellipse, 245 ms.

### 7.2 Linear map

Each label hangs right from its site. Its row is one above the highest row among the labels whose
sites its own box reaches, and the labels to its right whose boxes it meets. Labels at one position
are taken in order. A crowded run therefore rises in a staircase. The canvas keeps a right margin as
wide as the widest label.

**Why it cannot fail.** A label whose box reaches another's site always sits above it, so no leader
passes through a lower label, and two labels whose boxes meet never share a row. These are
Gemsa's requirements F1 to F3. A layout always exists, since at worst each label takes its own row.

Measured, the axis 900 pt wide:

| Map | Labels | Violations | Rows, staircase | Rows, two staircases per cluster |
| --- | --- | --- | --- | --- |
| pUC19, 99 unique 6+ cutters and 9 feature names | 108 | 0 | 59 | 47 |
| the same, names merged | 47 | 0 | 25 | 24 |
| lambda, 188 feature names | 188 | 0 | 162 | 159 |

The staircase took 0.3 to 5 ms. Trying every split of each cluster into a left-hanging and a
right-hanging staircase took up to 836 ms and saved little; Gemsa's dynamic program is the exact
answer if rows ever matter more. Lambda shows the limit of the model: a whole-genome map cannot show
every name above the line, so the linear map needs the row cap and notice of section 8.

### 7.3 Sequence view

The same rows, per 60-base line, with names but no positions, as SnapGene shows them. The line grows to hold its rows.

| Record, enzymes shown | Sites | Most rows on a line, names apart | Merged |
| --- | --- | --- | --- |
| pUC19, shipped enzymes | 20 | 4 | 3 |
| pUC19, unique 6+ cutters | 99 | 24 | 7 |
| pUC19, every commercial site | 4,719 | 202 | 49 |
| pUC19-GFP, shipped enzymes | 12 | 2 | 1 |
| pUC19-GFP, unique 6+ cutters | 84 | 12 | 4 |
| pUC19-GFP, every commercial site | 5,327 | 145 | 49 |

0 violations everywhere, and a whole record took 1 to 7 ms for the first two sets (measured). Up to
13 names per record ran past the end of a line, into the right margin.

## 8. When labels outnumber the space

In this order:

1. **Merge names at one cut position**, as SnapGene's map does: pUC19's 105 outside labels became
   44, and pUC19-GFP's 96 became 39 (measured). It loses nothing, since the enzymes cut at the same
   base.
2. **Grow the canvas**, up to a cap. A PNG or PDF can be taller; the measured costs are in section 7.
3. **Hide the lowest priority, and say so.** Keep features, then primers, then sites. Write the
   count on the drawing, as SnapGene writes `1 feature is hidden`, and in HTML list the hidden names
   in the notice.

**Abbreviating is not a fix for enzymes.** A shortened enzyme name no longer identifies the enzyme,
and a shorter label still needs its own row. DNA Features Viewer and Open Vector Editor truncate; SnapGene does not.

**Grouping into `name,+N` is an extra for HTML, not the rule.** SeqViz and Open Vector Editor expand
a group on hover, which a PNG cannot. The page and the PNG should share one layout, so hover can
reveal what the notice counts but should not change the layout.

## 9. How a test asserts it

The layout returns every label's box and leader, and the drawing's extent, in points. A test on a
crowded fixture asserts:

1. **No two label boxes overlap**: pairwise, touching allowed.
2. **No label overlaps the drawing**: on the circular map, each box's nearest point lies outside the
   outer radius; on a linear map or sequence line, each box lies above the line.
3. **No leader crosses a label other than its own**: segment against box.
4. **No leader enters the drawing**: each segment's closest approach to the centre stays outside the
   outer radius.

Leader crossings are counted rather than asserted, since section 7 measured but did not prove them
absent. Inline names get one more check, that each box lies inside its own arrow.

For fixtures, pUC19 with its unique 6+ cutters is crowded enough, plus a few seeded random maps. The
four checks took under 30 ms on 39 to 252 labels (measured), so they cost the gate nothing. Checking
pixels would be slower and could only find an overlap, not explain it. Width needs no pixels either:
measure with the pinned font in the test.

If matplotlib draws, the same checks run on `Text.get_window_extent`, the text alone, never on an
annotation's extent (section 5).

## 10. Judgement calls for the interface ticket (#81)

- **The caps**: how tall a circular map may grow, how many rows a linear map may take, and whether
  the sequence view ever caps.
- **Hiding priority**, and the notice's wording and place: SnapGene's bottom-right count with
  details, or a line under the map.
- **Merging in the sequence view.** SnapGene stacks names at one cut there but merges them on the
  map.
- **Feature names on a linear map**: below the line, as SnapGene, or in the rows above with the
  sites.
- **Names on arrows**: straight, as prototyped, or curved along the arc, as SnapGene; and whether
  arrows widen with the font.
- **The font and its size**: which file is pinned, and whether it ships as package data or comes
  from a required dependency such as matplotlib.
- **Names near a line's end**: a right margin, as prototyped, or hanging left.
- **Two staircases per cluster or Gemsa's exact rows**: worth it only if rows are tight.
- **Hover in HTML**: whether hovering a site or the notice reveals hidden names.

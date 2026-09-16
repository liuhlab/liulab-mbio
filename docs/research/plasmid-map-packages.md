---
search:
  exclude: true
---

# Plasmid map packages: which Python libraries draw a map, and how close do they get to SnapGene?

Research note for issue #77. It surveys the Python packages that draw an annotated sequence
as a circular or linear map, and asks how close each comes to the reference picture on map
issue #76: arrows along the strand, a label on the arrow when it fits and outside with a leader
line when it does not, enzyme sites outside as `name (position)`, primers outside as
`name (start .. end)`, switches to show or hide each layer, PNG and PDF, and one HTML file that
opens offline. It also asks which of them draws a base-by-base sequence view like SnapGene
Viewer's Sequence tab.

Everything below was read or measured on **2026-09-16**. Each section says which facts were
**measured** — drawn, installed or counted — and which were **read** from a package's own
source, docs or metadata.

## 1. The verdict, first

| Package | Draws | Nearest to the reference in | Leaves us to build | Carry into #80? |
| --- | --- | --- | --- | --- |
| **gbdraw** | circular, linear, a linear region | Closest picture out of the box: strand arrows, a label on the arc when it fits and outside with a leader line when it does not, name and length in the centre, an offline interactive SVG with popups | Sites and primers passed in as features; overlap is not prevented; no show/hide in the file; not on PyPI; its Python API is still being renamed | **Yes**, as the benchmark picture |
| **pyCirclize + pyGenomeViz** | circular (pyCirclize), linear and region (pyGenomeViz) | Lightest maintained pair on conda-forge and PyPI, MIT, one author, both on matplotlib; primitives that take our own spans; pyGenomeViz writes an offline HTML viewer | The label rule on the circle, splitting a span at the origin, show/hide in HTML, any circular HTML | **Yes**, as the matplotlib route |
| DNA Features Viewer | circular, linear, a cropped region | Label on the feature when it fits, otherwise outside; the one sequence view with a translation | Circle labels stack in a column above the ring; segments collapse to one span; no commits in a year; bioconda only | No — read its label and translation code |
| Biotite | circular, linear | Curved labels on the arc; a feature across the origin drawn as one arc | No outside labels at all: a label that does not fit is dropped | No |
| Biopython GenomeDiagram | circular, linear, a region | Already inside Biopython | No label placement; ReportLab adds 36 packages; untouched since 2024 | No |
| QUEEN | circular, linear | The only package that draws both strands base by base | The released version does not import on Python 3.13; PyPI only; no commits since 2024 | No — a reference for the sequence view |
| pLannotate | circular, in Bokeh | Hover details, a legend by feature type | GPL-3.0; takes its own BLAST hit table, not a record; installs BLAST, DIAMOND and Perl | No |
| DNAplotlib, plasmidcanvas, others | see section 9 | — | — | No |

**No package keeps labels apart on a circle with a guarantee.** Every circular candidate measured
that places outside labels overlapped some on the crowded map (section 4). The rule is ours to write,
which is issue #83's question.

**No package's HTML has show/hide switches.** The offline files that exist give hover details
and zoom; none lets a reader hide a layer or a feature type (section 6). That is issue #78's.

**No package draws SnapGene's Sequence tab.** The pieces exist separately — one strand with a
translation, both strands without one, enzyme names over both strands as text — but no package
joins them (section 7). The sequence view is ours to build.

**No package knows what an enzyme site or a primer is.** Every candidate takes them as generic
features or free-text annotations, so the `name (position)` and `name (start .. end)` labels,
and the purple, are ours in any route.

## 2. How this was measured

**The record.** `docs/examples/pUC19-GFP/product.dna`, written to GenBank: 3,347 bp, circular,
12 features, two of them joined from segments (lacZα and AmpR), 4 primers. Enzyme sites came
from `sites.find_sites` over the shipped enzymes, keeping those that cut once: NcoI, BtgZI,
BspQI and SapI.

**Three variants of it.**

- *plain* — the record as it is.
- *crowded* — plain, plus pUC19's own 14 multiple-cloning-site enzymes packed into 60 bp, so
  the map carries 30 labels, 14 of them within under 2% of the circle.
- *rotated* — the record's origin moved 2,700 bp, so AmpR crosses the origin as
  `complement(join(2934..3347,1..378,379..447))`.

A fourth, *harsh*, added 24 more site labels in 48 bp near position 2,000, for gbdraw only.

**The environment.** Every conda candidate was installed together from conda-forge and bioconda
on Python 3.13.15 and drew the variants; the PyPI-only candidates went into a second throwaway
environment. Dependency weight was measured by solving each candidate beside the
repository's own runtime set (Python 3.13, typer, biopython, primer3-py) for both platforms and
counting the packages the lock gained.

**Counting overlaps.** For the matplotlib packages, the figure was drawn and every pair of label
boxes compared in display coordinates — the text's own extent, or the rounded box around it where
the package draws one. gbdraw writes SVG, so its overlaps were read off the rendering.

## 3. The map views

Measured on the three variants unless marked *read*.

| | Circular | Linear | Region | Strand arrows | Per-feature colour | Show/hide | Input it takes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **gbdraw** 0.13.0 | yes | yes | yes, `--region`, 1-based | yes | yes, TSV rules by qualifier | at draw time: feature types, a visibility table, label whitelist | GenBank, or GFF3 + FASTA |
| **pyCirclize** 1.10.1 | yes | no | a sector may cover part of a record, still circular | yes | yes | no | its GenBank/GFF parsers, SeqFeatures, or plain coordinates |
| **pyGenomeViz** 1.7.0 | no | yes | yes, segments | yes | yes | no | its GenBank/GFF parsers, SeqFeatures, or plain coordinates |
| DNA Features Viewer 3.1.5 | yes | yes | yes, `crop()` | yes | yes | by filtering features in a translator | a Biopython record, a GenBank or GFF path, or plain features |
| Biotite 1.7.1 | yes | yes | yes, `loc_range` | yes | yes, a formatter function | no | Biotite's own `Annotation` |
| GenomeDiagram (Biopython 1.88) | yes | yes | yes, `start`/`end` | yes, `ARROW`, `BIGARROW` sigils | yes | by what is added to a track | SeqFeatures |
| QUEEN 1.2.0 (unreleased) | yes | yes | yes, `start`/`end` | yes | yes, by qualifier | a feature list | its own `QUEEN` object, read from GenBank |
| pLannotate 2.0.0 (*read*) | yes | no — `linear=True` only marks a break | no | yes | by feature type | no | its BLAST hit table |

Static formats: every matplotlib package writes PNG, PDF and SVG. GenomeDiagram writes PDF, SVG and
EPS, and PNG through ReportLab's bitmap renderer
([Biopython tutorial](https://biopython.org/docs/latest/Tutorial/chapter_graphics.html)). gbdraw
writes SVG, and PNG, PDF, EPS and PS through CairoSVG
([gbdraw README](https://github.com/satoshikawato/gbdraw/blob/8675f33bb2f2ac8f174f65e9a36529df942c5333/README.md)).
All were produced in the measurement.

### 3.1 How each one looks against the reference

- **gbdraw** is the nearest out of the box. On the crowded map it drew each feature as an
  arrow on its strand ring, put `lacZα`, `GFP`, `AmpR` and `ori` on their arcs, fanned the 18
  site labels outside with leader lines, and wrote the name, length and GC content in the centre.
  Its default feature types are `CDS`, `rRNA`, `tRNA`, `tmRNA`, `ncRNA`, `misc_RNA` and
  `repeat_region`
  ([Python API reference](https://github.com/satoshikawato/gbdraw/blob/8675f33bb2f2ac8f174f65e9a36529df942c5333/docs/REFERENCE/python-api.md)),
  so promoters, sites and primers appear only when named. It chose `/product` or `/note` over
  `/label` until given a qualifier priority table, and drew primers as plain boxes; the shape is
  set per feature type.
- **pyCirclize** gives rings, arrows, ticks, centre text and `annotate`, a label with a leader
  line. It draws what it is told, so the picture is as close as the calling code makes it.
- **pyGenomeViz** gives the same primitives on a line, and stacks annotations upward until they
  clear each other. On the region view, the axis counts from the region's start rather than the
  record's.
- **DNA Features Viewer** puts a label inside its feature when it fits and outside with a line
  otherwise. On a circle the outside labels do not fan around the ring: they stack in a vertical
  column above it, with leader lines running down to the features.
- **Biotite** writes each label as curved text on its arc and drops the label when it is longer
  than the arc
  ([`plot_plasmid_map`](https://www.biotite-python.org/latest/apidoc/biotite.sequence.graphics.plot_plasmid_map.html)).
  A site or primer, being short, got no label at all. Its linear map overlapped long labels.
- **GenomeDiagram** places a label at a fixed position and angle per feature. The crowded site
  labels piled into one black mass on both views.
- **QUEEN** draws features on stacked rings inside the circle with radial labels, the reverse of
  the reference's outside labels.

## 4. Labels: who keeps them apart

One fact each — issue #83 goes deep.

| Package | What it does | A guarantee? | Measured on the crowded map |
| --- | --- | --- | --- |
| gbdraw | fans outside labels with leader lines; with `--labels both`, some go inside the ring | **No.** Its FAQ answers "My labels overlap" with: reduce the label size, filter the label set, or move the track ([FAQ](https://github.com/satoshikawato/gbdraw/blob/8675f33bb2f2ac8f174f65e9a36529df942c5333/docs/FAQ.md#L102-L108)) | the site labels clear each other; one primer label laid over the centre GC text. On *harsh*, at least four site-label pairs overlapped as well |
| pyCirclize | shifts each annotation along the circle, then outward, up to 1,000 steps, and skips adjustment past 200 annotations ([`adjust_annotation`](https://github.com/moshi4/pyCirclize/blob/5a0f36111a4bbfab3e3d765e7365a1108f891dcb/src/pycirclize/annotation.py#L19), [limits](https://github.com/moshi4/pyCirclize/blob/5a0f36111a4bbfab3e3d765e7365a1108f891dcb/src/pycirclize/config.py#L54-L70)) | **No** — capped iterations, and an open request for bent leader lines ([#66](https://github.com/moshi4/pyCirclize/issues/66)) | **18** overlapping pairs among 28 labels |
| pyGenomeViz | moves each annotation up until its box clears every earlier one, with no iteration cap, up to 200 per track ([`_adjust_annotation`](https://github.com/moshi4/pyGenomeViz/blob/fb137ded78d56b1d55229309d500aa54c80e3abe/src/pygenomeviz/genomeviz.py#L732-L767), [limit](https://github.com/moshi4/pyGenomeViz/blob/fb137ded78d56b1d55229309d500aa54c80e3abe/src/pygenomeviz/config.py#L4-L16)) | **Text boxes, yes, on a line**; leader lines still cross other labels | 0 of 32 |
| DNA Features Viewer | assigns labels to levels by graph colouring so overlapping spans never share one ([`compute_features_levels`](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/blob/049bbe4e3063e90ae9b0f88ac7e92f47b735d38a/dna_features_viewer/compute_features_levels.py#L30)) | **Linear, in practice; circular, no.** Labels overlap when records share an axis ([#64](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/issues/64)) | linear 0 of 30 on all four maps; circular 1 on *crowded*, 1 on *rotated*, 5 on *rotated crowded* |
| Biotite | drops a label that does not fit its arc ([source](https://github.com/biotite-dev/biotite/blob/1aaf784a826a828fde989ada69e37e05013ac3fc/src/biotite/sequence/graphics/plasmid.py#L91-L96)) | **Yes, by omission**; a request to move labels outward instead is open since 2023 ([#459](https://github.com/biotite-dev/biotite/issues/459)) | no outside labels drawn |
| GenomeDiagram | none | **No** | labels piled on both views |
| QUEEN | README says labels are "automatically adjusted to prevent overlaps" ([README](https://github.com/yachielab/QUEEN/blob/9798bf0a88f08c51f5327575d5f1707f9ce198ad/README.md#L1459-L1461)) | not measured on *crowded* | — |
| pLannotate | a fixed offset per feature, by quadrant ([source](https://github.com/mmcguffi/pLannotate/blob/eb043a4fbf16ddca9206451c490c9dbc38dcdb0e/plannotate/bokeh_plot.py#L275-L300)) | **No** (*read*) | — |

**pyGenomeViz's loop is the one sufficient rule found, and it works because a line has a free
direction.** Stacking upward always ends. A circle has no such free direction without pushing
labels arbitrarily far out, which is why the circular candidates cap their search and give up.

## 5. Across the origin, and joined features

SnapGene files hold a feature as segments, and a circular record can carry one across the
origin. Both matter to the reference picture, and the packages differ sharply. Measured on
*plain* (lacZα and AmpR joined) and *rotated* (AmpR also across the origin).

| Package | Joined feature | Feature across the origin |
| --- | --- | --- |
| gbdraw | segments in place, joined by a thin grey arc | correct, continuous through position 1 |
| Biotite | each segment in place | correct, continuous through position 1 |
| QUEEN | one outer span | correct, continuous through position 1 |
| GenomeDiagram | each segment its own arrow | correct place, but each piece gets its own arrowhead |
| pyGenomeViz | segments joined by an intron line | the two ends joined by a line running the whole map |
| DNA Features Viewer | collapsed to one outer span, because its translator reads only `location.start` and `location.end` ([source](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/blob/049bbe4e3063e90ae9b0f88ac7e92f47b735d38a/dna_features_viewer/BiopythonTranslator/BiopythonTranslatorBase.py#L32-L55)) | **a full ring** — open as [#29](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/issues/29) |
| pyCirclize `genomic_features` | **wrong on the reverse strand**: it takes the first part's start and the last part's end ([source](https://github.com/moshi4/pyCirclize/blob/5a0f36111a4bbfab3e3d765e7365a1108f891dcb/src/pycirclize/track.py#L1440-L1447)), but Biopython lists a reverse feature's parts in biological order, so AmpR drew with zero length and lacZα from the wrong ends | a wrong arc |

For the matplotlib route this matters less than it looks: `liulab_mbio` already holds each
feature as segments with the coordinate rule of ADR 0001, so a `plot` module would call
pyCirclize's `arrow` and pyGenomeViz's `add_feature` per piece and split a span at the origin
itself, never going through their SeqFeature readers.

## 6. HTML output, as far as each package offers it

Issue #78 compares HTML routes in depth. Here, only what each package writes by itself.

| Package | File | Opens offline | Details | Zoom | Show/hide |
| --- | --- | --- | --- | --- | --- |
| gbdraw | interactive SVG | **yes** — self-contained style and script ([docs](https://github.com/satoshikawato/gbdraw/blob/8675f33bb2f2ac8f174f65e9a36529df942c5333/docs/REFERENCE/output-formats-and-export.md#L36-L48)); measured: no external reference | popups | yes | no |
| pyGenomeViz | HTML viewer, linear only ([`savefig_html`](https://github.com/moshi4/pyGenomeViz/blob/fb137ded78d56b1d55229309d500aa54c80e3abe/src/pygenomeviz/genomeviz.py#L617)) | **yes** — its JavaScript libraries are embedded; measured: no external reference | a tooltip table on hover | yes | no — recolour, relabel, a filterable feature table, PNG/SVG download |
| DNA Features Viewer | Bokeh, **linear only**: the patches are straight arrows even for a circular record ([source](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/blob/049bbe4e3063e90ae9b0f88ac7e92f47b735d38a/dna_features_viewer/GraphicRecord/BokehPlottableMixin.py#L23-L59)); needs Bokeh and pandas, neither declared | **yes** with Bokeh's inline resources; measured: the one CDN address inside is a lazy MathJax loader used only for LaTeX labels | the label on hover | x-axis only | no |
| pLannotate | Bokeh, circular | only with `--htmlfull`; `--html` loads Bokeh from a CDN ([source](https://github.com/mmcguffi/pLannotate/blob/eb043a4fbf16ddca9206451c490c9dbc38dcdb0e/plannotate/models.py#L477-L496)) | name, type, identity on hover | yes | no |
| pyCirclize | none; tooltips only in a live Jupyter kernel ([README](https://github.com/moshi4/pyCirclize/blob/5a0f36111a4bbfab3e3d765e7365a1108f891dcb/README.md#L268-L280)) | — | — | — | — |
| Biotite, GenomeDiagram, QUEEN | none | — | — | — | — |

## 7. The sequence view

The reference is SnapGene Viewer's Sequence tab: both strands base by base, a position ruler,
features as bars under the bases, primers, enzyme names above their cut, and the translation of
an annotated CDS. **No Python package draws it.** Three draw a part.

| Package | Both strands | Ruler | Features | Primers | Enzyme names at the cut | Translation | Output |
| --- | --- | --- | --- | --- | --- | --- | --- |
| DNA Features Viewer | **top strand only** ([`plot_sequence`](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/blob/049bbe4e3063e90ae9b0f88ac7e92f47b735d38a/dna_features_viewer/GraphicRecord/SequenceAndTranslationMixin.py#L5)) | yes | arrows above the bases | as features | as features, a box over the site | **yes**, any range and strand ([`plot_translation`](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/blob/049bbe4e3063e90ae9b0f88ac7e92f47b735d38a/dna_features_viewer/GraphicRecord/SequenceAndTranslationMixin.py#L63)) | PNG, PDF, SVG; wrapped across pages of a PDF ([`plot_on_multiple_pages`](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/blob/049bbe4e3063e90ae9b0f88ac7e92f47b735d38a/dna_features_viewer/GraphicRecord/MultilinePlottableMixin.py#L103)) |
| QUEEN | **yes**, as coloured letter cells, `seq` and `rcseq` ([README](https://github.com/yachielab/QUEEN/blob/9798bf0a88f08c51f5327575d5f1707f9ce198ad/README.md#L1459-L1520)) | yes | arrows above the bases | as features | no | no | PNG, PDF, SVG, wrapped by `linebreak` |
| Biopython `Restriction` | **yes**, as text | start and end of each line | no | no | **yes**, name and position over the cut ([`print_as("map")`](https://github.com/biopython/biopython/blob/5bbc6c12c505301f2d681f932c30fdb8fcbe9a6e/Bio/Restriction/PrintFormat.py#L96-L106)) | no | plain text |

All three were measured on the pUC19 multiple cloning site. DNA Features Viewer and QUEEN put
the features above the sequence rather than as bars under it. The sequence view's own ticket
should take DNA Features Viewer's translation placement and Biopython's enzyme-name stacking as
worked examples, not as dependencies.

## 8. Licence, maintenance and install

Releases and activity read from PyPI, anaconda.org and the GitHub API; commit counts are for the
12 months to 2026-09-16, on the file or folder named. "Added" is measured: packages the lock gains
over the repository's runtime set, which is 39 on osx-arm64 and 42 on linux-64.

| Package | Licence | Last release | Commits, 12 months | Python 3.13 | conda | PyPI | Added (osx-arm64 / linux-64) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gbdraw | MIT | 0.13.0, 2026-07-05; source at 0.14.0 | **1,439** | measured | bioconda, noarch | **no** | 53 / 56: CairoSVG (LGPL-3.0), cairo, pandas, bcbio-gff, svgwrite |
| pyCirclize | MIT | 1.10.1, 2025-10-03 | 19 | measured | conda-forge, noarch | yes | 53 / 58: matplotlib, pandas |
| pyGenomeViz | MIT | 1.7.0, 2026-06-13 | 40 | measured | conda-forge, noarch | yes | 52 / 57: matplotlib |
| DNA Features Viewer | MIT | 3.1.5, 2025-05-09 | **0** | measured | **bioconda only** | yes | 52 / 57: matplotlib |
| Biotite | BSD-3-Clause | 1.7.1, 2026-06-21 | 4 in `sequence/graphics` | measured | conda-forge | yes | 67 / 72 with matplotlib: networkx, requests, biotraj |
| GenomeDiagram | Biopython licence | Biopython 1.88, 2026-08-06 | **0** in `GenomeDiagram`; newest 2024-06-26 | measured | conda-forge | yes | 36 / 39: ReportLab, pycairo, cairo |
| QUEEN | MIT | `python-queen` 1.1.0, 2022-08-09 | **0**; newest 2024-10-09 | **1.1.0 fails**: pins `requests~=2.23`, whose urllib3 does not import; unreleased 1.2.0 runs once current requests is forced | no | yes | — |
| pLannotate | **GPL-3.0** | 2.0.0, 2026-07-05 | 3 on `bokeh_plot.py` | not measured | bioconda | no | 76 / 76: BLAST, DIAMOND, Infernal, Perl, Bokeh, pandas |

Matplotlib alone adds 51 / 56 of those packages, so pyCirclize and pyGenomeViz together cost
matplotlib plus pandas. Bokeh alone adds 27 / 26.

Sources: [pyCirclize](https://pypi.org/project/pycirclize/),
[pyGenomeViz](https://pypi.org/project/pygenomeviz/),
[DNA Features Viewer](https://pypi.org/project/dna-features-viewer/) and its
[bioconda recipe](https://anaconda.org/bioconda/dna_features_viewer),
[Biotite](https://anaconda.org/conda-forge/biotite),
[gbdraw on bioconda](https://anaconda.org/bioconda/gbdraw),
[GenomeDiagram's history](https://github.com/biopython/biopython/commits/master/Bio/Graphics/GenomeDiagram),
[QUEEN's `setup.py`](https://github.com/yachielab/QUEEN/blob/9798bf0a88f08c51f5327575d5f1707f9ce198ad/setup.py#L15-L24),
[pLannotate's licence](https://github.com/mmcguffi/pLannotate/blob/eb043a4fbf16ddca9206451c490c9dbc38dcdb0e/LICENSE).

**Three install facts bear on the choice.**

- **gbdraw is not on PyPI.** `liulab-mbio` lists its runtime dependencies in `pyproject.toml`
  for its wheel, so a bioconda-only library would be a pixi-only dependency, as `ipcr` is.
- **gbdraw's API is moving.** 1,439 commits in a year, and its reference already rejects
  "retired flat label names"
  ([Python API](https://github.com/satoshikawato/gbdraw/blob/8675f33bb2f2ac8f174f65e9a36529df942c5333/docs/REFERENCE/python-api.md#L130-L134)).
  The package-root `draw_circular` in the 0.14.0 docs does not exist in the 0.13.0 build on
  bioconda (measured).
- **pLannotate's licence rules it out as code or dependency.** GPL-3.0 code cannot be copied into
  an MIT package.

## 9. Also found, and why each is out

| Package | What it is | Why it is out |
| --- | --- | --- |
| [DNAplotlib](https://github.com/VoigtLab/dnaplotlib) | SBOL Visual diagrams of a genetic design: glyphs in order, not a position-scaled map; linear only | No circular map; PyPI 1.0 from 2017; no commits since 2024-03; not on conda. Its README says OSL-3.0 while its `LICENSE.txt` says MIT ([README](https://github.com/VoigtLab/dnaplotlib/blob/7bed9fac3f2ea495e44086c2defdf0fe91a1a376/README.md#L3)). Imports on 3.13 (measured) |
| [plasmidcanvas](https://github.com/th0mr/plasmidcanvas) | a matplotlib circular map; features added by hand; moves overlapping features inward | No record input, no label collision handling (read), no linear view; 1.0.0 in 2024, 4 stars; PyPI only. Imports on 3.13 (measured) |
| [SpliceCraft](https://github.com/Binomica-Labs/SpliceCraft) | a terminal plasmid workbench that exports PNG or SVG | An application, not a library; pins `biopython<1.87`, which conflicts with this repository's `>=1.88` ([PyPI](https://pypi.org/project/splicecraft/)) |
| [plasmidviewer](https://github.com/ponnhide/plasmidviewer) | a circular GenBank map | GPL-3.0; not on PyPI; last push 2021 |
| [pyCircos](https://github.com/ponnhide/pyCircos) | circos plots | GPL-3.0; last release 2022; superseded for this use by pyCirclize |
| [LoVis4u](https://github.com/art-egorov/lovis4u) | linear comparisons of genomic loci | Comparative, not a plasmid map; [bioconda](https://anaconda.org/bioconda/lovis4u) lists its licence as BSD-3-Clause AND GPL-3.0 AND WTFPL |
| [SeqViz](https://github.com/Lattice-Automation/seqviz), Open Vector Editor | JavaScript viewers with a map and a sequence view | Not Python; issue #78's to weigh as an HTML route |

## 10. What #80 should prototype

**Carry two routes.**

1. **pyCirclize for the circle and pyGenomeViz for the line and the region.** It is the lightest
   maintained stack that reaches both views: matplotlib plus pandas, MIT, on conda-forge and PyPI,
   written by one author to one style. The prototype has to supply what they do not: labels that
   never overlap on the circle (#83), splitting a span at the origin, sites and primers as labelled
   layers, and show/hide in HTML (#78). The sequence view could be drawn on the same matplotlib
   figure, so one plotting library serves both views.
2. **gbdraw, as the benchmark for the picture.** It is the closest drawing out of the box, and the
   prototype should put its crowded map beside the matplotlib route's. Adopting it would mean a
   bioconda-only, pixi-only dependency whose API is still changing, CairoSVG for PNG and PDF, and
   the same missing label guarantee and show/hide.

**Drop the rest as dependencies.** DNA Features Viewer has stalled and cannot draw a joined or
origin-crossing feature; Biotite drops the outside labels the reference needs; GenomeDiagram has
no label placement; QUEEN does not install as released; pLannotate is GPL. DNA Features Viewer's
inline-or-outside rule and translation layout, and QUEEN's two-strand rows, are worth reading
before the sequence view is designed.

## Open gaps

- **gbdraw's overlaps were read off a rendering, not counted.** The count for the matplotlib
  packages compared boxes in display coordinates; gbdraw's SVG was not measured that way.
- **QUEEN's label claim was not tested** on the crowded map.
- **pLannotate's map was read, not drawn.** Its environment pulls BLAST and DIAMOND, and its
  input is a hit table rather than a record, so no like-for-like drawing was possible.

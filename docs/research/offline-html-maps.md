---
search:
  exclude: true
---

# Offline HTML maps: the ways to write one file, and what each costs

Research note for issue #78. Everything below was read or measured on **2026-09-16**.

It compares the ways the planned `liulab_mbio.plot` module could write a sequence record as one
HTML file that opens offline: a map, the sequence view below it, details on hover, layers and
feature types shown or hidden, and two switches at the top right, as in SnapGene Viewer — circular
or linear, and the sequence view shown or hidden. It also asks what the same drawing gives as PNG
and PDF. Label placement (#83), colours (#79) and the Python map packages (#77) appear only where a
route decides them.

## 1. The verdict, first

| Route | One offline HTML file | Hover, show/hide | The two switches | Sequence view | PNG and PDF | Carry to #80 |
| --- | --- | --- | --- | --- | --- | --- |
| **SVG written by hand, with our script** | **100 KB**: both maps, the sequence view, every switch | ours, 0.8 KB of script | ours: both maps in the page, 0.4 KB of script | ours; a click on a feature in either map highlights its bases | the same SVG through CairoSVG or vl-convert, no browser | **Yes — the route** |
| SeqViz, inlined | 594 KB | hover gives the name only; show/hide ours, 0.8 KB | no linear map; the sequence view goes below the map only as a second viewer, which the map no longer scrolls | free, closest to SnapGene's | none without a browser | **As the reference** for the sequence view, not a dependency |
| matplotlib SVG, with our script | 69 KB, map only | ours, as the first row | ours, as the first row | ours, in matplotlib | native | No |
| Altair (Vega-Lite) and vl-convert | 906 KB, map only | free: tooltips, checkboxes, legend | a second chart, switched by a parameter | ours, awkward in Vega-Lite | vl-convert, no browser | No |
| Bokeh, inline | 1.5 MB, map only | free: hover, legend hide, zoom | a second figure, in tabs or by a callback | ours | needs Playwright or Selenium and a browser; no PDF | No |
| plotly, inline | 4.4 MB, map only | free: hover, legend toggle, zoom | a second figure, or buttons that switch traces | ours | kaleido needs Chrome installed | No |
| CGView.js with d3, inlined | 528 KB, map only | free popovers and zoom; show/hide through its API | circular/linear is one call; no sequence view to show | none | from a browser only | No |
| Open Vector Editor, inlined | 4.9 MB, fonts still outside | free: tooltips with type and span, a visibility menu | free panels: circular, linear, sequence | free, inside an editor | none without a browser | No |

Sizes are for one plasmid, pUC19-GFP (section 2). *Ours* means code this package would write;
*free* means the tool does it.

**Only the hand-written SVG gives one drawing for all three outputs.** Python writes the SVG. The
page wraps it with a stylesheet and a script, as `protocol/render.py` already wraps its simulated
gel, and CairoSVG or vl-convert turns the same SVG into PNG and PDF without a browser. Whatever rule
keeps labels apart (#83) therefore runs once, in Python, and the PNG shows what the page shows. The
page holding a circular map, a linear map, the sequence view below and every switch was 99,857
bytes, and every behaviour asked for — hover, show/hide, both switches, and a click on a feature
highlighting its bases — worked offline on 1.9 KB of script. The cost is building: both maps, the
sequence view, label placement, and zoom if #80 wants it.

**SeqViz is the sequence view to match, not to ship.** It gives, free, most of what SnapGene's
Sequence tab shows, and a click on its map moves its sequence view. But it draws no primers on its
circular map, its hover shows only a name, and it resolved 12 of the 30 enzyme names asked for and
dropped the rest without a message. It has no linear map, and it sets the sequence view below the
map only as a second viewer, whose scroll the map then no longer drives. It keeps 7 of 52 rows in
the page at once, and nothing turns it into PNG or PDF without a browser. Shipping it would still
leave a static sequence view to build: two drawings of one view.

**The plotting libraries buy hover and zoom with megabytes, and none keeps labels apart.** plotly,
Bokeh and Altair drew the enzyme names over each other exactly as the hand-written SVG did, from the
same coordinates. From Python, neither plotly nor Bokeh writes a PNG without a browser: kaleido stopped with
"Kaleido requires Google Chrome to be installed" once Chrome was hidden, and Bokeh documents
Playwright or Selenium. A lab compute node cannot count on a browser.

**The JavaScript viewers each break a settled rule.** Open Vector Editor draws its own enzyme list
and amino-acid rows the file does not hold, SeqViz searches the sequence for sites itself, and
CGView.js has no sequence view. All three need a browser for any static file.

## 2. How this was measured

- **One record**, pUC19-GFP: 3,347 bp, circular, 12 features, 4 primer binding sites, and 30
  enzymes that cut once — the shipped enzymes that cut once, topped up from Biopython's list of
  commercially available enzymes to crowd the map.
- **One geometry.** Every route drew the same shapes from the same coordinates — feature arrows on a
  ring, ticks with `name (position)` outside, primer arrows — so only the output layer differed.
- **Throwaway environments**, one per route, from conda-forge and bioconda: Python 3.13 on
  osx-arm64. A weight is the disk space a route adds over Python and NumPy alone.
- **Offline in a browser.** Each page was opened in Chrome with the network emulated offline. Hover,
  switches and links were driven by script and read back from the page, and every request the page
  made was listed.
- **Measured or read.** Every byte count, time, weight and behaviour described as seen was
  measured. Requirements, licences and behaviour quoted from docs, package metadata or source were
  read, with a link. Times are one run on one laptop, and only rank the routes.

## 3. SVG written by hand, with a little script

`protocol/render.py` already writes one page with its stylesheet and script inline, 16.4 KB
together, and draws the simulated gel as inline SVG. This route does the same for a map.

| Part (measured) | Bytes |
| --- | --- |
| circular map SVG | 25,981 |
| linear map SVG | 11,675 |
| sequence view SVG, 56 rows of 60 bp | 58,631 |
| map page with hover and show/hide | 28,115; the script 805 |
| map and sequence view side by side, linked | 86,999; the link script 677 more |
| both maps, the sequence view below, both switches | 99,857; the switch script 448 more |

**Interaction, all ours, all measured working offline.** Each drawn group carries a `data-info`
attribute that a hover card shows. A checkbox per layer — features, enzyme sites, primers — and per
feature type hides the matching groups. A click toggles a highlight. A click on a feature in either
map highlighted its span in every sequence row and scrolled there.

**The two switches.** Python writes both maps into the page. Radio buttons at the top right show
one; a checkbox hides the sequence view below. The script was 448 bytes, and the page grew by the
second map, 11.7 KB here. Zoom was not built: pyGenomeViz inlines Panzoom (MIT, 10,125 bytes for
4.6.2, [npm](https://www.npmjs.com/package/@panzoom/panzoom)) for exactly this, untried here.

**The sequence view is ours.** The prototype drew both strands, feature bars, primers, enzyme names
above their cut and the GFP translation; not yet a tick ruler, a three-letter translation or the
overview strip. Each row of bases is one `text` element whose `textLength` pins it to the column
grid whatever font draws it, and CairoSVG and vl-convert both honoured that: in both PNGs the bases
end where the feature bars end. The view also split into 7 A4-shaped pages for a PDF. It grows by
about 17.5 bytes per base (58,631 / 3,347), so a 10 kb plasmid would be about 175 KB — computed,
not measured.

**Labels: nothing free.** The first rendering printed enzyme names over each other, in the page and
the PNG alike. Here the rule is Python writing coordinates, so the page, PNG and PDF agree, provided
text is the same width in each (section 10).

**Python drives it** with strings, as `render.py` does. The page needs no dependency; PNG and PDF
need one converter (section 9).

## 4. matplotlib SVG, with hover added

matplotlib writes the SVG, `set_gid` names each artist's group, and matplotlib's own gallery adds
tooltips and click-to-hide to such a file with a script
([SVG tooltip](https://matplotlib.org/stable/gallery/user_interfaces/svg_tooltip_sgskip.html),
[SVG histogram](https://matplotlib.org/stable/gallery/user_interfaces/svg_histogram_sgskip.html)).

Measured: drawing and saving PNG, PDF and SVG took 0.3 to 0.4 s. The SVG was 59,318 bytes with text
kept as text (`svg.fonttype: none`) and 115,459 with each glyph drawn as a path, the default. The
page with section 3's script was 68,806 bytes; the PNG 194 KB, the PDF 39 KB.

It adds two things over section 3: matplotlib measures text itself, and PNG and PDF need no
converter. It costs 134 MB and 51 packages, and the ids reach the SVG through a regular expression
over matplotlib's output. With `none`, matplotlib will "assume fonts are installed on the machine
where the SVG will be viewed"
([matplotlibrc](https://matplotlib.org/stable/users/explain/customizing.html)), so the page draws
text at a width matplotlib did not measure (section 10); drawing glyphs as paths keeps the width
but doubles the file. matplotlib has no label placement of its own; adjustText (MIT,
[conda-forge](https://anaconda.org/conda-forge/adjusttext)) is the usual add-on, for #83 to judge.
The switches cost what they cost in section 3, and the sequence view would be drawn in matplotlib.

## 5. plotly, with plotly.js inlined

`write_html(include_plotlyjs=True)` puts "the plotly.js source code (~3MB)" in the file
([write_html](https://plotly.com/python-api-reference/generated/plotly.io.write_html.html)).
Measured: 4,365,275 bytes. The same page with a CDN script tag is 72,078, so the figure itself is
about 72 KB. plotly 7.0.0 inlined plotly.js 4.0.0, and the page made no request offline.

**Free:** hover on filled shapes (measured: AmpR's arrow showed name, type, span and strand); zoom,
pan and a PNG download button in the page's toolbar (measured); and a legend whose entries hide
their traces, since "users may show or hide traces by clicking or double-clicking on their
associated legend item" ([legend](https://plotly.com/python/legend/)). `write_html` inlines only the
full bundle. The smaller `plotly-basic.min.js` is 1,191,749 bytes (4.1.1,
[npm](https://www.npmjs.com/package/plotly.js-basic-dist-min)), but using it means writing the page
ourselves.

**The switches, read:** circular or linear is a second figure in the page, about 72 KB more and a
few lines of script, or `updatemenus` buttons whose `"update"` method sets each trace's `visible`
([custom buttons](https://plotly.com/python/custom-buttons/)). The sequence view would be another
figure drawn from text traces, shown or hidden the same way; not attempted.

**Static:** "Starting with Kaleido version 1.0.0, Chrome is not included in the package itself;
instead, Kaleido looks for a compatible version of Chrome (or Chromium) already installed on the
machine" ([static image export](https://plotly.com/python/static-image-export/)). Measured: PNG
and PDF in 2.9 s with Chrome installed. With the browser path pointed at nothing, `write_image`
raised "Kaleido requires Google Chrome to be installed."

**Labels:** none; the enzyme names overlapped.

## 6. Bokeh, with inline resources

`INLINE` resources "provide minified BokehJS from library static directory"
([resources](https://docs.bokeh.org/en/latest/docs/reference/resources.html)). Measured: 1,462,596
bytes, and 74,287 with CDN resources instead; BokehJS 3.10.0; no request offline.

**Free:** hover tooltips, pan, wheel zoom and save (measured in the page's toolbar), and a legend
set to `click_policy = "hide"` hides a feature type's glyphs when clicked
([legends](https://docs.bokeh.org/en/latest/docs/user_guide/interaction/legends.html)). Labels
overlapped, and those past the plot frame were cut off (measured).

**The switches, read:** "Tab panes allow multiple plots or layouts to be shown in selectable tabs",
and a `CustomJS` callback on a widget can change another model
([widgets](https://docs.bokeh.org/en/latest/docs/user_guide/interaction/widgets.html)). Either way
the linear map and the sequence view are further figures.

**Static:** "you need a headless browser backend. Bokeh supports two options: Playwright and
Selenium", for PNG and SVG
([export](https://docs.bokeh.org/en/latest/docs/user_guide/output/export.html)); the page names no
PDF export.

## 7. Altair and vl-convert

One Vega-Lite spec gives the page and, through vl-convert, PNG and PDF. vl-convert "embeds the
official Vega and Vega-Lite JavaScript libraries, so conversions do not require a browser or
Node.js" ([vl-convert](https://github.com/vega/vl-convert)).

Measured: the inline page 906,012 bytes, no request offline; the spec 7.7 KB; PNG and PDF 2.3 s.
**Free:** tooltips, a checkbox per layer bound to a parameter, a legend click selecting a feature
type. **The switches:** the linear map is a second chart, shown by a parameter; the sequence view
likewise.

**Labels.** Vega's `label` transform "positions text marks so that they do not overlap with each
other or with other marks in the chart", trying anchor positions round each mark
([label transform](https://vega.github.io/vega/docs/transforms/label/)), and a label that fits
nowhere gets opacity 0
([source](https://github.com/vega/vega/blob/main/packages/vega-label/src/LabelLayout.js)). It is a
Vega transform, so the spec would be written in Vega rather than Vega-Lite, and it sets labels
beside marks rather than in columns with leader lines. Because vl-convert runs Vega itself, a PNG
would place them exactly as the page does.

**Sequence view:** thousands of positioned letters in wrapped rows is not a chart Vega-Lite
describes naturally; not attempted.

## 8. JavaScript plasmid viewers, inlined

Each ships as a minified bundle. Inlining one means keeping the bundle as package data, rebuilt by a
script and sourced in a note; writing the record into the page as JSON; and calling the viewer. The
page must escape `</script` inside the bundle, and must carry the bundle's licence notice. None of
them needs anything installed.

### 8.1 SeqViz

`seqviz.min.js` 3.10.24 is 587,072 bytes, MIT ([npm](https://www.npmjs.com/package/seqviz)). Its
header points to a separate 1,605-byte notice file that also covers the React it bundles, so the
page would carry that file too. Python writes `seqviz.Viewer(element, props).render()` with `seq`,
`annotations` and `primers` — "0-based start (inclusive) and end (exclusive)", this package's own
rule — plus `enzymes`, `translations` and `viewer`, one of `"linear"`, `"circular"`, `"both"` or
`"both_flip"` ([README](https://github.com/Lattice-Automation/seqviz)). A span across the origin is
`end < start`
([Annotations.tsx](https://github.com/Lattice-Automation/seqviz/blob/develop/src/Circular/Annotations.tsx)).

Measured:

- **Size.** 592,262 bytes, of which the record 4,967. With our show/hide script, which calls
  `setState`, 594,238; the script was 794 bytes.
- **Offline.** By default the page asked Google Fonts for Roboto Mono.
  `disableExternalFonts: true` removed the request.
- **Hover.** An SVG `title` with the feature's name only.
- **Link.** In one viewer showing both, a click on the map moved the sequence view's cursor and
  scrolled to it.
- **Sequence view.** Both strands, a position ruler under each row, features and primers as arrows,
  enzyme names above boxed recognition sites, a one-letter translation. No overview strip.
- **Enzymes.** Names resolve against its own list of 236. It resolved 12 of our 30 and dropped 18
  without a message. A custom enzyme is `{name, rseq, fcut, rcut}`
  ([README](https://github.com/Lattice-Automation/seqviz)), so ours could be passed, but SeqViz
  still finds the sites itself.
- **Primers** appear only in the sequence view; the circular viewer's source has no primer layer
  ([Circular.tsx](https://github.com/Lattice-Automation/seqviz/blob/develop/src/Circular/Circular.tsx)).
- **The switches.** `viewer: "both"` sets map and sequence view side by side; custom placement uses
  `refs` and `children`, a React-only option
  ([README](https://github.com/Lattice-Automation/seqviz)). Two viewers stacked, map above, took
  1,247 bytes of script. The checkbox hid the lower viewer. A click on the map reached the sequence
  view through its `selection` prop, which placed the cursor but did not scroll to it. The
  `"linear"` choice is not a linear map: at `zoom: {linear: 0}` it is the sequence view wrapped at
  about 580 bp per row, still scrolling, without enzyme names.
- **Labels.** Colliding labels on the circle are grouped under the first name with `+N`, spread a
  little when a group of up to four fits, and opened on hover
  ([Labels.tsx](https://github.com/Lattice-Automation/seqviz/blob/develop/src/Circular/Labels.tsx)).
  The page showed `MCS,+1`. A static capture would show the group closed.
- **Static.** None without a browser. The sequence view "Renders only the seqBlocks that are
  visible"
  ([InfiniteScroll.tsx](https://github.com/Lattice-Automation/seqviz/blob/develop/src/Linear/InfiniteScroll.tsx)):
  7 of 52 rows were in the page at once, so printing the page does not print the sequence.
  `renderToString()` exists, but needs a JavaScript runtime.

### 8.2 Open Vector Editor

The `openVectorEditor` repository is archived and says "THIS REPO HAS MOVED" to `TeselaGen/tg-oss`
([old repository](https://github.com/TeselaGen/openVectorEditor)). Its npm package stopped at
18.3.6 in July 2023; the editor is now `@teselagen/ove`, 0.8.42 in April 2026
([npm](https://www.npmjs.com/package/@teselagen/ove)). The repository is MIT
([LICENSE](https://github.com/TeselaGen/tg-oss/blob/main/LICENSE)). Python writes
`createVectorEditor(node, {readOnly: true})` and `updateEditor({sequenceData, panelsShown})`, with
annotations "0-based inclusive"
([README](https://github.com/TeselaGen/tg-oss/blob/main/packages/ove/README.md)).

Measured on 18.3.6:

- **Size.** 4,858,400 bytes, and the page still asked for three font files beside it, which one
  file would have to inline. The current release's bundle and stylesheet are 7,398,558 and
  1,962,246 bytes ([jsDelivr](https://www.jsdelivr.com/package/npm/@teselagen/ove)) — read, not
  built.
- **Free.** Tooltips such as "Primer - pUC19 backbone reverse - Start: 378 End: 395"; a toolbar
  with a visibility menu; zoom and rotate sliders. The switches are its panels — `circular`, `rail`
  (a linear map) and `sequence` — chosen through `panelsShown`; side by side here, and not tried
  stacked.
- **Against the settled rules.** It drew its own enzymes (EcoO109I, SspI, ScaI and more that were not
  asked for) and amino-acid rows across the sequence that the file does not annotate. Both would
  need turning off.
- **Labels.** Spread by angle, the rest grouped as `+3,M13 fwd` and opened on hover
  ([Labels](https://github.com/TeselaGen/tg-oss/blob/main/packages/ove/src/CircularView/Labels/index.js)).
- **Notices.** The bundle holds code under Apache-2.0 alongside React and other MIT notices, with no
  collected notice file, and the npm package declares no licence.

### 8.3 CGView.js

CGView.js 1.8.2 is Apache-2.0 and needs d3 (ISC) loaded first
([npm](https://www.npmjs.com/package/cgview)). Python writes its JSON — features with start and
stop, tracks, a legend — and calls `cgv.io.loadJSON(json)` then `cgv.draw()`.

Measured: 527,929 bytes — the bundle 222,830, d3 279,706, the stylesheet 14,980 — and no request
offline.

- **Free.** Popovers on hover, giving type, name, length and track (AmpR: CDS, 861 bp), smooth
  zoom, and the linear map from one call, `settings.update({format: "linear"})`.
- **Not free.** Show/hide, through each feature's `visible`
  ([Feature](https://js.cgview.ca/api/Feature.html)). The first draw sat off-centre until a
  `resize`.
- **Sequence view.** None; bases appear on the ring only at deep zoom, so the second switch has
  nothing to show.
- **Labels.** Leader lines to labels outside; `priorityMax` is "the number of priority labels that
  will be drawn for sure. If they overlap the label will be moved until they no longer overlap"
  ([Annotation](https://js.cgview.ca/api/Annotation.html)).
- **Static.** `downloadImage` and `getSVG` give "the currently visible map", in the browser
  ([IO](https://js.cgview.ca/api/IO.html)).

### 8.4 Others

angularplasmid (ISC, 20,020 bytes) was last published in February 2015 and depends on AngularJS
([npm](https://www.npmjs.com/package/angularplasmid)). pyGenomeViz writes one offline HTML file
from matplotlib, measured at 712,285 bytes, but "is a genome visualization python package for
comparative genomics" drawing linear tracks ([docs](https://moshi4.github.io/pyGenomeViz/)); #77
covers it.

## 9. From one SVG to PNG and PDF

| Converter (measured) | PNG | PDF | Circular map, PNG and PDF | Sequence view pages | Needs | Weight | Licence (read) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CairoSVG 2.9.1 | yes | yes | 0.37 s | 7 page PDFs and a PNG, 0.37 s | cairo, from conda-forge | +75 MB, 38 packages | LGPL-3.0 |
| vl-convert-python 1.9.0 | yes | yes | 1.18 s | 7 page PDFs, one tall PDF and a PNG, 14.4 s | nothing else | +70 MB, 1 package | BSD-3-Clause |
| resvg-py 0.5.0 | yes | no | PNG only, 0.33 s | — | nothing else | +2 MB, 1 package | MIT |
| kaleido 1.3.0, plotly only | yes | yes | 2.9 s | — | Chrome installed | +47 MB, and a browser | MIT |

CairoSVG "can export SVG files to PDF, EPS, PS, and PNG files"
([CairoSVG](https://github.com/Kozea/CairoSVG)); LGPL-3.0 covers CairoSVG itself, which the
package would import, not copy. resvg-py exposed one function, `svg_to_bytes`, which returns a PNG.
The newest vl-convert-python builds on conda-forge are 2.0.0 release candidates, and 1.9.0 is the
newest stable ([conda-forge](https://anaconda.org/conda-forge/vl-convert-python)). CairoSVG and
vl-convert both honoured `textLength`.

## 10. Text width, one fact for #83

A rule that keeps labels apart must know how wide each label is. Where Python lays labels out and a
browser draws them, the two widths must agree.

Measured on this map's 42 labels: Chrome's default sans-serif, Helvetica on macOS, drew them between
15% narrower and 5% wider than matplotlib measured them in its bundled DejaVu Sans. A monospaced font
removes the question: Menlo in Chrome and DejaVu Sans Mono in matplotlib both advance 0.602 of the
font size per character, so a label's width is its length times one number. For proportional text,
the page and the converter must draw with the same font file, which #83 decides how to ensure.

Only the hand-written SVG and Vega draw the page and the static file from one layout. plotly, Bokeh
and the JavaScript viewers leave text to the browser.

## 11. Licences for JavaScript inside an MIT package's output

| JavaScript | Licence (read) | What putting it in a page asks |
| --- | --- | --- |
| our own script | MIT, this package | nothing |
| plotly.js | MIT | written by plotly; the bundle keeps its licence comments (measured) |
| BokehJS | BSD-3-Clause ([npm](https://www.npmjs.com/package/@bokeh/bokehjs)) | written by Bokeh; the bundle keeps its copyright line (measured) |
| Vega, Vega-Lite | BSD-3-Clause | written by Altair |
| SeqViz, and the React it bundles | MIT | inline its 1,605-byte notice file |
| Open Vector Editor | MIT repository; the npm package declares none; the bundle also holds Apache-2.0 code | collect the notices ourselves |
| CGView.js; d3 | Apache-2.0; ISC ([npm](https://www.npmjs.com/package/d3)) | a copy of the Apache licence goes with the page |

The MIT licence asks that "the above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software" ([MIT](https://opensource.org/license/mit)).
Apache-2.0 asks that recipients get "a copy of this License"
([Apache-2.0, section 4](https://www.apache.org/licenses/LICENSE-2.0)). A vendored bundle is package
data, so under this repository's rule it is rebuilt by a script and sourced in a note; a route that
writes only our own script has none.

## 12. Dependencies

| Route | Adds to the environment (measured on osx-arm64) | conda-forge, linux-64 and osx-arm64 (read) |
| --- | --- | --- |
| hand-written SVG, page only | nothing | — |
| hand-written SVG, with PNG and PDF | CairoSVG +75 MB, or vl-convert-python +70 MB | CairoSVG noarch; vl-convert-python built for both |
| matplotlib | matplotlib-base +134 MB, 51 packages | built for both |
| plotly | +45 MB, 3 packages; kaleido +47 MB and Chrome | noarch |
| Bokeh | +61 MB, 27 packages; export adds Playwright or Selenium and a browser | noarch |
| Altair | altair and vl-convert-python, +88 MB | altair noarch; vl-convert-python built for both |
| a JavaScript viewer | nothing; the bundle is package data | — |

Platforms and licences are from anaconda.org's package pages:
[plotly](https://anaconda.org/conda-forge/plotly),
[python-kaleido](https://anaconda.org/conda-forge/python-kaleido),
[bokeh](https://anaconda.org/conda-forge/bokeh),
[matplotlib-base](https://anaconda.org/conda-forge/matplotlib-base),
[cairosvg](https://anaconda.org/conda-forge/cairosvg),
[resvg-py](https://anaconda.org/conda-forge/resvg-py),
[vl-convert-python](https://anaconda.org/conda-forge/vl-convert-python),
[altair](https://anaconda.org/conda-forge/altair) and
[playwright-python](https://anaconda.org/conda-forge/playwright-python).

## 13. What #80 should take from this

- **Prototype one route:** both maps and the sequence view as SVG from Python, in one page with the
  switches of section 3, on a record crowded enough to test #83's rule.
- **Pick one converter:** CairoSVG, faster over many pages but LGPL and pulling in cairo; or
  vl-convert, one BSD package but slower over the sequence view's pages. Both passed on this record.
- **Hold the sequence view against SeqViz:** its ruler, strand layout and boxed cut sites are the
  nearest open example of SnapGene's Sequence tab.
- **Still open:** zoom on this route, untried; and how long a record the page must carry. At about
  17.5 bytes a base every row is in the page, which prints whole but grows with the record; drawing
  rows on demand, as SeqViz does, would split the page from the PNG.

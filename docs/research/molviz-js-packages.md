---
search:
  exclude: true
---

# Which JavaScript packages liulab-molviz can build on

Research note for issue #575, a sub-issue of map #572. liulab-molviz will draw in the browser from
data mbio writes: a circular and a linear plasmid map, a sequence view with both strands and
translations, Sanger trace tracks, a simulated gel, a plate grid, and later 3D protein structures.
A command line must also write SVG, PNG and PDF with no browser. The reader promises it inherits:
no label overlaps another and a hidden label is counted; bases are real text that find reaches and
copy takes; one file opens offline. The map also decided that the viewer draws only what it is
given: mbio finds cut sites, translates and places features.

Every claim below was **read** (from a package's source, docs, LICENSE or npm metadata) or
**measured** (installed, run, downloaded and sized, or counted through the GitHub API) on
**2026-10-17**, unless its line gives another date. Sizes are the package's own browser build,
minified as shipped, then `gzip -9`. Commit counts are commits on the default branch since
2025-10-17, from `gh api .../commits?since=2025-10-17T00:00:00Z --paginate`. Claims from memory are
marked **unverified**.

## 1. The verdict, first

**SeqViz is a look to copy, not a package to build on** (measured and read, 2026-10-17). It does
render to an SVG string under Node, but only through a test hook, and four of its traits break a
need:

- **The sequence view is virtualised.** It mounts five rows at a time. On the server the linear
  view returned the first 450 of 3,347 bases and nothing else; in a browser, rows off screen are not
  in the page, so find cannot reach them (read: `InfiniteScroll`, `visibleBlocks` starts at
  `Math.min(5, …)`).
- **Bases are text but not selectable.** Each base is a `<tspan>` inside one `<text>`, styled
  `user-select: none`; copy goes through SeqViz's own Ctrl+C handler (measured in the server
  output; read: `copyEvent`).
- **It computes what mbio computes.** It takes enzyme names or `{rseq, fcut, rcut}` and digests the
  sequence itself, and takes translation ranges and translates them itself. There is no prop for a
  cut site at a given position (read: `SeqVizProps`).
- **Its element kinds are fixed**: annotations, primers (linear viewer only), translations, enzymes,
  highlights, search. `children` and `refs` rearrange its two viewers; they add no track kind
  (read: `SeqVizProps`).

It also needs React (bundled in its 188 KB build) and fetches Roboto Mono from Google Fonts unless
`disableExternalFonts` is set (read).

**The strongest alternative for the prototype is no viewer package: a framework-free SVG core of
molviz's own on d3-path, d3-shape and d3-scale, with rbush for collisions and opentype.js for text
width** (measured: all five run under Node with no DOM, section 5). It is the only route that meets
every need at once, about 84 KB gzipped together (sum of the dist files measured, d3-scale's own
dependencies excluded). It can draw SeqViz's look, which is flat SVG. Among packaged viewers
**CGView.js** is the runner-up for maps: Apache-2.0, SVG export, a plugin system, and the nearest
label rule to ours: priority labels pushed outward until free, then every other label drawn only
where it overlaps nothing. It drops the rest without counting them, though the count is one
subtraction away. It draws to a canvas: its labels and bases are pixels on screen, its SVG export
needs a DOM (`d3.select('body')`) and svgcanvas, and it has no sequence view. Its label code is
worth reading, and Apache-2.0 lets molviz port it with the notice kept.

| View | Build on | Read or copy from | Reject |
| --- | --- | --- | --- |
| Circular and linear map | own core on d3-path, d3-shape, d3-scale and rbush | SeqViz (look, `name,+N` groups), CGView.js (label islands), @carabennemsi/plasmid (pure SVG function) | OVE (React, Redux, 2.2 MB gzip), nitro-bio (React) |
| Sequence view | own core; one `<text>` per strand per row, all rows in the page | SeqViz (layout of strands, primers, translation boxes) | SeqViz's virtualised view |
| Sanger trace | own core | teal (approach only: GPL-3.0), sanger-viewer (MIT, pure `exportSvg`) | OVE's trace (canvas, React) |
| Simulated gel | own core | hatchlings `GelViewer` (Svelte) | nothing else exists on npm |
| Plate grid | own core | hatchlings `PlateHeatmap` (Svelte) | plate-map (AGPL-3.0), @nitro-bio/platemap (GPL-3.0) |
| 3D structure, later | **Mol\*** through `PluginContext` without its React UI, fed MolViewSpec | 3Dmol.js as the light fallback | iCn3D (licence, network), NGL (maintenance only), LiteMol (dead) |
| Label placement | own 1D de-overlap on rbush that returns the hidden labels | d3fc-label-layout `removeOverlaps` (marks `hidden = true`) | labella (no hidden output, unmaintained) |
| Text width with no browser | opentype.js, with the same font file inlined in the page | fontkit (same widths) | canvas `measureText` (needs a canvas) |

## 2. The decisive check: SeqViz under Node

### 2.1 What was run

In a fresh empty directory outside the repo (measured, 2026-10-17; Node 26.10.0, pnpm 12.4.1):

```sh
echo '{"name":"seqviz-env","private":true}' > package.json
pixi exec -s nodejs -s pnpm -- pnpm add seqviz react react-dom
# → seqviz 3.10.25, react 19.3.1, react-dom 19.3.1
```

The record was `docs/examples/pUC19-GFP/product.dna` read with Biopython's SnapGene reader: 3,347
bp, 21 features, 11 of them `primer_bind`, passed as SeqViz `annotations` and `primers`, with
twelve enzymes (`EcoRI`, `SacI`, `KpnI`, `SmaI`, `BamHI`, `XbaI`, `SalI`, `PstI`, `SphI`,
`HindIII`, `NcoI`, `BsaI`) and the GFP span as a translation. The script lived in a separate
directory and loaded the packages by path:

```js
// node seqviz-ssr.cjs ENV_DIR RECORD_JSON OUT_DIR
const req = createRequire(path.join(envDir, "package.json"));
const React = req("react");
const { renderToStaticMarkup } = req("react-dom/server");
const seqviz = req("seqviz");
const base = {
  name: rec.name, seq: rec.seq, annotations, primers, enzymes,
  translations: [{ start: gfp.start, end: gfp.end, direction: gfp.direction }],
  disableExternalFonts: true,
  style: { width: "1000px", height: "800px" },
};
if (process.env.TEST_SIZE) base.testSize = { width: 1000, height: Number(process.env.TEST_SIZE) };
for (const viewer of ["circular", "linear", "both"]) {
  renderToStaticMarkup(React.createElement(seqviz.SeqViz, { ...base, viewer }));  // route 1
  seqviz.Viewer(undefined, { ...base, viewer }).renderToString();                 // route 2
}
```

Each output's `<svg>` was pulled out and rasterised with `vl_convert.svg_to_png` to look at it.

### 2.2 What came out

| Route | Size given | Result |
| --- | --- | --- |
| `Viewer(…).renderToString()`, the README's "server-side rendering" call | any | **Throws** `ReferenceError: document is not defined`, for all three viewers |
| `renderToStaticMarkup(<SeqViz>)` | `style` only | Empty containers, no `<svg>`: width and height come from a resize detector, which reports 0 on the server |
| `renderToStaticMarkup(<SeqViz>)` with `testSize` 1000 × 800 | test hook | **Circular: one complete `<svg>`**, 24,774 bytes, 24 `<text>`, 51 `<path>`: ring, ticks, feature arrows, labels on the arc through `<textPath>` (lacZα, ori, AmpR, GFP), outside labels with leader lines (AmpR promoter, ATGA junction, NcoI, BsaI ×2, lac promoter, lac operator) and one group `MCS,+1`. No primers: the circular viewer does not draw them |
| same | test hook | **Linear: five `<svg>` rows**, both strands, primers above and below, the GFP translation as boxed residues, annotations, index — **bases 1 to 450 only** |
| same, `testSize` height 60,000 | test hook | Still five rows: the window is fixed at five blocks, not sized to the height |

`testSize` is documented only in `SeqViewerContainer.d.ts`: "a forced height/width that overwrites
anything from sizeMe. For testing" (read). `SeqViz` passes every prop through to the container, so
it reaches it (read). The layout itself needs no DOM: label widths are estimated at a fixed 7.2 px
per character (`CHAR_WIDTH`), and nothing is measured (read).

**So**: a SeqViz circular map can become SVG text on a command line through an undocumented hook;
a whole sequence view cannot.

## 3. Comparison, one row per package per view

Size is min / gzip of the package's own browser build (measured). Activity is the latest release
and its date, then releases and commits since 2025-10-17. "No browser" asks whether its layout
can emit SVG text under Node. All rows: 2026-10-17.

| Package | View | Licence | Size | Activity | Framework | No browser | Crowded labels | Bases as text | Own tracks and kinds |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SeqViz 3.10.25 | circular map | MIT | 631,230 / 187,615 B (`seqviz.min.js`, React bundled); 309,663 / 65,502 B without React (`dist/index.js`, unminified) | 2026-09-23; 9 releases, 126 commits | React 16.8 to 19 | **Yes with `testSize`** (measured) | Labels within 15 px on one side group as `name,+N`, names on hover; groups of up to 4 spread if neighbours allow; horizontal overlap not checked (read) | Only for records of 250 bp or less (`RENDER_SEQ_LENGTH_CUTOFF`, read) | No: fixed kinds |
| SeqViz 3.10.25 | linear map and sequence view | MIT | same | same | React | **First 5 rows only** (measured) | Primer and annotation names truncated with `…` to fit (measured: "Sequencing forwa…") | `<tspan>` per base, `user-select: none`, virtualised (measured, read) | No: fixed kinds |
| CGView.js 1.8.2 | circular map | Apache-2.0 | 222,830 / 59,959 B, plus d3 7.0.1 (pinned) 274,269 / 89,829 B and css 3,095 B gzip | 2026-05-04; 2 releases, 64 commits | none (d3) | No: draws to canvas; `getSVG()` needs svgcanvas and `d3.select('body')`; widths by canvas `measureText` (read) | Labels sorted by priority; the first `priorityMax` are pushed outward until their rect is free, or set in angled "islands"; each remaining label is drawn only if it overlaps nothing, and the rest are **dropped without a count** (read: `Annotation`, `LabelPlacementDefault`) | No: letters of both strands drawn on the canvas when zoomed in (read) | Yes: its own tracks, and a plugin system since 1.8.0 (read: CHANGELOG) |
| CGView.js 1.8.2 | linear map | Apache-2.0 | same | same | none | No | same | No | Yes |
| Open Vector Editor (@teselagen/ove 0.8.42) | circular map | MIT (repo; npm field missing) | 7,398,558 / 1,479,935 B (`index.umd.js`, as shipped, unminified) plus css 1,962,246 / 729,986 B | 2026-04-16; 37 releases, 84 commits under `packages/ove` | React 18, Redux, Blueprint 3 | No (unverified: not tried; built on react-measure and html2canvas) | Groups crowded labels (read 2026-09-16 in `plasmid-map-packages.md`) | n/a | Unverified |
| Open Vector Editor 0.8.42 | linear and sequence view | MIT | same | same | React, Redux | No (unverified) | unverified | unverified | unverified |
| Open Vector Editor 0.8.42 | Sanger trace | MIT | same | same | React, Blueprint | No: canvas `getContext("2d")` (read) | n/a | No | No |
| @nitro-bio/sequence-viewers 2.4.0 | circular, linear, sequence | MIT | 1,110,009 B unpacked (npm) | 2026-10-01; 6 releases; 31 stars | React 18 or 19, Radix | Not tried | Not read | Not read | Not read |
| @carabennemsi/plasmid 0.1.3 | circular map | MIT | 115,172 B unpacked (npm) | 2026-07-24; 2 releases; repo created 2026-07-11, 0 stars | none | **Yes**: `renderPlasmidSVG(record)` returns an SVG string (measured) | Narrow features without room get **no label, silently**: of the 17 features not labelled on their arc, 3 get an outside label and 14 none, with no count (measured, pUC19-GFP) | No sequence view | Features and cut sites only |
| @molbiohive/hatchlings 0.8.4 | map, sequence, trace, gel, plate | MIT | 521,179 B unpacked; ships `.svelte` source | 2026-03-31; repo created 2026-03-04, 2 stars, no commit since | Svelte 5 | Not tried (needs the Svelte compiler) | Not read | Not read | Its own fixed components |
| teal (gear-genomics) | Sanger trace with reference row | GPL-3.0 | not on npm | last commit 2026-09-03 | none (Custom Element) | Its SVG is built as strings; the element wrapper needs a DOM (read by subagent) | n/a | Unverified | Takes tracy JSON: four channels, calls, qualities, gapped reference |
| sanger-viewer (animeshkundu) | Sanger trace | MIT | not on npm | last commit 2026-09-10; 1 star | none | `exportSvg` documented as pure, no DOM (read by subagent) | n/a | Unverified | No |
| plate-map 3.0.1 (NEB) | plate grid | AGPL-3.0-only | not sized | last commit 2026-10-06 | jQuery, jquery-ui, select2, svg.js | No | n/a | n/a | Editor, not a drawer |
| @nitro-bio/platemap 2.0.3 | plate grid | GPL-3.0-or-later | not sized | not checked | React | No | n/a | n/a | n/a |
| Mol\* (molstar 5.13.1) | 3D | MIT | 5,268,685 / 1,504,319 B (`build/viewer/molstar.js`) | 2026-10-08; 20 releases, 912 commits; master is a v6 prototype | React bundled; `PluginContext` and `Canvas3D` usable without its React UI (read) | PNG or JPEG through `HeadlessPluginContext` and the native `gl` module; **no SVG** (read) | Labels are 3D primitives | Residue labels, not sequence text | Yes: MolViewSpec colours, labels, tooltips, primitives, annotation files; custom themes and representations |
| PDBe Mol\* 3.12.0 | 3D | Apache-2.0 | 6,098,173 / 1,755,677 B | 2026-04-24; 5 releases, 46 commits; pins molstar 5.8.0 | React bundled, Lit | Inherits Mol\* | — | — | PDBe helpers |
| 3Dmol.js 2.5.5 | 3D | BSD-3-Clause | 537,792 / 156,450 B (`3Dmol-min.js`, re-measured) | 2026-05-22; 2 releases, 52 commits | none | No: WebGL in a browser only (read) | — | — | Yes: `setStyle` by selection, `addLabel`, `addShape`, `addSurface` |
| NGL 2.5.0 | 3D | MIT | 1,292,474 / 356,820 B | 2026-09-03; 1 release, 5 commits | none (three.js) | No (unverified) | — | — | Selection language, Shape API |
| iCn3D 3.54.5 | 3D | US public domain, except DelPhi code: no commercial use | 2,313,310 / 567,881 B | 2026-10-14; 33 releases, 83 commits | jQuery, full app | Node scripts analyse, do not render | — | — | Own command language |
| d3-path 3.1.0, d3-shape 3.2.0, d3-scale 4.0.2 | geometry for any 2D view | ISC | 1,107 B; 7,614 B; 4,890 B gzip (d3-scale's own dependencies excluded) | last releases 2022-12 to 2023-04; last commit 2023-10-06 | none | **Yes**: arcs, lines and scales give path strings under Node (measured) | — | — | — |
| rbush 4.0.1 | label collisions | MIT | 2,340 B gzip | 2024-08; last commit 2026-04-03 | none | **Yes**, `collides()` (measured) | Building block | — | — |
| @d3fc/d3fc-label-layout 7.0.4 | label placement | MIT | — | 2024-06 | none in `src` | npm tarball fails to import: `main` points at a missing build (measured); `src` has no `document` | `removeOverlaps` hides the most-overlapped label until none overlap, marking it `hidden = true` (read, re-read) | — | — |
| labella 1.1.4 | 1D label placement | Apache-2.0 | 6,802 B gzip | last change 2022-06; last commit 2021-10-21 | none | Yes (measured) | Moves overflow to more rows; reports nothing as unplaced (inferred) | — | — |
| opentype.js 2.0.0 | text width | MIT | 67,605 B gzip | 2026-05-06; last commit 2026-05-19 | none | **Yes**: advance width of `BsaI` at 12 px = 24.012 (measured) | — | — | — |
| fontkit 2.0.4 | text width with kerning | MIT on npm | 110,540 B gzip | 2026-09-30 | none | Yes, same width (measured) | — | — | — |

Rows for teal, sanger-viewer, the plate packages, the five 3D viewers and the layout helpers come
from two subagent reports; the rows marked "re-measured" or "re-read" were checked again by hand,
and so were teal's GPL-3.0 and sanger-viewer's MIT licences.

## 4. Licences

molviz is MIT and ships one built file, so every bundled package's notice must travel inside it.

| Package | Licence | Read from | For an MIT single-file build |
| --- | --- | --- | --- |
| SeqViz | MIT | `package.json`, GitHub | Fine; its own build already carries `seqviz.min.js.LICENSE.txt` listing React |
| Open Vector Editor | MIT | tg-oss GitHub licence; `@teselagen/ove` npm has no `license` field | Fine, once the field gap is accepted |
| @nitro-bio/sequence-viewers, @carabennemsi/plasmid, @molbiohive/hatchlings, sanger-viewer | MIT | npm, GitHub | Fine |
| Mol\*, NGL, rbush, d3fc-label-layout, opentype.js | MIT | npm, LICENSE | Fine |
| fontkit | MIT on npm; GitHub detects no licence file | npm | Check the tarball's LICENSE before bundling |
| d3-path, d3-shape, d3-scale | ISC | npm | Fine |
| nightingale-sequence | ISC on npm; repo MIT | npm, GitHub | Fine |
| CGView.js, PDBe Mol\*, labella | Apache-2.0 | npm, file header | Fine; keep the licence and any NOTICE text; porting CGView's label code means keeping its copyright line |
| 3Dmol.js | BSD-3-Clause, incorporating GLmol, three.js and jQuery code | LICENSE | Fine; keep all three notices |
| teal | GPL-3.0 | GitHub | **Copying its code would make molviz GPL**: read the approach only |
| plate-map | AGPL-3.0-only | npm | Out |
| @nitro-bio/platemap | GPL-3.0-or-later | npm | Out |
| iCn3D | US Government work, public domain, except DelPhi (Columbia; no commercial use without contact) | LICENSE.md | Out |

## 5. Building blocks with no DOM, measured

The layout helpers were checked in a scratch directory with `typeof document` and
`typeof window` both `undefined` (subagent, measured 2026-10-17):

```text
d3-shape arc: M0,-100A100,100,0,0,1,100,0L90,0A90,90,0,0,0,0,-90Z
d3-shape line monotoneX: M0,0C0.333,2.5,0.667,5,1,5C1.333,5,1.667,2,2,2C2.333,2,2.667,5,3,8
d3-path: M0,0L50,0A50,50,0,1,1,-50,6.123233995736766e-15
d3-scale: 300 0,1000,2000,3000,4000,5000
labella positions: 10.0,33.0,17.0,40.0,63.0 layers: 1,1,0,0,0
rbush collides: true false
opentype.js advance 'BsaI' @12px: 24.012
fontkit advance 'BsaI' @12px: 24.012
```

`@d3fc/d3fc-label-layout` failed with `ERR_MODULE_NOT_FOUND` from its own `index.js`, and the d3fc
umbrella threw "d3fc-element depends on Custom Elements (v1)" at load.

A text width from opentype.js matches the browser only if the page draws with the same font, so
the font file is inlined in the page and read by the command line. SeqViz avoids the question by
estimating 7.2 px per character; CGView.js asks a canvas.

**The pure-function route has one working example.** `@carabennemsi/plasmid` renders the pUC19-GFP
record to a standalone SVG string under Node (measured: 12,328 bytes, 13 `<text>`). Its labels show
why that is not enough alone: lacZα, GFP, AmpR and ori sit on their arcs, three narrow features
and the NcoI site get outside labels, and the other 14 features, 10 primers among them, get none,
with no count.

## 6. Also found

| Package | What it is | Why it is not in the table |
| --- | --- | --- |
| igv.js 3.8.10 | MIT genome browser, released 2026-10-14 | Linear genome tracks, no plasmid circle, no trace view (subagent: 0 hits for "chromatogram"); not assessed further |
| Nightingale 5.11.0 (EBI) | Custom elements for protein feature tracks, d3 7.9.0 | Needs a DOM; a candidate for protein tracks later, not assessed |
| @omicsos/circular-genome-viewer 0.1.4 | React wrapper around CGView.js | Adds React to CGView |
| seqviz-plus, dash_seqviz, forks of SeqViz | wrappers | Inherit SeqViz |
| AngularPlasmid 1.0.5, @yapv/svg 0.2.2, @metrichor/plasmid 0.2.0 | plasmid drawers | Last releases 2015, 2020, 2022 |
| open-vector-editor 18.3.6 | the old OVE package | Repo archived, marked "DEPRECATED"; peers React 16 |
| LiteMol | 3D | Last release 2019-06-08; replaced by Mol\* |
| JSmol | 3D | Not checked (unverified: LGPL, Java port) |
| Trace parsers: @teselagen/bio-parsers, @shandley/abif-ts, @fishka/seqio | `.ab1` readers | molviz is handed data, so it reads no `.ab1` |

## 7. The earlier findings, checked

The claims carried into this ticket from 2026-10-16:

| Claim | Now |
| --- | --- |
| SeqViz MIT, 3.10.25, active | Right: 3.10.25 on 2026-09-23, 9 releases and 126 commits in a year, last push 2026-10-13 |
| SeqViz 188 KB min+gz with React bundled | Right: 187,615 B gzip of a 631,230 B build |
| SeqViz has a map and a sequence view with both strands and translations | Right. Primers show in the linear viewer only |
| SeqViz has select and copy | Right, through its own handler: the base text is `user-select: none` |
| SeqViz has no export | **Half wrong**: no export call, but its React component renders a circular map to SVG under Node with the `testSize` hook. Its documented `renderToString` throws |
| CGView.js Apache-2.0, label placement, SVG and PNG export | Right, but the SVG export needs svgcanvas and a DOM, and labels that do not fit are dropped uncounted |
| CGView.js has no base-level view | **Half wrong**: it draws both strands' letters on the backbone when zoomed in, as canvas pixels; it has no row-by-row sequence view or translations |
| OVE MIT, React and Redux | Right; it now lives in `TeselaGen/tg-oss` as `@teselagen/ove`, the old repo archived; last release 2026-04-16 |
| OVE about 1.5 MB JS and CSS | **Wrong**: the shipped build is 7.4 MB JS and 2.0 MB CSS, 2.2 MB gzipped together; 1.5 MB is the JS alone, gzipped |
| Mol\* for 3D (from memory) | Confirmed: MIT, 5.13.1 on 2026-10-08, 912 commits in a year |

## Open gaps

- **SeqViz in a real browser, offline.** Whether Ctrl+F finds bases is inferred from the source
  (`user-select: none`, five mounted rows), not tried. The SeqViz prototype should open its page
  with the network off and search for a base run outside the first rows.
- **SeqViz on a crowded map.** The 30-label variant of `plasmid-map-packages.md` was not drawn, so
  how often `name,+N` groups or overlaps there is unmeasured.
- **Whether SeqViz would take a documented size prop** in place of `testSize`, and an option to
  mount every row. Not asked upstream.
- **CGView.js under jsdom with svgcanvas** was not tried, so its command-line cost is read, not
  measured.
- **@molbiohive/hatchlings** was not compiled and rendered with `svelte/server`; it is the one
  package covering trace, gel and plate, but is seven months old with two stars and no commit since
  March.
- **@nitro-bio/sequence-viewers** was not installed: only its README was read.
- **Mol\* headless** was not run; whether the native `gl` module builds under the pinned Node on
  osx-arm64 and linux-64 is unknown, and the v6 reorganisation on master may move the API.
- **sage's trace code** was not read, nor whether sanger-viewer's `exportSvg` draws a reference row.
- **Bundlephobia** was not used; sizes are the packages' own builds. The own-core total is a sum of
  dist files, not a bundled build.
- **OVE's headless behaviour and label rule at 0.8.42** were not re-read.
- **JSmol, vega-label, tagmap.js** were not checked.
- **This note is not in `docs/research/index.md`**: the branch carries only this file.

## Sources

All retrieved 2026-10-17. Downloads are kept under `reference_docs/molviz-packages/`.

- SeqViz 3.10.25: npm tarball (`package.json`, `README.md`, `dist/index.js`, `dist/seqviz.min.js`,
  `dist/SeqViewerContainer.d.ts`, `dist/SeqViz.d.ts`); <https://github.com/Lattice-Automation/seqviz>
- CGView.js 1.8.2: npm tarball (`README.md`, `CHANGELOG.md`, `dist/cgview.esm.js`,
  `dist/cgview.min.js`); <https://github.com/sciguy/cgview-js>; d3 7.0.1 from
  <https://cdn.jsdelivr.net/npm/d3@7.0.1/dist/d3.min.js>
- Open Vector Editor 0.8.42: <https://cdn.jsdelivr.net/npm/@teselagen/ove@0.8.42/index.umd.js>,
  `ove.css`, the jsDelivr file listing; <https://github.com/TeselaGen/tg-oss>;
  <https://github.com/TeselaGen/openVectorEditor> (archived);
  `packages/ove/src/RowItem/Chromatograms/Chromatogram.js`
- @carabennemsi/plasmid 0.1.3: npm tarball, README; <https://github.com/cara/plasmid>
- @nitro-bio/sequence-viewers 2.4.0: npm README; <https://github.com/nitro-bio/sequence-viewers>
- @molbiohive/hatchlings 0.8.4: <https://github.com/molbiohive/hatchlings> README,
  `src/lib/index.ts`, `docs/viewers/trace-viewer.md`
- teal: <https://github.com/gear-genomics/teal>, `client/src/static/js/traceView.js`
- sanger-viewer: <https://github.com/animeshkundu/sanger-viewer>, `src/export/svg.ts`
- plate-map: <https://github.com/nebiolabs/plate-map>
- Mol\* 5.13.1: <https://github.com/molstar/molstar>, `docs/docs/plugin/instance.md` and
  `docs/docs/plugin/file-formats.md` at v5.13.0, `build/viewer/molstar.js` from jsDelivr;
  molrender <https://github.com/molstar/molrender>; MolViewSpec
  <https://github.com/molstar/mol-view-spec> (`docs/docs/annotations.md`, `primitives.md`,
  `selectors.md`); PDBe Mol\* <https://github.com/molstar/pdbe-molstar>
- 3Dmol.js 2.5.5: <https://github.com/3dmol/3Dmol.js>, `src/GLViewer.ts`, LICENSE,
  `build/3Dmol-min.js` from jsDelivr
- NGL 2.5.0: <https://github.com/nglviewer/ngl>, `dist/ngl.js` from jsDelivr
- iCn3D 3.54.5: <https://github.com/ncbi/icn3d>, LICENSE.md, `icn3dnode/README.md`
- Layout helpers: npm metadata and dist files of d3-path 3.1.0, d3-shape 3.2.0, d3-scale 4.0.2,
  rbush 4.0.1, labella 1.1.4, @d3fc/d3fc-label-layout 7.0.4 (and
  `packages/d3fc-label-layout/src/removeOverlaps.js` on d3fc master), opentype.js 2.0.0,
  fontkit 2.0.4
- Earlier reading: `docs/research/plasmid-map-packages.md` (2026-09-16) for SeqViz 3.10.24 and
  OVE label grouping

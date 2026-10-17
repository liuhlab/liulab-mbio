---
search:
  exclude: true
---

# What mbio draws today, and what each drawing needs

Research note for issue #574, on the map that founds `liulab-molviz`. It is the inventory that
viewer's first roadmap starts from: every drawing mbio writes, what each shows, what it is drawn
from, what a reader can do with it, what it promises, and every call mbio makes to get it.

It describes features and calls only. How `mbio.plot` lays a drawing out today is left out on
purpose, because the viewer starts from first principles. Where a fact below names a file, that
is where the drawing is written, not a design to port.

**Every claim was read or measured on 2026-10-17, at commit `10e8785`.** Each carries its
evidence: *read* (the code or a page, at the path given), *measured* (a script over
`docs/examples/`, byte counts exact), or *inferred*.

## 1. The drawings

Six drawings, one row each. "Example size" is the HTML a reader opens, measured.

| Drawing | What it shows | Drawn from | What mbio computes before drawing | What a reader can do | Written by | Example size |
| --- | --- | --- | --- | --- | --- | --- |
| **Map page** (`.html`) | One record. A circular record drawn whole: as a circle, and as a line. A region, or any linear record: a line only, numbered as the record is. Its features, primers, cut sites and insertions, each labelled; a count of the labels left out. The sequence view of the same stretch, both strands, up to 100,000 bases | A `SequenceRecord` (name, topology, bases, features with segments, strand, type, colour and qualifiers; primers with sequence and binding sites), plus switches, enzymes, region, highlight, insertions and tracks | Cut sites, translations, primer tails and mismatches, default colours, label text, hover text, which labels hide | Hover, switch layers and feature types, flip circle and line, zoom, pan, click to select in both views, select bases and copy them (section 4) | `mbio plot map -o x.html`; every pipeline whose protocol has a figure (`*-map.html`, one per record) | 12 pages, 186,987 to 1,053,579 B; 7,402,877 B in all |
| **Map as PNG or PDF** | The same map, with only what is switched on. The sequence view under it when asked. A PDF puts the map on page 1 and the sequence view's rows on the pages after | As the map page | As the map page | Nothing: a file to print or paste | `mbio plot map -o x.png` or `-o x.pdf`, and `Drawing.write` from Python. No pipeline writes one | 7 PDFs in `docs/examples/`, 15,356 to 33,871 B |
| **Protocol figure** | One record, or several stacked as rows (one molecule becoming the next), drawn to one span, beside the step that makes it. Every unlit item, or row, dimmed. A caption, and a citation where the figure stands for a published one | A `Figure` spec in `protocol.json`: record paths, `span`, `linear`, `sequence_view`, `enzymes`, `highlight`, `caption`, `citation` (`protocol/model.py:787-834`) | As the map page, for each row; which row answers to each highlighted name | Nothing on the drawing. "Explore the map" opens that record's map page in place, and "Open it in a new tab" | The figure library (`protocol/figures.py`) and each method's steps; rendered by `protocol/render.py` `_figure` | 16 figures across 9 protocol pages; 205,850 B of SVG in all |
| **Verification page** | One clone's verdict and its table, then the expected record drawn whole with each sequencing result as an arrow over the bases it trusts, and each disagreement coloured by its region's verdict. Then one close-up per disagreement: 30 bases either side, at base level, each result's own disagreements marked on it, and a Sanger trace's four channels under the bases | The expected record, the results, the `Verification` (checks, placements, disagreements) and each trace's channels and peaks | A stand-in record: disagreements as `misc_difference` features, insertions as insertion marks, each result as a primer whose bases carry its disagreements (a deleted base as `N`) (`verification/page.py:245-266`). Traces placed base by base | Nothing: no script on the page | `mbio sequence-verify ... --out DIR` (`verification.html`) | 1 page, 94,927 B |
| **Plate wells** | A plate's wells by row label and column number, each well filled by what sits in it, a legend of what each fill means. Past 10 kinds, no legend: every used well one colour, and the caption counts the kinds. A transfer drawn as two plates, "From" and "Into", every well it touches filled, with its moves listed under a closed toggle | A `Plate` (format, seating, labels, catalogue, holds, note, rules) or a `Transfer` and its `Stamp` (`protocol/model.py`: `Plate` at 984, `Stamp` at 1103, `Transfer` at 1152) | Row labels and well names from the format; which wells a stamp touches | Nothing on the drawing; open the list of moves | `protocol/render.py` `_plate`, `_transfer` (through `draw_plate`); iGGA and DMX protocols | 19 plate figures across 6 pages; 1,133,447 B of SVG in all |
| **Simulated gel** | A ladder lane and sample lanes, each band at its size; ladder sizes labelled; a lane legend giving each lane's band sizes, or "no band" | A `Gel`: ladder name and band sizes, lanes with label and band sizes, title (`protocol/model.py:741-784`) | Where each band runs, linear in the log of its size (`Gel.migration`) | Hover a band to see its size (an SVG `<title>`) | `protocol/render.py` `_gel`, outside `mbio.plot`; every cloning method and iGGA | 11 gels across 5 pages; 32,275 B of SVG in all |

Evidence: rows read from `plot/drawing.py`, `plot/page.py`, `protocol/render.py:2058-2253` and
`:2754-2798`, `verification/page.py`; counts and sizes measured.

## 2. What each view of a map shows

A map page holds up to three views of one record; a figure and the verification page use the
same three. Read from the module docstrings of `plot/circular.py`, `plot/linear.py` and
`plot/sequence_view.py`, and `plot/layers.py`.

| View | Shows | Shown when |
| --- | --- | --- |
| Circle | The record's name and length, a bp scale, each feature as an arrow along its strand (a box where it has none), each primer at each binding site, each cut site labelled with its enzymes and position, each insertion as a mark labelled with how many bases it adds, and the hidden-label notice | A circular record drawn whole |
| Line | The same items along a stretch, numbered as the record numbers it, with `• • •` at an end the molecule carries on past. A feature the stretch cuts is drawn cut, with no arrowhead at the cut | Always on a page; the only shape for a region or a linear record |
| Sequence view | A position ruler, the top strand, a tick rail, the bottom strand, each row's last position. Each primer along its binding site, its 5' tail off the template and each base it does not pair with marked red. Each cut site through each strand where its enzyme cuts that strand, so the overhang shows. Each feature as a bar, with a thin line across the gap a joined feature leaves; a CDS's amino acids in three letters over their codons, a stop in red. Each insertion through both strands. Each track's curves under the rows | Up to 100,000 bases; a page leaves it and its switch out past that |

## 3. The annotation kinds

What each item carries, and which drawings need it. "Faked" marks something mbio passes as
another kind to get it drawn. Read from `plot/layers.py`, `verification/page.py` and
`protocol/render.py`.

| Kind | What it carries | Needed by |
| --- | --- | --- |
| Feature | Name, type, strand or none, one or more segments (a join), the file's colour and a segment's own where it differs; a default colour by type group (Paul Tol's light scheme), `primer_bind` purple, `regulatory` by its `/regulatory_class`; hover name, type, span and length, 1-based | Every map view |
| Joined span | Segments in top-strand order; the hover lists each run; a translation reads across the joins from `/codon_start` | Every map view |
| Span across the origin | A segment, a region or a figure span ending past the record's length, drawn as one piece and numbered as the record numbers it | Map page, figures (ligation figures can run across it), verification close-ups |
| CDS translation | One codon per amino acid, standard code, stops included, `X` where a base is not ACGT | Sequence view only |
| Primer | One item per binding site; its 5' tail length; the positions where it does not pair; label with its span, as `M13 fwd (378 .. 394)` | Every map view |
| Cut site | Every enzyme cutting at one position, merged and named alphabetically, as `BanII - SacI (406)`; each enzyme's cut count over the whole record (a unique cutter is bold) and its stagger, so the overhang shows. Positioned by the base after which the top strand is cut, as SnapGene numbers it. The shipped unique cutters by default, or every site of the enzymes named | Every map view |
| Insertion | Name, the point between two bases, how many bases it adds, colour; hover gives the point as `before^after` | Every map view; verification (a read's insertion) |
| Highlight | Names lit (a feature, a primer, or an enzyme, which lights every site it cuts); every other item dims to one grey in every view and format, and nothing moves; a lit label is the last to hide. Across a stacked figure a name lights the row that has it and dims the other rows whole | Map page, figures, verification |
| Track | Named curves sampled alike, and which sample lies over which base | Verification close-ups: a Sanger trace's four channels |
| Sequencing result | Faked as a primer over its trusted span, on its strand, its bases those it read | Verification page |
| Disagreement | Faked as a `misc_difference` feature coloured by its region's verdict (pass, warn, fail, none), or as an insertion | Verification page |
| Hidden-label notice | A count by kind, as "3 enzyme sites and 1 primer are hidden"; each hidden item's kind, name and spans | Circle and line, on every page and file |
| Source feature | Off unless asked for | Every map view |
| Well | Its name, what it holds; a fill per kind up to 10 kinds | Plate drawings |
| Band | A size in bp, per lane | Gels |

## 4. What a reader can do

**On a map page** (read from `plot/plot.js`, 668 lines, and `plot/page.py`):

1. Hover any item to see its name, type, span and length, 1-based. Hover the notice to list the
   hidden labels that show.
2. Switch each kind on or off in place (Features, Primers, Cut sites, Insertions), and each
   feature type. The notice recounts what hid among what shows.
3. Flip between circle and line, where both exist, and show or hide the sequence view. Each
   shape keeps its own zoom and scroll position.
4. Zoom the map from 1× to 8× by a slider and a reset button, or by Ctrl or ⌘ with the wheel or
   a trackpad pinch, about the pointer. The circle zooms freely and the line in four steps; a
   label hidden at 1× can show on a longer line (`tests/plot/test_linear.py:342`).
5. Drag to pan the map. A touch scrolls.
6. Click an item to select it in both views; the other view scrolls to it. Click again to clear.
7. See the sequence view beside the map when the window has room, and under it otherwise, in
   rows half as long where full ones would not fit.
8. Toggle the bottom strand.
9. Hover a base to see its position. Drag across bases to select them, the view scrolling while
   the pointer is past its edge; the caption says `first .. last (N bp)`. Copy puts the top
   strand's selected bases on the clipboard, by the button or the copy keys.
10. Print: the switches and zoom caption drop, the full rows and the line at 1× print.

There is no find box, no reverse-complement copy, no keyboard control and no export from the
page (read: `plot.js` has no `keydown` handler). The browser's own find reaches labels and the
bases within one row, since both are text (inferred from the page's markup).

**On a protocol page** (read from `protocol/render.py:2105-2186`, `protocol/protocol.js:189-198`
and `protocol/protocol.css:247-288`): a figure is a still drawing. "Explore the map" loads the
record's map page into the protocol page the first time it opens, and "Close the map" puts the
figure back; "Open it in a new tab" opens the same file. A figure prints small whether its map
was open or not. Plates and gels have no control but a gel band's tooltip and a transfer's
closed list of moves.

**On the verification page**: nothing beyond reading and the browser's own find. It has no
script (read: `verification/page.py:136-182`).

## 5. The promises a reader relies on

The map's decisions (issue #572) keep the first four as needs. Each is pinned by a test where
one is named.

1. **No label overlaps another label or the drawing** — circle, line and sequence view
   (`tests/plot/test_circular.py:316,321`, `test_linear.py:222`, `test_sequence_view.py:363`).
2. **A hidden label is counted, by kind**, on the drawing, recounted as switches change, listed
   on hover, and printed by the command after the files it wrote (`test_cli.py:69`). The
   sequence view hides nothing. A name on a feature never hides; labels hide cut sites first,
   the most frequent cutter first, then primers, features and insertions, and unlit before lit
   (`test_drawing.py:1300,1455`).
3. **Bases and labels are real text** on a page, so they can be found, selected and copied. A
   PNG and a PDF draw every letter as an outline instead, so they look the same anywhere and
   need no font (`test_drawing.py:274`).
4. **One file opens offline**: a page loads nothing over the network and embeds the fonts its
   text was measured in (`test_drawing.py:387,397`). A protocol page is one file too, but its
   "Explore the map" opens a sibling `*-map.html`, so a protocol is a folder.
5. **Positions read as the record numbers them**: 1-based and inclusive in every label and
   tooltip, a region keeps the record's numbering, across the origin too
   (`test_drawing.py:884,1156`).
6. **A highlight moves nothing** and dims the same way in every view and format
   (`test_drawing.py:1393,1446`).
7. **A page carries everything; a file draws only what is on.** A page keeps every item behind
   its switches; a PNG or PDF leaves out what is switched off (`test_drawing.py:256,649`).
8. **The same input writes the same bytes.** `scripts/check_examples.py` regenerates every
   example and compares bytes, in the docs CI job. The viewer's output, or that check, has to
   allow for this.
9. **A drawing keeps a white ground whatever the colour scheme.** A map page and the
   verification page declare a light scheme; a protocol page has a dark one and paints each
   drawing's ground white (`protocol/protocol.css:247-250`).

## 6. Every call mbio makes into drawing

Read from the call sites; `mbio.plot` is imported only by `mbio/cli.py`, `protocol/render.py`,
`protocol/figures.py` and `verification/page.py`. `synbio` imports none of it: it hands over
`Figure` specs only.

| Caller | Call | Hands over | Takes back |
| --- | --- | --- | --- |
| `mbio plot map` (`plot/cli.py:97-113`) | `draw_map(...)`, `Drawing.write(path, dpi=)` per `-o`, `layers.notice(drawing.hidden)` | A record file and 13 options: `--output/-o` (repeatable, `.html`, `.png`, `.pdf`), `--region` (a feature name, or `START..END` 1-based, END before START across the origin), `--linear`, `--sequence-view`, `--no-features`, `--no-primers`, `--no-cut-sites`, `--enzyme` (repeatable), `--hide-type` (repeatable), `--source`, `--bases-per-row` (60), `--one-strand`, `--dpi` (300) | Files written; the hidden-label line printed |
| Python users, the `plot-map` skill | `draw_map(record or path, region=, highlight=, linear=, sequence_view=, features=, primers=, cut_sites=, enzymes=, hide_types=, source=, bases_per_row=, both_strands=, insertions=, tracks=)` | As above, plus highlight, insertions and tracks | A `Drawing`: `.write`, `.element(outlines=)`, `.hidden`, `.answering` |
| Protocol render, one map page per record (`protocol/render.py:479-515`) | `draw_map(path, enzymes=).write(dir / "NAME-map.html")` | Each record any figure draws; the union of the enzymes its figures name, the unique cutters where none names any | A page; its name goes into each opener |
| Protocol render, a figure (`render.py:2105-2186`) | `draw_map(path, region=span, linear=, sequence_view=, enzymes=)`, then the highlight set by `dataclasses.replace`, `.answering`, `.element()` | A `Figure` spec, one row per record | An inline SVG per row; names each row answers to |
| Protocol render, plates (`render.py:2075-2103`, `2245-2251`) | `draw_plate(name, rows, columns, row_labels, seating=).element()`, `.layout.legend`, `.layout.kinds` | A plate's shape and seating (its labels and seating merged), or a stamp's touched wells | An inline SVG; the legend and the kind count, written as HTML |
| Protocol render, fonts (`render.py:1480-1488`) | `page.font_face(SANS, BOLD, MONO)` | Whenever the body holds a figure marked `drawing` | Three `@font-face` rules |
| Protocol render, gels (`render.py:2754-2798`) | None into `mbio.plot`: the SVG is written here | A `Gel` | An inline SVG and legend |
| Figure library (`protocol/figures.py`) | `sequence_view.LIMIT` | — | The 100,000-base cap a ligation figure checks |
| Verification page (`verification/page.py:245-281`) | `draw_map(copy, cut_sites=False, insertions=, highlight=)` for the map; `draw_map(copy, region=, sequence_view=True, bases_per_row=window, cut_sites=False, insertions=, tracks=, highlight=)` per close-up; `layers.Insertion`, `sequence_view.Track`, `sequence_view.Curve`; reads `drawing.sequence_view.rows[].strips`; `page.font_face` | The stand-in record, insertion marks and trace tracks | Inline SVGs, and whether any trace was drawn (for its colour legend) |
| PNG and PDF (`plot/convert.py`) | `vl_convert.svg_to_png`, `svg_to_pdf`, joined by `pypdf` | An SVG with letters as outlines | Bytes |

`PlateDrawing.write` (`.html`, `.png`, `.pdf`) has no caller in `src/` (read).

## 7. Who writes each drawing, and how large

Measured over `docs/examples/`: 30 HTML pages, 11,048,492 B.

| Example (writer) | Map pages | Protocol figures | Plates | Gels | Other |
| --- | --- | --- | --- | --- | --- |
| pUC19-GFP Golden Gate (`mbio cloning goldengate plan`) | product, GFP amplicon, backbone amplicon | 3: two amplicon ends at base level, the product with its junctions lit | — | 2 | `product-map.html` is also written by `mbio plot map ... --enzyme BbsI`, byte for byte |
| pUC19-GFP Gibson | product | 1: the product, both overlaps lit | — | 2 | — |
| pUC19-GFP restriction and ligation | product | 1: the product, both junctions lit | — | 4 | — |
| Gateway BP and LR | entry clone, product | 2: each product, its att sites lit | — | 1 | — |
| AP-1 library (`synbio igga plan`) | oligo, block vector 1, round 1, round 2, product | 9 on 5 pages, 3 of them stacked rows of two records | 14 on 5 pages, up to 1,536 wells | 2 | — |
| AP-1 read-back (`synbio dmx plan`) | — | — | 5 drawings | — | — |
| Maps (`mbio plot map`) | — | — | — | — | 7 PDFs: pUC19, the GFP region, AP-1 vector, rounds, product, barcode block |
| Verification (`scripts/build_example_verification.py`) | — | — | — | — | `verification.html`: a map and 2 close-ups |

Where the bytes go (measured):

- **`pUC19-GFP/product-map.html`, 656,401 B**, for a 3,347 bp record with 12 features, 9 primers
  and no cut site: 565,907 B of SVG in 7 drawings (the circle 17,893; the line four times,
  86,902; the sequence view at 60 and at 30 bases a row, 461,112), 53,717 B of fonts, 5,605 B
  of styles and 28,903 B of script. The record as compact JSON (name, topology, bases, features
  with segments, strand and colour, primers with their sites) is 5,507 B: the page is 119 times
  that.
- **AP-1 `round-1-map.html`, 977,314 B**, for 6,161 bp, 20 features and 4 cut sites (BbsI and BsaI): 886,735 B
  of SVG, 780,877 of it the sequence view twice. The record as compact JSON is 8,140 B. The
  largest page is AP-1 `product-map.html`, 1,053,579 B.
- **Fonts**: the same three faces, 53,717 B as `@font-face` rules (40,020 B of WOFF2), are
  embedded again by each of the 25 pages that show a drawing: 1,342,925 B across the examples.
- **The script**: 28,903 B, in each of the 12 map pages.
- **Plates outweigh maps on protocol pages**: 1,133,447 B of plate SVG against 205,850 B of
  figure SVG. AP-1 `05-cargo-validation-index-pcr.html`, 574,024 B, holds 414,012 B of plates.

## 8. The pointers re-measured

| Pointer (from the session that drafted the map) | Re-measured on 2026-10-17 |
| --- | --- |
| `plot/drawing.py` `_html`, about 562-595 | 562-595 |
| `plot/page.py` `render` | line 54 |
| `Figure` spec, about 817-825 | the class at 787-834, its fields at 817-825 |
| `render.py` `_figure` about 2105, `_opener` about 2150, `_write_maps` about 479-514 | 2105, 2151, 479-515 |
| `verification/page.py` stand-in record, about 255-266 | 245-266 (features from 245, the `draw_map` call at 266) |
| Plate wells, `drawing.py` about 686-795 | 685-795 |
| The hidden-label notice worded in 5 places: `layers.py:485`, `circular.py:430`, `linear.py:669`, `drawing.py:521`, `plot.js:66-97` | Worded in 2: `layers.py:485` (`notice`) and `plot.js:66-83` (`said`). The other three place what `notice` returns, and `plot/cli.py:113` prints it: 6 sites touch it |
| Fonts about 53 KB a page | 53,717 B, on 25 of the 30 example pages |
| `pUC19-GFP/product-map.html` 656,401 B, about 566 KB SVG, 7 layouts | 656,401 B; 565,907 B; 7 |
| AP-1 `round-1-map.html` 977 KB | 977,314 B |
| The record as compact JSON about 7 KB | 5,507 B for the fields listed in section 7; qualifiers and cut sites would add to it |

## Open gaps

1. **The page behaviours in section 4 were read from `plot.js`, not run in a browser.** The gate
   parses each page's HTML but runs no script (read: no browser or JavaScript runner in `tests/`
   or `pyproject.toml`), so nothing pins them. The viewer's roadmap needs its own tests of each.
2. **An explore page's cut sites follow its figures' enzymes.** A figure naming an empty list
   writes a map page with no cut sites and no switch for them: the Gibson and Gateway product
   pages have none (measured: 0 cut-site items). Whether the viewer keeps that rule, or lets the
   reader add enzymes, is the roadmap's call.
3. **A plate past 10 kinds tells the reader nothing about its wells.** The legend drops, every
   used well takes one colour, and a well's contents ride only as data attributes nothing shows
   (read: `plot/plate.py`, and no hover on protocol pages). AP-1's picked plate carries 288 well
   labels.
4. **Drawings the figure-sources note wants and nothing draws**: lines joining one stacked row's
   cut to the next row's end (iGGA Figure S1A); a plate coloured by a per-well verdict (DMX panel
   a, step 10); a Sankey of a counted run and a cost chart (DMX panels c and d); a workflow track
   (`docs/research/figure-sources.md`). That note says nothing lays out more than one record;
   stacked figure rows have since shipped, so only the joining lines are missing.
5. **The verification page fakes two kinds.** Results ride as primers and disagreements as
   `misc_difference` features. The viewer needs a read and a difference as kinds of their own;
   until then the page inherits primer labels and feature colours for them.
6. **A reader cannot find a sequence across rows**, copy the bottom strand, or move by keyboard
   on any page.
7. **Not measured**: PNG output (no example writes one), the time a page takes to open, and
   how a page behaves on a phone.
8. **The index.** `docs/research/index.md` does not list this note: the branch carries the note
   alone and is never merged.

---
search:
  exclude: true
---

# What mbio draws that molviz already does, and what it cannot

Research note for issue #583, under the map issue #573. mbio is to stop drawing on its own and
draw through liulab-molviz. This note takes every caller of drawing in mbio and sorts each
option, flag, drawing feature and supporting job into one of three classes:

- **molviz**: molviz already does it. mbio deletes its own version and calls molviz.
- **mbio**: it stays in mbio as data handed over, or as the input mbio reads and the page that
  holds the drawing. The row names the molviz field it fills.
- **gap**: molviz cannot do it. Section 6 drafts the molviz issue.

Everything was **read** from source on 2026-10-31: mbio at `8e0d608` (main), molviz at `a59ba40`
(main). Nothing was drawn or measured. Where a row says two things behave alike, both sides'
code was read; no drawing was compared.

## 1. The answer

| Class | Rows |
| --- | --- |
| molviz | 62 |
| mbio | 19 |
| gap | 1 |
| all | 82 |

The one gap is small: a stacked figure in a protocol puts an "Explore the map" control under
each row, and molviz draws a figure as one SVG. mbio can instead put the controls under the
whole figure, and then there is no gap at all (section 6).

molviz covers nearly everything because it was ported from `mbio.plot`: the type colours, the
trace colours, the plate palette, the gel's band rule, the label text and the hover text are the
same values and rules on both sides. What it does differently, and visibly:

- The sequence view writes each amino acid as one letter in a chip; mbio writes three letters.
- A drag on a page's map selects bases; in mbio it pans the map. molviz turns the circle with the
  wheel and zooms with Ctrl or ⌘ and the wheel.
- A plate's legend is drawn inside the SVG, so the protocol drops its own HTML legend.

Things the move costs mbio, none of them a molviz issue:

1. **mbio needs a public converter**, `SequenceRecord` to a molviz record. molviz's
   `scripts/write_test_records.py` is the model, but it imports `mbio.plot.layers.items` and
   `translation`, which go away with `plot/`, and four private names of
   `mbio.verification.page`: `_name`, `_verdict`, `_VERDICT_COLORS`, `_tracks`. A CDS's
   translation and a primer site's mismatches have to move into mbio's biology modules
   (`translate`, `sequence` or `primers`) before `plot/` is deleted.
2. **The converter misses one colour rule.** mbio colours a `regulatory` feature whose
   `/regulatory_class` is `terminator` or `polyA_signal_sequence` as a terminator, pink
   (`layers.default_color`). molviz sees no qualifiers, so it draws it orange as a promoter. The
   converter must set `features[].color` for it.
3. **Every shipped enzyme's cuts go into `cuts`**, not only the ones drawn. molviz counts unique
   cutters among the cuts it is handed, and lists every enzyme it is handed in the page's
   switches. The converter already does this.
4. **molviz refuses two reads of one name.** mbio verification lets two results share a name
   (`verification/cli.py`, `page.py`'s `draw_page`), so mbio must make read names unique before
   handing them over.
5. **mbio's map pages may stop passing `bases_per_row=60`.** Left out, molviz fits the rows to
   the page's width; given, it keeps them.
6. mbio can drop `vl-convert-python` and `pypdf` from its own dependencies: only
   `plot/convert.py` imports them. `liulab-molviz` brings both, and `mini-racer`, which runs the
   viewer's JavaScript in the process. The fonts in `src/mbio/data/fonts/` and
   `scripts/build_fonts.py` go too.

## 2. The command `mbio plot map`

molviz's command line has one verb, `version`. It reads no `.dna`, GenBank or FASTA file, so the
verb stays mbio's: read the file, convert it, call `molviz.draw`, write each file.

| # | Flag or job | Class | Where it goes |
| --- | --- | --- | --- |
| 1 | `RECORD`, a `.dna`, GenBank or FASTA file | mbio | `mbio.io.read_record`, then the converter: `name`, `circular`, `bases`, `features`, `primers`, `cuts` |
| 2 | The verb itself: print each file written, `error: ...` and exit 1 | mbio | stays in `plot/cli.py`, or wherever the verb moves |
| 3 | `--output`/`-o`, repeated, the suffix picking the format | molviz | `Drawing.write`, once per path; molviz also writes `.svg` |
| 4 | `--region NAME` | mbio | mbio finds the feature, case-blind, first-segment start to last-segment end, and passes `span` |
| 5 | `--region START..END`, 1-based, END before START across the origin | mbio | mbio parses it and passes `span`, 0-based, `end` past the length across the origin |
| 6 | `--linear` | molviz | `map="line"` |
| 7 | `--sequence-view`, at most 100,000 bases | molviz | `sequence_view=True`; molviz raises past `SEQUENCE_VIEW_LIMIT` |
| 8 | `--no-features`, `--no-primers`, `--no-cut-sites` | molviz | `kinds_off` |
| 9 | `--enzyme`, once per enzyme | molviz | `enzymes` |
| 10 | An `--enzyme` name no shipped enzyme answers to raises `KeyError`; an alias resolves | mbio | mbio checks names against `mbio.enzymes` before the call; molviz ignores a name that cuts nothing |
| 11 | `--hide-type` | molviz | `types_off` |
| 12 | `--source` | molviz | `types_off` without `"source"`; molviz's default leaves `source` out too |
| 13 | `--bases-per-row`, default 60, half on a narrow page | molviz | `bases_per_row`, passed only when given (section 1, item 5) |
| 14 | `--one-strand` | molviz | `bottom_strand=False` |
| 15 | `--dpi` | molviz | `Drawing.write(path, dpi=...)` |
| 16 | Print the notice of hidden labels | molviz | `Drawing.notice` |

## 3. What `mbio.plot` draws, and the jobs behind it

| # | Feature or job | Class | Where it goes |
| --- | --- | --- | --- |
| 17 | Circle layout (`circular.py`) | molviz | `map/circle.ts` |
| 18 | Line layout, a region numbered as the record is (`linear.py`) | molviz | `map/line.ts`, the `span` option |
| 19 | Sequence view: rows, both strands, numbering (`sequence_view.py`) | molviz | `sequence/view.ts` |
| 20 | Translation drawn over its codons | molviz | molviz draws one-letter chips, mbio three letters (section 1) |
| 21 | Working out a CDS's translation across joined segments, the origin and `/codon_start` (`layers.translation`) | mbio | `features[].translation`, `codon_start` |
| 22 | A primer's 5' tail bent off the template | molviz | from `primers[].bases` longer than the site |
| 23 | Finding a primer site's mismatches (`layers._mismatches`) | mbio | `primers[].sites[].mismatches` |
| 24 | Marking a mismatch on the primer's arrow | molviz | `sequence/kinds.ts` |
| 25 | Finding cut sites; dropping a cut that falls off a linear record (`find_sites`) | mbio | `cuts[]`: `enzyme`, `start`, `end`, `at`, `stagger` |
| 26 | Cuts drawn through both strands, staggered | molviz | `cuts.ts`, from `at` and `stagger` |
| 27 | Default to the unique cutters; a unique cutter's name in bold | molviz | counted among the `cuts` handed over (section 1, item 3) |
| 28 | Several enzymes at one point merged into one label, `A - B` (`merge_cuts`) | molviz | `cuts.ts` |
| 29 | Labels: `name (position)` for a cut, `name (span)` for a primer | molviz | `labels.ts` |
| 30 | Hover: name, type, span, length | molviz | `labels.tipOf`, the same four lines |
| 31 | Label placement: nothing overlaps, the map grows before it hides (`labels.py`) | molviz | `labels.ts`, `staircase.ts`, `map/columns.ts` |
| 32 | Hiding order and the hidden list (`layers.hiding`, `Drawing.hidden`) | molviz | `Drawing.hidden`, as dicts, not `layers.Item` |
| 33 | The notice's wording (`layers.notice`) | molviz | `Drawing.notice` |
| 34 | Feature colour by type, Paul Tol's light scheme; `primer_bind` purple | molviz | the same table, `reference/record.md` "Colours by type" |
| 35 | A file's own colour for a feature or a segment, SnapGene's `noColor` read as none | mbio | `features[].color`, `segments[].color` |
| 36 | A `regulatory` terminator coloured as a terminator, from `/regulatory_class` | mbio | `features[].color` (section 1, item 2) |
| 37 | Text colour on a fill, outline colour (`text_color`, `outline_color`) | molviz | `paint.ts` |
| 38 | Highlight: lit names keep their colour, all else one grey, nothing moves; a lit label hides last | molviz | `highlight` |
| 39 | A highlight naming nothing raises | molviz | `ValueError`, `options.highlight[0]` |
| 40 | Insertions, a mark at a point (`layers.Insertion`) | molviz | `insertions[]`, or a disagreement with `at` |
| 41 | Tracks: four curves under each row their peaks reach (`view.Track`, `view.Curve`) | molviz | `tracks[]`; the curve colours are molviz's, the same four hex values |
| 42 | Fonts: three DejaVu faces, measured widths (`fonts.py`, `data/fonts/`, `scripts/build_fonts.py`) | molviz | molviz's faces; `font_rules()` for a page |
| 43 | SVG writing, each letter an outline in a file (`svg.py`) | molviz | `Drawing.write(".svg")`, `still_svg()` |
| 44 | PNG at a dpi, refusing one too large (`convert.png`) | molviz | `Drawing.write(".png", dpi=...)`, the same limits |
| 45 | PDF, the map on page one and the rows after (`convert.pdf`, `drawing._pages`) | molviz | `Drawing.write(".pdf")` |
| 46 | `Drawing.write` makes the missing parent directory and returns the path | molviz | `write` returns `None` and makes no directory; the caller makes it |
| 47 | Offline page: switches for kinds, feature types and the view (`page.py`, `plot.js`, `plot.css`) | molviz | the record page `Drawing.write(".html")` writes; molviz adds a switch per enzyme |
| 48 | Page: circle and line flip | molviz | the `map` switch |
| 49 | Page: zoom, the line laid out again longer at each step | molviz | `page/zoom.ts`, 1× to 8×; `page/turn.ts` turns the circle |
| 50 | Page: drag to pan the map | molviz | molviz's drag selects bases instead; the box scrolls (section 1) |
| 51 | Page: hover details, hover on the notice lists hidden labels | molviz | `page/hover.ts`, the notice |
| 52 | Page: a click finds an item in both views; a drag selects bases to copy | molviz | `page/link.ts`, `page/caption.ts` |
| 53 | Page: sequence view beside the map or under it, narrower rows on a narrow page | molviz | `page/arrange.ts`, rows fitted to the width |
| 54 | Embedding a drawing in another page (`Drawing.element`, `page.font_face`) | molviz | `still_svg()`, and `font_rules()` once per page |
| 55 | `Drawing.answering`, which names a record draws, to check a highlight across records | molviz | `draw_figure` checks a highlight against every record |
| 56 | `Drawing.layout`, `Drawing.sequence_view`, `PlateDrawing.layout` | molviz | no equivalent needed: only tests and the two uses in rows 61 and 71 read them |

## 4. The protocol page: `protocol/render.py` and `protocol/figures.py`

| # | Feature or job | Class | Where it goes |
| --- | --- | --- | --- |
| 57 | `_write_maps`: one map page per record a figure draws, with the enzymes its figures name | molviz | `molviz.draw(record, enzymes=...).write("<stem>-map.html")`; the enzyme list stays the figure's |
| 58 | `_fonts`: the faces once, on a page that holds a drawing | molviz | `font_rules()` |
| 59 | `_figure` and `_row`: a figure's `span`, `linear`, `sequence_view`, `enzymes` | molviz | the `span`, `map`, `sequence_view` and `enzymes` options |
| 60 | A figure's records stacked, each row named | molviz | `draw_figure`, which writes each name over its record |
| 61 | A figure's highlight across rows: a row it names lights, the others dim whole; a name no row draws raises | molviz | `draw_figure(..., highlight=...)`; mbio wraps the `ValueError` with the step's name |
| 62 | An "Explore the map" control under each stacked row | gap | section 6; or mbio puts the controls under the figure |
| 63 | `ligation_figure`, `tail_figure`: the span around a junction or a primer's end, what to light | mbio | the figure spec's `span` and `highlight`, handed over as options |
| 64 | `figures.py` imports `sequence_view.LIMIT` | molviz | `molviz.SEQUENCE_VIEW_LIMIT` |
| 65 | `_plate`: a plate's wells, title, legend; the protocol's labels and seating | molviz | `draw_plate`; the legend is drawn in the SVG, and past ten kinds a count, as mbio does |
| 66 | What each well holds, `{**labels, **seating}` | mbio | plate `wells`, `rows`, `columns`, `title` |
| 67 | `row_labels` passed to `plot.draw_plate` | molviz | molviz letters rows A to Z, then AA, as `Plate.row_labels` does; the parameter goes |
| 68 | The kind count in the plate's caption | mbio | mbio counts its own wells' contents; molviz's legend lists them |
| 69 | `_stamped`: a stamp as two plates, every touched well filled with the transfer's title | mbio | two plate inputs, "From" and "Into"; molviz's plate how-to builds a stamp the same way |
| 70 | `_gel`: lane marks, ladder sizes beside the gel skipped when crowded, the legend | molviz | `draw_gel`; the same 10-unit spacing rule |
| 71 | `Gel.migration`: log of size, the largest band 0, the smallest 1, all one size 0.5 | molviz | molviz's band rule is the same; `migration` and its test go |
| 72 | Band sizes and the ladder | mbio | gel `ladder.name`, `ladder.bands`, `lanes[].name`, `lanes[].bands`, `title` |
| 73 | Gel styling in `protocol.css` (`.gel-image`, `.well`, `.band`) | molviz | molviz's look; the rules go, the print rule can stay on the wrapper |

## 5. The verification page: `verification/page.py`

The page's verdict, its table and its HTML stay mbio's: they are not drawings. Its drawings change
as follows. molviz's "Check a clone against its reads" how-to draws the same map from the same
data.

| # | Feature or job | Class | Where it goes |
| --- | --- | --- | --- |
| 74 | Each placed result drawn as a primer over the bases it trusts (`_primer`) | molviz | `reads[]`: `name`, `start`, `end`, `strand` |
| 75 | Each disagreement written into its result's bases, so the mark shows on that arrow only | molviz | `disagreements[].reads`; molviz marks those reads |
| 76 | Each disagreement drawn as a `misc_difference` feature, an insertion as a `layers.Insertion` | molviz | `disagreements[]` over a span, or at a point with `at` |
| 77 | A disagreement's name, `696 A>G` (`_name`) | mbio | `disagreements[].name`, `bases` |
| 78 | Its colour by its region's verdict (`_verdict`, `_VERDICT_COLORS`) | mbio | `disagreements[].color` |
| 79 | The map without cut sites, every read and disagreement lit | molviz | `kinds_off=["cuts"]`, `highlight` |
| 80 | A close-up: 30 bases either side, one row (`_window`, `FLANK`) | mbio | `span`, `bases_per_row`, `sequence_view=True` |
| 81 | A trace: peaks from the alignment, a bottom-strand read's channels complemented (`_tracks`) | mbio | `tracks[]`: `name`, `curves`, `peaks` |
| 82 | The trace colour legend, written when a close-up holds a trace | molviz | each strip draws its own key |

The counts in section 1 are of the 82 numbered rows in sections 2 to 5. A row that is part
kept and part deleted, such as 2 or 50, is counted by where the job goes.

## 6. Missing in molviz: the one draft issue

**Hand back each record of a figure on its own.** mbio's protocol draws a stacked figure as rows,
each row a record, and puts a control under each row that opens that record's own map page in
place. `draw_figure` returns one SVG for the whole stack, so a page cannot put anything between
its rows. Drawing each record with `draw` instead loses the figure's highlight rule: a row that
names nothing in `highlight` should dim whole, but `draw` refuses a highlight that names nothing
in its record. Either would do: a figure `Drawing` that also returns each record's still SVG,
laid out as in the stack; or a `draw` option under which a highlight that names nothing in the
record dims all of it. If mbio instead lists the controls under the whole figure, this issue is
not needed; that is the lighter choice and the recommended one.

## 7. mbio's public names under `mbio.plot`, and what each becomes

| Name | Used by | Becomes |
| --- | --- | --- |
| `draw_map` | `plot/cli.py`, `protocol/render.py`, `verification/page.py`, the plot-map skill, `docs/reference/python.md` | the converter, then `molviz.draw` |
| `Drawing` | `protocol/render.py` (type), `verification/page.py`, `docs/reference/python.md` | `molviz.Drawing` |
| `Drawing.write` | the verb, `_write_maps`, the skill | `molviz.Drawing.write` |
| `Drawing.element` | `render.py`, `verification/page.py` | `molviz.Drawing.still_svg` |
| `Drawing.hidden` | the skill | `molviz.Drawing.hidden` |
| `Drawing.answering` | `render._figure` | `draw_figure`'s highlight check |
| `draw_plate`, `PlateDrawing` | `render.py` | `molviz.draw_plate`, `molviz.Drawing` |
| `PlateDrawing.layout.legend`, `.kinds` | `render._plate` | molviz's legend; mbio counts the kinds |
| `plot.page.font_face` | `render.py`, `verification/page.py` | `molviz.font_rules` |
| `plot.fonts.SANS`, `BOLD`, `MONO`, `INVISIBLE` | `render.py`, `verification/page.py`, `scripts/build_fonts.py` | `molviz.font_rules`; the script goes |
| `plot.layers.notice` | `plot/cli.py` | `molviz.Drawing.notice` |
| `plot.layers.Insertion` | `verification/page.py` | a disagreement with `at` |
| `plot.layers.items` | `tests/methods/igga/test_ap1_demo.py`, molviz's converter script | `mbio.sites.find_sites` for which enzymes cut; the converter's own code for primers |
| `plot.layers.translation` | molviz's converter script | a biology module (section 1, item 1) |
| `plot.layers.DIM` | `tests/protocol/test_render.py` | gone; a test asserts the `highlight` handed over |
| `plot.sequence_view.LIMIT` | `protocol/figures.py`, its test | `molviz.SEQUENCE_VIEW_LIMIT` |
| `plot.sequence_view.Track`, `Curve`, `CELL` | `verification/page.py`, its test | molviz `tracks[]`; the test asserts the data handed over |
| `plot.cli.app` | `mbio/cli.py`, mounted as `plot` | stays, over molviz |

## 8. The skill, the docs and the scripts

None of these draws; each describes or runs the drawing, so each changes with it.

- `skills/plot-map/SKILL.md` keeps the verb and its flags. Its account of the page, the label
  hiding order and the errors becomes molviz's, and `from mbio.plot import draw_map` changes
  with section 7.
- `docs/methods/maps.md` keeps its commands. "Where the numbers come from" and "References"
  cite the DejaVu faces mbio ships; they should cite molviz instead. `docs/reference/cli.md`
  keeps `mbio plot map`. `docs/reference/python.md` drops `mbio.plot.draw_map` and `Drawing`.
- The method pages (`golden-gate.md`, `gibson.md`, `restriction-ligation.md`, `gateway.md`,
  `projects/igga.md`, `index.md`) link maps and protocols under `docs/examples/`; those are
  written again, not edited.
- `scripts/check_examples.py` runs `mbio plot map` and the pipelines and lists the files they
  write. The commands stay; the files are rebuilt through molviz.
- `docs/adr/0006-drawings-as-geometry.md` is superseded. The research notes under "Drawing a
  record" in the index record why `plot/` was built; they stay as history.
- `tests/plot/` (ten files) goes, and the tests above that read `plot` internals assert the data
  mbio hands molviz instead.

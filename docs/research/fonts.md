---
search:
  exclude: true
---

# Fonts: the DejaVu release, its licence, and what the shipped tables are

Research note for issue #90, carrying out the decision in #86. Read and measured on
**2026-09-16**. It records where `src/liulab_mbio/data/fonts/` comes from, what the licence asks
of it, which characters and kerning it holds, and how closely its numbers match what a browser
draws.

## 1. What ships

`scripts/build_fonts.py`, run as `pixi run build-fonts`, writes seven files, which
`liulab_mbio.plot.fonts` reads:

| File | What it is |
| --- | --- |
| `sans.json`, `sans-bold.json`, `sans-mono.json` | a face's table: advance widths, kerning pairs, glyph outlines |
| `sans.woff2`, `sans-bold.woff2`, `sans-mono.woff2` | the same face cut down to those characters, for a page to embed |
| `LICENSE` | the release's licence file, unchanged |

| Face | Cut from | Renamed to | Characters | Kerning pairs | Table | WOFF2 |
| --- | --- | --- | --- | --- | --- | --- |
| `sans` | `DejaVuSans.ttf` | Liulab Mbio Sans, Regular | 411 | 1,348 | 123,943 B | 17,248 B |
| `sans-bold` | `DejaVuSans-Bold.ttf` | Liulab Mbio Sans, Bold | 411 | 748 | 113,131 B | 16,848 B |
| `sans-mono` | `DejaVuSansMono.ttf` | Liulab Mbio Sans Mono, Regular | 96 | 0 | 21,296 B | 5,924 B |

A page embedding all three subsets carries 40,020 bytes of font, 53,360 as base64. Loading a
table the first time a face is used took about 3 ms.

A table is JSON:

| Field | Meaning |
| --- | --- |
| `family`, `style` | the names the subset carries |
| `source` | the release, its address and SHA-256, the file inside it, and where the licence is |
| `units_per_em` | 2048 for all three: the unit of every other number |
| `ascender`, `descender` | `hhea`'s, 1901 and -483 for all three |
| `advances` | each character shipped, and how far its glyph advances |
| `kerning` | each pair of characters kerned, keyed by the two, and how far the first one's advance moves |
| `outlines` | each character's glyph as SVG path data, y up from the baseline, empty for a space |

Rebuilding gave byte-identical files, from the archive on disk and from a fresh download: no
date is written, and each subset keeps the release's timestamp.

## 2. The release, pinned

The archive is `dejavu-fonts-ttf-2.37.zip` from the GitHub release `version_2_37`
([release](https://github.com/dejavu-fonts/dejavu-fonts/releases/tag/version_2_37)), 5,522,795
bytes, SHA-256 `7576310b219e04159d35ff61dd4a4ec4cdba4f35c00e002a136f00e96a908b0a`. 2.37 is the
newest release; matplotlib bundles 2.35.

GitHub publishes no digest for that asset (its API returns `"digest": null`). The project's
SourceForge file list publishes two for the same file name
([files](https://sourceforge.net/projects/dejavu/files/dejavu/2.37/)): MD5
`33e1e61fab06a547851ed308b4ffef42` and SHA-1 `bdc4350b771c6c23e8a40b2413595d1055cff2ab`. The
archive downloaded from each site matches both, and its SHA-256 is the one pinned.

The script refuses any other archive before reading it: one with a byte appended exits 1,
naming both digests, and writes nothing.

## 3. The licence

DejaVu's `LICENSE` ([text](https://dejavu-fonts.github.io/License.html)) has three parts: the
glyphs from Bitstream Vera are Bitstream's, those imported from Arev are Tavmjong Bah's, and
DejaVu's own changes are public domain. The Vera and Arev terms say the same things:

> Permission is hereby granted, free of charge, to any person obtaining a copy of the fonts
> accompanying this license ("Fonts") and associated documentation files (the "Font Software"),
> to reproduce and distribute the Font Software, including without limitation the rights to use,
> copy, merge, publish, distribute, and/or sell copies of the Font Software [...]
>
> The above copyright and trademark notices and this permission notice shall be included in all
> copies of one or more of the Font Software typefaces.
>
> The Font Software may be modified, altered, or added to [...] only if the fonts are renamed to
> names not containing either the words "Bitstream" or the word "Vera".
>
> The Font Software may be sold as part of a larger software package but no copy of one or more
> of the Font Software typefaces may be sold by itself.

Arev's terms name "Tavmjong Bah" and "Arev" instead. What is done about each:

- **Renamed.** Every name record a subset keeps that names the font (IDs 1, 2, 3, 4 and 6) is
  rewritten to Liulab Mbio Sans or Liulab Mbio Sans Mono. None holds a reserved word, nor
  "DejaVu", so a DejaVu installed on the reader's machine cannot answer to the name either.
- **Notices kept.** Each subset keeps name ID 0, the copyright, and IDs 13 and 14, the licence
  text and its address. `LICENSE` ships beside the tables and subsets, since the outlines are
  copies of the glyphs too.
- **Not sold alone.** The fonts ship inside an MIT package and the pages it writes.
- **Embedding.** `OS/2.fsType` is 0 in all three, which puts no restriction on embedding.

## 4. Characters

Sans and Sans Bold ship every character of these ranges the font has; Sans Mono ships printable
ASCII and U+FFFD.

| Range | Block | Characters shipped |
| --- | --- | --- |
| U+0020–U+007E | printable ASCII | 95 |
| U+00A0–U+00FF | Latin-1 Supplement | 95 |
| U+0370–U+03FF | Greek and Coptic | 135 |
| U+2000–U+206F | General Punctuation | 84 |
| U+2122 | ™ | 1 |
| U+FFFD | � | 1 |

Characters `plot.fonts` drops are left out: Unicode's control and format characters (Cc, Cf),
and the line and paragraph separators (Zl, Zp). In these blocks that is the soft hyphen U+00AD,
U+200B–U+200F, U+2028–U+202E, U+2060–U+2064 and U+206A–U+206F. The separators joined the dropped
set after Chrome drew U+2028 as wide as a space where the font's glyph advances by zero
(section 7). Anything else a face lacks is measured and drawn as U+FFFD.

**Why 411 and not the 284 #86 counted.** #86 and the spec name these blocks and count 284; the
count's script is not recorded. The only set found that gives exactly 284 is narrower than the
blocks: Latin-1 without the no-break space and soft hyphen, the modern Greek letters
U+0386–U+03CE, and U+2010–U+2027 alone. It leaves out the prime and double prime (U+2032,
U+2033, the typographic way to write 5′ and 3′), ‰, ‹ › and ⁄, which would then draw as U+FFFD.
Shipping the whole blocks costs 4,372 bytes of WOFF2 in Sans and 4,256 in Bold (the 284 give
12,876 and 12,592 bytes), and a page still carries less font than #86's estimate of 59 KB.

## 5. Kerning

The table's kerning is the GPOS `kern` feature as HarfBuzz applies it to Latin text.

- **Which lookups.** DejaVu Sans's `latn` script applies lookups 14 and 15 under `kern`; `DFLT`
  and `grek` apply lookup 15 alone (13 and 14, and 14, in Bold). Each is one
  pair adjustment subtable of format 2: glyphs grouped into classes, one advance change for each
  pair of classes, on the first glyph only (`ValueFormat1` 4, `ValueFormat2` 0). The script
  refuses a kerning lookup of any other shape.
- **What they hold among the characters shipped.** Lookup 15 kerns none of their pairs. Every
  pair the `latn` lookups kern holds a Latin letter: in Sans 856 are Latin–Latin, 287 Latin then
  another character, 205 another character then Latin; in Bold 544, 161 and 43. So one table of
  pairs serves every label, whatever else it spells.
- **Class 0.** Within one lookup, the first subtable whose coverage holds the first glyph
  decides the pair even when its value is zero: HarfBuzz's format 2 `apply` returns true once
  the first glyph is covered and both classes are in range
  ([source](https://github.com/harfbuzz/harfbuzz/blob/main/src/OT/Layout/GPOS/PairPosFormat2.hh)).
  Each lookup here has one subtable, so the rule changes no number today.
- **The legacy `kern` table** (2,727 pairs in Sans, 1,538 in Bold) is not read. HarfBuzz falls
  back to it only for a font without GPOS kerning, and the subsetter drops it.
- **Sans Mono** has no `kern` feature.

## 6. The subsets

Each subset is cut with fontTools' subsetter to the characters its table ships:

- **Only `kern`.** `layout_features` is `["kern"]` and GSUB is dropped, so no ligature, case form
  or contextual form can change what the table measured. Read back, the GPOS of Sans and Bold
  holds the `kern` feature alone, and Mono's none.
- **No hinting.** TrueType instructions and the tables holding them are dropped. They fit
  outlines to pixels and change no advance width: the subsets shape to the table's widths
  (section 7). Kept, they made the Sans subset 30,752 bytes rather than 17,248.
- **Also dropped:** `FFTM`, a FontForge timestamp, and `MATH`, which only Sans carries.
- **Written as WOFF2**, compressed by brotli.

## 7. Agreement with HarfBuzz and with Chrome

Chrome and Firefox shape text with HarfBuzz ([HarfBuzz](https://github.com/harfbuzz/harfbuzz)).
Each string below was measured with `Font.width` and against two references:

- **HarfBuzz 14.4.0** (uharfbuzz 0.56.1), shaping the string in one buffer with its script
  guessed from the text and ligatures off, at 2048 units to the em, where both give whole font
  units. Once with the release's TTF, once with the shipped WOFF2.
- **Chrome 153** on macOS, with the three WOFF2 subsets embedded in a page as `@font-face`, an SVG
  `<text>` at 100 px, read with `getComputedTextLength()`.

| Strings | HarfBuzz, release and subset | Chrome |
| --- | --- | --- |
| every ordered pair of characters shipped: 168,921 in Sans and in Bold, 9,216 in Mono | 0 differ | not run |
| 86 labels, as `Font.drawn` gives them | 0 differ | all within 0.03% |
| every kerned pair (1,348 Sans, 748 Bold) and every single character | 0 differ | all within 0.03% |
| `α-T` | 0 differ | 6% wider in Sans, 9% in Bold |

The labels are every record, feature and primer name in `tests/data/GFP.dna`,
`tests/data/pUC19.dna` and `docs/examples/pUC19-GFP/product.dna`, every shipped enzyme name
alone and as `name (1234)`, and a few written to break it: `lacZα`, `β-lactamase`, `5′ LTR`,
`Δ-crystallin ™`, `‰ ‘q’`, `AVATAR`, `Tet™ ×2`, a leading byte-order mark, a decomposed `é`
followed by a CJK character, a tab, U+0374, and U+2028 inside a word. Chrome's largest gap in
all of them was a thin space, 19.9707 px measured and 19.9766 px drawn, a 1/64-pixel rounding.

**Where a browser still differs.** A browser shapes a label in runs of one script and does not
kern across two. In `α-T` the hyphen joins the Greek run, so Chrome leaves `-T` unkerned while
the table kerns it, as HarfBuzz does with the whole string in one buffer. Only punctuation
between a letter of another script and a Latin letter reaches this, and no label above does
except the one written for it. A page pinning each label's width with `textLength` absorbs the
difference.

The outlines were also drawn letter by letter at `Font.letters`' positions and converted with
vl-convert: each line filled the box its width and `ascender` and `descender` gave, with the
kerning visible in `AVATAR`.

## 8. Regeneration

```sh
pixi run build-fonts                         # downloads the release
pixi run build-fonts --release PATH.zip      # reads one already downloaded
```

`build-fonts` is a task of the `fonts` feature, so pixi runs it in the `fonts` environment, the
only one with fontTools and brotli-python. That environment is outside the default solve group,
so installing it moves no version the gate runs, and no CI job installs it. The builder's tests,
`tests/scripts/test_build_fonts.py`, read no font: they check the checksum refusal, the
characters each face ships, the renaming and the kerning rule, on plain values.

## Sources

- DejaVu 2.37 release: <https://github.com/dejavu-fonts/dejavu-fonts/releases/tag/version_2_37>
- DejaVu on SourceForge, with MD5 and SHA-1: <https://sourceforge.net/projects/dejavu/files/dejavu/2.37/>
- DejaVu licence: <https://dejavu-fonts.github.io/License.html>
- OpenType GPOS, pair adjustment: <https://learn.microsoft.com/en-us/typography/opentype/spec/gpos#lookup-type-2-pair-adjustment-positioning-subtable>
- OpenType `name` IDs: <https://learn.microsoft.com/en-us/typography/opentype/spec/name#name-ids>
- HarfBuzz, `PairPosFormat2`: <https://github.com/harfbuzz/harfbuzz/blob/main/src/OT/Layout/GPOS/PairPosFormat2.hh>
- fontTools subsetter: <https://fonttools.readthedocs.io/en/latest/subset/>

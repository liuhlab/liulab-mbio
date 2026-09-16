---
search:
  exclude: true
---

# Feature colours: where a default colour for each feature type can come from

Research note for issue #79. Everything below was read or measured on **2026-09-16**. It settles
where `liulab_mbio.plot` takes the colour of a feature whose file gives it none, under a licence
this MIT package can ship, and which feature types that default has to cover. The rule it serves
is already settled on the map (#76): a feature keeps the colour its file gives it, and only a
feature without one takes a colour chosen by its type.

Each claim says whether it was **measured** (a script run over files) or **read** (a document or
source file, linked).

## 1. The verdict, first

| Question | Answer | Evidence |
| --- | --- | --- |
| **Source of the automatic colours** | **Paul Tol's light scheme**: 8 of its 9 colours, one per group of feature types | Designed for black labels on coloured fills; measured distinct in simulated colour-blind vision (section 6) |
| **Licence** | **BSD-3-Clause**, in Tol's `tol_colors.py`. Ship the copyright notice beside the values | Read, section 3 |
| **SnapGene's colours** | **Do not copy them.** There is also little to copy: its per-type defaults are nearly all one grey, and its colourful maps come from per-feature colours the `.dna` file already stores | Read and measured, section 4 |
| **Feature types to cover** | **Seven groups plus "other"**, which take every current INSDC key, the 12 keys INSDC retired in 2014, and SnapGene's `LTR` | Measured, section 5 |
| **Label text on a fill** | **Black or white, whichever contrasts more.** That never falls below 4.58:1, so it passes WCAG's 4.5:1 on any fill, a file's own included | Computed, section 7 |
| **Arrow outline** | **Always draw one.** No light-scheme colour, and 77% of the colours SnapGene files give features, reach 3:1 against a white page | Measured, section 7 |
| **Primers** | **One colour, Tol bright purple `#AA3377`**: 6.1:1 on white | Measured, section 8 |
| **Enzyme sites** | **Black text**, as SnapGene draws them | Read, section 8 |

**The reader, not the palette, is the larger finding.** The GenBank reader in `io.py` never sets
`Feature.color`. A GenBank file that carries a colour in any of the three conventions found —
ApE's `ApEinfo_fwdcolor`, DNA Features Viewer's `color`, SnapGene's `/note="color: #…"` — keeps it
as a plain qualifier, so under the settled rule every GenBank feature would take the automatic
colour even when its file names one (section 2).

**Why Tol's light scheme and not the others.** Of the twelve palettes measured, it is the only one
that does both things a labelled arrow needs: black text passes 4.5:1 on every colour, and its
closest pair stays 7.5 CIEDE2000 units apart or more in normal vision and in each simulated
colour-blind vision. ColorBrewer's Set2 and Paired, SeqViz's palette and the colours SnapGene
files carry each fall to between 0.1 and 1.8 in at least one simulation. Okabe-Ito separates further (10.9) but fails black text on its blue,
and would serve equally well with the black-or-white rule; that choice is a judgement (section 9).

## 2. What the files already hold

Read from `src/liulab_mbio/sequence.py`, `snapgene.py` and `io.py`; measured with
`read_record` on `docs/examples/pUC19-GFP/product.dna` and on a three-feature GenBank test file.

| Where | What the model keeps |
| --- | --- |
| `.dna` feature | `Feature.color` from the first segment's `color` attribute; `Segment.color` only where a segment differs. All 12 features of `pUC19-GFP` carry one |
| `.dna` write | a feature with no colour is written with SnapGene's grey `#a6acb3` (`snapgene._DEFAULT_COLOR`) |
| GenBank feature | `Feature.color` is always `None`. Colour qualifiers survive only as qualifiers |
| Primer | no colour field at all |

The GenBank test file measured the three conventions this note found in use:

| Convention | Qualifier | `Feature.color` after `read_record` |
| --- | --- | --- |
| ApE, and tools writing its format ([ApE manual](https://jorgensen.biology.utah.edu/wayned/ape/Documentation/ApE_Manual_current.pdf)) | `/ApEinfo_fwdcolor="#346ee0"`, `/ApEinfo_revcolor=…` | `None` |
| DNA Features Viewer ([source](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/blob/master/dna_features_viewer/BiopythonTranslator/BiopythonTranslator.py)) | `/color="#ffd700"` | `None` |
| SnapGene's GenBank export ([format page](https://support.snapgene.com/hc/en-us/articles/10242682237588-What-is-Genbank-SnapGene-Format)) | last `/note="color: #9d1b1c"`; a multi-segment feature lists `/ #rrggbb /` per segment | `None` |

SnapGene's format page also defines a primer colour, as one of seven names (`black red orange
green blue purple gray`) in the last `/note` of a `primer_bind` feature. None of the 1,024 `.dna`
files of section 5 stores a colour in its primer packet (**measured**), so SnapGene's primer
colour is a display default in practice.

## 3. Licences, source by source

A list of hex values may not be protectable at all. This note does not rely on that: it picks
a source whose licence says yes, and reports the rest as facts about the software.

| Source | What it offers | Licence or terms | Can we ship it? |
| --- | --- | --- | --- |
| SnapGene | per-type defaults; per-feature colours in its feature database and files | EULA forbids derivative works and reverse engineering; website forbids republishing its material | **No** |
| Addgene | maps drawn by SnapGene | noncommercial use only | **No** |
| pLannotate | 7-row colour table by type | GPL-3.0 | **No**: do not copy the table |
| DNA Features Viewer | one default, `#7245dc`; reads `/color` | MIT | Yes, but nothing per type |
| pyGenomeViz | orange for every feature; `CDS="orange"` in its GUI | MIT | Yes, but nothing per type |
| pyCirclize | reads a `/facecolor` qualifier or a caller's function | MIT | Yes, but nothing per type |
| Open Vector Editor (TeselaGen) | 90 types, 86 distinct colours | MIT | Yes; not chosen (section 6) |
| SeqViz | 12 colours, by index or at random | MIT | Yes, but nothing per type |
| Benchling | 16 colours the user picks from | proprietary | **No** |
| ApE | one default per strand; colours in its feature library | free download; no licence stated | **No** |
| SBOL Visual | glyph shapes; colour left to the author | spec CC BY 4.0; glyphs CC0 | Glyphs yes; no colours to take |
| Okabe-Ito | 7 colours and black | "feel free to use", crediting the authors | Yes, with credit |
| Paul Tol | bright 7, vibrant 7, muted 10, light 9, and more | BSD-3-Clause (`tol_colors.py`) | **Yes**, keeping the notice |
| ColorBrewer | qualitative schemes up to 12 | Apache-2.0 | Yes; not colour-blind safe past 3-4 colours |

The rows, with their sources:

- **SnapGene.** Its [terms page](https://www.snapgene.com/legal/terms-of-service) splits by
  purchase date. The [SnapGene Terms](https://www.snapgene.com/downloads/SnapGene-Terms.pdf), last
  updated 15 January 2026, cover SnapGene and the free SnapGene Viewer alike. Section 4.5: the
  user "may not alter, merge, modify, adapt, translate, decompile, reverse engineer, disassemble,
  or otherwise reduce the Software to a human-perceivable form"; 4.7: "may not create derivative
  works based upon the Software"; 5.1: "All rights not expressly granted in this Agreement are
  reserved by SnapGene." Purchases from 12 August 2026 fall under the
  [Dotmatics terms](https://www.dotmatics.com/terms-and-conditions), whose Software License
  Attachment, section 2.1, forbids the same: "modify, adapt, translate or create derivative works
  based upon the Software" and "reverse engineer, decompile, disassemble or otherwise attempt to
  derive the source code". The [legal disclaimers](https://www.snapgene.com/legal-disclaimers) add
  "you must not republish material from this website". Images are the one stated freedom: "We
  place no restriction on the publication of images generated by SnapGene or SnapGene Viewer"
  ([support article](https://support.snapgene.com/hc/en-us/articles/14097213072788-Can-I-publish-images-I-have-generated-with-SnapGene)),
  which covers a picture, not a colour table.
- **Addgene.** Its [terms of use](https://www.addgene.org/terms-of-use/) (updated 24 January
  2023): "Content may not be reproduced, duplicated, copied, sold, resold, or otherwise exploited
  for any commercial purpose without the express written consent". Its plasmid maps are SnapGene
  drawings: the pUC19 map on [Addgene #50005](https://www.addgene.org/50005/sequences/) is stamped
  "Created by SnapGene" (read from the image).
- **pLannotate.** [`colors.csv`](https://github.com/mmcguffi/pLannotate/blob/master/plannotate/data/data/colors.csv)
  and [`bokeh_plot.py`](https://github.com/mmcguffi/pLannotate/blob/master/plannotate/bokeh_plot.py)
  in a GPL-3.0 repository. Its feature database is SnapGene's: "we obtained the feature database
  used by SnapGene, which was first constructed as the GenoLIB biological part database"
  ([McGuffie and Barrick 2021](https://doi.org/10.1093/nar/gkab374), CC BY 4.0). Reading the table
  to report what it holds is not copying it; no value of it goes into `src/`.
- **DNA Features Viewer**, **pyGenomeViz**, **pyCirclize**: MIT per the GitHub licence API.
  `default_feature_color = "#7245dc"` in
  [`BiopythonTranslator.py`](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/blob/master/dna_features_viewer/BiopythonTranslator/BiopythonTranslator.py);
  `kwargs.setdefault("fc", "orange")` in pyGenomeViz's
  [`patches.py`](https://github.com/moshi4/pyGenomeViz/blob/main/src/pygenomeviz/patches.py) and
  `type2color … CDS="orange"` in its
  [`gui/config.py`](https://github.com/moshi4/pyGenomeViz/blob/main/src/pygenomeviz/gui/config.py);
  the `facecolor` qualifier in pyCirclize's
  [`track.py`](https://github.com/moshi4/pyCirclize/blob/main/src/pycirclize/track.py).
- **Open Vector Editor**:
  [`featureTypesAndColors.js`](https://github.com/TeselaGen/tg-oss/blob/master/packages/sequence-utils/src/featureTypesAndColors.js),
  MIT. **SeqViz**: [`colors.ts`](https://github.com/Lattice-Automation/seqviz/blob/develop/src/colors.ts),
  MIT.
- **Benchling.** "There are 16 colors to choose from when assigning colors to annotations. Hex
  color codes assigned to features on an import file are mapped to the closest one of these
  colors" ([help article](https://help.benchling.com/hc/en-us/articles/39768268305805-How-to-characterize-sequences-with-annotations-and-translations)).
  The page does not list them. Benchling's GenBank export has an "Export Annotation Colors"
  option ([release notes, October 2020](https://help.benchling.com/hc/en-us/articles/9684262210957-October-2020-Product-Release-Notes))
  but no public page names the qualifier it writes.
- **ApE.** The [home page](https://jorgensen.biology.utah.edu/wayned/ape/) (v3.1.10.1) offers a
  [Default Features](https://jorgensen.biology.utah.edu/wayned/ape/Download/Default_Features_11_10_23.txt)
  library and states no licence for it. The manual: "Feature default (fwd) / (rev) — default
  colors for forward and reverse strand features that have no individually assigned color."
- **SBOL Visual.** [`LICENSE.html`](https://github.com/SynBioDex/SBOL-visual/blob/master/LICENSE.html):
  "this work is licensed under a Creative Commons Attribution 4.0 International License.
  Exception: … the SBOL developers have waived all copyright and related or neighboring rights to
  SBOL Visual glyphs."
- **Okabe-Ito.** [Color Universal Design](https://jfly.uni-koeln.de/color/) (Okabe and Ito,
  2002, modified 2008): "Please feel free to use these items for your classes and seminars.
  (Please don't forget, however, to mention Masataka Okabe and Kei Ito for reference.)" The
  values are given as RGB in the page's
  [palette image](https://jfly.uni-koeln.de/color/image/pallete.jpg).
- **Paul Tol.** [`tol_colors.py`](https://sronpersonalpages.nl/~pault/data/tol_colors.py):
  "Copyright (c) 2022, Paul Tol. All rights reserved. License: Standard 3-clause BSD". The light
  set in that file is `'#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99',
  '#BBCC33', '#AAAA00', '#DDDDDD'`, matching the [scheme page](https://sronpersonalpages.nl/~pault/).
  The `tol-colors` package on PyPI is BSD-3-Clause too.
- **ColorBrewer.** Apache-2.0 per its
  [repository](https://github.com/axismaps/colorbrewer); its own colour-blind flags are read in
  section 6.

## 4. SnapGene's defaults, and why they are not a palette

**Read.** SnapGene keeps a default colour per feature type, shown in Features → Manage Feature
Types ([Standard Feature Types in SnapGene](https://support.snapgene.com/hc/en-us/articles/10383975167636-Standard-Feature-Types-in-SnapGene)).
The article's screenshots show grey for nearly every type — `-10_signal`, `3'UTR`, `enhancer`,
`gene`, `misc_feature`, `variation` among them — with CDS dark magenta, the peptide types pink,
`source` white, `exon` slate and `gap` near-white. Changing a type's default "will only apply to
new features of that type"; existing features "retain their original set color".

**Measured.** Over 1,024 files of SnapGene's public plasmid collection (section 5), the colours
files carry do not follow those type defaults:

| Type | Colours seen | Most common |
| --- | --- | --- |
| CDS | 16 | `#993366` 2,030; `#ccffcc` 1,302; `#cc99b2` 1,003 |
| promoter | 1 | `#ffffff`, all 2,708 |
| enhancer | 2 | `#ffffff` 324 of 337, though the type default is grey |
| rep_origin | 2 | `#ffff00` 1,592 of 1,594 |
| terminator | 2 | `#ffffff` 689; `#ff0000` 82 |
| misc_feature | 16 | `#ffe4c4` 595; `#99ccff` 506; `#a6acb3` 267 |
| protein_bind | 1 | `#31849b`, all 644 |
| primer_bind | 1 | `#a020f0`, all 317 |
| polyA_signal | 1 | `#a6acb3`, all 407 |

So the colour of a SnapGene map is set feature by feature, from the colour each entry in its
feature database carries ([Edit Common Features](https://support.snapgene.com/hc/en-us/articles/10383924910612-Edit-Common-Features)
lists colour among the fields of an entry), and the `.dna` reader already keeps it. What SnapGene
would give a feature **without** a colour is its type default — grey for most types. There is
no SnapGene per-type palette worth copying, and the terms forbid it anyway (section 3).

What is worth keeping is the **look**, which is not anyone's property: pale fills, a dark
outline, purple primer labels, black enzyme names.

## 5. Which feature types a default must cover

### 5.1 What real plasmid files carry

**Measured** with `read_record` over 1,024 `.dna` files of SnapGene's public collection
(`snapgene.com/local/fetch.php`): all 308 basic cloning vectors and 286 CRISPR plasmids, 3
coronavirus resources, 67 fluorescent protein plasmids, and an evenly spaced 60 from each of the
insect cell, mammalian expression, pET and Duet, plant, viral and yeast sets. 14,227 features;
14,116 carry a colour, and the 111 without (104 introns, 7 misc_features) are stored as `noColor`.

| Type | Features | Share | Cumulative | Files carrying it |
| --- | --- | --- | --- | --- |
| CDS | 4,578 | 32.2% | 32.2% | 98.4% |
| promoter | 2,708 | 19.0% | 51.2% | 93.2% |
| misc_feature | 1,695 | 11.9% | 63.1% | 70.1% |
| rep_origin | 1,594 | 11.2% | 74.3% | 93.1% |
| terminator | 771 | 5.4% | 79.7% | 45.2% |
| protein_bind | 644 | 4.5% | 84.3% | 41.4% |
| polyA_signal | 407 | 2.9% | 87.1% | 30.3% |
| enhancer | 337 | 2.4% | 89.5% | 29.6% |
| LTR | 329 | 2.3% | 91.8% | 16.1% |
| primer_bind | 317 | 2.2% | 94.0% | 15.3% |
| misc_RNA | 140 | 1.0% | 95.0% | 12.9% |
| 23 more types | 707 | 5.0% | 100% | each under 11% |

The 23: `RBS`, `intron`, `variation`, `repeat_region`, `gene`, `misc_recomb`, `mat_peptide`,
`mRNA`, `oriT`, `mobile_element`, `5'UTR`, `misc_signal`, `source`, `3'UTR`, `misc_binding`,
`exon`, `polyA_site`, `stem_loop`, `TATA_signal`, `sig_peptide`, `ncRNA`, `unsure`, `regulatory`.
RNA is where SnapGene's collection and other annotators differ most: `ncRNA` is 2 features here,
but 18 features in 9 of the 10 GenBank files pLannotate keeps as annotation controls
([`tests/test_data/annotation_controls/default/`](https://github.com/mmcguffi/pLannotate/tree/master/tests/test_data/annotation_controls/default)),
where Rfam finds them.

**Ten types cover 94.0% of features.** The long tail is short in count but long in names.

### 5.2 Two spellings of the same thing

**Read.** The INSDC feature table, [version 11.4, April 2026](https://www.insdc.org/submitting-standards/feature-table/),
defines 52 keys. Its `regulatory` key "has replaced the following Feature Keys on 15-DEC-2014:
enhancer, promoter, CAAT_signal, TATA_signal, -35_signal, -10_signal, RBS, GC_signal,
polyA_signal, attenuator, terminator, misc_signal", with the kind carried by
`/regulatory_class` from a [controlled vocabulary](https://www.insdc.org/submitting-standards/controlled-vocabulary-regulatoryclass/)
(`promoter`, `enhancer`, `terminator`, `polyA_signal_sequence`, `ribosome_binding_site`,
`riboswitch`, `insulator`, `silencer`, … and `other`). `LTR` is not a current key either.

**Measured.** SnapGene still writes the retired keys — `promoter` in 93.2% of files,
`terminator`, `enhancer`, `polyA_signal`, `RBS`, `LTR` — and only one of 14,227 features used
`regulatory`. A default therefore has to accept both spellings, and colour a `regulatory` feature
by its `/regulatory_class`.

### 5.3 Groups

Three sources split plasmid parts along the same lines: pLannotate's seven colour
rows (origin of replication, promoter, CDS, misc feature, primer bind, terminator, ncRNA), SBOL
Visual's sequence-feature glyphs (`promoter`, `cds`, `terminator`, `origin-of-replication`,
`primer-binding-site`, `non-coding-rna`, `specific-recombination-site`, `operator`, …;
[glyph directory](https://github.com/SynBioDex/SBOL-visual/tree/master/Glyphs/SequenceFeatures)),
and INSDC's regulatory classes. The proposal follows them:

| Group | Keys | Features measured | Tol light colour |
| --- | --- | --- | --- |
| Coding | `CDS gene mat_peptide sig_peptide transit_peptide propeptide proprotein exon`, and the immunoglobulin segments and regions | 33.1% | light blue `#77AADD` |
| Promoters and enhancers | `promoter enhancer RBS -10_signal -35_signal CAAT_signal TATA_signal GC_signal attenuator misc_signal 5'UTR`, and `regulatory` by class | 22.4% | orange `#EE8866` |
| Terminators and polyA | `terminator polyA_signal polyA_site 3'UTR` | 8.4% | pink `#FFAABB` |
| Replication | `rep_origin oriT` | 11.3% | light yellow `#EEDD88` |
| Binding site | `protein_bind misc_binding primer_bind` | 6.8% | light cyan `#99DDFF` |
| RNA | `ncRNA misc_RNA mRNA tRNA rRNA tmRNA precursor_RNA prim_transcript` | 1.2% | mint `#44BB99` |
| Repeat and recombination | `LTR repeat_region mobile_element misc_recomb stem_loop D-loop` | 3.4% | pear `#BBCC33` |
| Other | everything else: `misc_feature source variation intron unsure`, custom types | 13.2% | pale grey `#DDDDDD` |

Olive `#AAAA00` is left out: with pear it is the scheme's closest pair in normal, protan and
deutan vision (**measured**). Every current
INSDC key, every retired one and every key seen in the survey lands in exactly one group
(**measured**); the ones left to "other" are `misc_feature`, `source`, `variation`,
`intron`, `unsure`, `misc_difference`, `misc_structure`, `modified_base`, `old_sequence`,
`operon`, `STS`, `assembly_gap`, `gap`, `centromere`, `telomere` and `iDNA`.

With these eight colours the closest pair, in CIEDE2000 units, is 12.4 in normal vision
(light yellow and pear), 9.2 protan and 8.4 deutan (pink and pale grey), and 7.5 tritan (light
blue and mint, which is coding against RNA, 1.2% of features) (**measured**).

## 6. Colour-blind-safe palettes, measured

**Method.** For each palette: WCAG 2 contrast of black text on each colour; and the smallest
CIEDE2000 difference between any two of its colours, in normal vision and in protanopia,
deuteranopia and tritanopia simulated by [Machado, Oliveira and Fernandes 2009](https://www.inf.ufrgs.br/~oliveira/pubs_files/CVD_Simulation/CVD_Simulation.html)
(severity 1.0, applied in linear RGB). The CIEDE2000 code was checked against the first test pair
of [Sharma, Wu and Dalal 2005](https://doi.org/10.1002/col.20070) (2.0425). No
source gives a threshold for "distinct", so the numbers compare palettes rather than pass them.

| Palette | Colours | Black text ≥ 4.5:1 | Normal | Protan | Deutan | Tritan |
| --- | --- | --- | --- | --- | --- | --- |
| **Tol light** | 9 | **9** | 9.3 | 8.8 | 7.9 | 7.5 |
| Tol light, 8 hues | 8 | 8 | 9.3 | 8.8 | 7.9 | 7.5 |
| Okabe-Ito, without black | 7 | 6 | 21.7 | 12.2 | 11.6 | 10.9 |
| Tol bright | 7 | 5 | 20.5 | 14.2 | 14.8 | 8.6 |
| Tol vibrant | 7 | 5 | 18.3 | 13.7 | 14.4 | 8.6 |
| Tol muted | 10 | 6 | 15.0 | 11.8 | 14.3 | 11.5 |
| Tol pale | 7 | 7 | 7.2 | 3.6 | 4.4 | 3.3 |
| ColorBrewer Set2 | 8 | 8 | 16.8 | 1.6 | 3.6 | 4.6 |
| ColorBrewer Paired | 12 | 8 | 13.8 | 1.2 | 2.9 | 8.2 |
| pLannotate | 7 | 7 | 16.4 | 8.6 | 6.8 | 12.5 |
| SeqViz | 12 | 12 | 9.3 | 0.1 | 1.8 | 3.4 |
| SnapGene files, 10 most common colours | 10 | 8 | 12.5 | 1.9 | 1.8 | 5.8 |
| Open Vector Editor | 90 | 64 | 0.0 | — | 0.0 | — |

**Read.** Tol designed the light scheme "to fill labelled cells with more and lighter colours
than contained in the bright scheme …, but keeping black labels clearly readable. However, it can
also be used for general qualitative maps" ([scheme page](https://sronpersonalpages.nl/~pault/),
Figure 8). An arrow with its name printed on it is a labelled cell. ColorBrewer flags its own
schemes: in `colorbrewer_schemes.js` a `1` is "ok", `2` "maybe" and `0` not colour-blind friendly
([`colorbrewer.js`](https://github.com/axismaps/colorbrewer/blob/master/colorbrewer.js)); Set2 and
Dark2 are "ok" at 3 classes only, Paired at 3-4, and Set1 is "maybe" at every size.

**What the table says.**

- Only the palettes designed for colour-blind vision — Tol's and Okabe-Ito — keep every pair
  apart under every simulation. ColorBrewer's, SeqViz's and SnapGene's each fall to between 0.1
  and 1.8 in at least one. Tol's pale scheme, which Tol says is "not meant for lines or maps", is
  the exception among his.
- **One colour per type does not work.** Open Vector Editor's table gives 90 types 86 colours;
  13 pairs are under 5 units in normal vision and 73 in simulated deuteranopia, and 26 of its
  colours fail black text. Grouping is what makes a palette of 7-9 colours enough.
- **Tol light trades separation for readable labels**: its closest pairs sit nearer than
  Okabe-Ito's, but black text reads at 8.3:1 or better on every colour, and none is dark enough to
  need white text.

## 7. Label contrast and outlines

**Read.** WCAG 2.2 asks text for "a contrast ratio of at least 4.5:1" (1.4.3) and graphical
objects for "at least 3:1 against adjacent color(s)" (1.4.11), with contrast
`(L1 + 0.05) / (L2 + 0.05)` over relative luminance ([WCAG 2.2](https://www.w3.org/TR/WCAG22/)).
ApE meets the same problem with a command that "adjusts the shade of each feature color using a
luminance calculation to ensure sufficient contrast between the feature background and overlaid
text" (manual, 4.12).

**Computed.** Over every luminance an sRGB colour can have, the better of black and white text
never falls below **4.583:1**; the worst case is a mid-grey where the two are equal. The rule
"black or white, whichever contrasts more" therefore passes 4.5:1 on any fill, whatever colour a
file gives.

**Measured** over the 14,116 coloured features of section 5:

- Black text falls under 4.5:1 on **17.0%** of them — nearly all SnapGene's CDS magenta
  `#993366` (3.0:1, 2,030 features) and primer-binding purple `#a020f0` (4.0:1, 317) — and white
  contrasts more on exactly those. A black-only label rule would fail one feature in six of real
  files.
- **76.8%** of the fills are under 3:1 against a white page, the white promoters and terminators
  first among them. A pale arrow on a white page is visible only by its outline, so the outline
  is not optional. The same holds for every Tol light colour: none reaches 3:1 against white.

## 8. Primers and enzyme sites

These are drawn as text on the page rather than as filled arrows, so what matters is contrast
against white.

**Read.** SnapGene draws primer labels in purple and enzyme names in black, the unique cutters in
bold by default ([Don't Highlight Unique Cutters in Bold](https://support.snapgene.com/hc/en-us/articles/10383729806356-Don-t-Highlight-Unique-Cutters-in-Bold)),
with blunt cutters optionally green ([Highlight Blunt Cutters in Green](https://support.snapgene.com/hc/en-us/articles/10383756333716-Highlight-Blunt-Cutters-in-Green)).
The purple is visible on the Addgene pUC19 map; SnapGene's `primer_bind` features carry `#a020f0`.

**Measured** contrast on white: SnapGene's `#a020f0` 5.3:1; **Tol bright purple `#AA3377` 6.1:1**,
in `tol_colors.py` and so under the same BSD notice; Okabe-Ito's reddish purple `#CC79A7` 3.1:1,
too faint for text; black 21:1. Tol's newer dark scheme, whose colours all reach 4.5:1 on white,
has a purple `#882288` at 8.0:1, but it is on the scheme page only and not in `tol_colors.py`.

## 9. What remains a judgement for the interface ticket

- **Which hue goes to which group.** Section 5.3 is one assignment. It keeps SnapGene's yellow
  origin and puts the most common group on the calmest colour; nothing measured prefers it over
  another.
- **Tol light or Okabe-Ito.** Tol light reads with black labels throughout; Okabe-Ito separates
  further, and its seven hues fit the seven groups with a grey borrowed for "other". With the
  black-or-white rule both pass.
- **Where `primer_bind` goes.** Here it is a binding site. It could take the primer purple instead,
  so a primer and a feature describing one look alike.
- **Whether `source` is drawn at all.** It usually spans the whole record.
- **Whether a colour qualifier in a GenBank file counts as "the file's colour".** If it does, the
  reader has to lift it into `Feature.color`, in some order among the three conventions of
  section 2; that is a change to `io.py`, not to plotting.
- **Whether a primer keeps a colour from its file.** The model has no field for one, and no
  measured file stored one.
- **Dark backgrounds.** Every measurement here assumes a white page.

## Open gaps

- **The survey is SnapGene's own collection**: vendor and repository vectors, annotated by
  SnapGene. The lab's own files, and GenBank files from other tools, were not surveyed. Nine of
  its nineteen sets were not sampled, and six were sampled at 60 files each after a faster
  download was refused (HTTP 403).
- **Addgene's GenBank and `.dna` downloads returned 404** to a plain request, so no Addgene file
  was measured, only its map image.
- **Benchling's GenBank colour qualifier is undocumented** and no Benchling export was measured.
- **CIEDE2000 has no agreed "distinct" threshold**, and the colour-blind simulation is a model of
  full dichromacy; milder forms sit between normal vision and it.
- **The licence reading is not legal advice.** It takes the conservative position: ship only what
  a stated licence permits.

## Sources

All read on 2026-09-16.

- INSDC. *The DDBJ/ENA/GenBank Feature Table Definition*, version 11.4, April 2026.
  <https://www.insdc.org/submitting-standards/feature-table/>; and the `/regulatory_class`
  vocabulary, <https://www.insdc.org/submitting-standards/controlled-vocabulary-regulatoryclass/>.
- SnapGene. Terms and Conditions page, <https://www.snapgene.com/legal/terms-of-service>, leading
  to the SnapGene Terms of 15 January 2026, <https://www.snapgene.com/downloads/SnapGene-Terms.pdf>,
  and to the Dotmatics terms, <https://www.dotmatics.com/terms-and-conditions>, with their Software
  License Attachment of 7 February 2024; Legal Disclaimers, <https://www.snapgene.com/legal-disclaimers>; support articles
  10383975167636, 10242682237588, 10383924910612, 14097213072788, 10383729806356 and
  10383756333716 at <https://support.snapgene.com/>; plasmid collection, <https://www.snapgene.com/plasmids>.
- Addgene. Terms of Use, 24 January 2023, <https://www.addgene.org/terms-of-use/>; pUC19 map,
  <https://www.addgene.org/50005/sequences/>.
- McGuffie, M.J. and Barrick, J.E. (2021) pLannotate: engineered plasmid annotation.
  *Nucleic Acids Res.* 49, W516-W522. [doi:10.1093/nar/gkab374](https://doi.org/10.1093/nar/gkab374).
  CC BY 4.0. Source: <https://github.com/mmcguffi/pLannotate>, GPL-3.0.
- DNA Features Viewer, <https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer>, MIT.
- pyGenomeViz, <https://github.com/moshi4/pyGenomeViz>, and pyCirclize,
  <https://github.com/moshi4/pyCirclize>, MIT.
- TeselaGen tg-oss (Open Vector Editor), <https://github.com/TeselaGen/tg-oss>, MIT.
- SeqViz, <https://github.com/Lattice-Automation/seqviz>, MIT.
- Benchling help centre articles 39768268305805 and 9684262210957, <https://help.benchling.com/>.
- Davis, M.W. ApE, A plasmid Editor, v3.1.10.1, <https://jorgensen.biology.utah.edu/wayned/ape/>,
  with its user manual and Default Features library.
- SBOL Visual 3.0 (2021), <https://github.com/SynBioDex/SBOL-visual>: section 4.2, "The interior of
  a glyph MAY be given any fill color, as long as the choice of fill does not interfere with
  recognizing the glyph"; spec CC BY 4.0, glyphs CC0.
- Okabe, M. and Ito, K. Color Universal Design (2002, modified 2008),
  <https://jfly.uni-koeln.de/color/>.
- Tol, P. Colour schemes, <https://sronpersonalpages.nl/~pault/>; `tol_colors.py`,
  <https://sronpersonalpages.nl/~pault/data/tol_colors.py>, BSD-3-Clause; technical note
  SRON/EPS/TN/09-002, issue 3.2, <https://sronpersonalpages.nl/~pault/data/colourschemes.pdf>.
- Brewer, C.A. ColorBrewer, <https://github.com/axismaps/colorbrewer>, Apache-2.0.
- Machado, G.M., Oliveira, M.M. and Fernandes, L.A.F. (2009) A physiologically-based model for
  simulation of color vision deficiency. *IEEE Trans. Vis. Comput. Graph.* 15, 1291-1298.
  [doi:10.1109/TVCG.2009.113](https://doi.org/10.1109/TVCG.2009.113).
- Sharma, G., Wu, W. and Dalal, E.N. (2005) The CIEDE2000 color-difference formula:
  implementation notes, supplementary test data, and mathematical observations. *Color Res. Appl.*
  30, 21-30. [doi:10.1002/col.20070](https://doi.org/10.1002/col.20070).
- W3C. Web Content Accessibility Guidelines 2.2, <https://www.w3.org/TR/WCAG22/>.

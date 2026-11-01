---
name: plot-map
description: >-
  Draw a sequence record as a map with `mbio`, as SnapGene Viewer's Map and Sequence tabs
  show it: a circle or a line, the whole record or one region, its features, primers and enzyme
  cut sites labelled so no label overlaps another, and the bases with it on request. Writes one
  offline HTML page to explore, or an SVG, PNG or PDF figure. Use whenever someone wants a record
  drawn or shown: a plasmid map, a map of a `.dna`, GenBank or FASTA file, a region or insert
  drawn as a line, the sequence view of a construct, where a feature, primer or cut site lies, a
  figure of a plasmid as PNG or PDF, or a page to look through a record without opening SnapGene.
---

# Plot a map

`mbio.plot` reads the record, works out its biology (every enzyme's cuts, each CDS's
translation, each primer's mismatches) and hands it to liulab-molviz, which lays it out and
draws it. Your job is to choose the record, the stretch, what shows and the files to write,
then hand those files over. Never draw a map yourself or state a position from your own
knowledge: the drawing places every feature, primer and cut site from the record, and nothing
checks one you placed.

## Run it

```bash
pixi run mbio plot map plasmid.dna -o map.html
pixi run mbio plot map plasmid.gb -o map.png -o map.pdf --sequence-view
pixi run mbio plot map plasmid.dna --region insert -o insert.html
pixi run mbio plot map plasmid.dna --region 2680..10 --enzyme BsaI --enzyme BsmBI -o sites.pdf
```

The suffix picks the format: `.html`, `.svg`, `.png` or `.pdf`. `-o` repeats to write several
from one drawing. It prints each file written, then a notice if the map hid any label.
`--help` lists every switch.

- **`.html`** is one page that opens offline and lays the map out in the browser: details on
  hover, switches for each kind, feature type and enzyme, the circle or the line, zoom, and the
  bases in rows that fit the page. The command's switches set only what it shows first. Hand it
  to someone who wants to explore.
- **`.svg`**, **`.png`** and **`.pdf`** look the same on any machine and leave out what is
  switched off, the bases unless `--sequence-view`. The PDF puts the bases' rows on the pages
  after the map.

To look at a map yourself, write a PNG without `--sequence-view` and read the image.

## What it draws

Only what the file holds, plus cut sites: by default the shipped enzymes that cut the whole
record once, in bold. `--enzyme` names others instead, and draws every site each one cuts. A
region, and a linear record, is drawn as a line. There is no colour, size or font option: a
feature draws in its file's colour, so recolour the record rather than the map.

## Regions and positions

`--region` takes a feature's name, whatever its case, or `START..END` counted from 1 with both
ends included, as the map prints positions; END before START runs across the origin. `draw_map`
takes the span 0-based and half-open (`docs/adr/0001-coordinates.md`). Quote positions to the
user as the drawing prints them. To find a feature's exact name:

```bash
pixi run python -c "from mbio.io import read_record; print([f.name for f in read_record('plasmid.dna').features])"
```

## From Python

```python
from mbio.plot import draw_map

drawing = draw_map("plasmid.dna", region="insert", sequence_view=True, highlight="GFP")
drawing.write("map.html")
drawing.notice  # the sentence counting hidden labels by kind, empty when none hid
drawing.hidden  # each hidden label: its kind, name and place
```

`record` is a path or a `SequenceRecord`, such as a pipeline's product. `draw_map` returns
molviz's own `Drawing`. Read the docstring rather than reconstructing a call:

```bash
pixi run python -c "from mbio.plot import draw_map; help(draw_map)"
```

## A crowded map

A label that would overlap another hides: cut sites first, the enzyme cutting most often before
the rest, then primers, then feature names. A name on its own arrow never hides, and the
sequence view hides nothing. Where the user needs a hidden label, draw a region, name only the
enzymes wanted, or switch a layer off, rather than saying where it lies yourself.

## When it refuses

The command prints `error: ...` and exits 1; `draw_map` raises. A `KeyError` names an enzyme
nothing ships. A `ValueError` comes back for another suffix, a region naming no feature or lying
off the record, a highlight naming nothing drawn, a sequence view past 100,000 bases, and a PNG
too large at its `--dpi`. For the last two, draw a shorter region, or take the highest dpi or
the PDF the error names.

## Before you hand it over

Give the user each path written, and the notice if the map hid anything.

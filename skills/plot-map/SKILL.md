---
name: plot-map
description: >-
  Draw a sequence record as a map with `liulab_mbio`, as SnapGene Viewer's Map and Sequence tabs
  show it: a circle or a line, the whole record or one region, its features, primers and enzyme
  cut sites labelled so no label overlaps another, and the bases beside it on request. Writes one
  offline HTML page to explore, or a PNG or PDF figure. Use whenever someone wants a record drawn
  or shown: a plasmid map, a map of a `.dna`, GenBank or FASTA file, a region or insert drawn as
  a line, the sequence view of a construct, where a feature, primer or cut site lies, a figure of
  a plasmid as PNG or PDF, or a page to look through a record without opening SnapGene.
---

# Plot a map

`liulab_mbio.plot` lays the record out and draws it. Your job is to choose the record, the
stretch, what shows and the files to write, then hand those files over. Never draw a map
yourself or state a position from your own knowledge: the drawing places every feature, primer
and cut site from the record, and nothing checks one you placed.

## Run it

```bash
pixi run liulab_mbio plot map plasmid.dna -o map.html
pixi run liulab_mbio plot map plasmid.gb -o map.png -o map.pdf --sequence-view
pixi run liulab_mbio plot map plasmid.dna --region insert -o insert.html
pixi run liulab_mbio plot map plasmid.dna --region 2680..10 --enzyme BsaI --enzyme BsmBI -o sites.pdf
```

The suffix picks the format, and `-o` repeats to write several from one layout. It prints each
file written, then a notice if the map hid any label. `--help` lists every switch.

- **`.html`** is one page that opens offline: details on hover, a switch for each layer and
  feature type, a circle or line flip, the sequence view beside the map, a click that finds an
  item in both views, and a drag that selects bases to copy. The command's switches set only
  what it shows first. Hand it to someone who wants to explore.
- **`.png`** and **`.pdf`** look the same on any machine and leave out what is switched off, the
  sequence view unless `--sequence-view`. The PNG is one image, the sequence view under the map,
  at `--dpi`. The PDF puts the sequence view's rows on the pages after the map.

Write a figure and a page in separate runs when they should show different things. To look at a
map yourself, write a PNG without `--sequence-view` and read the image.

## What it draws

Only what the file holds, plus cut sites: by default the shipped enzymes that cut the whole
record once, in bold, each as `name (position)`, the base the top strand is cut after.
`--enzyme` names others instead, and draws every site each one cuts. The sequence view sets the
bases out in rows, each coding sequence's translation under them. A region, and a linear
record, is drawn as a line. There is no colour, size or font option: a
feature draws in its file's colour, so recolour the record rather than the map.

## Regions and positions

`--region` takes a feature's name, whatever its case, or `START..END` counted from 1 with both
ends included, as the map prints positions; END before START runs across the origin. `draw_map`
takes the span 0-based and half-open (`docs/adr/0001-coordinates.md`). Quote positions to the
user as the drawing prints them. To find a feature's exact name:

```bash
pixi run python -c "from liulab_mbio.io import read_record; print([f.name for f in read_record('plasmid.dna').features])"
```

## From Python

```python
from liulab_mbio.plot import draw_map

drawing = draw_map("plasmid.dna", region="insert", sequence_view=True)
drawing.write("map.html")
drawing.write("map.pdf")
drawing.hidden  # each label the map left out, in the order it hid
```

`record` is a path or a `SequenceRecord`, such as a pipeline's product. Read the docstrings
rather than reconstructing a call:

```bash
pixi run python -c "from liulab_mbio.plot import draw_map; help(draw_map)"
```

## A crowded map

The map grows before it hides a label, then hides enzyme sites first, then primers, then boxed
feature names; a name on its arrow never hides. The notice says how many, and the page lists
them when you hover over it. Where the user needs a hidden label, draw a region, name only the
enzymes wanted, or switch a layer off, rather than saying where it lies yourself.

## When it refuses

The command prints `error: ...` and exits 1; `draw_map` raises. A `KeyError` names an enzyme
nothing ships. A `ValueError` comes back for another suffix, a region naming no feature or lying
off the record, a sequence view past 100,000 bases, and a PNG too large to make at its `--dpi`.
The last two say what would draw instead — a region, another dpi or a PDF: pass that on.

## Before you hand it over

Give the user each path written, and the notice if the map hid anything.

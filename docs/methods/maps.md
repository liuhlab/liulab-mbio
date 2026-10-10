# Maps and figures

A map draws a record the way you would read it on paper: the molecule as a circle or a line, its
features as coloured arrows, its primers, and the places an enzyme cuts, with every label set so
none covers another. One command draws any record you hold, as a figure for a paper or as a page
to look through, and can set the bases out under it.

## When to use it

- You want a figure of a plasmid for a paper, a slide or a lab book.
- You want to see what a plan built: which features carried over, and where the joins fell.
- You want to read a stretch of bases with the features, primers and cut sites on them.

Every plan here ends in a record, so the product of [Golden Gate](golden-gate.md),
[Gibson](gibson.md), [restriction and ligation](restriction-ligation.md) and
[Gateway](gateway.md) is drawn by the command below.

## What you need

One sequence file. SnapGene `.dna`, GenBank and FASTA all read, and a `.dna` file brings its own
feature colours with it. The runs below use pUC19, which ships with the repo, and the Golden Gate
product published on this site.

## Run it

```bash
pixi run mbio plot map tests/data/pUC19.dna -o puc19-map.pdf
```

It prints the file it wrote, and nothing else:

```text
puc19-map.pdf
```

The circle holds everything the file carries: 2,686 bp, the features in their own colours, both
M13 primers, and the cut sites. No `--enzyme` was named, so it drew the shipped enzymes that cut
this record exactly once, BsaI and EcoRI among them: the ones worth planning around. Name
`--enzyme` once per enzyme to draw every site of those you name instead.

Nothing was printed after the file name, which is the second thing to read. Where a label will
not fit, the map leaves it out and prints one more line counting what hid, by kind. A hidden
label leaves you a figure, but not the one you meant.

## What it wrote

| File | What it is |
| --- | --- |
| [puc19-map.pdf](../examples/pUC19-GFP/puc19-map.pdf) | pUC19 drawn whole, as a circle |
| [product-insert.pdf](../examples/pUC19-GFP/product-insert.pdf) | the GFP of the finished plasmid, as a line |
| [product-map.html](../examples/pUC19-GFP/product-map.html) | the finished plasmid as a page to explore. The Golden Gate plan writes this same page, and its protocol's figure opens it |

The other two came from the finished Golden Gate plasmid, 3,347 bp:

```bash
pixi run mbio plot map docs/examples/pUC19-GFP/product.dna --region GFP -o product-insert.pdf
pixi run mbio plot map docs/examples/pUC19-GFP/product.dna --enzyme BbsI -o product-map.html
```

Each name links to what those commands wrote, published here unedited, and the same file always
draws the same figure. `-o` repeats, and the suffix picks the format: `.html`, `.png` or `.pdf`.

**A page keeps every item behind its own switch, whatever you asked for.** A PNG or a PDF draws
only what is switched on.

## How to read the map

The circle is the whole molecule, its name and length in the centre and a scale in bp round the
inside. Each feature is an arrow along the strand it reads on, in the colour its file gives it;
one with no colour of its own takes one by its type, and overlapping features sit on rings
further in. Primers are thin arrows outside the backbone, and each cut site is labelled with the
enzymes cutting there and the base they cut at.

`--region` draws a stretch as a line: a feature's name, as `--region GFP` does here, or
`START..END` counted from 1 with both ends included, END before START running across the origin.
The line keeps the record's own numbering, so the GFP line is captioned 396 to 1112, 717 bp.
`• • •` at each end says the molecule carries on.

A page loads nothing over the network, so it opens anywhere. It carries
the circle, the record opened as a line, that line again at two, four and eight times its length,
and the bases. Switches above the drawing turn each kind of item and each type of feature on and
off in place, and a click lights an item in both views.
`--sequence-view` sets the bases on in a PNG or a PDF too, 60 to a row
(`--bases-per-row`) and up to 100,000 of them; `--one-strand` drops the bottom strand.

## Before you publish

- Check that every feature you meant to show is drawn.
- Read the last line printed: it counts the labels left out.
- Decide whether the enzymes drawn are the ones you want, and name `--enzyme` where they are not.
- Make sure the numbering is the record's own, and not a count from 1.
- Open the file where you are sending it, at the size it will be printed.

## Tips and troubleshooting

**The figure is crowded and a label is missing.** Switch off what the figure does not need, with
`--no-primers`, `--no-cut-sites` or `--hide-type` once per feature type. One region drawn as a
line gives a label far more room than the whole circle.

**The PDF has no bases in it.** A PNG and a PDF draw only what is switched on, so ask for the
bases with `--sequence-view`. A page needs no such flag, and carries them already.

**The enzyme you need is not on the map.** With no `--enzyme`, only the enzymes that cut the whole
record once are drawn, and that count is over the whole record even where you draw one region.
Name `--enzyme` once per enzyme to draw every site of each.

## Where the numbers come from

A map judges nothing: no threshold, no verdict, none of the badges a plan or a barcode set
carries. The one number it reports is how many labels it hid. Every position and length is read
off the record itself.

The cut sites come from the enzyme data that ships with the package, the same data the cloning
commands plan with, so a site on a map and a site in a protocol are one site. Text is measured in
font tables that ship too, and a PNG and a PDF draw each letter as an outline, so the figure
needs no font installed where it is opened.

`pixi run mbio plot map --help` lists the rest: opening a circular record as a line, drawing the
source feature, and a PNG's resolution.

## References

No literature stands behind a drawing. Three things do, and the package ships all three:

- DejaVu Fonts, release 2.37. Every label is measured and drawn in three of its faces, each
  shipped as a subset with the release's licence.
- Paul Tol (2022) `tol_colors`, the light qualitative scheme, under BSD-3-Clause. A feature with
  no colour of its own takes one from it.
- The formats. A page is HTML with the drawing written in as SVG; a PNG and a PDF are drawn from
  that SVG, by vl-convert and pypdf.

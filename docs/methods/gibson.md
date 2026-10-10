# Gibson assembly

An exonuclease chews one strand back at each end, so two fragments that share the same bases
there anneal and are sealed in one tube, with no enzyme site needed anywhere. One command picks
each overlap, designs the primers that carry it, builds the plasmid you should get, and writes
the bench protocol.

## When to use it

- Two to about five fragments, joined in a defined order, in one reaction.
- You want a join that adds nothing. The shared bases are already in one of the two parts.
- No part has to be free of an enzyme site, because nothing is cut.

Use [Golden Gate](golden-gate.md) instead when a dozen fragments go in at once and an enzyme is
free to cut them. Use [restriction and ligation](restriction-ligation.md) when the two plasmids
already share a usable pair of sites.

## What you need

One file for the vector, and one for each insert, named in the order they go round the
product. SnapGene `.dna`, GenBank and FASTA all read. The run below uses the same two files as the
Golden Gate page.

## Run it

```bash
pixi run mbio cloning gibson plan tests/data/pUC19.dna tests/data/GFP.dna --out plan/
```

It prints the design in one line, then the five files it wrote:

```text
pUC19-GFP: 3346 bp, 2 fragments, overlaps 16 bp, 16 bp, NEBuilder HiFi DNA Assembly Master Mix, checks warn
plan/product.dna
plan/primers.tsv
plan/protocol.json
plan/product-map.html
plan/protocol.html
```

Read that line before anything else. Both overlaps came out at **16 bp**, inside the 15 to
20 bp the NEBuilder manual documents, where 12 bp is the shortest that works at all. Both are
taken from the pUC19 backbone and carried by the other fragment as a primer tail, so the opened
backbone belongs to no insert in particular. The kit named in that line is the assembly
product: the mix you pipette, not the plasmid you get.

`checks warn` means nothing failed, but something wants a look. Two of the overlap checks come
back with no verdict. A check with no verdict is one nothing sourced judges: it reports a number
and stops.

## What it wrote

| File | What it is |
| --- | --- |
| [product.dna](../examples/pUC19-GFP/gibson/product.dna) | the finished plasmid, 3,346 bp, features carried over and both joins marked. Opens in SnapGene |
| [primers.tsv](../examples/pUC19-GFP/gibson/primers.tsv) | the nine oligos to order, with length and melting temperature |
| [protocol.json](../examples/pUC19-GFP/gibson/protocol.json) | the same protocol as data |
| [product-map.html](../examples/pUC19-GFP/gibson/product-map.html) | the finished plasmid as a map you can explore, which the figure on the page opens |
| [protocol.html](../examples/pUC19-GFP/gibson/protocol.html) | the protocol as one page: no network, nothing to install |

Each name links to what that run wrote, published here unedited. Nothing in it is typed by
hand, so the same two input files always give the same five. To change the design, change the
command and run it again. To change what the page says, edit `protocol.json` and make the page
from it again:

```bash
pixi run mbio protocol render plan/protocol.json
```

## How to read the protocol

Open `protocol.html`. The top states the whole design in seven cards and six sentences. Here
they say where the two shared stretches sit, at 380 and at 1113, that GFP goes in on the
opposite strand from the lac promoter, that nothing is annotated ahead of GFP to start
translation, so the clone is not expected to glow, and that the insertion breaks lacZα, so
correct clones are white. All six are read off the finished plasmid, so they cannot disagree
with the map.

Under them sit the badges. `primers` is amber here: of nine oligos, three carry a warning and
none fails, on GC in two and melting temperature in one. The two without a verdict are GC, and
how alike the two overlaps are. Nobody publishes a figure for either, so read them as
unmeasured, never as fine.

Then eleven numbered steps, from the two PCRs to sequencing. Reaction tables rescale when you
change the number of reactions, and the assembly sits at 50 °C for 15 min. Gel steps draw the
bands to expect: at the screen, a correct clone gives 160 and 837 bp and an empty vector gives
177 bp, so read the big band first.

## Before you order

- Read the protocol once, top to bottom.
- Look at every badge that is not green, and decide about it.
- Read the two badges with no verdict, and answer those yourself.
- Check the product map holds what you meant to clone.
- Check the screening lanes can still be told apart.

## Tips and troubleshooting

**Every colony reads as empty vector.** The plasmid template the backbone was amplified from
survived, and it transforms far better than an assembly does. Run the DpnI step on every
amplicon made from a plasmid, and check the PCR gel gives one band at 2629 bp and one at
749 bp.

**A warning sits on a primer with a tail.** The overlap fixes where that primer binds. The
design tries every binding site the overlap allows and keeps the one that warns least, so a
warning left there is one nothing nearby could clear.

**Few colonies, or none.** Read each overlap's melting temperature and how far it folds back
on itself: here both melt at 48 °C, against a floor of 48 °C, and fold back over 4 bp,
where 10 is the limit. A highly palindromic overlap costs NEB up to tenfold fewer colonies.

## Where the numbers come from

The overlap rules are the kit's own, and so is the DNA that goes in: 0.03 to 0.2 pmol at this
fragment count, against the 0.0927 pmol this run puts in. Name another kit with `--product`
and its own rules come with it, never borrowed from one supplier to fill a gap in another.

Volumes, units, times and temperatures are the supplier's, cited step by step inside the
protocol.

`pixi run mbio cloning gibson plan --help` lists the rest: where the inserts go, which way
round each one sits, whether an insert is amplified or built from oligos, a junction joined by
one bridging oligo, and which kit, polymerase and strain the protocol names.

## References

- NEB, NEBuilder HiFi DNA Assembly Master Mix / Cloning Kit instruction manual, NEB
  #E2621S/L/X and #E5520S, version 6.0_1/26.
- New England Biolabs (2022) NEBuilder HiFi DNA Assembly Reaction (E2621), protocols.io, under
  CC BY.
- Gibson, D.G. et al. (2009) Enzymatic assembly of DNA molecules up to several hundred
  kilobases. *Nature Methods* 6, 343–345.
- NEB, Agarose Gel Resolution.

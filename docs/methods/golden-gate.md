# Golden Gate assembly

A Type IIS enzyme cuts outside its own site, so each fragment can be cut to a four-base
overhang you choose, and every fragment joined in one tube in a fixed order. One command picks
the enzyme, places the overhangs, designs the primers that carry the tails, builds the plasmid
you should get, and writes the bench protocol.

## When to use it

- Two to about a dozen fragments, joined in a defined order, in one reaction.
- You want to choose what the junction reads, down to the base, and leave no scar you did not
  pick.
- No part holds a site for one of the Type IIS enzymes, or you will let a silent codon change
  clear the sites that are there.

Use [Gibson](gibson.md) instead when the junctions may gain a few bases and you would rather
not think about enzymes. Use [restriction and ligation](restriction-ligation.md) when the two
plasmids already share a usable pair of sites.

## What you need

One file for the vector, and one for each insert, named in the order they go round the
product. SnapGene `.dna`, GenBank and FASTA all read. The run below uses two files that ship
with the repo.

## Run it

```bash
pixi run mbio cloning goldengate plan tests/data/pUC19.dna tests/data/GFP.dna --out plan/
```

It prints the design in one line, then the four files it wrote:

```text
pUC19-GFP: 3347 bp, BbsI, 2 fragments, overhangs ATGA, TGGC, fidelity 100% (measured), checks warn
plan/product.dna
plan/primers.tsv
plan/protocol.json
plan/protocol.html
```

Read that line before anything else. It chose **BbsI** because BsaI reads a site in both pUC19
and GFP, and BsmBI reads two in pUC19; cutting with either would cut the parts as well as the
ends, and BbsI cuts neither. Where no enzyme is free, the command stops and says so. It never
changes your sequence unless you ask.

`fidelity 100% (measured)` is how well those two overhangs should pair, from ligation data
measured on BbsI itself. `checks warn` means nothing failed, but something wants a look.

## What it wrote

| File | What it is |
| --- | --- |
| [product.dna](../examples/pUC19-GFP/product.dna) | the finished plasmid, 3,347 bp, features carried over and both joins marked. Opens in SnapGene |
| [primers.tsv](../examples/pUC19-GFP/primers.tsv) | the nine oligos to order, with length and melting temperature |
| [protocol.json](../examples/pUC19-GFP/protocol.json) | the same protocol as data |
| [protocol.html](../examples/pUC19-GFP/protocol.html) | the protocol as one page: no network, nothing to install |

Each name links to what that run wrote, published here unedited. Nothing in it is typed by
hand, so the same two input files always give the same four. To change the design, change the
command and run it again. To change what the page says, edit `protocol.json` and make the page
from it again:

```bash
pixi run mbio protocol render plan/protocol.json
```

## How to read the protocol

Open `protocol.html`. The top states the whole design in eight cards and four sentences. Here
they say GFP goes in on the opposite strand from the lac promoter, that nothing is annotated
ahead of GFP to start translation, so the clone is not expected to glow, and that the insertion
breaks lacZα, so correct clones are white. All four are read off the finished plasmid, so they
cannot disagree with the map.

Under them sit the badges. `primers` is amber here: of nine oligos, two carry a warning and
none fails. The oligo sheet gives each oligo its own verdict in words, and a line under the
sheet opens to show which check fired on which oligo.

Then eleven numbered steps, from the two PCRs to sequencing. Reaction tables rescale when you
change the number of reactions. Gel steps draw the bands to expect: at the screen, a correct
clone gives 160 and 838 bp and an empty vector gives 177 bp, so read the big band first.

## Before you order

- Read the protocol once, top to bottom.
- Look at every badge that is not green, and decide about it.
- Open the line under the oligo sheet and decide about each warning.
- Check the product map holds what you meant to clone.
- Check the screening lanes can still be told apart.

## Tips and troubleshooting

**Every colony reads as empty vector.** Either the plasmid template survived the DpnI digest,
or the vector closed again. Check the 60 °C soak at the end of the program ran: it is a digest,
not heat inactivation, and it cuts vector that never opened.

**A warning sits on a primer with a tail.** The junction fixes where that primer binds. The
design tries every binding site the junction allows and keeps the one that warns least, so a
warning left there is one nothing nearby could clear.

**More fragments than about a dozen.** A larger overhang set pairs less cleanly. The fidelity
in the first line is the number to read before you order.

## Where the numbers come from

Fidelity comes from ligation data measured for five enzymes. For any other enzyme one of those
five stands in, and the protocol says which it used and that it is not specific. A ligase
matrix you hold can be used instead, with `--ligase-matrix`; see
[ligation fidelity](../reference/ligation-fidelity.md).

Volumes, units, times and temperatures are NEB's, cited step by step inside the protocol.

`pixi run mbio cloning goldengate plan --help` lists the rest: where the inserts go, which way
round each one sits, holding a junction in frame, forcing an enzyme, and which polymerase and
strain the protocol names.

## References

- Engler, C., Kandzia, R. and Marillonnet, S. (2008) A one pot, one step, precision cloning
  method with high throughput capability. *PLoS ONE* 3(11): e3647.
- Pryor, J.M., Potapov, V., Kucera, R.B., Bilotti, K., Cantor, E.J. and Lohman, G.J.S. (2020)
  Enabling one-pot Golden Gate assemblies of unprecedented complexity using data-optimized
  assembly design. *PLoS ONE* 15(9): e0238592.
- NEB, NEBridge Golden Gate Assembly Kit (BsaI-HFv2) instruction manual, NEB #E1601S/L,
  version 5.0_6/26.
- NEB, Protocol for NEBridge Ligase Master Mix (NEB #M1100).

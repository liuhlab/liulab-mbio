# Gateway cloning

A phage enzyme swaps a pair of 25-base att sites between two plasmids, so an insert moves from
one backbone into the next with no cut and no ligase. One command plans the BP reaction that
makes the entry clone and the LR reaction that moves it on, designs the primers, builds both
clones, and writes the bench protocol.

## When to use it

- One insert bound for several destination vectors: make the entry clone once, and each
  destination after it is one reaction and one plate.
- The vector you want is sold as a Gateway destination and nothing else.
- Your construct can carry a whole 25-base att site at each end of the insert.

Use [Golden Gate](golden-gate.md) or [Gibson](gibson.md) when the join has to read what you
chose, and [restriction and ligation](restriction-ligation.md) when the plasmids already share
a usable pair of sites.

## What it cannot do

It looks for att sites in your own files and refuses a record that has none. An insert without
one needs `--amplify` first; a destination vector without one is refused outright, and the run
names the count it found.

## What you need

Three files: the insert or the entry clone, a donor vector and a destination vector. No record
shipped with the repo carries an att site, so this example brings its own:
[pDONR-demo](../examples/gateway/donor.gb), 661 bp with two attP sites, and
[pDEST-demo](../examples/gateway/destination.gb), 473 bp with two attR sites. Both are built
from the published att sequences, everything outside those sites is filler, and neither is a
plasmid to order. Your own run takes a real pDONR and a real pDEST.

## Run it

```bash
pixi run mbio cloning gateway plan tests/data/GFP.dna docs/examples/gateway/destination.gb \
  --donor docs/examples/gateway/donor.gb --amplify --out plan/
```

It prints the design in one line, then the seven files it wrote:

```text
pDEST-demo-GFP: 915 bp, BP then LR, 751 bp insert, junctions attB1 at 89, attB2 at 831, checks warn
plan/entry-clone.dna
plan/product.dna
plan/primers.tsv
plan/protocol.json
plan/entry-clone-map.html
plan/product-map.html
plan/protocol.html
```

Read that line before anything else. GFP carries no att site, so `--amplify` put one on each
end by PCR: each primer carries 4 G residues, the 25 bp att site and the frame bases the
fusion needs, 29 bases either way. **BP then LR** is two reactions on two days: BP makes the
entry clone, and LR moves the same 751 bp on into the destination.

`checks warn` means nothing failed, but something wants a look: two oligos out of six, on GC,
hairpin and longest run. A junction here is the whole 25-base att site, not a few bases of
scar. attB1 spells ACAAGTTTGTACAAAAAAGCAGGCT, attB2 spells ACCACTTTGTACAAGAAAGCTGGGT, and a
fusion reads through both.

## What it wrote

| File | What it is |
| --- | --- |
| [entry-clone.dna](../examples/gateway/entry-clone.dna) | the entry clone BP made, 1,037 bp |
| [product.dna](../examples/gateway/product.dna) | the expression clone, 915 bp, features carried over and both joins marked. Opens in SnapGene |
| [primers.tsv](../examples/gateway/primers.tsv) | the six oligos to order, with length and melting temperature |
| [protocol.json](../examples/gateway/protocol.json) | the same protocol as data |
| [entry-clone-map.html](../examples/gateway/entry-clone-map.html), [product-map.html](../examples/gateway/product-map.html) | each clone as a map you can explore, which the figure of the step making it opens |
| [protocol.html](../examples/gateway/protocol.html) | the protocol as one page: no network, nothing to install |

Each name links to what that run wrote, published here unedited; seven files rather than five,
because BP ran and wrote the entry clone and its map. Nothing is typed by hand, so the same three
inputs always give the same seven.
To change the design, change the command and run it again. To change what the page says, edit
`protocol.json` and make the page from it again:

```bash
pixi run mbio protocol render plan/protocol.json
```

## How to read the protocol

Open `protocol.html`. The top states the whole design in ten cards and six sentences. The 4 G
residues leave with the BP by-product, so this entry clone is the same record as one made from
an insert that already had attB ends. LR drops 751 bp into the 164 bp of destination vector
outside its own att sites, and the ccdB cassette leaves with it. GFP sits 8 bp downstream of
the T7 promoter with a ribosome binding site ahead of it, so this clone may make its protein.

Under them sit the badges, and each reaction is judged in its own right: `BP junctions` reads
attL1 at 136 and attL2 at 878, `LR junctions` attB1 at 89 and attB2 at 831. One badge carries
no verdict, because nothing sourced judges it: grow the donor and destination in One Shot
ccdB Survival 2 T1R, the only strain that grows one.

Then thirteen numbered steps, with a miniprep between the two reactions. Each is set up, run
at 25 °C, stopped with proteinase K, transformed and plated. CcdB kills TOP10, so unreacted
vector does not grow. Name an F′ strain with `--host` and that badge turns red: F′ brings
ccdA and cancels the selection. At the screen a correct clone gives 871 bp and an empty
vector 429 bp.

## Before you order

- Read the protocol once, top to bottom.
- Look at every badge that is not green, and decide about it.
- Answer the badge with no verdict: which strain grows your vectors.
- Check the product map holds what you meant to clone.
- Check the two screening lanes can still be told apart.

## Tips and troubleshooting

**The fusion has to stay in frame.** The attB tail adds frame bases only at an end you name
with `--fusion`, and none by default. Nothing in a map says which coding sequence is a tag, so
it asks rather than guessing.

**The entry clone will not grow on your plate.** Pick it on kanamycin, which comes from the
donor's own KanR. The expression clone after LR goes on ampicillin or carbenicillin instead.

**Your insert already carries an att site.** The plan looks inside the DNA that moves and
fails that check when it finds one. Recombination would then happen where nobody meant it to.

## Where the numbers come from

The att sites are the published ones, and the plan looks for them in your own file rather than
holding a copy of one vendor's vector. The 7 bases at the centre must match exactly; the rest
of the 25 may differ by one, because the vendor says these sequences drift between products.

Volumes, units, times and temperatures are Invitrogen's for the two Clonase mixes and NEB's
for the PCR and the gels, cited step by step inside the protocol.

`pixi run mbio cloning gateway plan --help` lists the rest: naming the donor, amplifying an
insert with no att ends, which end a tag sits at, and which polymerase and strain the protocol
names.

## References

- Hartley, J.L., Temple, G.F. and Brasch, M.A. (2000) DNA cloning using in vitro site-specific
  recombination. *Genome Research* 10, 1788–1795.
- Brasch, M., Cheo, D., Hartley, J. and Temple, G., US 7,670,823 B1, FIG. 9, for the att site
  sequences.
- Invitrogen, Gateway Technology with Clonase II, MAN0000470, revised 2 April 2012.
- Invitrogen, Gateway BP Clonase II Enzyme Mix, product sheet 11789.II.pps, revision
  31 October 2010.
- Invitrogen, Gateway LR Clonase II Enzyme Mix, MAN0001032 revision A.0.

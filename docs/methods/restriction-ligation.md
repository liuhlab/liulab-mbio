# Restriction and ligation

Two enzymes cut the vector and the insert to matching sticky ends, and a ligase seals them.
The pair has to cut those ends and nothing else, which is the whole of the design problem. One
command chooses the pair, designs the primers that add the sites, builds the plasmid you
should get, and writes the bench protocol.

## When to use it

- One insert into one vector, at a pair of sites.
- The two plasmids already carry a usable pair, or the insert can be amplified onto one.
- You can live with a recognition site left at each join.
- The enzymes and T4 ligase are in the freezer already, and no kit has to be bought.

Use [Golden Gate](golden-gate.md) or [Gibson](gibson.md) instead when the join has to read
what you chose, or when more than one insert goes in at once.

## What you need

One file for the vector, and one for the insert, which may be the plasmid it is cut out of.
SnapGene `.dna`, GenBank and FASTA all read. The run below uses two files that ship with the
repo.

## Run it

```bash
pixi run mbio cloning restriction plan tests/data/pUC19.dna tests/data/GFP.dna --out plan/
```

It prints the design in one line, then the four files it wrote:

```text
pUC19-GFP: 3403 bp, HindIII, SphI chosen over 89 refused pairs, 719 bp insert into a 2684 bp backbone, junctions GCATGC at 442, AAGCTT at 1165, checks warn
plan/product.dna
plan/primers.tsv
plan/protocol.json
plan/protocol.html
```

Read that line before anything else. No enzyme was named, so the pair was chosen: **HindIII
and SphI**, over 89 pairs refused, 75 on the vector site rule, 12 on the insert site rule and
2 on the backbone rule. The pair that got furthest was KpnI and XmaI, whose two backbone ends
anneal to each other, so the backbone would close on itself and the insert could go in either
way round. Name it yourself and the plan puts a dephosphorylation in instead.

`checks warn` means nothing failed, but something wants a look: two oligos out of six, both on
GC. Then read the junctions. Each end came from a recognition site, the sequence the enzyme
binds rather than the point where it cuts, and ligating two such ends puts that site back. The
product reads GCATGC at 442 and AAGCTT at 1165, and no option takes them out.

## What it wrote

| File | What it is |
| --- | --- |
| [product.dna](../examples/pUC19-GFP/restriction/product.dna) | the finished plasmid, 3,403 bp, features carried over and both joins marked. Opens in SnapGene |
| [primers.tsv](../examples/pUC19-GFP/restriction/primers.tsv) | the six oligos to order, with length and melting temperature |
| [protocol.json](../examples/pUC19-GFP/restriction/protocol.json) | the same protocol as data |
| [protocol.html](../examples/pUC19-GFP/restriction/protocol.html) | the protocol as one page: no network, nothing to install |

Each name links to what that run wrote, published here unedited. Nothing in it is typed by
hand, so the same two input files always give the same four. To change the design, change the
command and run it again. To change what the page says, edit `protocol.json` and make the page
from it again:

```bash
pixi run mbio protocol render plan/protocol.json
```

## How to read the protocol

Open `protocol.html`. The top states the whole design in seven cards and seven sentences. The
one to read first says the insert is amplified rather than cut out, so its ends come from its
primers: 6 spacer bases then the SphI site on the forward primer, 6 then HindIII on the
reverse. The spacer is what lets an enzyme cut a site that close to a fragment's end, and the
741 bp amplicon is digested in place of a plasmid.

Under them sit the badges, and the digest ones settle how the day runs. Both enzymes are
supplied in rCutSmart Buffer and both cut at 37 °C, so one tube takes both. Heat stops both,
HindIII-HF at 80 °C and SphI-HF at 65 °C, 20 min each. The backbone is cut to a 5' AGCT and a
3' CATG, which do not anneal, so it cannot close on itself and needs no phosphatase.

Then twelve numbered steps, from the PCR to sequencing. The preparative gel in the middle
recovers 2684 bp of backbone and the 719 bp insert, and each slice gives back 70% to 90% of
the DNA in it. At the screen a correct clone gives 828 bp and an empty vector 111 bp, and a
miniprep cut with both enzymes gives 2684 bp and 719 bp.

## Before you order

- Read the protocol once, top to bottom.
- Look at every badge that is not green, and decide about it.
- Check neither restored site sits where you will want to cut later.
- Check the product map holds what you meant to clone.
- Check the screening lanes can still be told apart.

## Tips and troubleshooting

**Every colony reads as empty vector.** This backbone cannot close on itself, so suspect a
digest that did not finish: a vector cut at one site still transforms. Take both digests to
completion, and cut the backbone out of the gel before you ligate.

**pUC19 gives only one band on the preparative gel.** The second piece is 2 bp and runs off
the bottom, past the smallest ladder band. Cut the 2684 bp band out and ignore the rest of
the lane.

**The restored site is in your way.** Name a different pair with `--enzyme` and read the
junctions it leaves instead. Or use [Golden Gate](golden-gate.md) or [Gibson](gibson.md),
where the join reads what you choose.

## Where the numbers come from

The pair is chosen by rules rather than by a score. An enzyme has to cut the vector where the
insert goes and nowhere else, leave the insert whole, and leave two backbone ends that cannot
anneal to each other.

Buffer, temperature, the heat that stops a digest, what Dam and Dcm methylate, the 3:1 ratio and
what a gel slice gives back are NEB's, cited step by step inside the protocol.

`pixi run mbio cloning restriction plan --help` lists the rest: naming the enzymes yourself,
and which polymerase, strain and product name the protocol carries.

## References

- NEB, Optimizing Restriction Endonuclease Reactions, for the digest and star-activity
  conditions.
- NEB, Double Digests, for cutting with two enzymes in one buffer.
- NEB, Ligation Protocol with T4 DNA Ligase (NEB #M0202), for the ligation, its ratios and
  both incubations.
- NEB, Restriction Endonuclease Technical Guide, version 5.0 – 7/17, for what Dam and Dcm
  methylate.
- NEB, Cleavage Close to the End of DNA Fragments, for the bases a site needs 5' of it to be
  cut on a PCR product.

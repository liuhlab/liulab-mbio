# Clone GFP into pUC19 with two enzymes

This page follows one cloning job from end to end. The vector is pUC19. The insert is the GFP
coding sequence. Two enzymes cut both, and a ligase joins them.

This is the oldest way to clone, and the cheapest. There is no kit. Every plasmid map in the
lab is drawn around it.

Both sequence files ship with the repo, under `tests/data/`.

## Run it

```bash
pixi run liulab_mbio cloning restriction plan tests/data/pUC19.dna tests/data/GFP.dna --out plan/ --enzyme EcoRI --enzyme BamHI
```

It prints one summary line, then the four paths it wrote:

```text
pUC19-GFP: 3388 bp, EcoRI, BamHI, 723 bp insert into a 2665 bp backbone, junctions GAATTC at 396, GGATCC at 1119, checks warn
```

Read that line first. It is the whole design in one sentence.

| Part of the line | What it says |
| --- | --- |
| `3388 bp` | how long the finished plasmid is |
| `EcoRI, BamHI` | the two enzymes that cut both pieces |
| `723 bp insert into a 2665 bp backbone` | the two bands to cut out of the gel |
| `junctions GAATTC at 396, GGATCC at 1119` | what each join now reads, and where it sits |
| `checks warn` | something is worth a look before you order |

### Why name the enzymes

You do not have to. Leave `--enzyme` off and the pair is chosen for you.

The package ships 18 enzymes that cut inside their own site, which makes 153 pairs. On these
two files it refuses 89 of them and ranks the other 64. A pair is refused when it cuts nothing
in pUC19, when it cuts inside the insert, when the cut backbone can close on itself, or when
the digest throws away more of pUC19 than it keeps. Each refusal is one sentence, and the plan
keeps them all.

What is left is ranked. First on what the bench says: the buffer both enzymes come in, the
temperature they want, whether heat stops them, and whether the host's own methylation blocks
them. Then on what the cloning costs: how much of the vector it gives up, and how much of the
insert it carries whole.

Name the pair when you have a reason the sequences cannot show. Here the reason is the freezer.
The pUC19 multiple cloning site was built for EcoRI and BamHI, and most labs stock both.

### Why `checks warn`

Six oligos were designed and two carry a warning, both for GC content. Neither fails. The
badge names the kind, and the oligo sheet gives every oligo a verdict of its own.

Both warnings sit on the two primers that amplify GFP. Those primers carry a tail with the
recognition site in it, so the design cannot move them far. A warning left there is one that
nothing nearby could clear.

Every other check passes. Two are worth reading anyway:

- **The reading frame.** GFP lands inside lacZα, and it adds 723 bases, which is a whole number
  of codons. So the frame holds. Had it not, the check would say so here rather than after
  sequencing.
- **The clean-up.** The supplier states no heat step for BamHI-HF, so heat cannot stop it. The
  gel takes it away instead. Do not swap that step for a heat step.

## What comes out

| File | What it is |
| --- | --- |
| `product.dna` | the finished plasmid, with every feature carried over and both joins marked. Opens in SnapGene |
| `primers.tsv` | the six oligos, with length and melting temperature. This is your order |
| `protocol.json` | the bench protocol written as data. The page is made from this file |
| `protocol.html` | the bench protocol: one page, no network, nothing to install |

Nothing here is written by hand, so the same two input files always give the same four files.
To change the design, run the command again.

To change what the page says without changing the design, edit `protocol.json` and turn it back
into a page:

```bash
pixi run liulab_mbio protocol render protocol.json
```

## How to read the protocol

Open `protocol.html` in a browser. It is built to be followed at the bench, not filed away.

**The top of the page** carries seven cards: the vector, the insert, the enzymes, the backbone,
the junctions, the product and what to select on. Under them sit a few sentences about the
clone, then the badges. A badge that is not green spells out why beneath it.

**The middle** is twelve numbered steps. GFP is amplified first, because the file holds the
coding sequence on its own and nothing else carries the two sites. Then both pieces are cut,
run on a gel and cut out of it. Then the ligation, the transformation, and three ways of
checking the result.

Reaction tables rescale when you change the number of reactions at the top of the table. Gel
steps draw the bands to expect.

**The ligation table asks for picomoles, not nanograms.** It wants three insert molecules for
each backbone molecule. GFP is short, so three of it weigh less than one backbone: 26.72 ng of
insert against 32.83 ng of backbone. Weighing them out equal would starve the reaction of
insert.

**The screening step** is the one to read twice. A correct clone gives one band at 829 bp. An
empty vector gives one at 127 bp. There is no third lane, because the two ends of this insert
cannot pair with each other. GFP can only go in one way round.

**The last check before sequencing** cuts a miniprep with the same two enzymes. A correct clone
falls into 2665 bp and 723 bp. A colony that carried plain pUC19 gives 2665 bp and 21 bp, and
that 21 bp piece runs off the bottom of the gel. So read the 723 bp band: it is there or it is
not.

## What the junction costs

Every join here reads the site the two ends came from. The product gains those bases. It cannot
be otherwise: the enzyme cut inside its own site, so putting the ends back together writes the
site again.

That is the trade. The junction is not clean, and in return nothing has to be designed for it.
The plan names what each junction reads, and it checks whether those extra bases hold your
reading frame. If you need a join that adds nothing, use Golden Gate instead — see
[Put GFP into pUC19](golden-gate.md).

## When one enzyme does both ends

Name one enzyme and the plan still works, and two things change.

The cut vector can now close on itself with nothing in it, which would fill your plate with
empty vector. So the plan adds a step: a phosphatase takes the phosphates off the cut vector
before the ligation. Without a phosphate the ligase cannot seal the vector back to itself.

The insert can also go in either way round. So the screening gains a third lane, and it uses a
primer that reads out of the insert itself. Two primers sitting in the vector would give the
same band both ways round and tell you nothing.

## When the insert already sits in another plasmid

Then it is cut out rather than amplified, and you order no primers for it at all.

The command takes that plasmid as its second argument, and reads what the record carries:

| What that record holds | What the plan does |
| --- | --- |
| a site for each enzyme | cuts the insert out and takes it off a gel |
| neither site | amplifies it, with a spacer and the site on each primer tail |
| any other count | stops, and names every site it found |

This is where the method is cheapest. Two enzymes, a ligase, and nothing to order.

## A check with no verdict

One check here may come back neither green nor red.

Two enzymes can only share a tube when both work in one buffer. The package knows which buffer
each enzyme is sold in, and it says so when both are sold in the same one. When they are not,
answering the question needs figures for how much activity each keeps in the other's buffer,
and nobody publishes those under a licence this package may ship.

So that check carries no verdict. The page shows it as `not judged`, and it names the
supplier's own table to look the pair up in. It is never shown as a pass. A number nobody
measured is worse than no number at all.

For EcoRI and BamHI the question does not arise: both are sold in the same buffer, so the check
passes and the reaction table names it.

## When it refuses

It stops rather than guessing, and the message says why. Ask for BsaI and EcoRI on these files
and it answers:

```text
error: GFP reads 1 BsaI site (at 643) and no EcoRI site; this method cuts the insert out between two cuts -- one site of each enzyme, or two of one -- and a record reading none of them is amplified with the sites on its primer tails instead, so neither route reaches two cut ends; a synonymous codon change reaches every one of them, leaving the insert to be amplified instead: BsaI at 643 by GAC to GAT in GFP
```

The last clause is the useful part. That BsaI site sits inside a coding sequence, so one codon
could be swapped for another that reads the same amino acid, and the site would be gone. The
protein would not change. Whether to do it is your call, not the tool's.

## Before you order

- Read the whole protocol once, top to bottom.
- Look at every badge that is not green, and decide about it.
- Look for any badge reading `not judged`, and answer it yourself.
- Check that the product map holds what you meant to clone.
- Check that the screening bands can still be told apart.

`pixi run liulab_mbio cloning restriction plan --help` lists the rest of the options: which
polymerase to amplify with, which strain the protocol names, and what to call the product.

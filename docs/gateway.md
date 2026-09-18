# Move a gene into an expression vector with Gateway

This page follows one cloning job from end to end. Nothing is cut here, and nothing is joined by
a ligase. Two short sites swap their halves, and the insert moves from one plasmid into the
next.

Gateway is what a shelf of expression vectors is built for. You make the entry clone once. Every
destination after that is one reaction and one plate.

## What you bring

Three files, and all three are yours.

| File | What it is |
| --- | --- |
| your gene | the insert. A plain coding sequence, a fragment that already has att ends, or an entry clone you made last year |
| a donor vector | where the entry clone comes from, such as a pDONR plasmid |
| a destination vector | where the gene ends up, carrying the promoter and the tag you want |

No vector ships with this package, and no list of vectors either. Export your own maps from
SnapGene, or download them from the supplier.

That is on purpose. The supplier's own att sites drift between products. A tool holding a copy
of a site would match the old one and quietly miss yours. So the plan reads your file and looks
for the sites in it, allowing a base to differ outside the core of each one. The reasoning is in
[the att sites decision](adr/0008-att-sites-and-vectors.md).

## Run it

Three commands, and the files you hold decide which one you want.

An entry clone already in hand, so only the LR reaction is needed:

```bash
pixi run liulab_mbio cloning gateway plan entry.dna destination.dna --out plan/
```

A fragment that already carries att ends, so BP runs first:

```bash
pixi run liulab_mbio cloning gateway plan insert.dna destination.dna --donor donor.dna --out plan/
```

A plain gene, which has to be amplified onto att ends first:

```bash
pixi run liulab_mbio cloning gateway plan gene.dna destination.dna --donor donor.dna --amplify --out plan/
```

## What it prints

One summary line, then every path it wrote. Read the line first. It is the whole design in one
sentence:

```text
NAME: LENGTH bp, BP then LR, INSERT bp insert, junctions attB1 at START, attB2 at START, checks pass
```

| Part of the line | What it says |
| --- | --- |
| `LENGTH bp` | how long the finished plasmid is |
| `BP then LR` | which reactions it planned. It reads `LR` alone when you brought the entry clone |
| `INSERT bp insert` | how much DNA moves between the att sites |
| `junctions attB1 at ...` | what each join now reads, and where it sits |
| `checks pass` | the worst verdict of every check, so read it before you order |

## What comes out

| File | What it is |
| --- | --- |
| `entry-clone.dna` | the entry clone, written only when the run made one |
| `product.dna` | the expression clone, with every feature carried over and both joins marked. Opens in SnapGene |
| `primers.tsv` | every oligo it designed, with length and melting temperature. This is your order |
| `protocol.json` | the bench protocol written as data. The page is made from this file |
| `protocol.html` | the bench protocol: one page, no network, nothing to install |

So five files where the BP reaction ran, and four where you brought the entry clone. Nothing
here is written by hand, so the same files in always give the same files out. To change the
design, run the command again.

To change what the page says without changing the design, edit `protocol.json` and turn it back
into a page:

```bash
pixi run liulab_mbio protocol render protocol.json
```

## How to read the protocol

Open `protocol.html` in a browser. It is built to be followed at the bench.

**The top of the page** carries a card for each vector, each clone, the reaction, the insert,
the junctions, the product and what to select on. Under them sit a few sentences about the
clone, then a badge for each check.

One badge is worth knowing about before you start. An att site inside your own insert would
recombine where nobody meant it to, and you would learn it from a plate of wrong clones. The
plan looks for one on both strands, and fails that check when it finds one.

**The middle** is the bench work. A run from a plain gene has thirteen numbered steps, in five
groups:

- amplify the gene onto att ends, then clean the product up
- set up BP, run it for an hour at 25 °C, stop it with proteinase K, transform and plate
- pick a colony, grow it overnight and prep the entry clone
- the same four steps again for LR
- screen the colonies, then sequence the clone

**The miniprep is a step of its own, and that is the point.** LR takes clean entry clone DNA,
weighed out. A stopped BP reaction is not that. No plan can weigh a yield nobody has measured,
so the page says what the step is for and leaves the weight to your own reading.

**The tables count molecules, not just weight.** BP takes the same number of molecules of each
piece of DNA, and a short fragment weighs less than a vector at the same count. LR asks for 50
to 150 ng of the entry clone. Both units are printed so that you do not do the arithmetic at the
bench.

**The screening step** has two lanes, and the page draws both. A correct clone gives one band. A
destination vector that never reacted gives a band of another size, because the cassette it
still carries is not the insert. There is no third lane: the two att sites differ from each
other, so the insert cannot go in the wrong way round.

**The sequencing primers are designed against your own product**, reading in from outside each
join. The supplier's named primers are not used. One pair suits a single kit's vector, and the
others cross well over a hundred bases of vector before they reach the insert.

## The junction is not clean

The insert arrives with a whole att site at each end, so the clone gains 25 bases there. Where
your protein is meant to run into a tag, the ribosome reads through those bases, and the frame
has to hold.

The plan checks that frame only at an end you name:

| `--fusion` | What it changes |
| --- | --- |
| `none` | nothing is added, and neither frame is judged |
| `N-terminal` | two more bases on the forward primer, so the tag in front stays in frame |
| `C-terminal` | one more base on the reverse primer, so the tag behind stays in frame |
| `both` | both ends |

Why name it rather than let the tool work it out? Nothing in a plasmid map says which coding
sequence is a tag. The one nearest the join is as likely to be the resistance gene. A tool
guessing there would fail plenty of correct work, so it asks instead.

## Two strains, and they are not the same one

A donor and a destination vector both carry ccdB, a gene that kills ordinary *E. coli*. That is
what keeps your plate clean: after the reaction, leftover vector cannot grow.

It also means two different strains:

- Grow the vectors themselves in the resistant strain the page names. An ordinary strain will
  not take them.
- Select the clone in an ordinary strain, which ccdB kills.
- Never use a strain carrying the F′ episome. It brings *ccdA*, which cancels ccdB, and your
  plate fills with background. The plan fails that check by name and tells you why.

## A check with no verdict

One check here comes back neither green nor red.

The manuals say which strain grows a ccdB vector, and they forbid the F′ one. None of them says
what happens if you plate the reaction itself on a resistant strain. It should lose the
selection, but that is reasoning rather than a measurement.

So the page shows the check as `not judged`, with the guidance beside it. It is never shown as a
pass. A number nobody measured is worse than no number at all.

## When it refuses

It stops rather than guessing, and the message says what it looked for. Hand it a plasmid with
no att sites and it answers:

```text
error: pUC19 carries 0 attR1 sites, needing one: looked for attR1 and attR2 on either strand, matching the 7 bp overlap exactly and the rest of the 25 bp region within 1 base(s)
```

Read that as: this is not a destination vector, or it is not the file you meant. The same
message names attL sites when it was given something that should be an entry clone.

## Before you order

- Read the whole protocol once, top to bottom.
- Look at every badge that is not green, and decide about it.
- Look for the badge reading `not judged`, and answer it yourself.
- Check that the product map holds what you meant to clone.
- Check that the two screening bands can still be told apart.

`pixi run liulab_mbio cloning gateway plan --help` lists the rest of the options: which
polymerase to amplify with, which strain the protocol names, and what to call the product.

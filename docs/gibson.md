# Join GFP into pUC19 by Gibson

This page follows one cloning job from end to end. The vector is pUC19. The insert is the GFP
coding sequence. Nothing is cut, so no enzyme has to be free. The two fragments share a short
stretch of bases at each join, and one kit does the rest.

Both sequence files ship with the repo, under `tests/data/`.

## Run it

```bash
pixi run liulab_mbio cloning gibson plan tests/data/pUC19.dna tests/data/GFP.dna --out plan/
```

It prints one summary line, then the four paths it wrote:

```text
pUC19-GFP: 3346 bp, 2 fragments, overlaps 16 bp, 16 bp, NEBuilder HiFi DNA Assembly Master Mix, checks warn
```

Read that line first. It is the whole design in one sentence.

| Part of the line | What it says |
| --- | --- |
| `3346 bp` | how long the finished plasmid is |
| `2 fragments` | the opened vector, plus GFP |
| `overlaps 16 bp, 16 bp` | how many bases the parts share at each join |
| `NEBuilder HiFi DNA Assembly Master Mix` | the kit the numbers came from |
| `checks warn` | something is worth a look before you order |

### Why nothing was cut

BsaI reads a site in both pUC19 and GFP, and BsmBI reads two in pUC19. That rules both enzymes
out of a Golden Gate plan. Here it does not matter at all. The vector is opened by
PCR rather than by a digest, and the join is made by bases the two fragments share. So the parts
go in as they are, and nothing has to be rewritten to make an enzyme fit.

The shared bases are the overlap. Both joins take theirs from pUC19, and GFP's two primers carry
them as tails. That way the opened backbone belongs to no insert in particular, and you can make
it once and reuse it.

### Why `checks warn`

Nine oligos were designed and three carry a warning: GC content on two of them, melting
temperature on one. None fails. A tail fixes where such a primer binds, so a warning left on one
is a warning nothing nearby could clear.

Two checks say `not judged` rather than passing. Nobody publishes a GC band for an overlap, and
nobody says how alike two overlaps may be before they join in the wrong order. Where no
published number exists, the page says so. Read that as unmeasured, never as fine.

## What comes out

| File | What it is |
| --- | --- |
| `product.dna` | the finished plasmid, with every feature carried over and both joins marked. Opens in SnapGene |
| `primers.tsv` | the nine oligos, with length and melting temperature. This is your order |
| `protocol.json` | the bench protocol written as data. The page is made from this file |
| `protocol.html` | the bench protocol: one page, no network, nothing to install |

Nothing here is written by hand, so the same two input files always give the same four files.
To change the design, run the command again.

## How to read the protocol

Open `protocol.html` in a browser. It is built to be followed at the bench.

**The top of the page** carries seven cards: the vector, the insert, the fragment count, the
overlaps, the kit and its incubation, the product, and what to select on. Under them sit six
sentences about what the clone should do, then a badge for each check.

Two of those sentences matter before you start:

- GFP goes in on the opposite strand from the lac promoter, so that promoter does not read it.
- Nothing is annotated ahead of GFP to start translation, so the clone is not expected to glow.

**The middle** is eleven numbered steps. Two PCRs, a gel, a DpnI digest, a column cleanup, a
reading on the spectrophotometer, the assembly, the incubation, the plate, the colony PCR and
the sequencing. Each step says what to do, what a good result looks like, and what to do when it
is wrong.

**The assembly table** gives each fragment in picomoles and in nanograms. The vector goes in at
50 ng. GFP goes in at twice the moles of the vector, which is what the kit asks for with two
fragments. The count is in moles rather than in mass because a short insert weighs less at the
same number of molecules. Pipette from whatever your own concentrations are; the two units are
there so you do not do that arithmetic yourself.

**The incubation** is 50 °C for 15 minutes with two fragments. It grows to an hour with four or
more. The number is the kit's own, and it changes when you name another kit.

**The screening step** is the one to read twice. A correct clone gives two bands, 160 and 879 bp.
An empty vector gives one band at 219 bp. A clone with GFP the wrong way round gives 220 and
879 bp. Read the big band first: only the two lanes that hold GFP carry it. The empty vector and
the flipped insert are 1 bp apart, so the small band alone will not tell them apart.

## Name the kit on your bench

Three kits are shipped as rows of data. Each carries its own overlap rule, reaction, incubation,
molar ratio and fragment count, each read from that supplier's own manual.

```bash
pixi run liulab_mbio cloning gibson plan vector.dna insert.dna --out plan/ --product in-fusion
```

The start of a name is enough: `nebuilder`, `gibson` or `in-fusion`. Naming one changes the
design and not just a label. Ask for In-Fusion and the overlaps come out at 15 bp rather than
16, the reaction is half the size, and the melting temperature check stops judging, because
Takara publishes no melting temperature for an overlap. Nothing is borrowed from one supplier to
fill a gap in another.

## Two jobs that need no PCR

A short part, such as a linker or a tag, can be built from oligos that overlap each other and
tile both strands. They go into the same tube as everything else:

```bash
pixi run liulab_mbio cloning gibson plan vector.dna linker.fasta --out plan/ --route stitch
```

Twelve oligos are the most the method allows, which tile 500 bases. Above that the command stops
and tells you to amplify instead. It also warns outside 60 to 150 bp, the window this route is
worth using in: below it a primer tail can carry the bases, and above it a PCR is cheaper.

Two fragments that share nothing can be joined by a single oligo carrying bases from each end.
Neither fragment then needs a tail, so an amplicon you already have goes in as it is:

```bash
pixi run liulab_mbio cloning gibson plan vector.dna insert.dna --out plan/ --bridge "GFP:pUC19 backbone"
```

Name the two parts in the order they go round the product. Both kinds of oligo land on the same
order sheet as the primers, each with the job it does. Neither one primes anything, so no
threshold judges it and its row carries no verdict.

## Before you order

- Read the whole protocol once, top to bottom.
- Look at every badge that is not green, and decide about it.
- Read the checks that say `not judged`, and decide about those too.
- Check that the product map holds what you meant to clone.
- Check the three screening lanes can still be told apart.

`pixi run liulab_mbio cloning gibson plan --help` lists the rest of the options: where the
inserts go, which way round each one goes, which kit and polymerase to use, and which strain the
protocol should name.

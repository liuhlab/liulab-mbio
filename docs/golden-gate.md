# Put GFP into pUC19

This page follows one cloning job from end to end. The vector is pUC19. The insert is the GFP
coding sequence. One command does the design and writes the bench protocol.

Both sequence files ship with the repo, under `tests/data/`.

## Run it

```bash
pixi run liulab_mbio goldengate plan tests/data/pUC19.dna tests/data/GFP.dna --out docs/examples/pUC19-GFP
```

It prints one summary line, then the three paths it wrote:

```text
pUC19-GFP: 3347 bp, BbsI, 2 fragments, overhangs ATGA, TGGC, fidelity 100% (measured), checks warn
```

Read that line first. It is the whole design in one sentence.

| Part of the line | What it says |
| --- | --- |
| `3347 bp` | how long the finished plasmid is |
| `BbsI` | the enzyme it chose |
| `2 fragments` | the opened vector, plus GFP |
| `overhangs ATGA, TGGC` | the four bases each join is cut to |
| `fidelity 100% (measured)` | how well those two joins should pair up |
| `checks warn` | something is worth a look before you order |

### Why BbsI

The tool ranks enzymes by how many sites they read in the parts. BsaI reads a site in both
pUC19 and GFP, and BsmBI reads two in pUC19. Cutting either would cut the parts as well as the
ends. BbsI reads none, so nothing has to be changed to use it.

If no enzyme were free, the command would stop and say so. It never quietly edits your
sequence to make one fit.

### Why `checks warn`

Nine oligos were designed and five carry a warning: two on GC content, two on the GC clamp and
one on melting temperature. None fails. The badge names those kinds, and the oligo
sheet gives every oligo a verdict of its own.

Each oligo is judged by what it is for. A sequencing primer may be as short as 16 bases, since
one of the free primers Genewiz itself offers is that long. A PCR primer needs 18. The GC clamp
asks for one to three G or C bases among the last five. No published rule sets that band, so
the sheet calls it proposed.

A warning is yours to judge. The design tries every length a primer could have and keeps the
one that warns least, so each of these comes from where that primer has to sit.

## What comes out

| File | What it is |
| --- | --- |
| `product.dna` | the finished plasmid, with every feature carried over and both joins marked. Opens in SnapGene |
| `primers.tsv` | the nine oligos, with length and melting temperature. This is your order |
| `protocol.html` | the bench protocol: one page, no network, nothing to install |

Nothing here is written by hand, so the same two input files always give the same three files.
That means you can run it again instead of editing the output.

## How to read the protocol

Open `protocol.html` in a browser. It is built to be followed at the bench, not filed away.

**The top of the page** carries eight cards: the vector, the insert, the enzyme, the fragment
count, the overhangs, the fidelity, the product and what to select on. Each is a few words.
Under them sit four sentences about what the clone should do, then the badges — `sites`,
`junctions`, `primers` and one per part. A badge that is not green spells out why beneath it.

Two of those sentences matter before you start:

- GFP goes in on the opposite strand from the lac promoter, so that promoter does not read it.
- Nothing is annotated ahead of GFP to start translation, so the clone is not expected to glow.

That is the design telling you what it is and is not. It reads both from the finished plasmid,
so it cannot disagree with the map.

**The oligo sheet** gives each oligo its own verdict, written as a word rather than a colour:
`pass`, `warn` or `fail`. One line under the sheet opens to show which check fired on which
oligo, what it measured, and the band it missed. It stays closed, so the sheet you order from
is still a list of oligos.

**The middle** is eleven numbered steps, from the first PCR to sequencing. Each step has a
check mark, what to do, what a good result looks like, and what to do when it is wrong.
Reaction tables rescale when you change the number of reactions at the top of the table, so
you do not do that arithmetic yourself. Gel steps draw the bands to expect.

**The screening step** is the one to read twice. A correct clone gives two bands, 160 and
897 bp. An empty vector gives one band at 236 bp. A clone with GFP the wrong way round gives
220 and 897 bp. Those three lanes are different on purpose, so the gel alone tells you which
you have.

## The finished example

The plan above is published here, exactly as the command wrote it:

- [the bench protocol](examples/pUC19-GFP/protocol.html)
- [everything the run wrote](examples/pUC19-GFP/index.md)

## Point it at your own ligase data

Fidelity says how likely the joins are to pair correctly. For five enzymes that number is
measured, and the measurements ship with the package. For the rest, it is worked out from
rules.

You can do better if you hold a copy of the T4 ligase data from Potapov and colleagues. That
archive is not ours to ship, so nothing from it is in this package. Point the command at your
own copy instead:

```bash
pixi run liulab_mbio goldengate plan vector.dna insert.dna --out plan/ \
    --ligase-matrix ~/potapov/FileS03_T4_18h_25C.xlsx
```

An environment variable does the same for every run:

```bash
export LIULAB_MBIO_LIGASE_MATRIX=~/potapov/FileS03_T4_18h_25C.xlsx
```

Reach for `FileS03` first: 25 °C for 18 hours is what NEB's own viewer assumes. Both `.xlsx`
and `.csv` are read.

The report then says where its number came from, and keeps the three kinds apart:

| Scored against | What the protocol prints |
| --- | --- |
| the enzyme's own measurements | `measured` |
| your ligase file | `measured ligase profile, not specific to PaqCI` |
| the rules | `rule-based estimate` |

The enzyme's own measurements win where they exist. A ligase file stands in for a measurement
nobody has made; it does not replace one somebody has. Without such a file nothing changes.

## Before you order

- Read the whole protocol once, top to bottom.
- Look at every badge that is not green, and decide about it.
- Open the line under the oligo sheet and decide about each warning.
- Check that the product map holds what you meant to clone.
- Check the three screening bands are different from each other.

`pixi run liulab_mbio goldengate plan --help` lists the rest of the options: where to put the
insert, which way round it goes, holding a join in frame, forcing an enzyme, and which
polymerase and strain the protocol should name.

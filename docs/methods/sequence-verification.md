# Sequence verification

A sequencing result tells you what a clone is. This command tells you whether it is the clone
you planned. It lines up one clone's results against the plasmid it should be, and gives each
junction and each insert a verdict. Every base that differs is placed on the plasmid's map.

## When to use it

- Your Sanger reads of a miniprep came back, and you want to know which colony to keep.
- A sequencing service sent a whole-plasmid consensus, and you want every difference placed.
- A cloning protocol from this package reached its last step, which checks each read against
  the planned plasmid.

## What you need

The plasmid the clone should be. Each of the four cloning plans writes it as `product.dna`, and
marks each junction on it. The command reads those marks, and judges each junction and each
stretch between two of them as an insert. The vector is left out, and
[the command line reference](../reference/cli.md#check-a-clones-sequencing) says how it is found.

Every result one sample gave. A Sanger read is its `.ab1` file. A whole-plasmid consensus is the
FASTA or GenBank file the service sent. Give a read from each side of the insert together, in
one run: a base counts as read when any result covers it.

## Run it

The example is the Golden Gate plan's GFP in pUC19. No real clone of it was sequenced, so the
result below was made from the plan's own plasmid with two bases changed: one in GFP, and one
in AmpR, the ampicillin resistance gene in the backbone.

```bash
pixi run mbio sequence-verify docs/examples/pUC19-GFP/product.dna \
  docs/examples/pUC19-GFP/verification/consensus.fasta
```

It prints:

```text
ATGA junction    pass        every base read, and none disagrees
GFP insert       fail        substitution at 696
TGGC junction    pass        every base read, and none disagrees
substitution at 2687, in AmpR: outside every region, so no verdict
consensus.fasta  pass        reads 1 .. 3347 on the forward strand
consensus.fasta: a consensus alone cannot show a mixed sample, because the commonest plasmid becomes the consensus
not verified
```

The first three lines judge the two junctions and GFP between them. Both junctions read true.
GFP fails at base 696, where the consensus reads another base than the plan. Positions count
from 1, as a map shows them.

The fourth line is the AmpR change at base 2687. It lies outside every junction and insert, so
it is listed with the feature it falls in and judges nothing.

The next two lines are the result's own. The first says where it landed: the whole plasmid, on
the forward strand. The second says what a consensus cannot show.

The last line is the clone's. GFP failed, so the clone is not verified, and the command exits
with 1.

| File | What it is |
| --- | --- |
| [product.dna](../examples/pUC19-GFP/product.dna) | the plasmid the Golden Gate plan says to make |
| [consensus.fasta](../examples/pUC19-GFP/verification/consensus.fasta) | the same bases with two changed, standing in for what a sequencing service sends |

## What each verdict means

| Verdict | On a junction or an insert |
| --- | --- |
| pass | every base was read, and none differs from the plan |
| fail | a base differs and every read covering it agrees, or the base or the whole sample is mixed |
| warn | two reads disagree with each other at one base |
| no verdict | part of it went unread, and nothing failed |

A clone is verified only when every junction and every insert passes. Nothing failing is not
enough. A junction no read reached has no verdict, and the clone stays unverified until a read
covers it.

In the example, the GFP change alone decides that.

## A difference outside the junctions and inserts

The AmpR change is listed with the feature it falls in, and it gets no verdict. Whether it
matters depends on what that feature does, and the command cannot know that.

Read it three ways before you discard the clone:

- The parent plasmid may never have matched its map. Addgene finds errors in 30% of the
  plasmids deposited with it, and often a few mismatches in the origin of replication or other
  backbone parts.
- The sequencer may have made it. On a nanopore consensus, a methylated site or a run of more
  than nine of one base is a known error.
- It may be real and harmless. Addgene counts a change as a problem only when it may affect
  what the plasmid does.

Sequencing the parent beside the clone settles the first. A change the parent shares was not
made by your cloning.

## Sanger reads and a whole-plasmid consensus

**A Sanger read** is trusted only in its middle. The first bases after the primer and the tail of
the read are poor, so each read is trimmed by its own quality scores. Inside what is left, a base
counts at quality 20 or better, a 1 in 100 chance of a wrong call. A base whose second peak
reaches a quarter of the first is mixed, and a mixed base fails.

**A whole-plasmid consensus** is trusted whole. The service has already lined up its reads and
called each base. But a consensus cannot show a mixed sample, because the commonest plasmid
becomes the consensus. The command says so each time it reads one. Check the service's
per-base table or trace when a mixture matters.

## A plasmid no plan wrote

A record you made yourself has no junctions marked. Name what to judge with `--feature`, once
for each feature:

```bash
pixi run mbio sequence-verify my-plasmid.gb read-F.ab1 read-R.ab1 --feature GFP
```

With no marked junction and no `--feature`, every difference is listed and nothing is judged.

## Tips and troubleshooting

**A junction has no verdict.** No trusted base reached it. The first bases after a primer are
unreadable, so sequence again from a primer further out.

**A result does not read as this plasmid.** More than a tenth of the bases it trusts disagree
or fail to line up. It is a different plasmid, or another sample's file, and it judges nothing.

**Two reads disagree.** Neither outvotes the other, so the verdict is warn. Read across that
base again.

`pixi run mbio sequence-verify --help` lists the rest.

## Where the numbers come from

- Quality 20 is the line Applied Biosystems' analysis software marks acceptable, and the one
  GENEWIZ counts a read's length by.
- A second peak at a quarter of the first is the top of the 15 to 25% range Applied Biosystems
  recommends for calling a mixed base, so a clean trace calls the fewest.
- Where a result gives the number of reads behind each base, a base needs more than 10.
  GENEWIZ fails a sample at 10 or fewer. Below 20 it warns, because Plasmidsaurus and Eurofins
  call about 20 enough for an accurate consensus.

## References

- Applied Biosystems, *Sequencing Analysis Software v5.4 Quick Reference Card*, PN 4401738
  Rev. B, April 2009.
- Applied Biosystems, *DNA Sequencing by Capillary Electrophoresis Chemistry Guide*, PN 4305080,
  chapter 8.
- GENEWIZ, Sanger sequencing service page, undated.
- Azenta/GENEWIZ, *Plasmid-EZ* FAQ, undated.
- Plasmidsaurus, *Technical Documentation: Plasmid*, PS-0002-E v1.2, 22 September 2026.
- Plasmidsaurus, technical note on Oxford Nanopore error modes, 12 September 2025.
- Eurofins Genomics, *Results Interpretation Guide* for nanopore sequencing, undated.
- Addgene blog, A. Hazen, *Addgene moves to NGS verification*, 3 August 2017, updated 6 April
  2021.
- Addgene blog, A. Shepard, *A Look at Addgene's QC Process*, 6 May 2025.

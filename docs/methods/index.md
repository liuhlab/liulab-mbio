# Choosing a method

Four ways of joining DNA are planned here. They differ first in what the junction costs you.
Golden Gate and Gibson leave the join reading what you chose, base for base. Restriction and
ligation puts back the site the enzyme binds, at both ends, because that is where the two ends
came from. Gateway leaves a whole att site at each end, which a fusion protein has to read
through. Decide that first, then count your fragments.

| Method | Fragments | What the junction costs | What the parts must carry | Bench cost |
| --- | --- | --- | --- | --- |
| [Golden Gate assembly](golden-gate.md) | two to about a dozen, in one tube | nothing: the overhang is yours to choose | no site for the Type IIS enzyme chosen, or a silent codon change to clear one | 9 oligos, 11 steps |
| [Gibson assembly](gibson.md) | up to about five, in one tube | nothing: the overlap is homology, not added bases | nothing; a primer tail carries the overlap | 9 oligos, 11 steps |
| [Restriction and ligation](restriction-ligation.md) | one insert into a vector | a recognition site at each end, GCATGC and AAGCTT in this example | a usable pair of sites, or primers that add them | 6 oligos, 12 steps |
| [Gateway cloning](gateway.md) | one insert into a destination | a whole 25 bp att site at each end | att sites, or a PCR that adds them | 6 oligos, 13 steps |

Those counts all come from one job planned four ways: GFP into pUC19, from the same two files, so
the four pages compare directly. Golden Gate and Gibson are one tube and one incubation after the
PCRs. Restriction and ligation adds a preparative gel to recover the cut fragments, and a second
digest to check the miniprep. Gateway runs two reactions on two days with a miniprep between
them, and needs both vendor enzyme mixes and a ccdB-surviving strain to grow the vectors in.

Each of the four writes the plasmid you should get, the oligos to order, and the bench protocol
as one page you can follow at the bench. Gateway also writes the entry clone that the BP reaction
makes.

## The other five pages

[Primer design](primers.md) designs or judges a pair on its own, for a fragment you only want to
amplify, and can check that pair against a whole genome before you order it.

[Barcode sets](barcodes.md) draws a set of barcodes far enough apart that no two read as one
another, free of the enzyme sites you forbid.

[Codon optimisation](codon-optimisation.md) writes a protein as DNA in a host's commonest codons,
and moves the codons that would spell an enzyme site you named.

[Sequence verification](sequence-verification.md) checks a clone's sequencing against the
plasmid it should be, once the reads come back, and gives each junction and insert a verdict.

[Maps and figures](maps.md) draws any record as a figure or as a page to look through, including
every product the four methods above write.

# Example: an AP-1 domain library

A combinatorial library over three positions, built from 72 protein sequences and nothing else.
Every output file below was written by the command, from the three input files beside it.
Nothing here is edited by hand.

```bash
pixi run liulab_synbio library plan docs/examples/ap1-library/parts.fasta \
  --scheme docs/examples/ap1-library/scheme.json \
  --vector docs/examples/ap1-library/vector.gb \
  --host human --coverage 300 --name 'AP-1 DESynR' --out docs/examples/ap1-library
```

To change what these files say, change the code and run that command again. The same inputs
write the same bytes, so a run that changes nothing leaves them alone.

## What goes in

| File | What it is |
| --- | --- |
| [parts.fasta](parts.fasta) | 72 proteins: 24 for each of the N, DBD and C positions |
| [scheme.json](scheme.json) | the architecture: which enzyme does which job, and the stuffer each part carries |
| [vector.gb](vector.gb) | the destination the first round opens |

The names say the position: `N_JUN` fills N, `DBD_JUN` fills DBD, `C_JUN` fills C. No overhang,
stuffer, barcode or codon is given. The planner chooses all four.

## What comes out

| File | What it is |
| --- | --- |
| [parts.tsv](parts.tsv) | the synthesis order sheet: 72 blocks, 133 to 1,149 bases |
| [barcodes.tsv](barcodes.tsv) | which barcode names which part, and where it sits |
| [changes.tsv](changes.tsv) | every amino acid the overhang standard moved, wild type beside synthesised |
| [round-1.dna](round-1.dna), [round-2.dna](round-2.dna) | one annotated record a round |
| [product.dna](product.dna) | one member of the finished library, 2,276 bases |
| [protocol.json](protocol.json) | the bench protocol as data |
| [protocol.html](protocol.html) | the same protocol as a page to work from |

The three part lists make 24 x 24 x 24 = 13,824 members. The product file holds one of them,
with the rest differing only in which protein and which barcode sits at each position.

## Where the proteins come from

The 72 sequences are the AP-1 domains of Takacsi-Nagy, O. et al. (2026) Synthetic transcription
factors designed by domain recombination enhance CAR T cell antitumor function. *Cell* 189,
1-20, [doi:10.1016/j.cell.2026.07.054](https://doi.org/10.1016/j.cell.2026.07.054), used under
[CC BY 4.0](http://creativecommons.org/licenses/by/4.0/). They were read from that paper's
Table S1 and translated back to protein; the DNA here is this package's own.

[Highly parallel DNA synthesis and assembly](../../synthesis-and-assembly.md) explains the
method these files follow.

# Barcode sets

A barcode is a stretch of DNA that names one part of a library. The set has to hold three
things: no two barcodes may read as one another, none may spell a site the cloning enzyme
reads, and none may put a stop codon in the frame the construct is read in. One
command designs a set that holds all three, and says how far apart it stands.

## When to use it

- You are naming the parts of a pooled or combinatorial library, one barcode each.
- The parts will be cut and joined, so no barcode and no join may spell a site the enzyme
  reads.
- The barcode is translated, so it must not put a stop in the frame.

[Golden Gate assembly](golden-gate.md) joins the parts that carry these barcodes, and
[codon optimisation](codon-optimisation.md) writes the parts themselves free of the same
sites.

## What you need

No input file. You give the numbers: how many barcodes, how long each one is, the cloning scar
that joins one to the next, and every enzyme to keep out. The run below asks for the set a
small barcoded library would use.

## Run it

```bash
pixi run mbio barcode-design 24 --length 11 --scar AGCG --forbid BsaI --forbid BbsI --out plan/
```

It prints the set, then the two files it wrote:

```text
24 barcode(s) of 11 bases, scar AGCG, 3 apart by sequence-levenshtein, free of BsaI, BbsI, seed 0, checks pass
plan/barcodes.tsv
plan/checks.tsv
```

Candidates are drawn from the whole space of 11-base sequences in a shuffled order and kept
when they pass every rule. Walking that space in order instead opens a set on a run of one
base, which no later filter repairs. Each candidate is read with the scar on both sides of it,
so a BsaI or BbsI site that only shows up across a join is caught as well.

`3 apart by sequence-levenshtein` is the rule the set holds to. Sequence-Levenshtein is not
plain Levenshtein: it counts substitutions, insertions and deletions, and charges nothing for
the bases a read gains or loses past the barcode's end, which is what a lost base does to a
read. `checks pass` means nothing failed.

## What it wrote

| File | What it is |
| --- | --- |
| [barcodes.tsv](../examples/design/barcodes.tsv) | the 24 barcodes, numbered in the order they were accepted |
| [checks.tsv](../examples/design/checks.tsv) | the verdict on the set, one row per check |

Both names link to what that run wrote, published here unedited. The same count, the same
rules and the same seed always give the same 24 barcodes, so you can store the command rather
than the set. To change the set, change the command and run it again.

## How to read the sheet

`barcodes.tsv` holds the number and the bases. The first five here are TACCATTGTCT,
TCCTCCCCAGG, ACCAGTGGCCT, GACAGCAAGTT and TTGATAGTGAA. Keep that numbering, because it is what
pairs a barcode with the part it names.

`checks.tsv` holds a row for each check: its name, its verdict, the number it measured and a
line of detail. `rules` found nothing broken. `separation` is 3, which is how far apart the
set really stands, against a rule of 3.

`deletion_ambiguity` has no verdict, on purpose. It measures the share of one-base deletions
that another barcode could have left too, which is zero here, and then stops, because no
sourced threshold says what share is too much. A check nothing sourced judges reports its
number and leaves the decision with you.

## Before you order

- Check the count covers the parts you have and the ones you may add.
- Check `--scar` is the four bases your own cloning leaves behind.
- Check `--phase` is the frame the construct reads the barcode in, or pass `--untranslated`.
- Name every enzyme you will cut with, including one you only use in a later round.
- Read `deletion_ambiguity` and decide whether that share is one you can live with.

## Tips and troubleshooting

**The design gives up before it has your set.** The space ran out under the rules you asked
for, and the message says how many it found and which rule rejected the most candidates. Add a
base to `--length`, or drop `--distance` by one.

**Nothing translates your barcode.** Pass `--untranslated`, and the stop-codon rule and the
whole-codon rule both go away. Otherwise `--phase` says how many bases of the barcode's first
codon are read before it: 0, 1 or 2.

**You want the same set again next year.** One seed gives one set, so keep the command.
Another `--seed` gives a different set of the same quality.

## Where the numbers come from

The default distance of 3 is the distance a published barcode set was measured to hold.
Sequence-Levenshtein is the default metric because the long reads a barcode block is read with
fail mostly by insertion and deletion. Hamming counts only the positions two barcodes of one
length differ in, so it cannot see an indel at any distance.

The homopolymer cap of 5 bases is on for the same reason, since a run is where an indel
lands. That number is a judgement rather than a measurement. A GC band is off, and
`--gc-band` turns one on: every band in the barcode literature traces back to one uncited
sentence in Hamady et al. 2008, and Xu et al. 2009, the one study that imposed a band and
then measured the result, found GC was not what separated its probes.

A set that is ligated has its junctions judged the way a Golden Gate overhang set is; see
[ligation fidelity](../reference/ligation-fidelity.md). `pixi run mbio barcode-design --help`
lists the rest.

## References

- Buschmann, T. and Bystrykh, L.V. (2013) Levenshtein error-correcting barcodes for
  multiplexed DNA sequencing. *BMC Bioinformatics* 14: 272.
- Hamady, M., Walker, J.J., Harris, J.K., Gold, N.J. and Knight, R. (2008) Error-correcting
  barcoded primers for pyrosequencing hundreds of samples in multiplex. *Nature Methods* 5,
  235–237.
- Xu, Q., Schlabach, M.R., Hannon, G.J. and Elledge, S.J. (2009) Design of 240,000 orthogonal
  25mer DNA barcode probes. *PNAS* 106, 2289–2294.
- Hawkins, J.A. et al. (2018) Indel-correcting DNA barcodes for high-throughput sequencing.
  *PNAS* 115, E6217–E6226.

# Primer design

A PCR primer has to bind one place on the template, melt near its partner, and end on a base
the polymerase can extend. You say how far each primer may move; the package judges every
binding site in that room and keeps the best pair. One command does it and writes the sheet
you order from.

## When to use it

- You want every check on a pair before you order it.
- A primer has to carry a tail: a restriction site, a Type IIS site, or bases that hold a
  fusion in frame.
- You are amplifying genomic DNA, where a pair clean on a plasmid may prime in twenty places.

You do not run this for the cloning pages: [Golden Gate](golden-gate.md),
[Gibson](gibson.md), [restriction and ligation](restriction-ligation.md) and
[Gateway](gateway.md) design their own primers.

## What you need

One record to amplify from, as SnapGene `.dna`, GenBank or FASTA, and the region you want; on
a genome, that FASTA and its index. The run below uses a file that ships with the repo.

## Run it

```bash
pixi run mbio primers design tests/data/pUC19.dna --region 20..340 --out plan/
```

It prints the pair, then the file it wrote:

```text
pUC19: CGGTGAAAACCTCTGACA / AATCGCCTTGCAGCACAT, 321 bp, Ta 63 °C, checks warn
plan/primers.tsv
```

`--region 20..340` counts from 1 and takes both ends in; every verb here reads a span that
way. Told nothing more, the design pins each primer's 5' end to one end of the span and varies
only its length, so the amplicon is the 321 bp you asked for. Tell it more and the primer
may move: a tail to anchor it, a junction to read across, a flank to sit inside. It judges
every site that room allows and keeps the best pair, so a warning left on a primer is one
nothing nearby could clear.

`Ta 63 °C` is the temperature to anneal at. `checks warn` means nothing failed and something
wants a look.

## What it wrote

| File | What it is |
| --- | --- |
| [primers.tsv](../examples/design/primers.tsv) | the two oligos to order, with length and melting temperature |

That name links to what the run wrote, published here unedited. The same record and region
always give the same pair, so changing the design means changing the command:
`--target-tm` moves the temperature aimed for, 62 °C, and `--polymerase` picks whose rules
judge the pair, Q5.

## Check the pair on a genome

Give it a genome FASTA and the assembly's name, and it designs against the whole genome:

```bash
pixi run mbio primers design GENOME.fa --assembly ce11 --region chrI:401..700 --out plan/
```

The region now names the sequence it lies on, as `NAME:START..END`. The forward primer goes in
the left flank and the reverse in the right, so the amplicon covers the whole region;
`--flank` says how far out either may sit. The best-ranked pairs go into one genome search,
and the first making only the amplicon you asked for wins. Where none does, the primers that
primed elsewhere are ruled out and it searches again, and after a fixed number of searches it
hands back the best pair it saw, marked not specific.

The summary line gains the genome and the region in front; a second line follows it, with
whether the pair is specific or how many other amplicons it makes, how many searches that
took, and the worst of its checks. Beside the sheet it writes `genome.tsv`, a row for each
amplicon on the genome: its place in the spelling you typed, its length, whether the pair or
one primer alone makes it, the mismatches at each end, and whether it is yours.

A large FASTA gets a perfect-match screen: above a size limit a site with one mismatch goes
unseen, and the report says so. `ipcr` runs the search. It comes from bioconda, not PyPI, so
`pixi install` brings it and a pip install does not.

No genome FASTA ships here and none is downloaded; the lab's genomes live in liulab-genome, so
ask it for a path. It is the one section here with no committed output behind it, because a
genome cannot be committed.

## How to read the sheet

Open `primers.tsv` in a spreadsheet or an editor. One row per oligo: the name you order it
under, the bases 5' to 3', how many there are, and the melting temperature. Both oligos here
are 18 bases, at 61.7 and 65.3 °C.

That temperature is the binding part's, not the whole oligo's. The annealing temperature on
the summary line comes off the lower of the two, by the polymerase maker's own rule.

The sheet carries no verdicts; the summary line's last words hold the worst of every check.
Two checks are never judged: the melting temperature of the whole oligo, and how tightly its
last bases hold. Nothing sourced says what a good value is.

## Before you order

- Read the summary line: both sequences, the amplicon, the annealing temperature.
- Check the amplicon is the region you meant, both ends counted in.
- Check each tail reads the way the enzyme needs it.
- Set the block to that annealing temperature.
- On genomic DNA, read `genome.tsv`: the only amplicon should be yours.

## Tips and troubleshooting

**The pair warns and you want a cleaner one.** More room is the only thing that helps. Move
the ends of `--region` a few bases, or aim at another temperature with `--target-tm`.

**Your primer needs a tail.** `--forward-tail` and `--reverse-tail` join bases to a 5' end.
The tail stays out of the search for a binding site, and it still lengthens the amplicon.

**No index sits beside the genome FASTA.** The region is read through that index, so the run
stops without it. Ask liulab-genome for an indexed copy.

## Where the numbers come from

Every check is held to a band with a source: primer3's defaults for the hairpin and dimer
limits, IDT's ranges for length and melting temperature, NEB's calculator for the gap
between a pair's two temperatures. Where no source states a band, the package either proposes
one and says so, or judges nothing.

Melting temperature is the nearest-neighbour model, computed with primer3 and salt-corrected
for the polymerase's own buffer, so it agrees with NEB's Tm Calculator. The annealing rule and
the extension rate are that polymerase's, and the off-target rules are Primer-BLAST's
defaults.

`pixi run mbio primers design --help` lists the rest.

## References

- SantaLucia, J. (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA
  nearest-neighbor thermodynamics. *PNAS* 95, 1460–1465.
- Ye, J. et al. (2012) Primer-BLAST: a tool to design target-specific primers for polymerase
  chain reaction. *BMC Bioinformatics* 13: 134.
- NEB, Tm Calculator, and the Q5 High-Fidelity DNA Polymerase protocol (NEB #M0491).

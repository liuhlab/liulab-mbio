# Check primers on a genome

The package checks every primer pair against the template you hand it. For a plasmid or a
fragment that is the whole story, because the PCR never sees anything else.

Genomic DNA is not like that. The primers meet the whole genome, and a pair that is clean on a
3 kb plasmid can prime in twenty other places on a 3 Gb one. So check a pair against the genome
before you order it.

There is no command for this yet. It is Python.

## Find the FASTA

The lab's genomes live in liulab-genome. Ask it where one is:

```bash
genome assembly files hg38 --json
```

That prints one line of JSON and never downloads anything. Three fields matter:

| Field | What it is |
| --- | --- |
| `genome_files.fasta` | the path to the FASTA |
| `genome_files.fai` | the path to its index, which sits beside the FASTA |
| `assembly` | the name to pass along, so a report says which genome it read |

An assembly that is not registered, or not trusted, prints nothing, exits 1, and names
`genome assembly register` on the error stream.

This command is merged on liulab-genome's main branch. It is in no tagged release yet.

## The three calls

| What you have | Call |
| --- | --- |
| one pair | `evaluate_pair_on_genome` |
| several pairs | `evaluate_on_genome` |
| a region you want a pair for | `design_pair_on_genome` |

One search over many pairs costs little more than one search over a single pair. So pass your
pairs together where you can. Each pair is still its own reaction, and primers of different
pairs are never mixed.

`design_pair_on_genome` is the one that designs. It puts the forward primer in the left flank
of your region and the reverse in the right, then checks the best-ranked pairs in one search.
The first pair that makes only the amplicon you asked for wins. Where none does, it rules out
the primers that primed elsewhere and searches again. After `Thresholds.genome_rounds` searches
it gives up, and hands back the best pair it saw, marked not specific.

The docstrings carry the worked calls and every argument. They are on the
[API reference](api.md) page.

## Read the result

A `GenomeReport` says what a pair makes:

- `off_target` lists every amplicon but the one you meant. Each names its sequence, where it
  starts and ends, how long it is, and whether the pair made it or one primer made it alone.
- `status` is `pass`, `warn` or `fail`, whichever is worst among the checks.
- `near_matches_checked` is false where the search looked for perfect matches only.

A `GenomeDesign` adds `specific`. False means the search gave up. The pair it hands back is
then the best one it checked, not a clean one, so read its off-target amplicons before you
order.

## A large genome gets a perfect-match screen

Looking for near matches is quick on a small genome and slow on a large one. So the size of the
FASTA on disk decides what the search does:

| FASTA size | What the search finds |
| --- | --- |
| up to `Thresholds.genome_size_limit` | sites with up to `Thresholds.genome_mismatches` mismatches, none of them in the last `Thresholds.genome_terminal_window` bases |
| above it | perfect matches only |

Above that limit the check is a screen, not a proof. A site with one mismatch goes unseen, so
every report from such a search says near-match sites were not checked.

## The search tool

`ipcr` does the searching. It comes from bioconda and is not on PyPI, so `pixi install` brings
it and a pip install does not. Where it is missing, the check says so and tells you to run
`pixi install`.

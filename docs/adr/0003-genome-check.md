---
search:
  exclude: true
---

# The genome check runs a local in-silico PCR tool on a FASTA path

A primer pair is checked against a genome by `ipcr`, run on a FASTA the lab already holds,
behind one function in `primers.genome`. The tool streams the FASTA with no index, pairs binding
sites into amplicons, counts those one primer makes alone, and holds a 3'-terminal window to a
perfect match. The package's own Tm rule then drops a site too weak to prime. The caller passes
the FASTA path and the assembly name, both of which liulab-genome gives.

On hg38, `ipcr` checked one pair in 7 s and 35 pairs in 10 s with no mismatches, but one
mismatch took over 3 minutes for one pair (#44). On the 100 Mb worm genome, 3 mismatches with a
3-base window took 0.3 s (#47). So a genome up to the size limit is searched with mismatches,
and a larger one with perfect matches only.

## Considered options

- **NCBI Primer-BLAST.** It documents no API. The BLAST web API is rate-limited, sends the
  sequences off the machine, and never pairs hits into amplicons.
- **Local BLAST, pairing its hits ourselves.** Another engine and another index per genome, and
  strands, pairing and the origin rewritten here.
- **`dicey`**, an FM-index search. Its hg38 index build passed 12 minutes and 22 GB before it
  was stopped (#44), and it is GPL-3.
- **Depending on liulab-genome** to find the FASTA. It would tie sequence design to a download
  cache, where a path is all the check needs.

## Consequences

Over the size limit, the check is a perfect-match screen: a site with one mismatch goes unseen,
and each report says near-match sites were not checked. `ipcr` keeps at most 10,000 sites per
primer on one sequence, so a primer in a common repeat lists only some of its off-target
amplicons. It ignores SIGTERM, so a run past its time limit is killed outright. It misses a
binding site spanning the origin, so the check doubles each circular record itself. It is a
bioconda package absent from PyPI: a pip install lacks it, and the check says how to get it.

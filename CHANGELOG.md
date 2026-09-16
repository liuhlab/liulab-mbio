# Changelog

Every change worth knowing about, newest first. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Version numbers are CalVer tags
of the form `vYYYY.M.PATCH`, and the tag is where the version comes from — nothing here
sets one.

## [Unreleased]

### Added

- A Golden Gate plan picks the codons of the host the part is for: `plan_assembly` takes
  `codon_table`, `goldengate plan` takes `--codon-table`, and a plan that has to propose
  domestication names the codons it would move to. The strain the protocol transforms is still
  `--host`, and still *E. coli*.
- A barcoded library of protein combinations, built in rounds rather than one pot.
  `liulab_mbio.library.plan_library` takes lists of proteins, or coding DNA, with a scheme and a
  destination vector. It picks the overhang standard that costs the proteins fewest changed
  residues, writes each part's synthesis sequence with its own barcode, fits a vector that cannot
  be opened yet, simulates every round, and counts the colonies a round needs for the coverage
  asked for. `LibraryPlan.write` puts the synthesis order sheet, the barcode table, the
  amino-acid change table, a record for each round, the assembled product, the protocol as JSON
  and the page rendered from it in one directory.
- `liulab_mbio library plan` on the command line, and a repo-local `protein-assembly` skill that
  calls it.
- Writing DNA for a protein, and choosing its codons for a host while clearing sites it must not
  spell (`liulab_mbio.translate`); and barcode sets held a set distance apart
  (`liulab_mbio.barcodes`). Each has its own skill, `codon-optimize` and `barcode-design`,
  because both are wanted outside a library build.
- SrfI and PmeI join the shipped enzymes, rebuilt through the enzyme-data builder.
- Two research notes: the library method with every number sourced to the paper it comes from,
  and whether a barcode set needs limits on GC and on repeated bases. The second measured the
  evidence rather than following custom, and the answer is a cap on repeated bases and no GC
  band at all.
- Golden Gate cloning, end to end. `liulab_mbio.goldengate.plan_assembly` takes a vector and
  any number of inserts, picks a Type IIS enzyme with no site in the parts, designs the whole
  overhang set, checks every primer, simulates the assembly, and designs the colony PCR and
  sequencing that confirm the clone. `Plan.write` writes the annotated product, the primer
  order sheet, the protocol as JSON data (`write_protocol`), and an interactive HTML bench
  protocol rendered from that data, whose order sheet carries each oligo's verdict, and says
  which check fired and what it measured.
- `liulab_mbio goldengate plan` and `liulab_mbio protocol render` on the command line, and a
  repo-local `golden-gate-assembly` skill that calls them.
- SnapGene `.dna` read and write, with editing that carries features and primer binding sites,
  under one 0-based half-open coordinate model shared by every module.
- Restriction enzyme, codon usage and ligation fidelity data, each rebuilt by a script in
  `scripts/` and sourced in a note under `docs/research/`.
- Codon usage for human and mouse, beside *E. coli* K-12: `codon_usage("human")` and
  `codon_usage("mouse")`. Each counts one coding sequence per protein-coding gene, the one
  GENCODE marks as canonical, on hg38 and mm39.
- Primer design and evaluation against NEB's published rules, per polymerase.
- Fidelity scoring against a ligase matrix the user holds, by `--ligase-matrix` or
  `LIULAB_MBIO_LIGASE_MATRIX`. Nothing from that archive ships here.
- Docs: a walkthrough of GFP into pUC19, the plan it writes published as a live example, and
  an API reference covering every public module.

### Changed

- Barcode sets are designed on an indel-aware distance by default. `liulab_mbio.barcodes` takes
  the metric as a dial: `sequence-levenshtein`, which counts a lost or gained base, or `hamming`,
  which counts mismatches and cannot see one. Measured over the 11-mer space, up to 46% of the
  single deletions of a Hamming set read as another barcode of the same set, against none of an
  indel-aware one; the indel-aware rule costs about one base of barcode length at a given set
  size, rejects nothing in the published set, and every set it designs is a Hamming set too.
  `deletion_ambiguity` reports that share for any set, and a library protocol now prints it where
  the barcode block is read back. A library plan draws different barcodes than it did, the same
  seed still returning the same ones. `docs/research/barcode-design.md` section 10 is the
  measurement and the reasoning.
- Each designed oligo is judged by what it is for, and `plan_assembly` takes its thresholds for
  each role. A sequencing primer passes from 16 bases, the length of Genewiz's own M13F primer,
  while a PCR primer still needs 18. A design picks the primer length that warns least on
  length, GC, GC clamp and Tm. A band no published rule sets, such as the GC clamp, reads as
  proposed on the order sheet.

### Removed

- The template's placeholder `greet()`, its CLI verb and its test.

# Changelog

Every change worth knowing about, newest first. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Version numbers are CalVer tags
of the form `vYYYY.M.PATCH`, and the tag is where the version comes from — nothing here
sets one.

## [Unreleased]

### Added

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
- Primer design and evaluation against NEB's published rules, per polymerase.
- Fidelity scoring against a ligase matrix the user holds, by `--ligase-matrix` or
  `LIULAB_MBIO_LIGASE_MATRIX`. Nothing from that archive ships here.
- Docs: a walkthrough of GFP into pUC19, the plan it writes published as a live example, and
  an API reference covering every public module.

### Changed

- Each designed oligo is judged by what it is for, and `plan_assembly` takes its thresholds for
  each role. A sequencing primer passes from 16 bases, the length of Genewiz's own M13F primer,
  while a PCR primer still needs 18. A design picks the primer length that warns least on
  length, GC, GC clamp and Tm. A band no published rule sets, such as the GC clamp, reads as
  proposed on the order sheet.

### Removed

- The template's placeholder `greet()`, its CLI verb and its test.

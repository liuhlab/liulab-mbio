---
name: golden-gate-assembly
description: >-
  Plan a Golden Gate cloning experiment end to end from a vector and an insert sequence file
  (SnapGene .dna, GenBank or FASTA): pick a Type IIS enzyme, remove internal sites, design and
  check PCR primers, simulate the assembly, design colony PCR validation, and write an
  interactive HTML bench protocol with expected results. Use when the user wants to clone,
  insert or subclone a sequence into a plasmid by Golden Gate, asks for primers with BsaI,
  BsmBI, BbsI, PaqCI or SapI tails, or wants a ready-to-run cloning protocol from two
  sequence files.
---

# Golden Gate assembly

> **Scaffold.** The `liulab_mbio` functions this skill needs do not exist yet; issue #1
> tracks them. Until they land, tell the user the skill is not usable and point to #1.
> Do not hand-design primers or protocols in its place.

## Inputs

- **Vector**: a sequence file, usually circular.
- **Insert**: a sequence file.
- **Optional**: where to insert (a feature name or coordinates), orientation, in-frame
  fusion, preferred enzyme, polymerase, host strain. Ask only for what the files do not
  settle.

## Outputs

One output directory holding:

- the annotated product as `.dna`;
- a primer order sheet;
- one self-contained interactive HTML protocol, with expected results at every step.

## Workflow

Each step binds to package functions once #1 lands.

1. Read both files; report topology, length and features.
2. Scan both for Type IIS sites. Choose an enzyme with no site in the parts kept in the
   product; domesticate only when none is free.
3. Choose junctions and overhangs; score the overhang set.
4. Design vector and insert primers with tails; evaluate every primer and pair.
5. Simulate PCR, digestion and ligation. Confirm the product has no remaining site and
   carries the insert intact.
6. Design colony PCR and sequencing primers; compute expected bands for correct and empty
   clones.
7. Render the protocol: reagents, reaction tables, thermocycler programs, checks and
   troubleshooting per step.

## Smoke test

`tests/data/pUC19.dna` and `tests/data/GFP.dna`: insert GFP into the pUC19 multiple cloning
site and validate the insertion by colony PCR.

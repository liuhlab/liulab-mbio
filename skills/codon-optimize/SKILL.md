---
name: codon-optimize
description: >-
  Write a protein as a synthesis-ready coding sequence for a named host with `liulab_mbio`:
  reverse-translate it with the codons that host counts most often, take out every restriction
  or Type IIS site you name by synonymous codon change, and report every codon it moved. A
  sequence already coded is checked rather than written again. Use whenever someone wants a gene
  codon-optimised for a host, a protein turned into DNA to order from a vendor, an enzyme site
  removed from a coding sequence without changing the protein, or asks which codons a construct
  must change to stay free of BsaI, BbsI or another enzyme.
---

# Codon optimise

`liulab_mbio.translate` writes the DNA and checks it. Your job is to name the host and the
enzymes, and to read back what it changed. Never write a codon, a sequence or a site position
from your own knowledge: the package works each one out from a shipped codon usage table, and
nothing checks a base you invented (`docs/adr/0002-editable-protocols.md`).

## Name the host — there is no default

Each shipped table is counted over one genome. Ask the user which host the construct expresses
in, and pass that table's name. Never guess a host, and never write a codon table: a table
comes from `liulab_mbio.codons` or it does not exist.

```bash
pixi run python -c "from liulab_mbio.codons import codon_tables; print(codon_tables())"
```

## Run it

```bash
pixi run liulab_mbio codon-optimize PROTEIN --kind protein --host HOST --forbid BsaI --forbid BbsI
```

`--kind` says what you passed, `protein` or `dna`, and it is required: A, C, G and T are all
amino-acid letters too, so nothing can tell the two apart for you. `--forbid` goes once per
enzyme. It prints a summary line, one line per codon it changed, and the sequence last.

## From Python

```python
from liulab_mbio.translate import optimize_coding_sequence, optimize_protein, reverse_translate

optimize_protein(protein, host="e-coli-k12", forbidden=["BsaI"], name="part A")
optimize_coding_sequence(dna, host="e-coli-k12", forbidden=["BsaI"], name="part A")
```

Both return a `CodingSequence`: `dna`, the `protein` it spells, the `host` it was written for,
the `forbidden` enzymes it is free of, and `changes` — one entry per codon moved, carrying the
old codon, the new codon, the amino acid both spell, and the site it took away. `recoded` is
false when a coded sequence was checked and left alone. `reverse_translate` is the plain
protein-to-DNA step on its own, with no site removal.

Read the docstrings rather than reconstructing a call:

```bash
pixi run python -c "from liulab_mbio.translate import optimize_protein; help(optimize_protein)"
```

## A coded sequence is checked, not written again

Hand DNA to `optimize_coding_sequence` and every codon that spells no forbidden site stays
exactly as it is. Use it whenever the user already has a coding sequence: re-coding one in
silence would throw away choices they may have reasons for.

## Every site, both strands

A site is found on either strand, because an enzyme reads it on either. A change that would
clear one site by spelling another forbidden one is passed over — so pass every enzyme the
scheme uses in one `forbidden` list, rather than one enzyme at a time.

## When it refuses

`SiteNotRemovableError`, a `ValueError`, when a forbidden site is left that no synonymous change can
reach. The message names the sequence, the site, its strand, the codons it covers and why: the
amino acids there have one codon each, or every synonym spells another forbidden site. Pass that
cause on — the choice is the user's, either changing an amino acid or dropping an enzyme from the
list. Never hand over a part carrying a site it was told to avoid.

A `ValueError` also comes back for a sequence that is not a coding sequence — a length that is
not whole codons, a base outside ACGT, or a stop before the end — and a `KeyError` for a host or
an enzyme name nothing ships.

## Before you hand it over

Tell the user the host, the enzymes the sequence is free of, and every codon that changed: old
codon, new codon, and the amino acid, which is the same either side by construction. Someone
paying a vendor for a gene should see what was done to it.

`golden-gate-assembly` plans a cloning reaction from sequences that already exist; come here to
write the sequence itself.

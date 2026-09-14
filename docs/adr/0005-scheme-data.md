---
search:
  exclude: true
---

# The scheme is data the user supplies, and the paper's is a worked example

The positions, the enzymes that cut internally, externally and bluntly, the stuffers, the cloning
scar and the barcode length are one scheme object, loaded from JSON the user supplies and validated
as it is loaded. The scheme the paper used ships as a worked example under `docs/`, with its values
sourced in `docs/research/protein-library-assembly.md`. The package ships no scheme as package
data.

This is surprising without the reason, because every other biology-layer table here is package
data, rebuilt by a script in `scripts/` and sourced in a note. The ligase profile is the precedent
for the other direction: the package ships none and reads a copy the user holds.

The trade-off is real and it is paid by the user. Supplying a file is friction, and a shipped
default would make the paper's library one command. Against that, the paper's scheme is one
instance of a general pattern — any number of positions, any four enzymes — and a shipped default
is the instance everyone would run, so the general path would stop being exercised. The stuffers
and barcodes are also the paper's own sequences, whose licence the package would then answer for.
That is hard to unwind: once they sit in `src/liulab_mbio/data/` they ship in every wheel, and
every signature above them has assumed one scheme.

## Considered options

- **Ship the paper's scheme as package data**, as the enzyme and codon tables ship. It takes on
  the licence question and makes one instance the default.
- **Hard-code it in `library/`.** The parameters become constants, so a second scheme is a rewrite
  and not a file.
- **Take the parts of the scheme as separate arguments.** They would cross every boundary one at a
  time, and the rules between them — a stuffer prefix ending in the next position's entry
  overhang, barcode plus scar a whole number of codons — would have nowhere to be checked.

## Consequences

The pipeline cannot run without a scheme file, and the worked example is what a user copies and
edits. An invalid scheme fails where it is loaded and not at the bench. An agent can read the
scheme and show the user what a design assumed. Being documentation rather than package data, the
example owes no builder in `scripts/`.

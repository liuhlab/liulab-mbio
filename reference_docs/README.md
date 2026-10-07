# reference_docs

Papers, manuals and web pages downloaded while researching a method. Everything here but
this file is git-ignored: vendors do not let us redistribute their documents, so what ships
is what we learned from them, cited in a note under `docs/research/`.

Downloads live here. What the lab writes is tracked under `docs/research/`, however much of
it was drafted while reading these files. The one thing of ours that stays is a script whose
whole job is to rebuild a dump of a file here, because the file it reads cannot move.

## Using it

Read anything here. Cite the source itself in the research note — its URL, version and page
— never a path here, which a fresh clone does not have.

## What belongs here, and what does not

A file here is reference material: downloaded, read, cited in a note under `docs/research/`,
and replaceable by that citation. A file a **shipped artefact is built from** is an input, not
a reference document. It is tracked — converted to a text format, in `tests/data/`, with its
provenance in the research note.

The test is whether deleting the file breaks a build. A paper does not; a sequence a shipped
record is rebuilt from does. `tests/data/dmx0001.gb` and `tests/data/pcr-blunt-ii-topo.gb` are
that second kind: `scripts/build_dmx_vector.py` rebuilds the iGGA destination from them, so CI
can check the destination anywhere rather than only on a laptop holding this tree.

## Adding to it

One directory per method or topic, named for it (`gibson/`, `gateway/`). Inside, sort into
subdirectories that say what they hold — by kind (`papers/`, `manuals/`, `web/`), by source
(`twist/`, `Lund2024/`), or however the material divides.

- Keep the source and one copy you can grep: a PDF with its `pdftotext` dump beside it, a
  web page as `.txt` or `.md`, never its raw HTML.
- Name a file for what it is: author and year for a paper, the vendor's filename for a
  manual, the snapshot date for an archived page.
- A method directory has a `README.md`: the research note it serves, what each subdirectory
  holds, and what was left out with where it is online. A source worth an account of its own
  gets a `README.md` in its own directory, listed from the method's.
- Per-file provenance — URL, version, retrieval date — belongs in the research note, which is
  tracked and outlives the download. A README here need not repeat it.

## Keeping it clean

Before you finish, delete what was only a step: raw HTML once its text is out, page renders,
duplicate or superseded downloads, failed fetches, `.DS_Store`. A source that stays online at
a fixed address — a patent, an open-access paper, a GenBank accession — need not be kept; put
its URL in the README instead.

# reference_docs

Papers, manuals and web pages downloaded while researching a method. Everything here but
this file is git-ignored: vendors do not let us redistribute their documents, so what ships
is what we learned from them, cited in a note under `docs/research/`.

## Using it

Read anything here. Cite the source itself in the research note — its URL, version and page
— never a path here, which a fresh clone does not have.

## Adding to it

One directory per method or topic, named for it (`gibson/`, `gateway/`). Inside, sort by kind:

| Subdirectory | What goes in it |
| --- | --- |
| `papers/` | Journal articles, preprints, patents |
| `manuals/` | Vendor manuals, protocols, product sheets, app notes, brochures; one subdirectory per vendor when there are several |
| `web/` | Web pages as text; `live/` and `archived/` when both were read |
| `sequences/` | GenBank, SnapGene or FASTA records |

- Keep the source and one copy you can grep: a PDF with its `pdftotext` dump beside it, a
  web page as `.txt` or `.md`, never its raw HTML.
- Name a file for what it is: author and year for a paper, the vendor's filename for a
  manual, the snapshot date for an archived page.
- Each method directory has a `README.md`. It names the research note it serves, lists every
  file with its source URL and version, and lists what was left out and where it is online.

## Keeping it clean

Before you finish, delete what was only a step: raw HTML once its text is out, page renders,
duplicate or superseded downloads, failed fetches, `.DS_Store`. A source that stays online at
a fixed address — a patent, an open-access paper, a GenBank accession — need not be kept; put
its URL in the README instead.

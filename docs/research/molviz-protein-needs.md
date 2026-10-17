---
search:
  exclude: true
---

# What liulab-protein will need a viewer to draw

Research note for issue #576, on map #572. It lists the views `liulab-protein`, the viewer's
second user, would ask `liulab-molviz` for. Needs only: no design, no package choice.

**What was read, on 2026-10-17:** `liulab-protein`'s default branch at commit `eab571e`
(committed 2026-09-04): its code, `README.md`, `CONTEXT.md`, its eight ADRs, its guides and
its research notes; and its whole tracker, 79 issues (all closed) and 16 pull requests (one
open, a dependabot bump), with their 180 comments. Read-only: nothing was run against it and
nothing was written to it.

Every claim below is **read** from a named file or thread, or marked **inferred**. "Today"
means 2026-10-17.

## 1. The answer first

| # | View | Starts from | What a reader does with it | Drawn in liulab-protein today? |
| --- | --- | --- | --- | --- |
| 1 | A structure in 3D, whole entry or one chain | `Structure.atoms` or `Chain.atoms`, a biotite `AtomArray` | turns and zooms the fold; looks at one chain alone | **Yes**: `Structure.view()`, `Chain.view()` over py3Dmol |
| 2 | A prediction coloured by its own confidence | a `Structure` from `ESMFold2.fold` or `Structure.from_file`; per-residue confidence in the B-factor column | sees which parts of the fold the model trusts | **Yes**, by a 3Dmol.js style passed to `view()` |
| 3 | The pairwise confidence matrix as a heatmap, chain boundaries marked | `Structure.confidence.pairwise()`, an `(L, L)` float32 array | judges whether an interface in a complex is real | **No**. Named as the one uncovered gap |
| 4 | Per-residue confidence as a line along the chain | the B-factor column | the same as view 2, flat | **No**, and declined: view 2 covers it |
| 5 | An alignment: its rows against the query's columns, and its depth by position | an `MSA` | sees whether the alignment is deep where it must be | **No**. Deferred, not declined |
| 6 | Sequence search hits | the `DataFrame` from `Protein.search` | ranks hits; sees where on the query each one lands | **No**. A table only |
| 7 | Structure search hits | the `DataFrame` from `Structure.search` or `Chain.search` | the same, plus how well the shapes agree | **No**. A table only |
| 8 | Which stretches of a protein have a solved structure | the `DataFrame` from `Protein.structures` | picks the entry and chain to open | **No**. A table only |
| 9 | Sparse features along a sequence | a `SaeActivation`, joined to the feature descriptions | sees which named features fire where | **No**. A lookup only |

Views 1 and 2 exist, so for them the need is replacement, if `liulab-protein` ever moves off
py3Dmol. View 3 is the only one its own tracker names as missing. Views 4 to 9 are what its
types hold and no thread asks to see, except view 5, deferred, and view 4, declined.

**No protein sequence with domains.** A `Protein` carries no feature or domain annotation at
all: its sequence, an accession, an entry name, a description and the FASTA header fields
`organism`, `taxon_id`, `gene`, `protein_existence` and `sequence_version` (read,
`src/protein/core.py` and `src/protein/db/swissprot.py`, 2026-10-17). Nothing in the package
or its tracker fetches UniProt features or InterPro domains. A sequence view drawn for
`liulab-protein` today would carry the per-residue spans its own types hold instead: the SIFTS
ranges (view 8), search-hit ranges (views 6 and 7), per-residue confidence (view 2) and SAE
features (view 9).

## 2. Each view

### 2.1 A structure in 3D

**Starts from.** `Structure` is one set of coordinates, a PDB entry's asymmetric unit or any
coordinate file; `Chain` is one polymer in it, reached as `structure["A"]`. `.atoms` is every
atom of the first model, ligands, water and nucleic acids included; `.models` holds every
model of an NMR entry (38 for `1L2Y`, per its docstring); `.chain_ids` lists the labels in file
order; `Chain.kind` is `"protein"`, `"nucleic"` or `"other"` (read, `CONTEXT.md`,
`src/protein/structure/structure.py`, `src/protein/structure/chain.py`, 2026-10-17).

**What is drawn today.** `protein.structure.view.view(atoms, name=, width=, height=, style=)`
serialises the atoms as mmCIF text and returns py3Dmol's own `view`. The default style is one
cartoon per chain coloured N to C (`{"cartoon": {"color": "spectrum"}}`); `style={}` hands the
drawing to the caller. A chain's view writes that chain's atoms, so it shows that chain and
not the entry around it (read, `src/protein/structure/view.py`; merged in PR #96 on
2026-09-04).

**What the guides ask the viewer for** (read, `docs/guides/structures.md`, 2026-10-17):

- turn by dragging, zoom by scrolling;
- a cartoon coloured N to C, by chain;
- one chain alone, protein or nucleic (`1BNA` chain A, one DNA strand);
- sticks with a carbon colour scheme, and a van der Waals surface at 0.6 opacity;
- colour by a per-residue number read from the file (view 2).

**Where it is shown.** Three places (read, ADR-0008, `docs/guides/structures.md`,
`docs/start.md`, 2026-10-17):

- a notebook, by `.show()`;
- a static docs page, by `.write_html()` returning HTML that `markdown-exec` writes into the
  page while the site builds, with no network and no `LIULAB_DATA`;
- a page of its own, by `.write_html(open_file)`.

**Constraints a replacement inherits** (read, ADR-0008 and the module docstrings,
2026-10-17):

- **A view is asked for.** Neither class has `_repr_html_`, because displaying a
  `Structure` must not parse its file or reach RCSB.
- **Labels are case-sensitive and may be longer than one character.** `10EG` holds both
  `A` and `a`; mmCIF is the format because it carries such labels.
- **The file carries no secondary structure.** Only `atom_site` is written, and 3Dmol.js
  computes the helices and sheets itself.
- **The page holds the coordinates inline**, so a large entry makes a large page.
- **Today's page needs a network.** py3Dmol loads 3Dmol.js from a CDN. That falls short of
  the "one file opens offline" promise map #572 lists for mbio (decided 2026-10-17).
  **Inferred**: `liulab-protein` would gain offline pages from molviz, though it has not
  asked for them.

### 2.2 A prediction coloured by its confidence

**Starts from.** A `Structure` returned by `ESMFold2.fold`, or one reopened with
`Structure.from_file`. Per-residue confidence is the mmCIF's B-factor column, "where every
viewer looks" (read, `CONTEXT.md`, entry Prediction, 2026-10-17). The scalars are on
`Structure.confidence`, a `Confidence(plddt, ptm, iptm, pairwise_file)`, each scalar 0 to 1
(read, `src/protein/fold/predictions.py`, 2026-10-17).

**Two scales.** The B-factor column runs 0 to 100, while `Confidence.plddt` runs 0 to 1. The
docs once coloured with `min: 0, max: 1` and painted every residue one colour; PR #97 fixed
it (read, PR #97, merged 2026-09-04). Both guides now pass
`{"prop": "b", "gradient": "roygb", "min": 50, "max": 100}`, blue confident and red not
(read, `docs/start.md`, `docs/guides/structures.md`, 2026-10-17).

**A reopened prediction carries no `Confidence`.** On a cache hit, or after `from_file`,
`confidence` is `None` while the B-factor column still holds pLDDT. On a crystal structure
the same column is a temperature factor, and no value range tells the two apart (read,
issue #93 comment of 2026-09-04; `docs/guides/folding.md`, 2026-10-17). **Inferred**: a
viewer colouring by confidence cannot work out on its own whether the column means
confidence; the caller has to say so.

**What a prediction holds.** Protein, DNA and RNA chains fold together. The getting-started
page shows AP-1: FOS and JUN bZIP domains with a DNA duplex, 4 chains and 164 residues (read,
PR #97, 2026-09-04). Chain labels come from position in the request, so a caller reads them
off `chain_ids` (read, `docs/start.md`, 2026-10-17).

### 2.3 The pairwise matrix

**Starts from.** `Structure.confidence.pairwise()` reads `<name>.pairwise.npy` beside the
mmCIF and returns `(L, L)`, one entry per pair of residues; lower is surer. It raises
`FileNotFoundError` where the model reported no matrix. The file is written only for a
prediction that reported one, and the folding guide calls it the thing that "says whether two
chains sit correctly against each other" (read, `src/protein/fold/predictions.py`,
`docs/guides/folding.md`, 2026-10-17).

**Asked for.** Issue #93 proposed `Structure.plot_pairwise(ax=None)`, "the `(L, L)` matrix
as a heatmap", with chain boundaries marked (opened 2026-09-04). Its closing comment
(2026-09-04) dropped the confidence curve and kept this one: "The pairwise matrix has no
substitute … there is no way to look at it … if plotting comes back, it comes back for the
heatmap." It closed without filing a follow-up, pending somebody hitting the gap in practice.

**What a reader needs on it** (read, issue #93 body, 2026-09-04): the matrix, and where one
chain ends and the next begins. **Inferred**: the boundaries come from the prediction's
chains in `chain_ids` order and each chain's residue count; nothing in the code states that
the matrix's row order is that order.

### 2.4 Per-residue confidence as a line

Proposed as `Structure.plot_confidence(ax=None)` in issue #93 (2026-09-04), then dropped in
its closing comment the same day: the 3D view coloured by the B-factor column "reads better
than a curve". Listed here so the roadmap does not count it as an open need.

### 2.5 An alignment

**Starts from.** An `MSA`: `(header, row)` pairs of A3M text, row 0 the query.
`.query`, `.query_header`, `.depth` (rows), `.match_states` (columns), `.to_a3m()` and
`.compress()`. Case is the match state: an uppercase residue or `-` holds a column, and a
lowercase residue is an insertion that holds none. Headers are carried byte for byte,
because a `key=<taxon>` field pairs the chains of a complex, and so is a leading `#` line
(read, `CONTEXT.md` entry MSA, `src/protein/msa/msa.py`, `docs/guides/alignments.md`,
2026-10-17). Made by `Protein.msa(db)` or `align(sequences, query=...)`.

**Asked for.** Issue #93, under "Not in this issue" (2026-09-04): "An MSA coverage plot, rows
against position, is the standard way to see whether an alignment is deep where it needs to
be." It was deferred until the alignments guide existed. The guide shipped and points at no
picture (read, `docs/guides/alignments.md`, 2026-10-17).

**What a reader needs on it.** Depth at each query column, and **inferred** the rows
themselves laid against the query's columns, insertions hidden or marked. The guide's own
depth rule, a few thousand rows as a floor, says what a reader checks (read,
`docs/guides/alignments.md`, 2026-10-17).

### 2.6 Sequence search hits

**Starts from.** `Protein.search(db)` runs MMseqs2 and returns a `pandas.DataFrame`, one row
per hit: `query`, `target`, `pident` (a percentage), `alnlen`, `mismatch`, `gapopen`,
`qstart`, `qend`, `tstart`, `tend`, `evalue`, `bits`. No aligned strings: the package asks
the tool for these columns and no others (read, `src/protein/external.py`,
`docs/guides/search.md`, 2026-10-17).

**Inferred** need: hits as spans along the query, from `qstart` and `qend`, ranked by `evalue`
or `bits`. Nothing in the code, docs or tracker asks for a picture of a hit table.

### 2.7 Structure search hits

**Starts from.** `Structure.search(db)` or `Chain.search(db)` runs Foldseek. The table has the
same columns, with `fident` (a fraction, not a percentage) for identity, plus `alntmscore`
and `lddt`, which score how well the shapes agree. `query` is `<entry>_<chain>` (read,
`src/protein/external.py`, `docs/guides/search.md`, 2026-10-17).

**A target is not a `Structure` key.** In Foldseek's `pdb100` a target is a biological
assembly chain such as `1ubq-assembly1_A`, and a label may be longer than one character
(`6dwu-assembly44_DU`), while `Structure("1UBQ")` holds the deposited asymmetric unit (read,
`docs/research/foldseek-pdb100-contents.md` and `CONTEXT.md` entry Structure, 2026-10-17).
`pdb100` holds C-alpha traces only and renumbers residues, so coordinates for a hit come from
RCSB (read, `CONTEXT.md`, 2026-10-17). **Inferred**: a "query and hit side by side in 3D"
view needs the target named as an entry and a chain the asymmetric unit carries, which the
hit table does not give directly.

No thread asks for a picture of these hits.

### 2.8 Which stretches of a protein have a structure

**Starts from.** `Protein.structures`, a `DataFrame` with one row per SIFTS segment: `pdb`,
`chain`, `accession`, `res_beg`, `res_end` (the range in the chain) and `sp_beg`, `sp_end`
(the range in the UniProt sequence). One chain and accession may carry several rows, and the
two ranges are not always the same length, so both come back verbatim and no offset is
computed (read, `CONTEXT.md` entry SIFTS, `src/protein/sifts.py`, `docs/concepts.md`,
2026-10-17). `Chain.uniprot` is the other direction, a tuple of none, one or several
accessions.

**What a reader does.** The concepts page walks from an accession to a row to
`Structure(row["pdb"])[row["chain"]]` (read, `docs/concepts.md`, 2026-10-17). **Inferred**
need: the segments drawn along the UniProt sequence, so a reader picks the entry that covers
the residues they care about. Nothing asks for it.

**A trap for any view joining sequence to structure.** `Chain.sequence` is the observed
sequence: a residue with no coordinates is absent, so a disordered loop is missing rather
than filled from `SEQRES` (read, `src/protein/structure/chain.py`, 2026-10-17). A chain's
residue index, its author numbering and the UniProt position are three different numbers.

### 2.9 Sparse features along a sequence

**Starts from.** A `SaeActivation`: `indices` and `values`, both `(L, k)`, naming which
codebook features fired at each residue and how hard; `reconstruction_loss`, `(L,)`; `.max()`
for the per-sequence vector; and the facts naming the SAE, its parent, layer and codebook.
`features.descriptions("6b-layer60-k64-cb16384")` returns a `label` and a `description` per
feature index, for that one SAE only (read, `src/protein/embed/esm/sae.py`,
`src/protein/embed/esm/features.py`, `docs/guides/embedding.md`, 2026-10-17).

**What a reader does today.** Joins by index: `frame.loc[activation.indices[0]]`, "what fired
at the first residue" (read, `docs/guides/embedding.md`, 2026-10-17). **Inferred** need: a
track along the sequence showing where a chosen feature fires, labelled from the
descriptions. This is the nearest thing `liulab-protein` holds to a domain track. Nothing
asks for it.

**Embeddings themselves** (`Embedding`, `(L, d_model)`) have no view asked for anywhere, and
none is implied beyond what the SAE gives.

## 3. Where liulab-protein's own stance on drawing came from

Read from the tracker, 2026-10-17, in order:

1. **2026-09-04, issue #92.** The docs-site plan noted "Nothing here plots anything" and
   proposed matplotlib in the docs environment only, pointing readers at ChimeraX or PyMOL
   for the structure: "A 3D cartoon is not something this package should grow."
2. **2026-09-04, PR #96.** That was overtaken the same day: `view()` over py3Dmol landed,
   with ADR-0008. py3Dmol is held rather than wrapped, so a caller gets every 3Dmol.js call.
   It is the one runtime dependency taken from PyPI rather than conda-forge, because the
   conda-forge recipe makes IPython a hard dependency.
3. **2026-09-04, PR #97.** The getting-started page ends the fold with the AP-1 prediction in
   3D, coloured by confidence.
4. **2026-09-04, issue #93.** The flat plots were closed as not worth their cost, except the
   pairwise heatmap, left uncovered.

Nothing after that touches drawing. The PRs merged after #97 are a menu change (#98), a
fold-request change (#99) and its docs (#100), all on 2026-09-04, the date of the newest
commit on the default branch. The one open PR is a dependency bump (read, 2026-10-17).

## Open gaps

- **Whether liulab-protein moves to molviz at all.** ADR-0008 chose py3Dmol on purpose, for
  the whole 3Dmol.js API at the cost of one small module. Replacing it is that repo's
  decision. Nothing here asks it.
- **Aligned strings for a hit.** The hit tables ask for no aligned-sequence columns. Whether
  MMseqs2's and Foldseek's `--format-output` can add them, and whether `search`'s `extra`
  pass-through can change the format, was not checked. A hit drawn as an alignment depends on
  both.
- **Superposing two structures.** No code superposes a prediction on its crystal, or a query
  on a Foldseek hit. The getting-started page's "Where next" says to "open the prediction and
  the crystal, and look at them in 3D", but the structures guide shows them one at a time.
  Whether a superposed view is a need is open, and whether Foldseek reports the transform was
  not checked.
- **Assembly against asymmetric unit.** A Foldseek target names an assembly chain, and a
  `Structure` holds the asymmetric unit (section 2.7). How a hit is opened in 3D is
  unresolved in `liulab-protein`.
- **The pairwise matrix's order.** That its rows follow `chain_ids` and each chain's residues
  is inferred, not read (section 2.3).
- **Domain annotations.** No source exists in `liulab-protein`. A domain track would need
  data that package does not fetch today.
- **The published site was not opened.** The two viewers on it were read from the page
  sources and the PRs, not looked at.

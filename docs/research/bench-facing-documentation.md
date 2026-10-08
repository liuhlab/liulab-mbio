---
search:
  exclude: true
---

# How comparable tools document themselves for a wet-lab reader

Research note for issue #489, a sub-issue of map #487. It surveys five comparable tools and asks
one question of each: how fast does an experimental biologist — strong wet-lab knowledge, basic
programming — learn what the tool ships and what it does for them at the bench? Findings feed the
page list (#490) and the model page (#488).

Everything below was read on **2026-10-08** from each tool's own documentation. No blog posts or
third-party tutorials were used as evidence.

**One access note.** `support.snapgene.com` returns HTTP 403 to a scripted fetch of its HTML. Its
section tree and article bodies came from the Zendesk Help Center JSON API on the same host, which
answers plainly: `/api/v2/help_center/en-us/categories.json`, `/sections.json`, and
`/api/v2/help_center/articles/search.json`. Those are the same articles the site serves.

## 1. The verdict, first

| Tool | Reader it actually writes for | Time to "what can this do for me" | Takes | Leaves |
| --- | --- | --- | --- | --- |
| **pydna** | a biologist who already decided to script | Fast — the top nav is four items and one is a cloning notebook list | A section tree named after bench operations: restriction and ligation, PCR, primer design, Gibson, CRISPR | Toy sequences; nothing says what comes out at the bench; no protocol, no order sheet |
| **Biopython** | a bioinformatician | Slow — 26 chapters, and chapter 2 is the only orientation | Chapter 2 is the best single page in the survey: one real dataset carried the whole way | A flat chapter list sorted by module, not by task; cloning is absent |
| **Benchling SDK** | a developer integrating a platform | Never, for our reader — the docs answer "how do I call the API", not "what does this do" | Honest scoping: it says up front the SDK adds nothing the API does not already have | Biology appears only as sample data; every heading is a software noun |
| **SnapGene** | exactly our reader | Immediate — the user guide's sections are named after methods and menus a biologist already knows | A section per cloning method, each holding three to eight short task articles; six short videos as the front door | Nothing reusable in the writing itself — it is screenshot-led and tied to one GUI |
| **Addgene protocols** | the person running the bench step | Immediate — a flat list of ~50 procedures grouped by what you are trying to do | Page shape worth copying outright: Summary, Why this method, Procedure, Tips and Troubleshooting, Resources and References | No software at all, so nothing about where a number came from |

**The strongest single finding.** The two tools written for our reader — SnapGene and Addgene —
both organise by **what the reader is trying to do**, and both put a short orienting page in front
of the detail. The two written for programmers — Biopython and Benchling — both organise by **what
the software contains**. pydna sits between: its section tree is bench-named, but its pages are
written as library tutorials.

## 2. pydna — the closest analogue

Sources: [docs home](https://pydna.readthedocs.io/),
[`index.rst`](https://pydna.readthedocs.io/stable/_sources/index.rst),
[`getting_started.rst`](https://pydna.readthedocs.io/latest/_sources/getting_started.rst),
[example gallery](https://pydna.readthedocs.io/stable/example_gallery.html),
[how pydna works](https://pydna.readthedocs.io/stable/how_pydna_works.html),
[README](https://github.com/BjornFJohansson/pydna/blob/master/README.md).

### Section tree

Four top-level pages, plus an API reference:

```text
installation
getting_started
example_gallery
how_pydna_works
reference/pydna
```

`getting_started` is itself a grouped list of nine notebooks:

| Group | Notebooks |
| --- | --- |
| First steps | `Dseq`, `Dseq_Features`, `Importing_Seqs` |
| Cloning | `Restrict_Ligate_Cloning`, `PCR`, `primer_design`, `Gibson`, `CRISPR` |
| Cloning history | `history` |
| More examples | a pointer to the gallery |

That is the whole nav. A reader sees the entire surface on one screen.

### Page shape

Every teaching page is a Jupyter notebook rendered to HTML. `Restrict_Ligate_Cloning` runs:
title, then *Cutting with one or more restriction enzymes*, *Ligating fragments*, *Circularizing
fragments*, and a closing *Extra Notes: what happens to features when cutting and ligating?* Prose
and code are roughly balanced, about six code cells, each followed by its real printed output —
`Dseqrecord` summaries showing size, circularity, feature count and sticky or blunt ends.

### Examples

Inline, executed, with outputs shown. The first cell of every notebook is a `%%capture` block that
installs pydna when running on Colab, and each notebook carries a Colab badge, so a reader can run
it in a browser with nothing installed. That is a real strength.

The sequences are **toy**: `Restrict_Ligate_Cloning` works on a 338 bp circular `sample_seq.gb`
with one gene feature. The example gallery is where real material lives — Gibson assembly of
*R. cellulolyticum* fragments reproducing the original Gibson paper, CRISPR deletion in
*K. phaffi*, promoter libraries in *S. cerevisiae* — but the gallery is a link index with
one-line descriptions and no code on the page.

### Time to orientation

Fast. The README's opening sentence states the scope ("a human-readable formal description of
cloning and genetic assembly strategies", for simulation and verification) and then lists the
operations: primer design, PCR, restriction digestion, ligation, gel image generation, homologous
recombination, Gibson assembly, with Golden Gate marked in progress. The first README code example
builds a 60 bp record, labels a feature, and prints GenBank.

### Right for our reader

- **The section tree is named in bench words.** "Restriction and ligation", "PCR", "Gibson",
  "CRISPR" — a biologist scanning the nav knows immediately whether their job is covered.
- **Capability list on the landing page**, as a plain list of operations, not of modules.
- **Runnable in a browser.** The Colab badge removes the install from the first contact.
- **Outputs are shown.** The reader sees what a result looks like before deciding to install.

### Wrong for our reader

- **Toy sequences teach nothing about the bench.** A 338 bp `sample_seq.gb` is a software fixture.
  The biologist cannot map it onto their own plasmid, and the gallery that would — real genomes,
  a reproduced paper — is one click away and carries no code on the page.
- **The output stops at the screen.** Nothing in the teaching pages produces a file anyone takes
  to a bench: no primer sheet, no order list, no protocol. The reader learns the library can
  simulate a ligation, not what they then do on Monday.
- **No page says what is shipped as a whole.** There is no catalogue. `how_pydna_works` is one
  section long and covers only the assembly model.
- **The API reference is docstrings.** The project's own guidance is to read them with Python's
  `help()`. Our reader will not.

## 3. Biopython — the genre's reference point

Sources: [Tutorial contents](https://biopython.org/docs/latest/Tutorial/index.html),
[chapter 2, Quick Start](https://biopython.org/docs/latest/Tutorial/chapter_quick_start.html).

### Section tree

26 chapters, flat, in this order: Introduction; Quick Start; Sequence objects; Sequence annotation
objects; Sequence Input/Output; Sequence alignments; Pairwise sequence alignment; Multiple Sequence
Alignment objects; pairwise2; BLAST (new); BLAST and other search tools; Entrez; Swiss-Prot and
ExPASy; the PDB module; PopGen; Phylo; motifs; Cluster analysis; Graphics including GenomeDiagram;
KEGG; phenotype; Cookbook; the testing framework; contributing; a Python appendix; bibliography.

The ordering axis is the package's own module layout. A reader with a task in hand has to already
know which module owns it.

### Page shape

A chapter is a long linear essay with numbered subsections, prose alternating with interactive
sessions. Chapter 2 runs: *General overview of what Biopython provides*; *Working with sequences*;
*A usage example*; *Parsing sequence file formats* (with FASTA, GenBank, and a subsection titled
"I love parsing — please don't stop talking about it!"); *Connecting with biological databases*;
*What to do next*.

### Examples

Interactive `>>>` sessions with `...` continuation and the printed result below. The first example
is three lines — `from Bio.Seq import Seq`, then `Seq("AGTACACTGGT")` — self-contained and runnable
in a bare Python session.

The data is **real and carried through**: `ls_orchid.fasta` and `ls_orchid.gbk`, a genuine NCBI
nucleotide search for Cypripedioideae returning 94 records, with real accessions (Z78533.1,
Z78439.1) quoted in the output. The orchid dataset runs through the entire chapter, so every later
example builds on something the reader has already seen.

### Time to orientation

Slow by structure, fast if the reader lands on chapter 2. Chapter 2 is doing all the work: it
opens with a prose overview of what the package provides, then gets the reader to a running
example within a screen. Chapters 3 to 26 are reference material that happens to be written as
prose.

### Right for our reader

- **One real dataset, carried through a whole chapter.** This is the single most transferable idea
  in the survey. The reader builds a mental model of one record and reuses it; no example costs
  them a fresh setup.
- **A deliberate orientation chapter** that exists only to answer "what can you do with this".
  It is even titled with the question.
- **Output shown under every snippet**, so the page can be read without running anything.

### Wrong for our reader

- **The chapter list is the module list.** A biologist wanting to clone has no entry point —
  cloning is not a chapter, because it is not a module.
- **26 chapters is not a surface anyone scans.** There is no page saying what the package ships.
- **Chapter length.** A chapter is an essay; the reader cannot tell where their answer sits.
- **No bench artefact anywhere.** It is a library tutorial throughout, which is correct for
  Biopython and wrong for us.

## 4. Benchling developer docs — a commercial tool, the same reader

Sources: [developer platform access](https://help.benchling.com/hc/en-us/articles/9714802977805),
[docs index](https://docs.benchling.com/docs),
[Getting Started with the SDK](https://docs.benchling.com/docs/getting-started-with-the-sdk),
[Common SDK Interactions and Examples](https://docs.benchling.com/docs/common-sdk-interactions-and-examples).

### Section tree

Titled "Warehouse & API Guides". Sections: General (Developer Platform Overview, Limits,
Stability); Benchling Apps (nine pages); Webhooks; App Canvas; Python SDK (one page: Common SDK
Interactions and Examples); Warehouse (an overview, then fifteen warehouse-table pages including
Molecular Biology, Inventory, Workflows, Registry); Events; Analyses; Documents; Technical
Accelerators; Examples (four, including "Interacting with the API using R"); V3 APIs.

Every heading names a part of the platform.

### Page shape

*Getting Started with the SDK* runs: About; Why Use the SDK; Getting Started (Installation —
stable and preview; Using the SDK — API-key auth, OAuth auth, then one sample call listing DNA
sequences); Next Steps; Support. Callouts carry caveats — one notes the models come from
Benchling's OpenAPI spec, another states plainly that the SDK adds no functionality beyond the
underlying APIs.

*Common SDK Interactions and Examples* is a long single page of software concerns: entity
creation, updates, registration, pagination, schema fields, async tasks, retries, error handling,
the `returning` parameter, unsupported endpoints, client customisation, self-signed certificates,
the `Unset` type, and forward compatibility in enums and polymorphic types.

### Examples

Inline Python with placeholder values — `https://my.benchling.com`, an API-key stand-in. Biology
appears only as sample payload: a DNA sequence entity built from a base string with `species` set
to `"mouse"`. The task framing is entirely software: client, endpoint, serialisation, timeout,
enum.

### Time to orientation

Our reader never gets there, and before that they cannot even open the door: developer-platform
access is gated behind a tenant admin setting it to Full in the Tenant Admin Console. The docs
assume the reader has decided to integrate and now needs call mechanics.

### Right for our reader

- **Scope stated up front, including what the thing is not.** "The SDK adds no functionality
  beyond the underlying APIs" saves a reader an hour. We should be as blunt about what the package
  does not do.
- **Callouts carry the caveat, not the body text.** The prose stays an instruction; the warning
  sits beside it. That maps cleanly onto our protocol `cautions` discipline.
- **Installation split stable and preview**, so nobody installs an alpha by accident.

### Wrong for our reader

- **Every heading is a software noun.** "Pagination", "Forward-Compatibility in Enums". A
  biologist scanning this tree learns nothing about what they could accomplish.
- **No task is ever posed in biology terms.** Nothing is titled "register a plasmid" or "pull a
  sequence into a notebook entry"; the Molecular Biology page is a warehouse-table schema.
- **Access gating comes before orientation.** The reader must be entitled before they can read
  what they would be entitled to.

## 5. SnapGene help — written for exactly our audience

Sources: [old guide redirect](https://help.snapgene.com),
[help-centre categories](https://support.snapgene.com/api/v2/help_center/en-us/categories.json),
[sections](https://support.snapgene.com/api/v2/help_center/en-us/sections.json),
[article search](https://support.snapgene.com/api/v2/help_center/articles/search.json),
[Get Started series](https://snapgene.com/series/getting-started).

### Section tree

Three categories: **SnapGene User Guide**, **SnapGene FAQ**, **SnapGene Server Guide**.

The User Guide holds 46 sections, in this order: What's New; Installation and License Management;
Preferences; Projects; Display Options for DNA Sequences; Display Options for Protein Sequences;
Other Display Options; Info Panel; Enzymes; Features; Custom Feature Types; Primers; Primer
Database; Translations; Colors; Searching; Focus on Region; History; RNA; Proteins; CRISPR;
**Restriction Cloning and Linear Ligation; Gateway Cloning; Gibson Assembly; Golden Gate Assembly;
NEBuilder HiFi DNA Assembly; In-Fusion Cloning; TA or GC Cloning; TOPO Cloning**; PCR and
Mutagenesis; Agarose Gel Simulation; Collections; Batch Operations; Importing and Exporting;
Editing Sequences; Assembly and Alignment; Pairwise and Multiple Alignment; Single Stranded
Sequences; Nucleic Acid Secondary Structure; Designing Sequences; Codon Usage; Dotmatics
Bioregister; Electronic Laboratory Notebooks; Command Line Interface; Blast; Localization;
Shortcuts.

The FAQ is a parallel 15-section tree — General, Licenses, SSO, Printing and Exporting, Importing
Files, Files, Primers, Restriction Enzymes, Features, Translations, Agarose Gels, Alignment and
Assembly, Sequence Traces, History, Troubleshooting, Actions.

**Three things to notice.** First, each cloning method is its own section — eight of them, by
name, in the order a reader would meet them. Second, display and preferences come *before* the
methods: orientation first. Third, the FAQ is not a dumping ground; it mirrors the guide's topics,
so "how do I" and "why did it" live in parallel trees.

### Page shape

A section holds three to eight short task articles. Golden Gate Assembly holds, among others: *Set
the Enzymes for Golden Gate Assembly*; *Golden Gate Cloning of a Single Fragment into a Vector*;
*Golden Gate Cloning of Multiple Fragments into a Vector with Type IIS Sites*; *…with no Type IIS
sites*; *Golden Gate Cloning into a Vector Cut with non-Type IIS Enzymes*; *How is Golden Gate
Fidelity Predicted in SnapGene?*; *Golden Gate Assembly in Older Versions of SnapGene*.

Note the shape of that set: **one article per situation the reader is actually in**, not one
article per feature. "Does my vector have Type IIS sites?" is a question the biologist can answer
about their own plasmid, and it is the title.

Note also the one explanatory article, *How is Golden Gate Fidelity Predicted?* — a single page
answering "where does this number come from", separate from the task articles.

An article is short, imperative, and screenshot-led. *Set the Enzymes for Golden Gate Assembly*
states the default list (nine Type IIS enzymes, BsaI the default), gives the menu path, says what
each checkbox does, and says what the reader will then see in the dropdown. It is four screens
long and answers one question.

### Examples

There are no code examples — it is a GUI. The equivalent is the screenshot plus the menu path. The
material is real throughout: real enzyme names, real menus, real defaults.

### Time to orientation

Immediate, by two routes. Six videos, 1m30 to 3m each: *A Brief Tour of the Interface*; *Create,
Open and Import a File*; *Using the Project Panel*; *Annotating Sequences*; *Use Actions to Plan
and Predict*; *Use Tools to Analyze and Verify*. About fourteen minutes end to end, and the arc
goes orientation → daily file work → planning → verification. The written guide then covers
everything in detail.

### Right for our reader

- **A section per method, named after the method.** This is the clearest answer in the survey to
  "what does this ship": the reader counts eight cloning methods in the nav. We ship four plus
  two named methods, and should be as countable.
- **One article per situation, not per feature.** Titles phrased as the reader's own case.
- **A separate "where does this number come from" page per method**, so the task article stays
  imperative and the explanation still exists.
- **Orientation before methods**, and a short video arc as the front door.
- **A FAQ mirroring the guide's topics**, so troubleshooting has a predictable address.

### Wrong for our reader

- **46 sections is a lot of nav** — it works only because the method names are self-evident.
- **Screenshot-led writing does not survive.** Every article is tied to one version of one GUI.
  Nothing in it transfers to a package whose surface is a CLI and a protocol page.
- **It stops at the screen.** SnapGene tells the reader what to click, never what to pipette. The
  bench step is assumed. That gap is exactly what we fill, and the survey's two halves —
  SnapGene's method tree and Addgene's procedure page — do not exist together anywhere.

## 6. Addgene protocols — bench procedures for people who run them

Sources: [protocols index](https://www.addgene.org/protocols/),
[Gibson Assembly](https://www.addgene.org/protocols/gibson-assembly/),
[Plasmid Cloning by Restriction Enzyme Digest](https://www.addgene.org/protocols/subcloning/).

### Section tree

Five groups, flat inside each, roughly fifty protocols:

| Group | Count | What it holds |
| --- | --- | --- |
| Intro to the Lab Bench | 7 | PPE, biosafety, water baths, pipetting, centrifugation, microscopy, weighing |
| Basic Molecular Biology | 18 | plates, streaking, cultures, glycerol stocks, DNA purification and quantification, restriction digests, gels, ligation, transformation, PCR, gel purification, *How to Design a Primer*, *Sequence Analysis* |
| Plasmid Cloning | 6 | Restriction Cloning, Cloning by PCR, Annealed Oligo Cloning, Gibson Assembly, Ligation Independent Cloning, pLKO.1 |
| Virus | 11 | transfection, lentivirus, AAV production, purification, titration |
| Antibodies | 6 | transfection, purification, stains, western, ICC, ELISA |

The grouping axis is **how far along the reader is**, not what Addgene's software contains. The
first group assumes nothing. The methods group assumes the first.

### Page shape

Consistent and short. The Gibson page runs:

1. Title
2. **You may also like…** — three cross-links, above the content
3. **Summary** — about 85 words
4. **Why Gibson Cloning?** — four bullets, about 35 words
5. **Procedure** — six numbered steps, about 400 words, the bulk of the page
6. **Tips and Troubleshooting** — about 100 words, with named sub-cases ("Stitching fragments
   together using oligos"; "Number of fragments assembled simultaneously")
7. **Resources and References** — about 150 words: Gibson 2009, Gibson 2010, Rabe and Cepko 2020,
   plus an NEB resource

Three diagrams, four "Pro-Tip" callouts, no tables, no video. **The whole page is about 800
words.**

The subcloning page is the same skeleton with a *Background* and a *Design (Choosing enzymes)*
section added, and its procedure split into five named stages rather than numbered steps: Digest
your DNA; Isolate your insert and vector by gel purification; Ligate your insert into your vector;
Transformation; Isolate the Finished Plasmid.

### Examples

There is no code and no worked numeric example. The equivalent is the **expected result stated
inline**: the transformation step says the insert plate should carry significantly more colonies
than the vector-only plate; the diagnostic digest step says expect two bands, one vector-sized and
one insert-sized; gel purification asks for crisp bands. The reader is told what success looks
like at the step where they would doubt it.

Instructions are plain imperatives that name the number: "Incubate the mix for 1 hour at 50 °C or
follow manufacturer's instructions." Advice that is not an instruction goes to a Pro-Tip: "Avoid
strong secondary structures in the homology region."

### Time to orientation

Immediate. Fifty titles on one page, each one naming a thing the reader has either done or needs
to do.

### Right for our reader

- **The page skeleton is directly copyable**: Summary → Why this method → Procedure → Tips and
  Troubleshooting → Resources and References. We already render most of these fields in a protocol;
  the method pages should use the same spine so a reader meets one shape twice.
- **Why this method, in four bullets.** The reader's real first question is "should I use this",
  and it is answered before any procedure. No tool in the survey other than Addgene answers it.
- **Expected result stated at the step.** This is our `expected` field, and Addgene shows it
  belongs in the narrative page too, not only in the rendered protocol.
- **Pro-Tips keep the imperative clean.** Same separation Benchling gets with callouts.
- **References named and dated**, three or four per page — a bar our method pages can meet.
- **800 words is enough for a complete method.** Our method pages should be judged against that,
  not against a tutorial's length.
- **Cross-links above the content.** "You may also like" sits at the top, where a reader who
  landed on the wrong page can leave at once.

### Wrong for our reader

- **No design work shown.** The page says "design primers with 15–40 bp of homology" and leaves
  the reader to do it. That is the gap the package closes, and our pages must show the design
  being done, not just prescribed.
- **No provenance per number.** Numbers appear as settled fact with a reference list at the end,
  not a citation at the number. ADR 0019 already holds us to the stricter rule.
- **No materials list.** Reagents are named inside prose. Our protocol renders a materials table
  and an order sheet, which is better; the method page should link to it.

## 7. What recurs across all five

- **A landing page that lists operations, not modules.** pydna's README and SnapGene's video arc
  both do it; Biopython's chapter 2 does it in prose. Benchling does not, and is the hardest to
  orient in.
- **One short orienting unit in front of the detail.** Chapter 2; the six videos; the Summary
  block. In all three cases it is under fifteen minutes of the reader's time.
- **Real material beats toy material, and carried-through material beats both.** Biopython's
  orchids are the best example in the survey; pydna's `sample_seq.gb` the weakest.
- **Explanation is separated from instruction**, always by the same device: a Pro-Tip, a callout,
  or a separate "how is this predicted" article.
- **Every tool that works for our reader is organised by reader intent.** Every tool organised by
  software structure fails them. Map #487 already settled a reader-intent nav axis; this survey is
  evidence for it, not against.

## 8. What this suggests for us

For the page list (#490) and the model page (#488).

### Section tree

- **Keep the settled reader-intent axis** — Home / Methods / Protocols / Projects / Reference. The
  survey supports it: every bench-facing tool in it is organised that way.
- **Under Methods, one page per method, titled by the method.** SnapGene's eight cloning sections
  are the model. A reader should be able to count what we ship off the nav: Golden Gate, Gibson,
  restriction and ligation, Gateway, plus primer design, barcodes, codon optimisation and maps.
  Resist a "Cloning / Design / Maps" sub-layer unless the page list grows past what one screen
  holds — map #487 already flags that subsection split as unsettled, and this survey says collapse
  it.
- **One "what this ships" page on Home**, listing operations in bench words, as a list a reader
  scans in under a minute. pydna's four-item nav and SnapGene's video arc both achieve this; a
  97-module API dump achieves the opposite.
- **Protocols as a flat catalogue grouped by how far along the reader is**, after Addgene: the
  procedure comes first, the protocol system beneath it. Map #487 settled that order; Addgene's
  five groups show the grouping axis to use.
- **A FAQ or troubleshooting tree that mirrors the Methods tree**, after SnapGene, rather than one
  undifferentiated troubleshooting page.

### Page shape — the model page

Adopt Addgene's spine, with two additions the survey says are missing everywhere:

| Block | Length | Content |
| --- | --- | --- |
| Summary | ~80 words | what this method does, in one paragraph |
| When to use this | 4 bullets | the reader's real first question. Addgene's "Why Gibson Cloning?" |
| What you need | short | the inputs: which files, which sequences — our analogue of a materials list |
| Worked run | the bulk | one command or one function call, its real printed output, and the files it writes, named |
| What it produced | short | each output file, what it is for, with a link to the rendered protocol |
| Tips and troubleshooting | ~100 words | named sub-cases |
| Where the numbers come from | short | the research note or ADR behind the design choices — SnapGene's "How is Golden Gate Fidelity Predicted?", made a standing block |
| References | 3–5 | named and dated |

Target **under 1000 words** for a method page. Addgene covers a complete method in 800.

### Example policy

Map #487 settled a full generated, committed corpus per page. The survey sharpens three choices
inside that:

- **One real record, carried across pages.** Biopython's orchids are the pattern. We already have
  `docs/examples/pUC19-GFP/`; reuse it everywhere a method page needs an input, so a reader builds
  one mental model and never pays a setup cost twice. A method that cannot use it gets its own
  named record, not a toy.
- **Show the output, not just the call.** Every example on a method page carries its real printed
  result, and names the files written. pydna shows outputs and is better for it; its failing is
  that no output is a bench artefact.
- **Name the bench artefact in the example.** This is where we beat every tool in the survey. The
  worked run ends by naming `protocol.html`, the primer sheet and the product, and linking the
  rendered protocol. SnapGene stops at the screen; Addgene starts at the bench; nothing joins them.

### Two things to avoid

- **Do not organise anything by module, package or import name.** The split between `mbio` and
  `synbio` is stated once on the "two packages" page and shown by the CLI prefix — already settled
  on #487, and the survey is strong evidence for it. Benchling's tree is what happens otherwise.
- **Do not let the orientation page grow.** Its job is to be scanned in a minute. Every tool that
  kept it short is easy to enter.

## Open gaps

- **The SnapGene article bodies were read through the help-centre search API**, not the rendered
  pages, so screenshot density and page length were inferred from the text rather than measured on
  screen. SnapGene Viewer is installed on this machine and was not opened.
- **Benchling's help centre (`help.benchling.com`)**, which is the non-developer side and may be
  closer to our reader than `docs.benchling.com`, was read only through the one access article.
- **Article counts per SnapGene section were not tallied**; "three to eight" comes from the Golden
  Gate search result and the section list, not a full enumeration.
- **No reader was timed.** "Time to orientation" is read off structure — nav width, position of
  the first runnable example — not measured on anyone.

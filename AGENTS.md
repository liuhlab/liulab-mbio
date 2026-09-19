# liulab-mbio

Molecular biology design tools for DNA sequences, enzymes, primers and cloning. Two pipelines
plan an experiment end to end: a cloning job from a vector and its inserts, by any method under
`cloning/`, and a barcoded combinatorial library built from lists of proteins in rounds. Each
picks its enzymes, designs the DNA, simulates the product, and writes a bench protocol someone
can follow. Repo-local skills call them; the lab uses both. Distribution name
**`liulab-mbio`**, import name **`liulab_mbio`**.

**Easiest thing to get wrong: coordinates.** Every module is 0-based and half-open, and a span
across the origin of a circular record ends past the record's length. Read
`docs/adr/0001-coordinates.md` before touching a span. File formats convert at their own
boundary and nowhere else.

## Architecture

One direction, bottom to top — nothing lower imports anything higher.

| Layer | Modules | What it owns |
| --- | --- | --- |
| model | `sequence`, `checks` | `SequenceRecord`, `Feature`, `Segment`, `Primer`, the coordinate rule; `Check`, its status and the worst-of rule |
| files | `io`, `snapgene`, `edits` | reading and writing records; editing spans, and carrying what a record annotates into another one |
| biology | `enzymes`, `sites`, `codons`, `translate`, `barcodes`, `overhangs`, `ligase` | shipped enzyme data, cut sites, domestication; reverse translation and whole-sequence codon choice; distance-separated barcode sets; whether two cut ends anneal, the rules an overhang set is held to, its ligation fidelity, and a ligase profile the user holds |
| plot | `plot/` | a record drawn as a map: `drawing` is the way in, `layers` resolves items, `circular`, `linear` and `sequence_view` lay them out, `labels` keeps labels apart, `fonts` measures, `svg` and `page` write, `convert` makes a PNG or PDF |
| primers | `primers/` | `polymerase`: Tm, Ta and its PCR profile; `thresholds` and their wording; `placement`, `evaluation`, `design`; `genome`, which runs `ipcr` |
| protocol | `protocol/` | `model`, read from and written to JSON, and `render`, its self-contained HTML page |
| bench | `bench/` | what any pipeline shares: `amounts`, `reactions`, `pcr`, `gels`, `validation`, `inactivation`, `phenotype`, `oligos`, `steps` |
| pipeline | `cloning/` | `plan`, what every cloning plan writes and how it is judged; `goldengate/`: `design`, `assembly`, `bench` (its reaction and cycling), `oligos`, `steps`, joined by its own `plan`; `gibson/`: the same modules, where `design` chooses each junction's overlap and lays out a stitched part's and a bridging oligo, and `bench` holds each assembly product's own numbers; `restriction/`: those modules again, plus `digest`, `amplify`, `ligation` and `verdicts`, where `design` chooses the enzyme pair; `gateway/`: `att`, the site sequences and the arithmetic a junction follows, then `design` for the attB tail and its PCR, `recombination` for one reaction on two records, `checks`, `oligos`, `bench` and `steps` |
| pipeline | `library/` | `scheme`, `standard`, `parts`, `vector`, `rounds`, `coverage`, `bench`, `steps`, joined by `plan` |
| command line | `cli`, and each feature's own `cli` | the verbs: the root app mounts one sub-app per feature, `cloning/cli` one per method and the spine they share |

Each pipeline has one way in. `cloning.goldengate.plan_assembly` writes four files: the
product, the primer sheet, `protocol.json` and the `protocol.html` rendered from it.
`cloning.gibson.plan_gibson` and `cloning.restriction.plan_restriction` write the same four.
`cloning.gateway.plan_gateway` writes those four, and the entry clone as a fifth where it
planned a BP reaction.
`library.plan_library` writes the synthesis order sheet, the barcode and amino-acid change
tables, a record per round, the product, and those same two protocol files.
A pipeline's protocol is data an agent may edit and render again, never a place to invent a
number the package computes: `build-protocol` says how, `docs/adr/0002-editable-protocols.md` why.
A subpackage re-exports its own way in, for callers outside it. Inside the package, import a
name from the module that owns it; the top-level `__init__.py` re-exports only `__version__`,
and `Check`, `Junction`, `Part` and `Files` each mean different things in every module that
defines one — a `Junction` and a `Files` belong to the cloning method that defines them.

Package data is in `src/liulab_mbio/data/`. Each file is rebuilt by a script in `scripts/` and
sourced in a note under `docs/research/`. Never hand-edit one, and ship nothing whose licence
forbids it — a ligase matrix is read from a copy the user holds.

## Restraint

Generalizable, lightweight, uncustomized — in that order, and ahead of thorough.

Every gate, lint rule and cap is paid by everyone who works here afterwards. They arrive one at
a time, each reasonable alone, and nothing measures the sum. Before adding one: does it
generalize, what is the total someone must already satisfy, and would a narrower rule do?

**A measurement outranks a hypothesis.** When a rule is shown to fire on correct work, that is
evidence; a defect it might also catch is not. Removing one that misfires is a contribution.

## Toolchain

- **pixi** only. Never bare pip, uv, or conda. `pyproject.toml` holds the dependencies,
  environments and tasks.
- **Python 3.13**, declared by `requires-python` and the pixi pin. Write it anywhere else and
  `conformance` holds that copy to the floor.
- **hatchling + hatch-vcs**. The version comes from the newest git tag, CalVer `vYYYY.M.PATCH`.
  Never hand-edit a version.
- Platforms: `osx-arm64` and `linux-64`. Dependencies: `typer`, `biopython`, `primer3-py`,
  `vl-convert-python`, `pypdf`, and `ipcr` from bioconda, a pixi dependency only.

## Commands

| What | Command |
| --- | --- |
| the gate | `pixi run check` |
| one test | `pixi run test -- tests/test_sites.py::test_name` |
| the docs | `pixi install -e docs`, then `pixi run docs-build` |
| a Golden Gate plan | `pixi run liulab_mbio cloning goldengate plan VECTOR INSERT --out DIR` |
| a Gibson plan | `pixi run liulab_mbio cloning gibson plan VECTOR INSERT --out DIR` |
| a restriction and ligation plan | `pixi run liulab_mbio cloning restriction plan VECTOR INSERT --out DIR` |
| a Gateway plan | `pixi run liulab_mbio cloning gateway plan CARRIER DESTINATION --out DIR` |
| a library plan | `pixi run liulab_mbio library plan PARTS --scheme S --vector V --out DIR` |
| the skills | `python skills/install.py --target all`, and `--check` |

## Layout

```text
src/liulab_mbio/  the package
tests/            pytest, mirroring src/
docs/             the published site; docs/adr/, docs/agents/ and docs/research/ are agent-facing
skills/           repo-local agent skills
scripts/          the gate runner, and the package-data builders
reference_docs/   reference material downloaded for a method — papers, manuals; git-ignored
CONTEXT.md        the glossary — the words this repo uses
```

## Gates

`pixi run check` must be green before you commit: `lint` and `fmt-check` (ruff), `typecheck`
(pyright, `standard`, plus annotated parameters outside `tests/`), `vale` and `markdownlint`,
`conformance`, and `test`. It reports **all** failures, not just the first — read to the bottom
before fixing anything. The docs build is not part of it and runs as its own CI job.

Work on a branch and merge through a pull request.

Four traps:

- **`--doctest-modules` runs every `Examples` block in `src/` as a test.** An example that
  downloads, shells out or needs a large file wants `# doctest: +SKIP`, and the marker covers
  only the line it sits on.
- **`filterwarnings = ["error"]`.** Each tolerated warning gets a targeted entry in
  `pyproject.toml` with a comment saying why. Never a blanket ignore.
- **`src/` is collected**, so every module is imported at test time. Keep a heavy import inside
  the method body that needs it, not at module scope.
- **A substring of `--help` output is not a substring of the help.** Typer prints it through
  rich, which styles the first dash of an option on its own. Colour is off under `ssh` and on in
  CI; reproduce CI with `FORCE_COLOR=1`, and strip the styling before you assert.

## Writing rules

Three rules, all enforced by `vale`: be concise; agent-facing documents have word caps;
human-facing prose avoids jargon and stays readable. Read `docs/agents/writing.md` before
writing either kind — the caps are lower than you expect, and this file is subject to one.

### Comments and docstrings are short

Nothing checks these, so the rule is on you. A docstring says what a thing does and what it
promises. A comment says why a line is surprising. Neither is a notebook.

**Do not record** measurements (sizes, timings, percentages — true once, on one machine), the
environment (hostnames, paths, versions), implementation detail a reader can see, or history
(what a previous version did, what was tried, which ticket decided it).

Where a fact has to survive, it has a home that is checked:

| The fact | Where it lives |
| --- | --- |
| What the code does | a test pins it |
| A decision and its trade-off | an ADR |
| Vocabulary | `CONTEXT.md` |
| What is still open | an issue |
| A measurement | nowhere — delete it |

## Read next

| When | Read |
| --- | --- |
| Before changing code | `CONTEXT.md`, then any ADR covering the area |
| Writing or moving a test | `docs/agents/testing.md` |
| Recording vocabulary or a decision | `docs/agents/domain.md` |
| Filing or working an issue | `docs/agents/issue-tracker.md` |
| Labelling someone else's issue | `docs/agents/triage-labels.md` |
| Writing anything | `docs/agents/writing.md` |

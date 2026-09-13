# liulab-mbio

Molecular biology design tools for DNA sequences, enzymes, primers and cloning. The package
plans a Golden Gate experiment from a vector and its inserts: it picks a Type IIS enzyme,
designs the overhangs and the primers, simulates the product, and writes a bench protocol
someone can follow. A repo-local skill calls it; the lab uses both. Distribution name
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
| files | `io`, `snapgene`, `edits` | reading and writing records; editing spans |
| biology | `enzymes`, `sites`, `codons` | shipped enzyme data, cut sites, domestication |
| primers | `primers/` | `polymerase`, its Tm and Ta; `thresholds` and their wording; `placement`, `evaluation`, `design` |
| protocol | `protocol/` | the protocol model, and its self-contained HTML render |
| pipeline | `goldengate/` | `design`, `assembly`, `bench`, `ligase`, `steps`, joined by `plan` |

`goldengate.plan_assembly` is the one way in, and `Plan.write` puts the product, the primer
sheet and the protocol in a directory. `__init__.py` re-exports only `__version__`: `Check`,
`Fragment` and `Junction` each mean different things in two modules, so import by module path.
`Check` is the judged check in `checks` and its displayed form in `protocol`.

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
- Platforms: `osx-arm64` and `linux-64`. Dependencies: `typer`, `biopython`, `primer3-py`.

## Commands

| What | Command |
| --- | --- |
| the gate | `pixi run check` |
| one test | `pixi run test -- tests/test_sites.py::test_name` |
| the docs | `pixi install -e docs`, then `pixi run docs-build` |
| the pipeline | `pixi run liulab_mbio goldengate plan VECTOR INSERT --out DIR` |
| the skills | `python skills/install.py --target all`, and `--check` |

## Layout

```text
src/liulab_mbio/  the package
tests/            pytest, mirroring src/
docs/             the published site; docs/adr/, docs/agents/ and docs/research/ are agent-facing
skills/           repo-local agent skills
scripts/          the gate runner, and the package-data builders
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
| Recording vocabulary or a decision | `docs/agents/domain.md` |
| Filing or working an issue | `docs/agents/issue-tracker.md` |
| Labelling someone else's issue | `docs/agents/triage-labels.md` |
| Writing anything | `docs/agents/writing.md` |

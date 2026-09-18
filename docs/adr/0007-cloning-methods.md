---
search:
  exclude: true
---

# A cloning method is a directory under a shared roof, and the overhang rules sit below them all

Golden Gate was the only method, so its package held everything it used and its skill held
everything an agent had to read. #74 adds Gibson, Gateway and restriction and ligation cloning.
Each is a mechanism of its own, and each wants part of what Golden Gate owned.

So one shape, in two places. In the package a method is a package under `liulab_mbio.cloning`,
and `cloning.plan` and `cloning.cli` own what every plan is: the files it writes, its status,
and the spine each plan verb runs. In the skills, `skills/molecular-cloning/SKILL.md` chooses
the method and a method's own `METHOD.md` carries the detail, read only once the choice is made.

What a method is *held to* rather than what it decides — the ligation matrix, the fidelity
score, the rejection rules and the distance two overhangs must keep — sits at the biology layer
in `overhangs`, with the ligase profile reader beside it, below every pipeline. Golden Gate
keeps its own chooser: which Type IIS enzyme is free, and which overhangs the junctions take.

## Considered options

- **Leaving the overhang rules inside Golden Gate**, where they were written. `library` needed
  them and imported `cloning.goldengate.design` — the only sideways pipeline-to-pipeline import
  in the package, which the layer table forbids. Restriction and ligation cloning would be the
  next to reach across. A rule every method is held to, owned by one method, makes every other
  method import that one.
- **One skill per method**, as `golden-gate-assembly` was. An agent then carries every method's
  description on every task, and the judgement a bench scientist actually makes — which method
  fits this job — is written nowhere. Two other skills already linked to that one by name, so a
  user arriving from them landed in one branch instead of at the choice.

## Consequences

- `liulab_mbio cloning --help` is the supported set, read off the tool rather than off prose.
  The group's membership is mounted from the root `cli.py`: the shared spine and the group
  cannot both sit in `cloning/cli.py` without an import cycle.
- A method file is never `SKILL.md` — conformance holds every one of those to
  `skills/<name>/SKILL.md` — so it is `METHOD.md`, and the single Vale section covering skills
  widened to every Markdown file under `skills/`. A second section naming `skills/` fails
  another rule, so widening was the only move, not a preference.
- `cloning.plan` is a layer and not a method, so the library pipeline may import it. The guard
  keeping a pipeline out of `liulab_mbio.cloning` exempts that module by name.

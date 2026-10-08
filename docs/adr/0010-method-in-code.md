---
search:
  exclude: true
---

# The method is code and a project is a file, which reverses ADR 0005

`liulab_synbio.igga.method` holds `IGGA`: one `Scheme` instance built in code, whose
invariants are checked when the module is imported. What one build chooses — its positions, its
parts, its vector, its host, its oligo length, its batch size, its completeness, its seed, its
barcode length and distance, and any further enzyme it needs kept clear — is read from one
`project.json` and checked where it is read. `read_scheme()`, the scheme JSON format and the
CLI's `--scheme` are gone.

The test is what the shelf already carries. The four ccdB cassettes, the DMX vector, the part
carrier and the barcode kit are ordered molecules. A parameter their DNA spells is not a knob,
because turning it would not change the tube: `TATG`, `AGGA`, `TTCC` and `CTAA`, the four enzyme
roles, the 34 bp internal stuffer and the shared external stuffers are therefore constants.
Everything a second build could rationally choose is an input.

ADR 0005 argued the opposite, that shipping a scheme would make the general path stop being
exercised. Under this package's boundary synbio **is** the method, so there is no general path:
a second scheme is a second method, not a second file. The reason inverts, and 0005 is
superseded by this.

**Nothing overrides.** Every parameter sits on one side. Where a build legitimately varies
something the method also constrains, it composes rather than replaces: the enzymes a block is
kept clear of are the method's unioned with the build's `reserved_extra`. That is the guard
against the method becoming a spec language.

Two further consequences. `Position` collapses to a name, because every position carried
identical stuffers and a record of one field is not a record. And the barcode frame rule is
checked rather than stated: a build's barcode length plus the method's cloning scar is a whole
number of codons, which refuses 12 and admits 14.

## Considered options

- **Keep the scheme a file and ship the method's as a default.** The default is the instance
  everyone runs, so the file is friction with no reader.
- **Let a build override a method constant.** It turns the method into a configuration
  language, and no reviewer can then say what a build assumed.

## Consequences

An invalid method fails at import rather than at the bench, and an invalid build fails where it
is read. A designer reads one constant instead of threading a scheme through its signature. A
second build is a second `project.json` and no code change.

---
search:
  exclude: true
---

# A lab's habit goes in the skill that runs the command, not in a package default

A default is what every user of a verb gets without asking. **A package default is the value
whose numbers the package has sourced.** What this lab happens to keep on its shelf is a habit,
and it goes on the command line of the repo-local skill that runs the verb, as a flag the verb
already takes. CLAUDE.md's package boundary sorts `mbio` against `synbio`; this sorts the
package against `skills/`.

The worked case is the clean-up kit. The lab uses Zymo's `D4003`. The package default stays
NEB's `T1130`, because that is the kit whose manual `COLUMN_RECOVERY` is read from. A run naming
another kit gets no recovery sentence and no citation of that manual, rather than NEB's yield
under Zymo's name. Moving the default would strip that cited number from every restriction
page that cleans up a PCR, and hand every other user one lab's supplier. Each method's
`METHOD.md` instead runs its plan with `--cleanup-kit D4003`, so an agent following the lab's
skill names Zymo without anyone typing it.

The same question recurs for a host, a polymerase, a kit vendor and a supplier. Ask whether the
package has sourced a number to that value. If it has, that value may be the default; if the
only reason is that the lab uses it, the skill names it.

This is not ADR 0010's line. A method's own choice is code, and what one run chooses is its
build file. A habit is neither: the lab keeps it whichever method it runs.

## Considered options

- **Make `D4003` the package default.** A restriction page cleaning up a PCR loses its cited
  recovery, and every user outside the lab inherits the lab's supplier.
- **A default set in `synbio`.** Never read: the lab runs `mbio cloning goldengate plan …`,
  which imports nothing from `synbio`; the dependency runs the other way.
- **A config file the user holds.** ADR 0020's first store is for a fact a licence forbids
  shipping or one that goes stale. A preference is neither, and a config mechanism for one
  string is surface everyone pays for.

## Consequences

Nothing is built. Someone running the verb by hand, without the skill, gets `T1130` and, on a
restriction page, the recovery it cites. Nothing tests the skill's line; the verbs' own tests
pass `--cleanup-kit D4003`, so the flag and the key stay alive.

A habit with no flag to carry it is not expressed. The gel extraction kit is the restriction
method's own `T1120`, its manual being where the gel numbers come from, so the skill cannot name
Zymo's `D4007` there. A flag for it would follow this rule: `T1120` stays the default, and the
skill names `D4007`.

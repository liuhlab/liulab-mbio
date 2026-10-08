---
search:
  exclude: true
---

# Writing rules

Four rules. Two are checked by `vale` in `pixi run check`; the other two are on you.

## 1. Be concise

Shorter beats longer — in documents, issues, commit messages, replies, and code. If a
sentence survives deletion without loss, delete it. The caps below are ceilings, not targets.

Code's share of this rule is its comments and docstrings, which no cap reaches. `AGENTS.md`
says what belongs in one and where each fact lives when it does not.

## 2. Agent-facing documents have word caps

| Files | Cap |
| --- | --- |
| `AGENTS.md`, `docs/agents/*`, every `.md` under `skills/` | 1000 (`Lab.LengthDoc`) |
| `docs/adr/*` | 500 (`Lab.LengthAdr`) |
| `CONTEXT.md` | 200 words per glossary entry, checked by `conformance` |
| `docs/research/*` | none — research notes are long by nature |

`CONTEXT.md` has no file cap. A glossary grows one term at a time and never shrinks, so a
file-level count measures how big the domain is, not how well the entries are written — a cap
set for fourteen terms fires again at sixteen. The per-entry cap is the unit that matches, and
it is the only one. What can go wrong with a long glossary is the number of entries, which no
word cap was going to catch.

**Measure with the gate, not with `wc`.** `wc -w` counts table pipes, link targets and
shell flags as words; vale does not, and vale is what enforces the cap. Vale's `words` metric
counts a table's content as zero rather than as less, so a table row costs the cap nothing
while the same words in a paragraph cost it every one, and `wc` will tell you a reference
page is near a cap it is nowhere near. Run `pixi run vale` and believe it.

The caps are dials with tight defaults. Raising one is a one-line diff in `styles/Lab/` —
do that deliberately, and say why in the commit message. Do not raise a cap because a
document ran long; that is the cap working.

## 3. Human-facing prose avoids jargon and stays readable

`README.md` and every `docs/` page a human browses: no terms from the lab jargon list
(`Lab.Jargon` — architecture-speak and Latinate verbs), and reading grade 11 or below
(`Lab.Readability`). No length cap — a tutorial is as long as the task.

**When `Lab.Readability` fails, match a passing exemplar — do not attack the number.** The
grade is a ratio over structure, so shaving words inside sentences you already committed to
moves it by hundredths. Find text that passes — an older section, a sibling page — and rewrite
toward how it is segmented. Measured on a 21,000-word changelog: optimising the grade directly
moved it 12.46 to 12.41; matching a passing section in the same file moved it under 11, while
cutting only 3% of the words and *raising* the entry count.

**This rule also covers what you say, not only what you write.** Nothing checks a chat
reply, so it is on you: when a human asks, answer in plain language and explain the term
you would otherwise reach for. An unexplained term in conversation is the same failure as
one in the README.

## 4. Protocol text is written for the bench

Every field a protocol renders is read by someone standing at the bench with gloves on, who
cannot ask you anything. Nothing checks these, so the rule is on you.

| Field | What it holds |
| --- | --- |
| `instructions` | one imperative sentence an arm can execute: the action, the thing, the amount, the vessel, the time or temperature. No clause explaining the choice |
| `cautions` | what would hurt the person or the material, written as the action to take |
| `expected` | what the reader sees if it worked, written as an observation |
| `notes` | the only place a *why* may stand, and only with a citation. An uncited note is deleted |
| `figures` | one, on a step that changes a molecule: a PCR, a digest, a ligation, an assembly, a recombination. It draws what the step makes, with what the step changed lit. A step that moves, waits, incubates or reads carries none, and a plate draws itself. A figure that replaces no explanation is decoration |

Five things are forbidden in every rendered field: a name the reader cannot look up, such as a
route named by a letter, or a file the page does not link — `Protocol.files` is how it links
one; markup, which renders as the characters it spells and not as emphasis; a code identifier
or a module name; a tracker or code reference; and this package's own vocabulary where a bench
word exists.

The explanation you delete has a home already. A decision true for this run goes in the
protocol's own `overview` and `highlights`; a fact about the method goes to `docs/research/` or
the method's `METHOD.md`; a trade-off goes to an ADR.

## Which is which

Agent-facing means written for a machine that has to act: capped, exempt from the jargon
and readability rules, kept out of the site navigation. Human-facing means written for a
person reading the published site: checked for jargon and readability, uncapped. A file is
one or the other — if you are adding a document and cannot tell, it is human-facing.
Rendered protocol text is neither. It is product, `vale` reads none of it, and rule 4 governs
it.

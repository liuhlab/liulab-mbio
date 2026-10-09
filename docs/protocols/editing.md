# Change a protocol

The protocol is the JSON file; the page is made from it. So to skip a step, add one of your own,
or reword one, edit that file and make the page again. One rule holds the whole thing together:
never write in a number the package worked out. Change the wording and the order freely, and
leave the volumes, masses, counts and temperatures alone.

## Find the file

A run of one sitting leaves `protocol.json` beside `protocol.html`, in the directory you named
with `--out`: that is what a cloning run and a read-back run each write. A library run is
several sittings, so it leaves one `project.json` in the `protocol/` folder, holding the whole
chain of pages.

The file is plain text, and the words in it are the words on the page. Search it for a sentence
you can see and you have found the field holding it.

## Change the words, not the numbers

| Yours to change | Leave to the package |
| --- | --- |
| the wording of any instruction | a volume, a mass or a concentration |
| the order of the steps | a count of tubes, wells, colonies or cycles |
| a step of your own, with its own instructions | a temperature, or a time the design fixed |
| a note or a caution | a band size on a gel, an overhang, a primer or its melting temperature |
| the title, and the sentences at the top | a verdict on a badge or an oligo |

Everything in the right column came out of the design. Typing a different number in does not
change the design. It only makes the page disagree with the plasmid the command built, and
nobody at the bench can tell which of the two is right.

Dropping a step leaves mentions of it behind: a row in the materials, a line in another step's
troubleshooting, an oligo nothing now uses. Nothing checks for those, so read the whole page
again after an edit.

## Make the page again

```bash
pixi run mbio protocol render plan/protocol.json
```

It prints the file it wrote:

```text
plan/protocol.html
```

Hand it a folder instead and the whole chain is made again, one page a protocol:

```bash
pixi run mbio protocol render library/protocol
```

```text
library/protocol/project.json
library/protocol/01-part-carrier.html
library/protocol/02-primer-plates.html
library/protocol/03-cargo-ordering-and-pool-preparation.html
library/protocol/04-cargo-creation.html
library/protocol/05-cargo-validation-barcode-ligation.html
library/protocol/05-cargo-validation-index-pcr.html
library/protocol/06-library-assembly-in-rounds.html
library/protocol/07-final-cargo-ligation.html
library/protocol/index.html
library/protocol/reagents.html
library/protocol/references.html
```

Each page is written over. Add `--output` to send a single page somewhere else and leave the
original where it is.

## What it refuses

Every value is checked for its type, and a wrong one is named where it sits:

```text
error: plan/protocol.json: protocol.steps[3].instructions: expected a list, got a string
```

Here a list of instructions had been replaced by one string. The render stops and writes
nothing, so a slip shows up now rather than on the page you carry to the bench.

## Keep it where a rerun will not land on it

Planning again into the same `--out` directory writes over `protocol.json` and the page with it,
and your edits are gone. The same inputs always write the same bytes, so a rerun puts back
exactly what the command wrote the first time.

Copy the edited directory somewhere else, or plan into a fresh one. Either way, the edited
protocol is yours to keep and nothing here writes over it.

## When the change is a design change

Swapping the enzyme, the polymerase, the inserts or which way round one sits is not an edit to
the page. Those go back through the plan, so everything resting on them is worked out again: the
primers, the overhangs, the product, the bands on the gel. Run the command again with the option
you want. Its `--help` lists them.

## Why the rule is worth keeping

The protocol is data so that it can be adapted. A lab asks for things no pipeline anticipated,
and you, or an agent working for you, should be able to answer without waiting for the package
to grow an option for each one.

That only works while one thing computes the numbers. If a volume on the page could have come
from either the package or a person, no reader can tell which, and the page can drift away from
the design it was built from without anyone noticing. Leaving every computed number in the
package's hands is what makes the page safe to follow.

[How to read a protocol page](reading.md) says what each of those numbers is, and
[Protocols](index.md) lists the pages you can open now.

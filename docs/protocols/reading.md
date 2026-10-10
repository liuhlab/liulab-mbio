# How to read a protocol page

One page is the whole of a run, read top to bottom. This walk-through uses the Golden Gate page
the repo ships: [GFP into pUC19](../examples/pUC19-GFP/protocol.html), eleven steps. Open it
beside this page and follow along. Nothing on it was typed by hand, so every number on it came
from the design the command built.

## The top of the page

The title and one sentence say what the run does. Under them sit eight cards: vector, insert,
enzyme, fragments, overhangs, fidelity, product and selection. Here they read `pUC19, 2686 bp`,
`GFP, 717 bp`, `BbsI-HF (R3539) at 37 °C`, `2 in one reaction`, `ATGA and TGGC`,
`100%, measured`, `pUC19-GFP, 3347 bp` and `ampicillin or carbenicillin`.

Below a thin line, under the heading **Keep in mind**, four sentences follow, about this clone in
particular. GFP reads on the opposite strand from the lac promoter. No ribosome binding site is
annotated ahead of it, so the clone is not expected to glow. The insertion interrupts lacZα, so
correct colonies are white. All four are read off the finished plasmid, so they cannot disagree
with the map.

## The badges

Another thin line, then a row of badges: `sites`, `pUC19 backbone`, `GFP`, `junctions` and
`primers`. Four are green. `primers` is amber, and only a badge that is not green prints its line.
That line says nine oligos were designed, two carry a warning and none fails.

A badge can also read **not judged**. That is not a pass and not a failure. It means nothing
judged it, because no sourced threshold fits what this run does.
[The library page](../examples/ap1-library/protocol/index.html) carries one on
`ligation fidelity`: the measurement that exists was made with a different ligase in a different
buffer, so the package will not score the round from it. Treat it as a call you have to make.

## Materials and the oligo sheet

**Materials** is everything the run consumes: name, supplier, catalogue number, storage, how
much per run, and a note. Fifteen rows here, with the equipment named on a line under them.

**Oligos** is the order sheet. Nine rows, each with its sequence, length, melting temperature,
what it is for, its working stock and its own verdict. A copy button sits on every sequence, and
one at the foot copies all nine. `Order sheet: primers.tsv` names the file those nine also sit
in, ready for a vendor's form. A fold under the table reads `Checks on 2 of 9 oligos`. Open it
and each warning is named: `Tm 64.1 °C (band 60-64)` on one, `GC 35% (band 40-60)` on the other.

**A `band` here is a range, not anything on a gel.** `band 60-64` is the melting temperature the
check wants to see. Further down this page a gel's bands are DNA. On a library page a third kind
turns up in the bill: a price band, the quantity bracket a quote is charged in. Say which you
mean.

## The steps

Eleven numbered steps run from the two PCRs to sequencing. A column on the right lists them,
grouped by section, and stays in view as you scroll: click a step to jump to it. Each section
counts its own ticks, and the one you are in opens when you come back. On a narrow screen the list
sits above the title instead.

Each step and each instruction has a tick box. Tick every instruction of a step and the step ticks
itself; clear one and it clears. A step partly done shows a dash in its box. Tick the step and all
its instructions tick with it. The bar at the top counts steps: `0 of 11 steps done`. The page
remembers your ticks, your reaction counts, the concentrations you type and your running timers,
so closing the tab does not lose the sitting. `Reset page` puts them all back to what was written.
`Print` lays the page out for paper.

A step holds some of these, always in this order: cautions, instructions, a reaction table, a
thermocycler program, timers, a map, the expected result, troubleshooting, and notes. A note is the one
place the page says why, and it sits last so the doing comes first.

## Reaction tables

A reaction table has a `Reactions` box. Type how many tubes you are running and every volume
rescales. `1 rxn (µL)` is one tube; `Mix for N (µL)` is the master mix, which already includes
10% extra. A line under the table says how to split it, such as
`Put 49 µL of mix in each tube, then add 1 µL Template DNA`, because the template is the one
thing that differs tube to tube.

Step 6 measures each fragment, and step 7's table takes what you measured. Its two DNA rows
have a field in the `Stock` column: `82 ng/µL` for the pUC19 backbone. That is the
concentration the planned 1 µL assumes. Type the one you measured and the row gives the volume
that carries the same picomoles. The water takes up the difference, so the reaction stays 15
µL. The page keeps your number, and a button beside it offers the planned one back. Once you
tick anything in step 7, the field holds still, so no volume moves while you pipette. Clear the
ticks to change it. If the DNA no longer fits, a line under the table says so in step 6's own
troubleshooting words.

**Two senses of `reaction` sit on one page.** In a reaction table, one reaction is one tube of
that mix, so a step done in 72 wells says 72. In the cards at the top, `2 in one reaction` means
two fragments joined in a single tube — and in that sense one reaction spread over a 96-well
array is still one reaction, not 96.

## Thermocycler programs

A program is a table of step, temperature, time and cycles, with the whole run length in its
caption: `32 min 30 s plus ramps`. A small bracketed number beside a row says where that row came
from, such as `[1]` on the cycle count of the first PCR. It jumps to the source with that number
in the list at the foot of the page. Hold the pointer over it to see the document and the part
cited, here FAQ 11 of the kit manual.

## Timers

An incubation of fixed length gets a timer, such as `DpnI digest 1:00:00`. Press `Start` to
start it. A thermocycler program gets one under its table, set to the run time in its caption.
That time leaves out the ramps, so the machine ends a little later. A time written only in
words, such as an overnight growth, gets no timer. Step 9 carries five: 10 minutes to thaw, 30
minutes on ice, a 30-second heat shock, 5 minutes back on ice, and an hour of outgrowth.

The time is yours to change. Click it, type a new one and press Enter: `1:30:00`, `90:00`, or
`90` for 90 minutes. On a paused timer, what you type is the time left. The page keeps your
time, and a button beside it offers the planned time back. `Reset page` puts it back as well.

A running timer keeps its deadline and a paused one keeps the seconds left, so turning the page
does not lose an incubation. When it runs out it beeps, vibrates where the device can, flashes,
and puts `Time up` in the tab's title. A browser lets a page make a sound only after you have
touched it. After you reload a page with a timer running, click anywhere on it once.

## Maps

A step that makes a molecule shows it as a map, in a box of its own with its caption inside.
Step 8 draws pUC19-GFP as the reaction closes it: the two junctions in colour, the rest grey.
Press `Explore the map` and the whole plasmid opens in its place. Switch features, primers and
cut sites on and off, zoom in, or show the bases. `Close the map` puts the small one back, and
`Open it in a new tab` gives the map a window of its own. Paper prints the small map.

## Expected result and troubleshooting

Every step ends with what you should see, written as an observation: `One band at 2662 bp.`
after the first PCR, `Clean DNA, free of polymerase, primers and dNTPs.` after the spin columns.
A step that shows you nothing at the time says so, and says when it will show. Troubleshooting
follows: a short list of a problem and what to do about it. Read it before the step, not after
it has gone wrong.

## The simulated gel

Where a step is read on a gel, the expected result draws the gel the design predicts: a ladder
lane, then one lane a sample, each band at the size it should run at. At the colony PCR, a
correct clone gives bands at 160 and 838 bp and an empty vector gives one at 177 bp, against the
NEB 100 bp ladder. Read the big band first, since that is the one telling the two apart. The
drawing is a prediction, not a photograph.

## Where the numbers came from

The page closes with two lists. **Sources** names each document a number was read from, with its
edition and the day it was read. Each is numbered as the bracketed numbers in the text name it,
and says which parts of the document the page cites, so a number still leads somewhere on paper.
**References** is the wider reading behind the method. Where the package cannot source a number, it says so rather
than printing a guess, and a chain collects those in a `Holes` list on its way-in page.

[Change a protocol](editing.md) says what you may edit here, and what you must leave to the
package.

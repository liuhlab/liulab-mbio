---
search:
  exclude: true
---

# A plan is finished with the holes no source can fill still standing

A number nobody published is left as a hole rather than guessed, so "no hole left" cannot be what
finishes a design. Some holes close only at a bench, and some only when a lab names its own stock.

The AP-1 demo's final assembly is where this first bites. Two holes stand there. **H24** asks the
working vector's mass and the ratio it meets the cargo at; departure D12 purifies nothing between
the release and the assembly, so the cargo is never pipetted into anything and the pot already
holds 50 µL. NEB's 15 µL two-fragment reaction, Qian's column-cleaned amplicon assembly and the
paper's own room-temperature T7 ligation each size a different tube, and none of them this one.
Only a pilot settles it. **H31** asks which backbone the library ends in, which is a stock a lab
holds and an application chooses.

So a plan is finished when every hole left is of that sort: `undecided` where only a bench
settles it, `lab` where only the user's own stock does. A hole a source would close is unfinished
work, and the step that carries it is where to look.

The demo therefore names no working vector. Naming one closes H31 and leaves H24 standing, which
`test_a_named_working_vector_fills_the_enzyme_and_its_cycling_in` pins, so it buys no finished
step. It also wants a backbone `sites.domesticate` will not produce, the one site it must clear
being non-coding, and it moves four of the five files `tests/synbio/igga/test_gate.py` reads as
its corpus — changing what the gate is proved against, for a demo customised to one lab's stock.

## Considered options

- **Hold a plan open until no hole remains.** It asks for a number nobody published, which is the
  defect a hole exists to prevent.
- **Name a working vector in the demo anyway.** The measurement above prices it: a hand-edited
  backbone the package cannot write, a CLI option that does not exist, and most of the gate's
  corpus moved.
- **Count only `undecided` holes, leaving `lab` ones unmentioned.** The two differ in who fills
  them, not in whether the step is runnable as shipped, and a reader at the bench meets both.

## Amendment: a build may state what only a bench or a shelf settles

A lab that has run the pilot is not held to a hole. A build may state the assembly's mass and
ratio, the two pool cycle counts, the linkage pass mark and the working vector, each optional
and each leaving its hole standing by default. Stated, the step prints it as that run's own
measurement, never a published figure. The measurement and the five kinds stand.

## Consequences

A pipeline's acceptance reads each hole's kind rather than counting them. The AP-1 demo ships two
on its final assembly and is finished; a third of either kind would be a finding, and a hole
naming a source nobody read would be a defect.

A project of the user's own that names a working vector drops H31 and keeps H24 until it states
the assembly too, and a demo stays what a reader copies rather than a record of one lab's shelf.

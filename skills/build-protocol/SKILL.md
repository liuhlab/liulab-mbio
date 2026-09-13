---
name: build-protocol
description: >-
  Write a bench protocol as a JSON data file and render it to one self-contained interactive
  HTML page: overview, materials, numbered steps with check marks, reaction tables that rescale
  by reaction count, thermocycler programs, timers, expected results, troubleshooting and a
  simulated agarose gel. Use whenever someone wants a protocol, a bench sheet, a wet-lab
  procedure, a printable step-by-step for an experiment, or a page someone can follow at the
  bench — and whenever another skill has to hand back a protocol rather than prose.
---

# Build a protocol

The model and the renderer are `liulab_mbio.protocol`; they hold no task-specific content. Your
job is what goes in the protocol. Write the data, render it, read the page back.

## Render

```bash
pixi run liulab_mbio protocol render protocol.json -o protocol.html
```

From Python, when another module computes the numbers:

```python
from liulab_mbio.protocol import Protocol, read_protocol, render_html, write_html

write_html(read_protocol("protocol.json"), "protocol.html")   # from a data file
write_html(Protocol(title="...", steps=(...)), "protocol.html")  # built in code
```

The page needs no network, keeps check marks in the browser, rescales its tables and prints.
Never hand-edit the HTML: change the data and render again.

## The data format

JSON keys are the field names of the classes in `liulab_mbio.protocol.model`, lists become
tuples, and only fields without a default are required. Two places to look:

```bash
pixi run python -c "import liulab_mbio.protocol as p; print(p.__doc__)"   # the outline
```

`src/liulab_mbio/protocol/model.py` — each class docstring says what its fields mean.
`tests/data/pcr-protocol.json` — a worked example.

Text is escaped, so write plain sentences; markup will show up as characters.

## What a good protocol contains

Picture the reader: standing at the bench, gloves on, your page open, unable to ask you
anything. Everything they need to decide what to do next is on the page, or it does not exist.

- **Summary and overview.** One paragraph, then the facts they check before starting: what goes
  in, what comes out, how long it takes.
- **Materials.** Every reagent, kit, oligo and piece of equipment, with storage. Oligos carry
  their sequence, which gives the reader a copy button.
- **Steps in the order they happen**, each one action per instruction, in the imperative. A
  step the reader could split in two usually should be.
- **Volumes in a reaction table**, never buried in a sentence: the reader scales it to their
  reaction count and pipettes from it. Mark anything added to each tube separately, such as
  template, as not part of the master mix.
- **Thermocycler programs as programs**, with stages, cycles and the lid temperature.
- **Timers** for any wait the reader has to track.
- **Expected result at every step**, so the reader knows whether to carry on. Where the answer
  is a gel, give the gel: band sizes per lane and a named ladder.
- **Troubleshooting** for what actually goes wrong at that step, each with what to do.
- **Cautions** before a step, **notes** after it.
- **References** for numbers you took from a vendor protocol or a paper.

## Before you hand it over

Render the file and read the HTML back. Check that every step has an expected result, that
volumes and times are numbers rather than "some" or "a while", and that nothing on the page
asks the reader a question you should have answered.

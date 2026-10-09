# Projects

A project is a named method this lab runs end to end, built on the methods this site already
describes. It is a long run, not one reaction: one file goes in, and a folder of bench pages
comes out, one page for each sitting. The method fixes the molecules, the enzymes and the
order; a run chooses only its own numbers.

## What earns a page

A project earns a page when it has a pipeline entry that writes files someone takes to the
bench. Nothing else belongs here. A technique with no entry of its own is a
[method](../methods/index.md), and a run with no files to hand over is a lab notebook entry.
So this section holds two pages and will not grow quietly.

## The two

| Project | What it is for | What runs it |
| --- | --- | --- |
| [iGGA](igga.md) | joining lists of proteins into every combination, as one plasmid library, each part barcoded so sequencing says what a member carries | `pixi run synbio igga plan PROJECT --out DIR` |
| [DMX](dmx.md) | reading a plate of picked wells back, so every well is named rather than assumed | `pixi run synbio dmx plan BUILD --out DIR` |

Both commands come with the same install as `mbio`, under the second command name `synbio`.
Each page runs its own worked example and links every file that run wrote.

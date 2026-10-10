# Projects

A project is protocols run in order, each handed what the ones before it produced. It is built
on the methods this site already describes. One file goes in and the bench pages come out: a
folder of them for a library built in rounds, one page where the whole run is a single
procedure. The molecules, the enzymes and the order are fixed in code; the file holds only what
this run chooses.

## What earns a page

A project earns a page when it has a pipeline entry that writes files someone takes to the
bench. Nothing else belongs here. A technique with no entry of its own belongs under
[Methods](../methods/index.md), and a run with no files to hand over is a lab notebook entry.
So this section holds two pages and will not grow quietly.

## The two

| Project | What it is for | What runs it |
| --- | --- | --- |
| [iGGA](igga.md) | joining lists of proteins into every combination, as one plasmid library, each part barcoded so sequencing says what a member carries | `pixi run synbio igga plan BUILD --out DIR` |
| [DMX](dmx.md) | reading a plate of picked wells back, so every well is named rather than assumed | `pixi run synbio dmx plan BUILD --out DIR` |

Both run under `synbio`, which comes with the same install as `mbio`. Each page runs its own
worked example and links every file that run wrote.

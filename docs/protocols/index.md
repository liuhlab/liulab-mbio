# Protocols

A protocol page is one HTML file. It needs no network and nothing installed: open it in a
browser at the bench and work down it. Tick a step and the page remembers. Change how many
reactions you are running and every table rescales. Start an incubation and the timer keeps its
deadline. A command wrote every page below, and nothing on them was typed by hand.

The catalogue runs in the order a job does: design and order, build, validate.

## Design and order

| Page | What it covers |
| --- | --- |
| [Domesticating the working vector](../examples/ap1-library/working-vector-domestication.html) | taking six BsaI and two BsmBI sites out of a backbone, then confirming what arrives |
| [Cargo ordering and pool preparation](../examples/ap1-library/protocol/03-cargo-ordering-and-pool-preparation.html) | ordering a pool of 131 oligos and the 74 primers that pull batches out of it |
| [Primer plates](../examples/ap1-library/protocol/02-primer-plates.html) | ordering 74 primers as one plate, then splitting a working copy to pipette from |
| [Reagents and equipment](../examples/ap1-library/protocol/reagents.html) | what a whole run orders, priced from a list the lab holds |

## Build

| Page | What it covers |
| --- | --- |
| [Golden Gate assembly: GFP into pUC19](../examples/pUC19-GFP/protocol.html) | one cloning run end to end: eleven steps, from two PCRs to sequencing |
| [Part carrier](../examples/ap1-library/protocol/01-part-carrier.html) | putting each of 72 parts in its own plasmid, one part a well, so any can be released again |
| [Cargo creation](../examples/ap1-library/protocol/04-cargo-creation.html) | pulling each block out of the pool and cloning it into its position's vector |
| [Library assembly in rounds](../examples/ap1-library/protocol/06-library-assembly-in-rounds.html) | three Golden Gate rounds, one part list a round, over 27 steps |
| [Final cargo ligation](../examples/ap1-library/protocol/07-final-cargo-ligation.html) | moving the finished library into the working vector an application needs |

## Validate

| Page | What it covers |
| --- | --- |
| [Cargo validation: barcode ligation](../examples/ap1-library/protocol/05-cargo-validation-barcode-ligation.html) | taking archived designs to clonal wells and marking each so sequencing names it: one of two routes |
| [Cargo validation: index PCR](../examples/ap1-library/protocol/05-cargo-validation-index-pcr.html) | the same job marked by index PCR, the other route. A run picks one; nobody does both |
| [AP-1 cargo read-back](../examples/ap1-readback/protocol/index.html) | the same job as a run of its own, over a plate of frozen stock and a list of names |
| [Design read-back: index PCR](../examples/ap1-readback/protocol/01-design-read-back-index-pcr.html) | its bench page: array and pick, mark every well, then call each one |

## Work through a run of many sittings

A cloning run is one page. A library run is a chain of them: one page a sitting, in the order
someone does them, each saying what it is handed and what it leaves behind. Mail one page to
whoever runs that sitting, or send the folder and they have the run.

[The way in](../examples/ap1-library/protocol/index.html) draws the chain, lists the checks on
the design, and collects every number the run has no source for.
[The references](../examples/ap1-library/protocol/references.html) say where each number was
read from, and every page in the chain shares them.

## Three more methods write one

Gibson assembly, restriction and ligation, and Gateway cloning each write a protocol of the same
shape: the same cards, badges, numbered steps, reaction tables and drawn gels.
[Gibson](../methods/gibson.md), [restriction and ligation](../methods/restriction-ligation.md)
and [Gateway](../methods/gateway.md) each run the command and link the page it wrote.

## Read one, or change one

[How to read a protocol page](reading.md) walks the Golden Gate page part by part, and names the
two words a biologist is most likely to read the wrong way.

[Change a protocol](editing.md) covers the rest: how to skip a step, add one or reword one, what
you must leave alone, and how to make the page again afterwards.

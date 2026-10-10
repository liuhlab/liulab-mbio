---
search:
  exclude: true
---

# Research notes

A note is added to this index in the commit that writes the note. Nothing checks that, so it is
on you.

Every note here is agent-facing: out of the site navigation, excluded from search, and exempt
from the jargon and readability rules a published page is held to. A human-facing page never
links to one. Where a bench reader needs a fact a note holds, the fact is restated on a page
they can reach.

The directory is flat and stays flat. A note is cited by full path from `src/`, `tests/`,
`scripts/` and the committed example corpus, so moving one is expensive; and the obvious tree,
one directory per method, would file the evidence `mbio.ligase` and `mbio.overhangs` cite under
a method's name, which inverts the package boundary. The grouping belongs in this index, where a
note can sit under two headings at once. Within the directory, a note in a family takes the
family's name as a prefix, so `ls` clusters it at no cost.

Thirty-nine notes, under nine headings. Four of them earn a second heading, listed at the end.

## Shipped data and where it came from

| Note | What it settled |
| --- | --- |
| [`codon-usage.md`](codon-usage.md) | where the shipped codon tables come from, which published tables were refused, and why these are counted from sequence |
| [`restriction-enzyme-data.md`](restriction-enzyme-data.md) | where the shipped enzyme records come from, what each source allows, and which value is read from whom |
| [`ligation-fidelity.md`](ligation-fidelity.md) | the fidelity matrices' source and licence, which way round a matrix reads, and the rules applied where no measurement covers an enzyme |
| [`fonts.md`](fonts.md) | the DejaVu release behind the shipped faces, what its licence asks, and how closely the measured tables match a browser |
| [`dmx-destination.md`](dmx-destination.md) | what the build script does to DMX0001, why each step is there, and what the result measures |
| [`working-vector-plvx-tetone.md`](working-vector-plvx-tetone.md) | pLVX-TetOne-Puro-GFP measured against the four requirements a lentiviral parent has to meet |
| [`synthesis-and-assembly-barcode-kit.md`](synthesis-and-assembly-barcode-kit.md) | the 96 barcodes, their four groups, the overhangs they chain on, the two universal primers, and the accession range |

## Cloning methods

| Note | What it settled |
| --- | --- |
| [`golden-gate-assembly.md`](golden-gate-assembly.md) | NEBridge's reaction and cycling, overhang fidelity, and how the other assembly methods compare |
| [`gibson-assembly.md`](gibson-assembly.md) | the assembly products, the two oligo routes, and every number the Gibson pipeline carries |
| [`restriction-ligation.md`](restriction-ligation.md) | the digests, the buffers, what each source allows, and what nobody stated |
| [`gateway-cloning.md`](gateway-cloning.md) | the att sites, the BP and LR reactions, and which sequences may ship |
| [`dad-split.md`](dad-split.md) | what DAD is, what the two existing implementations do, and where the line falls between them and this package |

## Primers, PCR and barcodes

| Note | What it settled |
| --- | --- |
| [`primer-design-and-pcr.md`](primer-design-and-pcr.md) | primer design, how Tm is computed, and what validates a colony by PCR |
| [`oligo-pool-pcr-cycles.md`](oligo-pool-pcr-cycles.md) | the cycle counts for the two oligo-pool PCRs, which had looked sourced and were not |
| [`bench-numbers.md`](bench-numbers.md) | every number the generated protocol carries and cannot compute, stage by stage |
| [`barcode-design.md`](barcode-design.md) | whether a barcode set needs a GC band and a homopolymer cap beyond distance and site freedom |
| [`heat-inactivation-order.md`](heat-inactivation-order.md) | whether one 80 °C hold inactivates an enzyme the supplier lists at 65 °C, or the tube needs both |

## What a vector tolerates, and whether it survived

| Note | What it settled |
| --- | --- |
| [`domestication-methods.md`](domestication-methods.md) | which route domesticates a working vector, and when each one is the right one |
| [`lentiviral-tolerance.md`](lentiviral-tolerance.md) | which parts of a third-generation lentiviral transfer plasmid may be changed at all, and what each change costs |
| [`promoter-readthrough-lentivirus.md`](promoter-readthrough-lentivirus.md) | what a primary source measures when a constitutive promoter sits upstream of an inducible one inside the LTRs |
| [`vector-qc-panel.md`](vector-qc-panel.md) | how you find out whether a domestication broke the vector: which assay per failure mode, and how small a change it catches |

## The library method

| Note | What it settled |
| --- | --- |
| [`synthesis-and-assembly.md`](synthesis-and-assembly.md) | the evidence behind the iGGA method, every source and when it was read, and the method as specified |
| [`synthesis-and-assembly-departures.md`](synthesis-and-assembly-departures.md) | every way our design differs from the paper, what forced each change, and whether the reason survives checking |
| [`synthesis-and-assembly-materials.md`](synthesis-and-assembly-materials.md) | what you need for each step, by the process that uses it |
| [`protein-library-assembly.md`](protein-library-assembly.md) | the published method behind the pipeline, its numbers measured from the supplement, and its licence |
| [`ap1-demo-project.md`](ap1-demo-project.md) | one real build specified end to end, so the pipeline had something concrete to build toward |

## Reading a library back

| Note | What it settled |
| --- | --- |
| [`route-choice.md`](route-choice.md) | how a reader picks between barcode ligation and index PCR, or neither, and what a run that reads nothing back holds |
| [`route-b-index-pcr.md`](route-b-index-pcr.md) | the per-well numbers of LevSeq's index PCR, the one stage the bench-numbers inventory left out |
| [`route-b-index-primers.md`](route-b-index-primers.md) | where the index primers bind and what marks a well, against a premise that turned out wrong |

## Drawing a record

| Note | What it settled |
| --- | --- |
| [`plasmid-map-packages.md`](plasmid-map-packages.md) | which Python packages draw a map, and how close each comes to SnapGene. The answer was none, which is why `plot/` exists |
| [`label-placement.md`](label-placement.md) | the one rule that keeps every label clear of every other label and of the drawing, on all three layouts |
| [`feature-colours.md`](feature-colours.md) | where a default colour for a feature type can come from under a licence this package may ship |
| [`figure-sources.md`](figure-sources.md) | what four published figures show, and what ours has to vary |
| [`offline-html-maps.md`](offline-html-maps.md) | the ways to write a record as one HTML file that opens offline, and what each costs |
| [`offline-navigation.md`](offline-navigation.md) | how an offline document navigates with no server, measured in a browser |

## Writing for the bench

| Note | What it settled |
| --- | --- |
| [`bench-facing-documentation.md`](bench-facing-documentation.md) | how five comparable tools document themselves for a wet-lab reader, and which page shape to take from each |
| [`docs-site-diagrams.md`](docs-site-diagrams.md) | whether the docs site renders a mermaid fence, what the theme gives unconfigured, and what shape the orientation schematic takes |

## Surveys whose answer was no

| Note | What it settled |
| --- | --- |
| [`dnachisel-evaluation.md`](dnachisel-evaluation.md) | what DnaChisel and the Edinburgh Genome Foundry stack would replace here, and why little of it is worth adopting |
| [`gantt-conventions.md`](gantt-conventions.md) | whether a ready convention charts a protocol whose durations are mostly unknown |

## Four notes that belong under two headings

This is the reason the directory has no subdirectories: a tree makes you pick one, and an index
lists a note twice for nothing.

- [`fonts.md`](fonts.md) is shipped data, and it is also how a map draws text.
- [`ligation-fidelity.md`](ligation-fidelity.md) is shipped data, and it is also what makes a
  Golden Gate overhang set choosable.
- [`synthesis-and-assembly-barcode-kit.md`](synthesis-and-assembly-barcode-kit.md) is shipped
  data, and it is also how a well is marked for reading back.
- [`figure-sources.md`](figure-sources.md) is about drawing, and it is also about the library.

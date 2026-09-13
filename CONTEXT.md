# Context

## Glossary

### Sequence record

One DNA sequence with its topology, features, primers and notes: what a sequence file holds,
whatever its format.
_Avoid_: SeqRecord, file, map

### Feature

A named, typed annotation of a sequence record, such as a CDS or a promoter, lying on a strand
over one or more segments.
_Avoid_: annotation, region

### Segment

One contiguous span of a feature. Most features have one; a promoter drawn with its -35 and -10
boxes has several.
_Avoid_: range, interval, part

### Primer

A named oligonucleotide, written 5' to 3', that primes DNA synthesis. A sequence record
annotates each of its primers with its binding sites.
_Avoid_: oligo

### Binding site

Where a primer's 3' part anneals to a sequence record, and on which strand. A 5' tail, such as
an enzyme site added for cloning, lies outside it.
_Avoid_: primer site, hybridization site

### Across the origin

Of a span on a circular sequence record: running through the last base and on into the first.
_Avoid_: wrapping, wraparound

### Protocol

A bench procedure someone can follow without asking anything: summary, materials, numbered
steps, expected results and references. Rendered as one self-contained HTML page.
_Avoid_: SOP, recipe, method

### Step

One numbered unit of a protocol: what to do, what it needs, what a successful result looks
like, and what to do when the result is wrong.
_Avoid_: task, procedure

### Reaction table

The volumes of one reaction, scaled to a master mix for several reactions plus an overage. A
component marked per-tube, such as template, is added to each tube instead of to the mix.
_Avoid_: mix table, recipe

### Thermocycler program

Stages run in order, each a list of incubations repeated for a number of cycles.
_Avoid_: cycling conditions, PCR conditions

### Simulated gel

A drawing of the agarose gel a step should produce: a ladder lane and sample lanes, each band
placed by the logarithm of its size in base pairs.
_Avoid_: virtual gel, gel image

### Ladder

The size marker lane of a simulated gel: a name and the band sizes it carries.
_Avoid_: marker, standard

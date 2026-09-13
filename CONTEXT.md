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

### Recognition site

The sequence a restriction enzyme binds, written 5' to 3' in IUPAC codes.
_Avoid_: recognition sequence, restriction site

### Cut offset

Where an enzyme cuts one strand, counted from the first base of its recognition site: a 0-based
boundary, so offset 7 severs the strand between the seventh and eighth base. An offset outside
the site is what makes an enzyme Type IIS. The same two cuts placed in a record's own
coordinates are its cut positions.
_Avoid_: cleavage coordinate, cut index

### Overhang

The bases two staggered cuts leave single-stranded. The end is 5' when the bottom strand is cut
further along than the top, 3' when it is cut nearer, and blunt when neither strand overhangs.
Golden Gate joins two fragments by matching overhangs.
_Avoid_: sticky end, extension

### Type IIS

A restriction enzyme that cuts outside its recognition site, so the overhang it leaves is
whatever the DNA holds there rather than part of the site. BsaI, BsmBI, BbsI, SapI, PaqCI and
BtgZI are the ones this package ships.
_Avoid_: type 2S, IIs

### Isoschizomer

An enzyme a supplier sells that reads the same site as another and cuts it in the same place. A
neoschizomer shares the site but cuts elsewhere — XmaI to SmaI — and is not an isoschizomer.
_Avoid_: equivalent enzyme, alias

### Commercial name

The name a supplier sells an enzyme under, such as BsaI-HFv2 for BsaI. It carries the
supplier's own incubation, heat-inactivation and methylation answers, which differ between
products of one enzyme.
_Avoid_: brand, product name

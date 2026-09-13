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

### Annealing region

The 3' part of a primer that pairs with its binding site. Length, GC content and Tm are read
from it, not from the whole primer.
_Avoid_: priming region, homology arm

### Tail

The 5' bases of a primer that pair with nothing on the template, such as a spacer, an enzyme
site and an overhang. Copies made from the primer carry them, so later cycles anneal it whole.
_Avoid_: overhang, flap, adapter

### Melting temperature

Tm: the temperature at which half of a primer is paired with its complement, in a named
polymerase's buffer. Of the annealing region unless said otherwise.
_Avoid_: melting point, annealing temperature

### Annealing temperature

Ta: the temperature of a PCR's annealing step, set from the lower Tm of the pair by a rule
belonging to the polymerase.
_Avoid_: hybridization temperature, Tm

### 3'-anchored dimer

A self-dimer or heterodimer that pairs a primer's 3' end, which the polymerase can extend into
primer dimer. Worse than one pairing only inside the primer.
_Avoid_: end dimer, self-complementarity

### Off-target site

A place other than its binding site where a primer's annealing region can prime: its 3' end
pairs, and it melts near the temperature a perfect match would.
_Avoid_: mispriming site, secondary binding site

### Amplicon

What a primer pair copies, from one binding site to the other, tails included.
_Avoid_: PCR product, band

### Check

One pass, warn or fail verdict on a primer or a pair, carrying the value it judged.
_Avoid_: validation, rule, test

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

### Cut site

Where an enzyme reads its recognition site in a record: the strand carrying it, where each
strand is severed, and the overhang that leaves. A site is reported wherever the enzyme may
cut, so an IUPAC code in the template that could spell the site counts as one.
_Avoid_: hit, match, restriction site

### Digest fragment

One piece a digest leaves, bounded by two cuts in the top strand and carrying an overhang at
each end. Its length is the top strand's, which is what a gel measures. An uncut circular
record leaves none, nothing having been cut; an uncut linear one is a single fragment already.
_Avoid_: band, piece

### Domestication

Taking a recognition site out of a coding sequence by changing one codon for another spelling
the same amino acid, so the reading frame and the protein survive. The replacement is the one
the host uses most often among those that create no new site. A site outside a coding sequence
cannot be domesticated, and is reported for someone to decide about.
_Avoid_: silent mutation, site removal

### Codon usage

How often a host spells each of the 64 codons, counted over every complete coding sequence of
its genome. A whole-genome table is the background the genome itself uses; a highly expressed
table counts a small reference set and is a different thing.
_Avoid_: codon bias, codon frequency table

### Packet

One unit of a SnapGene `.dna` file: a type byte, a length, then the payload. A reader keeps the
packets the sequence record does not hold, so that a writer can put them back unchanged.
_Avoid_: block, chunk, section

### Stale packet

A packet describing bases that have since changed, such as SnapGene's cut-site cache or its
history. The writer drops it rather than write it back, leaving SnapGene to rebuild it.
_Avoid_: invalid packet, dirty cache

### Edit report

What an edit did to the features and primer binding sites it did not simply shift: **trimmed**
ones lost bases, **dropped** ones lost all of them, and **changed** ones now span the new bases
because the edit fell inside them.
_Avoid_: diff, changelog, summary

### Molar ratio

How many molecules of insert an assembly gets for each molecule of vector. Counted in moles and
not in mass, because a short insert weighs less at the same ratio: a protocol prints pmol beside
ng for that reason.
_Avoid_: insert ratio, stoichiometry

### Junction

Where two fragments meet in an assembled product: the four bases, three for SapI, that one
overhang paired with its match. Validation reads across it.
_Avoid_: joint, seam, fusion site

### End soak

The 60 °C step every Golden Gate program ends with. It is a digest and not heat inactivation: it
cuts destination plasmid that never opened or has closed again, so fewer empty colonies grow.
_Avoid_: final incubation, kill step

### Colony PCR

A PCR templated by a colony lifted off the plate with a toothpick, which the long hot step at
the start lyses. It says which clone a colony carries without a plasmid prep.
_Avoid_: screening PCR, colony screen

### Reversed insert

A clone carrying its insert the other way round. Two vector primers flanking the insert give it
the same band as the correct clone, so a gel separates the two only with a third primer inside
the insert.
_Avoid_: flipped clone, wrong orientation

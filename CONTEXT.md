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
annotates each of its primers with its binding sites. Every primer is an **Oligo**; say primer
only for one that primes something.
_Avoid_: oligo (the wider word, and its own entry)

### Oligo

A short synthetic DNA a protocol orders: a row of the order sheet, with the name it is ordered
under, the sequence, the length and the Tm. It is what a bench buys rather than what a reaction
does, so a protocol may order one that primes nothing.
_Avoid_: primer (unless it primes), synthetic DNA

### Binding site

Where a primer's 3' part anneals to a sequence record, and on which strand. A 5' tail, such as
an enzyme site added for cloning, lies outside it.
_Avoid_: primer site, hybridization site

### Placement

Where a primer's binding site may lie: the positions each of its two ends may take. A tail fixes
the 5' end, a target to read across bounds each end's distance from it, and a region bounds the
amplicon.
_Avoid_: window, flexibility, allowance

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

### Off-target amplicon

An amplicon a pair, or one of its primers alone, makes on a genome other than the intended one.
Where no amplicon is named as intended, every amplicon is one.
_Avoid_: non-specific product, secondary product

### Check

One pass, warn or fail verdict on a primer, a pair or an assembled product, carrying the value
it judged. A check no sourced threshold judges carries no verdict at all rather than a pass,
and says so wherever it is shown.
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

How often a host spells each of the 64 codons, counted over the complete coding sequences of
its genome, one transcript per gene where a gene has several. A whole-genome table is the background the genome itself uses; a highly expressed
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

### Cloning method

How an experiment joins its fragments into one plasmid: the mechanism, the enzymes it needs, and
what it leaves at each junction. It is chosen for a job — how many fragments join, whether the
junction may gain bases, what the parts already carry, speed and cost — rather than preferred in
general. Golden Gate is the one this package plans; a library built in rounds is a pipeline over
a method, not a method of its own.
_Avoid_: cloning strategy, technique, approach

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

### Overhang set

The overhangs an assembly's junctions take, together. Chosen together and not one at a time,
because every rule and the fidelity score are about the whole set.
_Avoid_: fusion sites, overhang list

### Ligation fidelity

The chance that every junction of an assembly joins the two ends meant for it, counted as the
product over the junctions of correct ligations over all ligations. It is a measurement of one
ligase at one temperature, so it belongs to the enzyme it was measured with.
_Avoid_: accuracy, specificity, efficiency

### Ligase profile

A count matrix measured with the ligase alone, over every overhang pair, at one temperature for
one incubation. Because the joining is the ligase's work, it scores the overhangs of a Type IIS
enzyme nobody has measured. It is still not a measurement of that enzyme, and a report scored
against one says so. This package ships none: it reads a file the user already holds.
_Avoid_: ligase data, T4 table

### Near-duplicate

Of two overhangs in one set: differing in fewer bases than the set's distance rule allows. Two
bases is what modular cloning standards require; a measurement can show a particular pair is
safer than that, which is why the rule is a dial.
_Avoid_: similar overhang, close pair

### Scarless junction

A junction whose overhang is the bases the part already spells there, so the product reads as
the parts do and gains nothing. Moving one changes where the cut falls and not what the product
spells.
_Avoid_: seamless, native junction

### In-frame junction

A scarless junction inside a coding sequence that has to read through it: it sits on a codon
boundary and moves only by whole codons, so the protein either side survives.
_Avoid_: fusion junction, translational junction

### Free enzyme

An enzyme with no recognition site in any of the parts, so an assembly can use it without
changing a base. A design prefers one; only when none is free does it propose domestication.
_Avoid_: clean enzyme, available enzyme

### Domestication candidate

A site an enzyme reads inside a coding sequence, which a synonymous codon change could take
away. A site outside a coding sequence is not one: removing it changes what the part spells,
which is a decision for a human.
_Avoid_: removable site, fixable site

### Part

One piece going into an assembly: what it contributes to the product, its ends arranged so that
the enzyme leaves the overhangs the design chose. How it is made is a separate question — a span
amplified from a template with a tail at each end, or one synthesised block ordered whole with a
stuffer at each end. The vector is a part like any other, and so is the product of an earlier
round.
_Avoid_: module, component, element

### Linearised vector

A circular vector opened by PCR with its primers facing outward from the span the assembly
replaces, so the whole backbone amplifies. The plasmid that templated it is taken away with
DpnI, which cuts only the methylated sites a Dam-positive host leaves, so the template cannot
transform as itself and grow as background.
_Avoid_: cut vector, destination, opened plasmid

### Ligation

Joining two fragment ends whose overhangs match. Golden Gate cuts and ligates in one tube, so a
join that puts the recognition site back is cut again and only the intended ones last.
_Avoid_: joining, annealing, sealing

### Product

The circular plasmid an assembly makes: every part's features carried to their new coordinates,
its primers annotated where they anneal, and each junction marked. It keeps the vector's origin,
so the vector's own coordinates still read true and no junction sits at base zero.
_Avoid_: construct, output, final plasmid

### Assembly plan

The whole design for one assembly: the enzyme, the overhang set, the parts and their primers, the
simulated product, the bench quantities, and the colony PCR and sequencing that confirm it. It
holds every number a protocol prints, and writes four files: the product, the primer order sheet,
the protocol as data, and the page rendered from that data.
_Avoid_: design, run, recipe

### Primer order sheet

Every oligo one design asks for, as a table to order from: the name each is ordered under, the
sequence 5' to 3', its length, and the Tm of the part that anneals.
_Avoid_: oligo list, primer table

### Phenotype

What a clone is expected to show rather than what it holds: the colour it grows on an indicator
plate, the antibiotic it resists, and whether the insert is transcribed and translated. Read off
the product's own features, so a protocol states it instead of a person asserting it.
_Avoid_: behaviour, readout, construct

### Blue/white screening

Reading a plate by colour: an insertion interrupting the vector's lacZ fragment leaves a colony
white where an empty vector is blue. It reports nothing unless the host supplies the rest of
lacZ, so the strain is part of the screen.
_Avoid_: colour screen, X-gal screen

### Fragment count

How many parts one reaction joins, the vector counted among them. NEB's tables tier on it — the
reaction volume, the enzyme dose and the cycling all step up as it rises — and the ligation
fidelity of a set falls with it, every junction added being another chance to mis-ligate.
_Avoid_: number of parts, complexity

### Insert order

The order the inserts are given in, which is the order they go round the product. Each insert's
junction takes the bases that insert begins with, so the order chooses the overhangs and the
overhangs hold the order in the tube.
_Avoid_: part order, arrangement

### Junction primer

A colony PCR primer annealing inside one insert rather than in the vector either side. It gives
its own junction a band, and it is what separates a reversed insert from a correct clone, which
two flanking vector primers cannot do.
_Avoid_: internal primer, screening primer

### Part list

The members of one position of a scheme: the protein or coding sequences that may fill it, each
named so that it says which position it belongs to. One round joins one part list to the library,
so every member of a list carries the same entry overhangs and differs only in what it codes for
and in its barcode.
_Avoid_: pool, insert list, position list

### Scheme

What a library design is given rather than works out: the positions, the enzymes that cut
internally, externally and bluntly, the stuffers, the cloning scar and the barcode length. It is
one object the user supplies, checked as it is read, and any one scheme is an instance of the
pattern rather than the only one. A part list is what fills one of its positions.
_Avoid_: config, standard, design, layout

### Entry overhang

The overhang that admits a part at one position, so that a part enters only where the scheme meant
it to. Each position has its own: a part's 5' external stuffer begins with its own entry overhang
and its internal stuffer begins with the next position's, which is how a part carries its place.
Every member of a part list shares them, which is what lets one round take a whole list.
_Avoid_: fusion site, position tag, adapter

### Internal stuffer

The piece a part carries where the next part will go, which the internal enzyme excises to open
it. It holds that enzyme's two sites facing inward and a blunt enzyme's site in its core, so the
excised piece is cut again and cannot ligate back. Its first bases are the next position's entry
overhang.
_Avoid_: filler, spacer, placeholder, dummy insert

### External stuffer

The piece at each end of a synthesised part, outside what the part contributes to the product. The
external enzyme cuts inside it to release the part as a digest fragment, and a blunt enzyme cuts
further out so that what is left of the block cannot ligate back. The 5' one begins with the
part's own entry overhang.
_Avoid_: adapter, arm, flank, tail

### Barcode

A short stretch of DNA naming one part, so that sequencing a product says which member of each
part list it carries. Its length is the scheme's, chosen with the cloning scar so that the two
together are a whole number of codons, and the barcodes of one part list stand far enough apart,
by the set's distance metric, that no two read as one. It lies in the product's reading frame, so
it spells no stop.
_Avoid_: index, tag, UMI, identifier

### Distance metric

How far apart two barcodes of one part list stand, counted one of two ways. **Hamming** counts the
positions they differ in, and sees no insertion or deletion at any distance. **Sequence-Levenshtein**
counts substitutions, insertions and deletions, charging nothing for the bases a read gains or
loses past the barcode's end, which is what a lost base really does to a read. Every set at a
Sequence-Levenshtein distance is at that Hamming distance as well.
_Avoid_: edit distance, Levenshtein distance (the plain one, which is neither), similarity

### Cloning scar

The bases a ligation leaves at a junction that neither part spells there — the opposite of a
scarless junction, which gains nothing. A scheme has one, shared at every part's 3' end, and it is
what joins each round's barcode to the barcodes already there.
_Avoid_: linker, junction sequence, spacer

### Barcode block

The run of barcodes a finished product carries, each separated from the last by the cloning scar.
Every round inserts its barcode upstream of the ones already there, so the block reads in the
reverse of the order the rounds added them. One primer pair reads the whole block, which is what
links a product back to its parts.
_Avoid_: barcode region, index block, tag array

### Round

One stage of an iterative assembly: the library built so far is cut internally to open it, that
round's parts are cut externally to release them, the two ligate, and the product is transformed,
grown and prepped to become the next round's destination. One round appends one part list and its
barcode to every member of the library at once. The product keeps the internal enzyme's sites,
which is what lets the next round open it.
_Avoid_: cycle, iteration, step

### Synthesis order sheet

Every part one library design asks for, as a table to order synthesis from: the name each is
ordered under, the whole synthesised block 5' to 3', its length, the position it fills and its
barcode. The barcode stands on the same row, so the sheet ordered from is also what decodes the
sequencing afterwards.
_Avoid_: gene list, construct table, primer order sheet (the oligo one, and its own entry)

### Library coverage

How many times over a round's colonies hold every distinct product the round could make: the
colonies counted against the number of those products. A round short of the coverage asked for
loses members no later round can put back, so it is counted for each round and not once at the
end.
_Avoid_: complexity, depth, diversity, representation

### Map

A drawing of a sequence record, whole or one region of it, as a circle or a line: features as
arrows along their strand, cut sites and primers labelled outside. A region is always drawn as a
line, numbered as the record is. No label overlaps another; one that cannot be placed is hidden,
and the map says how many.
_Avoid_: plasmid map (a linear record has one too), figure, plot

### Sequence view

A drawing of a sequence record base by base, in rows: both strands, the translation of every CDS,
features as bars, primers as arrows and enzyme names above their cut. It is drawn only beside a
map, which shows where each row lies.
_Avoid_: sequence panel, text view

### Unique cutter

An enzyme with one cut site in a sequence record, counted over the whole record even when a map
draws only a region of it, since that one cut is what opens the record for cloning. A map shows
the shipped unique cutters unless told which enzymes to show.
_Avoid_: single cutter, unique site

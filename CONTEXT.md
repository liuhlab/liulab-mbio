# Context

## Glossary

One glossary for both import packages. An entry is `liulab_mbio`'s unless it ends with a
`_Package_` line naming another.

### mbio

`liulab_mbio`, the general import package: everything a different method could use unchanged —
the sequence model, the file formats, enzymes and sites, codons, barcodes, overhangs, maps,
primers, the protocol model and the bench. A module is here when a different method could use it
unchanged. It imports **synbio** nowhere.
_Avoid_: core, common, the base package

### synbio

`liulab_synbio`, the import package holding this lab's own named methods, one subpackage each:
**iGGA** and **DMX**. A module is here when it encodes one method's choices. It imports
**mbio**, and one distribution ships both at one version.
_Avoid_: the library package, extension, plugin

### Pipeline

One way in that plans a whole experiment from the records it is given: it picks the enzymes,
designs the DNA, simulates the product and writes what the bench follows — one **protocol** for
each of the four cloning methods, and a whole **project** of them for a library. Five ship:
`plan_assembly`, `plan_gibson`, `plan_restriction`, `plan_gateway` and `plan_igga`, each the
single entry point of its own subpackage. A module below one decides a detail of the design.
_Avoid_: workflow, driver, orchestrator

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

### Reading order

The order a feature lists its segments in: the order the top strand reads them, whatever the
feature's strand, so one on the reverse strand reads them last to first. On a circular record, a
segment after the origin starts past the record's length, so the starts ascend and a sort keeps
the order. A feature cut apart at a linear record's two ends lists the higher start first.
_Avoid_: segment order, sorted order

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
steps, expected results and references. Rendered as one self-contained HTML page. It declares
what it consumes and what it produces, and nothing else about the run it belongs to: everything
else it needs it holds itself.
_Avoid_: SOP, recipe, method

### Item

One thing passed from a protocol to the one after it: a plasmid, a plate of colonies, a pooled
library. The name is the contract a **project** chains by; what it is, what it has to meet and
where it waits are prose for the bench, which nothing parses.
_Avoid_: input, output, artefact, deliverable

### Project

Protocols run in order, each handed what the ones before it produced. It holds what no single
protocol owns: the background a reader is told before the first protocol, the bill for the whole
run, the **source** its rows cite, and the checks that judge the design rather than one bench
procedure. What a build chooses
is a file it is read from, and `docs/adr/0010-method-in-code.md` draws that line; what a build
writes is this chain of protocols.
_Avoid_: workflow, pipeline, campaign

### Build

What one run of a **pipeline** chooses, as against what its method fixes: for a library, the
positions and their part lists, the vector, the host, the oligo length, the batch size, the
completeness, the seed, the barcode length and distance, and any further enzyme to keep clear.
It is read from one `project.json` and checked where it is read. It overrides nothing the method
states: where both have a say, the two compose. What a build writes is a **project**, which is
the chain of protocols and a different thing.
_Avoid_: project (the chain of protocols, and its own entry), configuration, run
_Package_: liulab_synbio

### Step

One numbered unit of a protocol: what to do, what it needs, what a successful result looks
like, and what to do when the result is wrong.
_Avoid_: task, procedure

### Section

Where in a protocol a step belongs, such as `Day 1` or `Round 2`. A label, not a container: the
steps stay one list and the numbering runs through them. Both navigation lists group by it, one
collapsible group each, counted from the bench's own marks.
_Avoid_: stage, phase, part

### Key

What a step, a protocol or a **project** is addressed by, assigned by its builder and never
derived from its wording, so rewording a title moves no tick the bench has made. A page's
anchors and its marks are built from it, and a key holds no dot, which is what keeps one step's
marks out of another's.
_Avoid_: id, name, slug (what a key is spelled as)

### Mint

To give a protocol or a **project** the **key** it was built without, from a digest of its own
content, so that one run's protocols are addressed apart from every other run's. Two protocols
of one run that mint the same key are a builder defect, and the run is refused.
_Avoid_: generate, hash, assign

### Wait

Time a step spends waiting on someone else, which nobody attends: a vendor's turnaround on an
oligo pool, a plate sent away to be sequenced. It is neither a thermocycler stage nor a
countdown someone starts, and most of a real project's calendar is this. How long is written as
whoever states it does, so `10-15 working days` stands as it is, and an empty one is an admitted
unknown. What a step holds beside it is its hands-on time: how much of it someone stands over,
unknown until a source states it and never written as a zero.
_Avoid_: delay, downtime, lead time, incubation

### Reaction table

The volumes of one reaction, scaled to a master mix for several reactions plus an overage. A
component marked per-tube, such as template, is added to each tube instead of to the mix.
_Avoid_: mix table, recipe

### Thermocycler program

Stages run in order, each a list of incubations repeated for a number of cycles. An incubation
whose temperature steps a set amount each cycle makes its stage a **touchdown**.
_Avoid_: cycling conditions, PCR conditions

### Figure

What a step shows rather than describes: the record to draw, the stretch of it, which view, and
what the map is pointed at. It is a spec and never a drawing, so it cannot disagree with the
design it names, and the page lays it out when it renders. A record is named by a path, relative
to the directory the protocol was read from, since a protocol holds no bases of its own. A map is
what a figure becomes, not another word for one.
_Avoid_: image, illustration, panel, picture

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

### Compatible ends

Two cut ends a ligase can join: both blunt, or overhangs spelling the same bases on the same
strand. The bases alone do not decide it, because a 5' overhang and a 3' overhang spelling the
same bases run the wrong way for each other. So an end is its overhang and its end type
together, and the end type comes from the enzyme that made the cut.
_Avoid_: matching ends, compatible overhangs

### Dephosphorylation

Taking the 5' phosphates off a cut vector so it cannot close on itself. A ligase seals a join
only where one of the two ends carries a phosphate, so a backbone whose own two ends are
compatible religates empty unless its phosphates are taken away; the insert keeps its own, and
those are what the ligase seals. It is needed only where those two ends can meet.
_Avoid_: phosphatase treatment, vector prep, CIP

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

### Double digest

Cutting one record with two enzymes in the same tube rather than in two rounds, which saves an
afternoon and a clean-up between them. It holds only where both enzymes keep enough activity in
one buffer and want the same temperature. This package judges the buffer from the one each
supplier sells that product in: the same buffer for both is an answer, and anything else carries
no verdict rather than a pass, because how much activity an enzyme keeps in another's buffer is
a measurement nothing shippable states.
_Avoid_: dual digest, two-enzyme digest

### Diagnostic digest

A digest run on a miniprep to say whether it carries the insert, read as band sizes on a gel
rather than as sequence. It is the check before sequencing rather than a substitute for it, and
it is worth running only where a correct clone and the plasmid it could have come from give
bands a gel can tell apart.
_Avoid_: test digest, check digest, analytical digest

### Domestication

Taking a recognition site out of a coding sequence by changing one codon for another spelling
the same amino acid, so the reading frame and the protein survive. The replacement is the one
the host uses most often among those that create no new site. A site outside a coding sequence
cannot be domesticated, and is reported for someone to decide about. What counts as a coding
sequence is a feature the record types `CDS`, or one whose own bases spell a whole protein —
`ATG`, whole codons, one closing stop — because most files type an open reading frame loosely.
`docs/adr/0013-coding-by-bases.md` holds that rule.
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
general. Golden Gate, Gibson assembly, classical restriction and ligation and Gateway are the
ones this package plans, the last of them joining nothing — its att sites recombine; a library
built in rounds is a pipeline over a method, not a method of its own.
_Avoid_: cloning strategy, technique, approach

### Molar ratio

How many molecules of insert an assembly gets for each molecule of vector. Counted in moles and
not in mass, because a short insert weighs less at the same ratio: a protocol prints pmol beside
ng for that reason.
_Avoid_: insert ratio, stoichiometry

### Junction

Where two fragments meet in an assembled product: in Golden Gate the four bases, three for SapI,
that one overhang paired with its match; in Gibson assembly the overlap the two share; in
restriction and ligation the recognition site the two cut ends came from, which the join puts
back; in Gateway the whole att site the recombination wrote, whose 25 bases straddle the
boundary, so where the insert starts is not where the junction does. Validation reads across it.
_Avoid_: joint, seam, fusion site

### Overlap

The bases two fragments both spell at a Gibson junction. One of the two already spells them and
the other carries them as a primer tail, so the junction adds nothing. Its length is set per
assembly product and per fragment count, never per fragment length, and it is lengthened inside
that band until its melting temperature reaches the product's floor.
_Avoid_: homology arm, overhang, homology region

### Assembly product

The kit a Gibson reaction is run with, such as NEBuilder HiFi. It is the row of data that
carries the overlap band, the reaction, the incubation, the molar ratio and the fragment-count
limit, each with the citation it came from. Not the **product**, which is the plasmid.
_Avoid_: master mix, kit, enzyme mix

### Oligo stitching

Making a short part out of overlapping oligos that tile both its strands, assembled in the same
reaction rather than amplified or ordered from a synthesis vendor. It suits a linker, a tag or a
short promoter. The oligo count, the oligo length and the overlap between neighbours are the
research note's, and a part too long for the oligos the note allows is refused rather than
stitched out of more.
_Avoid_: oligo annealing, oligo assembly, gene synthesis

### Bridging oligo

One oligo joining two fragments that share no homology at all, carrying homology to the end of
each, so neither needs a tailed primer. An amplicon made for something else goes in as it is.
It primes nothing, so no primer threshold judges it and its row on the order sheet carries no
verdict.
_Avoid_: splint, stitching oligonucleotide (a vendor's word for this, not for **oligo
stitching**)

### End soak

The 60 °C step every Golden Gate program ends with. It is a digest and not heat inactivation: it
cuts destination plasmid that never opened or has closed again, so fewer empty colonies grow.
_Avoid_: final incubation, kill step

### Colony PCR

A PCR templated by a colony lifted off the plate with a toothpick, which the long hot step at
the start lyses. It says which clone a colony carries without a plasmid prep.
_Avoid_: screening PCR, colony screen

### Reversed insert

A clone carrying its insert the other way round. Only a join whose two ends pair either way
round makes one — one enzyme cutting both, two leaving compatible ends, or blunt ends — and
the method that joins them is what builds the plasmid to measure. Two vector primers flanking
the insert give it the same band as the correct clone, so a gel separates the two only with a
third primer inside the insert.
_Avoid_: flipped clone, wrong orientation

### att site

One of the eight sites Gateway recombines: attB, attP, attL or attR, each numbered 1 or 2.
Twenty-five bases recombine, and the 7 bp overlap inside them is what decides which site pairs
with which. A reaction rewrites the pair it consumes — attB with attP gives attL and attR, attL
with attR gives back attB and attP — so the sites are the mechanism and not a scar left by it.
One is found in a record by searching for it, never by matching a shipped sequence, because a
vendor's own sites drift between products.
_Avoid_: recombination site, att sequence, recombination region

### Entry clone

The plasmid carrying an insert between attL1 and attL2: what a BP reaction makes, and what an LR
reaction spends. It is the record a Gateway job reuses, because one entry clone feeds every
destination vector after it. A plan given one plans no BP reaction; a plan that makes one writes
it as a file of its own.
_Avoid_: donor clone, middle vector, pENTR

### Destination vector

The plasmid an LR reaction moves an insert into: attR1 and attR2 around a ccdB cassette, with
the promoter, tag and marker the finished clone is wanted for. It is the user's own file. This
package ships no vector catalogue, and reads a vector's sites out of the bases it is handed.
_Avoid_: expression vector, target vector, pDEST

### Donor vector

The plasmid a BP reaction moves an insert into: attP1 and attP2 around a ccdB cassette, and the
backbone the entry clone then carries. It is the user's own file, read for its sites the way a
destination vector is. A plan takes one only where it was handed an insert rather than an entry
clone, because a plan given an entry clone runs no BP reaction.
_Avoid_: entry vector, BP vector, pDONR

### BP reaction

The reaction making an entry clone: an insert with an att site on each end recombines with a
donor vector, the plasmid holding attP1 and attP2 around a ccdB cassette. The insert gains attL
sites and the cassette leaves in the by-product. Its entry clone is grown up and purified before
an LR reaction takes it, never chained straight on.
_Avoid_: BP cloning, entry reaction, donor reaction

### LR reaction

The reaction making the expression clone: an entry clone recombines with a destination vector,
the insert gains attB sites, and that vector's ccdB cassette leaves in the by-product. Every
Gateway plan runs one, and the BP reaction before it only where no entry clone was handed in.
_Avoid_: LR cloning, expression reaction, destination reaction

### ccdB counter-selection

What keeps a Gateway plate clean. A donor or destination vector carries ccdB, which kills an
ordinary strain, so only a clone that has traded the cassette away grows. It decides two strains
rather than one: the vectors are grown in a strain resistant to CcdB, and the clone is selected
in a strain CcdB kills. A strain carrying F′ cancels the whole thing, because the episome's
ccdA neutralises CcdB.
_Avoid_: negative selection, ccdB selection, suicide gene

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

### Reserved overhang

An overhang a step outside a selection already spends, so the set being designed leaves it
alone: the cargo junction a finished library goes on to, for one. It binds a candidate the way
a junction already taken does — neither that overhang nor a near-duplicate of it — but it is
not a junction of this set, so taking one is refused as reserved rather than as a repeat.
_Avoid_: excluded overhang, forbidden overhang, blacklist

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

### Site out of reach

A site carried more than once in a record, with so much identical sequence on both sides that no
oligo placed there is unique to one copy. It is counted before anything else a domestication
route decides, because it is not a bench job at all: the enzyme changes, or the plasmid is
bought whole.
_Avoid_: repeat site, untouchable site

### Cargo enzyme

The enzyme that admits cargo to a working vector. Every other enzyme a method names is fixed by
DNA already on the shelf; this one is searched for per vector, the vector being the user's own.
_Avoid_: insertion enzyme, final enzyme

### Cassette

The piece a vector gives up to admit a part, bounded by one enzyme's two cuts: the entry overhang
at one end and the cloning scar at the other. A round's destination gives up the method's
internal stuffer; a working vector gives up a ccdB cassette, which kills anything that keeps it.
What a cassette bars outside itself is whatever shares that vector's tube, and the two vectors do
not share a list.
_Avoid_: insert site, cloning site, landing pad

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
so the vector's own coordinates still read true and no junction sits at base zero. A Gibson plan
holds it as `Plan.plasmid`, because `Plan.product` there is the **assembly product**.
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

### DAD-GGA-DMX

The experiment that makes cargo and reads it back, one design at a time. A gene is split into
fragments, the fragments come out of an oligo pool by nested PCR and assemble into the DMX
vector, and one design is archived per well. Reading it back is **DMX**: it takes that archive
to clonal wells, marks each well by the DMX barcode kit or by index PCR, sequences it, and calls
a pass per well, so identity is read per member and stays with well position. It is the only
place a member is picked or read on its own. Whether a design is read back at all is a **build**
choice, and the cargo an **iGGA** round consumes is made here.
_Avoid_: cargo pipeline, the validation pipeline
_Package_: liulab_synbio

### DMX

A protocol of its own for multiplexed validation: a design sits one per well, the well is
marked, sequenced and called on its own, and identity stays with well position throughout. It
takes any cargo, and what **iGGA** builds is one kind of cargo among others, so DMX stands
beside iGGA rather than downstream of it — a caller chains it. Two routes mark a well, barcode
ligation or index PCR, and one judgement reads them; the picking, the pass rule and the reformat
are shared, while the marking step, the plate and the depth floor are the route's own, because
each floor was measured on its own library prep.
_Avoid_: the validation pipeline, read-back pipeline, QC
_Package_: liulab_synbio

### Validation floor

The **fragment count** at or above which a design is read back one well at a time. A **build**
states it or leaves it out: left out, nothing is read and the cargo stays polyclonal, and zero
reads every design. No floor ships, because the measured curve gives a design's chance of a clean
colony and not the chance worth paying to check. A build that states one also names which of
the two marking routes reads its wells: **barcode ligation** or **index PCR**.
_Avoid_: validation threshold, QC cutoff, validation level, Route A, Route B
_Package_: liulab_synbio

### iGGA

The pooled experiment built from the parts DAD-GGA-DMX synthesised: **rounds** of Golden Gate,
each appending one part list and its barcode to every member of the library at once. Its product
is a pool that goes straight to a pooled screen or to a downstream experiment. No member of it is
ever picked, and nothing in it is re-validated for identity. What judges it are representation and
linkage reads over the whole pool — which combinations are there and how evenly, and which barcode
combination goes with which cargo — never a read of one member. A bound taken from a gate whose
purpose was pickable clones does not transfer to it.
_Avoid_: iterative assembly, library build, the library pipeline
_Package_: liulab_synbio

### Part list

The members of one position of a scheme: the protein or coding sequences that may fill it, each
named so that it says which position it belongs to. One round joins one part list to the library,
so every member of a list carries the same entry overhangs and differs only in what it codes for
and in its barcode.
_Avoid_: pool, insert list, position list

### Scheme

What a method fixes rather than works out: the positions, the enzymes that cut internally,
externally and bluntly, the stuffers, the cloning scar and the barcode length. It is one object
built in code, whose invariants are checked when its module is imported, and not a file anyone
supplies — `docs/adr/0010-method-in-code.md` reversed that, so a second scheme is a second
method rather than a second file. A **part list** fills one of its positions, and what one
**build** chooses sits beside it without overriding any of it.
_Avoid_: config, standard, design, layout
_Package_: liulab_synbio

### Entry overhang

The overhang that admits a part at one position, so that a part enters only where the scheme meant
it to. Each position has its own: a part's 5' external stuffer begins with its own entry overhang
and its internal stuffer begins with the next position's, which is how a part carries its place.
Every member of a part list shares them, which is what lets one round take a whole list.
_Avoid_: fusion site, position tag, adapter
_Package_: liulab_synbio

### Internal stuffer

The piece a part carries where the next part will go, which the internal enzyme excises to open
it. It holds that enzyme's two sites facing inward and a blunt enzyme's site in its core, so the
excised piece is cut again and cannot ligate back. Its first bases are the next position's entry
overhang.
_Avoid_: filler, spacer, placeholder, dummy insert
_Package_: liulab_synbio

### External stuffer

The piece at each end of a synthesised part, outside what the part contributes to the product. The
external enzyme cuts inside it to release the part as a digest fragment, and a blunt enzyme cuts
further out so that what is left of the block cannot ligate back. The 5' one ends with the
part's own entry overhang.
_Avoid_: adapter, arm, flank, tail
_Package_: liulab_synbio

### Block vector

The vector one position's synthesised cargo closes into, which is then the donor a round
releases that part from. It supplies the external stuffers the cargo is not synthesised with,
so it has to be opened on the overhang that position's parts enter on: one backbone offers one
pair, and a build over several positions needs one block vector each. They are the build's own
destination vector, differing only in the bases of that overhang.
_Avoid_: block backbone, donor plasmid, part vector
_Package_: liulab_synbio

### Barcode

A short stretch of DNA naming one part, so that sequencing a product says which member of each
part list it carries. Its length is the scheme's, chosen with the cloning scar so that the two
together are a whole number of codons, and the barcodes of one part list stand far enough apart,
by the set's distance metric, that no two read as one. It lies in the product's reading frame, so
it spells no stop.
_Avoid_: index, tag, UMI, identifier

### Barcode kit

The 96 plasmids **DMX** marks wells with, in four groups of 24, used as supplied. One member is
a plasmid carrying a group, an index and a UMI, and a well's marks are arithmetic from its
address rather than a recorded draw. The sequences are not
shipped: they are read from a copy the user holds. A member of this kit is not a **barcode**,
which names one part of a library.
_Avoid_: index set, tag kit, barcode plate
_Package_: liulab_synbio

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
which is what lets the next round open it. The first round opens the destination vector, which is
the library before any part list has been appended, so a build over _n_ positions runs _n_ rounds.
A step that seats one part in a carrier makes no library and is not a round.
_Avoid_: cycle, iteration, step
_Package_: liulab_synbio

### Synthesis order sheet

Every part one library design asks for, as a table to order synthesis from: the name each is
ordered under, the whole synthesised block 5' to 3', its length, the position it fills and its
barcode. The barcode stands on the same row, so the sheet ordered from is also what decodes the
sequencing afterwards.
_Avoid_: gene list, construct table, primer order sheet (the oligo one, and its own entry)
_Package_: liulab_synbio

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
features as bars, primers as arrows and enzyme names above their cut. It is drawn only with a
map, which shows where each row lies.
_Avoid_: sequence panel, text view

### Highlight

What a map is pointed at: the names it lights. A lit item keeps its colours, and every other
item, its label with it, paints one pale grey. Nothing moves, so each label keeps where it sits,
and a lit label is the last to hide where labels crowd.
_Avoid_: callout, emphasis, focus, selection

### Unique cutter

An enzyme with one cut site in a sequence record, counted over the whole record even when a map
draws only a region of it, since that one cut is what opens the record for cloning. A map shows
the shipped unique cutters unless told which enzymes to show.
_Avoid_: single cutter, unique site

### Reaction

One tube: the molecules in it, each bound to the role it plays there, and the enzymes acting on
them. A reaction is a digest, a ligation, or a one-pot that does both, and it carries how many
vessels differ only in which species fills a role -- one reaction over a 96-well array is one
reaction, not 96. One reaction's product is the next one's input; a sequence of reactions records
that order and no constraint of its own.
_Avoid_: step, tube, assembly, mix

### Role

What a molecule is in a reaction for: the destination opened, the donor cut out of its block, the
insert going in, or the carrier holding a part until it is wanted. A role binds to one species or
to a pool of them, and never to two pools in one tube.
_Avoid_: part type, slot, component

### Gate

What judges a finished design: the predicates `liulab_mbio` holds, called with one method's
parameters, over the molecules of each reaction. It reads reactions rather than bare records, so
it is blind to how a design was reached and still says which tube a failure belongs to. A design
an agent composed and one a pipeline wrote are judged the same way.
_Avoid_: validator, linter, QC
_Package_: liulab_synbio

### Judgement

One check the gate made, what it judged, and the findings behind it. A finding is whatever domain
object `liulab_mbio` already returns where the gate looked -- a cut site, a span -- so a failure
can be drawn on the record it occurred in. A failing judgement names what is wrong and need not
name a remedy.
_Avoid_: violation, error, issue
_Package_: liulab_synbio

### Vessel

Something the bench holds material in whose contents have no positions: a tube, a flask, a
reservoir, a 25 cm BioAssay plate, a flow cell. A vessel whose positions form an array is a
**plate**.
_Avoid_: container, tube (unless it is one)

### Plate

A vessel whose positions form an array, named by its well count: 12, 24, 96, 384 or 1536.
Format is one parameter, not a kind. A plate carries a **seating** and says nothing about what
differs between the reactions it holds.
_Avoid_: microplate, microtitre plate, array

### Stock plate

A plate a set of oligos is resuspended in once and kept at that concentration: the copy nothing
runs from. Every **working plate** is split from it, so it is thawed only to make one.
_Avoid_: master plate, mother plate, source plate

### Working plate

A copy of a **stock plate** at the concentration a reaction takes, seated the same way, so a
well's address names the same oligo on both. One is thawed for one run and thrown away, never
returned to the freezer; that is a **rule** on the plate as a material, not a field a protocol
declares. A run splits as many as it has runs ahead of it.
_Avoid_: daughter plate, aliquot plate, dilution plate

### Seating

Where each thing sits in a plate: a map from a well to the name of what it holds. A name has to
resolve against the protocol's own materials, oligos, vessels and plates, so a dangling one is
reported. A seating never restates what a reaction already says differs between its vessels.
_Avoid_: layout, plate map, well assignment

### Label

What a well is described as, as against what sits in it: a grouping such as `quarter 3`, a
printed address, or the name of something the protocol does not declare. A label resolves
against nothing and is drawn as it stands, so it is never a dangling reference. A well a step
has to name carries a **seating**; a well that only has something to say about itself carries
this.
_Avoid_: annotation, caption, tag

### Transfer

Material moved from one well to another, with a volume and the instrument that moves it. One
shape covers every move: one source to one destination is a plain or an acoustic transfer, many
sources to one destination is a **pool**, and a dense re-layout that leaves out the wells that
failed is a **compaction**. A protocol step says where a thing is by holding the transfer rather
than describing it.

A transfer whose moves take every destination well from one source plate at one stride, at one
volume, is a **stamp**: a stride and the source well it starts at stand for the whole move list,
and a page draws the two plates rather than printing a row each.
_Avoid_: dispense, reformat

### Source

A document a number was read from: what it is, its edition, where it was read, how, when, and
the research note under `docs/research/` it was read into, where it has one. A protocol names
each one once, and a run names the ones its own shared pages cite, such as the record pricing
its bill; a **citation** of a source and a locator hangs on the row that carries the number.
_Avoid_: reference (the protocol's own bibliography entry), provenance

### Store

Where a fact is kept. There are five: a file the user holds, `src/liulab_mbio/data/`, a
catalogue-keyed table on a **material**, the protocol JSON, and a sourced constant in the
method's own bench module. A new fact goes in the first that fits, asked in that order, so two
facts of one kind cannot end up apart.
`docs/adr/0020-a-fact-goes-in-the-first-store-that-fits.md` holds the questions.
_Avoid_: location, home, storage

### Material

Something a protocol consumes, named as the bench names it, with its supplier and catalogue
number where it has one. It carries what it brings into the tube, its own parameters, the
**rules** that follow it into every step that uses it, and its cautions — what would hurt the
person or the material, written as the action to take. A fact keyed by a catalogue number hangs
on the material and never on a step that mentions it, so it cannot be edited out of one step and
left in another. How much of it a step takes is the step's, since the same material is pipetted
at different amounts.
_Avoid_: reagent, consumable

### Rule

A prohibition or a requirement attached to a **material**, which follows it into every step that
uses it: `forbids` keeps something out of the tube, `requires` keeps it in, and each may be
conditional on what else the tube holds. A rule that computes a number is not one of these; that
stays a function. Nor is a caution: that is a separate thing a material carries beside its
rules, and it names the action that keeps someone or something safe.
_Avoid_: constraint, warning

### Hole

A number nobody sourced, standing where the number would be. The field it belongs to stays
empty and the hole stands beside it, so a hole is never read as a value and never judged. It
says what is missing, why and what would fill it, and names the ticket it is routed to where one
owns it. A hole naming a gap no source closes is what a finished plan keeps; a hole waiting on a
source nobody has read fails the plan. A guess is a defect.
_Avoid_: missing value, TODO, placeholder

### Price record

A banded tariff the user holds and the package never ships: rows of a key, its **bands**, a
charge and a **basis**, in one currency. The key is a catalogue number where the item has one,
so a price and a shipped parameter meet at the same key.
_Avoid_: price list, tariff, cost table

### Band

A quantity a price row is bounded by, inclusive at both ends, such as 101 to 200 oligos. A row
may carry more than one, since a vendor may price a pool by count and by length together, and
writes them in one cell: `count 101-200; length_nt 1-200`. A tier with no top leaves its high
end empty, so nobody writes a sentinel. A band bounds a quantity and is not a span in a
sequence, so the coordinate rule does not reach it.
_Avoid_: tier, bracket, range

### Basis

What a price row's charge is for: `per order` is flat inside the band, and `per unit` is
amortised over a pack. It is a property of the row and not a policy, which is what lets one
shape price a pool and an enzyme.
_Avoid_: pricing model, unit

### Bill

What a run consumes, and what it costs where a price record prices it. A quantity comes from the
design and is always there; money comes only from a record, and a row nothing prices carries a
hole instead. No figure is ever estimated, and price steers no part of a design.
_Avoid_: quote, invoice, budget

### Headroom

How far a quantity sits from the edge of the band pricing it, reported beside the row, such as
153 oligos with 47 below the next band. It is the whole of what a price cliff gets, because the reader
decides and the design does not.
_Avoid_: slack, margin, buffer

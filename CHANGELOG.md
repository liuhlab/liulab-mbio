# Changelog

Every change worth knowing about, newest first. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Version numbers are CalVer tags
of the form `vYYYY.M.PATCH`, and the tag is where the version comes from — nothing here
sets one.

## [Unreleased]

### Added

- Gateway cloning, end to end. `liulab_mbio.cloning.gateway.plan_gateway` and `liulab_mbio
  cloning gateway plan` take an insert and a destination vector and plan both reactions. Nothing
  is cut and nothing is ligated here: two att sites recombine, and the reaction rewrites the
  sites themselves. So there is no enzyme to pick and no overhang to score. There are three
  routes, and your files choose one. An entry clone you already hold plans the LR reaction
  alone. A fragment that already carries att ends plans the BP reaction first, and its entry
  clone feeds LR. `--amplify` designs the two attB primers for a plain gene and amplifies it
  onto those ends. `Plan.write` writes the four files every cloning plan writes, plus the entry
  clone as a fifth where the BP reaction ran.
- Both vectors are your own files. No vector list ships with this package, and none is planned.
  A vector is recognised by searching its bases for att sites, allowing one base to differ
  outside the core of a site and none inside the seven that decide which site it is. The
  supplier's own sites drift between products, so a stored copy would turn a working vector
  away. A record carrying no site is refused, and the message names every site it looked for.
- The junction is not clean, and the plan says what it costs. A whole att site lands at each end
  of the insert, so the clone gains 25 bases there. Where a tag is read through one of them,
  `--fusion` says which end, and the frame is judged there: two more bases at the front end, one
  at the back. The tail an attB primer carries is four G residues, the 25 bp site and those
  frame bases. Those G residues leave with the BP by-product, so the entry clone reads the same
  whichever route made it.
- The protocol names both strains, and it screens what it made. A donor and a destination vector
  carry ccdB, so each is grown in a strain that resists it, while the clone is selected in a
  strain it kills. A strain carrying F′ fails that check by name, because its `ccdA` cancels the
  selection. Whether a reaction plated on a resistant strain still counter-selects is stated
  nowhere, so that check carries no verdict rather than a pass. Every plan then screens colonies
  across both joins and designs sequencing primers reading in from outside each. There are two
  lanes and not three: the two att sites differ, so the insert cannot go in backwards.
- A `gateway/METHOD.md` behind the `molecular-cloning` skill, which now chooses between four
  methods on the same five things, with a `variants.md` beside it for the published routes this
  command does not plan. A docs page follows one Gateway job end to end.
- Classical restriction and ligation cloning, end to end.
  `liulab_mbio.cloning.restriction.plan_restriction` takes a vector and an insert, plans both
  digests and the gel that separates their fragments, works the ligation out in picomoles, and
  designs the colony PCR, the diagnostic digest and the sequencing that confirm the clone. It
  says what each junction now reads, because this method's join puts the enzyme's own site back,
  and it says whether those extra bases hold your reading frame. `Plan.write` writes the same
  four files every cloning plan writes. The second record may be the insert itself or the
  plasmid it is cut out of. One already carrying both sites is cut out and taken off a gel. One
  carrying neither is amplified first, with a spacer and the site on each primer tail.
- `liulab_mbio cloning restriction plan` on the command line, and a `restriction-ligation`
  method file behind the `molecular-cloning` skill, which now sends a job to this method instead
  of calling it unsupported.
- The enzyme pair is chosen for you when you name none. Every pair of the 18 shipped enzymes
  that cut inside their own site is weighed, and each rejected pair is kept with the one
  sentence that rejected it: no site in the vector, a site inside the insert, a cut vector that
  closes on itself, or a digest that throws away more of the vector than it keeps. What survives
  is ranked on the buffer, the temperature, the heat step and the host's own methylation first,
  then on how much of the vector it gives up and how much of the insert it carries whole. Name
  one or two enzymes and that pair is planned, or refused in the same words.
- One enzyme at both ends, and blunt ends, are planned rather than refused. The cut vector can
  then close on itself, so the plan adds a phosphatase step before the ligation. The insert can
  also go in either way round, so the colony PCR reads out of the insert itself and draws a lane
  for the reversed clone. Where the two ends cannot pair, that lane is left off: it was a lane
  for a plasmid that cannot exist.
- Gibson assembly, end to end. `liulab_mbio.cloning.gibson.plan_gibson` and `liulab_mbio cloning
  gibson plan` take a vector and up to five inserts, in the order they go round the product and
  either way round. Nothing is cut, so a part that reads a site for every Type IIS enzyme still
  goes in as it is. The vector is opened by PCR across the span the inserts replace, or handed in
  already linear. Each junction gets an overlap of the length and melting temperature the kit
  documents, carried as a primer tail by the part on the other side, so the junction spells only
  what the parts spell. `Plan.write` writes the same four files a Golden Gate plan writes. The
  protocol runs from the two PCRs through the DpnI digest, a column cleanup, a reading of each
  concentration, the assembly and the plate, and ends at the sequencing rather than at the
  transformation.
- The kit is yours to name, and it sets the numbers. `--product` takes NEBuilder HiFi, the Gibson
  Assembly Master Mix or In-Fusion, and the start of a name is enough. Each carries its own
  overlap rule, reaction, incubation, molar ratio and documented fragment count, read from that
  supplier's own manual. Nothing is borrowed from one supplier to fill a gap in another: where
  Takara publishes no melting temperature for an overlap and no picomole band, those checks say
  so instead of passing. So do an overlap's GC and how alike two overlaps are, which nobody
  quantifies for any kit.
- A short part can be built from oligos instead of amplified. `--route stitch` lays a linker, a
  tag or a short promoter out as oligos that overlap each other and tile both strands, assembled
  in the same tube. Two ceilings, and they differ: above 500 bases the plan refuses, because
  twelve 60-base oligos are the most the method allows; outside 60 to 150 bp it warns, the window
  the route is worth using in.
- Two fragments that share nothing can be joined by one oligo. `--bridge BEFORE:AFTER` designs an
  oligo carrying bases from each end, so neither fragment needs a tailed primer and an amplicon
  made for something else goes in as it is. Both kinds of oligo are rows of the one order sheet,
  with the job each does. Neither primes anything, so no threshold judges it and its row carries
  no verdict rather than a pass nothing measured. In-Fusion takes no single-stranded oligo, so it
  refuses both routes and names the kits that can.
- A `gibson/METHOD.md` behind the `molecular-cloning` skill, which now chooses between two
  methods on the same five things, and a docs page following one Gibson job end to end.
- Each shipped enzyme now names the buffer its supplier sells it in, as `Enzyme.supplied_buffer`.
  All sixteen multiple cloning site enzymes read `rCutSmart Buffer`, so any pair of them can be
  cut in one tube; BsmBI-v2 and BspQI read `NEBuffer r3.1`. The buffer belongs to the product
  rather than to the enzyme name, and each one was read from that product's own page.
- A map of a sequence record. `liulab_mbio.plot.draw_map` and `liulab_mbio plot map RECORD -o
  map.html` draw a `.dna`, GenBank or FASTA file as a circular map, written as one HTML page that
  opens offline. A feature keeps the colour its file gives it. One with no colour takes a
  default for its type that people with colour blindness can tell apart. A feature's name sits on
  its arrow when it fits there, curved along the circle and upright at the bottom. Any other name
  sits in a box outside the circle, and no two boxes overlap. Hovering over a feature shows its
  name, type, span and length.
- A map shows where each primer binds and where enzymes cut. A primer is drawn in purple and
  labelled with its span. By default the map names the shipped enzymes that cut the record once,
  in bold, each with the base it cuts after. `--enzyme` names others, and the map then shows every
  site each one cuts. `--no-features`, `--no-primers` and `--no-cut-sites` leave a layer off,
  `--hide-type` leaves out a feature type, and `--source` draws the `source` feature. `draw_map`
  takes the same choices.
- A map as a line. A linear record, every FASTA file included, is always drawn as one, and
  `--linear` opens a circular record into one, with `• • •` at each end. `--region` draws one
  stretch of a record as a line: a feature's name, or `START..END` as the map numbers bases,
  such as `2680..10` across the origin. The line keeps the record's numbering, and an enzyme is
  bold only if it cuts the whole record once. Features lie under the line, named inside or
  underneath, and the labels above rise in steps so none overlaps another. `draw_map` takes
  `linear` and `region`.
- A crowded map grows before it hides a label: pUC19 with all 99 of its unique 6+ cutters still
  shows every label. Past that size, enzyme sites hide first, then primers, then boxed feature
  names, and only where labels crowd one another. A name on its feature never hides. The map says
  what it hid at its bottom right, such as `3 enzyme sites and 1 primer are hidden`, and the page
  lists the names when you hover over it. `plot map` prints the same words after the files it
  wrote, and `Drawing.hidden` gives each hidden label.
- A map as a PNG or a PDF, which looks the same on any computer, with no font installed. Name
  the file `map.png` or `map.pdf`; `-o` takes several files in one run, and `--dpi` sets the
  PNG's resolution, 300 by default. The letters are drawn as shapes, so the text in them cannot
  be searched or copied.
- The bases beside the map, as SnapGene's Sequence tab shows them. `--sequence-view` sets the
  record out in rows of 60 bases, or as many as `--bases-per-row` says. Each row has a ruler, both
  strands with a tick at every base between them, and the number of its last base at its right.
  `--one-strand` leaves out the bottom strand. Features lie under the bases as bars, and every
  coding sequence shows its protein under its bases in three-letter codes, with each stop in red.
  A region keeps the record's numbering. The page shows this view up to 100,000 bases; past that,
  name a region. A PNG is one image with this view under the map; one too big to make stops with
  an error that asks for a region or a PDF. A PDF has the map on its first page and the rows on
  the pages after, each page as many rows as fit, so a long record prints. `draw_map` takes
  `sequence_view`, `bases_per_row` and `both_strands`.
- Switches on the page. Features, primers and cut sites each have one, and so does each feature
  type. A switch shows or hides them where they are, so no label moves. `--no-features`,
  `--hide-type`, `--source` and the others now set what the page shows first, and a PNG or PDF
  still leaves out what they switch off. A circular record drawn whole flips between a circle
  and a line at the top right, and `--linear` shows the line first. A click on a feature, primer
  or cut site highlights it. The notice of hidden labels counts only those whose switches are on.
- Primers and cut sites in the sequence view. Each primer is an arrow beside the strand it copies.
  Its 5' tail bends away from the bases, and each base it does not pair with is marked in red. Each
  cut is drawn through both strands, so you can see the overhang it leaves. The enzymes that cut
  there are named above it, one name a line, in bold if they cut the record once. The view shows
  the same enzymes, layers and feature types as the map.
- The map and the sequence view work together on the page. The page carries the sequence view of
  up to 100,000 bases behind a Sequence switch at the top right, and `--sequence-view` turns it on
  from the start. Each view scrolls on its own. A click on a feature, primer or cut site in one
  view highlights it in both, and scrolls the other view to it. Both strands shows or hides the
  bottom strand without moving anything else, and `--one-strand` hides it from the start. Hover
  over a base to see its position. Drag over bases to select them: the page shows the stretch as
  `start .. end (n bp)`, and copying, or the Copy button, puts the top strand's bases on the
  clipboard. Nothing on the page changes the record.
- A repo-local `plot-map` skill that calls `plot map` and `draw_map`, so a coding agent asked to
  show a record hands back a map.
- A Golden Gate plan picks the codons of the host the part is for: `plan_assembly` takes
  `codon_table`, `goldengate plan` takes `--codon-table`, and a plan that has to propose
  domestication names the codons it would move to. The strain the protocol transforms is still
  `--host`, and still *E. coli*.
- A barcoded library of protein combinations, built in rounds rather than one pot.
  `liulab_mbio.library.plan_library` takes lists of proteins, or coding DNA, with a scheme and a
  destination vector. It picks the overhang standard that costs the proteins fewest changed
  residues, writes each part's synthesis sequence with its own barcode, fits a vector that cannot
  be opened yet, simulates every round, and counts the colonies a round needs for the coverage
  asked for. `LibraryPlan.write` puts the synthesis order sheet, the barcode table, the
  amino-acid change table, a record for each round, the assembled product, the protocol as JSON
  and the page rendered from it in one directory.
- `liulab_mbio library plan` on the command line, and a repo-local `protein-assembly` skill that
  calls it.
- Writing DNA for a protein, and choosing its codons for a host while clearing sites it must not
  spell (`liulab_mbio.translate`); and barcode sets held a set distance apart
  (`liulab_mbio.barcodes`). Each has its own skill, `codon-optimize` and `barcode-design`,
  because both are wanted outside a library build.
- SrfI and PmeI join the shipped enzymes, rebuilt through the enzyme-data builder.
- Two research notes: the library method with every number sourced to the paper it comes from,
  and whether a barcode set needs limits on GC and on repeated bases. The second measured the
  evidence rather than following custom, and the answer is a cap on repeated bases and no GC
  band at all.
- Golden Gate cloning, end to end. `liulab_mbio.goldengate.plan_assembly` takes a vector and
  any number of inserts, picks a Type IIS enzyme with no site in the parts, designs the whole
  overhang set, checks every primer, simulates the assembly, and designs the colony PCR and
  sequencing that confirm the clone. `Plan.write` writes the annotated product, the primer
  order sheet, the protocol as JSON data (`write_protocol`), and an interactive HTML bench
  protocol rendered from that data, whose order sheet carries each oligo's verdict, and says
  which check fired and what it measured.
- `liulab_mbio goldengate plan` and `liulab_mbio protocol render` on the command line, and a
  repo-local `golden-gate-assembly` skill that calls them.
- SnapGene `.dna` read and write, with editing that carries features and primer binding sites,
  under one 0-based half-open coordinate model shared by every module.
- Restriction enzyme, codon usage and ligation fidelity data, each rebuilt by a script in
  `scripts/` and sourced in a note under `docs/research/`.
- Codon usage for human and mouse, beside *E. coli* K-12: `codon_usage("human")` and
  `codon_usage("mouse")`. Each counts one coding sequence per protein-coding gene, the one
  GENCODE marks as canonical, on hg38 and mm39.
- Primer design and evaluation against NEB's published rules, per polymerase.
- Fidelity scoring against a ligase matrix the user holds, by `--ligase-matrix` or
  `LIULAB_MBIO_LIGASE_MATRIX`. Nothing from that archive ships here.
- Docs: a walkthrough of GFP into pUC19, the plan it writes published as a live example, and
  an API reference covering every public module.

### Changed

- A plate step now names the medium the drug needs, and reads a marker the table does not know.
  `liulab_mbio.bench.phenotype` gained the markers Gateway's vectors carry: kanamycin, Zeocin
  and spectinomycin. Zeocin only works in low-salt medium, so `Phenotype.medium` sits beside the
  antibiotic and every method's plate step reads it. A resistance gene the table cannot name is
  now found anyway and named on the plate, rather than left out. A marker inside the piece a
  reaction throws away is skipped, because it is not what the plate selects.
- A check on a protocol page may now carry no verdict, and the page shows it as `not judged`.
  A check with no verdict used to be dropped from the page, which read as though nobody had
  asked the question. It now gets a badge of its own, in a muted colour, with the reason and
  what to go and look up beneath it. This holds for any protocol, not only a cloning one. The
  first check to use it asks whether two restriction enzymes can share one tube: the package
  knows the buffer each is sold in, and where those differ, answering needs figures for how much
  activity each keeps in the other's buffer, which nothing this package may ship states.
- An insert with no room for a junction primer no longer stops a colony PCR being designed. A
  junction primer anneals 100 bases inside the insert, so a linker or a tag is too short to hold
  one. `liulab_mbio.bench.validation.colony_pcr_check` used to refuse; it now leaves that primer
  off, keeps the flanking pair that reads across the junction, and says the gel cannot tell an
  insert that short from one the wrong way round. Golden Gate plans of a short insert are
  designed rather than refused for the same reason.
- The cloning methods are grouped under one verb. `liulab_mbio goldengate plan` is now
  `liulab_mbio cloning goldengate plan`, and `liulab_mbio cloning --help` lists the methods this
  package plans. `liulab_mbio library plan` is unchanged. The import path moves with the verb:
  `liulab_mbio.cloning.goldengate.plan_assembly`. The plan writes the same four files, with the
  same names and the same bytes. The `golden-gate-assembly` skill is now `molecular-cloning`,
  which picks the method for the job and reads `golden-gate/METHOD.md` once it has.
- Barcode sets are designed on an indel-aware distance by default. `liulab_mbio.barcodes` takes
  the metric as a dial: `sequence-levenshtein`, which counts a lost or gained base, or `hamming`,
  which counts mismatches and cannot see one. Measured over the 11-mer space, between 19% and
  46% of the single deletions of a Hamming set read as another barcode of the same set, against
  none of an indel-aware one; the indel-aware rule costs about one base of barcode length at a given set
  size, rejects nothing in the published set, and every set it designs is a Hamming set too.
  `deletion_ambiguity` reports that share for any set, and a library protocol now prints it where
  the barcode block is read back. A library plan draws different barcodes than it did, the same
  seed still returning the same ones. `docs/research/barcode-design.md` section 10 is the
  measurement and the reasoning.
- Each designed oligo is judged by what it is for, and `plan_assembly` takes its thresholds for
  each role. A sequencing primer passes from 16 bases, the length of Genewiz's own M13F primer,
  while a PCR primer still needs 18. A design picks the primer length that warns least on
  length, GC, GC clamp and Tm. A band no published rule sets, such as the GC clamp, reads as
  proposed on the order sheet.

### Removed

- The template's placeholder `greet()`, its CLI verb and its test.

### Fixed

- The reversed-insert lane of a colony PCR is now read off a plasmid that can exist. The lane
  used to come from the product with its insert turned over in place. At an end with an
  that put the overhang's bases on the wrong side of the join. The band was then off by the
  overhang, 231 bp where the real clone gives 227 with EcoRI. Now restriction and ligation turns
  the cut insert over and ligates it in, as the bench does.
- Golden Gate and Gibson no longer draw a lane for an insert the wrong way round. Neither can
  make one: no overhang in a Golden Gate set pairs with another backwards, and a Gibson insert
  shares its bases with the vector one way round only. With no orientation to tell, the colony
  PCR's reverse primer sits as close to its junction as the forward one, so both plans' bands
  are shorter than before.
- A plan no longer refuses a record that has a feature or primer across its origin when it has
  to read that record from the other strand. SnapGene often writes such a feature. Restriction
  and ligation refused any plasmid like this that an insert was cut from. Golden Gate and Gibson
  refused one given as an insert in reverse, and Gateway refused one whose first att site reads
  backwards. On the other strand, the feature still runs across the origin.
- A feature in pieces across the origin of a circular plasmid now keeps them in the order it
  reads. The GenBank reader sorted them by position, and so did turning a record over and
  building a plan's product. A feature that reads bases 91 to 98, then 2 to 8, came back reading
  2 to 8 first, so its bases came out of order. SnapGene files already kept the order.

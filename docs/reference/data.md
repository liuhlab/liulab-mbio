# What ships with it

Install the package and four sets of data come with it, in its own `data/` directory. Each one
is built by a script in `scripts/` from a named source. None is edited by hand.

Three other files the package reads are not included. You keep your own copy of each. They are
at the bottom of this page.

The calls below are Python, and [the Python reference](python.md) covers the rest. The same
data backs [the command line](cli.md).

## Enzymes

`enzymes.json` holds the restriction and Type IIS enzymes the package knows. There are 28.
Count them yourself:

```python
from mbio.enzymes import enzymes

len(enzymes())
```

Eighteen are plain restriction enzymes. Ten are Type IIS, which cut to one side of the site
they bind, so you choose the overhang they leave.

Each record holds:

- the recognition site, and where the enzyme cuts each strand
- the overhang that leaves, and whether the end is 5', 3' or blunt
- a supplier, a catalogue number and the commercial name
- the buffer that supplier sells the product in, and the incubation temperature
- the heat inactivation temperature and time, where the supplier states one
- sensitivity to Dam, Dcm and CpG methylation
- isoschizomers, the other enzymes that read the same site

Two kinds of source feed the file. The site, the cut offsets and the isoschizomers come from
REBASE, the open restriction enzyme database run by Richard Roberts and colleagues at New
England Biolabs. Version 609 is the one in the file.

The rest is what a supplier says about its own product: the name, the catalogue number, the
buffer, the temperatures, the methylation sensitivity. Each of those was read from that
supplier's catalogue page and typed in by hand, with the page cited. NEB's terms let you read
their pages for your own reference but not copy them, so no NEB file or table is reproduced
here. Twenty-six of the 28 are NEB products and two are Thermo Fisher's.

Every value names the source it came from. Where no source stated a value, the record says so
rather than filling the gap.

## Codon tables

`codon_usage.json` holds a codon table for three hosts. List them:

```python
from mbio.codons import codon_tables

codon_tables()
```

That returns `('e-coli-k12', 'human', 'mouse')`. The first is the default.

| Name | Organism | What was counted | Coding sequences | Codons |
| --- | --- | --- | --- | --- |
| `e-coli-k12` | *E. coli* K-12 MG1655 | `U00096.3`, every complete CDS | 4,317 | 1,342,016 |
| `human` | *Homo sapiens* | hg38, GENCODE v50, one canonical CDS per gene | 19,597 | 11,327,553 |
| `mouse` | *Mus musculus* | mm39, GENCODE vM39, one canonical CDS per gene | 21,479 | 11,894,511 |

All three are whole-genome counts. Every coding sequence counts once, whether the cell makes a
lot of its protein or none at all. A highly expressed reference set is a different thing. It
counts only the genes a cell makes most of, and it weights rare codons far lower. Each record
states which kind it is, and all three say whole genome.

The counts were made here, from sequence. No published codon usage table is copied. NCBI and
EMBL-EBI both state that they place no restrictions on passing on the records that were
counted. The bacterium came from the GenBank record. Human and mouse came from the assemblies
and GENCODE annotations the lab already holds.

Counting also gives the better table. The Kazusa entry usually reached for, *E. coli* K-12, rests
on 14 coding sequences, and one of its cells is zero. These tables have no zero cell at all. A
fraction or a log ratio taken from a zero is undefined.

## Ligation fidelity matrices

`ligation_fidelity.json` holds five matrices, one per enzyme: BsaI-HFv2, BsmBI-v2, Esp3I,
BbsI-HF and SapI. A matrix records how often each overhang was seen joining each other
overhang in one measured reaction. The first four cover four-base overhangs. SapI's are three
bases.

The numbers are supplementary tables S1 to S5 of Pryor et al., *PLOS One* 15(9):e0238592
(2020). The paper is open access under CC BY 4.0, so the tables ship with the credit the
licence asks for.

What the package does with a matrix is on [ligation fidelity](ligation-fidelity.md).

## Fonts

`fonts/` holds three typefaces cut down from the DejaVu 2.37 release: DejaVu Sans,
DejaVu Sans Bold and DejaVu Sans Mono. Labels on a map are set in the first two, and bases and
translations in the mono face. Each face ships twice. One file is a JSON table of advance
widths, kerning pairs and glyph outlines. The other is a WOFF2 subset for a page to embed. The
release's own `LICENSE` sits beside them, unchanged.

A figure carries the typefaces because of how a PNG or PDF is written. The package hands the
renderer an SVG in which every letter has already become its outline, a shape rather than a
character. The renderer draws nothing for text in a font it was not given, and it raises no
error while doing it, so a map would come back with its labels quietly missing. Outlines
remove the lookup, and the figure then draws the same on any machine. An HTML map is the other
way round: it keeps its letters as text you can search and copy, and embeds the WOFF2 subsets
so the text still lands where it was measured.

DejaVu's own changes are public domain. The glyphs it inherited are Bitstream's, from Bitstream
Vera, and Tavmjong Bah's, from Arev. Both grant free use and redistribution. Both ask that the
copyright and permission notice travel with the fonts, which is what the `LICENSE` beside them
is for.

## What does not ship

Three files the package reads are left out. Each is read from a copy you hold.

### The T4 ligase fidelity matrix

Potapov et al. measured T4 DNA ligase across all 256 four-base overhangs in 2018, which covers
the Type IIS enzymes nobody has published a matrix for. That archive is CC BY-NC 4.0, so none
of it is passed on here. Point at your own copy with `--ligase-matrix`, or with the
`LIULAB_MBIO_LIGASE_MATRIX` environment variable.
[Ligation fidelity](ligation-fidelity.md) says what it is used for.

### The 96 DMX barcode sequences

The kit's 96 barcodes are published with the kit, not here. Name your copy of that table with
the `LIULAB_SYNBIO_DMX_BARCODES` environment variable.
The [DMX barcode kit](dmx-barcode-kit.md) page covers what the file has to hold.

### A price record

A price record is a CSV file of what you pay: one row per item, with a charge, a currency, and
the quantity bands the charge holds over. A bench protocol turns it into a bill.

None ships. A tariff belongs to the lab that negotiated it, it goes stale, and a figure quoted
here would be ours rather than your supplier's. Point at yours with the `LIULAB_MBIO_PRICES`
environment variable. Without one a protocol still lists every item and quantity, and leaves
each money cell as a hole naming what went unpriced. Nothing is ever estimated.

## Licences cover files, not sequences

Four licences are named above. All of them sit on a file: a data archive, a paper's
supplementary tables, a set of typefaces. None of them reaches a DNA sequence. Bases carry no
licence, and an MTA or a UBMTA governs a plasmid someone shipped you rather than the sequence
of it. Nothing in this package weighs a licence when it picks a sequence.

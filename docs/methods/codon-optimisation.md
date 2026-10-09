# Codon optimisation

A protein can be spelled many ways in DNA, and a host does not use those spellings equally.
One command writes a protein as the codons one host counts most often, then takes out every
restriction site you name by swapping a codon for a synonym, so the protein itself is
unchanged. It reports every codon it moved.

## When to use it

- You are ordering a gene as synthetic DNA and want your host to read it well.
- A part must carry no BsaI, BbsI or other site, because the assembly cuts with it.
- You hold a coding sequence already and want to know which codons have to change.

A part written here goes straight into [Golden Gate assembly](golden-gate.md), whose enzyme
must not cut it, and [barcode sets](barcodes.md) are kept clear of the same sites.

## What you need

The protein itself, as one-letter amino acids typed on the command line rather than read from
a file. A host with a shipped codon table: `e-coli-k12`, `human` or `mouse`. The protein below
is GFP, read off a record that ships with the repo.

## Run it

```bash
pixi run mbio codon-optimize --kind protein --host e-coli-k12 \
  --forbid BsaI --forbid NdeI --forbid NcoI --name GFP --out plan/ \
  MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK
```

Every option comes first, because the sequence is the last thing on the line. It prints the
summary, then the two files it wrote:

```text
GFP: 714 bp, 238 aa, host e-coli-k12, 3 codon change(s), free of BsaI, NdeI, NcoI
plan/coding-sequence.dna
plan/codon-changes.tsv
```

Three codons moved. Writing GFP with the codons E. coli counts most often spelled two NdeI
sites and one NcoI site, and each went away when a histidine codon changed from CAT to CAC.
BsaI was forbidden too and needed nothing, because those same codons spell no BsaI site here.
Where no synonymous change can reach a site, the run stops rather than leave it in: a part
carrying it would cut itself.

`238 aa` and `714 bp` count no stop codon, because the protein handed over carries none. End
the protein with `*` and a stop codon is written.

## What it wrote

| File | What it is |
| --- | --- |
| [coding-sequence.dna](../examples/design/coding-sequence.dna) | the 714 bases to order, as a record. Opens in SnapGene |
| [codon-changes.tsv](../examples/design/codon-changes.tsv) | every codon it moved, and the site that moved it |

Both names link to what that run wrote, published here unedited. The same protein, host and
enzymes always give the same bases, so the file is a rerun away rather than something to keep
safe. To change the design, change the command and run it again.

## How to read the sheet

`codon-changes.tsv` holds one row for each codon moved: where it sits in the protein, the
codon before and after, the amino acid both of them spell, and the enzyme whose site the
change cleared. There are three rows here, codons 76 and 216 for NdeI and codon 230 for NcoI.

All three are histidine, CAT to CAC. Histidine has two codons, so there was one synonym to
move to. Read the amino acid column as the guarantee: the protein that came out is the
protein that went in.

`coding-sequence.dna` is the sequence to send a vendor, and it opens in SnapGene like any
other record. Check the name on it, which `--name` set, is the name you want on the order.

## Before you order

- Check the host is the one the protein will be made in, not the one you clone in.
- Name every enzyme the part will meet, including one you use only in a later round.
- Check the protein you pasted is the whole protein, both ends included.
- Decide whether you need a stop codon, and end the protein with `*` if you do.
- Check it starts where translation should, because nothing adds a start codon for you.

## Tips and troubleshooting

**A site cannot be removed.** The run stops and names the sequence, the site and the reason.
No synonymous change reached it, so either change the protein over that site or cut with a
different enzyme.

**You already have a coding sequence.** `--kind dna` checks what you hand it instead of
writing it again. Its codons are left as they are, apart from any that has to move to clear a
forbidden site.

**You want the bases on the screen.** Leave `--out` off, and the verb prints the summary, each
codon it moved, and then the sequence, so it pipes into whatever comes next.

## Where the numbers come from

A codon table here is a count of real coding sequences rather than a copy of a published
compilation. E. coli K-12 MG1655 is counted over every complete coding sequence in its
genome, fetched from NCBI. Human and mouse are counted over one canonical coding sequence per
protein-coding gene, from GENCODE on the assemblies the lab holds.

So each table is the background usage of a whole genome. It is not a set of highly expressed
genes, which is what a classical optimisation is built on. The genetic code that groups the
codons is the standard one.

`pixi run mbio codon-optimize --help` lists the rest.

## References

- NCBI E-utilities, spliced coding sequences for *Escherichia coli* str. K-12 substr. MG1655,
  GenBank U00096.3, read 14 September 2026.
- GENCODE v50 on GRCh38, one Ensembl canonical coding sequence per protein-coding gene, read
  14 September 2026.
- GENCODE vM39 on GRCm39, the same for mouse, read 14 September 2026.
- The standard genetic code, NCBI translation table 1.

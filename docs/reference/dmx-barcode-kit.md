# The DMX barcode kit

DMX is the lab's own read-back method. One design sits in one well, and each well is marked,
sequenced and called on its own. One of its two marking routes ligates four barcode plasmids
into the construct, in lysate, straight from the well. Those plasmids are this kit.

A barcode here marks a well. That is not the other sense of the word in this package, where a
barcode names a part in a library.

The package ships no barcode sequence. The 96 are read from a copy you hold. So this page is
the shape of the kit, and it prints no barcode either.

## Four groups of 24

The kit is 96 plasmids in four groups of 24. One barcode from each group goes into a well, so
four of them together name that well. Three of the groups tell the wells of one plate apart.
The fourth tells the plates apart.

All 24 members of a group carry the same two overhangs. What tells them apart is the UMI, 25
bases that differ between every plasmid in the kit.

## The four chain head to tail

Each group has its own pair of overhangs, read on the strand the cargo reads on:

| Group | Opens on | Closes on |
| --- | --- | --- |
| 1 | `AGGA` | `GTTC` |
| 2 | `GTTC` | `CCTT` |
| 3 | `CCTT` | `TCAG` |
| 4 | `TCAG` | `TTCC` |

Each group closes on the overhang the next one opens on. So the four can assemble in one order
only, and none of them can be left out. The chain closes through the design: group 4 ends on
`TTCC`, the design follows, and `AGGA` opens group 1 again.

That is five overhangs for four groups, because they are the boundaries and not the members.
Four things in a row have five ends: one before the first, one between each neighbouring pair,
and one after the last.

## One barcode

Each one is 97 bp. Only the UMI differs:

```text
GGTCTCA  [5' overhang]  [25 nt]  [25 nt UMI]  [25 nt]  [3' overhang]  CGAGACC
  BsaI                 constant               constant                  BsaI
```

That is 7 + 4 + 25 + 25 + 25 + 4 + 7 bases. The two BsaI sites face inwards, so one cut
releases the barcode on its two overhangs.

## The two universal primers

| Name | Sequence | Reads |
| --- | --- | --- |
| `dmx7` | `ATCGGTGACGGCGATTCTCACATTT` | forward, into the design |
| `dmx0` | `TTTGATACGCAGGAAGATGGCCCAC` | reverse, into the design |

`dmx7` is group 4's 3' constant region. `dmx0` is the reverse complement of group 1's 5'
constant. Both flank the design, and both exist only once the well is barcoded, because they
come from the barcodes. These two are the only sequences the package ships.

## Your own copy of the sequences

`read_kit` reads the 96 from a tab-separated file you name. Name no file and it reads the path
in `LIULAB_SYNBIO_DMX_BARCODES`.

| Column | What it holds |
| --- | --- |
| `name` | what the kit calls it, such as `DMX_1_1` |
| `group` | which of the four groups it is in, counting from one |
| `index` | which of the 24 in that group, counting from one |
| `overhang5` | the overhang its group opens on |
| `overhang3` | the overhang its group closes on |
| `umi` | the bases that tell it from the rest of its group |
| `final_seq` | the whole barcode, as it is ordered |

The columns may be in any order, and a missing one is named back to you.

The file is read only if all of this holds:

- 96 rows, one for each plasmid.
- No two rows share a UMI.
- Each group is indexed 1 to 24, nothing missing and nothing twice.
- All 24 members of a group carry one pair of overhangs.
- Each group closes on the overhang the next one opens on.

Anything else is refused, and the refusal says which of these failed.

Write both overhangs as they read on the strand the cargo reads on, the strand of the table
above. A plasmid map may show the other strand. That strand spells the chain backwards and
complemented, and a file written from it is refused.

## The calls

`read_kit` hands back a `Kit`. `Kit.group` gives one group's 24 barcodes in index order, and
`Kit.at` gives one barcode by its group and index. Both count from one. The arguments are in
the [Python reference](python.md).

## Ordering

The kit is deposited with Addgene, and four of its maps are held here:

| Barcode | Addgene |
| --- | --- |
| `DMX_1_1` | 255161 |
| `DMX_2_1` | 255185 |
| `DMX_3_1` | 255209 |
| `DMX_4_24` | 255256 |

Those four accessions are spaced exactly 24 apart, which implies
`accession = 255161 + 24 × (group − 1) + (index − 1)`, a formula inferred from four numbers and
not confirmed. Check it against Addgene's own kit listing before you order, and do not trust an
accession worked out from it.

The barcodes are used as supplied. Nothing is domesticated and nothing is changed.

The sequences come from Qian, Z. et al., *Nature Communications* (2026),
[doi:10.1038/s41467-026-76740-9](https://doi.org/10.1038/s41467-026-76740-9). Supplementary
Table 3 holds the UMIs, and Table 4 the primers.

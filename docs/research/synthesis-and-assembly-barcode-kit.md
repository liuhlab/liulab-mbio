---
search:
  exclude: true
---

# The DMX barcode kit

96 plasmids in four groups of 24. One barcode from each group marks a well, so a well is named
by a combination of four, and 24⁴ combinations is far more than any plate set needs. Used as
supplied — no domestication, no modification.

Source: Qian, Z. et al. *Nat. Commun.* (2026),
[doi:10.1038/s41467-026-76740-9](https://doi.org/10.1038/s41467-026-76740-9), Supplementary
Table 3 (the UMIs) and Table 4 (the primers). Deposited with Addgene.

## All 96 sequences

Held as a table extracted from Supplementary Table 3: name, group, index, both overhangs, the
25 nt UMI, and the full 97 bp barcode. All 96 UMIs are distinct.

Each row is 97 bp with the same layout, only the UMI varying:

```text
GGTCTCA  [oh5]  [25 nt constant]  [25 nt UMI]  [25 nt constant]  [oh3]  CGAGACC
  BsaI                 Seq5                           Seq3                BsaI
```

## How the four groups chain

The groups carry different overhangs, and those overhangs chain head to tail into the same
`AGGA`/`TTCC` site the cargo uses. So the four barcodes assemble in a fixed order:

```text
AGGA ── BC1 ── GTTC ── BC2 ── CCTT ── BC3 ── TCAG ── BC4 ── TTCC
```

A barcode cannot land in the wrong position, and none can be left out.

The diagram above is the strand the cargo reads on. Measured against the four depositor maps on
2026-10-06, the map strand runs the other way, `TTCC`-BC4-`TCAG`-BC3-`CCTT`-BC2-`GTTC`-BC1-`AGGA`,
and the extracted table is the reverse complement of the map. Reverse-complementing reverses the
order as well as the bases, so any statement of the chain has to say which strand it means. The
adjacency itself is the same either way, and it is derivable from all 96 rows of the table, one
overhang pair per group across all 24 members, not only from the four maps.

## Each group primes a different read

All eight constant regions either side of the UMI are primer landing sites. Six carry a
published DMX primer; the other two, group 1's 5' constant and group 4's 3' constant, are
`dmx0` and `dmx7`, for which the paper publishes a name and no sequence. The mapping below is
ours, derived from the published sequences, not stated by the paper:

| Group | 5' constant is | 3' constant is |
| --- | --- | --- |
| 1 | `dmx0` | DMX 1_rv |
| 2 | DMX 2_fw | DMX 3_rv |
| 3 | DMX 4_fw | DMX 5_rv |
| 4 | DMX 6_fw | `dmx7` |

Each region is one sequence across all 24 members of its group, and each published primer
matches its region exactly, forward or reverse-complemented. Measured over all 96 rows of
`dmx-barcodes.tsv` on 2026-10-06, then again on 2026-10-06 for the two unmatched.

## The two universal primers

`dmx0` and `dmx7` flank the design. The chain closes through it, so reading round the circle
gives `BC4`'s 3' constant, `GGAA`, the design, `TCCT`, then `BC1`'s 5' constant — which is the
`dmx7-design-dmx0` the amplicon table names. They exist only after barcoding, because they come
from the barcodes.

| Primer | Sequence | Reads |
| --- | --- | --- |
| `dmx7` | `ATCGGTGACGGCGATTCTCACATTT` | forward, into the design |
| `dmx0` | `TTTGATACGCAGGAAGATGGCCCAC` | reverse, into the design |

`dmx7` is group 4's 3' constant as the table holds it; `dmx0` is the reverse complement of
group 1's 5' constant. Both are 25 nt, like all six published primers, at 48% and 52% GC and
within 0.8 °C of each other — 70.7 and 71.5 °C in Q5 buffer.

**The names are the paper's, the sequences are ours by derivation.** The paper publishes no
sequence for either, so a mismatch against the authors' intent would show as a failed
amplification, not a wrong result.

That is why there are three primer pairs rather than three separate assays. Each pair enters
the same circle at a different rotation, so every amplicon carries the design plus all four
UMIs:

| Pair | Amplicon |
| --- | --- |
| DMX1 / DMX2 | `2-3-4-dmx7-design-dmx0-1` |
| DMX3 / DMX4 | `3-4-dmx7-design-dmx0-1-2` |
| DMX5 / DMX6 | `4-dmx7-design-dmx0-1-2-3` |

`dmx0` and `dmx7` are universal primers flanking the design. The three reactions are run
separately and pooled.

## Ordering

We hold four plasmid maps: DMX_1_1 (Addgene 255161), DMX_2_1 (255185), DMX_3_1 (255209) and
DMX_4_24 (255256). Their accessions are spaced exactly 24 apart, which implies:

```text
accession = 255161 + 24 × (group − 1) + (index − 1)
```

giving group 1 = 255161–255184, group 2 = 255185–255208, group 3 = 255209–255232,
group 4 = 255233–255256.

**Inferred from four data points, not confirmed.** Check it against Addgene's own kit listing
before ordering. The sequences themselves do not depend on it — those are measured from
Supplementary Table 3.

## What the kit is not

Barcoding happens in heat-lysed cell lysate on a throwaway aliquot, and the product is
PCR-amplified for sequencing, so a barcoded molecule never has to grow. The culture in the well
is the stock and is untouched. BsaI therefore does three separate jobs on the DMX vector —
plate barcoding, releasing cargo as an iGGA donor, releasing cargo into the working vector —
without those uses conflicting.

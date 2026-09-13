# Example: GFP into pUC19

Every file on this page was written by the command below, from the two sequence files in
`tests/data/`. Nothing here is edited by hand.

```bash
pixi run liulab_mbio goldengate plan tests/data/pUC19.dna tests/data/GFP.dna --out docs/examples/pUC19-GFP
```

To change what these files say, change the code and run that command again. The same inputs
write the same bytes, so a run that changes nothing leaves them alone.

| File | What it is |
| --- | --- |
| [protocol.html](protocol.html) | the bench protocol: eleven steps, reaction tables, expected bands |
| [product.dna](product.dna) | the finished plasmid, for SnapGene |
| [primers.tsv](primers.tsv) | the nine oligos to order |

[Put GFP into pUC19](../../golden-gate.md) explains how to read them.

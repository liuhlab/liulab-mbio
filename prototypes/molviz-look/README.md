# PROTOTYPE, throw away: what molviz copies from SeqViz's look

Throwaway branch, never merged. It answers [SeqViz beside SnapGene: build on it, or lay out from scratch?](https://github.com/liuhlab/liulab-mbio/issues/579).

One offline page with six views, flipped by the arrow keys: pUC19-GFP and the AP-1 library product (with mbio's 28 enzymes), each drawn by SeqViz 3.10.25 live, by mbio today, and as SnapGene Viewer showed it (window captures, kept out of git).

Rebuild, from the repo root, with SeqViz installed by `pixi exec -s nodejs -s pnpm -- pnpm add seqviz react react-dom` in an empty directory:

```bash
pixi run python prototypes/molviz-look/props.py docs/examples/pUC19-GFP/product.dna docs/examples/ap1-library/product.dna BUILD/props.json
pixi run mbio plot map docs/examples/ap1-library/product.dna --enzyme ... -o BUILD/mbio-crowded.html   # all 28 enzymes
python3 prototypes/molviz-look/build.py BUILD OUT.html SEQVIZ/node_modules/seqviz/dist/seqviz.min.js docs/examples/pUC19-GFP/product-map.html
```

`BUILD/` also needs `snapgene-map.jpg` and `snapgene-ap1.jpg`, captures of SnapGene Viewer's Map tab for the two files.

**Verdict, 2026-10-17:** copy SeqViz's sequence view and its interactions, keep SnapGene's map labels, and keep hiding and counting crowded labels. The full answer is on the issue.

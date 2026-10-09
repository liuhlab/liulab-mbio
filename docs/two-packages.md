# Two commands

One install gives you two commands.

| Command | What it does |
| --- | --- |
| `mbio` | the toolkit: cloning, primers, barcodes, codon choice, plasmid maps, bench protocols |
| `synbio` | the methods this lab works by, built from that toolkit |

Type `mbio` for a job with a standard name — Golden Gate, Gibson, restriction and ligation,
Gateway, a primer pair, a map. Type `synbio` for a method this lab named: iGGA, which builds
a library, and DMX, which reads one back.

```bash
pixi run mbio cloning goldengate plan vector.dna insert.dna --out plan/
pixi run synbio igga plan project.json --out library/
pixi run synbio dmx plan build.json --out readback/
```

Every command on this site starts with one of the two, so you never have to work out which.
From Python the names match:

```python
from mbio.cloning.goldengate import plan_assembly
from synbio.igga import plan_igga
from synbio.dmx import plan_dmx
```

Both come from the one install, `liulab-mbio`, and share its version number.

## Writing code for either one

Code that another method could use as it stands belongs in `mbio`. Code that holds one
method's own choices belongs in `synbio`. `CLAUDE.md` in the repo gives the test and
worked cases.

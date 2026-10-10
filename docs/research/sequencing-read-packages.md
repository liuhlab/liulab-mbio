---
search:
  exclude: true
---

# The packages that could align a sequencing read and judge it

Research note for issue #552, a sub-issue of map #551. It asks whether an existing package already
aligns a sequencing result against the product it should be, and judges it, so that this package
writes less. A result is a Sanger trace (`.ab1`), a whole-plasmid consensus, or a per-well
consensus an in-house route hands over. The baseline to beat is what is already installed:
Biopython ≥1.88, with `Bio.Align.PairwiseAligner` and the `abi` reader in `Bio.SeqIO`. File
formats, and what confirms a clone, are #553's.

Sources were read on **2026-10-09 and 2026-10-10** from each tool's own repository, docs, package
feed and changelog. Every measurement below was run on 2026-10-10 on one Apple-silicon laptop
(`osx-arm64`), through pixi: Biopython in this repo's environment, every other package with
`pixi exec` from conda-forge and bioconda. A time or a memory figure is true of that machine
once; the ratios between them are what carry.

## 1. The verdict, first

**No package does the job. Biopython already does the parts a package can do, so add no
dependency.** Every candidate stops at an alignment or a variant list; none gives a verdict per
junction and insert. Of the candidates that could be a dependency, only Tracy does a part
Biopython cannot do at all: it re-calls bases from the trace and splits a mixed trace into two
alleles. It is a command-line program with a linear reference, and it gives no verdict.

| Rank | Choice | What it adds over the baseline | Cost | Take it when |
| --- | --- | --- | --- | --- |
| 1 | **Biopython alone** (installed) | Nothing to add: it reads the trace, trims it (Mott), aligns either strand, and hands over the four dye channels for mixed peaks | None. `from Bio import Align, SeqIO` imports in 182 ms | Now |
| 2 | **mappy** (minimap2's Python binding) | Strand and rotation offset of a whole-plasmid consensus in 2 ms, and a difference string (`cs`) that lists each substitution and indel | bioconda only, compiled, MIT, 5 ms import | Only if a product too long for `PairwiseAligner`'s memory turns up (§3.4). Not now |
| 3 | **Tracy** (CLI) | Re-calls bases from the raw trace, splits a mixed trace into two alleles, writes variants as BCF | bioconda binary on both platforms, BSD-3; called as a subprocess; its BCF needs a reader | Only if mixed peaks must be split rather than flagged. Then it is an external tool, as `ipcr` is |
| 4 | edlib | Fast edit-distance alignment with a CIGAR | Unit costs only, no affine gaps; last GitHub release 2021 | No |
| 5 | parasail | SIMD Smith-Waterman | Same answer as `PairwiseAligner`, and no less memory (§4.2); last PyPI release 2023 | No |
| — | pysam | Reads SAM, BAM and VCF; aligns nothing | Compiled, bioconda | Only if a route hands over a BAM or VCF |
| — | wf-clone-validation, Sequeduct with Ediacara, sangerseqR, sangeranalyseR, TraceTrack, pydna | Each is out for a stated reason (§4) | — | Never as a dependency. wf-clone-validation's way of turning a circle is the one to copy (§5) |

**What stays ours whichever is chosen** (§5): turning the circle, choosing the strand, trimming
a trace that carries no qualities, the mixed-peak threshold, placing each disagreement on the
product's junctions and inserts, coverage of each, and the verdict and its thresholds.

## 2. What the job is, and who does which part

The job, from #552: read the file; align it against a circular reference, on either strand, a
whole-plasmid consensus starting anywhere; trim the unreadable ends; call substitutions, indels and
mixed peaks; give a verdict.

| Tool | Reads `.ab1` | Reads a consensus | Circular reference | Either strand | Trims ends | Substitutions and indels | Mixed peaks | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **Biopython 1.88** | yes | yes (FASTA, GenBank, SnapGene) | no; double the reference (measured, §3.3) | yes, align both and compare (measured) | Mott's modified algorithm, `abi-trim` | as alignment blocks; listing them is ours | channel data only (`DATA9`–`DATA12`, `PLOC2`); the call is ours | no |
| **mappy 2.31** | no | yes | no; double the reference (measured) | yes, reported | no | yes, `cs` string and `NM` | no | no |
| **edlib 1.3.9** | no | strings | no; infix mode on a doubled reference works (measured) | no, run both | no | CIGAR, unit costs | no | no |
| **parasail 1.3.4** | no | strings | no | no, run both | no | CIGAR | no | no |
| **Tracy 0.9.1** | yes, and SCF | no, traces only | **no** (measured: a rotated reference aligns half the read) | yes, keeps the better | sliding window, or fixed 50 + 50 bases by default | BCF from `decompose -v` | yes, `--pratio` 0.33 | no |
| **wf-clone-validation 1.8.4** | no | no, raw FASTQ or BAM | yes, by doubling the assembly and rotating it | yes | n/a | bcftools on a minimap2 alignment | no | ✓ or ✗ at ≥99 % identity and ≥95 % coverage, both defaults |
| **Sequeduct 0.4.6 + Ediacara** | no | no, raw nanopore FASTQ | doubles the reference; purpose unconfirmed | via minimap2 | NanoFilt | freebayes | no | PASS, FAIL, WARNING, LOW_COVERAGE |
| **sangerseqR 1.48.0** | yes | no | no | n/a | fixed bases | no table | `makeBaseCalls(ratio = 0.33)` | no |
| **sangeranalyseR 1.22.0** | yes | no | no | n/a | modified Mott, or a sliding window | contigs from reads, not against a reference | flags at 0.33 | no |
| **TraceTrack 1.0.0** | yes | no | no | yes | until 3 consecutive good bases | Clustal Omega alignment | peak-area ratio | coverage and identity, no pass or fail |
| **pydna 5.5.16** | yes, through Biopython | n/a | no | yes | no | no; keeps reads that contain the target exactly | no | no |

## 3. What the baseline does, measured

Inputs, all public at fixed addresses (§7): four traces from Biopython's own test suite; the
example trace and its reference from gear-genomics' SAGE, the one public trace found with a
reference it was read against; and a whole-plasmid nanopore consensus from wf-clone-validation's
test data. The script ran each through Biopython 1.88 in this repo's environment.

### 3.1 The reader

| Trace | Bases | Read time | Mean quality | `abi-trim` keeps | Non-ACGT calls | Interior peaks with second/first ≥0.3 |
| --- | --- | --- | --- | --- | --- | --- |
| `3730.ab1` | 1165 | 8 ms | 44.8 | 1075 | 7 | 3 |
| `3100.ab1` | 795 | 4 ms | 46.8 | 678 | 0 | 30 |
| `A6_1-DB3.ab1` | 839 | 5 ms | 52.0 | 820 | 21 | 14 |
| `310.ab1` | 868 | 6 ms | 0.0 | **0** | 265 | 401 |
| SAGE `sample.abi` | 1451 | — | 0.0 | **0** | 247 (`N`) | — |

- **Every trace carries the four channels and the peak positions** (`DATA9`–`DATA12`,
  `PLOC2`, dye order in `FWO_1`). A second-to-first peak ratio at each called base took ten
  lines. The 0.3 used here is a probe, not a threshold. Tracy, sangerseqR and sangeranalyseR all
  default to 0.33, and this survey found no measurement behind it.
- **A trace without qualities trims to nothing.** Two of five carry all-zero qualities, and
  `abi-trim` returns an empty record rather than falling back. A policy for that trace is ours.
- **The vendor's own calls can be poor.** SAGE's trace carries 247 `N` of 1451 calls. Aligned
  as called, it disagrees with its reference at 383 positions. Tracy, which re-calls bases from
  the channels, wrote a read with no `N` at all.

### 3.2 Choosing the strand

SAGE's trace against its own reference, local mode, `blastn` scoring: score 1113 on the forward
strand and 26 on the reverse, 0.02 s each. The consensus against itself scored 12722 and 42.
Biopython's own traces, against pUC19, which they were not read from, scored 20 to 24 on both
strands: that is what no match looks like. So aligning both and keeping the better is enough
where the read belongs to the product. Nobody publishes a margin below which the choice is
unsafe; that is a hole.

### 3.3 Turning the circle

| Case | As given | Reference doubled |
| --- | --- | --- |
| SAGE trace, reference rotated so the read crosses its origin (Biopython, local) | read 2–721 of 1451 aligned; the rest lost | read 2–1449 aligned, reference 791–2234 of a 1509 bp sequence |
| Same, Tracy `align` | 722 reference bases aligned; the read's last 730 face end gaps | not tried; Tracy takes one reference |
| Same, mappy `map-ont` | one hit, read 2–721 | read 2–1451, reference 791–2236 |
| 6361 bp consensus, rotated by 4000, reverse-complemented, one substitution and one 3 bp deletion planted (Biopython, local) | — | reverse strand, reference 4003–10364, 1 mismatch, 3 gap bases; 0.72 s a strand |
| Same, global mode with free end gaps | — | right placement, but `counts()` adds the end gaps: 6364 gap bases. Use local mode |
| Same consensus against `tests/data/pUC19.dna`, global, unrotated | identities 1855, gaps 4995: no answer | — |

**Doubling the reference gives ADR 0001's form for free.** The span 791–2234 on a 1509 bp
sequence starts inside the record and ends past its length, which is exactly how
`docs/adr/0001-coordinates.md` writes a span across the origin. The conversion at the boundary
is one modulo on the start.

The real consensus against pUC19 gives a partial local hit: 862–2613, with 1582 identities, 123
mismatches and 94 gap columns. It shares a pUC-like backbone and little else. So a verdict has to
ask whether each expected span is covered, not whether the best score is high.

### 3.4 Memory, which is the baseline's one real limit

`PairwiseAligner` keeps a full traceback matrix. One local alignment of an n bp query against a
2n bp doubled reference:

| n | Time | Peak memory |
| --- | --- | --- |
| 3000 | 0.16 s | 68 MB |
| 6000 | 0.65 s | 179 MB |
| 10000 | 1.85 s | 443 MB |
| 15000 | 4.25 s | 958 MB |
| 15000, `score()` only | 2.55 s | 31 MB |

Memory grows with the square of the plasmid. `score()` keeps linear memory, so the strand can
be chosen by score and only one strand aligned. Anchoring the rotation on a shared word, then
aligning undoubled, would cut the matrix fourfold; that was not measured. A Sanger read is a few
hundred bases against at most the product, so this limit bites only on whole-plasmid consensus
sequences of large products.

## 4. The candidates

### 4.1 Install cost and whether the package may depend on it

| Package | Channel, both `osx-arm64` and `linux-64` | Kind | Warm import | Licence | Latest release | Last commit |
| --- | --- | --- | --- | --- | --- | --- |
| Biopython | conda-forge, yes | compiled, installed | 182 ms (`Align`, `SeqIO`) | Biopython licence | 1.88, PyPI 2026-08-12 | 2026-10-09 |
| mappy | bioconda, yes | compiled (minimap2) | 5 ms | MIT | 2.31, 2026-05-19 | 2026-05-19 |
| edlib (`python-edlib`) | bioconda, yes; conda-forge `edlib` has no `osx-arm64` | compiled | 13 ms | MIT | 1.3.9.post1 PyPI 2024-09-04; GitHub release v1.2.7, 2021 | 2025-05-13 |
| parasail (`parasail-python`) | bioconda, yes | compiled | 220 ms | BSD | 1.3.4, PyPI 2023-02-17 | 2024-09-04 |
| pysam | bioconda, yes | compiled (htslib) | 64 ms | MIT | 0.24.1, 2026-09-07 | 2026-09-07 |
| Tracy | bioconda, yes; ran on `osx-arm64` | C++ binary | n/a, subprocess | BSD-3-Clause | 0.9.1, 2026-08-18 | 2026-10-02 |
| sangerseqR, sangeranalyseR | bioconda `noarch`, R | R packages | n/a | GPL-2; sangeranalyseR's LICENSE says MIT, its DESCRIPTION GPL-2 | 1.48.0; 1.22.0 | sangeranalyseR 2026-05-06 |
| wf-clone-validation | not a package: Nextflow, Docker or Singularity | Nextflow program | n/a | Oxford Nanopore Technologies PLC Public License 1.0 | 1.8.4, 2025-12-15 | 2025-12-15 |
| Sequeduct, Ediacara | not on either channel; Nextflow | Nextflow program; Python side on PyPI | n/a | GPL-3.0 | Sequeduct tag 0.4.6; Ediacara PyPI 0.2.2, 2023-08-23 | 2025-12-12; 2025-02-15 |
| TraceTrack | not a package: Flask, Celery, Redis web app | web app | n/a | GPL-3.0 | 1.0.0, 2022-08-03 | 2023-07-26 |
| pydna | bioconda `noarch`, but at 5.2.0 (2024) against PyPI's 5.5.16 | pure Python | not measured | BSD | 5.5.16, 2026-06-09 | 2026-09-29 |

Import times are the best of three fresh interpreters. Cold first imports through `pixi exec` ran
far longer (parasail 4.4 s, pysam 6.7 s), which matters because `src/` is imported at test time;
any of these would be imported inside the function that needs it.

### 4.2 The compiled aligners on the planted consensus

Same 6361 bp consensus, rotated, reverse-complemented, one substitution and a 3 bp deletion,
against the doubled original.

| Aligner | Placement found | Difference reported | Time | Process peak memory after |
| --- | --- | --- | --- | --- |
| Biopython `PairwiseAligner`, local | 4003–10364, reverse | 1 mismatch, 3 gap bases | 0.72 s a strand | 236 MB |
| edlib, infix mode (`HW`) | 4003–10363, reverse | `3358=1X2000=3D999=`, distance 4 | 0.002 s | no rise over the 127 MB held after importing all four |
| parasail `sw_trace_striped_16` | ends at 10363, reverse | `4003D3358=1X1999=3D1000=`, the offset as a leading deletion | 0.072 s | 442 MB |
| mappy, `asm5` | 4003–10364, strand −1 | `:3358*cg:1999-aat:1000`, `NM` 4 | 0.002 s | not separated |
| mappy, `asm5`, reference not doubled | split at the origin: 0–4003 and 4003–6361 | two hits to join | — | — |

All four find the same answer. Speed is the only difference, and at plasmid size Biopython's
0.72 s does not matter. Memory, at the sizes in §3.4, is the only reason to reach for mappy.

### 4.3 Why each Nextflow program, R package and app is out

- **wf-clone-validation.** It starts from raw nanopore reads, not the consensus a vendor returns,
  and runs only under Nextflow with containers. Its licence is Oxford Nanopore's own Public
  License 1.0. Its comparison is the useful part: with `full_reference` it doubles the assembly
  (`seqkit concat`), maps the reference to it with minimap2 `asm5`, reads the strand and offset,
  and rotates with `seqkit restart`. Then it aligns again and calls variants with bcftools
  (`main.nf` lines 146–208 and 444).
- **Sequeduct and Ediacara.** Raw nanopore reads again, through NanoFilt, minimap2, freebayes,
  bcftools and Canu. GPL-3.0. Ediacara gives the one ready verdict in the survey (PASS, FAIL,
  WARNING, LOW_COVERAGE, low-depth cutoff 30), but its last PyPI release is from 2023.
  `docs/research/dnachisel-evaluation.md` already read the stack as prior art, not as
  dependencies that will receive fixes, and nothing here changes that.
- **sangerseqR and sangeranalyseR.** R. sangeranalyseR assembles reads into a contig; its only
  reference is an amino-acid sequence for fixing frameshifts. sangerseqR's `setAllelePhase` takes
  a reference sequence but gives no variant table and no verdict.
- **TraceTrack.** A web app, effectively dormant since 2023, requiring a CDS in every reference.
- **pydna.** `Dseqrecord.map_trace_files` keeps the traces that contain the target exactly, on
  either strand, and adds each as a feature. It neither trims nor lists a difference.
- **Midnighter/plasmid-verification** (Apache-2.0, last commit 2025-08-14) turned up in
  discovery. Its README lists no features and it is not on PyPI, so it was not judged.

## 5. What stays ours, whichever is chosen

| Part | Why no candidate takes it |
| --- | --- |
| **Turning the circle.** Double the product (or rotate a consensus to it), align, and convert the span at the boundary by one modulo | Tracy, mappy, edlib and parasail all treat the reference as linear (§3.3, §4.2). wf-clone-validation does it, but only inside Nextflow |
| **Choosing the strand.** Align both, keep the better; `score()` first if memory matters | Biopython, edlib and parasail do not choose. Tracy and mappy report it, and the rest of this row still applies to them |
| **Trimming a trace with no qualities** | Biopython's `abi-trim` returns nothing (§3.1); Tracy's fixed 50 + 50 is a default, not a source |
| **The mixed-peak call and its threshold** | The channels are there (§3.1). The 0.33 every tool uses is unsourced, so a check on it carries no verdict until someone finds a source |
| **Placing each disagreement on the product.** A position, then the junction or insert it falls in, under the four cloning methods' spellings (map #551) | No tool reads features. Biopython's alignment blocks give the positions |
| **Coverage of every junction and insert** | A best score hides a missing span (§3.3). Only wf-clone-validation and Ediacara compute coverage, over the whole plasmid |
| **The verdict and its thresholds** | wf-clone-validation's 99 % identity and 95 % coverage are that program's defaults, with no source given. Ediacara's cutoffs are the same kind |

What the baseline already gives, and none of this has to rebuild: reading the trace, its
channels and its qualities; Mott trimming where qualities exist; reading FASTA, GenBank and
SnapGene; affine-gap alignment on either strand, with exact coordinates and `counts()`.

## 6. Holes

| Hole | What stood in |
| --- | --- |
| No public whole-plasmid consensus with the plasmid it should be. Addgene's sequence files now need a login (its blog, "Protecting the science you share: why sequence data requires a login"), and every scripted fetch returned 404; wf-clone-validation's test assemblies come with no full reference | The real consensus with a planted rotation, strand flip, substitution and deletion |
| No public trace that reads across a plasmid's origin | SAGE's trace against its reference rotated by 750 bp |
| Tracy's variant calls were run once, on SAGE's example (53 PASS and 383 LowQual records), never on a plasmid | — |
| Tracy's paper (Rausch et al. 2020) was not read; only its repository and `--help` | — |
| Why Sequeduct doubles its reference | The script name, `create_2x_fasta.py` |
| A safe margin between the two strands' scores; a mixed-peak ratio with a source; an identity or coverage threshold with a source | Nothing; each is a hole for #554 to carry as one |
| sangerseqR's last commit date: its Bioconductor git was unreachable | The Bioconductor release date |

## 7. Sources

All read 2026-10-09 or 2026-10-10.

### Measured inputs

- Biopython test traces `3730.ab1`, `3100.ab1`, `310.ab1`, `A6_1-DB3.ab1`:
  <https://github.com/biopython/biopython/tree/e38c64daa2/Tests/Abi>
- SAGE example trace and reference, `server/sample.abi` and `server/sample.fa`:
  <https://github.com/gear-genomics/sage/tree/459d339300/server> (GPL-3.0)
- wf-clone-validation test consensus
  `test_data/workflow_glue/find_inserts/assemblies/barcode01.final.fasta` (6361 bp):
  <https://github.com/epi2me-labs/wf-clone-validation/tree/b3bf4ee47f/test_data>
- Addgene's login notice:
  <https://blog.addgene.org/protecting-the-science-you-share-why-sequence-data-requires-a-login>
  (title read from the link on <https://www.addgene.org/browse/sequence/91118/>)
- pUC19, `tests/data/pUC19.dna` in this repo

### Biopython 1.88

- `Bio/SeqIO/AbiIO.py`, `_abi_trim`: "Richard Mott's modified trimming algorithm", read from the
  installed 1.88
- `Bio/Align/__init__.py`, `PairwiseAligner.score` and `Alignment.counts`, the installed 1.88
- <https://pypi.org/pypi/biopython/json>; <https://anaconda.org/conda-forge/biopython>

### Aligners

- mappy and minimap2: <https://github.com/lh3/minimap2>, v2.31 (2026-05-19);
  <https://anaconda.org/bioconda/mappy>
- edlib: <https://github.com/Martinsos/edlib>; <https://pypi.org/pypi/edlib/json>;
  <https://anaconda.org/bioconda/python-edlib>; <https://anaconda.org/conda-forge/edlib>
- parasail: <https://github.com/jeffdaily/parasail-python>; <https://pypi.org/pypi/parasail/json>;
  <https://anaconda.org/bioconda/parasail-python>
- pysam: <https://github.com/pysam-developers/pysam>, v0.24.1;
  <https://anaconda.org/bioconda/pysam>

### Trace tools

- Tracy: <https://github.com/gear-genomics/tracy> at 1365b0ebf7: README, LICENSE,
  `docs/cli/README.md`, `src/tracy.cpp`, `src/trim.h`, `src/sage.h`, `src/indigo.h`; v0.9.1 run
  from bioconda, `tracy align` and `tracy decompose --help`; paper
  <https://doi.org/10.1186/s12864-020-6635-8> (not read)
- sangerseqR: <https://bioconductor.org/packages/release/bioc/html/sangerseqR.html> and its
  1.48.0 reference manual
- sangeranalyseR: <https://bioconductor.org/packages/release/bioc/html/sangeranalyseR.html>;
  <https://github.com/roblanf/sangeranalyseR>: DESCRIPTION, LICENSE, `NEWS.md`,
  `docs/source/content/function_manual.rst`
- TraceTrack: <https://github.com/Merck/TraceTrack> at 7f69f179: README, LICENSE,
  `tracetrack/alignment_utils.py`, `tracetrack/entities/record.py`, `templates/help.html`
- pydna: <https://github.com/pydna-group/pydna> at 1685b755, `src/pydna/dseqrecord.py`,
  `LICENSE.txt`; <https://pypi.org/pypi/pydna/json>; <https://anaconda.org/bioconda/pydna>

### Whole-plasmid Nextflow programs

- wf-clone-validation: <https://github.com/epi2me-labs/wf-clone-validation> at b3bf4ee47f
  (v1.8.4): README, LICENSE, `main.nf`, `nextflow_schema.json` (`full_reference`,
  `insert_reference`, `expected_identity` default 99, `expected_coverage` default 95),
  `CHANGELOG.md`
- Sequeduct: <https://github.com/Edinburgh-Genome-Foundry/Sequeduct> at d3544eae: README,
  `install_sequeduct.sh`, `nextflow/sequeduct_analysis.nf`; Ediacara:
  <https://github.com/Edinburgh-Genome-Foundry/Ediacara>, README and `ediacara/Comparator.py`;
  <https://pypi.org/pypi/ediacara/json>; paper <https://doi.org/10.1021/acssynbio.3c00589>
- Midnighter/plasmid-verification: <https://github.com/Midnighter/plasmid-verification>, README
  and source tree

### Package feeds

`https://api.anaconda.org/package/<channel>/<name>` for every channel and
platform claim above.

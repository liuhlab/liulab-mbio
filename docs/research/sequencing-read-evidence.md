---
search:
  exclude: true
---

# What a sequencing result holds, and what confirms a clone

Research note for issue #553, a sub-issue of map #551. It asks four things of each kind of
sequencing result a construct comes back as — a Sanger trace, a whole-plasmid result, and what a
DMX well hands over: what the files hold, what makes one base in them trustworthy, what published
sources call a confirmed clone, and which public example files are small enough for a test suite
held to 30 s. Sibling #552 surveys which packages could do the aligning; nothing here weighs one.

Everything below was read on **2026-10-09 and 2026-10-10** from each owner's own documentation,
or measured on a file this session downloaded. Two `.ab1` files were opened with Biopython 1.88,
the version this repo pins, and their directories walked tag by tag.

## 1. The verdict, first

| | Sanger trace | Whole-plasmid result | DMX well |
| --- | --- | --- | --- |
| **What arrives** | one `.ab1` per read: calls, one quality value per call, peak positions, four raw and four processed channels; vendors add a text sequence and sometimes a quality file (section 2) | a folder per sample: consensus FASTA and GenBank, a per-base table, a coverage plot, a read-length histogram, the aligned raw reads as FASTQ, a summary, and an `.ab1` drawn from the reads (section 3) | a read count and the consensus sequences someone has already called, `judge_well(reads=..., called=...)` (section 4) |
| **Per-base trust** | a Phred quality value per call (section 2.4) | read depth and the base counts at each position; GENEWIZ's consensus FASTQ and `.ab1` carry a quality per base, on a scale no vendor states (section 3.4) | none: the consensus arrives with no per-base support |
| **The ends** | the first bases after the primer and the tail are low quality; Biopython trims both by Mott's rule, cutoff 0.05 | none: a circular consensus has no read ends | none |
| **A mixture shows as** | two peaks at one position, written as an IUPAC code by the KB basecaller | mixed base counts in the per-base table and the `.ab1`; a second species in the histogram; only the commonest species gets a consensus | more than one consensus, which `identity_check` already fails |
| **Public example small enough for tests** | yes: Biopython's test traces, 0.2–0.3 MB each (section 6) | no raw vendor result is public; one processed per-base table is (section 6) | not applicable: the input is two values |

**The strongest single finding.** A whole-plasmid result is already a judged consensus: the
vendor aligned the reads, called each base and reported what disagreed. A Sanger read is not. It
is one strand's calls with their own quality, and only its middle is trustworthy. So the two
kinds need different readers but can meet in one comparison — a called sequence, with a
trustworthy span, held against the expected product — and a DMX well is that same comparison with
the trustworthy span being the whole call.

## 2. A Sanger trace

### 2.1 What the format promises

The format is Applied Biosystems' ABIF, *Applied Biosystems Genetic Analysis Data File Format:
ABIF File Format Specification and Sample File Schema*, July 2006 (PDF metadata: created
2005-04-13, modified 2006-07-14; © 2006 Applied Biosystems). It is a tagged binary: every item is
a four-character tag and a number, "analogous to the keys in a (key, value) mapping". Two
statements in it bind any reader:

- "The ABIF file format by itself does not specify the schema … which tags are written and when.
  These schema are specific to the instrument and software version which created the file."
- "Applied Biosystems provides no guarantee that all tagged data elements will be present,
  consistent, or supported in future versions of the software."

The tags a reader of base calls needs, from the spec's Table 9 (SeqScape v2.5) and the
Sequencing Analysis v5.2 table:

| Tag | Number | What the spec says it holds |
| --- | --- | --- |
| `PBAS` | 1, 2 | "Basecalled sequence (edited)"; "Basecalled sequence" |
| `PCON` | 1, 2 | "Per-base quality values (edited)"; "Per-base quality values" |
| `PLOC` | 1, 2 | "Base locations (edited)"; "Base locations" — the scan at each call's peak |
| `DATA` | 1–4, 9–12 | channels 1–4 raw data; channels 1–4 analysed data |
| `FWO_` | 1 | "Base order" — which channel is which base |
| `phQL` | 1 | "Maximum quality value" |
| `phTR` | 1, 2 | the first and last base of a trim region, and the "Trim probability threshold used" |
| `SMPL` | 1 | "Sample name" |
| `S/N%` | 1 | "Signal level for each dye" |

Number 2 is what the basecaller wrote and number 1 is what a person edited. A reader judging a
clone wants the basecaller's. The spec defines no scale for `PCON`; the phd1 equivalence it
names for `phQL` places it on phred's scale (section 2.4).

### 2.2 What two real files hold

Both are Biopython test files (section 6), read with Biopython's `abi` and `abi-trim` readers and
a walk of the ABIF directory.

| | `3730.ab1` | `310.ab1` |
| --- | --- | --- |
| Size | 299,987 bytes | 222,099 bytes |
| Instrument (`MCHN1`, `MODL1`) | ABI-3730-XL, `3730` | ABI PRISM 310, `310` and a trailing space |
| Basecaller (`SVER2`, `APrN1`) | `KB 1.2`; `3730BDTv3-KB-DeNovo_v5.2` | not written; `BDTpop6Rapid` |
| Directory entries | 123 | 113 |
| Calls (`PBAS2`) | 1,165 | 868 |
| Quality (`PCON2`) | 0–55, mean 44.8; 1,089 calls at Q20 or above, 978 at Q30 | **every value 0** |
| Non-ACGT calls | 7 IUPAC codes (`K`, `Y`, `R`), all outside the trimmed span | `N` throughout |
| Trace length (`DATA9`) | 16,302 scans | 9,826 scans |
| `phTR` trim region | `(-1, -1)`, threshold `-1.0`: no trim written | `(-1, -1)`: no trim written |
| Biopython `abi-trim` keeps | 1,075 calls, `[14, 1089)` | 0 calls |
| What it reads | sample `226032_C-ME-18_pCAGseqF`; the calls, translated here, are mCherry, a GGSG linker, then EGFP | not identified |

Four things a reader has to handle, all seen here rather than read on a page:

- **A file may carry no quality.** The 310 trace has `PCON` present and every value zero. Any rule
  of the form "trust Q20 and above" trusts nothing in it.
- **The file does not trim itself.** Neither file wrote a trim region, so the trustworthy span is
  the reader's to compute.
- **Edited and original tags can both be present and identical.** `PBAS1` equals `PBAS2` in both.
- **No tag holds QV20+, CRL or a trace score** in either file, and the 2006 specification names
  none. Those are per-read numbers the vendor's software reports beside the file (section 2.3).

### 2.3 What a vendor sends with it

| Owner | What it says comes back |
| --- | --- |
| GENEWIZ DNA sequencing FAQ (<https://www.genewiz.com/public/resources/faqs/faqs-dna-sequencing>) | both the ".seq and .ab1 files"; an "Analysis Result" PDF that gives "tips and suggestions" |
| GENEWIZ Sanger service page (<https://www.genewiz.com/public/services/sanger-sequencing>, undated) | "Quality Score and Contiguous Read Length (CRL) provided in each read"; "read lengths up to ~1000 bases. A typical read will provide 800 bases Phred20" — the figure `domestication-methods.md` read on 2026-10-06, unchanged |
| GENEWIZ blog, 2025-02-17 (<https://blog.genewiz.com/analyzing-sanger-sequencing-data>) | "The quality score (QS) is the average QV for all peaks in the trace with an assigned base" |
| Applied Biosystems, *Sequencing Analysis Software v5.4 Quick Reference Card*, PN 4401738 Rev. B, 4/2009 (<https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_064539.pdf>) | what the vendor's software can export: the `.ab1`; a "Text file of the sequence … (.seq; all bases or only clear range bases)"; "Phred (.phd.1) files"; "Standard chromatogram format (.scf) files"; analysis reports as text, HTML, PDF or XML |
| Thermo Fisher *Behind the Bench* blog, 2015-07-31 (<https://www.thermofisher.com/blog/behindthebench/?p=9258>) | QV20+: "Total number of bases in the entire trace that have a basecaller Quality Value >=20"; Trace Score: "Average of basecall quality values for bases in the clear range" |

So the `.ab1` is the one file every vendor sends, and the only one that carries per-base quality.
A `.seq` is the calls alone, sometimes cut to the clear range. The summary numbers (QS, CRL,
QV20+, trace score) are per read, computed by the vendor's software, and **not in the `.ab1`
files inspected** (section 2.2). **CRL has no published definition** that this survey reached.

### 2.4 Where a read cannot be trusted

| Number | What it means | Owner |
| --- | --- | --- |
| QV = −10 log10(Pe); QV 20 is a 1% chance of a wrong call, 30 is 0.1% | the quality scale | AB v5.4 card, above; Richterich 1998, *Genome Res.* 8:251, doi:10.1101/gr.8.3.251, quoting phred: "q = −10 × log10(p)" |
| QV 0–14 "Are not acceptable"; 15–19 "Need manual review"; 20 and above "Are acceptable" | the software's default display bands | AB v5.4 card |
| "Typical high-quality pure bases have QVs of 20 to 50"; "mixed bases … 10 to 50"; the scale runs 1–99 | what a good call scores | AB v5.4 card |
| Clear range: what "remains after excluding the low-quality or error-prone sequence at both the 5´ and 3´ ends" | the trustworthy span, computed by the software | AB v5.4 card |
| "The first 20 to 40 bases are typically not well resolved"; primers "at least 60 bp, preferably 100 bp, away from key bases" | the start of a read | GENEWIZ blog, 2025-02-17 |
| "The quality of the first 25 to 35 bases of the sequence is expected to be poor"; "Sequence lengths of 800 to 900 bases are possible for high-quality templates" | the start and the length | MSU RTSF Genomics Core, *Sanger Sequencing Best and Worst Practices Guide*, 25 April 2024 (<https://rtsf.natsci.msu.edu/_assets/files/genomics/Sanger_Sequencing_Best_and_Worst_Practices_Guide_25April2024.pdf>) |
| Dye blobs: "Large broad peak normally seen at 85–90 bp or 125–130 bp", in severe cases also "~ 60–65 bp" | an artefact that hides calls inside an otherwise good start | Thermo Fisher user bulletin MAN0014435 Rev. A.0, 15 Jan 2016 (<https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0014435_Trbleshoot_Sanger_seq_data_UB.pdf>) |
| "Large peaks (blobs) in the first 120 bases" | the same | Applied Biosystems, *DNA Sequencing by Capillary Electrophoresis Chemistry Guide*, PN 4305080, ch. 8, read from the University of Notre Dame genomics core's copy (<https://genomics.conductor.nd.edu/assets/324286/troubleshootingdna_sequencing_by_capillary_electrophoresis_chemistry_guide_pn_4305080_.pdf>), revision not printed |
| Best region "typically from base 150 to 200"; errors "start to exceed 10% between bases 300 and 700" | error along a read, measured on 1998 instruments and gels | Richterich 1998 |
| Mott trimming with `cutoff = 0.05` and a 20-base minimum | Biopython's `abi-trim` | Biopython 1.88, `Bio/SeqIO/AbiIO.py`, `_abi_trim` |
| A homopolymer "typically observed in stretches >9 bases" leads to mixed sequence or baseline noise | where a Sanger read loses its place | MAN0014435 |

Three things follow. **The start and the end of a read are always low quality, but no owner gives
one number for how many bases**: 20–40 (GENEWIZ), 25–35 (MSU), blobs to 120 (Applied
Biosystems). A reader computes the trustworthy span from the qualities and does not assume a
fixed margin. **Q20 is the one threshold the owners agree on** — Applied Biosystems' "acceptable"
band, GENEWIZ's read-length unit. **Mott's 0.05 is a Biopython default**, not a published rule
this survey could trace further: phrap.org, which Biopython's docstring cites, refused the
connection, and Ewing and Green's full texts were behind a login.

### 2.5 What a mixed peak looks like

| Owner | Words |
| --- | --- |
| Chemistry Guide PN 4305080, ch. 8 | the "2nd highest peak threshold for mixed base identification … The recommended range is 15 to 25%" |
| MAN0014435 | a second sequence: "Mixed sequence content throughout the length of the trace"; off-target priming: "Mixed sequence content after the primer region"; a heterozygous insertion or deletion: "Mixed sequence content starting at a specific point" |
| MSU guide | "Mixed template means that there are multiple peaks (base calls) at a single position" |

In the file, a mixed position is an IUPAC code in `PBAS` with its own `PCON` value. The 3730
trace (section 2.2) carries seven, all within the first 14 calls or the last 76, where the trim
removes them. **A mixed call in the trimmed span is the Sanger sign of two templates**, and a
whole-trace mixture is a mixed colony or prep, not a base to call.

## 3. A whole-plasmid result

Plasmidsaurus and Eurofins say they sequence on Oxford Nanopore; GENEWIZ's Plasmid-EZ page says
only "long-read sequencing". Each returns a consensus, not reads to assemble. `docs/research/domestication-methods.md` section 7 covers price and turnaround; this
section covers the files.

### 3.1 Plasmidsaurus

From *Technical Documentation: Plasmid*, `PS-0002-E | v1.2 | Revised Sept. 22, 2026`
(<https://plasmidsaurus.com/technical-documentation/plasmid>).

| File | What it holds, in the vendor's words |
| --- | --- |
| `<ORDERCODE>_fasta-files/<OSID>_<SAMPLE>.fasta` | "Polished plasmid consensus sequence generated from the raw sequencing reads" |
| `<ORDERCODE>_genbank-files/….gbk` | the polished consensus, annotated |
| `<ORDERCODE>_per-base-data/<OSID>_<SAMPLE>.tsv` | "Position-level comparison between the raw reads and consensus sequence. Includes the consensus base, aligned-read depth, and the distribution of matches, mismatches, insertions, deletions, and individual nucleotide calls" |
| `<ORDERCODE>_ab1-files/<OSID>_<SAMPLE>.ab1` | "Chromatogram showing relative abundance of A, T, G, and C among raw reads aligned to the consensus at each sequence position" |
| `<ORDERCODE>_coverage-plots/….png` | relative coverage; "Large coverage gaps or abrupt changes may indicate an assembly issue or the presence of multiple plasmid species" |
| `<ORDERCODE>_histograms/….png` | the read-length histogram |
| `<ORDERCODE>_summary_files/….txt` | "consensus length, average coverage, relative molar and mass composition, total reads, total bases, and estimated E. coli genomic DNA contamination" |
| `<ORDERCODE>_interactive-map/….html` | an annotated map |
| `<ORDERCODE>_gel.png` | "Visualization of raw-read lengths across all samples" |
| `<ORDERCODE>.fastq.zip` | "Basecalled raw reads that align to the consensus sequence, including Phred quality scores. Reads that do not align to the consensus are excluded" |
| `<ORDERCODE>-summary-report.csv` | "Summary of sample processing status, consensus length, and estimated E. coli genomic DNA contamination" |

An archived 2024 *Results Interpretation Guide*
(<http://web.archive.org/web/20240617013345/https://plasmidsaurus.com/results_interpretation>)
also listed a `SAMPLE_multimer_analysis.txt`; v1.2 does not. **The file list changes between
versions**, so a reader keys on the folder suffix and the file extension, not on a fixed list.

### 3.2 Azenta/GENEWIZ Plasmid-EZ

From *Plasmid-EZ Quick Start Guide*, `11035-M&G-2 1025`, file dated 10-29-25
(<https://web.genewiz.com/hubfs/2023-04%20GEN%20NGS%20-%20Plasmid-EZ%20Launch/Quick%20Start%20Guide/11035-M%26G-2%2010-29-25%20Plasmid-EZ%20Quick%20Start%20Guide.pdf>),
and the Plasmid-EZ FAQ (<https://web.genewiz.com/faqs/ngs-plasmid-ez>).

| File | What it holds |
| --- | --- |
| `30-xxxxxxxxx_report.html`, `sample_reports/` | read counts, assembly status, contigs, a virtual gel |
| assembly: FASTA, FASTQ and `.ab1` | the consensus |
| annotation: GenBank, CSV | features |
| `{sample_name}_base_count.tsv` | "the coverage score by base positions … also … the spread of base calling in the assembled contig" |
| `{sample_name}_vc.tsv`, a VCF and its index | variant calls: "Any base position with more than 5% variability compared to the assembled contig will be highlighted" |
| `vs_ref.var.txt` | differences against a reference, where the customer gave one |
| raw reads | FASTQ |

Topology is reported per sample as Circular, Linear, Multiple circular, Mixed ("linear and
circular contigs assembled"), Multiple linear or Failed ("no assembly"). Mixed here is a topology,
not a mixed population. The guide's `.ab1` is "Sanger-like data": "hovering over each position
provides the called base and quality score", and "Mixed base positions identified in the ab1 file
can be further reviewed in the variation Excel file … to evaluate exact read counts".

### 3.3 Eurofins Genomics

From *Data Deliverables*
(<https://eurofinsgenomics.com/en/products/nanopore-sequencing/data-deliverables/>, no version
shown): an HTML analysis report; a NanoPlot summary; the consensus as FASTA and GenBank; an
annotated map; a read-length histogram; a coverage plot; a per-base table of "position, base,
matches, mismatches, insertions, deletions"; and raw reads that "only contain the sequences that
align with the consensus".

### 3.4 What makes a base trustworthy

| Claim | Owner | Words |
| --- | --- | --- |
| Coverage that supports an accurate consensus | Plasmidsaurus v1.2 | "coverage greater than approximately 20× generally supports a highly accurate consensus sequence" |
| The same | Eurofins *Results Interpretation Guide* | "coverage of approximately 20x or higher suggests a highly accurate consensus" |
| A failed sample | GENEWIZ FAQ | "failure of the sample to produce consensus sequence with 10x coverage or less" |
| A failed sample | Plasmidsaurus v1.2 | "considered unsuccessful when the sequencing data are insufficient to generate a consensus sequence" |
| Consensus accuracy | Plasmidsaurus v1.2 | "often greater than Q60, corresponding to an estimated accuracy of 99.9999%"; 2024 guide: "typically >99.99%" |
| Read accuracy | Eurofins product page, quoting ONT | "Oxford Nanopore's reports the raw read accuracy is 98.3% and the consensus accuracy for SNPs is 99.6% with 50X coverage"; its own guide says "raw read accuracy exceeds 99%" |
| Read counts | Plasmidsaurus v1.2 | "read counts ranging from tens to thousands" |
| Known error modes | Plasmidsaurus FAQ (<https://plasmidsaurus.com/faq/plasmid>, undated) | "methylation sites and homopolymers with a length >9" |
| The methylation motifs | Plasmidsaurus technical note, Sep 12, 2025 (<https://plasmidsaurus.com/technical-note/oxford-nanopore-specific-error-modes>) | "Dam motif Gm6ATC and the Dcm motif … C5mCTGG or C5mCAGG"; a homopolymer ">9" is "often truncated by a base or two"; any other mismatch is "almost certainly real" |
| The same, per vendor | GENEWIZ guide | "GATC, CCAGG and CCTGG sites are methylation sites in E. coli. These sites may show discrepancy"; at a methyl GATC "the variation file will show half A/half G"; "Sanger sequencing follow-up may be needed" |
| The same | Eurofins guide | "Dcm methylation sites (CC[A/T]GG)" and homopolymer deletions are "the most common error modes" |

**The earlier caution is settled.** `domestication-methods.md` flagged that its Plasmidsaurus and
Eurofins figures came back in overlapping wording. Read again at each owner: Q60 is
Plasmidsaurus's own (FAQ and v1.2), and 98.3% / 99.6% at 50X is on Eurofins' product page, quoted
from Oxford Nanopore rather than measured by Eurofins. Neither figure is on the other's pages.

### 3.5 What a mixture looks like

| Owner | Words |
| --- | --- |
| Plasmidsaurus FAQ | "This service is intended for a clonal population of molecules. If your species are very similar … the pipeline will most likely create a single consensus file, with mixed peaks observed in the .ab1 file where there are SNPs and indels. If your species are sufficiently distinct … a single consensus sequence for the molecular species that produces the largest amounts of total sequencing data" |
| Plasmidsaurus FAQ | multimers "are not considered different molecular species … you will only receive the monomer consensus"; "Sequencing is considered successful if the pipeline is able to generate a consensus, even if it is not your target" |
| Plasmidsaurus histogram note, Oct 8, 2025 (<https://plasmidsaurus.com/technical-note/how-to-interpret-plasmid-read-length-histograms>) | the histogram colours reads mapping to the consensus, to *E. coli*, and unmapped, which is how a second species shows |
| GENEWIZ FAQ | "This service is intended for a clonal population of plasmids. If plasmid pools or mixtures are submitted, the sequencing and analysis results cannot be guaranteed" |
| GENEWIZ guide | the 5% highlight in `_vc.tsv`; a worked example of 144 G against 46 A reported as "22% variation" |

So **a whole-plasmid consensus can be confidently wrong**: a second species that out-reads the
target becomes the consensus and the sample still counts as successful. The consensus alone does
not show a mixture. The per-base table, the `.ab1` and the histogram do.

### 3.6 The workflow some results resemble

Oxford Nanopore's open workflow `wf-clone-validation` (v1.8.4,
<https://github.com/epi2me-labs/wf-clone-validation>) writes a consensus with a per-base
quality, `{alias}.final.fastq`, "Sequence and quality score of the final assembly"; a
`sample_status.txt` whose values include "Completed successfully" and "Failed due to insufficient
reads"; annotations; and, given a reference, a BAM and variant calls. Its README says "The
workflow has no way of reporting contaminants". No vendor read here says it runs this workflow,
so it is evidence of a file shape, not of what any vendor sends. How it aligns and judges is
`sequencing-read-packages.md`'s.

## 4. A DMX well

What `synbio/dmx/method.py` hands over today, read from the code on `main`:

| | barcode ligation | index PCR |
| --- | --- | --- |
| Input to `judge_well` | `reads: int`, `called: Sequence[str]` | the same |
| Depth wanted | above 150 (Qian SI, "Generation of consensus sequences") | above 20, 10 tolerated (LevSeq SI) |
| How a consensus base is called | at ≥51% read support (Qian SI, per `route-choice.md` section 4) | not stated in the code |
| What the pass is | `identity_check`: an exact match, upper-cased, across the whole designed region; more than one consensus is mixed and fails | the same |

The code's own comment on `ROUTE_INDEX_PCR` names the gap: LevSeq's second criterion, a mean
error below 10%, "is a mean over per-position base counts, and a well reaches this package as a
read count and a consensus call". **So a DMX well hands over no per-base support.** It is the
whole-plasmid shape with the per-base table removed: a called sequence, how many reads stood
behind it, and nothing per position. Whether that changes is #551's "DMX's call through it", not
this note's.

## 5. What confirms a clone

| Owner | What it reads | Words |
| --- | --- | --- |
| Addgene, J. Taylor-Parker, "Addgene's Tips for Plasmid Quality Control", 2016-01-14 (<https://blog.addgene.org/addgenes-tips-for-plasmid-quality-control>) | Sanger, before 2017 | "All of our incoming plasmids are sequenced at least twice"; "We then select forward and reverse primers for each sample"; "Important features may include cloning junctions (where most sequence assembly errors occur)"; "If the insert is large, we may only check the 5’ and 3’ ends to verify that the sequence assembly and identity of the insert are correct" |
| Addgene, A. Hazen, 2017-08-03, updated 2021-04-06 (<https://blog.addgene.org/addgene-moves-to-ngs-verification-powered-by-seqwell>) | whole plasmid, MiSeq "2x251" | "we will often find a few mismatches in the origin of replication or other common backbone elements"; "these few minor mismatches usually don't affect the function of the plasmid" |
| Addgene, A. Shepard, "A Look at Addgene's QC Process", 2025-05-06, with a footnote dated March 2026 (<https://blog.addgene.org/a-look-at-addgenes-qc-process>) | whole plasmid, now mostly long-read | "a QC issue is defined as a discrepancy that we believe may affect the function of the plasmid"; otherwise a note that "it is not known to impact the function of the plasmid"; "the QC team finds errors in 30% of deposited plasmids!" |
| Twist Bioscience, *Clonal Genes* datasheet, DOC-001424 REV9, 2026-01 (<https://www.twistbioscience.com/content/dam/twistbioscience/resources/2026-01/DOC-001424%20Clonal%20Genes%20REV9%20singles.pdf>) | NGS | "Clonal Genes are made by cloning an insert (synthesized sequence) into a vector (non-synthesized sequence). After synthesis and cloning, we use NGS to verify that the insert is 100% sequence perfect"; its example shows "sufficient read depth across the entire plasmid" |
| Bio-Synthesis FAQ (<https://biosyn.com/faq/assurance-of-correct-gene-sequence.aspx>, undated) | Sanger | clones are verified "to at least single-strand depth across the synthesized insert" |
| This repo: `bench/steps.py`, the sequencing step a cloning method's protocol ends with | Sanger | "Check the read across every junction and the whole of each insert." |

**What the owners agree on.** A confirmed clone has **every junction read and the whole insert
covered**, and the read agrees with the expected sequence there. One strand is enough for a clone
(Bio-Synthesis's "at least single-strand depth"). Addgene reads two reactions from opposite sides,
but for reach, not for a second strand at every base. The one relaxation in print is Addgene's
for a large insert by Sanger, which reads only its two ends. Twist holds the insert to "100%
sequence perfect" and says nothing of the vector.

**No owner sets a standard.** This survey found no published rule for read depth, for both
strands, or for how far past a junction a read must run. "Sequence-verified" is each owner's
practice, stated in its own words.

### 5.1 A disagreement outside the insert

Three readings, each with an owner:

1. **The parent plasmid was never what its map says.** Addgene finds errors in 30% of deposited
   plasmids, and "a few mismatches in the origin of replication or other common backbone elements"
   often enough to say so. Bai et al. 2025, *Nucleic Acids Res.* 53(14):gkaf697,
   doi:10.1093/nar/gkaf697, PMID 40794862 (authors at VectorBuilder), abstract: of the plasmids
   they received "nearly half of them contained design and/or sequence errors", and of AAV
   transfer plasmids "about 40% carried mutations in the inverted terminal repeat (ITR) regions".
   The abstract does not split backbone from insert. So an expected product built from a vector's
   map inherits the map's errors. `vector-qc-panel.md` section 5 already asks for the parent to
   be sequenced beside the prep; a backbone disagreement the parent shares was not made by the
   cloning.
2. **A sequencing artefact.** On a nanopore consensus, a disagreement at a Dam or Dcm motif or in
   a homopolymer longer than 9 is a known error mode, which Plasmidsaurus marks "Likely Match";
   any other mismatch is "almost certainly real" (section 3.4). GENEWIZ sends such a position to
   Sanger.
3. **A real change that does not matter.** Addgene's rule is functional: a discrepancy is an
   issue when "we believe [it] may affect the function of the plasmid", and otherwise a note.

So no owner fails a clone for a disagreement outside the insert as such. Each judges it by the
feature it falls in, and first asks whether the sequencer or the reference made it.

## 6. Example files for a test suite

| File | Where | Size | Licence | What it held, opened here |
| --- | --- | --- | --- | --- |
| `3730.ab1` | Biopython, `Tests/Abi/` (<https://github.com/biopython/biopython/tree/master/Tests/Abi>) | 299,987 bytes | Biopython License Agreement (`LICENSE.rst`) | a 3730xl KB-called read with qualities and 7 IUPAC mixed calls at its ends; the calls read mCherry then EGFP (section 2.2) |
| `310.ab1` | the same | 222,099 bytes | the same | a 310 read with every quality zero and `N` calls: the no-quality case |
| `3100.ab1`, `A6_1-DB3.ab1`, `no_smpl1.ab1`, `nonascii_encoding.ab1`, `empty.ab1`, `fake.ab1` | the same | 209,224; 261,660; 253,428; 295,984; 261,522; 26 bytes | the same | not opened; the names say what each tests |
| `per_base_df.csv` | `odcambc/dimple-qc-app` (<https://github.com/odcambc/dimple-qc-app>) | 192,143 bytes | MIT | a **processed** copy of a Plasmidsaurus per-base table: the columns `pos, ref, reads_all, matches, mismatches, deletions, insertions, A, C, T, G, low_conf, homopolymer, methylation` as the app reads them, then columns the app adds; `pos` starts at 1 |
| `wf-clone-validation-demo.tar.gz` | Oxford Nanopore (<https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-clone-validation/wf-clone-validation-demo.tar.gz>) | 23,389,102 bytes | Oxford Nanopore Public License 1.0, research purposes only | raw reads for 11 barcodes and ~3 kb references: inputs to a consensus, not a result. Too large, and the wrong stage |
| Addgene sequencing results | each plasmid's `/sequences/` page | — | Addgene's terms | Sanger and full-plasmid results exist per plasmid, but downloading one needs a login, so CI cannot fetch them |

**Sanger: public files are enough to test a reader.** Biopython's traces are small, carry both
the quality and the no-quality case, and are under a licence that allows copying with the notice.
**None comes with the construct it was read from**, so a test of a junction verdict builds its
expected record around the read's own trimmed calls and plants the disagreement it wants.

**Tracked.** `3730.ab1` and `310.ab1` are copied unchanged into `tests/data/` from Biopython at
commit `e38c64daa2` (<https://github.com/biopython/biopython/tree/e38c64daa2/Tests/Abi>), fetched
2026-10-10. Biopython's `LICENSE.rst` from the same commit sits beside them as
`tests/data/biopython-LICENSE.rst`, the notice its licence asks a copy to carry. SHA-256:

| File | SHA-256 |
| --- | --- |
| `3730.ab1` | `e4663e4db40232576ccdda5b878dddb01ef80a3d1b032941ba053146ce53f77b` |
| `310.ab1` | `0a395d39520d24a89af4c77dd55fd7fb632e732a6da807ceb6ead0c2465cb887` |

**Whole-plasmid: no public raw vendor result exists.** No vendor publishes a demo folder, and the
one public per-base table has been processed. The lab supplies its own result — one Plasmidsaurus
folder for a plasmid it already holds — or a test writes a small file in the vendor's documented
shape. The header of Plasmidsaurus's per-base table is unpublished, so the lab's own file is the
only source for it.

## 7. Holes

- **CRL has no published definition** that this survey reached; GENEWIZ reports it per read.
- **Ewing and Green 1998 were not read in full** (*Genome Res.* 8:175–185 and 186–194; the
  publisher's pages asked for a login). The quality formula is cited from Richterich 1998 and
  Applied Biosystems instead. phrap.org, which documents phred's trimming, refused the
  connection, so **Mott's 0.05 rests on Biopython's code alone**.
- **How many bases at the start of a read are unreadable has no single number**: 20–40, 25–35,
  and blobs to 120, from three owners.
- **The KB basecaller's default mixed-base threshold** is not stated; only the recommended
  15–25%.
- **Plasmidsaurus publishes no column header for its per-base table**, and no vendor says whether
  its positions start at 0 or 1. The processed copy starts at 1. ADR 0001 says such a file
  converts at its own boundary.
- **No vendor defines a per-base confidence**, the scale of the qualities in a consensus FASTQ,
  or the smallest mixture it detects. GENEWIZ's 5% highlight is the only number.
- **Whether GENEWIZ's Plasmid-EZ runs on Oxford Nanopore** is not stated on its pages.
- **No standard defines a sequence-verified clone**: no depth, no strand rule, no distance past a
  junction.
- **Bai et al. 2025 was read only to its abstract**; figures from its body are left out.
- **NEB's sequencing-analysis page and an Addgene help article** returned 403.
- **A DMX well hands over no per-base support** (section 4), so nothing here judges LevSeq's mean
  error criterion.

### 7.1 The sibling note's three holes

`sequencing-read-packages.md` (#552) section 6 left three numbers without a source. What this
survey found for each:

| Its hole | What an owner states | Where here |
| --- | --- | --- |
| A margin between the two strands' scores | nothing | — |
| A mixed-peak ratio | Sanger: a "2nd highest peak threshold for mixed base identification", "recommended range is 15 to 25%" (Applied Biosystems Chemistry Guide). The 0.33 the tools use sits outside that range. Whole plasmid: GENEWIZ highlights a position at "more than 5% variability" | sections 2.5, 3.2 |
| An identity or coverage threshold | identity: the insert "100% sequence perfect" (Twist); coverage: every junction and the whole insert, "at least single-strand depth" (Bio-Synthesis, Addgene); whole-plasmid depth: "approximately 20×" supports an accurate consensus (Plasmidsaurus, Eurofins), "10x coverage or less" is a failure (GENEWIZ). None sets a share of the backbone, so wf-clone-validation's 99% identity and 95% coverage still have no source | sections 3.4, 5 |

## 8. Sources

Read 2026-10-09 and 2026-10-10 unless a date is given. URLs, editions and quotes are inline where
each is used; this list adds what is not.

| Source | Edition | How read |
| --- | --- | --- |
| Applied Biosystems, *Applied Biosystems Genetic Analysis Data File Format* | July 2006; © 2006 Applied Biosystems | PDF from a mirror at Moscow State University (<https://kodomo.fbb.msu.ru/~lewis/terms/term3/pr6/ABIF_File_Format.pdf>); the Applied Biosystems address other readers cite, `www6.appliedbiosystems.com/support/software_community/ABIF_File_Format.pdf`, was not tried |
| Biopython | 1.88, the version this repo pins | `Bio/SeqIO/AbiIO.py` in the pixi environment; test files from `master` |
| Plasmidsaurus *Technical Documentation: Plasmid* | PS-0002-E v1.2, revised Sept. 22, 2026 | page text |
| Plasmidsaurus FAQ, and two technical notes | FAQ undated; notes Sep 12, 2025 and Oct 8, 2025 | page text |
| Plasmidsaurus *Results Interpretation Guide* | archived 2024-06-17 | Wayback Machine |
| Eurofins Genomics: *Whole Plasmid Sequencing* (<https://eurofinsgenomics.com/en/products/nanopore-sequencing/whole-plasmid-sequencing>), *Data Deliverables*, *Results Interpretation Guide* (<https://eurofinsgenomics.com/en/products/nanopore-sequencing/results-interpretation-guide-for-nanopore-sequencing/>) | no version shown | page text |
| Azenta/GENEWIZ *Plasmid-EZ Quick Start Guide* and FAQ | 11035-M&G-2 1025 (10-29-25); FAQ undated | PDF and page text |
| Oxford Nanopore, `wf-clone-validation` (<https://github.com/epi2me-labs/wf-clone-validation>) | v1.8.4 | README, CHANGELOG and source |
| MGH DNA Core, sequencing troubleshooting (<https://dnacore.mgh.harvard.edu/new-cgi-bin/site/pages/sequencing_pages/seq_troubleshooting.jsp>) | undated | page text; agrees with section 2.5: "There is more than one template present in the sequence reaction" |
| Richterich, P. *Genome Res.* 8:251–259 (1998), doi:10.1101/gr.8.3.251 | PMC310698 | full text |

Text copies of the ABIF specification and of the whole-plasmid vendor pages are kept under
`reference_docs/sequencing-read-evidence/`.

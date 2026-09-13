---
search:
  exclude: true
---

# Codon usage tables: sources, licences and what the shipped table is

Research note for issue #10. Everything below was read on **2026-09-12**. It records where
`src/liulab_mbio/data/codon_usage.json` comes from, which published tables were rejected and
why, and what the shipped table is and is not.

## 1. The verdict, first

**The shipped table is counted from a genome, not copied from a codon usage database.**

Every published compilation that was checked either states no licence at all, reserves its
rights, or carries a non-commercial clause that an MIT package cannot take. The primary
sequence data does not have that problem: NCBI states it places no restrictions on the
distribution of the sequence records, and a count computed over those records is this
package's own measurement rather than a redistribution of anyone's table.

That is also the better answer on the merits. The obvious table to reach for — Kazusa's
*Escherichia coli* K-12 entry — turns out to be **14 coding sequences**, and one of its cells
is zero. Counting the genome gives 4,317.

## 2. What the data file holds

| Field | Meaning |
| --- | --- |
| `name` | the short name `codon_usage()` is asked for, `e-coli-k12` |
| `organism`, `taxid`, `accession` | *Escherichia coli* str. K-12 substr. MG1655, 511145, `U00096.3` |
| `cds_count`, `codon_count` | 4,317 coding sequences, 1,342,016 codons |
| `counts` | all 64 codons, stop codons included |
| `note` | that this is a whole-genome table and not a highly expressed set |

A coding sequence is skipped, rather than counted in part, when its length is not a whole
number of codons or when it carries a base outside `ACGT`: a partial or ambiguous record would
bias the table by the codons it does spell. One of the 4,318 records NCBI returns is skipped
for that reason.

**No cell is zero.** The rarest codon in the genome is `TAG` at 306, which matters because a
fraction, a log ratio or a codon adaptation weight computed from a zero is undefined.

## 3. Licences, source by source

### 3.1 NCBI GenBank: no restrictions, and an honest caveat

From NCBI's *Website and Data Usage Policies and Disclaimers*, under **Molecular Data Usage**
(`https://www.ncbi.nlm.nih.gov/home/about/policies/`, read 2026-09-12):

> Databases of molecular data on the NCBI Web site include such examples as nucleotide
> sequences (GenBank), protein sequences, macromolecular structures, molecular variation, gene
> expression, and mapping data. They are designed to provide and encourage access within the
> scientific community to sources of current and comprehensive information. Therefore, NCBI
> itself places no restrictions on the use or distribution of the data contained therein. Nor
> do we accept data when the submitter has requested restrictions on reuse or redistribution.
> However, some submitters of the original data (or the country of origin of such data) may
> claim patent, copyright, or other intellectual property rights in all or a portion of the
> data (that has been submitted). NCBI is not in a position to assess the validity of such
> claims and since there is no transfer of rights from submitters to NCBI, NCBI has no rights
> to transfer to a third party. Therefore, NCBI cannot provide comment or unrestricted
> permission concerning the use, copying, or distribution of the information contained in the
> molecular databases.

Read it as it is written. The first two sentences are strong — no restrictions, and NCBI
refuses submissions that carry any. The last sentence is a disclaimer: NCBI is **not** issuing
a licence. This is an absence of restrictions plus a refusal to accept restricted submissions,
which is the best position available for genome sequence, and it is a much stronger position
than redistributing a third party's compilation would be.

**Verdict: the genome is counted, and the counts ship.**

### 3.2 Kazusa CUTG: no licence statement anywhere, and the wrong table besides

The Codon Usage Database (`https://www.kazusa.or.jp/codon/`) states no licence, no terms of
use and no redistribution permission. The front page, the readme and the institute privacy
page were all read on 2026-09-12. The only normative statement is a citation request, from
`https://www.kazusa.or.jp/codon/readme_codon.html` under the heading **"Please cite"**:

> Codon usage tabulated from the international DNA sequence databases: status for the year
> 2000. Nakamura, Y., Gojobori, T. and Ikemura, T. (2000) Nucl. Acids Res. 28, 292.

The same page gives its input as:

> NCBI-GenBank Flat File Release 160.0 [June 15 2007].

The Integbio catalogue record (`https://catalog.integbio.jp/dbcatalog/en/record/nbdc00033`)
shows no information under "Link(s) to Terms of use"; the CC0 shown there applies to the
catalogue record, not to the database.

**A citation request is not a copyright licence, and silence is not permission.** Kazusa has
granted nothing, so an MIT package cannot rely on it.

**It is also the wrong table.** Its own header for the K-12 entry reads
`Escherichia coli K12 [gbbct]: 14 CDS's (5122 codons)` — fourteen coding sequences, about
0.3% of the genome's codons — and its `UAG` cell is `0.00 0.0 (0)`. Kazusa's other *E. coli*
entry, `species=37762`, counts 8,087 CDS and 2,330,943 codons, which *exceeds* one genome
because CUTG counts redundant and partial GenBank entries rather than one annotated genome.
No Kazusa entry corresponds to MG1655's coding sequences.

**Verdict: not shipped, and not trusted as a cross-check without saying how small it is.**

### 3.3 HIVE-CUTs: the article is CC BY, the database is not licensed

`https://hive.biochemistry.gwu.edu/review/codon` redirects to
`https://dnahive.fda.gov/dna.cgi?cmd=cuts_main` and is reachable (2026-09-12). The paper
(Athey et al. 2017, *BMC Bioinformatics* 18:391) is open access:

> Open Access This article is distributed under the terms of the Creative Commons Attribution
> 4.0 International License (`http://creativecommons.org/licenses/by/4.0/`), which permits
> unrestricted use, distribution, and reproduction in any medium, provided you give
> appropriate credit to the original author(s) and the source, provide a link to the Creative
> Commons license, and indicate if changes were made. The Creative Commons Public Domain
> Dedication waiver (`http://creativecommons.org/publicdomain/zero/1.0/`) applies to the data
> made available in this article, unless otherwise stated.

That CC BY covers the article, and the CC0 sentence is scoped to "the data made available in
this article" — what is printed in the paper, not the contents of the web resource. The
availability statement says only:

> This database is available for use by all, including both academics and non-academics.

That is an access policy, not a licence with defined redistribution terms. The site asserts no
licence over the tables, and no scriptable download endpoint was confirmed: a probe of a
plausible file URL returned the site's HTML shell rather than a data file.

One point in its favour, which is not enough on its own: HIVE-CUTs is FDA-produced and
FDA-hosted, and FDA's website policy
(`https://www.fda.gov/about-fda/about-website/website-policies`) says:

> Unless otherwise noted, the contents of the FDA website (`www.fda.gov`) — both text and
> graphics — are not copyrighted. They are in the public domain and may be republished,
> reprinted and otherwise used freely by anyone without the need to obtain permission from FDA.

Note the scope: that sentence names `www.fda.gov`, and HIVE-CUTs is served from
`dnahive.fda.gov`. Suggestive, not dispositive.

**Verdict: cite it, do not vendor it.**

### 3.4 GenScript: all rights reserved

The codon frequency table page (`https://www.genscript.com/tools/codon-frequency-table`)
carries `© 2002-2026 GenScript All rights reserved.` Its "Terms of Use and Privacy" link goes
to a privacy policy with no content-reuse clause, and the Standard Terms and Conditions behind
that are a services agreement about deliverables, not about website content. There is an
express reservation, no permission, and no stated provenance for the table.

**Verdict: not copied.** This is the weakest source of the set — worse than Kazusa, which at
least discloses where its numbers came from.

### 3.5 The rest

| Resource | What it says | May this package ship it? |
| --- | --- | --- |
| Codon Statistics Database (`codonstatsdb.unr.edu`) | article is CC BY-**NC**; the download page asserts no licence over the data | **No** — non-commercial clashes with MIT |
| CoCoPUTs (Alexaki et al. 2019, *JMB*) | not open access; no licence found on the article or the site | **No** — nothing to rely on |
| EMBOSS `Eecoli.cut`, `Eecoli_high.cut` | ship inside a **GPL** tarball, and their own headers name `CUTG146` and `TranstermMay1994` as their input | **No** — copyleft, and the provenance leads back to CUTG anyway |
| `python_codon_tables` | the package is **CC0 1.0**, but its README says "All tables are from kazusa.or.jp" | **No** — a downstream grant cannot convey rights the upstream never gave. Useful as precedent that others do this; not authority that it is allowed |

## 4. What the shipped table is, and what it is not

A **whole-genome** table counts every coding sequence and describes the background the genome
itself uses: mutational bias plus weak selection. A **highly expressed** table counts a small
reference set — ribosomal proteins, elongation factors, glycolytic enzymes — where
translational selection is strong and the bias is much sharper. They are different objects,
and a reader will assume whichever they need, so the data file says which it is in its `note`.

The distinction goes back to Ikemura's correlation of codon choice with tRNA abundance
(*J. Mol. Biol.* 146:1-21 and 151:389-409, 1981) and to the codon adaptation index of Sharp
and Li (*Nucleic Acids Res.* 15:1281-1295, 1987), whose abstract states:

> A simple, effective measure of synonymous codon usage bias, the Codon Adaptation Index, is
> detailed. The index uses a reference set of highly expressed genes from a species to assess
> the relative merits of each codon, and a score for a gene is calculated from the frequency of
> use of all codons in that gene.

EMBOSS quantifies the gap between the two for *E. coli*: its whole-genome file records
`#CdsCount: 5045` with GC3 55.80%, its highly expressed file `#CdsCount: 250` with GC3 57.41%,
and `GAA` at 0.76 there against roughly 0.68 genome-wide.

**The shipped table is the whole-genome one.** It is the right table for domestication, which
is what issue #10 needs it for: a synonymous swap should leave a codon the host uses readily,
and the genomic background answers that. A highly expressed set is what a codon adaptation
index would want, and shipping one is a separate decision with its own reference gene list.

Amino acids are grouped with the **bacterial** genetic code (NCBI table 11). It assigns the
same amino acids as the standard table and differs only in which codons may start a gene,
which is not something a codon usage table is asked about.

## 5. Rebuilding the file

```sh
pixi run python scripts/build_codon_usage.py                 # downloads from NCBI
pixi run python scripts/build_codon_usage.py --cds FILE      # from a local fasta_cds_na file
```

The script fetches the spliced coding sequences of one accession from NCBI's E-utilities:

```text
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
  ?db=nuccore&id=U00096.3&rettype=fasta_cds_na&retmode=text&tool=liulab-mbio
```

`fasta_cds_na` is used rather than the full GenBank record because NCBI has already spliced
the joins and reverse-complemented the coding sequences on the minus strand, which removes a
whole class of parsing error from the build.

Anonymous access works: no API key, no cookies, no browser. NCBI's usage guidance
(`https://www.ncbi.nlm.nih.gov/books/NBK25497/`) is:

> In order not to overload the E-utility servers, NCBI recommends that users post no more than
> three URL requests per second and limit large jobs to either weekends or between 9:00 PM and
> 5:00 AM Eastern time during weekdays. Failure to comply with this policy may result in an IP
> address being blocked from accessing NCBI.

and, from the policies page:

> Use the URL parameter email, and tool for distributed software, so that we can track your
> project and contact you if there is a problem.

The script sends `tool=liulab-mbio` and **no** address: a build script should not put whoever
runs it on a mailing list they did not ask for. One fetch per rebuild is far inside the rate
limit.

Tests never reach the network. `tests/scripts/test_build_codon_usage.py` counts a four-record excerpt.

**`U00096.3` is the accession, not `NC_000913.3`.** Both fetch identically and both report
4,318 coding sequences; U00096.3 is the INSDC record of record and RefSeq is derived from it.

## 6. Cross-checks

- **Two independent counts agree.** The table was counted twice from the same fetch by two
  separate pieces of code written without sight of each other, and both report 4,317 usable
  coding sequences and 1,342,016 codons, with the same per-codon numbers. Stop codons also
  sanity-check: `TAA` 2,763 + `TGA` 1,253 + `TAG` 306 = 4,322, or about one stop per coding
  sequence.
- **Against Kazusa's K-12 entry** (14 CDS, 5,122 codons), the two agree on the commonest
  codons and diverge widely elsewhere, which is the sampling noise five thousand codons buys:

  | Codon | Kazusa, 14 CDS | This table, 4,317 CDS | Difference |
  | --- | --- | --- | --- |
  | `AAA` | 33.2 | 33.7 | +1% |
  | `GAA` | 43.7 | 39.7 | -9% |
  | `CTG` | 46.9 | 53.1 | +13% |
  | `CAC` | 13.1 | 9.7 | -26% |
  | `AGG` | 1.6 | 1.1 | -32% |
  | `TGG` | 10.7 | 15.2 | +42% |
  | `TCT` | 5.7 | 8.4 | +48% |
  | `TAG` | **0.0, count 0** | 0.23, count 306 | undefined |

  Per thousand codons. The last row is the one that settles it: a table with a zero cell cannot
  be used to choose a codon.

## Sources

All read on 2026-09-12.

- NCBI, *Website and Data Usage Policies and Disclaimers*:
  [ncbi.nlm.nih.gov](https://www.ncbi.nlm.nih.gov/home/about/policies/)
- NCBI, *A General Introduction to the E-utilities* (usage guidelines and rate limits):
  [NBK25497](https://www.ncbi.nlm.nih.gov/books/NBK25497/)
- GenBank `U00096.3`, *Escherichia coli* str. K-12 substr. MG1655, complete genome
  (Blattner, F.R. et al. (1997) *Science* 277, 1453-1462), and RefSeq `NC_000913.3`
- Kazusa, *Codon Usage Database* and its readme:
  [kazusa.or.jp](https://www.kazusa.or.jp/codon/),
  [readme_codon.html](https://www.kazusa.or.jp/codon/readme_codon.html); K-12 table
  `cgi-bin/showcodon.cgi?species=83333`
- Nakamura, Y., Gojobori, T. and Ikemura, T. (2000) Codon usage tabulated from the
  international DNA sequence databases: status for the year 2000. *Nucleic Acids Res.* 28, 292
- Athey, J. et al. (2017) A new and updated resource for codon usage tables.
  *BMC Bioinformatics* 18, 391.
  [doi:10.1186/s12859-017-1793-7](https://doi.org/10.1186/s12859-017-1793-7) (CC BY 4.0), and
  the HIVE-CUTs resource at [dnahive.fda.gov](https://dnahive.fda.gov/dna.cgi?cmd=cuts_main)
- FDA, *Website Policies*:
  [fda.gov](https://www.fda.gov/about-fda/about-website/website-policies)
- GenScript, *Codon Frequency Table*:
  [genscript.com](https://www.genscript.com/tools/codon-frequency-table)
- Subramanian, A. et al. (2022) Codon Statistics Database, at `codonstatsdb.unr.edu` (CC BY-NC)
- Alexaki, A. et al. (2019) Codon and Codon-Pair Usage Tables (CoCoPUTs). *J. Mol. Biol.* 431,
  2434-2441
- EMBOSS, *Licence*: [emboss.sourceforge.net](https://emboss.sourceforge.net/licence/), and its
  `Eecoli.cut` and `Eecoli_high.cut` data files
- `python_codon_tables`, CC0 1.0, whose README names kazusa.or.jp as its input
- Ikemura, T. (1981) *J. Mol. Biol.* 146, 1-21 and 151, 389-409
- Sharp, P.M. and Li, W.H. (1987) The codon adaptation index. *Nucleic Acids Res.* 15,
  1281-1295
- `docs/research/restriction-enzyme-data.md` (issue #7) for the precedent this note follows:
  ship only what the licence allows, cite and re-enter the rest, and record what is unverified

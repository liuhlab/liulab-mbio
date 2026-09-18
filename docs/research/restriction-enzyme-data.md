---
search:
  exclude: true
---

# Restriction enzyme data: sources, licences and conventions

Research note for issue #7. Everything below was retrieved on **2026-09-12** unless a section
stamps its own date. It records where
`src/liulab_mbio/data/enzymes.json` comes from, what each source allows, and which values are
still unverified.

## 1. What the data file holds, and who supplied each field

| Field | Source |
| --- | --- |
| `site`, `top_cut`, `bottom_cut`, `isoschizomers`, `suppliers` | REBASE `withrefm` |
| `overhang_length`, `end`, `type` | derived from the cut offsets by the build script |
| `commercial_name`, `catalog_number`, `supplier` | supplier product page |
| `supplied_buffer` | supplier product page, under "Reaction Conditions" |
| `incubation_celsius`, `heat_inactivation_celsius`, `heat_inactivation_minutes` | supplier product page |
| `methylation` (Dam, Dcm, CpG) | supplier product page |

Every enzyme carries a `provenance` map naming the source of each value, and an `unverified`
list naming fields no source stated. The sources themselves are listed in the file's `sources`
table with their URL and the date read: `archived` is a snapshot, `retrieved` a live page. A
source carrying both was read twice, and section 6 records that the two readings agree.

## 2. Licences

### REBASE: open, with a citation

REBASE states its own terms on its citation page:

> REBASE is an open database provided to benefit the scientific community. Adaptation or
> distribution of information from REBASE should also remain open and accessible to the
> public. If you include or otherwise utilize REBASE information in or for your work, please
> cite the latest review.

**Verdict: REBASE-derived values ship.** This repository is open and MIT-licensed, which is
what "remain open and accessible" asks for, and the citation is carried in the data file and
here:

> Roberts, R.J., Vincze, T., Posfai, J., Macelis, D. REBASE: a database for DNA restriction
> and modification: enzymes, genes and genomes. *Nucleic Acids Res.* 51: D629–D630 (2023).
> [doi:10.1093/nar/gkac975](https://doi.org/10.1093/nar/gkac975)

The file header adds `Copyright (c) Dr. Richard J. Roberts, 2024. All rights reserved.` The
copyright notice and the open-database statement sit together, so the citation page is read as
the grant and the notice as the reservation of authorship. This closes the REBASE gap issue #4
left open.

### NEB: cite, never mirror

NEB's terms of use, on the `Copyright and Use` heading:

> As between you and us, all Information contained on the Sites is copyright of NEB unless
> specifically identified otherwise. You may not use or copy the Content other than for
> non-commercial individual reference purposes with all copyright or other proprietary notices
> retained. Except as expressly provided above … you may not otherwise copy, display, download,
> install, use, compile, distribute, mirror on another server, transfer, license, lease, loan,
> sell, modify, reproduce, republish, retransmit … or create any derivative work based on such
> Information.

**Verdict: no NEB file, page or table is redistributed.** What ships is a handful of facts per
enzyme — a temperature, a catalogue number, a sensitivity word — each read from a cited page
and re-entered by hand in `scripts/enzyme_properties.toml`. This is the same line issue #4 drew
("cite and paraphrase numbers, never redistribute the files"). It is also why the build script
fetches REBASE and nothing else: an automated merge of an NEB dataset would be mirroring it.

The same reading applies to `nc3.neb.com/NEBcutter/data/enzymes.json`, NEBcutter 3's enzyme
table, which is reachable without a 403. It was read to confirm product names and reaction
temperatures. It is **not** a build input and none of it is copied wholesale.

`supplied_buffer` was added for issue #134 under this same rule. Every product page prints the
buffer under "Reaction Conditions", beside the incubation temperature already re-entered from
it, so one word per product is read and typed in. NEB publishes the same values twice more, as
the *Performance Chart for Restriction Enzymes* and in NEBcutter's table; both are compiled
tables, neither is republished here and neither is a build input.
`docs/research/restriction-ligation.md` §4 argues that line for issue #120 and reaches the same
verdict.

### Thermo Fisher

Thermo's product pages were read for BpiI and AarI, the two enzymes in the set NEB does not
sell. The same rule applies: cited, re-entered, not mirrored.

## 3. Cut offsets: the convention, and how REBASE's notation maps onto it

A cut offset is a **0-based boundary counted from the first base of the recognition site as
written 5' to 3'**, like a Python slice point: offset `k` severs a strand between site bases
`k - 1` and `k`. `top_cut` is the cut in the strand carrying the site; `bottom_cut` is the cut
in the complementary strand, written in that same frame, so it says the bottom strand is cut
between the partners of bases `k - 1` and `k`. This matches `docs/adr/0001-coordinates.md`:
positions are 0-based, and a boundary is a slice point rather than a base.

An offset outside `0 .. len(site)` is a cut outside the site, which is what makes an enzyme
Type IIS. The build script derives `type` that way, rather than reading REBASE's subtype field.

REBASE writes cuts two ways, and both are read by `parse_site`:

| REBASE | Means | `top_cut` | `bottom_cut` |
| --- | --- | --- | --- |
| `G^AATTC` (EcoRI) | cut marked inside the site, same position on both strands | 1 | 5 |
| `GGTAC^C` (KpnI) | as above, leaving a 3' overhang | 5 | 1 |
| `CCC^GGG` (SmaI) | as above, blunt | 3 | 3 |
| `GGTCTC(1/5)` (BsaI) | offsets past the 3' end of the site | 7 | 11 |
| `GCTCTTC(1/4)` (SapI) | as above, three-base overhang | 8 | 11 |
| `AACTGG(-5/-1)` | negative offsets count back from the site's 3' end | 1 | 5 |

The two worked examples issue #7 asks for:

- **BsaI `GGTCTC(1/5)`**: site length 6, so `top_cut = 6 + 1 = 7` and `bottom_cut = 6 + 5 = 11`.
  On `GGTCTCN^NNNN`, the top strand is cut after one spacer base and the bottom strand four
  bases further on, leaving a four-base 5' overhang.
- **SapI `GCTCTTC(1/4)`**: site length 7, so `top_cut = 8` and `bottom_cut = 11`, a three-base
  5' overhang.

`end` is `5'` when `bottom_cut > top_cut`, `3'` when it is smaller, and `blunt` when they are
equal; `overhang_length` is the distance between them.

A site found on the reverse strand is handled by `Enzyme.cut_positions`, which mirrors both
offsets about the end of the matched span, so the returned pair is always (top strand, bottom
strand) in the record's own coordinates.

REBASE also writes a few enzymes with **two** pairs of offsets, one on each side of the site
(BcgI, `(10/12)GCANNNNNNTCG(12/10)`). None is in this set, and `parse_site` refuses them rather
than keeping half a record.

## 4. Isoschizomers: the rule used here

REBASE's `<2>` field lists every enzyme of the same **specificity**, which includes
neoschizomers — the same site cut in a different place — and hundreds of enzymes nobody sells.
`isoschizomers` in the data file is narrower: an enzyme is an isoschizomer here when it shares
both the site and both cut offsets, **and** REBASE lists at least one supplier for it.

The difference is not cosmetic. SmaI (`CCC^GGG`, blunt) and XmaI (`C^CCGGG`, four-base 5'
overhang) appear in each other's REBASE list; treating them as interchangeable would give a
blunt vector a sticky end. Under the rule above SmaI's isoschizomer list is empty, which is the
honest answer for the set we ship.

This also decides lookup. `get_enzyme` resolves a name, then a commercial name, then an
isoschizomer. `BstV2I` is an isoschizomer of both BbsI and BpiI — enzymes with different
suppliers, and different heat-inactivation answers — so it raises rather than picking one.

## 5. Enzymes shipped, and where each product's facts came from

Sixteen multiple-cloning-site enzymes (EcoRI, SacI, KpnI, SmaI, XmaI, BamHI, XbaI, SalI, PstI,
SbfI, SphI, HindIII, NdeI, NcoI, XhoI, NotI) and ten Type IIS enzymes (BsaI, BsmBI, Esp3I,
BbsI, BpiI, SapI, BspQI, PaqCI, AarI, BtgZI). Section 8 adds SrfI and PmeI, which are neither.

NEB returns HTTP 403 to scripts, so every NEB product page was first read through a dated Wayback
snapshot; the snapshot date is in the data file beside each source. The snapshots run from
2025-10-07 to 2026-08-14 — recorded per source rather than as one date, because that is what
each value's age actually is. All 28 pages were read again, live, on **2026-09-18** through the
`r.jina.ai` route `docs/research/restriction-ligation.md` records, which is where
`supplied_buffer` comes from; each source carries that date as `retrieved`.

Values worth flagging:

- **BsmBI-v2 `incubation_celsius` is 55**, the digestion optimum on NEB's product page. NEB's
  Golden Gate programs cycle BsmBI-v2 at 42 °C. That is a protocol temperature and belongs to
  the protocol, not to the enzyme record. BspQI is the same story: 50 °C here, 42 °C in NEB's
  Ligase Master Mix table.
- **BsmBI-v2 Dam and Dcm are both "not sensitive"**, and CpG is "blocked", read from the
  2026-08-14 snapshot. Issue #4 left these two cells unread; this closes that gap.
- **BtgZI**: a 2025-12-31 snapshot of the product page was found under
  `neb.com/en-us/products/r0703-btgzi`, which #4 recorded as having no snapshot. Incubation is
  60 °C, heat inactivation 80 °C for 20 minutes, CpG impaired.
- **KpnI-HF, BamHI-HF and PstI-HF** say `Heat Inactivation: No`. Their
  `heat_inactivation_celsius` is `null` and they are *not* listed as unverified: the supplier
  answered, and the answer was "heat does not stop it".
- **`supplied_buffer` is a property of the product, not of the enzyme name.** Twenty-four of the
  twenty-six NEB products ship in `rCutSmart Buffer`, including all sixteen multiple cloning site
  enzymes; BsmBI-v2 and BspQI ship in `NEBuffer r3.1`. The two Thermo products carry Thermo's own
  answer from the "Compatible Buffer" row of their pages — `Buffer G` for BpiI and the unique
  `Buffer AarI` for AarI — so no pair crossing the two suppliers reads as sharing a buffer.
  PaqCI's page gives `rCutSmart Buffer` and asks for the PaqCI Activator on top of it; the
  activator is not a buffer and `cloning.goldengate.bench` already adds it.
- **The concentration is not part of the name.** The pages write `1X rCutSmart™ Buffer` and
  `10X Buffer G`; the field holds the buffer.

### Unverified values

| Enzyme | Field | Why |
| --- | --- | --- |
| BpiI (ER1011) | `heat_inactivation_celsius`, `heat_inactivation_minutes` | Thermo's page says "Sensitive to Heat Inactivation: Yes" and gives no temperature |
| AarI (ER1581) | `heat_inactivation_celsius`, `heat_inactivation_minutes` | as above |

Thermo's usual 65 °C for 20 minutes is **not** recorded, because the page does not say it.

AarI's CpG value is Thermo's own coarser wording, `"sensitive"`, where NEB says
`"impaired by overlapping"` for PaqCI, the same specificity. Both are kept as their supplier
wrote them rather than normalised into one vocabulary that neither supplier uses.

## 6. Cross-checks

Three independent readings agree with the shipped cut offsets, and none of them is a build
input:

- REBASE's own EMBOSS export `emboss_e.609`, whose `c1`/`c2` columns are the same two
  boundaries this note defines. BsaI reads `GGTCTC 6 2 0 7 11 0 0`; SapI reads
  `GCTCTTC 7 2 0 8 11 0 0`.
- Biopython's `Bio.Restriction` dictionary, compiled from an older REBASE release:
  `BsaI` has `fst5 = 7`, `SapI` has `fst5 = 8`.
- NEBcutter 3's enzyme table: `BsaI {"ct1": 7, "cb1": 11}`, `SapI {"ct1": 8, "cb1": 11}`.

One more cross-check, on **2026-09-18**, is over the supplier values rather than the cut offsets.
Every one already shipped — incubation temperature, heat inactivation and all three methylation
words — was compared against the live page it is cited to, for all 28 records, and agreed. So the
snapshot dates above understate how fresh these values are rather than overstating it, and
nothing in the file needed correcting when `supplied_buffer` was read from the same pages. Worth
saying plainly, because two values look wrong and are not: **SmaI incubates at 37 °C**, not the
25 °C of older literature, and KpnI-HF's `Heat Inactivation: No` is still an answer.

Heat-inactivation temperatures were also checked against NEB's *Heat Inactivation* chart
(archived 2026-03-05), which agrees with every product page read: BsaI-HFv2 80 °C, Esp3I 65 °C,
BspQI 80 °C, BtgZI 80 °C, SapI 65 °C, and `No` for BamHI-HF, KpnI-HF and PstI-HF.

REBASE's methylation-sensitivity pages (`cgi-bin/msget` and `cgi-bin/damlist`) were read as
well. They answer a different question — whether an overlapping Dam, Dcm or CpG site *can*
occur and what was measured when it does — in a vocabulary of their own (`cut`, `blocked`,
`variable`, `some impaired`). The data file carries the supplier's wording instead, because the
supplier's wording is what a bench protocol quotes. REBASE's view is the place to look when a
supplier says nothing.

## 7. Regenerating the file

```sh
pixi run python scripts/build_enzymes.py                    # downloads REBASE withrefm
pixi run python scripts/build_enzymes.py --withrefm FILE    # from a local copy
```

The script fetches REBASE, merges `scripts/enzyme_properties.toml`, and writes
`src/liulab_mbio/data/enzymes.json`. Tests never reach the network: they parse excerpts.

The REBASE release is recorded in the data file as `rebase_version`; this build used **609**.
REBASE updates daily, so a later run may move a cut offset or add an isoschizomer. That is the
point of having the script.

## 8. SrfI and PmeI: two blunt eight-base cutters

Retrieved **2026-09-14**, for issue #60. The iterative library scheme names both by name, and
neither was in the set above.

| Enzyme | REBASE | `top_cut` | `bottom_cut` | `end` | Product page |
| --- | --- | --- | --- | --- | --- |
| SrfI | `GCCC^GGGC` | 4 | 4 | blunt | NEB #R0629, snapshot 2026-03-08 |
| PmeI | `GTTT^AAAC` | 4 | 4 | blunt | NEB #R0560, snapshot 2026-02-01 |

Each cuts the middle of its own eight-base palindrome, so both are Type II by the rule in
section 3 and neither belongs in `GOLDEN_GATE_ENZYMES`. Both pages give 37 °C, heat inactivation
at 65 °C for 20 minutes, and Dam and Dcm insensitivity. CpG is `blocked` for SrfI and
`blocked by some combinations of overlapping` for PmeI. Neither has an unverified field.

**Verdict: the rules in section 2 apply unchanged, and both enzymes ship.** REBASE supplies the
site and both cut offsets under the open-database grant. The two NEB pages are cited and their
handful of facts re-entered in `scripts/enzyme_properties.toml`, never mirrored.

**SmaI is not a substitute for SrfI.** `CCCGGG` lies inside `GCCCGGGC`, so SmaI cuts every SrfI
site and every other `CCCGGG` besides. An enzyme meant to destroy one fragment and spare the
rest has to be the rarer one.

PmeI's isoschizomer is MssI, which Thermo Fisher sells and which reads and cuts identically, so
`get_enzyme` answers to that name too. REBASE lists no other enzyme of SrfI's specificity.

Three readings agree with the offsets above, and none is a build input:

- REBASE's EMBOSS export `emboss_e.609`: `SrfI GCCCGGGC 8 2 1 4 4 0 0`, and `PmeI` the same.
- Biopython's `Bio.Restriction`: `SrfI` and `PmeI` both have `fst5 = 4` and report a blunt cut.
- NEBcutter 3's enzyme table: `SrfI {"ct1": 4, "cb1": 4}`, `PmeI {"ct1": 4, "cb1": 4}`.

## Sources

All read on 2026-09-12, except the two product pages section 8 names and the live re-read of
every product page on 2026-09-18.

- REBASE, *withrefm* — all enzymes with references and isoschizomers, release 609:
  [rebase.neb.com](https://rebase.neb.com/rebase/link_withrefm)
- REBASE, EMBOSS exports `emboss_e.609` and `emboss_r.609`:
  [rebase.neb.com](https://rebase.neb.com/rebase/rebase.f37.html)
- REBASE, *Citing REBASE*: [rebcit.html](https://rebase.neb.com/rebase/rebcit.html)
- REBASE per-enzyme pages, supplier pages (`cgi-bin/ecget`) and methylation sensitivity
  (`cgi-bin/msget`, `cgi-bin/damlist`)
- Roberts, R.J., Vincze, T., Posfai, J., Macelis, D. (2023) REBASE: a database for DNA
  restriction and modification: enzymes, genes and genomes. *Nucleic Acids Res.* 51: D629–D630.
  [doi:10.1093/nar/gkac975](https://doi.org/10.1093/nar/gkac975)
- NEB product pages for R3101, R3156, R3142, R0141, R0180, R3136, R0145, R3138, R3140, R3642,
  R3182, R3104, R0111, R3193, R0146, R3189, R3733, R0739, R0734, R3539, R0569, R0712, R0745,
  R0703, R0629 and R0560 — each via a dated `web.archive.org` snapshot, listed in the data file,
  and each read again live on 2026-09-18 through `r.jina.ai`, the route
  `docs/research/restriction-ligation.md` records for pages `neb.com` 403s
- NEB, *Heat Inactivation* chart — archived snapshot 2026-03-05:
  [web.archive.org](https://web.archive.org/web/20260305053645id_/https://www.neb.com/en-us/tools-and-resources/usage-guidelines/heat-inactivation)
- NEB, *Terms of Use* — archived snapshot 2026-05-23:
  [web.archive.org](https://web.archive.org/web/20260523091525id_/https://www.neb.com/en-us/terms-of-use)
- NEBcutter 3 enzyme table: [nc3.neb.com](https://nc3.neb.com/NEBcutter/data/enzymes.json)
- Thermo Scientific product pages for BpiI (#ER1011) and AarI (#ER1581) — archived snapshots
  2026-08-02, read again live on 2026-09-18
- `docs/research/golden-gate-assembly.md` (issue #4) for NEB's Golden Gate reaction
  temperatures and its per-enzyme table

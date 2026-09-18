---
search:
  exclude: true
---

# Restriction and ligation: the reactions, the buffers, and what may ship

Research note for issue #120, under spec #117. Everything below was retrieved on
**2026-09-18** unless a source carries its own date. It answers the twelve questions #120
asks, says what each source allows, and lists what nobody stated.

## How to read this note

`neb.com` returns HTTP 403 to `curl` and WebFetch for its HTML pages, and its current pages
are a JavaScript shell even when a snapshot of them exists. Five routes work and every fact
below came from one of them.

| Route | Works | Used for |
| --- | --- | --- |
| `neb.com/-/media/nebus/files/manuals/*.pdf` | yes (`curl`) | Monarch kit manuals T1120, T1130, T1020, T1030 |
| `neb.com/en/-/media/catalog/specifications/<letter>/<digit>/<catalog>_v1.pdf` | yes (`curl`) | product specifications (M0202, M0371) |
| `r.jina.ai/<neb url>` | yes (`curl`) | the **live** page, rendered to text: protocols, usage guidelines, product pages, charts |
| Wayback snapshot of an **old-style** `neb.com` URL (no `/en-us/`) | yes | the same pages as they were server-rendered, 2019–2023 |
| `nc3.neb.com/NEBcutter/data/enzymes.json` | yes | cross-check of per-product buffer, activity and temperature |

Three traps, each paid for once:

- **A 403 from `neb.com` can mean "no such file".** The body is the same size either way, so a
  misspelled filename reads as a block. `ps-m0202s-l.pdf` 403s; `m/0/m0202s_l_v1.pdf` under
  `/en/-/media/catalog/specifications/` returns the PDF. Try other spellings before recording
  a document as unreachable.
- **A Wayback snapshot of a current `/en-us/` URL returns HTML with no content.** The body is
  fetched by script. Snapshots of the pre-2024 URLs (`/protocols/0001/01/01/…`,
  `/tools-and-resources/…` without `/en-us/`) are server-rendered and do carry the text.
- **`web.archive.org` was reported offline earlier the same day.** Every number below that a
  plan would print comes from a live page read today; snapshots are used for pages NEB has
  since moved, and each carries its snapshot date.

Values taken from an archived snapshot carry the word **archived** and the date. Anything
still unverified is in [Open gaps](#15-open-gaps) and is never guessed.

## 1. Sources and licences

| Document or dataset | Licence | Can we ship it? |
| --- | --- | --- |
| NEB product pages, protocols, usage guidelines, troubleshooting guides | © NEB, all rights reserved | **No** — cite, and re-enter single facts by hand |
| NEB *Performance Chart for Restriction Enzymes* (the whole activity table) | as above | **No** — read for the shipped enzymes' rows, never republished or made a build input |
| NEB *Restriction Endonuclease Technical Guide* (PDF, v5.0 – 7/17) | as above | **No** — quoted here with citation only |
| NEB Monarch instruction manuals (T1120, T1130) and product specifications (M0202, M0371) | as above | **No** — same rule |
| NEBcutter 3 `enzymes.json` | as above | **No** — read as a cross-check, never a build input |
| REBASE enzyme records | open database, with a citation | **Yes**, as `docs/research/restriction-enzyme-data.md` §2 already records |
| Bauer et al. 2017, *PLoS ONE* 12(12): e0190062 | CC BY 4.0 | **Yes**, with attribution |
| Pheiffer & Zimmerman 1983, *Nucleic Acids Res.* 11, 7853 | © OUP; free to read on PMC | **No** — cite only |
| Addgene protocol pages | site terms: non-commercial use only | **No** — and no number here rests on one |

NEB's terms of use are quoted in full in `docs/research/restriction-enzyme-data.md` §2. The
operative sentence:

> You may not use or copy the Content other than for non-commercial individual reference
> purposes … you may not otherwise copy, display, download, install, use, compile, distribute,
> mirror on another server … or create any derivative work based on such Information.

That verdict — **cite, never mirror** — was set by issue #4 and confirmed by #7, and it holds
here unchanged. What ships is a handful of facts per product, each read from a cited page and
re-entered by hand in `scripts/enzyme_properties.toml`. Section 4 says where that line falls
for a buffer.

Addgene's terms, read 2026-09-18:

> Content may not be reproduced, duplicated, copied, sold, resold, or otherwise exploited for
> any commercial purpose … You are free to use Content for your informational, noncommercial
> purposes, so long as you retain all copyright, Marks … and other proprietary notices.

Non-commercial-only clashes with this repository's MIT licence, the same reason the Golden
Gate note refuses Potapov 2018's CC BY-NC data. Addgene is not used as a source above.

Bauer 2017 carries the PLoS licence line:

> This is an open access article distributed under the terms of the Creative Commons
> Attribution License, which permits unrestricted use, distribution, and reproduction in any
> medium, provided the original author and source are credited.

## 2. The digest reaction

All from NEB, *Optimizing Restriction Endonuclease Reactions* (live, read 2026-09-18), unless
marked.

> By definition, 1 unit of restriction enzyme will completely digest 1 μg of substrate DNA in a
> 50 μl reaction in 60 minutes … most researchers follow the "typical" reaction conditions
> listed, where a 5–10 fold overdigestion is recommended.

A "typical" digestion, as NEB tabulates it:

| Component | Amount |
| --- | --- |
| Restriction enzyme | 10 units, generally 1 µl |
| DNA | 1 µg |
| 10X NEBuffer | 5 µl (1X) |
| Total reaction volume | 50 µl |
| Incubation time | 1 hour |
| Incubation temperature | enzyme dependent |

Smaller volumes, from the same page: 25 µl takes 5 units, 0.5 µg and 2.5 µl of buffer; 10 µl
takes 1 unit, 0.1 µg and 1 µl, and "should not be incubated for longer than 1 hour to avoid
evaporation". The page's rules of thumb: 5–10 units per µg of DNA, 10–20 units per µg for
genomic DNA, and enzyme volume no more than 10% of the reaction.

The *Technical Guide* (v5.0 – 7/17) adds the other volume constraint, in its troubleshooting
table: spin-column purification leaves salt, so "DNA solution should be no more than 25% of
total reaction volume".

**How much DNA one digest takes.** 1 µg in 50 µl is the stated unit of work. For a vector
digest feeding one ligation that is already generous: section 7's ligation takes 50 ng of a
4 kb vector.

**Star activity.** The conditions and the countermeasures are NEB's own pairing, from the
*Technical Guide* (v5.0 – 7/17) and the *Optimizing* page:

| Condition that contributes | What NEB says to do |
| --- | --- |
| High glycerol, > 5% v/v | Enzyme added is at most 10% of the reaction; use the standard 50 µl volume |
| High enzyme per µg of DNA (varies by enzyme, usually 100 units/µg) | Use the fewest units that digest |
| Non-optimal buffer | Use the supplied buffer wherever possible |
| Prolonged reaction time | Use the minimum time for complete digestion |
| Organic solvents (DMSO, ethanol, ethylene glycol, dimethylacetamide, dimethylformamide, sulphalane) | Keep the DNA prep free of them |
| Mg²⁺ substituted by Mn²⁺, Cu²⁺, Co²⁺ or Zn²⁺ | Use Mg²⁺ |

NEB's own caveat, quoted: "The relative significance of each of these altered conditions will
vary from enzyme to enzyme," and the HF enzymes "have been engineered for reduced star
activity". Eleven of the sixteen multiple-cloning-site enzymes this package ships are HF
products already; NdeI, SmaI, XbaI, XhoI and XmaI are not.

**Cutting a PCR product.** NEB, *Cleavage Close to the End of DNA Fragments* (archived
2021-04-19): "As a general rule and for enzymes not listed below, 6 base pairs should be added
on either side of the recognition site to cleave efficiently." The per-enzyme table did render
in that snapshot — a gap the Golden Gate note left open — and it is measured as percent
cleavage of an annealed oligonucleotide (– is 0%, + is 0–20%, ++ is 20–50%, +++ is 50–100%).
For the products this package names, at 1 / 2 / 3 / 4 / 5 spacer bases:

| Product | 1 | 2 | 3 | 4 | 5 |
| --- | --- | --- | --- | --- | --- |
| BamHI-HF | + | + | +++ | +++ | +++ |
| EcoRI-HF | + | + | ++ | +++ | +++ |
| HindIII-HF | – | + | +++ | +++ | +++ |
| KpnI-HF | + | +++ | +++ | +++ | +++ |
| NcoI-HF | + | ++ | +++ | +++ | +++ |
| NdeI | + | + | +++ | +++ | +++ |
| NotI-HF | ++ | ++ | ++ | ++ | ++ |
| PstI-HF | ++ | +++ | +++ | +++ | +++ |
| SacI-HF | – | + | +++ | +++ | +++ |
| SalI-HF | – | ++ | +++ | +++ | +++ |
| SbfI-HF | ++ | +++ | +++ | +++ | +++ |
| SmaI | +++ | +++ | +++ | +++ | +++ |
| SphI-HF | ++ | ++ | +++ | +++ | +++ |
| XbaI | ++ | ++ | ++ | ++ | ++ |
| XhoI | ++ | ++ | ++ | +++ | +++ |
| XmaI | +++ | +++ | +++ | +++ | +++ |

Two readings matter for a tailed primer: only SmaI and XmaI cut well at a single base, and NotI
and XbaI never reach +++ even at five, so the six-base default is the right one and is what
`sites.primer_tail` already builds.

## 3. Double digest and buffer compatibility

NEB, *Double Digests* (live, read 2026-09-18), quoted:

> Double digests with NEB's restriction enzymes can be set up in rCutSmart Buffer. Otherwise,
> choose an NEBuffer that results in the most activity for both enzymes.
>
> The final concentration of glycerol in any reaction should be less than 5% … For example, in
> a 50 µl reaction, the total amount of enzyme added should not exceed 5 µl.
>
> If two different incubation temperatures are necessary, choose the optimal reaction buffer
> and set up reaction accordingly. Add the first enzyme and incubate at the desired
> temperature. Then, heat inactivate the first enzyme, add the second enzyme and incubate at
> the recommended temperature.
>
> NEB currently supplies two enzymes with unique buffers: EcoRI and DpnII … Note that EcoRI has
> an HF version which is supplied with rCutSmart Buffer.

**When a sequential digest is needed**, in NEB's own words on the same page:

> If there is no buffer in which the two enzymes exhibit > 50% activity, a sequential digest
> can be performed.

The sequential recipe: digest first with the enzyme whose recommended buffer has the lowest
salt, adjust the salt to suit the second, add the second enzyme; or put the DNA through a spin
column between the two.

The 2021 snapshot of this page named three unique-buffer enzymes (EcoRI, SspI, DpnII); the live
page names two. Nothing else in the text has changed.

**Where the compatibility numbers live, and what may be done with them.** NEB publishes them
twice: once as the *Performance Chart for Restriction Enzymes*, a single table covering more
than 200 enzymes with a "Supplied NEBuffer" column and four "% Activity in NEBuffer" columns;
and once per enzyme, on that enzyme's own product page, under **Reaction Conditions** and
**Activity in NEBuffers**. Read live on 2026-09-18, BsmBI-v2's page (NEB #R0739) gives
"1X NEBuffer r3.1", "Incubate at 55°C" and "NEBuffer r1.1: <10%, NEBuffer r2.1: 50%, NEBuffer
r3.1: 100%, rCutSmart Buffer: 25%"; SmaI's page (NEB #R0141) gives "1X rCutSmart Buffer",
"Incubate at 37°C" and "rCutSmart Buffer: 100%". The chart's own row for each agrees, read the
same day.

The chart is a compiled table and is covered by the terms in section 1: **it is read, never
republished, and never a build input**, here or in `scripts/build_enzymes.py`. The same holds for
NEBcutter 3's `enzymes.json`, which carries the same per-product fields (`buff_name`,
`buff_act`) and which #7 already ruled out as a build input. The per-enzyme values on a product
page are the same kind of fact as the incubation temperature and the methylation words this
package already ships from those pages. Section 4 turns that into a verdict.

## 4. A buffer field on the shipped records

**Verdict: `Enzyme` may gain a reaction buffer and its activity values, by the rule already in
force.** Each value is read from that enzyme's own product page — the page
`scripts/enzyme_properties.toml` already cites for `incubation_celsius`,
`heat_inactivation_celsius` and `methylation` — and re-entered by hand beside them. No NEB file
is fetched, no chart is copied, and `build_enzymes.py` goes on downloading REBASE and nothing
else.

Two fields, and what each buys:

| Field | Source per enzyme | What a plan can say with it |
| --- | --- | --- |
| `supplied_buffer` (a string) | product page, "Reaction Conditions" | Two enzymes supplied in the same buffer digest together in it |
| `buffer_activity` (buffer → percent) | product page, "Activity in NEBuffers" | NEB's own rule: pick the buffer where both are most active; sequential when no buffer gives both > 50% |

The narrower field answers the common case on its own. Of the twenty-six NEB products this
package ships, twenty-four are supplied in rCutSmart Buffer — including **all sixteen
multiple-cloning-site enzymes** — and only BsmBI-v2 and BspQI are supplied in NEBuffer r3.1.
AarI and BpiI are Thermo Fisher products and sit outside NEB's buffer system entirely. So every
pair drawn from the sixteen shares a buffer, and a plan that ships `supplied_buffer` alone can
already give a verdict for every multiple-cloning-site pair it will be asked about. (Surveyed
2026-09-18 across the shipped set, on the live Performance Chart and against NEBcutter 3's
table; confirmed on the R0739 and R0141 product pages. PaqCI's entry reads "rCutSmart Buffer +
Activator", which is the activator the Golden Gate note already records. No value from either
source is a build input.)

The wider field is what lets a plan judge a mixed pair — a multiple-cloning-site enzyme with
BsmBI-v2, say, where rCutSmart leaves BsmBI-v2 at 25% and NEBuffer r3.1 leaves most of the
others at 10%, so NEB's own > 50% rule calls for a sequential digest. That is a real answer the
narrower field cannot give.

**One thing neither field carries.** The Performance Chart flags, per cell, whether a
non-optimal buffer risks star activity; the product page's activity block does not. A plan that
proposes a buffer other than the supplied one cannot check that flag from a shippable source,
so it should prefer the supplied buffer, and say plainly when it is not using one.

**The conservative alternative**, if a maintainer reads NEB's terms more strictly than #7 did:
ship neither field, carry no verdict, and point the user at NEBcloner or the Double Digest
Finder. The spec already describes that fallback. This note recommends the first reading,
because it is the same reading that already put a temperature in the data file.

## 5. Which enzymes are missing

The package ships 28 enzymes, of which 16 are the multiple-cloning-site enzymes this method
uses. Measured against the repository's own vector fixture, `tests/data/pUC19.dna`:

- Twelve shipped enzymes cut once inside pUC19's multiple cloning site — EcoRI, SacI, KpnI,
  SmaI, XmaI, BamHI, XbaI, SalI, SbfI, PstI, SphI, HindIII — and each cuts nowhere else in the
  plasmid.
- Sixteen shipped enzymes cut pUC19 exactly once in total; BsmBI and Esp3I cut twice; ten do
  not cut it at all.
- Against the repository's own pinned list of commercial enzymes with a six-base-or-longer site
  that cut pUC19 once (`tests/data/pUC19-unique-6-cutters.tsv`), every unshipped enzyme in that
  region but three either reads the same site as one the package ships or reads a site
  containing it — Acc65I for KpnI, Cfr9I and TspMI for XmaI, Psp124BI for SacI, SdaI and
  Sse8387I for SbfI, PaeI for SphI, ApoI, AcsI and XapI for EcoRI. The three genuinely absent
  are **AccI** and **HincII**, which share pUC19's SalI position, and **BspMI**, a Type IIS
  enzyme whose site sits between PstI and SphI.

So for a pUC19-class vector the shipped set is already sufficient, and the enzyme set does not
have to widen for this method to land.

**What adding one costs.** `scripts/build_enzymes.py` says "the enzymes in the data file are
the ones `enzyme_properties.toml` names", and REBASE supplies the site, both cut offsets and
the isoschizomers. So one enzyme is: one `[sources.*]` block (title, URL, retrieval or snapshot
date) and one `[enzymes.*]` block of eight values — `commercial_name`, `catalog_number`,
`supplier`, `incubation_celsius`, `heat_inactivation_celsius`, `heat_inactivation_minutes`,
`methylation` (three words), `source` — each read from one cited product page, then
`pixi run python scripts/build_enzymes.py`. No code changes, and the data file is never
hand-edited.

**What this cannot answer.** The lab's own vector collection is not in this repository. pUC19
and a GFP coding sequence are the only records under `tests/data/`, and `docs/examples/` builds
on the same pair. The question "which enzymes do our vectors carry that we do not ship" needs
those vector files; it is in [Open gaps](#15-open-gaps) with what to ask for.

## 6. Dephosphorylation

**When.** When both vector ends can anneal to each other — one enzyme, or two enzymes leaving
compatible ends, or two blunt ends — the vector closes on itself and the plate fills with empty
vector. NEB's rSAP product page (archived 2021-06-24): "In cloning, dephosphorylation prevents
religation of linearized plasmid DNA. The enzyme acts on 5´ protruding, 5´ recessed and blunt
ends." A digest with two enzymes leaving incompatible ends needs no phosphatase, and NEB's
cloning controls treat those two cases as the same control (section 12).

**The reaction**, from NEB's live protocol *Dephosphorylation of 5´-ends of DNA using rSAP*
(NEB #M0371), read 2026-09-18:

| Component | 20 µl reaction |
| --- | --- |
| DNA | 1 pmol of DNA ends |
| rCutSmart Buffer (10X) | 2 µl |
| rSAP | 1 unit |
| Nuclease-free water | to 20 µl |

Incubate at 37°C for 30 minutes; stop by heat inactivation at **65°C for 5 minutes**. NEB's own
note: "1 pmol of DNA ends is about 1 µg of a 3 kb plasmid," and larger volumes scale
proportionally.

**In the digest, which is what a plan should do.** The same page:

- "The phosphatase can be added directly into the digestion reaction during or after DNA
  digestion."
- "rSAP is active in all NEB restriction enzyme buffers."
- "The restriction enzyme should be heat inactivated at the same time as the phosphatase after
  digest and dephosphorylation."
- "If restriction enzyme(s) cannot be heat inactivated, DNA purification is required before
  ligation."

The archived in-digest protocol (2020-09-29) gives the dose the same way — "Add 1 unit of rSAP
for every 1 pmol of DNA ends (about 1 μg of a 3 kb plasmid) and incubate at 37°C for 30–60
minutes" — after a 1–5 µg digest in 20 µl.

**How it is stopped, and why it matters.** rSAP is "completely and irreversibly inactivated by
heating at 65°C for 5 minutes, thereby making removal of rSAP prior to ligation or end-labeling
unnecessary" (rSAP product page, archived 2021-06-24). NEB's cloning troubleshooting guide
lists "inefficient dephosphorylation" and a live phosphatase as separate causes of background
and of failed ligation, and says to heat inactivate or remove the phosphatase before the
ligation.

**The alternative.** Antarctic Phosphatase (NEB #M0289), archived protocol 2020-08-05: 1 pmol of
DNA ends, 2 µl of its own 10X buffer, **5 units**, to 20 µl, 37°C for 30 minutes, stopped at
**80°C for 2 minutes**. It is "active in all NEB restriction enzyme buffers only when
supplemented with Antarctic Phosphatase Reaction Buffer, which provides Zn²⁺". One unit of rSAP
against five units of Antarctic Phosphatase for the same pmol of ends, and a lower stop
temperature, is why rSAP is the default here.

## 7. The ligation reaction

From NEB, *Ligation Protocol with T4 DNA Ligase (NEB #M0202)* (live, read 2026-09-18). NEB's
precaution comes first: "Heat-inactivate restriction enzymes or purify DNA using a spin column
(NEB #T1130) before ligation."

| Component | 20 µl reaction | Final amount |
| --- | --- | --- |
| T4 DNA Ligase Buffer (10X) | 2 µl | 1X |
| Vector DNA (4 kb) | X µl | 50 ng (0.020 pmol) |
| Insert DNA (1 kb) | X µl | 37.5 ng (0.060 pmol) |
| Nuclease-free water | to 20 µl | |
| T4 DNA Ligase (400,000 u/ml) | 1 µl | 400 units, or 2,000 units at high concentration |

The table is a 1:3 vector-to-insert molar ratio for those two sizes. The page's own guidance on
amounts:

> We recommend a vector amount of 20-30 fmol with an overall concentration of vector + insert
> between 1-10 ng/μl for efficient ligation. Concentration lower than 1 ng/μL may result in
> intramolecular ligation/circularization of the fragments. Vector: Insert molar ratios between
> 1:1 and 1:10 are optimal for single insertions (up to 1:20 for short adaptors).

Temperature and time, quoted:

| Ends | Incubation |
| --- | --- |
| Cohesive (sticky) | 16°C overnight, **or** room temperature for 10 minutes |
| Blunt, or a single-base overhang | 16°C overnight, **or** room temperature (25°C) for 2 hours, **or** high-concentration (2,000,000 u/ml) T4 DNA Ligase for 10 minutes at 25°C |

Then: heat-inactivate at 65°C for 10 minutes, chill on ice, and transform 1–5 µl of the
reaction into 50 µl of competent cells.

**The ligase and its unit.** NEB product specification PS-M0202S/L v1.0 (effective 23 Jun 2016,
fetched 2026-09-18): T4 DNA Ligase, 400,000 units/ml, "One unit is defined as the amount of
enzyme required to give 50% ligation of 6 µg of Lambda-HindIII DNA in 30 minutes at 16°C in a
total reaction volume of 20 µl." The product page (archived 2021-06-10) states the same unit
against a 5´ DNA termini concentration of 0.12 µM and both concentrations sold, 400,000 and
2,000,000 cohesive end units/ml. ATP is required: the live protocol notes ligation can run in
any of the four restriction NEBuffers or T4 PNK buffer "if they are supplemented with 1 mM of
ribo ATP … Deoxyribo ATP will not work."

**The five-minute alternative.** Quick Ligation Kit (NEB #M2200), archived protocol 2023-03-24:
10 µl of 2X Quick Ligase Reaction Buffer, the same 50 ng / 37.5 ng of vector and insert, 1 µl
Quick Ligase, to 20 µl; 5 minutes at 25°C; transform 1–5 µl. Two rules from that page are its
own: "Do not heat inactivate – heat inactivation dramatically reduces transformation
efficiency," and the ratio guidance, "Insert:vector ratios between 2 and 6 are optimal for
single insertions. Ratios below 2:1 result in lower ligation efficiency. Ratios above 6:1
promote multiple inserts." Its buffer contains PEG, which matters twice: it is why the reaction
is fast (section 8) and why the product must be cleaned up before electroporation.

## 8. Blunt against sticky

The only cost a supplier states is time, and it is stated in the same table for the same enzyme
and the same units: **10 minutes at room temperature for cohesive ends, 2 hours for blunt ends
or a single-base overhang** — a twelvefold longer incubation — or the same 10 minutes with
five times the ligase (high-concentration T4 at 2,000,000 u/ml). NEB states no ratio of
colonies, no yield and no fold difference anywhere read.

What is *measured*, and it cuts against the folklore: Bauer et al. 2017 (*PLoS ONE*, CC BY 4.0)
profiled T4 DNA ligase on defined blunt and cohesive oligonucleotide junctions by capillary
electrophoresis, 1 µM ligase, 100 nM substrate, 20 minutes at 25°C. Their finding, quoted:
"T4 DNA ligase, the most common enzyme utilized for in vitro ligation, had its greatest
activity on blunt- and 2-base overhangs, and poorest on 5′-single base overhangs," with "T4 DNA
ligase sealing over 60% of the A/T Blunt substrate in the presence or absence of PEG, and 50%
of the G/C blunt substrate in the presence of PEG". So the sealing chemistry is not what makes
blunt cloning hard; the hard part is holding two ends together long enough, which base-paired
overhangs do and blunt ends do not.

That is also what the compensations are for, and each is sourced:

- **More ligase**, or longer: the two rows of NEB's own table above.
- **A crowding agent.** Pheiffer & Zimmerman 1983 (*Nucleic Acids Res.* 11, 7853), abstract:
  "The rates of blunt-end and cohesive-end ligation of DNA by T4 DNA ligase are increased by
  orders of magnitude in the presence of high concentrations of a variety of nonspecific
  polymers such as polyethylene glycol, Ficoll, bovine plasma albumin, or glycogen." NEB's
  Quick Ligation Kit buffer is one: its own page tells you to reduce the PEG before
  electroporation, and the kit ligates "both blunt and cohesive (sticky) ends in 5 minutes at
  room temperature" (M0202 page, live). For what such a buffer holds, Bauer 2017's methods give
  the NEBNext Quick Ligation reaction buffer as 66 mM Tris pH 7.6, 10 mM MgCl₂, 1 mM DTT, 1 mM
  ATP and 6% PEG 6000.
- **Dephosphorylate the vector.** Two blunt ends are always compatible with each other, so
  section 6 applies to every blunt ligation, not only to the single-enzyme case.

One more measured caveat worth carrying: Bauer 2017 found end-joining "greatly reduced under
physiologically relevant ionic strengths" — 150 mM NaCl eliminated nearly all activity except
on the four-base overhang. A ligation set up in leftover digest buffer with high salt is a
ligation that will not work.

**Not sourced anywhere read:** a colonies-per-reaction ratio between a blunt and a cohesive
cloning of the same insert. See [Open gaps](#15-open-gaps).

## 9. Gel purification

**When it is needed rather than optional.** Three cases, each with a supplier line behind it:

- **A digest that leaves more than the backbone.** Cutting a vector with two enzymes whose
  sites are far apart releases a stuffer that religates into the backbone as readily as the
  insert does. NEB's heat-inactivation page offers three ways out of a finished digest — a
  column, a gel, or phenol/chloroform — and only the gel separates two fragments from each
  other, so where a stuffer comes out the choice is made for you.
- **Background from uncut or religated vector.** NEB FAQ, *How can I reduce the number of
  vector-only background colonies?* (archived 2019-11-23): "the vector should be a PCR product
  rather than a restriction fragment. If background continues to be a problem, the
  PCR-amplified vector can be treated with DpnI to remove the template carry-over … extracted
  from an agarose gel following electrophoresis."
- **The wrong PCR product.** NEB's cloning troubleshooting guide lists "Incorrect PCR amplicon
  was used during cloning" under colonies with the wrong construct, and answers "Gel purify the
  correct PCR fragment."

**What it costs.** NEB, *Monarch Spin DNA Gel Extraction Kit* instruction manual (NEB #T1120,
version 2.0 10.25, fetched 2026-09-18): typical recovery **70–90% for fragments below 10 kb and
50–70% at 10 kb or above**, binding capacity 5 µg, elution from 6 µl, purity A260/280 > 1.8.
The older Monarch DNA Gel Extraction Kit manual (NEB #T1020, version 2.1_4/21) states the same
band, 70–90% for 50 bp–10 kb and 50–70% for 11–23 kb, and warns that a high agarose percentage
"will decrease the recovery yield due to incomplete extraction".

So the working figure is: **plan on losing a fifth of the DNA, and up to half of a fragment
above 10 kb.** The digest amounts in section 2 are large enough to absorb that.

## 10. Between the digest and the ligation

**Heat inactivation**, from NEB's *Heat Inactivation* usage guideline (archived 2021-04-20):

> Incubation at 65°C for 20 minutes inactivates the majority of restriction endonucleases that
> have an optimal incubation temperature of 37°C. Enzymes that cannot be inactivated at 65°C
> can often be inactivated by incubation at 80°C for 20 minutes.
>
> For enzymes that cannot be heat-inactivated at 65°C or 80°C, we recommend using a column for
> cleanup … or running the reaction on an agarose gel and then extracting the DNA … or
> performing a phenol/chloroform extraction.

The per-enzyme answer already ships: `Enzyme.heat_inactivation_celsius` and
`heat_inactivation_minutes` come from each product page, and for **KpnI-HF, BamHI-HF and
PstI-HF the supplier's answer is "No"** — heat does not stop them, so a plan using any of the
three has to clean up on a column or a gel before ligating. Of the sixteen multiple-cloning-site
enzymes, those three are the only ones with no heat answer; the rest are 65°C for 20 minutes,
except HindIII, NcoI and SbfI at 80°C for 20 minutes.

**What each choice costs.**

| Route | Yield cost | What it leaves behind |
| --- | --- | --- |
| Heat, 65 or 80°C for 20 min | none | Salt, glycerol and the denatured enzyme stay in the tube, so the digest must stay a minority of the ligation volume |
| Spin column (NEB #T1130, manual v1.0 06.24, fetched 2026-09-18) | typical recovery **70–90%** | Nothing of the digest; but column eluates carry salt, so the *Technical Guide*'s 25%-of-volume rule applies to the next reaction |
| Gel extraction (section 9) | 70–90%, or 50–70% above 10 kb | Nothing, and it is the only route that also removes a stuffer |

NEB's own precaution on the ligation protocol reduces to one line: heat-inactivate, or spin
column, before the ligase. A phosphatase in the same tube is stopped by the same heat step
(section 6).

## 11. Dam and Dcm methylation

**What methylates what.** NEB *Technical Guide* (v5.0 – 7/17): Dam methyltransferase methylates
"the N6 position of the adenine in the sequence GATC"; Dcm methylates "the C5 position of
cytosine in the sequences" CCWGG. "Most laboratory strains of E. coli contain three site-specific
DNA methyltransferases," and "Mammalian and plant DNA that has been cloned into a methylating
E. coli strain will be Dam/Dcm methylated. Most commonly used laboratory E. coli strains
methylate DNA." Genomic DNA taken straight from a mammalian source is neither Dam nor Dcm
methylated, but is CpG methylated.

**Which shipped enzymes care.** Read from the shipped data file, whose values come from the
product pages `scripts/enzyme_properties.toml` cites. Among the sixteen multiple-cloning-site
enzymes:

- **Dam**: only **XbaI**, which its supplier calls "blocked by overlapping". Every other one is
  "not sensitive".
- **Dcm**: none.
- **CpG**, which matters only for DNA from a eukaryotic source: NotI, SalI and SmaI are
  "blocked"; XhoI and XmaI are "impaired"; EcoRI and SacI are "blocked by some combinations of
  overlapping".

**"Blocked by overlapping" is a property of the flanks, not of the site**, and it is
computable rather than assumed. A Dam site is GATC; XbaI reads TCTAGA. A GATC overlaps that
site only when the site is followed by TC (TCTAGATC) or preceded by GA (GATCTAGA). Nothing in
the enzyme record knows which; the record says the enzyme is sensitive and the sequence decides.
The package already holds both halves — `sites.DCM_SITE` is `CCWGG` and
`cloning.goldengate.assembly.DAM_SITE` is `GATC` — and `sites.primer_tail` already refuses a
tail that puts a Dcm site across the enzyme's site. The same check over a vector's own flanks is
what this method needs.

**Which strain.** *Technical Guide*: "If the enzyme is inhibited by Dam or Dcm methylation, grow
the plasmid in a dam⁻/dcm⁻ strain (NEB #C2925)." But that strain's own product page (archived
2021-05-07) adds the warning that decides *when*:

> Note that dam⁻ strains are not recommended as a host for primary cloning/ligation. The dam
> mutation can result in an increased mutation rate in the cell and a reduction in the
> transformation efficiency. DNA should be maintained in a dam⁺ strain unless there is a
> specific need for DNA free of Dam or Dcm methylation.

The same page gives its transformation efficiency as 1–3 × 10⁶ cfu/µg of pUC19, which is what
"a reduction in the transformation efficiency" means in practice. So the order is: transform the
ligation into the ordinary host,
then re-transform the finished plasmid into C2925 only to prepare DNA for a digest a
methylation-blocked enzyme has to make.

## 12. What the plate should look like

**No supplier states an absolute colony count for this method.** What NEB does state is a set of
four controls and what each should give relative to the others, from the *Troubleshooting Guide
for Cloning* (live, read 2026-09-18; identical in the archived 2021-05-16 copy):

1. Transform 100 pg–1 ng of **uncut vector** — cell viability, transformation efficiency, and
   that the antibiotic is right.
2. Transform the **cut vector** — background from undigested plasmid. "The number of colonies
   in this control should be **< 1%** of the number of colonies in the uncut plasmid control."
3. Transform a **vector-only ligation**. Its ends should not be able to religate, "because
   either they are incompatible … or the 5´ phosphate group has been removed". "This control
   transformation should yield the same number of colonies as control #2."
4. Digest the vector with a **single** enzyme, religate and transform. Its ends are compatible,
   so it "should result in approximately the same number of colonies as control #1."

Two absolute numbers do appear, and they are about the cells rather than the ligation: a
transformation efficiency below 10⁴ means the competent cells should be remade or replaced
(same guide), and NEB's dam⁻/dcm⁻ cells are 1–3 × 10⁶ cfu/µg of pUC19 (section 11).

So what a plan can promise is a *ratio*, not a count: the empty-vector background is measured on
its own plate, and controls 2 and 3 are what "an empty-vector background beside it" actually
means. The fraction of colonies carrying the insert is not stated by any source read; see
[Open gaps](#15-open-gaps).

## 13. Troubleshooting

From NEB's *Troubleshooting Guide for Cloning* (live, read 2026-09-18) and the *Technical
Guide*'s own table (v5.0 – 7/17). Only the rows this method can hit are kept.

| Problem | Cause NEB names | The fix NEB names |
| --- | --- | --- |
| Few or no transformants | Cells not viable | Transform uncut pUC19; below 10⁴ cfu/µg, replace the cells |
| | Inefficient ligation | At least one fragment must carry a 5´ phosphate; vary vector:insert from 1:1 to 1:10 (1:20 for short adaptors); purify away salt and EDTA; fresh buffer, because ATP degrades over freeze–thaws; heat-inactivate or remove the phosphatase |
| | Too much ligation mix transformed | Use < 5 µl of the ligation |
| | Enzyme did not cleave completely | Check methylation sensitivity; use the supplied buffer; clean up the DNA; leave at least 6 bases between the site and the end of a PCR fragment |
| | PEG in the mix, with electrocompetent cells | Drop dialysis, or a spin column, before electroporation |
| Too much background | Inefficient dephosphorylation | Heat-inactivate or remove the restriction enzymes before the phosphatase |
| | Kinase still active | Heat-inactivate it; "active kinase will re-phosphorylate the dephosphorylated vector" |
| | Enzyme did not cleave completely | Methylation, wrong buffer, contaminated DNA |
| No ligated product on a gel | Inefficient ligation | As above; single-base overhangs are "most difficult" and want Quick Ligase or concentrated T4; test the ligase on Lambda-HindIII DNA |
| Ligation ran as a smear | Ligase bound to the DNA | Proteinase K before loading |
| Digest ran as a smear | Enzyme bound to the DNA | Fewer units; SDS (0.1–0.5%) in the loading buffer, or Gel Loading Dye Purple; fresh gel and running buffer |
| | Nuclease contamination | Clean up the DNA |
| Incomplete digest | Salt inhibition | Clean up the DNA; keep it under 25% of the reaction volume |
| | Too few units | At least 3–5 units per µg, 1–2 hours |
| | Supercoiled substrate, or a slow site | More units, or a longer incubation |
| | Methylation | Check the enzyme's sensitivity; grow the plasmid in NEB #C2925 |
| Extra bands | Star activity | An HF enzyme; fewer units; enzyme under 10% of the volume; shorter incubation; the supplied buffer |
| Colonies without the insert | Recombination, toxicity, internal site | recA⁻ strain; plates at 25–30°C; check the insert for the enzyme's site |

## 14. What this settles for the spec

The three conditional decisions #117 left to this note:

1. **May a buffer or activity field ship?** Yes — `supplied_buffer`, and optionally the four
   activity percentages, re-entered by hand from each enzyme's own product page into
   `scripts/enzyme_properties.toml`, exactly as its temperature already is. The Performance
   Chart and NEBcutter's table stay out of the build (section 4).
2. **Does the enzyme set widen?** Not for this method to land. All sixteen multiple-cloning-site
   enzymes are supplied in one buffer and cover pUC19's whole site; the only absentees there are
   AccI, HincII and BspMI. Widening costs one `enzyme_properties.toml` entry per enzyme and a
   rebuild (section 5). The lab's own vectors may change this answer, and they are not in the
   repository.
3. **Can a compatibility verdict be given at all?** Yes, and by NEB's own rule rather than a
   guess: both enzymes supplied in the same buffer, digest together in it; otherwise the buffer
   where both are most active; sequential when no buffer leaves both above 50%; and always
   glycerol under 5%, which is enzyme volume under 10% of the reaction. Where the two enzymes
   want different temperatures, NEB's answer is sequential with a heat inactivation between,
   which for the sixteen is never needed — all sixteen incubate at 37°C.

What the method's own `bench` module can compute, each with a line above it: digest units and
volume (§2), star-activity conditions to print rather than predict (§2), spacer bases for a
tailed primer (§2), the buffer verdict (§3, §4), the phosphatase dose and its stop (§6), the
ligation table, ratio and both incubations (§7), the clean-up route and its yield (§9, §10), the
methylation check and the strain (§11), and the four transformation controls (§12).

## 15. Open gaps

Still unverified. None of these is guessed above.

| Item | Why it is missing | Workaround in use |
| --- | --- | --- |
| Colonies per ligation, in absolute numbers | No supplier or paper read states one for a generic restriction–ligation cloning | The four relative controls in §12 |
| Fraction of colonies carrying the insert | Same | Colony PCR distinguishes them; `bench.validation` already designs it |
| A measured blunt-against-cohesive ratio, in colonies | Searched NEB, PLoS, NAR and PMC; the measurements found are sealing yields on oligonucleotides (Bauer 2017), not transformants | §8: NEB's twelvefold incubation difference, and its own compensations |
| Lund, Duch & Pedersen 1996, *NAR* 24, 800, the transformant counts behind temperature-cycle ligation | The PMC copy is a scan; the PDF route returns HTML and PubMed holds no abstract | Not cited; nothing above rests on it |
| Whether NEB's protocols.io postings are CC BY | The rendered page shows no licence; the API needs a token | Every protocol number above is cited from `neb.com` under the cite-never-mirror rule |
| Star-activity flags per buffer cell | They are on the Performance Chart, which does not ship; the product page's activity block omits them | §4: prefer the supplied buffer and say when the plan does not |
| Whether DNA damage from a UV gel box costs ligation efficiency | Not stated in any manual read | §9 records recovery only |
| Which enzymes the lab's vectors carry | No lab vector file is in the repository | Ask for the vector records (`.dna` or GenBank); §5 answers it for pUC19 |

## Sources

Read 2026-09-18 unless a snapshot date is given. Nothing below is redistributed; the downloads
sit in `reference_docs/restriction-ligation/`, which is git-ignored.

- NEB, *Optimizing Restriction Endonuclease Reactions*:
  [usage guideline](https://www.neb.com/en-us/tools-and-resources/usage-guidelines/optimizing-restriction-endonuclease-reactions)
  (live; also archived 2021-04-22)
- NEB, *Double Digests*:
  [usage guideline](https://www.neb.com/en-us/tools-and-resources/usage-guidelines/double-digests)
  (live; also archived 2021-04-22)
- NEB, *NEBuffer Activity/Performance Chart with Restriction Enzymes* — read for the shipped
  enzymes' rows only, live and archived 2021-05-14:
  [usage guideline](https://www.neb.com/en-us/tools-and-resources/usage-guidelines/nebuffer-performance-chart-with-restriction-enzymes)
- NEB, *Heat Inactivation* — archived 2021-04-20:
  [usage guideline](https://www.neb.com/en-us/tools-and-resources/usage-guidelines/heat-inactivation)
- NEB, *Cleavage Close to the End of DNA Fragments* — archived 2021-04-19:
  [usage guideline](https://www.neb.com/en-us/tools-and-resources/usage-guidelines/cleavage-close-to-the-end-of-dna-fragments)
- NEB, *Ligation Protocol with T4 DNA Ligase (NEB #M0202)*:
  [protocol](https://www.neb.com/en-us/protocols/dna-ligation-with-t4-dna-ligase-m0202)
  (live; also archived 2023-01-26)
- NEB, *Quick Ligation Protocol (M2200)* — archived 2023-03-24:
  [protocol](https://www.neb.com/protocols/0001/01/01/quick-ligation-protocol)
- NEB, *Protocol for Dephosphorylation of 5´-ends of DNA using rSAP (NEB #M0371)*:
  [protocol](https://www.neb.com/en-us/protocols/protocol-for-dephosphorylation-of-5-ends-of-dna-neb-m0371)
  (live; the in-digest variant archived 2020-09-29)
- NEB, *Protocol for Dephosphorylation of 5´-ends of DNA using Antarctic Phosphatase (NEB
  #M0289)* — archived 2020-08-05:
  [protocol](https://www.neb.com/protocols/0001/01/01/vector-dephosphorylation-protocol)
- NEB, *Troubleshooting Guide for Cloning*:
  [guide](https://www.neb.com/en-us/tools-and-resources/troubleshooting-guides/troubleshooting-guide-for-cloning)
  (live; also archived 2021-05-16)
- NEB FAQ, *How can I reduce the number of vector-only background colonies?* — archived
  2019-11-23
- NEB product pages, read live: SmaI (#R0141), BsmBI-v2 (#R0739), EcoRI-HF (#R3101); archived:
  T4 DNA Ligase (#M0202, 2021-06-10), Shrimp Alkaline Phosphatase (#M0371, 2021-06-24), Quick
  Ligation Kit (#M2200, 2021-04-19), dam⁻/dcm⁻ Competent E. coli (#C2925, 2021-05-07)
- NEB product specifications PS-M0202S/L v1.0 and PS-M0371S/L v1.0, fetched under
  `/en/-/media/catalog/specifications/`
- NEB, *Restriction Endonuclease Technical Guide*, version 5.0 – 7/17, via the `neb-online.de`
  mirror:
  [PDF](https://www.neb-online.de/literatur/pdf/Restriction_Endonuclease_Technical_Guide.pdf)
- NEB instruction manuals, fetched under `/-/media/nebus/files/manuals/`: Monarch Spin DNA Gel
  Extraction Kit (#T1120, v2.0 10.25), Monarch Spin PCR & DNA Cleanup Kit (#T1130, v1.0 06.24),
  Monarch DNA Gel Extraction Kit (#T1020, v2.1_4/21), Monarch PCR & DNA Cleanup Kit (#T1030,
  v3.0_4/21)
- NEBcutter 3 enzyme table, cross-check only:
  [nc3.neb.com](https://nc3.neb.com/NEBcutter/data/enzymes.json)
- Bauer, R.J., Zhelkovsky, A., Bilotti, K., Crowell, L.E., Evans, T.C., McReynolds, L.A. and
  Lohman, G.J.S. (2017) Comparative analysis of the end-joining activity of several DNA ligases.
  *PLoS ONE* 12(12): e0190062.
  [doi:10.1371/journal.pone.0190062](https://doi.org/10.1371/journal.pone.0190062) (CC BY 4.0)
- Pheiffer, B.H. and Zimmerman, S.B. (1983) Polymer-stimulated ligation: enhanced blunt- or
  cohesive-end ligation of DNA or deoxyribooligonucleotides by T4 DNA ligase in polymer
  solutions. *Nucleic Acids Res.* 11, 7853–7871.
  [doi:10.1093/nar/11.22.7853](https://doi.org/10.1093/nar/11.22.7853)
- Addgene, *Terms of Use* — for the licence verdict only:
  [addgene.org](https://www.addgene.org/terms-of-use/)
- `docs/research/restriction-enzyme-data.md` (issue #7) for the REBASE and NEB licence verdicts
  this note reuses, and for the shipped enzyme set

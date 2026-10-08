---
search:
  exclude: true
---

# Does an 80 °C hold inactivate an enzyme listed at 65 °C?

Research note for issue #522. `liulab_mbio.bench.inactivation.heat_inactivations` groups the
enzymes in one tube by the condition each supplier gives, so a tube holding BsaI-HFv2 and PmeI
is held twice: 80 °C for 20 minutes, then 65 °C for 20 minutes. Forty minutes, two stages, one
tube.

**The question.** Does the 80 °C hold already inactivate PmeI, making the 65 °C hold redundant
and the two collapsible into one?

**The answer is no source says so.** The condition each enzyme carries is the condition New
England Biolabs assayed, and the assay held a reaction at one temperature, not at both. NEB's
only statement about the two temperatures runs the other way: 80 °C is offered as a fallback for
enzymes that 65 °C fails to stop. Nothing published says the hotter hold covers the cooler
enzyme, so on present evidence the two holds stay.

A quoted row below is verbatim from the cited document. A row marked **no source** means the
document was searched and states nothing; it is not an invitation to supply a figure.

## 1. What each source gives for the two enzymes

| Source | PmeI (NEB #R0560) | BsaI-HFv2 (NEB #R3733) |
| --- | --- | --- |
| *Heat Inactivation* usage guideline, table | `65°C` / `20minutes` | `80°C` / `20minutes` |
| product page, **Heat Inactivation** field | `65°C for 20 minutes` | `80°C for 20 minutes` |
| product page, property icon | `heat inactivation 65°` | `heat inactivation 80°` |
| product specification sheet | **no source** — the field is absent | **no source** — the field is absent |
| NEBcutter 3 enzyme table | **no source** — no such field | **no source** — no such field |
| NEBridge Golden Gate kit manual | not applicable | **no source** — the word *inactivation* is absent |

The two readings agree with each other and with what `src/liulab_mbio/data/enzymes.json`
already ships, which `docs/research/restriction-enzyme-data.md` section 6 records as checked
live on 2026-09-18.

**One trap worth naming.** The guideline table carries three BsaI rows, and they disagree:
`BsaI-HF®` is `65°C` / `20minutes` while `BsaI-HF®v2` is `80°C` / `20minutes`. Read the v2 row.
Reading the wrong one makes the question disappear by making both enzymes 65 °C.

## 2. The crux: what NEB says about the two temperatures

The *Heat Inactivation* usage guideline opens with this, verbatim:

> Heat inactivation is a convenient method for stopping a restriction endonuclease reaction.
> Incubation at 65°C for 20 minutes inactivates the majority of restriction endonucleases that
> have an optimal incubation temperature of 37°C. Enzymes that cannot be inactivated at 65°C can
> often be inactivated by incubation at 80°C for 20 minutes. The table below indicates whether or
> not an enzyme can be heat inactivated and the temperature needed to do so.

This is the only sentence on any source read here that puts 65 °C and 80 °C in one claim, and it
points the opposite way from the one the question needs. It says what to try when 65 °C is not
enough. It does not say that 80 °C does everything 65 °C does, and an enzyme listed at 65 °C is
by that same sentence one of the majority that 65 °C already handles — so it was never a
candidate for the 80 °C fallback and NEB has no occasion to say what 80 °C does to it.

The same page then states how the table was measured, verbatim:

> Heat inactivation was performed as follows to approximate a typical experiment. A 50 µl
> reaction mixture containing the appropriate NEBuffer, 0.5 µg of calf thymus DNA, and 5 or 10 µl
> of restriction endonuclease (at selling concentration) was incubated at 37°C for 60 minutes and
> then at 65°C or 80°C for 20 minutes. 0.5 µg of substrate DNA (usually lambda) was added to the
> reaction mixture and incubated at the optimal reaction temperature of the enzyme for 60 minutes.
> Any digestion (complete or partial) of the substrate DNA after the second incubation, as seen by
> agarose gel electrophoresis, was interpreted as incomplete heat inactivation.

**"65°C or 80°C"** is the load-bearing word. Each enzyme's reaction was held at one of the two,
and the table reports that one result. PmeI's entry is therefore evidence that PmeI stops at
65 °C and evidence about nothing else. It is not a measurement at 80 °C, and no column of that
table reports one.

The wording is stable. An archived copy of the same page from 2021-04-20 carries both paragraphs
word for word, so this is not a recent edit that a later page might undo.

## 3. The other two things the brief asked for

**Is either enzyme flagged as not heat inactivatable?** No. Both carry a temperature. The
guideline does use `No` as a value — `CviKI-1` and `CviQI` read `No` / `--` in the 2026 copy,
`AclI`, `AlwI` and `ApaLI` in the 2021 one — and `docs/research/bench-numbers.md` records T7
DNA Ligase the same way, with the
further caution that its 65 °C option is void in a buffer carrying PEG. So `No` is an answer NEB
gives when it applies, which makes the temperatures it does give meaningful rather than
boilerplate.

**Does any source give a combined condition for two enzymes in one tube?** No. The *Double
Digests* usage guideline is the nearest, and it describes the sequential case only:

> If two different incubation temperatures are necessary, choose the optimal reaction buffer and
> set up reaction accordingly. Add the first enzyme and incubate at the desired temperature. Then,
> heat inactivate the first enzyme, add the second enzyme and incubate at the recommended
> temperature.

In that procedure the second enzyme is not in the tube during the heat step, so the page never
has to say what one hold does to two enzymes, and it names no temperature for the step at all.
*Optimizing Restriction Endonuclease Reactions* mentions heat inactivation in one line and links
straight back to the guideline. Both NEBridge Golden Gate Assembly Kit manuals, E1601 and E1602,
contain no occurrence of *inactiv* in any form.

## 4. Sources consulted

Read on **2026-10-08** unless the date column says otherwise.

| Source | URL | How read | Date |
| --- | --- | --- | --- |
| NEB, *Heat Inactivation* usage guideline | `https://www.neb.com/en-us/tools-and-resources/usage-guidelines/heat-inactivation` | held under `reference_docs/synthesis_and_assembly/bench/neb/`, fetched through `r.jina.ai` | 2026-10-06 |
| the same page, re-read | `https://web.archive.org/web/20260305053645/https://www.neb.com/en-us/tools-and-resources/usage-guidelines/heat-inactivation` | Wayback snapshot, plain `curl` | 2026-10-08 |
| the same page, five years earlier | archived copy under `reference_docs/restriction-ligation/web/archived/` | already held | 2021-04-20 |
| NEB PmeI (#R0560) product page | `https://www.neb.com/en-us/products/r0560-pmei` | held under `reference_docs/synthesis_and_assembly/bench/neb/`, fetched through `r.jina.ai` | 2026-10-06 |
| NEB BsaI-HFv2 (#R3733) product page | `https://web.archive.org/web/20260221183023/https://www.neb.com/en-us/products/r3733-bsai-hf-v2` | Wayback snapshot, plain `curl` | 2026-10-08 |
| NEB product specification PS-R0560S/L v1.0, 30 Jul 2013 | `https://www.neb.com/en/-/media/catalog/specifications/r/0/r0560s_l_v1.pdf` | PDF held under `reference_docs/`, with its text dump; URL re-checked by `curl` | 2026-10-08 |
| NEB product specification PS-R3733S/L v1.0, 13 Dec 2017 | `https://www.neb.com/en/-/media/catalog/specifications/r/3/r3733s_l_v1.pdf` | PDF held under `reference_docs/`, with its text dump; URL re-checked by `curl` | 2026-10-08 |
| NEBridge Golden Gate Assembly Kit manual E1601 | `https://www.neb.com/-/media/nebus/files/manuals/manuale1601.pdf` | PDF held under `reference_docs/`, with its text dump | 2026-10-06 |
| NEBridge Golden Gate Assembly Kit manual E1602 | `https://www.neb.com/-/media/nebus/files/manuals/manuale1602.pdf` | PDF held under `reference_docs/`, with its text dump | 2026-10-06 |
| NEB, *Double Digests* usage guideline | `https://www.neb.com/en-us/tools-and-resources/usage-guidelines/double-digests` | held under `reference_docs/restriction-ligation/web/live/` | 2026-10-06 |
| NEB, *Optimizing Restriction Endonuclease Reactions* | `https://www.neb.com/en-us/tools-and-resources/usage-guidelines/optimizing-restriction-endonuclease-reactions` | held under `reference_docs/restriction-ligation/web/live/` | 2026-10-06 |
| NEBcutter 3 enzyme table | `https://nc3.neb.com/NEBcutter/data/enzymes.json` | plain `curl` | 2026-10-08 |
| NEB restriction enzyme technical guide | held under `reference_docs/restriction-ligation/manuals/` | already held | 2026-10-06 |

The BsaI-HFv2 product page needed a fresh read because the copy held under `reference_docs/`
is a `Page Not Found` shell: it was fetched from `r3733-bsai-hfv2`, and the live slug is
`r3733-bsai-hf-v2`. The shell carries the full site navigation, so it greps like a real page
and reads like one until you check its title.

## 5. What was searched and did not answer the question

- **NEB's two product specification sheets.** Both were read end to end. They give
  concentration, unit definition, shelf life, storage and the release assays, and no heat
  inactivation condition at all. The product page is the only place NEB states one.
- **The NEBcutter 3 enzyme table.** It carries the recognition sequence, both cut offsets, the
  incubation temperature, the supplied buffer, the Time-Saver flag and methylation notes. There
  is no heat-inactivation key in the record for either enzyme.
- **Both Golden Gate kit manuals.** A kit manual that ran BsaI-HFv2 and a second enzyme in one
  tube would be the natural place for a combined hold. Neither manual mentions inactivation.
- **Two NEB FAQ pages that a web search surfaced**, under `neb.com/en-us/faqs/`: the slugs
  `how-should-i-stop-my-restriction-digest` and
  `are-certain-restriction-enzymes-more-active-at-higher-incubation-temperatures-than-what-is-recommended`.
  Neither has a Wayback snapshot, by its CDX query, and `neb.com` would not serve either today.
  Their text is therefore unverified and nothing from them is quoted here. The first is worth
  another attempt whenever `neb.com` opens again; by the search engine's own summary it points
  the reader back to per-enzyme catalogue information, which would strengthen this note's answer
  rather than change it.
- **A forum thread reporting an e-mail from NEB.** A Protocol Online post says NEB replied by
  e-mail that the higher of two temperatures should be used. It is a user post quoting private
  correspondence, unverifiable, and it is not a published statement. It is recorded here so the
  next person does not spend the search finding it again.

## 6. How `neb.com` was reachable today, and how it was not

`docs/research/bench-numbers.md` records the routes past NEB's HTTP 403. One of them has
closed since, so this note updates that record.

| Route | Result on 2026-10-08 |
| --- | --- |
| `web.archive.org` snapshots, plain `curl` | **works** — both pages this note depends on came through it |
| `nc3.neb.com/NEBcutter/data/enzymes.json`, plain `curl` | **works** |
| the specification and manual PDFs, plain `curl` | **works**, and the spelling matters: `/r/0/r0560s_l_v1.pdf` returns the PDF where `r0560.pdf`, `r0560s.pdf`, `ps-r0560.pdf` and four other spellings all redirect to the 403 error page. Take `<catalog>s_l_v1.pdf` as the pattern for a product sold in two pack sizes |
| `r.jina.ai/<neb url>` | **closed.** Returns HTTP 401, `AuthenticationRequiredError`: `You have been blocked from performing anonymous queries due to bad network reputation (AS7018). Please authenticate.` This is not the rate limit the earlier note describes, and waiting does not clear it |
| plain `curl` to `neb.com` with a browser user agent | HTTP 403, as before |
| `WebFetch` on a `neb.com` HTML page | HTTP 403, as before |
| `prd-sccd00`, `prd-sccd01`, `prd-sccd02.neb.com`, `www.neb.cn` | HTTP 403 — the mirrors a search engine lists are behind the same block |
| headless Chrome | a Cloudflare `Verifying you're human` interstitial that does not clear |

So the live HTML routes are all shut, and the archive is what is left. A snapshot is dated,
which the table in section 4 records per source; the 2021 copy agreeing word for word with the
2026 one is the reason that matters less here than it usually would.

## 7. What would settle it

Nothing in this note forbids collapsing the two holds. It records that no published source
permits it. Three things would:

- **NEB saying so in writing.** A reply from technical support, quoted and dated, is a source;
  a reply a forum reports is not.
- **A measurement.** The guideline's own assay is the protocol to copy: hold a PmeI reaction at
  80 °C for 20 minutes, add fresh substrate DNA, incubate at 37 °C for an hour, and run a gel.
  Any digestion means incomplete inactivation.
- **Finding the FAQ page.** Section 5 names the one page that might already answer this and
  could not be read today.

Until one of those lands, the two holds are what the suppliers' own conditions add up to, and
`heat_inactivations` grouping by condition is reporting that honestly rather than padding the
protocol.

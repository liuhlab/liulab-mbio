---
search:
  exclude: true
---

# Bench numbers for the generated protocol, stage by stage

Research note for issue #233, under #225, superseding #196. Everything below was retrieved on
**2026-10-06** unless a source carries its own date. It is the inventory of every number the
generated protocol must carry and cannot compute, for the five stages the ticket fixes.

This is the inventory, not the wiring. Nothing here is a package default yet.

## The rule this note is held to

Every number comes from a protocol, a kit insert or a cited source. **Nothing is invented and
nothing comes from memory.** Where no source gives a number, the row says so and the hole is
left visible: a hole is a finding, a guess is a defect. Vendor and catalogue number are named
wherever identity is load-bearing.

Two conventions:

- A value in a table is quoted from the cited document. Where a document gives a range, the
  range is reproduced, never narrowed.
- **Ours** in a source column means no document states it — the method page asks for the step
  and nobody published a number for it. Those rows are the holes, collected again at the end.

PaqCI was falsified as the destination enzyme by #201 and #207. The design uses BbsI, for the
reason in departure D2. Nothing in this note uses PaqCI or its activator oligo.

## Where each source came from

| Source | What it is | How it was read | Date |
| --- | --- | --- | --- |
| Zero Blunt TOPO PCR Cloning Kit user guide, Pub. MAN0000062 Rev B.0 (2014) | the part carrier's own kit insert | `assets.thermofisher.com/TFS-Assets/LSG/manuals/zeroblunttopo_man.pdf`, plain `curl` | 2026-10-06 |
| One Shot TOP10 Competent Cells user guide, Pub. MAN0000633 Rev A.0 | the cells that kit ships with | `assets.thermofisher.com/TFS-Assets/LSG/manuals/oneshottop10_man.pdf`, plain `curl` | 2026-10-06 |
| Qian et al. 2026, Supplementary Information | the DMX vector's own bench protocol, day by day | already held under `reference_docs/synthesis_and_assembly/dmx/paper/` | 2026-10-06 |
| Takacsi-Nagy et al. 2026, STAR Methods and Key Resources Table | the iGGA round, as the source ran it | already held under `reference_docs/synthesis_and_assembly/prot-assembly/` | 2026-10-06 |
| Endura Competent Cells manual, MA133 26Feb2018 | the electroporation the iGGA round uses | Biosearch/LGC, plain `curl` | 2026-10-06 |
| Agencourt AMPure XP instructions for use, B37419AB (2016) and B37419AA (2013) | the SPRI clean-up both pipelines call for | Beckman Coulter, plain `curl` | 2026-10-06 |
| Oxford Nanopore kit and flow cell pages | the read-out's consumables | `nanoporetech.com`, plain `curl`, text extracted | 2026-10-06 |
| New England Biolabs manuals and product specifications | every enzyme, buffer and clean-up kit | `neb.com/-/media/nebus/files/manuals/*.pdf` and `neb.com/en/-/media/catalog/specifications/<letter>/<digit>/*.pdf`, both plain `curl`; `nc3.neb.com/NEBcutter/data/enzymes.json` for the cross-check | 2026-10-06 |
| New England Biolabs product pages, protocol pages and usage guidelines | heat inactivation, the two transformation protocols, the ligase buffer and the electroporation settings | `r.jina.ai/<neb url>` | 2026-10-06 |

`neb.com` returns HTTP 403 to `curl` and WebFetch for its HTML pages. The routes that work are
recorded in `docs/research/gibson-assembly.md` and `docs/research/restriction-ligation.md`, and
all three were used here. A wrong filename and a block return the same 403, so a missing file
reads as a block — compare response sizes before concluding a manual is unreachable.

Two traps on the `r.jina.ai` route, paid for here:

- **It rate-limits.** Three requests in a row returned HTTP 422 and one returned a Cloudflare
  "Just a moment" page of 548 bytes. Five seconds between requests, and a retry that rejects a
  response under 2 kB, fetched twelve of twelve.
- **A protocol's dated URL renders an empty shell.** The page NEB's own product page links to
  as `/protocols/2016/06/08/high-efficiency-transformation-protocol-...` comes back with
  navigation and no protocol. The undated slug, `/protocols/high-efficiency-transformation-
  protocol-c3040h`, returns the steps. Take the slug from the product page's protocol list.

---

## Stage 1 — Part cloning into pCR-Blunt II-TOPO, and the tailed PCR that retailors an overhang

The kit insert is the whole source for this stage. It is the best-sourced of the five.

### The blunt PCR product that goes in

| Number | Value | Source |
| --- | --- | --- |
| PCR reaction volume | 25 or 50 µL | Zero Blunt TOPO UG p. 10 |
| Polymerase | "a thermostable proofreading polymerase" — no catalogue number is required | Zero Blunt TOPO UG p. 10 |
| Polymerase the guide names as an example | Platinum Pfx (11708-013/-021/-039); AccuPrime Pfx (12344-024/-032); Pfx50 (12355-012) appear only in the ordering table | Zero Blunt TOPO UG pp. 10, 26 |
| Final extension | 7-30 minutes | Zero Blunt TOPO UG p. 10 |
| Primer 5' phosphate | must **not** be present; a phosphorylated product will not ligate | Zero Blunt TOPO UG p. 10 |
| Hold after cycling | ice, or -20 °C for up to 2 weeks | Zero Blunt TOPO UG p. 10 |
| Taq product rescue | remove the 3' A with a proofreading polymerase or T4 DNA polymerase plus dNTPs | Zero Blunt TOPO UG p. 21 |
| Insert size ceiling | over 3 kb, the guide sends you to the TOPO XL kit instead | Zero Blunt TOPO UG p. 21 |
| Per-cycle denature, anneal and extend times | **ours** — the guide defers to the polymerase | — |

The retailoring primer's own structure, `5'-[spacer]-CGTCTC-N-[new overhang]-[annealing]-3'`,
is the method page's and is designed per part, not sourced. Its Ta is the package's to compute
from `liulab_mbio.primers.polymerase`, once a polymerase is named.

### The TOPO reaction

| Number | Value | Source |
| --- | --- | --- |
| Fresh PCR product | 0.5-4 µL | Zero Blunt TOPO UG p. 12 |
| Salt Solution | 1 µL (stock 1.2 M NaCl, 0.06 M MgCl₂) | Zero Blunt TOPO UG pp. 6, 12 |
| Water | to 5 µL | Zero Blunt TOPO UG p. 12 |
| pCR-Blunt II-TOPO vector | 1 µL at 10 ng/µL | Zero Blunt TOPO UG pp. 6, 12 |
| Final volume | 6 µL | Zero Blunt TOPO UG p. 12 |
| Final salt in the reaction | 200 mM NaCl, 10 mM MgCl₂ | Zero Blunt TOPO UG p. 11 |
| Incubation | 5 minutes at room temperature, then on ice | Zero Blunt TOPO UG p. 12 |
| Incubation range allowed | 30 seconds to 30 minutes; 30 s suffices for routine subcloning, longer helps a large product | Zero Blunt TOPO UG p. 12 |
| What the salt buys | 2- to 3-fold more transformants, and tolerance of the longer incubation | Zero Blunt TOPO UG p. 11 |
| Storage of the reaction | -20 °C overnight | Zero Blunt TOPO UG p. 12 |
| Insert-to-vector molar ratio | **ours** — the guide gives volumes only, never a ratio, a ng of insert or a pmol | — |

### Transformation, selection and screening

| Number | Value | Source |
| --- | --- | --- |
| Reaction into chemically competent cells | 2 µL into one 50 µL vial, mixed gently, never by pipetting | Zero Blunt TOPO UG p. 14 |
| Ice | 5-30 minutes | Zero Blunt TOPO UG p. 14 |
| Heat shock | 30 seconds at 42 °C, no shaking, then straight to ice | Zero Blunt TOPO UG p. 14 |
| Outgrowth | 250 µL S.O.C., 37 °C, 200 rpm horizontal, 1 hour | Zero Blunt TOPO UG p. 14 |
| Plating | 10-50 µL on a prewarmed plate, two volumes, plus 20 µL S.O.C. to spread | Zero Blunt TOPO UG p. 14 |
| Plate prewarming | 37 °C for 30 minutes | Zero Blunt TOPO UG p. 13 |
| Plate incubation | overnight at 37 °C | Zero Blunt TOPO UG p. 14 |
| Selection | LB + 50 µg/mL kanamycin, **or** Low Salt LB + 25 µg/mL Zeocin | Zero Blunt TOPO UG p. 13 |
| Low Salt LB | 1% tryptone, 0.5% yeast extract, 0.5% NaCl, pH 7.5 | Zero Blunt TOPO UG p. 24 |
| Zeocin trap | regular LB kills Zeocin selection; Low Salt LB is required | Zero Blunt TOPO UG p. 13 |
| Expected colonies | "several hundred" from an efficient reaction | Zero Blunt TOPO UG p. 14 |
| Expected cloning efficiency | at least 95%, less for a large insert; the kit control gives >100 colonies, 95% carrying the 800 bp insert, and <5% background from vector alone | Zero Blunt TOPO UG pp. 20-21 |
| Colonies to pick | about 10 | Zero Blunt TOPO UG p. 14 |
| Screening culture | 2-6 colonies overnight in LB with the same antibiotic; patch the original colony | Zero Blunt TOPO UG p. 16 |
| Sequencing primers | M13 Forward (-20) `GTAAAACGACGGCCAG`, M13 Reverse `CAGGAAACAGCTATGAC`, supplied at 0.1 µg/µL | Zero Blunt TOPO UG p. 6 |
| Colony PCR | 1 µL each of 20 µM primers in 48 µL PCR SuperMix High Fidelity; 10 colonies into 50 µL; 94 °C 10 min; 20-30 cycles; 72 °C 10 min final | Zero Blunt TOPO UG pp. 16-17 |
| Colony PCR per-cycle times | **ours** — the guide gives the cycle count and the two bookends, and no denature, anneal or extend | — |
| Glycerol stock | grow to OD600 ≈ 0.5 in 1-2 mL selective LB; 0.85 mL culture + 0.15 mL sterile glycerol; -80 °C | Zero Blunt TOPO UG p. 17 |

Electroporation of the same reaction is an alternative the guide gives: dilute 6 µL with 18 µL
water, put 2 µL into a 50 µL vial, 0.1 cm cuvette, 250 µL S.O.C., shake at least 1 hour at
37 °C (UG p. 14). The dilution exists to bring the salt to 50 mM NaCl and 2.5 mM MgCl₂ and
stop arcing (UG p. 15). The electroporator settings are the user's own — the guide states none.

**Counter-selection.** The vector carries `ccdB` fused to the C-terminus of `lacZα`. A blunt
insert disrupts the fusion, so only recombinants grow and no blue/white screen is needed
(UG p. 8, citing Bernard et al. 1994).

### Catalogue

K2800-20 (25 reactions), K2800-40 (50), K2800-J10 (10), all with One Shot TOP10 chemically
competent cells; K2820-20/-40 with DH5α-T1R; K2830-20 with Mach1-T1R; K2860-20/-40 with TOP10
electrocompetent cells; K2800-02 adds the PureLink miniprep kit; 450245 ships the reagents
without cells (UG p. 5). Box 1 stores at -30 to -10 °C in a non-frost-free freezer; Box 2, the
cells, at -85 to -68 °C (UG pp. 5-6).

---

## Stage 2 — Working cassette assembly

The BsmBI digest of the input working vector, the gel purification, the one-pot assembly, the
white-colony pick and NEB Stable.

**This is the thinnest-sourced stage.** No published protocol runs this exact step: it is the
method page's own, assembled from the DMX vector's GGA protocol and NEB's kit. The enzyme and
clean-up numbers come from NEB; everything about the RFP stuffer and the white-colony pick is
ours.

### The digest and the gel

| Number | Value | Source |
| --- | --- | --- |
| Enzyme identity, concentration, buffer and temperature | see **The five enzymes** below | NEB product specifications |
| Mass of input working vector to digest | **ours** — no source digests this vector | — |
| Digest volume and time for a preparative digest | **ours** — the method page calls for a gel-purified backbone and nobody published the scale | — |
| Gel percentage for separating the cut backbone from the RFP stuffer | **ours** — it depends on the stuffer's length, which is "long enough to separate" and not a number | — |

Gel extraction has two sourced options. Zymo's is the one Qian runs; NEB's is the one the
method page's `## Reagents and equipment` implies by naming a column clean-up kit.

| Number | Value | Source |
| --- | --- | --- |
| Zymoclean Gel DNA Recovery Kit | Zymo #D4007/D4008, eluted into 8 µL (Day 1) or 9 µL (Day 4) | Qian SI Day 1.1 and Day 4.1 |
| Monarch Gel Extraction, NEB #T1020S/L | 5 µg binding capacity; ~50 bp to 25 kb; 70-90% recovery from 50 bp to 10 kb, 50-70% from 11 to 23 kb; A260/280 ≥1.8 | T1020 manual v2.1_4/21 p. 2 |
| T1020 dissolving | 4 volumes Gel Dissolving Buffer per mg of gel read as µL — 400 µL per 100 mg; 3-3.5 volumes if the slice is over 150 mg | T1020 manual p. 5 |
| T1020 dissolve conditions | 37-55 °C, typically 50 °C, 5-10 minutes | same |
| T1020 for fragments over 8 kb | add 1.5 volumes of water before loading | same |
| T1020 spins and washes | 16,000 × g for 1 minute (30 s allowed); two 200 µL washes | same |
| T1020 elution | ≥6 µL, typically 6-20 µL; wait 1 minute, spin 1 minute; warm the buffer to 50 °C for ≥10 kb | same |
| Monarch PCR & DNA Cleanup, NEB #T1030S/L | same capacity, size range and recovery; sample 20-100 µL | T1030 manual v3.0_4/21 pp. 2, 5 |
| T1030 binding ratios | 2:1 buffer to sample for dsDNA over 2 kb such as a plasmid; **5:1 under 2 kb**; 7:1 for ssDNA | T1030 manual p. 5 |
| Expected yield from the gel | ~100-200 ng per 25 µL PCR, for Qian's library amplicon — **not** for a preparative vector digest | Qian SI Day 1.1 |

### The one-pot assembly

Qian's own GGA recipes are the closest published source, and they use BsaI-HFv2 rather than
BsmBI. Both are quoted; neither is this step.

| Number | Value | Source |
| --- | --- | --- |
| Vector-to-insert molar ratio | 1:2 | Qian SI, "GGA cloning" |
| 1 µL reaction (acoustic) | 0.1 µL 10x T4 Ligase Buffer, 1.2 U BsaI-HFv2, 40 U T4 DNA Ligase, 2-4 fmol vector, 4-8 fmol insert, water to 1 µL | Qian SI, "GGA cloning" |
| 5 µL reaction (by hand) | 0.5 µL 10x T4 Ligase Buffer, 3 U BsaI-HFv2, 100 U T4 DNA Ligase, 20-30 fmol vector, 40-60 fmol insert, water to 5 µL | Qian SI, "GGA cloning" |
| Cycling for these two | 37 °C 15 min, then 60 °C 5 min "to reduce background" | Qian SI, "GGA cloning" |
| Plate handling | spin 1,000 × g for 30 s before incubation | Qian SI, "GGA cloning" |
| The same reaction with BsmBI instead of BsaI | not stated by Qian; their BsmBI reaction is the 40 µL library cloning in stage 3 | — |
| Vector stock a midiprep yields | 50-100 µg, enough for 2,400-4,800 one-µL reactions; 2 µg covers one 96-well plate | Qian SI, "Preparing GGA vector stocks" |

NEB's own kit is the other sourced route, and it is the one with a cycled program. The BsmBI
kit is **NEB #E1602**; the BsaI kit, which stage 5 uses, is **#E1601S/L**.

| Number | Value | Source |
| --- | --- | --- |
| Reaction, 20 µL | destination 0.05 pmol (75 ng of pGGAselect); each precloned insert 0.05 pmol; an amplicon insert 1:1 molar to the vector; 2 µL 10x T4 DNA Ligase Buffer; 1 µL NEBridge Golden Gate Enzyme Mix for ≤10 inserts, 2 µL above that; water to 20 µL | E1601 manual v5.0_6/26 p. 6; E1602 manual v3.0_6/26 |
| Scaling to 25 µL | allowed, with 0.5 µL more buffer | E1601 manual p. 6 |
| E1602 (BsmBI) cycling, 1 insert | 42 °C 5 min for cloning, or 42 °C 1 h for library prep; then 60 °C 5 min | E1602 manual |
| E1602 cycling, 2-10 inserts | (42 °C 1 min, 16 °C 1 min) × 30-60; then 60 °C 5 min | E1602 manual |
| E1602 cycling, 11-20+ inserts | (42 °C 5 min, 16 °C 5 min) × 30-60; then 60 °C 5 min | E1602 manual |
| What the final 60 °C step does | digests destination that stayed uncut or religated | E1601 manual, FAQ 10 |
| Cycle count | efficiency rises from 30 to 60-65 cycles with no loss of fidelity | E1601 manual, FAQ 14 |
| Scale-down | 2- to 3-fold, except for a very complex assembly | E1601 manual, FAQ 12 |
| Insert tiers | the current manuals band at 1, 2-10 and 11-20+ inserts | E1601 manual p. 6 |

### The pick, and NEB Stable

| Number | Value | Source |
| --- | --- | --- |
| Strain for a stock whose ccdB is expressed | DB3.1 or ccdB Survival 2 T1R | the method page. Qian SI, "GGA cloning", names NEB Stable, which is ccdB-sensitive; it works there because DMX0001's ccdB is T7-silent — #300 |
| NEB Stable culture | pick a single clone into 50 mL LB — **not** a rich medium such as Terrific Broth — in a 250 mL baffled flask, overnight | Qian SI, "Preparing GGA vector stocks" |
| Plasmid prep | midiprep; typical yield 50-100 µg | Qian SI, "Preparing GGA vector stocks" |
| NEB Stable, #C3040H/I | >1 × 10⁹ cfu/µg, measured on 50 µL of cells with 100 pg pUC19 on LB-amp at 37 °C; stored at -80 °C; 12-month shelf life | C3040 product specification PS-C3040H/I v1.0 |
| NEB Stable thaw | 10 minutes on ice | C3040H high-efficiency transformation protocol |
| DNA in | 1-5 µL containing 1 pg to 100 ng of plasmid; flick 4-5 times, never vortex | same |
| Ice | 30 minutes, without mixing | same |
| Heat shock | **exactly** 42 °C for **exactly** 30 seconds, then 5 minutes on ice, no mixing | same |
| Outgrowth | 950 µL room-temperature NEB 10-beta/Stable Outgrowth Medium (#B9035), **30 °C for 60 minutes**, horizontal at 250 rpm | same |
| Plating | 50-100 µL of cells or diluted cells; **24 h at 30 °C**, or overnight at 37 °C | same |
| When outgrowth can be skipped | an AmpR plasmid needs none; any other selection needs the 60 minutes at 30 °C | C3040 product page, product note 3 |
| **Why 30 °C** | "30 °C or 37 °C may be used for plate incubation, however 30 °C is recommended as some constructs may be unstable at elevated temperatures" | C3040 product page, product note 3 |
| Storage | -80 °C; -20 °C costs a significant amount of efficiency, and the cells lose efficiency whenever warmed above -80 °C even without thawing | C3040 product page, product note 1 |
| The kit's own transformation, as a comparable | 50 µL 10-beta cells thawed on ice 10 min; 2 µL of assembly; ice 30 min; 42 °C 30 s; ice 5 min; 950 µL room-temperature NEB 10-beta/Stable Outgrowth Medium; 37 °C 60 min at 250 rpm | E1601 manual p. 6 |
| The kit's own plating | 50 µL of a 1:5 dilution for a single insert, or 50-100 µL for a multi-insert assembly; overnight at 37 °C, **or 24-36 h at 30 °C**, or 48 h at 25 °C | E1601 manual p. 6 |
| What a white colony means here | **ours** — the RFP stuffer and the "a colony that lost it is the one to pick" rule are the method page's. No source gives a false-positive rate, a colonies-to-pick count, or an expected white fraction | — |
| Transformation efficiency to expect for this assembly | **ours.** The nearest sourced figure is the kit's own specification test: 75 ng pGGAselect plus 75 ng each of five plasmids, 2 µL into T7 Express, gives >250 colonies and >80% blue | E1601 manual p. 9 |
| Whether the carrier's own BsaI site and ccdB harm a part added as the carrier plasmid | **unresolved in the method itself** — `docs/research/synthesis-and-assembly.md` already carries this as undecided | — |

---

## Stage 3 — Cargo into the DMX vector, plate barcoding, and the ONT read-out

Qian et al. 2026's supplementary information is a day-by-day bench protocol and it covers this
stage end to end. It is the best-sourced stage after stage 1.

### Day 1.1 — library amplification from the oligo pool

| Number | Value | Source |
| --- | --- | --- |
| Reaction volume | 25 µL | Qian SI Day 1.1 |
| Master mix | 12.5 µL 2x KAPA HiFi HotStart Ready Mix | Qian SI Day 1.1 |
| Dye | 1.25 µL 20x EvaGreen | Qian SI Day 1.1 |
| Each primer | 0.75 µL at 10 µM | Qian SI Day 1.1 |
| Template | 2.5 ng oligo library DNA | Qian SI Day 1.1 |
| Cycling | 95 °C 3 min; **14 cycles** of 98 °C 20 s / 65 °C 15 s / 72 °C 40 s; 72 °C 1 min; 4 °C hold | Qian SI Day 1.1 |
| Gel | 2% agarose, band cut out | Qian SI Day 1.1 |
| Extraction | Zymoclean Gel DNA Recovery Kit (Zymo #D4007/D4008), eluted into 8 µL | Qian SI Day 1.1 |
| Quantification | 1 µL on Qubit dsDNA HS (Invitrogen #Q32851) | Qian SI Day 1.1 |
| Expected yield | ~100-200 ng per PCR reaction | Qian SI Day 1.1 |

The PCR1/PCR2 cycling of our own three-primer scheme is **not** this reaction. Qian amplifies a
whole library in one pass with the primer pair in Supplementary Table 4; the three-primer scheme
replaces it with a per-batch and a per-gene pair. **No source gives cycling for PCR1 or PCR2** —
`docs/research/synthesis-and-assembly.md` already records that nothing documents the scheme in
text.

### Day 1.2 — cloning into the DMX vector

| Number | Value | Source |
| --- | --- | --- |
| Reaction volume | 40 µL | Qian SI Day 1.2 |
| Enzyme | 40 units, BsmBI-v2 NEBridge Golden Gate Assembly kit | Qian SI Day 1.2 |
| Buffer | 4 µL 10x T4 DNA ligase buffer | Qian SI Day 1.2 |
| Vector | 0.1 pmol total | Qian SI Day 1.2 |
| Insert | 0.3 pmol total — a 3:1 library-to-vector molar ratio | Qian SI Day 1.2 |
| Cycling | 42 °C 60 min; 60 °C 5 min; 4 °C hold | Qian SI Day 1.2 |
| Clean-up | DNA Clean & Concentrator-5 (Zymo #D4003), eluted in 10 µL water, Qubit-quantified | Qian SI Day 1.2 |
| Electroporation input | 90 ng of cleaned GGA product into 50 µL E. cloni EXPRESS BL21(DE3) electrocompetent cells (Lucigen) | Qian SI Day 1.2 |
| Electroporator settings | "according to the manufacturer's protocol" — **not stated** | Qian SI Day 1.2 |
| Expected time constant | ~4.1 ms | Qian SI Day 1.2 |
| Recovery | 1 mL Recovery Medium, 37 °C, 1 hour, 250 rpm | Qian SI Day 1.2 |
| Titre plating | serial dilutions 1:12,500, 1:25,000, 1:50,000 on 10 cm LB-agar with 100 µg/mL carbenicillin, overnight at 37 °C | Qian SI Day 1.2 |
| What happens to the rest | stored at 4 °C overnight and plated the next day | Qian SI Day 1.2 |

### Day 2 — coverage and the colony plate

| Number | Value | Source |
| --- | --- | --- |
| Coverage rule | electroporation efficiency (CFU × library complexity) should be **>300** | Qian SI Day 2 |
| Colony density on a 25 cm BioAssay plate (Corning #431111) | ~2,500 colonies is optimal for picking | Qian SI Day 2 |
| Plating volume | the CFU-derived volume made up to 500 µL with SOC, spread with beads | Qian SI Day 2 |
| Plate selection and incubation | 100 µg/mL carbenicillin, overnight at 37 °C | Qian SI Day 2 |
| Home-made plate drying | ~30 minutes, to stop colonies smearing | Qian SI Day 2 |

**The drug does not carry over.** Qian's DMX plasmid was AmpR. Departure D11 rebuilds ours KanR
before first use, so carbenicillin on it selects nothing. The plate is **50 µg/mL kanamycin**
(Zero Blunt TOPO UG p. 13, the guide covering this method's own part carrier) — a cited
transfer, not a new measurement. The temperature, the colony density and the plating volume are
properties of the step, not of the plasmid, and carry over unchanged. The package quotes
neither drug: it reads the marker off the record and names the drug that marker selects, which
ADR 0016 settles.

### Day 3 — picking

| Number | Value | Source |
| --- | --- | --- |
| Picker | QPix XE Microbial Colony Picker (Molecular Devices) | Qian SI Day 3 |
| Destination | ECHO-qualified 384-well plates (Beckman Coulter #c74290) | Qian SI Day 3 |
| Medium per well | 60 µL low-salt LB + 100 µg/mL carbenicillin | Qian SI Day 3 |
| Medium per well, ours | 60 µL low-salt LB + the vector's own drug, 50 µg/mL kanamycin for the KanR backbone D11 gives it | Qian SI Day 3 for the volume and the medium; Zero Blunt TOPO UG p. 13 for the drug |
| Low-salt LB | 10 g tryptone, 5 g yeast extract, 0.5 g NaCl per litre | Qian SI Day 3 |
| Seal | Breathe Easier seals (Fisher #NC1664397) | Qian SI Day 3 |
| Growth | 37 °C overnight, shaking at 1,000 rpm | Qian SI Day 3 |
| Manual-picking ceiling | below 500 genes, hand-picking is offered as the alternative | Qian SI Day 3 |

### Day 4.1 — barcoding in lysate

| Number | Value | Source |
| --- | --- | --- |
| Transfer | 2 µL of culture per well, four 384-well plates compressed into one 1536-well plate (Greiner #782270) | Qian SI Day 4.1 |
| Instrument | ECHO 525 acoustic liquid handler (Beckman Coulter) | Qian SI Day 4.1 |
| Pre-transfer step | invert the 384-well plates for 30 minutes so cells gather at the meniscus | Qian SI Day 4.1 |
| Seal | Axygen foil (Fisher #PCR-AS-600), then a vacuum-sealed bag | Qian SI Day 4.1 |
| Lysis | 98 °C water bath, 30 minutes | Qian SI Day 4.1 |
| Barcode input | 2 ng of each of the four UMIs, 8 ng total, in 0.5 µL | Qian SI Day 4.1 |
| Water transfer | 0.5 µL minus the barcode volume | Qian SI Day 4.1 |
| GGA mastermix per well | 0.4 µL 10x T4 Ligase Buffer, 1.2 U BsaI-HFv2, 40 U Salt-T4 DNA Ligase, water to 0.5 µL | Qian SI Day 4.1 |
| Well composition | 2 µL lysate + 0.5 µL barcodes + 0.5 µL water + 0.5 µL mastermix = **3.5 µL** | Qian SI Day 4.1 |
| Cycling | 37 °C 60 min; 60 °C 5 min; 4 °C hold | Qian SI Day 4.1 |
| Pooling | invert the plate and spin at 200 × g into a reservoir | Qian SI Day 4.1 |
| Pool clean-up | Zyppy Plasmid Miniprep (Zymo #D4036), eluted in 20 µL; **two columns per 1536-well plate**, because one saturates | Qian SI Day 4.1 |

The mastermix volumes sum to 0.4 + enzyme + enzyme and are then filled to 0.5 µL, so the
buffer is 0.4 µL of a 0.5 µL mix — at 3.5 µL total that is a 1.14x buffer, not 1x. The source
states it as written and this note does not correct it.

### Day 4.1 — the three amplicon reactions

| Number | Value | Source |
| --- | --- | --- |
| Reaction volume | 25 µL, one per primer pair (DMX1/2, DMX3/4, DMX5/6), run separately | Qian SI Day 4.1 |
| Master mix | 12.5 µL 2x KAPA HiFi HotStart Ready Mix | Qian SI Day 4.1 |
| Each primer | 1.25 µL at 10 µM | Qian SI Day 4.1 |
| Template | 20 ng purified DNA | Qian SI Day 4.1 |
| Cycling | 95 °C 3 min; **10 cycles** of 98 °C 20 s / 65 °C 15 s / 72 °C 40 s; 72 °C 1 min; 4 °C hold | Qian SI Day 4.1 |
| Gel | 1% agarose; bands extracted with Zymoclean D4007/D4008, eluted in 9 µL, then pooled | Qian SI Day 4.1 |
| Quantification | Qubit dsDNA HS (#Q32851) | Qian SI Day 4.1 |

Note the change of gel percentage: 2% for the library amplicon on Day 1, 1% for the barcoded
amplicons on Day 4. Both are the source's.

### Day 4.2 — the ONT read-out

| Number | Value | Source |
| --- | --- | --- |
| Library input | 200 fmol of the pooled amplicons | Qian SI Day 4.2 |
| Kit | Ligation Sequencing Kit V14, SQK-LSK114 (ONT) | Qian SI Day 4.2 |
| Flow cell | MinION R10.4.1, FLO-MIN114 (ONT) | Qian SI Day 4.2 |
| Small-library alternative | a Flongle flow cell, FLO-FLG114 | Qian SI Day 4.2 |
| Everything between input and loading | "according to the manufacturer's protocol" — end prep, ligation, bead ratios, loading mass and priming volumes are **not restated by Qian** | Qian SI Day 4.2 |
| Basecalling | Dorado 0.9.1, super-accuracy model; samtools 1.21, chopper 0.9.0, nanoq 0.10.0, cutadapt 4.9, minimap2 2.28, python 3.11.3 | Qian SI, "Generation of consensus sequences" |
| Reads per well to call a pass, Route A | consensus called only where depth is above 150 — Qian's, a consensus-calling cutoff that the method page reads as the floor ("What counts as a pass") | Qian SI, "Generation of consensus sequences" |
| Reads per well to call a pass, Route B | above 20 wanted, above 10 tolerable | LevSeq SI checklist; `route-b-index-pcr.md` section 6 |
| Reads per flow cell to budget | **ours** — no public ONT page gives an output figure for a flow cell | — |

What Qian defers to is ONT's own protocol. **Our route does not need ONT's barcodes** — the
barcode is in the construct — so the kit is SQK-LSK114, not a barcoding kit. The nearest
public ONT protocol is the native-barcoding amplicon one, `SQK-NBD114.24`
(NBA_9168_v114_revU, 15 Jul 2026). Its end prep, adapter ligation, clean-ups and loading are
reproduced below because they are the only published numbers for an amplicon library of this
shape. **They are not Qian's kit.** The barcode-ligation rows are the part that does not apply.

| Number | Value | Source |
| --- | --- | --- |
| Amplicon input | 200 fmol per sample, which is 130 ng for a 1 kb amplicon | NBD114.24 protocol §2 |
| End prep | 11.5 µL amplicon + 1 µL diluted DCS + 1.75 µL Ultra II buffer + 0.75 µL enzyme = 15 µL | NBD114.24 protocol §3 |
| End-prep cycling | 20 °C 5 min, then 65 °C 5 min | NBD114.24 protocol §3 |
| End-prep clean-up | 15 µL AXP beads (1x), 5 min on a Hula mixer, two 200 µL 80% ethanol washes, 30 s dry, elute in 10 µL water for 2 min | NBD114.24 protocol §3 |
| Adapter ligation | 30 µL pool + 5 µL NA + 10 µL 5x Quick Ligation Buffer + 5 µL Quick T4 ligase = 50 µL, 20 min at room temperature | NBD114.24 protocol §5 |
| Adapter clean-up | 20 µL AXP (0.4x), 10 min Hula, two 125 µL washes with SFB (all sizes) or LFB (enriches ≥3 kb), elute in 15 µL EB, 10 min at 37 °C. **Ethanol must not be used here** | NBD114.24 protocol §5 |
| Loading mass | <1 kb, 100 fmol; 1-10 kb, 35-50 fmol; >10 kb, 300 ng. Below that, load the whole library | NBD114.24 protocol §5 |
| Priming mix | 1,170 µL FCF + 5 µL BSA at 50 mg/mL + 30 µL FCT = 1,205 µL | NBD114.24 protocol §6 |
| Loading | 800 µL through the priming port, wait 5 min, then 200 µL; library mix 37.5 µL SB + 25.5 µL LIB + 12 µL library = 75 µL, dropwise into SpotON | NBD114.24 protocol §6 |
| Flow cell warranty floor | 800 active pores for MinION/GridION; a Flongle is replaced below 50 | NBD114.24 protocol §2; ONT Flongle store page |
| Hands-on time | 140 minutes: end prep 20, barcode ligation 60, adapter ligation and clean-up 50, priming and loading 10 | NBD114.24 protocol §1 |
| Run length and expected output | **not stated on any public ONT page read.** The protocol says only to start MinKNOW at default settings | — |

The method page's Route B, index PCR, is Long et al. 2025's. Its plate maps and primer plates
are already inventoried in `docs/research/synthesis-and-assembly-materials.md`; its per-well
PCR numbers are in `docs/research/route-b-index-pcr.md`.

---

## Stage 4 — The iGGA round

Takacsi-Nagy et al. 2026's STAR Methods states this round as one paragraph, with catalogue
numbers in the Key Resources Table. It is the only published source for the split digest.

**One difference from our design, and it is already a recorded departure.** The source digests
the **donor** with BsaI then SrfI and the **destination** with BbsI and PmeI. Ours is the other
way round: destination BbsI + SrfI, donor BsaI + PmeI, because departure D5 moved PmeI from
each synthesised block's external stuffer into the vector. The volumes and times below carry
over unchanged; which tube holds which blunt cutter does not.

### The split digest

| Number | Value | Source |
| --- | --- | --- |
| DNA per digest | 1 µg of plasmid pool | Takacsi-Nagy STAR Methods, "Pooled Cloning — Comprehensive Libraries" |
| First enzyme | 2.5 µL BsaI-HFv2 (NEB #R3733L) | same; Key Resources Table |
| Volume | 50 µL total | same |
| Buffer | CutSmart (NEB) | same |
| First incubation | 1 hour at 37 °C | same |
| Second enzyme | 2.5 µL SrfI (NEB #R0629L) added to the same tube | same |
| Second incubation | a further 1 hour | same |
| The other tube | BbsI-HF (NEB #R3539L) and PmeI (NEB #R0560L), "using the same protocol" | same |
| Units per µL for any of the four | **not stated by the paper** — it gives volumes, never units. NEB's specifications below convert them | — |
| Heat inactivation | **not stated** — the source goes straight to SPRI | — |

### Pooling the part list

The step before the digest, and every number it states is the digest's own read backwards. The
pool's **total** is the 1 µg above; its members are pooled in **equal picomoles**, so each gives
that total over the member count, weighed at its own length; and its **concentration floor** is
what the digest leaves room for.

| Number | Value | Source |
| --- | --- | --- |
| Resuspension buffer | nuclease-free TE pH 8.0, or 10 mM Tris-HCl pH 8.0 | Twist, *How should Multiplexed Gene Fragments be resuspended?* |
| Resuspension, IDT's own | IDTE or molecular-grade water, vortex, 50 °C for 15-20 min, then verify | IDT, gBlocks resuspension |
| Vendor concentration | **a floor, not a value**: "at least 10 ng/µL is recommended for the stock dilution" (Twist); 10 ng/µL (IDT) | both, above |
| The pool's own floor | **computed, not quoted**: 1 µg over the 45 µL the 50 µL digest leaves after its two 2.5 µL enzymes, so **≥ 22.3 ng/µL** | the split-digest row above |

**The vendor default does not reach it, and that is the trap.** At 10 ng/µL an equimolar 1,000 ng
pool occupies 100 µL, more than twice the 45 µL the digest has for it. Following the vendor
figure literally makes the next step impossible, so the number the pool is brought to is set by
the digest and never typed into a step: `liulab_synbio.igga.bench.pool_floor_ng_ul` computes it
from the three constants already in that row, and the step states what it returns.

The two digests never share a tube. That is the method page's rule and the source's practice.

### The five enzymes

Every concentration, unit definition, buffer and temperature below is from that enzyme's NEB
product specification sheet. These are the numbers that turn "2.5 µL of BsaI-HFv2" into a unit
count, and they are shared with stages 2 and 5.

| Enzyme | Catalogue | Concentration | Unit definition | Buffer | Temp |
| --- | --- | --- | --- | --- | --- |
| BsaI-HFv2 | #R3733S/L | 20,000 U/mL | 1 µg pXba, 1 h, 37 °C, 50 µL | rCutSmart | 37 °C |
| BsmBI-v2 | #R0739S/L | 10,000 U/mL | 1 µg lambda, 1 h, 55 °C, 50 µL | NEBuffer r3.1 | 55 °C |
| BbsI-HF | #R3539S/L | 20,000 U/mL | 1 µg lambda, 1 h, 37 °C, 50 µL | rCutSmart | 37 °C |
| SrfI | #R0629S/L | 20,000 U/mL | 1 µg pNEB193-SrfI, 1 h, 37 °C, 50 µL | rCutSmart | 37 °C |
| PmeI | #R0560S/L | 10,000 U/mL | 1 µg lambda, 1 h, 37 °C, 50 µL | rCutSmart | 37 °C |

At those concentrations, Takacsi-Nagy's 2.5 µL is 50 units of BsaI-HFv2 or BbsI-HF and 25
units of SrfI or PmeI, in 50 µL on 1 µg of DNA. **That multiplication is this note's, not a
source's**, and it holds only if the lot matches the specification.

Two blunt-cutter figures bear on the design directly:

| Number | Value | Source |
| --- | --- | --- |
| SrfI re-ligation | after 20-fold over-digestion, ~75% of fragments ligate and >95% of those recut | R0629 specification |
| PmeI re-ligation | pNEB193 linearised with a 10-fold excess and religated gives <1% white colonies | R0560 specification |

Both say a blunt end these enzymes leave **can** religate under a DNA ligase. The method's
escape route is closed by the ligase refusing blunt ends, not by the cut site resisting them —
and that refusal is sourced under
**Does T7 DNA ligase refuse blunt ends?** below.

| Enzyme | Heat inactivation | Source |
| --- | --- | --- |
| BsaI-HFv2 | 80 °C for 20 minutes | NEB heat-inactivation usage guideline |
| BsmBI-v2 | 80 °C for 20 minutes | same |
| BbsI-HF | 65 °C for 20 minutes | same, and the R3539 product page |
| SrfI | 65 °C for 20 minutes | same, and the R0629 product page |
| PmeI | 65 °C for 20 minutes | same, and the R0560 product page |
| T7 DNA Ligase | **No** on the product page's property icons. Note 4 adds that 65 °C for 10 minutes inactivates it, **but only in a buffer without PEG** — "do not heat inactivate if there is PEG in the reaction buffer, as transformation will be inhibited" | M0318 product page |

The last row is the reason the iGGA round goes straight from ligation to SPRI: the ligation is
in StickTogether buffer, which carries PEG 6000, so heat inactivation is ruled out. The source
does not say this; the vendor's note and the source's practice agree.

| Number | Value | Source |
| --- | --- | --- |
| Time-Saver qualification | `enzymes.json` flags BsaI-HFv2, and a 15-minute functional digest test appears on the BsmBI-v2, BbsI-HF and SrfI sheets | `enzymes.json`; the three specifications |

### Clean-up, ligation and clean-up again

| Number | Value | Source |
| --- | --- | --- |
| Post-digest SPRI | 2x volume ratio, eluted in H₂O | Takacsi-Nagy STAR Methods |
| Ligase | T7 DNA ligase with StickTogether DNA Ligase Buffer (both NEB #M0318L) | same; Key Resources Table |
| Molar ratio | 1:1, insert to backbone | same |
| Backbone mass | 20 ng of digested backbone | same |
| Ligation volume | 200 µL | same |
| Ligation time and temperature | 30 minutes at room temperature | same |
| Post-ligation SPRI | 1x volume ratio, eluted in H₂O | same |
| Which SPRI beads | **not stated for the cloning steps.** The Key Resources Table names AMPure XP (Beckman Coulter) only for the ATAC/NGS size selection and AMPure PB for PacBio | — |
| T7 DNA Ligase, #M0318S/L | 3,000,000 U/mL (1 mg/mL); one unit gives 50% ligation of 100 ng lambda-HindIII in 30 min at 25 °C | M0318 specification PS-M0318S/L v2.0 |
| T7 ligase assay | 16 h at 37 °C in 1x StickTogether buffer gives >95% ligation of a lambda-HindIII digest — **cohesive ends only** | same |
| Ligase units per reaction | **not stated by the paper** | — |
| StickTogether DNA Ligase Buffer | **NEB #B0535S**, supplied as a 2x solution, 2 mL; it ships with M0318 as component B0535AVIAL, 1 × 1 mL with M0318S and 3 × 1 mL with M0318L | B0535 product page; M0318 product page, component table |
| 1x StickTogether buffer | 66 mM Tris-HCl, 10 mM MgCl₂, 1 mM ATP, 1 mM DTT, **7.5% PEG 6000**, pH 7.6 at 25 °C | M0318 product page, reaction conditions |
| T7 ligase incubation | 25 °C | M0318 product page |
| T7 ligase in a buffer without PEG | it still works — T4 DNA Ligase Buffer, NEBuffer r1.1-r3.1 and rCutSmart all serve, with 1 mM ATP added for the NEBuffers, at about **10-fold lower activity** | M0318 product page, note 3 |
| T4 DNA Ligase, #M0202S/L, for the one-pot reactions | 400,000 U/mL (0.4 mg/mL); one unit gives 50% ligation of 6 µg lambda-HindIII in 30 min at 16 °C in 20 µL; stored in 10 mM Tris-HCl, 50 mM KCl, 1 mM DTT, 0.1 mM EDTA, 50% glycerol | M0202 specification PS v1.0 |
| T4 ligase reaction setup, temperature, time, inactivation | **not stated** on the specification sheet; the manual and protocol page are blocked | — |
| Why T7 ligase | it will not join blunt ends, which is what removes escapees — the method page's reason, and departure D2's | `docs/research/synthesis-and-assembly-departures.md` |

### Does T7 DNA ligase refuse blunt ends?

**Yes, and NEB says so in its own words.** The method's escape-removal argument is sourced.

> "unlike T4 and T3 DNA Ligases, blunt end ligation is not efficiently catalyzed by T7 DNA
> Ligase. Addition of high concentrations of PEG 6000 [≥ 20% (w/v)] to the reaction can force
> T7 DNA Ligase to have measurable activity. However, under typical reaction conditions
> blunt-end DNA ligation does not occur in the presence of T7 DNA Ligase, making it a good
> choice for applications in which blunt and cohesive ends of DNA are present but only the
> cohesive ends are to be joined."
>
> — NEB #M0318 product page, read 2026-10-06

That last clause describes this method's reaction exactly. NEB also publishes the measurement
behind it: blunt fragments (ΦX174-HaeIII, #N3026) and sticky-end fragments (λ-HindIII, #N3012)
in one tube, 200 ng of each substrate, 1 µL of each ligase, 30 minutes at 25 °C in each
ligase's own buffer, resolved on a 1% agarose gel. T7 joins the sticky ends and not the blunt.

**One condition is worth carrying into the protocol.** The refusal is a PEG-concentration
effect, and it is forced at **≥20% (w/v) PEG 6000**. The buffer the method uses,
StickTogether, contains **7.5% PEG 6000** — well under that, and it is the buffer NEB ran its
own negative result in, so the margin is measured rather than assumed. The protocol must not
add PEG to this reaction, and nothing in the method does.

| Number | Value | Source |
| --- | --- | --- |
| Blunt-end ligation by T7 ligase | does not occur under typical conditions | M0318 product page |
| What forces measurable blunt activity | PEG 6000 at ≥20% (w/v) | same |
| PEG in the buffer the method uses | 7.5% PEG 6000 | same, reaction conditions |

The SPRI steps above are quoted as ratios. Beckman's own instructions for use give the
procedure those ratios sit inside:

| Number | Value | Source |
| --- | --- | --- |
| Default PCR-purification ratio | 1.8 µL AMPure XP per 1.0 µL of sample, which binds fragments 100 bp and larger | AMPure XP IFU B37419AB pp. 5-6 |
| Bind | pipette-mix, then 5 minutes at room temperature for maximum recovery | same |
| Magnet | 2 minutes to separate | same |
| Wash | 200 µL fresh 70% ethanol, 30 s at room temperature, aspirated; **two washes** | same |
| Dry | optional, and 2 minutes is the guide's conservative figure before an enzymatic reaction; do not over-dry a bead ring holding fragments 10 kb and larger — a cracked ring elutes badly | same |
| Elute | 40 µL of water, 10 mM Tris-acetate pH 8.0, or TE; pipette-mix 10 times, 2 minutes at room temperature | same |
| Below 40 µL | allowed, but needs extra mixing and may not elute the whole product. In 384-well format the guide uses 30 µL, minimum 15 µL | same |
| Bead equilibration time | **not stated** anywhere in the guide | — |
| Ethanol strength | Beckman says 70%; NEB's and ONT's bead steps say 80%. All three are quoted as published | AMPure XP IFU; NEB #E7103 manual v7.0_9/22; NBD114.24 protocol |
| Final magnet | 1 minute, then transfer the eluate | same |
| Ethanol trap | measure 70 mL ethanol and 30 mL water separately and combine; topping 70 mL up to 100 mL gives ~65% | same |
| Catalogue | A63880 (5 mL), A63881 (60 mL), A63882 (450 mL) | AMPure XP IFU B37419AB p. 4 |

Beckman's own ratio is 1.8x. Takacsi-Nagy uses 2x and then 1x. Both are quoted as published;
nothing here reconciles them, because a ratio chosen for size selection is not the IFU's
default and the source gives no reason for either figure.

### Electroporation, recovery and growth

| Number | Value | Source |
| --- | --- | --- |
| DNA into the cells | up to 100 ng of purified ligation product | Takacsi-Nagy STAR Methods |
| Cells | Endura electrocompetent (Lucigen #60242-2) | same; Key Resources Table |
| Instrument | Gene Pulser XCell (Bio-Rad) | same |
| Settings | **not stated by the paper** | — |
| Recovery | 1 hour at **30 °C** with shaking | same |
| Growth before prep | 12-16 hours at **30 °C** | same |
| Plasmid prep | ZymoPURE II (Zymo #D4201 midi, #D4203 maxi), eluted in H₂O | same; Key Resources Table |
| Colonies per round | **not stated** | — |
| Selection antibiotic and concentration for a round | **not stated** | — |
| Selection antibiotic for a round, ours | the destination vector's own marker, read off the record: **50 µg/mL kanamycin** for the KanR destination D11 gives this method, and the final transfer on the working vector's own | Zero Blunt TOPO UG p. 13 |

Endura's own manual supplies the settings the paper omits, at 37 °C:

| Number | Value | Source |
| --- | --- | --- |
| Cuvette | 0.1 cm gap, pre-chilled | Endura manual MA133 p. 4 |
| Cells per transformation | 25 µL | Endura manual MA133 p. 5 |
| Thaw | on wet ice, 10-20 minutes, completely | same |
| Optimal settings | 10 µF, 600 Ω, 1800 V | Endura manual MA133 p. 4 |
| Alternate settings | 25 µF, 200 Ω, 1400-1600 V — 20-50% lower efficiency | same |
| Expected time constant | 3.5-4.5 ms | same |
| DNA ceiling | more than 2 µL of ligation mix may arc; DNA must be in water or a low-ionic buffer such as TE | same |
| Recovery medium | 975 µL of the supplied Recovery Medium, added within 10 seconds of the pulse; SOC lowers efficiency | Endura manual MA133 p. 5 |
| Recovery | 1 hour at 37 °C, 250 rpm | same |
| Plating | up to 100 µL per plate, overnight at 37 °C | same |
| Plate medium | low-salt agar such as LB Lennox; LB Miller varies colony size | Endura manual MA133 p. 4 |
| Efficiency | >1 × 10¹⁰ cfu/µg pUC19, electrocompetent | Endura manual MA133 pp. 3, 5 |
| Catalogue | 60242-0 (2 × 50 µL), 60242-1 (6 × 50 µL), 60242-2 (12 × 50 µL) | Endura manual MA133 p. 3 |
| Recovery Medium supplied | 1 mL per transformation; the manual notes a library protocol recommending ~2 mL | Endura manual MA133 p. 3 |
| **Recovery at 30 °C** | **no vendor document gives one.** Endura's manual is 37 °C throughout; the 30 °C recovery and growth are Takacsi-Nagy's and come with no efficiency figure | — |

The 25 µL aliquot in Endura's manual and the 50 µL vials in its catalogue table disagree in
kind, not in fact: a 50 µL vial is two transformations. Both figures are the manual's.

### What a round is bounded by

**Decided in #259**, and recorded in `docs/research/synthesis-and-assembly.md` 6.9. In: the
masses above. Out: a plated dilution of the recovery and a no-donor control, gated on net
colonies.

| Number | Value | Source |
| --- | --- | --- |
| Colonies per round the source asks for | **none.** It plates nothing at any round | Takacsi-Nagy STAR Methods |
| What the source judges the library by instead | barcode sequencing of the plasmid library for representation, and a long-read amplicon over the gene and its barcodes for linkage | Takacsi-Nagy, Figures 1D and 1E |
| Linkage the source reached | ~95% of reads carried three valid barcodes; nearly 90% of the library correctly linked; >90% of it above 80% fidelity | same |
| Whole-plasmid sequencing | a **sample** of individual clones from the finished libraries, Plasmidsaurus | Takacsi-Nagy METHOD DETAILS p. e4 |
| The colony target for a round | **ours, and a project input.** It follows from the representation the screen downstream needs. Qian's >300 gated a transformation run to yield pickable clones and does not transfer | — |

---

## Stage 5 — Final assembly

The one-pot BsaI reaction, SPRI, and electroporation of a cargo pool.

No published source runs this step as we run it. **The paper runs its own last transfer,
though, and says it ran it on the round's numbers.** Takacsi-Nagy moves the assembled full-gene
libraries into the CAR vector and writes that they were "similarly cloned"; for the control
libraries, that "a single digestion and ligation into the CAR vector, bacterial transformation
and plasmid preparation was performed by the above methods" (METHOD DETAILS, pp. e4-e5). The
above methods are stage 4's: 1 µg of plasmid pool digested in 50 µL, 2× SPRI after the digest,
20 ng of backbone per 200 µL ligation at 1:1, 1× SPRI after it, up to 100 ng electroporated.

So the release that opens this stage and the clean-up that closes it are sourced, and what
departure D12 changes — the assembly itself — is not. Its components:

| Number | Value | Source |
| --- | --- | --- |
| One-pot BsaI reaction, 1 µL or 5 µL | the two recipes quoted in stage 2 | Qian SI, "GGA cloning" |
| Cycling for a two-part assembly | 37 °C 15 min; 60 °C 5 min | Qian SI, "GGA cloning" |
| NEB's kit reaction, 20 µL | #E1601S/L, the table in stage 2 | E1601 manual v5.0_6/26 p. 6 |
| E1601 cycling, 1 insert | 37 °C 5 min for cloning, or **37 °C 1 h for library prep**; then 60 °C 5 min | E1601 manual p. 6 |
| E1601 cycling, 2-10 inserts | (37 °C 1 min, 16 °C 1 min) × 30; then 60 °C 5 min | E1601 manual p. 6 |
| E1601 cycling, 11-20+ inserts | (37 °C 5 min, 16 °C 5 min) × 30; then 60 °C 5 min | E1601 manual p. 6 |
| SPRI | the AMPure XP procedure in stage 4 | AMPure XP IFU B37419AB |
| Mass into the release digest | 1 µg of the finished library in 50 µL with CutSmart, 2.5 µL of each enzyme | METHOD DETAILS p. e4, carried by "the above methods" |
| SPRI ratio after the assembly | **1x**, eluted in H2O, immediately before the electroporation | METHOD DETAILS pp. e4-e5, same carry |
| Electroporation | Endura's settings, or the BL21(DE3) route Qian uses for a library | Endura manual MA133; Qian SI Day 1.2 |

NEB publishes settings for its own electrocompetent strain, #C3020. They are not Endura's, and
the two disagree on every setting — which is the point: an electroporation program belongs to
the cells, not to the method.

| Number | Value | Source |
| --- | --- | --- |
| Cuvette | 1 mm, chilled | C3020 electroporation protocol |
| Cells and DNA | 25 µL of cells, 1 µL of DNA | same |
| Settings | 2.0 kV, 200 Ω, 25 µF, on a BTX ECM 630 or Bio-Rad GenePulser | same |
| Expected time constant | 4.8-5.1 ms | same |
| Recovery | 975 µL of 37 °C NEB 10-beta/Stable Outgrowth Medium added immediately; 1 hour at 37 °C, 250 rpm | same |
| Plates | pre-warmed at 37 °C for 1 hour; incubated overnight at 37 °C | same |
| Against Endura | Endura is 1800 V, 600 Ω, 10 µF with a 3.5-4.5 ms time constant, in a 0.1 cm cuvette | Endura manual MA133 p. 4 |
| Mass of cargo pool into the reaction | **ours** — it is not pipetted at all. D12 purifies nothing between the release and the assembly, so whatever the release freed is what the pot holds | — |
| Mass of working vector into the reaction | **ours** — and no published reaction fits this pot. NEB's two-fragment Golden Gate is 15 µL total, less than the release volume alone; Qian's pooled 0.1 pmol vector against 0.3 pmol library is a column-cleaned amplicon GGA; Takacsi-Nagy's 20 ng per 200 µL is a T7 ligation of purified fragments at room temperature | — |
| Colonies the library needs | the package computes this floor from the completeness asked for; Qian's own rule, CFU × complexity > 300, is the nearest sourced figure | Qian SI Day 2 |
| Expected representation | at least 99.5% of members seen, a 90th/10th percentile skew ratio under 10, judged at over 100 reads a member | Joung et al. 2017, quoted in `docs/research/vector-qc-panel.md` section 3.5 |

The method page's "Neither by-product needs a check" is a design claim, not a number.

---

## The holes, collected

Twenty, by stage. Each is a number the generated protocol would otherwise have to invent.

**The numbering has five gaps, and they are deliberate.** H11, H17, H19, H20 and H26 were holes
on the first pass and were closed on the second, when the `r.jina.ai` route reached the NEB
pages that a plain fetch cannot. The numbers are not reused, so a reference to a hole made
before that pass still resolves. What closed:

| Was | What it asked | Where the answer is now |
| --- | --- | --- |
| H11 | NEB Stable's transformation protocol, and why 30 °C | stage 2, and the C3040 product note that gives the reason |
| H17 | heat inactivation for the five enzymes | stage 4, **The five enzymes** |
| H19 | whether T7 DNA ligase refuses blunt ends | stage 4, **Does T7 DNA ligase refuse blunt ends?** — it does, in NEB's own words |
| H20 | the StickTogether buffer's catalogue number | stage 4: #B0535S, with its 1x composition |
| H26 | electroporation settings for the library strain | stage 5, NEB's #C3020 protocol |

### Stage 1 — part cloning

- **H1.** **No polymerase is named.** The kit insert says "a thermostable proofreading polymerase" and
  names Platinum Pfx only as an example in an ordering table. The method page names none at
  all. *Filled by:* the lab choosing one and a catalogue number, after which Tm, Ta and the
  PCR profile come from `liulab_mbio.primers.polymerase`.
- **H2.** **No insert-to-vector molar ratio for the TOPO reaction.** The guide gives 0.5-4 µL of PCR
  product and nothing else. *Filled by:* measuring the product's concentration and fixing a
  ratio, or accepting the volume range as the specification.
- **H3.** **No per-cycle times for the retailoring PCR.** *Filled by:* the chosen polymerase's own
  manual.
- **H4.** **No per-cycle times for the colony PCR.** The guide gives 20-30 cycles and the two bookends
  only. *Filled by:* the same.
- **H5.** **The guide never states which strains ccdB kills**, only that non-recombinants die in the
  cells supplied. *Filled by:* the ccdB literature, or the strain vendor.

### Stage 2 — working cassette

- **H6.** **No mass, volume or time for the preparative BsmBI digest of the input working vector.**
  *Filled by:* NEB's digest guidelines plus a decision about scale.
- **H7.** **No gel percentage, and none is derivable**, because the RFP stuffer's length is specified
  only as "long enough that the cut backbone separates". *Filled by:* fixing the stuffer's
  length, which the design can then compute from.
- **H8.** **Nothing quantifies the white-colony pick.** No expected white fraction, no false-positive
  rate, no colonies to pick. *Filled by:* a pilot, or by dropping the RFP route for the
  ccdB counter-selection the DMX vector already uses.
- **H9.** **Nothing names the strain, medium or antibiotic for this assembly's transformation.**
  NEB Stable is named; its volumes are the C3040H protocol's, above. What this step
  needs and nobody states is which of them applies here. *Filled by:* a decision.
- **H10.** **The carrier's own BsaI site and ccdB, in a part added without a digest**, is undecided in
  the method itself. *Filled by:* a decision, not a source.

### Stage 3 — cargo and read-out

- **H12.** **No cycling for PCR2 of the three-primer scheme.** Half closed by
  `oligo-pool-pcr-cycles.md`: **PCR1 now has a source.** Twist bands the count by the pool's
  length, 12-14 cycles at 151-350 nt, and the step prints the fewest of its band. **PCR2 still
  has none** — Twist's table amplifies the pool as it arrives, this reaction's template is
  PCR1's product, and Qian's 14-cycle and 10-cycle programs are for a different primer layout.
  *Filled by:* a pilot titrated against the heteroduplex hump, or a real-time run stopped before
  plateau.
- **H13.** **No electroporator settings in Qian**, only "according to the manufacturer's protocol" and
  an expected 4.1 ms time constant. *Filled by:* the Lucigen E. cloni EXPRESS BL21(DE3)
  manual, which was not read for this note.
- **H14.** **No reads-per-well pass mark for Route A that Qian states as one.** Qian's SI gives 150 only
  as the depth above which its pipeline writes a consensus; nobody there calls it a pass mark. The
  method page's "What counts as a pass" reads it as the Route A floor, which is our adaptation.
  Route B has none of this hole: LevSeq's SI checklist sets it (above 20, with 10 tolerable).
  No ONT document gives a per-sample depth target for plate-scale amplicon barcoding. *Filled by:*
  the lab, from a first Route A run.
- **H15.** **No run length and no expected output per flow cell** on any public ONT page read; the
  protocol says only to start MinKNOW at default settings. *Filled by:* an ONT flow cell
  specification sheet, or a first run.
- **H16.** **The SQK-LSK114 protocol itself was not read**, only the native-barcoding amplicon
  protocol that shares its end prep and loading. Where the two differ, this note does not
  know. *Filled by:* the LSK114 protocol page.

### Stage 4 — the iGGA round

- **H18.** **The SPRI beads are not identified for the cloning steps.** AMPure XP is named only for
  the NGS size selection. Beckman's own default is 1.8x; the paper uses 2x then 1x, with no
  reason given for either. *Filled by:* naming a bead product and a ratio.
- **H21.** **No electroporation settings in the paper**, and **no vendor document gives a recovery at
  30 °C** — Endura's manual is 37 °C throughout. The 30 °C recovery and the 12-16 hour 30 °C
  growth are the paper's, with no efficiency attached. *Filled by:* Endura's manual for the
  settings; the 30 °C penalty is unmeasured anywhere.
- **H22.** **No colony count and no selection antibiotic for a round.** Endura's manual gives cfu/µg
  and "plate up to 100 µL" and nothing about library complexity, dilution plating or titre
  plates. *Filled by:* the vector's own marker, and the colony floor the package computes.
  **The marker half is answered** (ADR 0016): a record annotating one is plated on the drug it
  selects, at the concentration cited above, and the hole is then not raised at all. Only a
  record annotating no marker still carries it.
- **H23.** **No NEB document describes a split digest.** Both kit manuals cover the one-pot reaction,
  destination and inserts together. The split digest is Takacsi-Nagy's alone. *Filled by:*
  nothing — it is the method's own and the paper is its only source.

### Stage 5 — final assembly

- **H24.** **Narrowed to the one-pot assembly.** The release that opens this stage and the SPRI that
  closes it are no longer ours: the paper runs its own last transfer "by the above methods"
  (METHOD DETAILS, pp. e4-e5), which gives 1 µg into the digest and 1× beads after the ligation,
  at exactly these two positions. What stays open is the reaction between them, which departure
  D12 rebuilt: the cargo is never pipetted, because nothing is purified after the release, so
  neither the working vector's mass nor the ratio it meets the cargo at has a source. No
  published reaction fits a pot that already holds 50 µL of release. *Filled by:* a pilot.
- **H25.** **Neither kit manual gives a DNA mass for a pooled-library one-pot assembly.** Both give
  0.05 pmol of vector for a defined assembly and a one-hour single-insert program "for
  library preparation", and no pooled-insert mass or colony target. *Filled by:* a pilot.

## What this note does not cover

- **Route B, index PCR.** Long et al. 2025's per-well PCR numbers are in
  `docs/research/route-b-index-pcr.md`. The plate maps and primer plates are inventoried in
  `docs/research/synthesis-and-assembly-materials.md`.
- **Live NEB product pages.** `neb.com` HTML 403s to `curl` and WebFetch, and the agent that
  fetched these documents chose not to route around the block. Everything NEB here is from a
  PDF manual, a product specification sheet or `enzymes.json`. No catalogue number was
  confirmed against a live page, and the specification sheets date from 2013 to 2020.
- **One file under `bench/neb/` is not a source.** `NEB_web-pages_blocked_search-snippets.md`
  holds search-engine summaries of the blocked pages. It is labelled unverified, nothing in
  this note cites it for a value. Both facts it carried — the StickTogether catalogue number
  and NEB Stable's transformation protocol — are now sourced from the pages themselves, so the
  file is kept only as a record of what was unreadable on the first pass.

## Sources

- Zero Blunt TOPO PCR Cloning Kit user guide, Thermo Fisher Scientific, Pub. No. MAN0000062
  Rev B.0 (2014).
  [assets.thermofisher.com/TFS-Assets/LSG/manuals/zeroblunttopo_man.pdf](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/zeroblunttopo_man.pdf)
- One Shot TOP10 Competent Cells user guide, Thermo Fisher Scientific, Pub. No. MAN0000633
  Rev A.0.
  [assets.thermofisher.com/TFS-Assets/LSG/manuals/oneshottop10_man.pdf](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/oneshottop10_man.pdf)
- Qian, Z. et al. (2026) Accelerating protein design by scaling experimental characterization.
  *Nat. Commun.* [doi:10.1038/s41467-026-76740-5](https://doi.org/10.1038/s41467-026-76740-5),
  Supplementary Information, Days 1-4 and the GGA cloning section.
- Takacsi-Nagy, O. et al. (2026) Synthetic transcription factors designed by domain
  recombination enhance CAR T cell antitumor function. *Cell* 189, 1-20.
  [doi:10.1016/j.cell.2026.07.054](https://doi.org/10.1016/j.cell.2026.07.054). CC BY 4.0.
  STAR Methods, "Pooled Cloning — Comprehensive Libraries", and the Key Resources Table.
- Endura Competent Cells manual, Lucigen / Biosearch Technologies, MA133 26Feb2018.
- Agencourt AMPure XP instructions for use, Beckman Coulter, B37419AB (2016), with B37419AA
  (2013) held beside it.
- Oxford Nanopore Technologies, "Ligation sequencing amplicons — native barcoding V14
  (SQK-NBD114.24)", NBA_9168_v114_revU, 15 Jul 2026, with the .96 version NBA_9170_v114_revS,
  and the store pages for SQK-LSK114, SQK-NBD114.24/.96 and the R10.4.1 Flongle flow cell.
- New England Biolabs instruction manuals: NEBridge Golden Gate Assembly Kit (BsaI-HFv2)
  #E1601S/L v5.0_6/26; (BsmBI-v2) #E1602 v3.0_6/26; Monarch DNA Gel Extraction Kit #T1020
  v2.1_4/21; Monarch PCR & DNA Cleanup Kit #T1030 v3.0_4/21.
- New England Biolabs product specification sheets for #R3733, #R0739, #R3539, #R0629, #R0560,
  #M0202, #M0318, #C3040 and #C3020, with `nc3.neb.com/NEBcutter/data/enzymes.json` as the
  cross-check on buffer and temperature.

Downloads are under `reference_docs/synthesis_and_assembly/bench/`, which is git-ignored. Cite
the source, never the path.

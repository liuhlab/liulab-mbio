---
search:
  exclude: true
---

# Primer design, Tm calculation and colony PCR validation

This note answers the questions in #5 (part of #1). It feeds primer design and evaluation (#8), bench
numbers and colony PCR validation (#13), and the protocol skill (#14). Everything was retrieved on
2026-09-12 unless a capture date is given. A value taken from an archived capture or a third-party
copy says so where it is used.

## Sources and how they were read

`www.neb.com` HTML pages, and most of its PDFs, answer scripts and a headless browser with HTTP 403 or
a bot check. NEB material therefore came by these routes:

| Route | What it gave |
| --- | --- |
| `tmcalculator.neb.com`, rendered in headless Chrome (app version 1.17.0) | Help page text; app script `scripts/main-690de58338.js` (Ta rules, salt code, default primer concentration) |
| NEB Tm API `https://tmapi.neb.com/tm` (API version 0.19.0), called with curl | Live Tm and Ta for test primers; product codes from `/docs/productcodes` |
| Wayback Machine raw captures | `tmcalculator.neb.com/data/tmcalculatordata.json` (captured 2023-10-09); 1 kb Plus DNA Ladder page (2026-03-07); 100 bp DNA Ladder page (2025-10-31); Q5 M0491 protocol page (2024-08-04); Agarose Gel Resolution table (2024-09-18); Nucleic Acid Data page (2024-04-20); pUC19 sequencing-primer FAQ (2025-09-10); NEBioCalculator script `scripts/main-6f3020c533.js` (2026-08-28) |
| NEB PDFs at URLs that search engines index | Q5 High-Fidelity PCR Kit manual E0555 (version 5.0, 7/25); Phusion High-Fidelity PCR Kit manual E0553; Gibson Assembly Master Mix manual E2611/E5510 (version 3.0, 1/26); application note "Robust Colony PCR from Multiple E. coli Strains using OneTaq Quick-Load Master Mixes" (Y. Xu, 11/13); application note "Universal Annealing Temperature in PCR" (06/26); product specification N0550 |
| NEB-authored protocols on protocols.io | OneTaq M0480, version 2 (2022-02-21, doi:10.17504/protocols.io.bd24i8gw); Taq M0273, version 1 (2015-01-29, doi:10.17504/protocols.io.ch7t9m) |
| NEB product specification PDFs, `https://www.neb.com/-/media/catalog/specifications/<r>/<n>/<cat>s_l_v<N>.pdf`, with plain curl | 1 kb Plus (PS-N3200S/L v2.0), 100 bp ladder (PS-N3231S/L v1.0), Q5 (PS-M0491S/L v2.0), Quick-Load Purple 1 kb Plus (PS-N0550S/L v2.0). These carry release criteria and storage conditions, not band tables |
| Third-party copy | 100 bp DNA Ladder datasheet, lot 0951209, hosted by genexpress.cl; it agrees with the archived NEB page |

primer3 and primer3-py were read at source: libnano/primer3-py tag `v2.3.1` (`primer3/bindings.py`,
`primer3/thermoanalysis.pyx`, `primer3/argdefaults.py` and the bundled `primer3/src/libprimer3`, which
its `ABOUT.txt` says is "derived from the Primer 2.6.1 source"), and the primer3 manual
`src/primer3_manual.htm` at primer3-org/primer3 tag `v2.6.1`. primer3-py 2.3.1 (bioconda, Python 3.13)
was also run, in a throwaway `pixi exec` environment that was not added to the project. Every
primer3-py number below comes from that run.

Not read: the Owczarzy 2004 and 2008 full texts (paywalled). Their abstracts come from PubMed, and
their coefficients from two independent implementations (primer3 C source, NEB calculator script)
that agree with each other.

## 1. Tm calculation

### 1.1 SantaLucia 1998 nearest-neighbour model

SantaLucia J Jr (1998) "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor
thermodynamics", PNAS 95:1460–1465, doi:10.1073/pnas.95.4.1460. Unified parameters (Table 2, 1 M
NaCl):

| Stack | ΔH° (kcal/mol) | ΔS° (cal/K·mol) | ΔG°37 (kcal/mol) |
| --- | --- | --- | --- |
| AA/TT | −7.9 | −22.2 | −1.00 |
| AT/TA | −7.2 | −20.4 | −0.88 |
| TA/AT | −7.2 | −21.3 | −0.58 |
| CA/GT | −8.5 | −22.7 | −1.45 |
| GT/CA | −8.4 | −22.4 | −1.44 |
| CT/GA | −7.8 | −21.0 | −1.28 |
| GA/CT | −8.2 | −22.2 | −1.30 |
| CG/GC | −10.6 | −27.2 | −2.17 |
| GC/CG | −9.8 | −24.4 | −2.24 |
| GG/CC | −8.0 | −19.9 | −1.84 |
| Initiation, terminal G·C | 0.1 | −2.8 | 0.98 |
| Initiation, terminal A·T | 2.3 | 4.1 | 1.03 |
| Symmetry correction | 0 | −1.4 | 0.43 (computed as −310.15 K × ΔS°) |

- Both primer3 (`oligotm.c`) and the NEB calculator script add the initiation term once for each end,
  choosing the A·T or G·C row by the terminal base.
- Duplex Tm: `Tm = ΔH° / (ΔS° + R ln(CT/x))`, R = 1.987 cal/(K·mol), CT = total strand concentration,
  x = 4 for two non-self-complementary strands at equal concentration, x = 1 for a self-complementary
  strand.
- Salt, for oligomers: `ΔS°[Na⁺] = ΔS°(1 M NaCl) + 0.368 × N × ln[Na⁺]`, N = total phosphates in the
  duplex divided by 2; the matching ΔG°37 coefficient is −0.114. primer3 codes N as length − 1:
  `delta_S + 0.368 * (len - 1) * log(K_mM / 1000.0)`.
- The paper covers Na⁺ only. It has no Mg²⁺ treatment.

### 1.2 Salt corrections for PCR buffers

Concentrations inside a logarithm are in mol/L unless marked mM.

**Divalent-to-monovalent equivalent.** von Ahsen N, Wittwer CT, Schütz E (2001) Clin Chem
47:1956–1961, as given in the primer3 manual under `PRIMER_SALT_DIVALENT` and coded in
`divalent_to_monovalent`:

```text
[Mon]eq (mM) = [Mon] (mM) + 120 × sqrt([Mg2+] − [dNTP])   (mM inside the root)
if [dNTP] ≥ [Mg2+], the divalent term is dropped
```

primer3 applies it in its `schildkraut` and `santalucia` salt modes, to primer Tm, product Tm,
hairpins, dimers and template mispriming.

**Schildkraut and Lifson (1965)**, Biopolymers 3:195–208: `Tm = Tm(1 M) + 16.6 × log10([Mon])`.

**Owczarzy et al. (2004)**, "Effects of sodium ions on DNA duplex oligomers: improved predictions of
melting temperatures", Biochemistry 43:3537–3554, doi:10.1021/bi034621r, PMID 15035624. The abstract
reports 92 duplexes at 69 mM–1.02 M Na⁺, a quadratic dependence on ln[Na⁺], an average error of
1.6 °C on an independent set, and a correction that depends on G·C fraction but not on length or
DNA concentration. Eq. 22, as coded in primer3 and in the NEB script (Tm in K):

```text
1/Tm = 1/Tm(1 M) + (4.29 × fGC − 3.95) × 1e-5 × ln[Na+] + 9.40e-6 × (ln[Na+])^2
```

**Owczarzy et al. (2008)**, "Predicting stability of DNA duplexes in solutions containing magnesium
and monovalent cations", Biochemistry 47:5336–5353, doi:10.1021/bi702363u, PMID 18422348. From the
abstract: "If the concentration ratio of [Mg (2+)] (0.5)/[Mon (+)] is less than 0.22 M (-1/2),
monovalent ions (K (+), Na (+)) are dominant. Effects of magnesium ions dominate and determine duplex
stability at higher ratios. Typical reaction conditions for PCR and DNA sequencing (1.5-5 mM magnesium
and 20-100 mM monovalent cations) fall within this range."

Eq. 16, as coded in primer3 `oligotm.c` and in the NEB script (the same coefficients in both):

```text
R = sqrt([Mg2+]) / [Mon+]
R < 0.22        -> use Eq. 22 above
R ≥ 0.22        -> 1/Tm = 1/Tm(1 M) + a + b·ln[Mg] + fGC·(c + d·ln[Mg])
                                   + (e + f·ln[Mg] + g·(ln[Mg])^2) / (2·(N − 1))
a = 3.92e-5   b = −9.11e-6   c = 6.26e-5   d = 1.42e-5
e = −4.82e-4  f = 5.25e-4    g = 8.31e-5
0.22 ≤ R < 6    -> a = 3.92e-5 × (0.843 − 0.352 · sqrt[Mon] · ln[Mon])
                   d = 1.42e-5 × (1.279 − 4.03e-3 · ln[Mon] − 8.03e-3 · (ln[Mon])^2)
                   g = 8.31e-5 × (0.486 − 0.258 · ln[Mon] + 5.25e-3 · (ln[Mon])^3)
```

primer3 uses free Mg²⁺ = [Mg²⁺] − [dNTP], set to 1e-11 M when dNTP ≥ Mg²⁺, and sets R = 6 when
[Mon] = 0.

**Finding: primer3 drops the length term of Eq. 16.** `oligotm.c` writes the term as
`(1/(2 * (len - 1))) * (e + (f * log(free_divalent)) + g * (pow((log(free_divalent)),2)))`. `len` is
an `int`, so `1/(2 * (len - 1))` is integer division and is 0 for every oligo. The same line is in
upstream primer3 `main` (commit 345cba9, 2026-08-02, `src/oligotm.c` line 518) and in the copy
bundled with primer3-py 2.3.1. The NEB script writes `.5/(this.wseq.length-1)` and keeps the term.
Checked with primer3-py 2.3.1 at 50 mM K⁺, 2 mM Mg²⁺, 0.2 mM dNTP, 500 nM, `owczarzy` (R = 0.85):

| Primer | Length | primer3-py | Eq. 16 with length term | Eq. 16 without |
| --- | --- | --- | --- | --- |
| GTAAAACGACGGCCAGT | 17 | 58.925 | 59.236 | 58.925 |
| CGCCAGGGTTTTCCCAGTCACGACG | 25 | 73.363 | 73.588 | 73.363 |
| ATGAGTAAAGGAGAAGAACTTTTC | 24 | 57.616 | 57.831 | 57.616 |

So primer3's `owczarzy` mode reads 0.2–0.3 °C low whenever R ≥ 0.22. It does not affect
`santalucia` or `schildkraut`, or `owczarzy` with `dv_conc=0`.

### 1.3 primer3-py 2.3.1 API

Signatures, introspected from the installed build:

```python
calc_tm(seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, dmso_conc=0.0,
        dmso_fact=0.6, formamide_conc=0.0, annealing_temp_c=-10.0, max_nn_length=60,
        tm_method='santalucia', salt_corrections_method='santalucia') -> float
calc_hairpin(seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0,
             max_loop=30, output_structure=False) -> ThermoResult
calc_homodimer(seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0,
               max_loop=30, output_structure=False) -> ThermoResult
calc_heterodimer(seq1, seq2, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0,
                 temp_c=37.0, max_loop=30, output_structure=False) -> ThermoResult
calc_end_stability(seq1, seq2, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0,
                   temp_c=37.0, max_loop=30) -> ThermoResult
```

| Keyword | Unit and meaning |
| --- | --- |
| `mv_conc` | mM monovalent cations (primer3 `PRIMER_SALT_MONOVALENT`, "usually KCl") |
| `dv_conc` | mM divalent cations (`PRIMER_SALT_DIVALENT`) |
| `dntp_conc` | mM, sum of all four dNTPs; the manual's example: 0.2 mM each = 0.8 |
| `dna_conc` | nM of each annealing oligo; primer3 uses CT/4 in the Tm equation |
| `temp_c` | °C at which ΔG is reported |
| `max_loop` | largest loop considered in structures |
| `dmso_conc`, `dmso_fact` | % DMSO, and °C subtracted per % (default 0.6) |
| `formamide_conc` | mol/L formamide |
| `annealing_temp_c` | °C; −10 means off. When positive, primer3 computes percent bound, which `calc_tm`'s float return does not expose |
| `max_nn_length` | longer sequences use `long_seq_tm`: `Tm = 81.5 + 16.6·log10([Mon]eq/1000) + 41·fGC − 600/N` |
| `tm_method` | `'santalucia'` (primer3 value 1, "RECOMMENDED") or `'breslauer'` (0) |
| `salt_corrections_method` | `'santalucia'` (1, "RECOMMENDED"), `'owczarzy'` (2), `'schildkraut'` (0) |

Traps:

- `dna_conc` is nM. The `_ThermoAnalysis` class docstring says "(mM)", which is wrong: `argdefaults.py`
  and `bindings.py` say nM, and the default 50 matches primer3's `PRIMER_DNA_CONC` "(nM)".
- The class keyword is `salt_correction_method` (singular). `bindings.calc_tm` spells it
  `salt_corrections_method`.
- `ThermoResult.dg` and `.dh` are in cal/mol, `.ds` is in cal/(K·mol) and `.tm` is in °C. Divide `dg` by
  1000 before comparing it with a kcal/mol threshold.
- `ThermoResult` has `structure_found`, `tm`, `dg`, `dh`, `ds`, `ascii_structure`,
  `ascii_structure_lines`, `check_exc()` and `todict()`. `ascii_structure` is filled only with
  `output_structure=True`.
- Alignment types (`thermoanalysis.pyx`): homodimer and heterodimer use `thal_alignment_any`, hairpin
  uses `thal_alignment_hairpin`, end stability uses `thal_alignment_end1` (3′ end of `seq1` anchored).
  `thal_alignment_end2` is not exposed. primer3 itself scores pair 3′ complementarity with both END1
  and END2 and keeps the larger (`libprimer3.c`, `align_thermod(s1, s2_rev, end1)` and `end2`).

Worked values:

| Call | Result |
| --- | --- |
| `calc_tm('GTAAAACGACGGCCAGT')` (M13 fwd), defaults | 54.70 °C |
| same, `salt_corrections_method='owczarzy'` | 55.20 °C |
| same, `salt_corrections_method='schildkraut'` | 50.45 °C |
| same, `dv_conc=0, dntp_conc=0, dna_conc=250` | 51.67 °C |
| same, `mv_conc=50, dv_conc=2.0, dntp_conc=0.2, dna_conc=500, salt_corrections_method='owczarzy'` | 58.93 °C |
| `calc_tm('CAGGAAACAGCTATGAC')` (M13 rev), defaults | 49.10 °C |

The default value checks out by hand. ΔH° = −133.6 kcal/mol and ΔS° = −360.7 cal/(K·mol): the ten
stacks, plus initiation for a 5′ G·C end and a 3′ A·T end. `[Mon]eq = 50 + 120 × sqrt(1.5 − 0.6) =
163.84 mM`, `ΔS°salt = ΔS° + 0.368 × 16 × ln(0.16384)`, and
`Tm = 1000·ΔH° / (ΔS°salt + 1.987 × ln(50e-9 / 4)) − 273.15 = 54.70 °C`.

3′-anchored dimers, measured at defaults:

| Primer | `calc_homodimer` | `calc_end_stability(p, p)` |
| --- | --- | --- |
| `ACTTACTGACTACGAATTCG` (GAATTC at the 3′ end) | tm 25.64 °C, dg −8.44 kcal/mol | tm 25.64 °C, dg −8.44 kcal/mol: same structure, 3′-anchored |
| `GAATTCGACTTACTGACTAC` (GAATTC at the 5′ end) | tm −19.25 °C, dg −2.53 kcal/mol | tm −102.27 °C, dg −0.54 kcal/mol: no 3′-anchored duplex |

`calc_end_stability(seq1, seq2)` takes `seq2` written 5′→3′ as the strand `seq1` would anneal to:

- self-dimer: `(p, p)`
- cross-dimer: `(fwd, rev)` and `(rev, fwd)`
- a primer on its template: `(p, reverse_complement(p))`, which gives tm 54.61 °C for M13 fwd
  (`calc_tm` gives 54.70 °C)

### 1.4 Reproducing the NEB Tm Calculator

Help page text (app version 1.17.0):

- "Tm values are calculated using thermodynamic data from Santa Lucia [1] and the salt correction of
  Owczarzy [2]. For Phusion® DNA Polymerases, the salt correction of Schildkraut [2] is used."
  References given: SantaLucia (1998) PNAS 95:1460-5; Owczarzy et al (2004) Biochem 43:3537-54;
  Schildkraut et al (1965) Biopolymers 3, 195-208.
- `Tm = (ΔH°i + ΔH°)·1000 / (ΔS°i + ΔS° + R·ln Cp) − 273.15`, where "the primer concentration Cp is
  assumed to be significantly greater (6x) than the target template concentration".
- Owczarzy: `Tm(corrected) = 1 / (1/Tm + [(4.29·fgc − 3.95)·ln(m) + 0.94·(ln(m))²]·10⁻⁵)`, m the
  monovalent cation concentration.
- Phusion: the page shows `Tm + 16.6·ln(m)`. The script computes `16.6*Math.log(saltc)/Math.LN10`, that
  is log10.
- "The calculations in the NEB Tm calculator were modified in 2016 to correct an error that was
  decreasing the effective primer concentration to one quarter of its input value. The correction
  raised Tm values by roughly 2°C for all polymerase reaction buffers except Phusion, which used a
  different algorithm." In the script, the Phusion path still divides the primer concentration by 4
  (`l=5,M/=4`).

Buffer values from `tmcalculatordata.json` (Wayback capture 2023-10-09). The script uses each value as
m in mM; no divalent term is set, so every product takes the Eq. 22 branch.

| Buffer key | m (mM) | Products (API prodcode) |
| --- | --- | --- |
| `q5` | 150 | Q5 M0491 (`q5-0`), Q5 Hot Start M0493 (`q5hs-0`) |
| `q5mm` | 150 | Q5 2X Master Mix M0492 (`q5-1`) |
| `q5u` | 170 | Q5U Hot Start |
| `q5bd` | 150 | Q5 Blood Direct |
| `phusion_hf`, `phusion_gc` | 222 | Phusion M0530 HF (`phusion-0`), Phusion Master Mix M0531 (`phusion-4`) |
| `onetaq_std` | 54 | OneTaq M0480 Standard (`onetaq-0`), OneTaq 2X Master Mix M0482 (`onetaq-2`) |
| `onetaq_gc` | 80 | OneTaq GC buffer (`onetaq-1`); the script also sets DMSO = 5 for this buffer, so Tm drops 3 °C |
| `standard_taq` | 55 | Taq M0273 (`taq-0`), Taq 2X Master Mix M0270 (`taqmaster-0`), Hot Start Taq M0495 (`hstaq-0`) |
| others | — | `thermopol` 40, `longamp` 100, `longamphs` 100, `multiplex` 90, `hemo_klentaq` 70, `crimson_taq` 55, `phusionflex_hf`/`_gc` 222 |

The capture is still current. Fitting m to the live API on 2026-09-12, without looking at the JSON,
returned 150 mM (Q5, 10 primers), 55 mM (Taq, 6), 54 mM (OneTaq, 4) and 222 mM (Phusion, 8), each
with a maximum residual of 0.005 °C or less.

Default primer concentration (script `setCt`): 500 nM for Q5 and Phusion, 400 nM for LongAmp, 200 nM
otherwise. The API requires `conc`, in µM: `GET
https://tmapi.neb.com/tm?seq1=…&seq2=…&prodcode=q5-0&conc=0.5&fmt=long`.

primer3-py settings that reproduce the NEB Tm:

| NEB product and primer concentration | primer3-py call |
| --- | --- |
| Q5, 500 nM | `calc_tm(s, mv_conc=150, dv_conc=0, dntp_conc=0, dna_conc=2000, salt_corrections_method='owczarzy')` |
| Taq, 200 nM | `calc_tm(s, mv_conc=55, dv_conc=0, dntp_conc=0, dna_conc=800, salt_corrections_method='owczarzy')` |
| OneTaq Standard, 200 nM | `calc_tm(s, mv_conc=54, dv_conc=0, dntp_conc=0, dna_conc=800, salt_corrections_method='owczarzy')` |
| Phusion HF, 500 nM | `calc_tm(s, mv_conc=222, dv_conc=0, dntp_conc=0, dna_conc=500, salt_corrections_method='schildkraut')` |

`dna_conc` is 4 × the primer concentration for the Owczarzy products: primer3 always applies CT/4,
while NEB uses Cp. For Phusion NEB divides by 4 itself, so pass the real concentration. With
`dv_conc=0`, the Eq. 16 defect in §1.2 cannot trigger.

Checked, primer3-py against the live API:

| Primer | Product | NEB API Tm | primer3-py |
| --- | --- | --- | --- |
| GTAAAACGACGGCCAGT | Q5, 0.5 µM | 62.27 | 62.27 |
| CAGGAAACAGCTATGAC | Q5, 0.5 µM | 56.31 | 56.31 |
| CGCCAGGGTTTTCCCAGTCACGACG | Q5, 0.5 µM | 76.60 | 76.60 |
| GTAAAACGACGGCCAGT | Taq, 0.2 µM | 53.96 | 53.96 |
| CAGGAAACAGCTATGAC | Taq, 0.2 µM | 47.99 | 47.99 |
| GTAAAACGACGGCCAGT | Phusion, 0.5 µM | 56.57 | 56.57 |
| CAGGAAACAGCTATGAC | Phusion, 0.5 µM | 50.89 | 50.89 |
| CGCCAGGGTTTTCCCAGTCACGACG | Phusion, 0.5 µM | 71.12 | 71.12 |

## 2. Annealing temperature

### 2.1 Rules in the NEB calculator

From `getAnnealTemp` in the script. s is the lower Tm of the pair and u the length of the shorter
primer. Ta uses the unrounded Tm and is shown to 0.1 °C (56.31 + 1 → 57.3).

| Product group | Ta | Cap |
| --- | --- | --- |
| Q5, Q5 Hot Start, Q5 Blood Direct; Q5 master mixes | s + 1 (when u > 7) | 72 °C |
| Q5U Hot Start | s + 2 | 72 °C |
| Phusion, Phusion Hot Start Flex; Phusion master mixes | 0.93 × s + 7.5 | 72 °C |
| Taq, Hot Start Taq, OneTaq, OneTaq Hot Start, Hemo KlenTaq, EpiMark Hot Start | s − 5 | 68 °C |
| LongAmp Taq, LongAmp Hot Start Taq | s − 5 | 65 °C |
| Vent, Deep Vent | s − 2 when u > 20 | 72 °C |
| Other master mixes | s − 5 | 68 °C |

Warnings (`validateTm`, `validateInput`):

- Tm difference greater than 5 °C: "Tm difference is greater than the recommended limit of 5 °C."
- Ta below 55 °C for Q5: "The minimum recommended annealing temperature for Q5 is 55 °C." Below 60 °C
  for Q5U. Below 45 °C otherwise: "Annealing temperature is lower than the recommended minimum of
  45 °C."
- A note when Ta reaches 72 °C for Q5, Phusion, Vent and Deep Vent ("should typically not exceed
  72°C"), or 65 °C for LongAmp.
- Primer concentration below 500 nM for Q5, Q5U or Phusion.
- Each primer must be longer than 7 nt. At most 3 ambiguous bases are allowed; Tm is then the minimum
  over the expanded sequences.

Live API results (2026-09-12):

| Pair (lengths) | Product, primer conc | Tm1 / Tm2 (°C) | Ta (°C) |
| --- | --- | --- | --- |
| GTAAAACGACGGCCAGT / CAGGAAACAGCTATGAC (17/17) | `q5-0`, 0.5 µM | 62.27 / 56.31 | 57.3 |
| same | `phusion-0`, 0.5 µM | 56.57 / 50.89 | 54.8 |
| same | `onetaq-0`, 0.2 µM | 53.82 / 47.85 | 42.8, below-45 warning |
| same | `taq-0`, 0.2 µM | 53.96 / 47.99 | 43.0, below-45 warning |
| CGCCAGGGTTTTCCCAGTCACGACG / AGCGGATAACAATTTCACACAGGAAAC (25/27) | `q5-0`, 0.5 µM | 76.60 / 67.31 | 68.3 |
| same | `phusion-0`, 0.5 µM | 71.12 / 63.86 | 66.9 |
| same | `taq-0`, 0.2 µM | 68.68 / 58.72 | 53.7 |
| CGCCAGGGTTTTCCCAGTCACGACGTTG / GCGGATAACAATTTCACACAGGAAACAGCTATGAC (28/35) | `q5-0`, 0.5 µM | 77.45 / 71.73 | 72 (capped) |
| same | `phusion-0`, 0.5 µM | 72.52 / 68.67 | 71.4 |
| CCCAGTCACGACGTTGTAAAACG / AGCGGATAACAATTTCACACAGG (23/23) | `onetaq-2`, 0.2 µM | 59.94 / 56.49 | 51.5 |
| same | `taq-0`, 0.2 µM | 60.09 / 56.64 | 51.6 |
| same | `q5-0`, 0.5 µM | 68.30 / 65.16 | 66.2 |

Unexplained: `taqmaster-0` (Taq 2X Master Mix) returned Ta 56.6 for the 23/23 pair, which is s with no
−5. `taq-0` and `onetaq-2` returned s − 5, as the script's rules predict.

### 2.2 What the manuals say

- **Q5** (M0491 protocol page, capture 2024-08-04; same wording in the E0555 manual): "Typically, use a
  10–30 second annealing step at 3°C above the Tm of the lower Tm primer." The Tm Calculator help page
  explains the gap with the calculator's own +1: "Tm values obtained from other calculators generally
  underestimate the Tm for use with Q5. If you obtain Tm values from another calculator, we suggest
  raising your annealing temperature to at least Tm + 3. However, for best results when using Q5, we
  still recommend using the annealing temperature provided by the NEB Tm Calculator." The Q5
  application note (06/26) finds a single Ta of 62 °C (M0491, M0493, M2499) or 60 °C (Q5 master mixes)
  works for many pairs, while the per-pair optimum still gives the best yield and specificity.
- **Phusion** (E0553 manual): "Typically, primers greater than 20 nucleotides in length anneal for
  10–30 seconds at 3°C above the Tm of the lower Tm primer. If the primer length is less than 20
  nucleotides, an annealing temperature equivalent to the Tm of the lower primer should be used." The
  calculator uses 0.93 × s + 7.5 with its Schildkraut Tm instead.
- **OneTaq** (protocols.io version 2) and **Taq** (version 1), same wording: "Annealing temperature is
  based on the Tm of the primer pair and is typically 45–68°C. Annealing temperatures can be optimized
  by doing a temperature gradient PCR starting 5°C below the calculated Tm."
- **IDT** ("A practical guide for PCR and qPCR primer design", updated 2024-04-08): Ta "should be no
  more than 5°C below the Tm of your primers".
- **primer3 manual** ("General thoughts on primer binding"): the annealing temperature is "usually
  chosen 6-10°C below the melting temperature of the primers", relative to primer3's own lower Tm
  (50 nM, CT/4).

### 2.3 5′ tails

- The calculator has no notion of a tail. It computes the Tm of whatever sequence is entered.
- NEB's Gibson Assembly manual (E2611/E5510) describes a primer as a 5′ overlap sequence plus a 3′
  gene-specific sequence, and says: "The Tm of the 3´ gene-specific sequence of the primer can be
  calculated using the Tm calculator found on the NEB website." Ta therefore comes from the
  template-matching 3′ part.
- In the first cycle only the 3′ part matches the template. Strands copied from a tailed primer carry
  the tail's complement, so later cycles can anneal the whole primer, whose Tm is higher. Measured:
  GFP start `ATGAGTAAAGGAGAAGAACTTTTC` has Tm 54.81 °C (primer3 defaults) and 59.59 °C (NEB Q5); with
  the 11-nt tail `TTGAAGACAAA` added, 63.05 °C and 66.24 °C.
- Not found in any source read: a program that runs the first cycles at the 3′-part Ta and later
  cycles at the full-primer Ta. It is common practice but unverified here.
- Hairpin and dimer checks should run on the full oligo, because the tail is physically present.
  This is inference from what the checks measure, not a quoted rule.

### 2.4 Two-step PCR

- **Q5**: "When primers with annealing temperatures ≥ 72°C are used, a 2-step thermocycling protocol
  (combining annealing and extension into one step) is possible." (M0491 protocol page, capture
  2024-08-04). The Q5 application note (06/26) agrees: "a combined annealing and extension step at
  72°C is recommended when primers with annealing temperatures ≥ 72°C are used". The E0555 manual
  prints "≤ 72°C" in the same sentence, taken here as a typo. E0555's two-step program: 98 °C 30 s;
  25–35 cycles of 98 °C 5–10 s and 72 °C 15–30 s/kb; 72 °C 2 min; hold 4–10 °C. The application
  note's two-step table uses 20–30 s/kb.
- **Phusion** (E0553): "When primers with annealing temperatures ≥ 72°C are used, a 2-step
  thermocycling protocol is recommended." The PDF encodes the sign as U+F0B3, which is ≥ in the Symbol
  font. Program: 98 °C 30 s; 25–35 cycles of 98 °C 5–10 s and 72 °C 15–30 s/kb; 72 °C 5–10 min; hold
  4 °C.
- **OneTaq**: "When primers with annealing temperatures of 68°C or above are used, a 2-step
  thermocycling protocol (combining annealing and extension into one step) is possible."
- **Taq**: "When primers with annealing temperatures above 65°C are used, a 2-step thermocycling
  protocol is possible", at 65–68 °C for 1 min/kb, final extension 65–68 °C for 5 min.
- Q5 application note (06/26), general: "two-step thermocycling protocols are commonly recommended when
  PCR primers have high Tm values (≥ 68°C)."

## 3. Cycling

### 3.1 Standard programs and extension rates

| | Q5 M0491 | Phusion (E0553 kit) | OneTaq M0480 | Taq M0273 |
| --- | --- | --- | --- | --- |
| Source | protocol page, capture 2024-08-04 | E0553 manual | protocols.io v2 | protocols.io v1 |
| Mg²⁺ in 1X buffer | 2.0 mM (5X Q5 Reaction Buffer) | 1.5 mM (HF or GC buffer) | 1.8 mM Standard; 2.0 mM GC | 1.5 mM Standard Taq buffer |
| dNTPs | 200 µM each | 200 µM each | 200 µM each | 200 µM each |
| Each primer | 0.5 µM | 0.5 µM (0.2–1) | 0.2 µM (0.05–1) | 0.2 µM (0.05–1) |
| Enzyme | 0.02 U/µl (1.0 U/50 µl) | 1.0 U/50 µl (0.5 U per 20 µl in the table) | 1.25 U/50 µl | 1.25 U/50 µl |
| Template, 50 µl | < 1,000 ng; genomic 1 ng–1 µg; plasmid or viral 1 pg–10 ng | < 250 ng; genomic 50–250 ng; plasmid or viral 1 pg–10 ng | < 1,000 ng; genomic 1 ng–1 µg; plasmid or viral 1 pg–10 ng | < 1,000 ng; genomic 1 ng–1 µg; plasmid or viral 1 pg–1 ng |
| Initial denaturation | 98 °C 30 s (up to 3 min) | 98 °C 30 s (up to 3 min) | 94 °C 30 s | 95 °C 30 s |
| Cycles | 25–35 | 30 in the table; 25–35 in the guidelines | 30 in the table; 25–35 in the guidelines | 30 |
| Denaturation | 98 °C 5–10 s | 98 °C 5–10 s | 94 °C 15–30 s | 95 °C 15–30 s |
| Annealing | 50–72 °C, 10–30 s | 45–72 °C, 10–30 s | 45–68 °C, 15–60 s | 45–68 °C, 15–60 s |
| Extension | 72 °C | 72 °C | 68 °C | 68 °C |
| Final extension | 72 °C 2 min | 72 °C 5–10 min | 68 °C 5 min | 68 °C 5 min |
| Hold | 4–10 °C | 4 °C | 4–10 °C | 4–10 °C |

Extension rates:

| Polymerase | Rate | Source |
| --- | --- | --- |
| Q5 | "20–30 seconds per kb for complex, genomic samples, but can be reduced to 10 seconds per kb for simple templates (plasmid, E. coli, etc.) or complex templates < 1 kb"; 40 s/kb for cDNA or long complex templates; 40–50 s/kb for products > 6 kb | M0491 protocol page (capture 2024-08-04), E0555 |
| Phusion | "Generally, an extension time of 15 seconds per kb can be used"; 30 s/kb for complex amplicons such as genomic DNA; 40 s/kb for cDNA | E0553 |
| OneTaq | "Extension times are generally 1 minute per kb" | protocols.io v2 |
| Taq | "Extension times are generally 1 minute per kb" | protocols.io v1 |

The Q5 2X Master Mix (E0555) contains 2 mM MgCl₂ and 200 µM each dNTP at 1X, with primers at 0.5 µM.

### 3.2 Colony PCR

NEB application note, "Robust Colony PCR from Multiple E. coli Strains using OneTaq Quick-Load Master
Mixes" (Y. Xu; dated 11/13, © 2018):

- Reaction (50 µl): 25 µl OneTaq Quick-Load 2X Master Mix with Standard Buffer (M0486) or the Hot Start
  version (M0488), each primer at 200 nM, water to 50 µl. Set up the non-hot-start mix on ice; the Hot
  Start mix can be set up at room temperature.
- Template: "Well-isolated bacterial colonies, ideally 1-2 mm in diameter". "Use a sterile toothpick to
  pick up individual colonies and dip into each reaction tube. As soon as the solution looks cloudy,
  remove the toothpick." Then keep each colony by dipping the toothpick into 3 ml medium with
  antibiotic, or by streaking it onto another plate.
- Program: initial denaturation 2 minutes, then 30 cycles of 94 °C 15–30 s, 45–68 °C 15–60 s and
  68 °C 1 min/kb, then 68 °C 5–10 min and a 10 °C hold. The initial-denaturation temperature cell is
  garbled in the PDF (it reads "PCR Primer"), so take the temperature from the protocols below.
- Load 4–6 µl directly on an agarose gel. The note reports success in 18 E. coli strains and with
  amplicons up to 10 kb.
- For inserts over 65% GC, the note points to the GC-buffer master mixes.

Initial denaturation for colony PCR in the protocols:

- OneTaq: "With colony PCR, an initial 2–5 minute incubation at 94°C is recommended to lyse cells."
- Taq: "With colony PCR, an initial 5 minute denaturation at 95°C is recommended."

Not found in any source read: an NEB recommendation for or against Q5 in colony PCR.

## 4. Primer checks and thresholds

| Check | primer3 2.6.1 default (manual) | NEB | IDT | Other |
| --- | --- | --- | --- | --- |
| Length | `PRIMER_MIN_SIZE` 18, `PRIMER_OPT_SIZE` 20, `PRIMER_MAX_SIZE` 27; max "cannot be larger than 35", the limit of its Tm formula | "generally 20–40 nucleotides in length", tails included (Q5, Phusion, Taq manuals) | 18–30 bases | Sanger primers 18–24 (Genewiz) |
| GC content | `PRIMER_MIN_GC` 20, `PRIMER_MAX_GC` 80, `PRIMER_OPT_GC_PERCENT` 50 | "ideally have a GC content of 40–60%" | 35–65%, ideal 50% (guide); 35–80% (FAQ) | Sanger 45–55% (Genewiz) |
| GC clamp | `PRIMER_GC_CLAMP` 0 (off); `PRIMER_MAX_END_GC` 5 = most G or C allowed in the last five 3′ bases | "Avoid GC-rich 3´ ends" (Q5 E0555 troubleshooting) | — | Sanger: "Have a G or C at 3′ end" (Genewiz) |
| Runs | `PRIMER_MAX_POLY_X` 5 = longest mononucleotide run; N counts as the worst case | — | "should not contain regions of 4 or more consecutive G residues" | No source read gives a dinucleotide-repeat limit |
| Tm | `PRIMER_MIN_TM` 57, `PRIMER_OPT_TM` 60, `PRIMER_MAX_TM` 63 | Use the NEB calculator | 60–64 °C, ideal 62 °C | — |
| Pair Tm difference | `PRIMER_PAIR_MAX_DIFF_TM` 100 (effectively off) | Calculator warns above 5 °C | "should not differ by more than 2°C" (guide); "less than 5°C" (FAQ) | Q5 application note: "a maximum difference of 5°C" |
| 3′ end stability | `PRIMER_MAX_END_STABILITY` 100 (off). ΔG in kcal/mol to disrupt the last five 3′ bases; with SantaLucia 1998 parameters, 6.86 for GCGCG (most stable) to 0.86 for TATAT (most labile). The manual's Primer3Plus settings list 9.0 | — | — | — |
| Hairpin | `PRIMER_MAX_HAIRPIN_TH` 47.0 °C, Tm of the structure; "10 degrees lower than the default value of PRIMER_MIN_TM" | — | ΔG weaker (more positive) than −9.0 kcal/mol | — |
| Self-dimer | `PRIMER_MAX_SELF_ANY_TH` 47.00 °C; `PRIMER_MAX_SELF_END_TH` 47.00 °C (3′-anchored). Alignment-score versions: `PRIMER_MAX_SELF_ANY` 8.00, `PRIMER_MAX_SELF_END` 3.00 | "Verify that primers are non-complementary, both internally and to each other" | ΔG weaker than −9.0 kcal/mol | — |
| Hetero-dimer | `PRIMER_PAIR_MAX_COMPL_ANY_TH` 47.00 °C; `PRIMER_PAIR_MAX_COMPL_END_TH` 47.00 °C; score version `PRIMER_PAIR_MAX_COMPL_END` 3.00 | as above | ΔG weaker than −9.0 kcal/mol; "Avoid 3' complementarity between the two primers" | — |
| Off-target on template | `PRIMER_MAX_TEMPLATE_MISPRIMING_TH` −1 (off); the manual: "47.0 would be a reasonable choice if PRIMER_MIN_TM is 57.0". Needs `PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT` = 1 (default 0). Pair version `PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH` sums both primers' ectopic Tm | "Verify that primers have no additional complementary regions within the template DNA" | Run BLAST for uniqueness | — |

Sources: primer3 manual (tag v2.6.1); NEB Q5 E0555, Phusion E0553 and Taq M0273 manuals; IDT FAQ "How
can I check my PCR primers using the OligoAnalyzer program…" and IDT guide (updated 2024-04-08);
Genewiz Technical Notes page; NEB Q5 application note (06/26).

Notes:

- primer3 explains why the 3′ end matters: `PRIMER_MAX_SELF_END` "tries to bind the 3'-END to a
  identical primer … This is critical for primer quality because it allows primers use itself as a
  target and amplify a short piece (forming a primer-dimer)."
- The primer3 `_TH` limits are structure Tms at primer3's own conditions (50 mM K⁺, 1.5 mM Mg²⁺,
  0.6 mM dNTP, 50 nM). They are set 10 °C below `PRIMER_MIN_TM` 57, so they only mean something next
  to a Tm computed the same way.
- IDT does not state the temperature or salt at which its −9.0 kcal/mol applies. primer3-py reports ΔG
  at `temp_c`, 37 °C by default. The threshold is therefore approximate.

## 5. Validation

### 5.1 Expected colony PCR bands

The general method: simulate PCR on each candidate plasmid and list every product. Primer sites are
found on both strands, the origin is joined for circular sequences, and the size is the distance from
the forward primer's 5′ end to the reverse primer's 5′ end. The candidates are:

1. the correct product,
2. the parental vector (uncut plasmid, or PCR template carried past DpnI),
3. the product with the insert reversed,
4. the linearized vector closed without insert, when that is possible.

Symbols: `E` = empty-vector amplicon between the two vector primers; `D` = vector bases between the
two junctions that the design removes (stuffer, dropout, or bases deleted by outward PCR); `I` = insert
length including any bases the design adds (scar, spacer, start or stop codon).

| Primer layout | Correct | Parental vector | Insert reversed | Notes |
| --- | --- | --- | --- | --- |
| Both primers in the vector, flanking the site (e.g. M13 pair) | `E − D + I` | `E` | `E − D + I`, the same as correct | Cannot tell orientation. Gives a band for every clone, so a failed PCR is visible |
| Vector primer + insert primer across one junction | vector primer to insert primer | none | none; the insert primer now points the wrong way | Tells orientation. No band from a failed PCR either, so include a positive control |
| Three primers: both vector primers + one insert primer | two bands: `E − D + I` and the junction band | one band `E` | `E − D + I` plus a different-size band from the insert primer pairing with the other vector primer | Tells orientation. Read sizes off the simulation rather than deriving them by hand |

Golden Gate points for #13:

- Non-palindromic, distinct overhangs make the reversed product unlikely; the simulation still lists it.
- The parental vector's `E` band is the main background when the vector is PCR-linearized. This is
  why the template is digested with DpnI.
- On pUC19, a lacZα insert also gives white rather than blue colonies (issue #1).

Gel sizing: keep every expected band at 100 bp or more and inside the ladder's range (§6.2). NEB's
agarose table has no row below 0.2 kb, and NEB's 100 bp ladder page gives an optimum of 2% agarose.

### 5.2 Smoke-test numbers (pUC19 fixture)

Primer sites located in `tests/data/pUC19.dna` (1-based, fixture features: M13 fwd 379–395, M13 rev
465–481, MCS 396–452):

| Primer | Sequence 5′→3′ | Named by | pUC19 site |
| --- | --- | --- | --- |
| M13 fwd | GTAAAACGACGGCCAGT | SnapGene feature in the fixture | 379–395, top strand |
| M13 Forward (−21) | TGTAAAACGACGGCCAGT | Addgene (content reviewed 2025-10-22) | 378–395, top |
| M13F | GTAAAACGACGGCCAG | Genewiz free universal primers | 379–394, top |
| M13 Forward (−40) / M13-40FOR | GTTTTCCCAGTCACGAC | Addgene / Genewiz | 359–375, top |
| M13/pUC Forward | CCCAGTCACGACGTTGTAAAACG | Addgene | 364–386, top |
| M13 rev / M13 Reverse / M13R | CAGGAAACAGCTATGAC | SnapGene feature / Addgene / Genewiz | 465–481, bottom |
| M13/pUC Reverse | AGCGGATAACAATTTCACACAGG | Addgene | 478–500, bottom |
| M13-48REV | CGGATAACAATTTCACACAG | Genewiz | 479–498, bottom |

The fixture also contains `CGCCAGGGTTTTCCCAGTCACGAC` (352–375, top) and `AACAGCTATGACCATG` (461–476,
bottom). A search-engine summary named these NEB's "(−47)" and "(−24)" sequencing primers, but no NEB
page could be read to confirm it. NEB's FAQ "What primers should I use to sequence an insert (pUC19,
pNEB193, LITMUS)?" (capture 2025-09-10) only says they are "what most labs refer to as 'universal M13
primers'" and points to the catalog.

Empty-vector amplicons on pUC19:

| Forward | Reverse | Empty amplicon |
| --- | --- | --- |
| M13 fwd | M13 rev | 103 bp (matches #1) |
| M13 Forward (−21) | M13 Reverse | 104 bp |
| M13 fwd | M13/pUC Reverse | 122 bp |
| M13 Forward (−40) | M13 Reverse | 123 bp |
| M13/pUC Forward | M13/pUC Reverse | 137 bp |
| M13 Forward (−40) | M13/pUC Reverse | 142 bp |

Worked example, for illustration only; the real `D` and `I` come from the #8/#9 design. If the 717-bp
GFP CDS replaced exactly MCS 396–452 (`D` = 57) with no added bases, the M13 amplicon would be
`103 − 57 + 717 = 763 bp` against 103 bp for the parental vector.

The M13 17-mers give Ta 42.8 °C in OneTaq (Tm 53.82/47.85), below NEB's 45 °C minimum. The 23-mer
M13/pUC Forward and Reverse give Ta 51.5 °C in the OneTaq 2X Master Mix, and a 137-bp empty band.

### 5.3 Sanger read geometry and primer placement

- Genewiz FAQ: "We recommend designing your sequencing primers in a region that is 100 bases upstream
  of your sequence of interest. If you do not have the luxury of having this buffer, the closest you
  want the primer to be to your area of interest is 50-60 bases."
- Genewiz Technical Notes: "Approximately 18-24 bases in length", "Melting temperature (Tm) between 50
  – 60 degrees", "GC content should be about 45 – 55%", "Have a G or C at 3′ end". The FAQ says Tm
  56–60.
- Read length: no page read here states a typical high-quality read length. The Genewiz FAQ mentions
  "a read length of between 500-1000bp" only for a low-concentration plasmid service. Unverified; ask
  the provider.
- Geometry on the fixture: the M13 fwd 3′ end is at 395, the base before MCS 396. The M13 rev 3′ end
  is 12 bases past MCS end 452, and M13-48REV's is 26 bases past. All are closer than the 50–60-base
  minimum, so each read starts too close to its own junction. Each junction is then covered only by
  the read from the far side, after it has crossed the whole insert. A Sanger primer needs its 3′ end
  at least 50–60 (ideally 100) bases outside each junction, and a read must reach the far junction:
  `distance to near junction + insert length + margin`.

## 6. Gels and quantities

### 6.1 Agarose percentage by size

NEB "Agarose Gel Resolution" (capture 2024-09-18), optimum resolution for linear DNA:

| % gel | Range (kb) |
| --- | --- |
| 0.5 | 1.0–30 |
| 0.7 | 0.8–12 |
| 1.0 | 0.5–10 |
| 1.2 | 0.4–7 |
| 1.5 | 0.2–3 |

Also:

- Azenta/Genewiz "Sanger Quick Tips", Table 1: under 250 bp "Not Recommended"; 250–1,500 bp 2.00%;
  300–4,000 bp 1.50%; 400–7,000 bp 1.20%; 500–9,000 bp 1.00%; 800–11,000 bp 0.75%; 1,000–30,000 bp
  0.50%. The same guide gives about 10 ng per band as the minimum and about 150 ng per band as the
  most.
- NEB ladder pages: 1 kb Plus "Recommended gel percentage range: 1-1.4%", "Optimum separation on 1.2%";
  100 bp ladder "Recommended gel percentage range: 1.2-3%", "Optimum separation on 2%".
- Lee PY et al. (2012) J Vis Exp 62:e3923, doi:10.3791/3923: "most gels ranging between 0.5%-2%".

### 6.2 NEB ladders

**1 kb Plus DNA Ladder, N3200** (NEB page, capture 2026-03-07). Formerly the 2-Log DNA Ladder, "no
changes to product formulation or specifications". 19 bands, 100 bp to 10 kb; "The 0.5, 1.0 and 3.0
kb bands have increased intensity to serve as reference points." Mass at a 1.0 µg load, image on 1.0%
TBE:

| Band | bp | ng |
| --- | --- | --- |
| 1 | 10,002 | 40 |
| 2 | 8,001 | 40 |
| 3 | 6,001 | 48 |
| 4 | 5,001 | 40 |
| 5 | 4,001 | 32 |
| 6 | 3,001 | 120 |
| 7 | 2,017 | 40 |
| 8 | 1,517 | 57 |
| 9 | 1,200 | 45 |
| 10 | 1,000 | 122 |
| 11 | 900 | 34 |
| 12 | 800 | 31 |
| 13 | 700 | 27 |
| 14 | 600 | 23 |
| 15 | 500/517 | 124 |
| 16 | 400 | 49 |
| 17 | 300 | 37 |
| 18 | 200 | 32 |
| 19 | 100 | 61 |

**100 bp DNA Ladder, N3231** (NEB page, capture 2025-10-31). The 2014-lot datasheet copy (lot
0951209) gives the same values. 12 bands, 100–1,517 bp; "The 500 and 1,000 base pair bands have
increased intensity to serve as reference points." Mass at a 0.5 µg load, image on 1.3% TAE:

| Band | bp | ng |
| --- | --- | --- |
| 1 | 1,517 | 45 |
| 2 | 1,200 | 35 |
| 3 | 1,000 | 95 |
| 4 | 900 | 27 |
| 5 | 800 | 24 |
| 6 | 700 | 21 |
| 7 | 600 | 18 |
| 8 | 500/517 | 97 |
| 9 | 400 | 38 |
| 10 | 300 | 29 |
| 11 | 200 | 25 |
| 12 | 100 | 48 |

NEB counts 500/517 as one band in both ladders. The Q5 kit (E0555) ships the Quick-Load Purple 1 kb
Plus DNA Ladder (N0550): "19 bands … from 100 bp to 10 kb", with the same reference bands.

The band tables above come from the archived product pages. The live specification PDFs for both
ladders carry no band table; they give release criteria, and both name a 1.2% agarose gel with
0.5 µg/ml ethidium bromide as the QC condition for the banding pattern.

### 6.3 ng ↔ pmol

- **NEB Nucleic Acid Data** (capture 2024-04-20):
  - "Average weight of a DNA basepair (sodium salt) = 650 daltons"
  - "MW of a double-stranded DNA molecule = (# of base pairs) X (650 daltons/base pair)"
  - "Moles of ends of a double-stranded DNA molecule = 2 X (grams of DNA) / (MW in daltons)"
  - "1 µg of 1000 bp DNA = 1.52 pmol"; "1 µg of pUC18/19 DNA (2686 bp) = 0.57 pmol"; "1 pmol of 1000 bp
    DNA = 0.66 µg"
  - 1.0 A260 unit: dsDNA 50 µg/ml, ssDNA 33 µg/ml, ssRNA 40 µg/ml
- **NEB assembly manuals** (Gibson E2611/E5510, version 3.0 1/26): "pmols = (weight in ng) x 1,000 /
  (base pairs x 650 daltons)"; "50 ng of 5000 bp dsDNA is about 0.015 pmols."
- **NEBioCalculator** (script `main-6f3020c533.js`, capture 2026-08-28):
  - dsDNA MW (g/mol) = 36.04 + bp × 615.94
  - dsDNA MW from sequence = 36.04 + (A+T) × (312.23 + 303.21) + (G+C) × (328.23 + 288.2)
  - ssDNA MW = 18.02 + nt × 307.97; from sequence 18.02 + 312.23·A + 303.21·T + 328.23·G + 288.2·C
  - ssRNA MW = 18.02 + nt × 320.47
  - moles of ends (dsDNA) = 2 × moles

The two NEB constants differ by about 5.5%. 1 µg of pUC19 is 0.573 pmol by the 650-Da rule and
0.604 pmol by NEBioCalculator (`1e-6 / (36.04 + 2686 × 615.94)` mol). A search-engine summary gave
617.96 g/mol per bp for NEBioCalculator; the script says 615.94.

## Implications for this package

### For #8, primer design and evaluation

1. **Tm and Ta, NEB-compatible.** Use SantaLucia 1998 through primer3-py, with a per-polymerase
   profile. The output should say, for example: "Tm: SantaLucia 1998 nearest-neighbour, Owczarzy 2004
   salt correction at the NEB buffer's monovalent equivalent (NEB Tm Calculator method); Ta: NEB rule
   for Q5".

   | Profile | `salt_corrections_method` | `mv_conc` | `dv_conc`, `dntp_conc` | `dna_conc` passed | Ta | Cap | Warn if Ta below |
   | --- | --- | --- | --- | --- | --- | --- | --- |
   | Q5 (default), 500 nM | `owczarzy` | 150 | 0, 0 | 2000 | lower Tm + 1 | 72 | 55 |
   | Phusion HF, 500 nM | `schildkraut` | 222 | 0, 0 | 500 | 0.93 × lower Tm + 7.5 | 72 | 45 |
   | OneTaq Standard, 200 nM | `owczarzy` | 54 | 0, 0 | 800 | lower Tm − 5 | 68 | 45 |
   | Taq Standard, 200 nM | `owczarzy` | 55 | 0, 0 | 800 | lower Tm − 5 | 68 | 45 |

   Round Ta to 0.1 °C and compute it from unrounded Tm. Warn when the pair's annealing-part Tm differ
   by more than 5 °C (NEB, IDT). Pin these in tests: M13 fwd/rev Q5 62.27/56.31, Ta 57.3; Phusion
   56.57/50.89, Ta 54.8; Taq 53.96/47.99, Ta 43.0; the 28/35-nt pair under Q5, Ta 72 (capped).
2. **Do not use primer3's `owczarzy` mode with `dv_conc` > 0** until the Eq. 16 defect (§1.2) is fixed
   upstream; the profiles above avoid it.
3. **Tails.** Take Tm and Ta from the 3′ annealing part only. Report the full-primer Tm for
   information. Run hairpin and dimer checks on the full oligo.
4. **Two-step.** When the computed Ta hits the cap (72 °C for Q5 and Phusion), emit a two-step program
   at 72 °C. For OneTaq the trigger is Ta ≥ 68 °C, for Taq Ta > 65 °C.
5. **Checks with sourced defaults.** Thresholds live in one table. Rows marked *proposal* are this
   note's choices, not a quoted rule.

   | Check | Pass | Warn | Fail | Basis |
   | --- | --- | --- | --- | --- |
   | Annealing length | 18–30 nt | 15–17 or 31–35 nt | < 15 or > 35 nt | primer3 18/27 and 35 cap; IDT 18–30; *proposal* for the warn band |
   | Annealing length, sequencing primer | 16–24 nt | 15 or 25–35 nt | < 15 or > 35 nt | Genewiz 18–24 (§4); Genewiz's own M13F universal primer is 16 nt (§5.2); warn band as above |
   | GC % (annealing part) | 40–60 | 20–40 or 60–80 | < 20 or > 80 | NEB 40–60; primer3 20–80 |
   | G/C in last 5 bases | 1–3 | 0 or 4–5 | — | primer3 `MAX_END_GC` 5; NEB "avoid GC-rich 3´ ends"; *proposal* |
   | Mononucleotide run | ≤ 4 (G ≤ 3) | 5, or G run ≥ 4 | ≥ 6 | primer3 `MAX_POLY_X` 5; IDT G run |
   | Dinucleotide repeat | < 4 repeats | ≥ 4 repeats | — | *proposal*, no source found |
   | Hairpin, dimers: structure Tm at primer3 default conditions | ≤ 47 °C | > 47 °C, not 3′-anchored | > 47 °C and 3′-anchored (END1 either way) | primer3 `_TH` 47; primer3 SELF_END rationale |
   | Hairpin, dimers: ΔG at 37 °C | > −9.0 kcal/mol | ≤ −9.0, not 3′-anchored | ≤ −9.0 and 3′-anchored | IDT −9.0; 3′ rule *proposal* |
   | 3′ end stability | report only | — | — | primer3 default off; no primary threshold found |
   | Off-target: other template sites, both strands, origin-joined, 3′ end anchored | Tm < annealing-part Tm − 10 °C | within 10 °C | would give an extra product with the partner primer inside the size range | primer3 manual (47 vs 57); *proposal* for fail |

   Every other row holds for a sequencing primer as well. Genewiz's own universal primers miss its
   45–55% GC (M13F 56%, M13-48REV 40%, §5.2), and its 50–60 °C Tm (§5.3) names no method, so it
   cannot be placed on a polymerase's scale. A colony PCR primer is a PCR primer: §3.2 gives it no
   primer rules of its own.

### For #13, bench numbers and colony PCR; #14, protocol skill

1. **PCR program from polymerase and amplicon.** Q5 default: 98 °C 30 s; 30 cycles of 98 °C 10 s,
   Ta 20 s, 72 °C at 20 s/kb (NEB allows 10 s/kb for plasmid templates); 72 °C 2 min; 4 °C. Primers
   0.5 µM, plasmid template 1 pg–10 ng per 50 µl. Round extension time up; the rounding rule is a
   *proposal*.
2. **Colony PCR default:** OneTaq Quick-Load 2X Master Mix, primers 200 nM, 25 µl (a *proposal*; NEB's
   note uses 50 µl). Touch a 1–2 mm colony with a toothpick until cloudy, then streak the same
   toothpick onto a numbered replica plate. Program: 94 °C 5 min (the upper end of NEB's 2–5 min);
   30 cycles of 94 °C 30 s, NEB OneTaq Ta for 30 s, 68 °C 1 min/kb (at least 30 s is a *proposal*);
   68 °C 5 min; 10 °C. Load 5 µl.
3. **Colony PCR primers for the smoke test:** M13/pUC Forward + M13/pUC Reverse (23-mers, Ta 51.5 °C
   in OneTaq, 137-bp empty band) instead of the 17-mer M13 pair (Ta 42.8 °C, below NEB's minimum).
   Report the correct size from the simulation (`E − D + I`) and 137 bp for empty vector. Add a
   junction primer when orientation matters.
4. **Gel:** choose from the largest and smallest expected band. Below 1 kb, including any empty-vector
   band: 2% agarose with the 100 bp ladder. 1–10 kb: 1% with the 1 kb Plus ladder. The simulated gel
   uses the ladder tables in §6.2, with 500/517 drawn as one band.
5. **ng ↔ pmol:** use NEBioCalculator's `36.04 + 615.94 × bp` as the NEB tool default. State the
   constant in the protocol, because NEB's manual formula (650 Da/bp) differs by about 5.5%.
6. **Sequencing:** design or pick Sanger primers whose 3′ end is at least 60, ideally 100, bases
   outside each junction. The M13 17-mers sit 0–26 bases from the pUC19 MCS, too close to read their
   own junction.

## Open questions

- Report the primer3 Eq. 16 integer-division defect upstream (primer3-org/primer3), and to
  libnano/primer3-py for its bundled copy.
- Typical Sanger read length and the size of the unreadable region after the primer: no primary
  source read.
- The conditions (temperature, salt, oligo concentration) behind IDT's −9.0 kcal/mol threshold.
- Why the NEB API gives Ta = lower Tm for `taqmaster-0` while the script's rule gives lower Tm − 5.
- The initial-denaturation temperature in NEB's colony PCR application note table, which is garbled.
- NEB's names for the sequencing primers (−20, −24, −47): no NEB page with sequences could be read.
- Whether NEB recommends Q5 for colony PCR.

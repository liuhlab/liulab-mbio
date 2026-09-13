"""Design PCR primers and judge them.

Tm is the SantaLucia (1998) nearest-neighbour model computed by primer3-py, salt-corrected as
NEB's Tm Calculator does it: Owczarzy (2004) at the monovalent equivalent NEB assigns the
buffer, or Schildkraut (1965) for Phusion. Each `Polymerase` carries that buffer, the rule
turning a pair's Tms into an annealing temperature (Ta), and its extension rate:

- `Q5` (the default), each primer 500 nM: Ta is the lower Tm + 1 °C, at most 72 °C, and NEB
  asks for at least 55 °C. Extension 72 °C, 20 s/kb.
- `PHUSION`, 500 nM: Ta is 0.93 * the lower Tm + 7.5 °C, at most 72 °C. Extension 72 °C,
  15 s/kb.
- `TAQ` and `ONETAQ`, 200 nM: Ta is the lower Tm - 5 °C, at most 68 °C, and NEB asks for at
  least 45 °C. Extension 68 °C, 60 s/kb.

Hairpins and dimers are structure Tms at primer3's own default conditions, which is where its
47 °C threshold comes from. Sources, and the values these reproduce:
``docs/research/primer-design-and-pcr.md``.
"""

import itertools
import math
from collections.abc import Iterable
from dataclasses import KW_ONLY, dataclass
from typing import Literal, TypedDict

from liulab_mbio.sequence import (
    BindingSite,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)

type Status = Literal["pass", "warn", "fail"]

_RANK: dict[Status, int] = {"pass": 0, "warn": 1, "fail": 2}

#: primer3 refuses a thermodynamic alignment on anything longer.
_THERMO_MAX = 60

#: SantaLucia (1998) nearest-neighbour ΔG°37, kcal/mol.
_NEAREST_NEIGHBOUR_DG = {
    "AA": -1.00, "AC": -1.44, "AG": -1.28, "AT": -0.88,
    "CA": -1.45, "CC": -1.84, "CG": -2.17, "CT": -1.28,
    "GA": -1.30, "GC": -2.24, "GG": -1.84, "GT": -1.44,
    "TA": -0.58, "TC": -1.30, "TG": -1.45, "TT": -1.00,
}  # fmt: skip

#: Duplex initiation and the penalty for each terminal A or T, as primer3 applies them.
_INITIATION_DG = 1.96
_TERMINAL_AT_DG = 0.05


@dataclass(frozen=True, slots=True)
class Polymerase:
    """A DNA polymerase, the buffer its Tm is computed in, and its cycling rules.

    Parameters
    ----------
    name
        As sold, such as ``"Q5"``.
    monovalent_mm
        The monovalent cation concentration NEB's calculator assigns this buffer, mM. It
        stands for the whole buffer, so Mg²⁺ and dNTPs are left out of the Tm: primer3's
        Owczarzy (2008) magnesium term divides by an integer and vanishes.
    salt_correction
        primer3's `salt_corrections_method`.
    primer_nm
        Each primer in the reaction, nanomolar.
    tm_dna_conc_nm
        What primer3 is given, which is four times `primer_nm` under the Owczarzy correction:
        primer3 always divides by four, where NEB divides only for Phusion.
    annealing_slope, annealing_offset
        Ta is `annealing_slope` times the lower Tm of the pair, plus `annealing_offset`, °C.
    annealing_max, annealing_min
        Ta is capped at `annealing_max`; NEB warns below `annealing_min`.
    extension_temperature
        °C.
    extension_seconds_per_kb
        Extension time per kilobase of amplicon.
    """

    name: str
    _: KW_ONLY
    monovalent_mm: float
    salt_correction: str
    primer_nm: float
    tm_dna_conc_nm: float
    annealing_slope: float
    annealing_offset: float
    annealing_max: float
    annealing_min: float
    extension_temperature: float
    extension_seconds_per_kb: int

    def annealing_temperature(self, tm: float, other_tm: float) -> float:
        """Return the annealing temperature for a pair with these Tms, °C to a tenth."""
        lower = min(tm, other_tm)
        return round(
            min(self.annealing_slope * lower + self.annealing_offset, self.annealing_max), 1
        )

    def extension_seconds(self, amplicon_length: int) -> int:
        """Return the extension time for an amplicon, rounded up to whole kilobases."""
        return max(1, math.ceil(amplicon_length / 1000)) * self.extension_seconds_per_kb


#: NEB Q5 High-Fidelity DNA Polymerase (M0491).
Q5 = Polymerase(
    "Q5",
    monovalent_mm=150.0,
    salt_correction="owczarzy",
    primer_nm=500.0,
    tm_dna_conc_nm=2000.0,
    annealing_slope=1.0,
    annealing_offset=1.0,
    annealing_max=72.0,
    annealing_min=55.0,
    extension_temperature=72.0,
    extension_seconds_per_kb=20,
)

#: NEB Phusion High-Fidelity DNA Polymerase (M0530). NEB corrects its salt the older way and
#: divides the primer concentration itself.
PHUSION = Polymerase(
    "Phusion",
    monovalent_mm=222.0,
    salt_correction="schildkraut",
    primer_nm=500.0,
    tm_dna_conc_nm=500.0,
    annealing_slope=0.93,
    annealing_offset=7.5,
    annealing_max=72.0,
    annealing_min=45.0,
    extension_temperature=72.0,
    extension_seconds_per_kb=15,
)

#: NEB Taq DNA Polymerase with Standard Taq Buffer (M0273).
TAQ = Polymerase(
    "Taq",
    monovalent_mm=55.0,
    salt_correction="owczarzy",
    primer_nm=200.0,
    tm_dna_conc_nm=800.0,
    annealing_slope=1.0,
    annealing_offset=-5.0,
    annealing_max=68.0,
    annealing_min=45.0,
    extension_temperature=68.0,
    extension_seconds_per_kb=60,
)

#: NEB OneTaq DNA Polymerase (M0480) in Standard Reaction Buffer, the colony PCR default.
ONETAQ = Polymerase(
    "OneTaq",
    monovalent_mm=54.0,
    salt_correction="owczarzy",
    primer_nm=200.0,
    tm_dna_conc_nm=800.0,
    annealing_slope=1.0,
    annealing_offset=-5.0,
    annealing_max=68.0,
    annealing_min=45.0,
    extension_temperature=68.0,
    extension_seconds_per_kb=60,
)


class _Conditions(TypedDict):
    mv_conc: float
    dv_conc: float
    dntp_conc: float
    dna_conc: float
    salt_corrections_method: str


def _conditions(polymerase: Polymerase) -> _Conditions:
    return {
        "mv_conc": polymerase.monovalent_mm,
        "dv_conc": 0.0,
        "dntp_conc": 0.0,
        "dna_conc": polymerase.tm_dna_conc_nm,
        "salt_corrections_method": polymerase.salt_correction,
    }


def melting_temperature(sequence: str, polymerase: Polymerase = Q5) -> float:
    """Return the Tm of a sequence in a polymerase's buffer, °C.

    The value NEB's Tm Calculator gives for that polymerase.

    Examples
    --------
    >>> round(melting_temperature("GTAAAACGACGGCCAGT"))
    62
    """
    import primer3

    return primer3.calc_tm(sequence.upper(), tm_method="santalucia", **_conditions(polymerase))


@dataclass(frozen=True, slots=True)
class Band:
    """The values a check passes on, and the wider band it only warns on.

    Parameters
    ----------
    low, high
        A value between them, inclusive, passes.
    warn_low, warn_high
        A value between them warns; outside them the check fails. Infinite by default, so a
        value outside the passing band only warns.
    """

    low: float
    high: float
    warn_low: float = -math.inf
    warn_high: float = math.inf

    def grade(self, value: float) -> Status:
        """Return the status of a value."""
        if self.low <= value <= self.high:
            return "pass"
        if self.warn_low <= value <= self.warn_high:
            return "warn"
        return "fail"


@dataclass(frozen=True, slots=True)
class Thresholds:
    """Every threshold the checks are judged by, with its source.

    Sources are in ``docs/research/primer-design-and-pcr.md``; where it found no published
    rule, the note's own proposal is marked below.

    Parameters
    ----------
    length
        Annealing region, nt. primer3's minimum of 18 and IDT's 18-30, warning out to the
        longest primer primer3's Tm formula takes. The warning band is the note's proposal.
    gc_percent
        Of the annealing region. NEB asks for 40-60%; primer3 allows 20-80.
    gc_clamp
        Gs and Cs in the last five 3' bases. primer3 allows five; one to three is the note's
        proposal, from NEB's "avoid GC-rich 3' ends".
    tm
        Of the annealing region, °C, on this polymerase's scale. IDT asks for 60-64 °C; the
        warning band is the note's proposal.
    tm_difference
        Between a pair's annealing regions, °C. NEB's calculator warns above five.
    mononucleotide_run, guanine_run
        Longest run of one base, and of G alone, anywhere in the primer. primer3 allows a run
        of five; IDT warns at four Gs.
    dinucleotide_repeat
        Most repeats of one dinucleotide anywhere in the primer. Four is the note's proposal;
        no source gives a limit.
    hairpin, dimer
        Melting temperature of the structure at primer3's default conditions, °C: primer3's
        `PRIMER_MAX_HAIRPIN_TH` and `PRIMER_MAX_SELF_ANY_TH`, both 47 °C.
    dimer_3prime
        As `dimer`, for a structure holding the 3' end, which the polymerase can extend:
        primer3's `PRIMER_MAX_SELF_END_TH`, failing rather than warning.
    binding_sites, products
        Places on a template where the annealing region matches, and amplicons the pair can
        make. A primer wants exactly one of each.
    off_target
        Other places where it can prime, which warn.
    binding_min_length
        Bases at the 3' end that must match for a primer carrying no binding site to be placed
        on a template. The note's proposal.
    off_target_mismatches, off_target_3prime_window, off_target_3prime_mismatches
        What is worth scoring at all: at most five mismatches, and at most one in the last five
        3' bases. Primer-BLAST's defaults (Ye et al. 2012, *BMC Bioinformatics* 13:134).
    off_target_margin
        How far under a perfect match's Tm a place can still prime, °C. primer3 sets its
        mispriming threshold 10 °C below its own minimum Tm.
    """

    length: Band = Band(18, 30, 15, 35)
    gc_percent: Band = Band(40.0, 60.0, 20.0, 80.0)
    gc_clamp: Band = Band(1, 3, 0, 5)
    tm: Band = Band(60.0, 64.0, 55.0, 70.0)
    tm_difference: Band = Band(0.0, 5.0)
    mononucleotide_run: Band = Band(0, 4, 0, 5)
    guanine_run: Band = Band(0, 3)
    dinucleotide_repeat: Band = Band(0, 3)
    hairpin: Band = Band(-math.inf, 47.0)
    dimer: Band = Band(-math.inf, 47.0)
    dimer_3prime: Band = Band(-math.inf, 47.0, -math.inf, 47.0)
    binding_sites: Band = Band(1, 1, 1, 1)
    products: Band = Band(1, 1, 1, 1)
    off_target: Band = Band(0, 0)
    binding_min_length: int = 15
    off_target_mismatches: int = 5
    off_target_3prime_window: int = 5
    off_target_3prime_mismatches: int = 1
    off_target_margin: float = 10.0


#: The thresholds every check uses unless a caller passes its own.
THRESHOLDS = Thresholds()

#: The Tm design aims for, °C: IDT's ideal, in the middle of `Thresholds.tm`.
TARGET_TM = 62.0


@dataclass(frozen=True, slots=True)
class Check:
    """One verdict on a primer or a pair, with the value it judged.

    Parameters
    ----------
    name
        What was measured, such as ``"gc_clamp"``.
    status
        ``"pass"``, ``"warn"`` or ``"fail"``, or ``None`` where no sourced threshold judges it.
    value
        The measurement, in the unit `Thresholds` documents for it.
    detail
        What a reader needs besides the number.
    """

    name: str
    status: Status | None
    value: float
    detail: str = ""


@dataclass(frozen=True, slots=True)
class _Wording:
    """How one check is written for a reader."""

    label: str
    unit: str = ""
    decimals: int = 0
    #: The `Thresholds` field holding its band, where that is not the check's own name.
    band: str = ""


#: The words every check is printed in. A name absent from here reads as itself.
_WORDING: dict[str, _Wording] = {
    "length": _Wording("length"),
    "gc_percent": _Wording("GC", "%"),
    "gc_clamp": _Wording("GC clamp"),
    "tm": _Wording("Tm", " °C", 1),
    "tm_full": _Wording("full-primer Tm", " °C", 1),
    "end_stability": _Wording("3' end stability", " kcal/mol", 1),
    "mononucleotide_run": _Wording("longest run"),
    "dinucleotide_repeat": _Wording("dinucleotide repeat"),
    "hairpin": _Wording("hairpin", " °C", 1),
    "self_dimer": _Wording("self-dimer", " °C", 1, "dimer"),
    "self_dimer_3prime": _Wording("3'-anchored self-dimer", " °C", 1, "dimer_3prime"),
    "tm_difference": _Wording("Tm difference", " °C", 1),
    "heterodimer": _Wording("heterodimer", " °C", 1, "dimer"),
    "heterodimer_3prime": _Wording("3'-anchored heterodimer", " °C", 1, "dimer_3prime"),
    "binding_sites": _Wording("binding sites"),
    "off_target": _Wording("off-target sites"),
    "products": _Wording("products"),
    "amplicon_size": _Wording("amplicon", " bp"),
}


@dataclass(frozen=True, slots=True)
class Reading:
    """One check in the few words a page prints.

    Parameters
    ----------
    label
        What a reader calls the check, such as ``"GC"``.
    value
        What it measured, with its unit, such as ``"39%"``.
    limit
        The band it was held to, such as ``"band 40-60"``, empty where nothing judged it.
    """

    label: str
    value: str
    limit: str = ""

    @property
    def detail(self) -> str:
        """The value and the band it was held to."""
        return f"{self.value} ({self.limit})" if self.limit else self.value


def reading(check: Check, thresholds: Thresholds = THRESHOLDS) -> Reading:
    """Return `check` in the words a page prints: its label, its value and its band.

    A check no sourced threshold judges carries no band, so it says only what it measured.

    Examples
    --------
    >>> reading(Check("gc_percent", "warn", 39.1)).detail
    '39% (band 40-60)'
    """
    wording = _WORDING.get(check.name, _Wording(check.name.replace("_", " ")))
    band = getattr(thresholds, wording.band or check.name, None)
    return Reading(
        wording.label,
        f"{check.value:.{wording.decimals}f}{wording.unit}",
        _band_text(band) if isinstance(band, Band) else "",
    )


def _band_text(band: Band) -> str:
    """Return the values a band passes, in a few words."""
    if band.low == band.high:
        return f"exactly {band.low:g}"
    if not math.isfinite(band.low):
        return f"max {band.high:g}"
    if not math.isfinite(band.high):
        return f"min {band.low:g}"
    return f"band {band.low:g}-{band.high:g}"


@dataclass(frozen=True, slots=True)
class PrimerReport:
    """Every check on one primer.

    Parameters
    ----------
    primer
        The primer judged.
    checks
        In the order they were run.
    """

    primer: Primer
    checks: tuple[Check, ...]

    @property
    def status(self) -> Status:
        """Return the worst status of any check that was judged."""
        return _worst(check.status for check in self.checks if check.status is not None)

    def __getitem__(self, name: str) -> Check:
        """Return the check of that name.

        Raises
        ------
        KeyError
            If no check has it.
        """
        return _named(self.checks, name)


def design_primer(
    template: SequenceRecord,
    position: int,
    strand: Strand,
    *,
    tail: str = "",
    name: str = "",
    target_tm: float = TARGET_TM,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
) -> Primer:
    """Design a primer annealing at a position on a template.

    A forward primer's annealing region starts at `position` and reads towards higher
    coordinates; a reverse primer's ends there and reads back. Either crosses the origin of a
    circular template. Of the lengths `Thresholds.length` allows, the chosen one grades best
    and then lands closest to `target_tm`. `tail` joins its 5' end and stays out of the
    binding site.

    Raises
    ------
    ValueError
        If no annealing region fits the template there.

    Examples
    --------
    >>> template = SequenceRecord("GGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCT")
    >>> design_primer(template, 0, Strand.FORWARD, name="MCS fwd").sequence
    'GGCGTAATCATGGTCATAGC'
    """
    options = _annealing_options(template, position, strand, polymerase, thresholds)
    if not options:
        raise ValueError(f"no annealing region fits this template at {position}")
    best = min(options, key=lambda option: _option_score(option, target_tm))
    return Primer(name, tail + best.sequence, binding_sites=(best.site,))


def design_pair(
    template: SequenceRecord,
    start: int,
    end: int,
    *,
    forward_tail: str = "",
    reverse_tail: str = "",
    forward_name: str = "",
    reverse_name: str = "",
    target_tm: float = TARGET_TM,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
) -> tuple[Primer, Primer]:
    """Design a pair amplifying `start` to `end`, with their Tms as near each other as they go.

    `end` passes the length of a circular template when the amplicon crosses the origin. The
    two lengths are chosen together: the pair grading best, then the Tms sitting closest to
    each other and to `target_tm`.

    Raises
    ------
    ValueError
        If no annealing region fits the template at either end.
    """
    forwards = _annealing_options(template, start, Strand.FORWARD, polymerase, thresholds)
    reverses = _annealing_options(template, end, Strand.REVERSE, polymerase, thresholds)
    if not forwards or not reverses:
        raise ValueError(f"no annealing region fits this template at {start} or at {end}")
    forward, reverse = min(
        ((one, other) for one in forwards for other in reverses),
        key=lambda pair: _pair_score(pair[0], pair[1], target_tm, thresholds),
    )
    return (
        Primer(forward_name, forward_tail + forward.sequence, binding_sites=(forward.site,)),
        Primer(reverse_name, reverse_tail + reverse.sequence, binding_sites=(reverse.site,)),
    )


def evaluate_primer(
    primer: Primer,
    template: SequenceRecord | None = None,
    *,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
) -> PrimerReport:
    """Judge one primer, and where it anneals when a template is given.

    The annealing region is the 3' part its binding site covers; a primer carrying none is
    placed on the template by `find_binding_sites`, and judged whole without one. Length, GC
    and Tm are of the annealing region, runs and structures of the whole primer. A primer
    longer than primer3 will align is judged on its 3'-terminal bases. The full-primer Tm and
    the 3'-end stability are reported without a verdict.

    Examples
    --------
    >>> report = evaluate_primer(Primer("M13 fwd", "GTAAAACGACGGCCAGT"))
    >>> report["gc_clamp"].value, report.status
    (3, 'warn')
    """
    import primer3

    sites = _placed(primer, template, thresholds)
    annealing = _annealing_region(primer, sites)
    thermo, note = _thermo_sequence(primer.sequence)
    gc = 100.0 * sum(annealing.count(base) for base in "GC") / len(annealing)
    checks = (
        _graded("length", len(annealing), thresholds.length),
        _graded("gc_percent", gc, thresholds.gc_percent),
        _graded("gc_clamp", sum(annealing[-5:].count(base) for base in "GC"), thresholds.gc_clamp),
        _graded("tm", melting_temperature(annealing, polymerase), thresholds.tm),
        Check("tm_full", None, melting_temperature(primer.sequence, polymerase)),
        Check("end_stability", None, _end_stability(annealing)),
        _run_check(primer.sequence, thresholds),
        _graded(
            "dinucleotide_repeat",
            _longest_dinucleotide_repeat(primer.sequence),
            thresholds.dinucleotide_repeat,
        ),
        _graded("hairpin", primer3.calc_hairpin(thermo).tm, thresholds.hairpin, note),
        _graded("self_dimer", primer3.calc_homodimer(thermo).tm, thresholds.dimer, note),
        _graded(
            "self_dimer_3prime",
            primer3.calc_end_stability(thermo, thermo).tm,
            thresholds.dimer_3prime,
            note,
        ),
    )
    if template is not None:
        checks += _template_checks(annealing, sites, template, thresholds)
    return PrimerReport(primer, checks)


@dataclass(frozen=True, slots=True)
class PairReport:
    """Every check on a primer pair, and the numbers its PCR needs.

    Parameters
    ----------
    forward, reverse
        What each primer scored on its own.
    checks
        What only a pair can be judged on.
    annealing_temperature
        °C, by the polymerase's rule over the two annealing regions.
    amplicon_length
        Bases between the primers' 5' ends, tails included, or ``None`` unless they make
        exactly one product.
    extension_seconds
        For that amplicon, or ``None``.
    """

    forward: PrimerReport
    reverse: PrimerReport
    checks: tuple[Check, ...]
    annealing_temperature: float
    amplicon_length: int | None
    extension_seconds: int | None

    @property
    def status(self) -> Status:
        """Return the worst status of either primer or of any pair check that was judged."""
        judged = [check.status for check in self.checks if check.status is not None]
        return _worst((self.forward.status, self.reverse.status, *judged))

    def __getitem__(self, name: str) -> Check:
        """Return the pair check of that name.

        Raises
        ------
        KeyError
            If no check has it.
        """
        return _named(self.checks, name)


def evaluate_pair(
    forward: Primer,
    reverse: Primer,
    template: SequenceRecord,
    *,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
) -> PairReport:
    """Judge a primer pair on a template, with the amplicon it makes.

    Each primer is judged as on its own, and the pair adds the gap between their Tms, their
    heterodimers, and every product their binding sites can make. An amplicon runs from one
    primer's 5' end to the other's, tails counted, across the origin where it must.

    Examples
    --------
    >>> template = SequenceRecord("GGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCT")
    >>> pair = design_pair(template, 0, len(template))
    >>> report = evaluate_pair(pair[0], pair[1], template)
    >>> report.amplicon_length, report["products"].status
    (48, 'pass')
    """
    import primer3

    reports = tuple(
        evaluate_primer(primer, template, polymerase=polymerase, thresholds=thresholds)
        for primer in (forward, reverse)
    )
    tms = tuple(report["tm"].value for report in reports)
    first, second = (_thermo_sequence(primer.sequence) for primer in (forward, reverse))
    note = first[1] or second[1]
    anchored = max(
        primer3.calc_end_stability(first[0], second[0]).tm,
        primer3.calc_end_stability(second[0], first[0]).tm,
    )
    products = _products(forward, reverse, template, thresholds)
    length = products[0] if len(products) == 1 else None
    checks = (
        _graded("tm_difference", abs(tms[0] - tms[1]), thresholds.tm_difference),
        _graded(
            "heterodimer", primer3.calc_heterodimer(first[0], second[0]).tm, thresholds.dimer, note
        ),
        _graded("heterodimer_3prime", anchored, thresholds.dimer_3prime, note),
        _graded(
            "products",
            len(products),
            thresholds.products,
            ", ".join(f"{size} bp" for size in products),
        ),
        Check("amplicon_size", None, length or 0),
    )
    return PairReport(
        reports[0],
        reports[1],
        checks,
        polymerase.annealing_temperature(*tms),
        length,
        None if length is None else polymerase.extension_seconds(length),
    )


@dataclass(frozen=True, slots=True)
class PrimingSite:
    """Where a primer's 3' end can anneal on a template, and how strongly.

    Parameters
    ----------
    site
        What the annealing region covers, running past the end of a circular template when it
        crosses the origin.
    tm
        Of the primer annealed there, °C at primer3's default conditions.
    mismatches
        Bases of the annealing region that do not pair.
    """

    site: BindingSite
    tm: float
    mismatches: int


def find_binding_sites(
    sequence: str, template: SequenceRecord, *, thresholds: Thresholds = THRESHOLDS
) -> tuple[BindingSite, ...]:
    """Return every place a primer's 3' end matches a template exactly, on either strand.

    The match runs from the 3' end back towards the 5' end, so a tail hangs off it; it must
    reach `Thresholds.binding_min_length`. A site crosses the origin of a circular template.

    Examples
    --------
    >>> template = SequenceRecord("CGGCGTAATCATGGTCATAGCTGTTTCC")
    >>> sites = find_binding_sites("TTGGTCTCAGGCGTAATCATGGTCATAGC", template)
    >>> [(site.start, site.end, site.strand.name) for site in sites]
    [(1, 21, 'FORWARD')]
    """
    dna = sequence.upper()
    length = len(template)
    reverse = reverse_complement(dna)
    sites = []
    for index in range(length):
        matched = _matched(template, index, dna[::-1], -1)
        if matched >= thresholds.binding_min_length:
            start = (index - matched + 1) % length
            sites.append(BindingSite(start, start + matched, Strand.FORWARD))
        matched = _matched(template, index, reverse, 1)
        if matched >= thresholds.binding_min_length:
            sites.append(BindingSite(index, index + matched, Strand.REVERSE))
    return tuple(sorted(sites, key=lambda site: (site.start, site.strand)))


def find_priming_sites(
    sequence: str, template: SequenceRecord, *, thresholds: Thresholds = THRESHOLDS
) -> tuple[PrimingSite, ...]:
    """Return every place an annealing region can prime a template, on either strand.

    A place counts when few enough of its bases mismatch, fewest of all at the 3' end, and
    when the primer annealed there melts within `Thresholds.off_target_margin` of a perfect
    match. Sites cross the origin of a circular template.
    """
    import primer3

    dna = sequence.upper()
    size = len(dna)
    length = len(template)
    if size > length:
        return ()
    circular = template.topology == "circular"
    top = template.sequence + (template.sequence[: size - 1] if circular else "")
    reverse = reverse_complement(dna)
    window = thresholds.off_target_3prime_window
    floor = primer3.calc_end_stability(dna, reverse).tm - thresholds.off_target_margin
    found = []
    for start in range(length if circular else length - size + 1):
        here = top[start : start + size]
        for strand, probe, anchor, annealed in (
            (Strand.FORWARD, dna, _mismatches(dna[-window:], here[-window:]), None),
            (Strand.REVERSE, reverse, _mismatches(reverse[:window], here[:window]), here),
        ):
            if anchor > thresholds.off_target_3prime_mismatches:
                continue
            mismatches = _mismatches(probe, here)
            if mismatches > thresholds.off_target_mismatches:
                continue
            tm = primer3.calc_end_stability(dna, annealed or reverse_complement(here)).tm
            if tm >= floor:
                found.append(PrimingSite(BindingSite(start, start + size, strand), tm, mismatches))
    return tuple(found)


def amplicon_sizes(
    forward: Primer,
    reverse: Primer,
    template: SequenceRecord,
    *,
    thresholds: Thresholds = THRESHOLDS,
) -> tuple[int, ...]:
    """Return the size of every amplicon a pair can make on a template, smallest first.

    Tails count, as they do in `PairReport.amplicon_length`. A primer carrying binding sites is
    taken to bind where they say, so asking about a template other than the one it was placed on
    means clearing them first.
    """
    return tuple(_products(forward, reverse, template, thresholds))


def _matched(template: SequenceRecord, index: int, probe: str, step: int) -> int:
    """Return how many bases of `probe` match the template from `index`, walking by `step`."""
    length = len(template)
    circular = template.topology == "circular"
    matched = 0
    while matched < min(len(probe), length):
        position = index + step * matched
        if not circular and not 0 <= position < length:
            break
        if template.sequence[position % length] != probe[matched]:
            break
        matched += 1
    return matched


def _mismatches(one: str, other: str) -> int:
    return sum(base != base_here for base, base_here in zip(one, other, strict=True))


def _placed(
    primer: Primer, template: SequenceRecord | None, thresholds: Thresholds
) -> tuple[BindingSite, ...]:
    if primer.binding_sites or template is None:
        return primer.binding_sites
    return find_binding_sites(primer.sequence, template, thresholds=thresholds)


def _template_checks(
    annealing: str,
    sites: tuple[BindingSite, ...],
    template: SequenceRecord,
    thresholds: Thresholds,
) -> tuple[Check, ...]:
    length = len(template)
    intended = {_three_prime_end(site, length) for site in sites}
    elsewhere = [
        site
        for site in find_priming_sites(annealing, template, thresholds=thresholds)
        if _three_prime_end(site.site, length) not in intended
    ]
    detail = ", ".join(
        f"{site.site.start} {site.site.strand.name.lower()}, {site.mismatches} mismatched"
        for site in elsewhere
    )
    return (
        _graded("binding_sites", len(sites), thresholds.binding_sites),
        _graded("off_target", len(elsewhere), thresholds.off_target, detail),
    )


def _three_prime_end(site: BindingSite, length: int) -> tuple[int, Strand]:
    end = site.start if site.strand is Strand.REVERSE else site.end - 1
    return end % length, site.strand


@dataclass(frozen=True, slots=True)
class _Option:
    site: BindingSite
    sequence: str
    tm: float
    grade: Status


def _annealing_options(
    template: SequenceRecord,
    position: int,
    strand: Strand,
    polymerase: Polymerase,
    thresholds: Thresholds,
) -> list[_Option]:
    length = len(template)
    circular = template.topology == "circular"
    if not circular and not 0 <= position <= length:
        return []
    shortest = int(_finite(thresholds.length.warn_low, thresholds.length.low))
    longest = min(int(_finite(thresholds.length.warn_high, thresholds.length.high)), length)
    options = []
    for size in range(shortest, longest + 1):
        start = position if strand is Strand.FORWARD else position - size
        if circular:
            start %= length
        elif not 0 <= start <= length - size:
            continue
        bases = template.extract(Segment(start, start + size))
        sequence = bases if strand is Strand.FORWARD else reverse_complement(bases)
        tm = melting_temperature(sequence, polymerase)
        options.append(
            _Option(
                BindingSite(start, start + size, strand),
                sequence,
                tm,
                _worst((thresholds.length.grade(size), thresholds.tm.grade(tm))),
            )
        )
    return options


def _option_score(option: _Option, target_tm: float) -> tuple[int, float, int]:
    return _RANK[option.grade], abs(option.tm - target_tm), len(option.sequence)


def _pair_score(
    forward: _Option, reverse: _Option, target_tm: float, thresholds: Thresholds
) -> tuple[int, float, int]:
    difference = abs(forward.tm - reverse.tm)
    grade = _worst((forward.grade, reverse.grade, thresholds.tm_difference.grade(difference)))
    drift = max(abs(forward.tm - target_tm), abs(reverse.tm - target_tm))
    return _RANK[grade], drift + difference, len(forward.sequence) + len(reverse.sequence)


def _finite(value: float, fallback: float) -> float:
    return value if math.isfinite(value) else fallback


def _named(checks: tuple[Check, ...], name: str) -> Check:
    for check in checks:
        if check.name == name:
            return check
    raise KeyError(name)


def _products(
    forward: Primer, reverse: Primer, template: SequenceRecord, thresholds: Thresholds
) -> list[int]:
    """Return the size of every amplicon the two primers' sites can make, tails included."""
    placed: list[tuple[BindingSite, int]] = []
    for primer in (forward, reverse):
        sites = _placed(primer, template, thresholds)
        annealing = max((site.end - site.start for site in sites), default=len(primer.sequence))
        placed.extend((site, len(primer.sequence) - annealing) for site in sites)
    products = []
    for site, tail in placed:
        if site.strand is not Strand.FORWARD:
            continue
        for other, other_tail in placed:
            if other.strand is not Strand.REVERSE:
                continue
            span = _span(site, other, template)
            if span is not None:
                products.append(span + tail + other_tail)
    return sorted(products)


def _span(forward: BindingSite, reverse: BindingSite, template: SequenceRecord) -> int | None:
    length = len(template)
    if template.topology == "circular":
        return (reverse.end - forward.start) % length or length
    span = reverse.end - forward.start
    return span if span > 0 else None


def _graded(name: str, value: float, band: Band, detail: str = "") -> Check:
    return Check(name, band.grade(value), value, detail)


def _worst(statuses: Iterable[Status]) -> Status:
    worst: Status = "pass"
    for status in statuses:
        if _RANK[status] > _RANK[worst]:
            worst = status
    return worst


def _run_check(sequence: str, thresholds: Thresholds) -> Check:
    longest = _longest_run(sequence)
    guanines = _longest_run(sequence, "G")
    guanine_status = thresholds.guanine_run.grade(guanines)
    return Check(
        "mononucleotide_run",
        _worst((thresholds.mononucleotide_run.grade(longest), guanine_status)),
        longest,
        f"{guanines} Gs in a row" if guanine_status != "pass" else "",
    )


def _annealing_region(primer: Primer, sites: tuple[BindingSite, ...] = ()) -> str:
    sites = sites or primer.binding_sites
    if not sites:
        return primer.sequence
    return primer.sequence[-max(site.end - site.start for site in sites) :]


def _thermo_sequence(sequence: str) -> tuple[str, str]:
    if len(sequence) <= _THERMO_MAX:
        return sequence, ""
    return sequence[-_THERMO_MAX:], f"judged on the 3'-terminal {_THERMO_MAX} bases"


def _end_stability(sequence: str) -> float:
    """Return ΔG°37 of the 3'-terminal pentamer, kcal/mol. Bases outside ACGT add nothing."""
    pentamer = sequence[-5:]
    total = _INITIATION_DG + sum(
        _NEAREST_NEIGHBOUR_DG.get(pair, 0.0)
        for pair in (pentamer[i : i + 2] for i in range(len(pentamer) - 1))
    )
    return total + _TERMINAL_AT_DG * sum(base in "AT" for base in (pentamer[0], pentamer[-1]))


def _longest_run(sequence: str, base: str = "") -> int:
    longest = run = 1 if not base else 0
    for previous, this in itertools.pairwise(sequence):
        run = run + 1 if this == previous else 1
        if not base or this == base == previous:
            longest = max(longest, run)
    return longest


def _longest_dinucleotide_repeat(sequence: str) -> int:
    longest = 1
    for start in range(len(sequence) - 1):
        unit = sequence[start : start + 2]
        if unit[0] == unit[1]:
            continue
        repeats = 1
        while sequence[start + 2 * repeats : start + 2 * repeats + 2] == unit:
            repeats += 1
        longest = max(longest, repeats)
    return longest

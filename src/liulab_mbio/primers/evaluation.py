"""Judge one primer, or a pair on a template, against the thresholds.

Hairpins and dimers are structure Tms at primer3's own default conditions, which is where its
47 °C threshold comes from. Sources, and the values these reproduce:
``docs/research/primer-design-and-pcr.md``.
"""

import itertools
from dataclasses import dataclass
from functools import lru_cache

from liulab_mbio.checks import Check, Status, worst, worst_of
from liulab_mbio.primers.placement import amplicon_sizes, find_binding_sites, find_priming_sites
from liulab_mbio.primers.polymerase import Q5, Polymerase, melting_temperature
from liulab_mbio.primers.thresholds import THRESHOLDS, Band, Thresholds
from liulab_mbio.sequence import BindingSite, Primer, SequenceRecord, Strand, Topology

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
        return worst_of(self.checks)

    def __getitem__(self, name: str) -> Check:
        """Return the check of that name.

        Raises
        ------
        KeyError
            If no check has it.
        """
        return _named(self.checks, name)


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
    bases = None if template is None else (template.sequence, template.topology)
    checks = _checks(primer.sequence, primer.binding_sites, bases, polymerase, thresholds)
    return PrimerReport(primer, checks)


@lru_cache(maxsize=1 << 13)
def _checks(
    sequence: str,
    sites: tuple[BindingSite, ...],
    bases: tuple[str, Topology] | None,
    polymerase: Polymerase,
    thresholds: Thresholds,
) -> tuple[Check, ...]:
    """Return every check `evaluate_primer` runs on one primer.

    A design judges the same candidates on the same template again and again, so these are
    remembered: they depend on nothing but the primer's bases and binding sites, the template's
    bases and topology, the polymerase and the thresholds.
    """
    template = None if bases is None else _template(*bases)
    if not sites and template is not None:
        sites = find_binding_sites(sequence, template, thresholds=thresholds)
    annealing = sequence[-max(site.end - site.start for site in sites) :] if sites else sequence
    thermo, note = _thermo_sequence(sequence)
    hairpin, self_dimer, anchored_dimer = _structure_tms(thermo)
    checks = (
        *annealing_checks(annealing, polymerase=polymerase, thresholds=thresholds),
        Check("tm_full", None, melting_temperature(sequence, polymerase)),
        Check("end_stability", None, _end_stability(annealing)),
        _run_check(sequence, thresholds),
        _graded(
            "dinucleotide_repeat",
            _longest_dinucleotide_repeat(sequence),
            thresholds.dinucleotide_repeat,
        ),
        _graded("hairpin", hairpin, thresholds.hairpin, note),
        _graded("self_dimer", self_dimer, thresholds.dimer, note),
        _graded("self_dimer_3prime", anchored_dimer, thresholds.dimer_3prime, note),
    )
    if template is not None:
        checks += _template_checks(annealing, sites, template, thresholds)
    return checks


@lru_cache(maxsize=8)
def _template(sequence: str, topology: Topology) -> SequenceRecord:
    """Return a bare record of a template's bases, which is all a primer is judged against."""
    return SequenceRecord(sequence, topology=topology)


def annealing_checks(
    annealing: str, *, polymerase: Polymerase = Q5, thresholds: Thresholds = THRESHOLDS
) -> tuple[Check, ...]:
    """Judge what an annealing region's length decides: its length, GC, GC clamp and Tm.

    A design judges every candidate by `evaluate_primer` instead; these are the four checks
    that its length alone moves.

    Examples
    --------
    >>> [check.status for check in annealing_checks("GTAAAACGACGGCCAGT")]
    ['warn', 'pass', 'pass', 'pass']
    """
    gc = 100.0 * sum(annealing.count(base) for base in "GC") / len(annealing)
    return (
        _graded("length", len(annealing), thresholds.length),
        _graded("gc_percent", gc, thresholds.gc_percent),
        _graded("gc_clamp", sum(annealing[-5:].count(base) for base in "GC"), thresholds.gc_clamp),
        _graded("tm", melting_temperature(annealing, polymerase), thresholds.tm),
    )


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
        return worst(
            (self.forward.status, self.reverse.status, *(one.status for one in self.checks))
        )

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
    >>> from liulab_mbio.primers.design import design_pair
    >>> template = SequenceRecord("GGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCT")
    >>> pair = design_pair(template, 0, len(template))
    >>> report = evaluate_pair(pair[0], pair[1], template)
    >>> report.amplicon_length, report["products"].status
    (48, 'pass')
    """
    reports = tuple(
        evaluate_primer(primer, template, polymerase=polymerase, thresholds=thresholds)
        for primer in (forward, reverse)
    )
    checks = pair_checks(reports[0], reports[1], template, thresholds=thresholds)
    tms = tuple(report["tm"].value for report in reports)
    length = int(checks[-1].value) or None
    return PairReport(
        reports[0],
        reports[1],
        checks,
        polymerase.annealing_temperature(*tms),
        length,
        None if length is None else polymerase.extension_seconds(length),
    )


def pair_checks(
    forward: PrimerReport,
    reverse: PrimerReport,
    template: SequenceRecord,
    *,
    thresholds: Thresholds = THRESHOLDS,
) -> tuple[Check, ...]:
    """Judge what only a pair can be judged on, from what each primer already scored.

    The gap between their annealing Tms, their heterodimers, and every product their binding
    sites can make. The last check, the amplicon size, carries no verdict and is zero unless
    the pair makes exactly one product. `evaluate_pair` is the whole judgement; a design that
    has judged each primer already calls this.
    """
    import primer3

    primers = (forward.primer, reverse.primer)
    tms = (forward["tm"].value, reverse["tm"].value)
    first, second = (_thermo_sequence(primer.sequence) for primer in primers)
    note = first[1] or second[1]
    anchored = max(
        primer3.calc_end_stability(first[0], second[0]).tm,
        primer3.calc_end_stability(second[0], first[0]).tm,
    )
    products = amplicon_sizes(*primers, template, thresholds=thresholds)
    return (
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
        Check("amplicon_size", None, products[0] if len(products) == 1 else 0),
    )


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


def _named(checks: tuple[Check, ...], name: str) -> Check:
    for check in checks:
        if check.name == name:
            return check
    raise KeyError(name)


def _graded(name: str, value: float, band: Band, detail: str = "") -> Check:
    return Check(name, band.grade(value), value, detail)


def _run_check(sequence: str, thresholds: Thresholds) -> Check:
    longest = _longest_run(sequence)
    guanines = _longest_run(sequence, "G")
    guanine_status = thresholds.guanine_run.grade(guanines)
    return Check(
        "mononucleotide_run",
        worst((thresholds.mononucleotide_run.grade(longest), guanine_status)),
        longest,
        f"{guanines} Gs in a row" if guanine_status != "pass" else "",
    )


def _thermo_sequence(sequence: str) -> tuple[str, str]:
    if len(sequence) <= _THERMO_MAX:
        return sequence, ""
    return sequence[-_THERMO_MAX:], f"judged on the 3'-terminal {_THERMO_MAX} bases"


@lru_cache(maxsize=1 << 15)
def _structure_tms(sequence: str) -> tuple[float, float, float]:
    """Return the hairpin, self-dimer and 3'-anchored self-dimer Tms of one sequence, °C.

    A design judges thousands of candidates, so these are remembered: each depends on nothing
    but the sequence.
    """
    import primer3

    return (
        primer3.calc_hairpin(sequence).tm,
        primer3.calc_homodimer(sequence).tm,
        primer3.calc_end_stability(sequence, sequence).tm,
    )


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

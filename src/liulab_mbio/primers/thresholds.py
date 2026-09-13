"""The bands every primer check is judged by, and the words a check is printed in."""

import math
from collections.abc import Mapping
from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import Literal

from liulab_mbio.checks import Check, Status


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
    proposed
        Whether the passing band is the research note's own proposal rather than a published
        rule. A page prints it as proposed.
    """

    low: float
    high: float
    warn_low: float = -math.inf
    warn_high: float = math.inf
    proposed: bool = False

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
    off_target_amplicons
        Amplicons a pair, or either primer alone, makes on a genome besides the intended one,
        which warn as an off-target site does.
    intended_amplicon
        The amplicon named as intended, which a genome must hold.
    genome_size_limit, genome_mismatches, genome_mismatches_above_limit, genome_terminal_window
        How near a match a genome search finds: up to 3 mismatches, none in the 3'-terminal 3
        bases, on a FASTA up to 110 MB on disk, and a perfect match only on a larger one. On
        hg38 one mismatch took minutes where a perfect match took seconds; the limit takes in
        the 100 Mb worm genome, searched with 3 mismatches in under a second, and an E. coli
        genome beside it (#44, #47).
    genome_max_amplicon
        The longest amplicon a genome search reports, bp: the length #44 measured with.
    """

    length: Band = Band(18, 30, 15, 35)
    gc_percent: Band = Band(40.0, 60.0, 20.0, 80.0)
    gc_clamp: Band = Band(1, 3, 0, 5, proposed=True)
    tm: Band = Band(60.0, 64.0, 55.0, 70.0)
    tm_difference: Band = Band(0.0, 5.0)
    mononucleotide_run: Band = Band(0, 4, 0, 5)
    guanine_run: Band = Band(0, 3)
    dinucleotide_repeat: Band = Band(0, 3, proposed=True)
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
    off_target_amplicons: Band = Band(0, 0)
    intended_amplicon: Band = Band(1, 1)
    genome_size_limit: int = 110_000_000
    genome_mismatches: int = 3
    genome_mismatches_above_limit: int = 0
    genome_terminal_window: int = 3
    genome_max_amplicon: int = 4000


#: The thresholds every check uses unless a caller passes its own: a PCR primer's.
THRESHOLDS = Thresholds()

#: What a primer is for, which chooses the thresholds it is designed and judged by.
type PrimerRole = Literal["amplification", "colony PCR", "sequencing"]

#: The thresholds for each role. A colony PCR primer is held to a PCR primer's. A sequencing
#: primer passes on 16-24 bases, warning beyond them as `Thresholds.length` does: Genewiz asks
#: for 18-24, and its own M13F universal primer is 16. The note's implications for primer
#: design say why its other bands stay a PCR primer's.
THRESHOLDS_FOR: Mapping[PrimerRole, Thresholds] = MappingProxyType(
    {
        "amplification": THRESHOLDS,
        "colony PCR": THRESHOLDS,
        "sequencing": replace(THRESHOLDS, length=Band(16, 24, 15, 35)),
    }
)

#: The Tm design aims for, °C: IDT's ideal, in the middle of `Thresholds.tm`.
TARGET_TM = 62.0


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
    "off_target_amplicons": _Wording("off-target amplicons"),
    "intended_amplicon": _Wording("intended amplicon"),
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
    """Return the values a band passes, in a few words, saying so where it is only proposed."""
    if band.low == band.high:
        text = f"exactly {band.low:g}"
    elif not math.isfinite(band.low):
        text = f"max {band.high:g}"
    elif not math.isfinite(band.high):
        text = f"min {band.low:g}"
    else:
        text = f"band {band.low:g}-{band.high:g}"
    return f"proposed {text}" if band.proposed else text

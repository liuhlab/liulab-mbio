"""A DNA polymerase: the buffer its Tm is computed in, its cycling rules, and NEB's PCR for it.

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

Each also carries its `PcrProfile`: the buffer, the enzyme and the cycling NEB's protocol sets.
Sources, and the values these reproduce: ``docs/research/primer-design-and-pcr.md``.
"""

import math
from dataclasses import KW_ONLY, dataclass
from typing import TypedDict


@dataclass(frozen=True, slots=True)
class PcrProfile:
    """What NEB's protocol puts in one polymerase's PCR, and the program around it.

    Parameters
    ----------
    buffer_name, buffer_fold
        The reaction buffer as supplied, so its volume is the reaction over `buffer_fold`.
    units_per_ul
        Polymerase in the reaction.
    stock_units_ul
        Polymerase in the tube it is pipetted from.
    initial_denaturation_c, initial_denaturation_seconds
        The step before the cycles.
    denaturation_c, denaturation_seconds, annealing_seconds
        Inside each cycle; the annealing temperature is the primer pair's.
    final_extension_seconds
        At the polymerase's own extension temperature.
    two_step_celsius
        The lowest annealing temperature that gets a two-step program, annealing and extension
        combined. Taq's rule is "above 65 °C", which at a tenth of a degree is 65.1.
    cycles, dntp_um_each, hold_c
        The rest of NEB's table.
    """

    buffer_name: str
    _: KW_ONLY
    buffer_fold: float
    units_per_ul: float
    stock_units_ul: float
    initial_denaturation_c: float
    denaturation_c: float
    denaturation_seconds: int
    annealing_seconds: int
    final_extension_seconds: int
    two_step_celsius: float
    cycles: int = 30
    initial_denaturation_seconds: int = 30
    dntp_um_each: float = 200.0
    hold_c: float = 4.0


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
    pcr
        Its buffer, enzyme dose and cycling, as NEB's protocol for it sets them.
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
    pcr: PcrProfile

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
    pcr=PcrProfile(
        "Q5 Reaction Buffer",
        buffer_fold=5.0,
        units_per_ul=0.02,
        stock_units_ul=2.0,
        initial_denaturation_c=98.0,
        denaturation_c=98.0,
        denaturation_seconds=10,
        annealing_seconds=20,
        final_extension_seconds=120,
        two_step_celsius=72.0,
    ),
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
    pcr=PcrProfile(
        "Phusion HF Buffer",
        buffer_fold=5.0,
        units_per_ul=0.02,
        stock_units_ul=2.0,
        initial_denaturation_c=98.0,
        denaturation_c=98.0,
        denaturation_seconds=10,
        annealing_seconds=20,
        final_extension_seconds=300,
        two_step_celsius=72.0,
    ),
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
    pcr=PcrProfile(
        "Standard Taq Reaction Buffer",
        buffer_fold=10.0,
        units_per_ul=0.025,
        stock_units_ul=5.0,
        initial_denaturation_c=95.0,
        denaturation_c=95.0,
        denaturation_seconds=30,
        annealing_seconds=30,
        final_extension_seconds=300,
        two_step_celsius=65.1,
    ),
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
    pcr=PcrProfile(
        "OneTaq Standard Reaction Buffer",
        buffer_fold=5.0,
        units_per_ul=0.025,
        stock_units_ul=5.0,
        initial_denaturation_c=94.0,
        denaturation_c=94.0,
        denaturation_seconds=30,
        annealing_seconds=30,
        final_extension_seconds=300,
        two_step_celsius=68.0,
    ),
)

#: Every polymerase this package ships.
POLYMERASES: tuple[Polymerase, ...] = (Q5, PHUSION, TAQ, ONETAQ)


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

"""A DNA polymerase: the buffer its Tm is computed in, and its annealing and extension rules.

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

Sources, and the values these reproduce: ``docs/research/primer-design-and-pcr.md``.
"""

import math
from dataclasses import KW_ONLY, dataclass
from typing import TypedDict


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

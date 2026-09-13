"""Design PCR primers and judge them.

Tm is computed by primer3-py with the SantaLucia (1998) nearest-neighbour parameters and salt
correction; Mg²⁺ and dNTPs enter as a monovalent equivalent (von Ahsen et al. 2001), as primer3
does. A `Polymerase` holds the buffer a Tm is computed in and the rule that turns a pair's Tms
into an annealing temperature (Ta). Buffers are monovalent, Mg²⁺, dNTPs and each primer:

- `Q5`: 50 mM, 2.0 mM, 0.8 mM, 500 nM. Ta is the lower Tm + 3 °C, at most 72 °C.
- `TAQ`: 50 mM, 1.5 mM, 0.8 mM, 200 nM. Ta is the lower Tm − 5 °C, at most 68 °C.
- `ONETAQ`: 44 mM, 1.8 mM, 0.8 mM, 200 nM. Ta as for `TAQ`.
"""

import math
from dataclasses import KW_ONLY, dataclass


@dataclass(frozen=True, slots=True)
class Polymerase:
    """A DNA polymerase, its reaction buffer and its cycling rules.

    Parameters
    ----------
    name
        As sold, such as ``"Q5"``.
    monovalent_mm, divalent_mm, dntp_mm
        Final concentrations, millimolar. dNTPs are summed over all four.
    primer_nm
        Final concentration of each primer, nanomolar.
    annealing_offset
        Added to the lower Tm of a pair to give the annealing temperature, °C.
    annealing_max
        The highest annealing temperature, °C.
    extension_temperature
        °C.
    extension_seconds_per_kb
        Extension time per kilobase of amplicon.
    """

    name: str
    _: KW_ONLY
    monovalent_mm: float
    divalent_mm: float
    dntp_mm: float
    primer_nm: float
    annealing_offset: float
    annealing_max: float
    extension_temperature: float
    extension_seconds_per_kb: int

    def annealing_temperature(self, tm: float, other_tm: float) -> float:
        """Return the annealing temperature for a primer pair with these Tms, °C."""
        return min(min(tm, other_tm) + self.annealing_offset, self.annealing_max)

    def extension_seconds(self, amplicon_length: int) -> int:
        """Return the extension time for an amplicon, rounded up to whole kilobases."""
        return max(1, math.ceil(amplicon_length / 1000)) * self.extension_seconds_per_kb


#: NEB Q5 High-Fidelity DNA Polymerase (M0491) protocol; 30 s/kb is its upper figure. NEB does
#: not publish the buffer's monovalent salt, so that is primer3's default.
Q5 = Polymerase(
    "Q5",
    monovalent_mm=50.0,
    divalent_mm=2.0,
    dntp_mm=0.8,
    primer_nm=500.0,
    annealing_offset=3.0,
    annealing_max=72.0,
    extension_temperature=72.0,
    extension_seconds_per_kb=30,
)

#: NEB Taq DNA Polymerase with Standard Taq Buffer (M0273): 50 mM KCl, 1.5 mM MgCl₂.
TAQ = Polymerase(
    "Taq",
    monovalent_mm=50.0,
    divalent_mm=1.5,
    dntp_mm=0.8,
    primer_nm=200.0,
    annealing_offset=-5.0,
    annealing_max=68.0,
    extension_temperature=68.0,
    extension_seconds_per_kb=60,
)

#: NEB OneTaq DNA Polymerase (M0480), Standard Reaction Buffer: 22 mM KCl, 22 mM NH₄Cl,
#: 1.8 mM MgCl₂.
ONETAQ = Polymerase(
    "OneTaq",
    monovalent_mm=44.0,
    divalent_mm=1.8,
    dntp_mm=0.8,
    primer_nm=200.0,
    annealing_offset=-5.0,
    annealing_max=68.0,
    extension_temperature=68.0,
    extension_seconds_per_kb=60,
)


def _conditions(polymerase: Polymerase) -> dict[str, float]:
    return {
        "mv_conc": polymerase.monovalent_mm,
        "dv_conc": polymerase.divalent_mm,
        "dntp_conc": polymerase.dntp_mm,
        "dna_conc": polymerase.primer_nm,
    }


def melting_temperature(sequence: str, polymerase: Polymerase = Q5) -> float:
    """Return the Tm of a sequence paired with its complement, in a polymerase's buffer, °C.

    Examples
    --------
    >>> round(melting_temperature("GTAAAACGACGGCCAGT"))
    59
    """
    import primer3

    return primer3.calc_tm(
        sequence.upper(),
        **_conditions(polymerase),
        tm_method="santalucia",
        salt_corrections_method="santalucia",
    )

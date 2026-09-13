"""Design PCR primers and judge them.

Five modules, and every public name in them imports from here as well:

- `polymerase`: a polymerase, the Tm in its buffer, and its annealing and extension rules.
- `thresholds`: the bands a check is judged by, and the words it is printed in.
- `placement`: where a primer anneals on a template, where else it primes, and amplicons.
- `evaluation`: judging a primer, or a pair on a template.
- `design`: choosing the annealing region of a primer or a pair.
"""

from liulab_mbio.primers.design import design_pair, design_primer
from liulab_mbio.primers.evaluation import PairReport, PrimerReport, evaluate_pair, evaluate_primer
from liulab_mbio.primers.placement import (
    PrimingSite,
    amplicon_sizes,
    find_binding_sites,
    find_priming_sites,
)
from liulab_mbio.primers.polymerase import (
    ONETAQ,
    PHUSION,
    Q5,
    TAQ,
    Polymerase,
    melting_temperature,
)
from liulab_mbio.primers.thresholds import (
    TARGET_TM,
    THRESHOLDS,
    Band,
    Reading,
    Thresholds,
    reading,
)

__all__ = [
    "ONETAQ",
    "PHUSION",
    "Q5",
    "TAQ",
    "TARGET_TM",
    "THRESHOLDS",
    "Band",
    "PairReport",
    "Polymerase",
    "PrimerReport",
    "PrimingSite",
    "Reading",
    "Thresholds",
    "amplicon_sizes",
    "design_pair",
    "design_primer",
    "evaluate_pair",
    "evaluate_primer",
    "find_binding_sites",
    "find_priming_sites",
    "melting_temperature",
    "reading",
]

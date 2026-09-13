"""Design PCR primers and judge them.

Every public name in its modules imports from here as well:

- `polymerase`: a polymerase, the Tm in its buffer, its annealing and extension rules, and
  NEB's PCR for it.
- `thresholds`: the bands a check is judged by, and the words it is printed in.
- `placement`: where a primer may anneal, where it does, where else it primes, and amplicons.
- `evaluation`: judging a primer, or a pair on a template, and what a primer's length decides.
- `design`: choosing the annealing region of a primer or a pair.
- `genome`: checking pairs against a genome FASTA, and the amplicons they make there.
"""

from liulab_mbio.primers.design import design_pair, design_primer
from liulab_mbio.primers.evaluation import (
    PairReport,
    PrimerReport,
    annealing_checks,
    evaluate_pair,
    evaluate_primer,
    pair_checks,
)
from liulab_mbio.primers.genome import (
    Amplicon,
    GenomeReport,
    Locus,
    MadeBy,
    evaluate_on_genome,
    evaluate_pair_on_genome,
)
from liulab_mbio.primers.placement import (
    Placement,
    PrimingSite,
    amplicon_sizes,
    find_binding_sites,
    find_priming_sites,
)
from liulab_mbio.primers.polymerase import (
    ONETAQ,
    PHUSION,
    POLYMERASES,
    Q5,
    TAQ,
    PcrProfile,
    Polymerase,
    melting_temperature,
)
from liulab_mbio.primers.thresholds import (
    TARGET_TM,
    THRESHOLDS,
    THRESHOLDS_FOR,
    Band,
    PrimerRole,
    Reading,
    Thresholds,
    reading,
)

__all__ = [
    "ONETAQ",
    "PHUSION",
    "POLYMERASES",
    "Q5",
    "TAQ",
    "TARGET_TM",
    "THRESHOLDS",
    "THRESHOLDS_FOR",
    "Amplicon",
    "Band",
    "GenomeReport",
    "Locus",
    "MadeBy",
    "PairReport",
    "PcrProfile",
    "Placement",
    "Polymerase",
    "PrimerReport",
    "PrimerRole",
    "PrimingSite",
    "Reading",
    "Thresholds",
    "amplicon_sizes",
    "annealing_checks",
    "design_pair",
    "design_primer",
    "evaluate_on_genome",
    "evaluate_pair",
    "evaluate_pair_on_genome",
    "evaluate_primer",
    "find_binding_sites",
    "find_priming_sites",
    "melting_temperature",
    "pair_checks",
    "reading",
]

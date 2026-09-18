"""Every oligo a restriction and ligation plan orders, and what it is for."""

from dataclasses import dataclass

from liulab_mbio.primers.evaluation import PrimerReport
from liulab_mbio.primers.thresholds import PrimerRole


@dataclass(frozen=True, slots=True)
class DesignedOligo:
    """One oligo a plan orders, and what it is for.

    Parameters
    ----------
    report
        What it scored, the primer included. An amplification primer carries a spacer and a
        recognition site as a 5' tail, which moves no annealing region: the melting temperature
        is of the annealing region and the structure checks are of the whole oligo.
    role
        What it is for: amplifying the insert, colony PCR, or sequencing the clone. It chooses
        the thresholds the oligo is designed and judged by, and the step its row points at.
    """

    report: PrimerReport
    role: PrimerRole

"""Every oligo a Gateway plan orders, and what it is for."""

from dataclasses import dataclass

from liulab_mbio.primers.evaluation import PrimerReport
from liulab_mbio.primers.thresholds import PrimerRole


@dataclass(frozen=True, slots=True)
class DesignedOligo:
    """One oligo a plan orders, and what it is for.

    Parameters
    ----------
    report
        What it scored, the primer included. An attB primer carries a whole att site as a 5'
        tail, which moves no annealing region: the melting temperature is of the annealing
        region and the structure checks are of the whole oligo.
    role
        What it is for: amplifying the insert onto attB ends, screening colonies, or sequencing
        the expression clone. It chooses the thresholds the oligo is designed and judged by, and
        the step its row points at.
    """

    report: PrimerReport
    role: PrimerRole

"""Every oligo a Golden Gate plan orders, and what it is for."""

from dataclasses import dataclass

from liulab_mbio.cloning.goldengate.assembly import Part
from liulab_mbio.primers.evaluation import PrimerReport
from liulab_mbio.primers.thresholds import PrimerRole


@dataclass(frozen=True, slots=True)
class DesignedOligo:
    """One oligo a plan orders, and what it is for.

    Parameters
    ----------
    report
        What it scored, the primer included.
    role
        What it is for: amplifying a part, colony PCR, or sequencing the clone. It chooses the
        thresholds the oligo is designed and judged by.
    part
        The part it amplifies, given for an amplification primer and for nothing else.

    Raises
    ------
    ValueError
        If a part is given for any role but amplification, or not given for amplification.
    """

    report: PrimerReport
    role: PrimerRole
    part: Part | None = None

    def __post_init__(self) -> None:
        """Refuse a part on anything but an amplification primer, and one missing from it."""
        if (self.role == "amplification") != (self.part is not None):
            raise ValueError(
                f"oligo {self.report.primer.name!r}: an amplification primer names its part, "
                "and no other role does"
            )

"""Reading a construct back one well at a time: how many colonies that takes, and what came back.

A design sits one per well and the well is sequenced and called on its own, so identity stays
with well position. The curve sizes the pick that fills those wells; the rest judges what they
read as. How a well is marked, and how deep its read has to be before anyone calls it, belong
to the method doing the marking and are not here.

The clean-colony curve is Lund et al. 2024's measurement, recorded in
``docs/research/synthesis-and-assembly.md``; the protocol that picks the colonies cites it.
"""

import itertools
from collections.abc import Sequence
from dataclasses import dataclass

from mbio.checks import Check, Status, worst
from mbio.protocol.model import Well

#: The fragment counts Lund et al. 2024 measured a clean colony at, and the share clean at each.
#: Six anchors and nothing between them: what an unmeasured count is worth is the caller's.
CLEAN_COLONY_CURVE: tuple[tuple[int, float], ...] = (
    (2, 1.0),
    (3, 0.938),
    (5, 0.846),
    (8, 0.667),
    (12, 0.40),
    (16, 0.0),
)


def identity_check(called: Sequence[str], designed: str) -> Check:
    """Return the verdict on whether a well's call is the design, base for base.

    A pass is an exact match across the whole designed region — both entry overhangs, the
    fragment, the stuffer and the barcode. A silent mismatch is not a pass: it leaves the
    protein right and the barcode wrong, and a barcode that no longer names its member cannot be
    put right by the linkage read. More than one consensus is mixed, and mixed fails. No
    consensus at all is nothing to judge, so it carries no verdict.

    Raises
    ------
    ValueError
        If `designed` is empty.

    Examples
    --------
    >>> identity_check(("ACGT",), "ACGT").status
    'pass'
    >>> identity_check(("ACGT", "ACGA"), "ACGT").status
    'fail'
    """
    if not designed:
        raise ValueError("a well is judged against a designed region, and this one is empty")
    consensus = [one.upper() for one in called]
    if not consensus:
        return Check("well_identity", None, 0.0, "no consensus was called from this well")
    if len(consensus) > 1:
        return Check(
            "well_identity",
            "fail",
            float(len(consensus)),
            f"{len(consensus)} consensus sequences: the well is mixed",
        )
    if consensus[0] != designed.upper():
        return Check("well_identity", "fail", 1.0, "the call is not the design, base for base")
    return Check("well_identity", "pass", 1.0, "the call is the design across its whole length")


@dataclass(frozen=True, slots=True)
class WellVerdict:
    """What one well's read came to, and whether reformatting carries it forward.

    Parameters
    ----------
    well
        Where the well is, as a plate name and a well name.
    checks
        The depth and the identity, in that order.
    """

    well: Well
    checks: tuple[Check, ...]

    @property
    def status(self) -> Status:
        """The worst verdict the well carries, by `mbio.checks.worst`."""
        return worst(one.status for one in self.checks)

    @property
    def called(self) -> bool:
        """Whether anything judged this well at all."""
        return any(one.status is not None for one in self.checks)


def reformat(verdicts: Sequence[WellVerdict]) -> tuple[WellVerdict, ...]:
    """Return the wells reformatting carries forward: everything that did not fail.

    A well nobody could call is kept. Compacting it out would throw away a design that may be
    clean and has only been read too thinly.

    Examples
    --------
    >>> reformat(())
    ()
    """
    return tuple(one for one in verdicts if one.status != "fail")


def clean_colony_chance(fragments: int) -> float:
    """Return the chance one picked colony carries a clean copy of a design of this many pieces.

    Linear interpolation between Lund's six measured anchors, `CLEAN_COLONY_CURVE`: a design in
    two pieces came out clean every time, one in sixteen never. Fewer pieces than the first
    anchor takes the first anchor's value and more than the last takes the last's, because
    nothing was measured outside them.

    Raises
    ------
    ValueError
        If `fragments` is not positive.

    Examples
    --------
    >>> clean_colony_chance(8), clean_colony_chance(16)
    (0.667, 0.0)
    """
    if fragments < 1:
        raise ValueError(f"a design is built from at least one fragment, got {fragments}")
    anchors = CLEAN_COLONY_CURVE
    if fragments <= anchors[0][0]:
        return anchors[0][1]
    for (low, left), (high, right) in itertools.pairwise(anchors):
        if fragments <= high:
            return left + (right - left) * (fragments - low) / (high - low)
    return anchors[-1][1]

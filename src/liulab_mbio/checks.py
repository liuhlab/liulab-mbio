"""A check: one verdict and the value it judged, and the worst of several verdicts.

A check no sourced threshold judges carries no verdict, ``None``, rather than a pass.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

#: One verdict a check can carry.
type Status = Literal["pass", "warn", "fail"]

#: The three verdicts, best first: the order `worst` ranks them in.
STATUSES: tuple[Status, ...] = ("pass", "warn", "fail")


@dataclass(frozen=True, slots=True)
class Check:
    """One verdict on a primer, a pair or an assembled product, with the value it judged.

    Parameters
    ----------
    name
        What was measured, such as ``"gc_clamp"``.
    status
        One of `STATUSES`, or ``None`` where no sourced threshold judges it.
    value
        The measurement, in the unit the module that judged it documents.
    detail
        What a reader needs besides the number.
    """

    name: str
    status: Status | None
    value: float
    detail: str = ""


def worst(statuses: Iterable[Status | None]) -> Status:
    """Return the worst of these verdicts, passing over any that is ``None``.

    Nothing judged at all is a pass.

    Examples
    --------
    >>> worst(("pass", None, "warn"))
    'warn'
    """
    result: Status = "pass"
    for status in statuses:
        if status is not None and STATUSES.index(status) > STATUSES.index(result):
            result = status
    return result

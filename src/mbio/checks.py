"""A check: one verdict and the value it judged, and the worst of several verdicts.

A check no sourced threshold judges carries no verdict, ``None``, rather than a pass.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal


def counted(count: int, noun: str, plural: str = "") -> str:
    """Return a count and its noun, the noun plural unless the count is one.

    Every field a reader follows says ``1 tube`` and ``3 tubes``; none prints ``tube(s)``.
    `plural` carries an ending the noun does not take an ``s`` for.

    Examples
    --------
    >>> counted(1, "tube"), counted(3, "tube"), counted(2, "colony", "colonies")
    ('1 tube', '3 tubes', '2 colonies')
    """
    return f"{count:,} {noun if count == 1 else plural or noun + 's'}"


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


def worst_of(checks: Iterable[Check]) -> Status:
    """Return the worst verdict these checks carry, by the rule `worst` states."""
    return worst(check.status for check in checks)

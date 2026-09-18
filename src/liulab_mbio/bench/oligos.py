"""The primer order sheet: every oligo a design asks for, as a table to order from.

`primer_sheet` writes it as a file, and `oligo_row` makes one row of a protocol's own sheet.
"""

from collections.abc import Sequence

from liulab_mbio.bench.pcr import PRIMER_STOCK_UM
from liulab_mbio.primers.evaluation import PrimerReport
from liulab_mbio.primers.thresholds import Thresholds, reading
from liulab_mbio.protocol.model import Check, Oligo

#: The columns of the primer order sheet.
SHEET_COLUMNS = ("name", "sequence", "length", "tm_c")


def primer_sheet(reports: Sequence[PrimerReport]) -> str:
    """Return these oligos as a tab-separated sheet, one row each, in the order given.

    The columns are `SHEET_COLUMNS`: the name to order it under, the sequence 5' to 3', its
    length, and the Tm of the part that anneals.
    """
    rows = ["\t".join(SHEET_COLUMNS)]
    for report in reports:
        primer = report.primer
        rows.append(
            "\t".join(
                (
                    primer.name,
                    primer.sequence,
                    str(len(primer.sequence)),
                    f"{report['tm'].value:.1f}",
                )
            )
        )
    return "\n".join(rows) + "\n"


def oligo_row(report: PrimerReport, *, purpose: str, thresholds: Thresholds) -> Oligo:
    """Return one oligo as a row of a protocol's order sheet, carrying its verdict.

    Parameters
    ----------
    report
        What the oligo scored, the primer included.
    purpose
        The title of the step that uses it.
    thresholds
        What it was judged by, so a row prints the band beside each value that missed it.
    """
    return Oligo(
        report.primer.name,
        report.primer.sequence,
        purpose=purpose,
        tm_c=round(report["tm"].value, 1),
        stock=f"{PRIMER_STOCK_UM:g} µM",
        status=report.status,
        checks=_fired(report, thresholds),
    )


def _fired(report: PrimerReport, thresholds: Thresholds) -> tuple[Check, ...]:
    """Return the checks that did not pass, each as its value and the band it missed."""
    fired = []
    for check in report.checks:
        status = check.status
        if status is None or status == "pass":
            continue
        word = reading(check, thresholds)
        fired.append(Check(word.label, status, detail=word.detail))
    return tuple(fired)

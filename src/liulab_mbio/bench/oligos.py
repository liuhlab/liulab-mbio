"""The primer order sheet: every oligo a design asks for, as a table to order from.

`primer_sheet` writes it as a file, and `oligo_row` and `ordered_row` make one row of a
protocol's own sheet. Not every oligo a design asks for is a primer: one that anneals to another
fragment rather than to a template primes nothing, so no threshold in
`liulab_mbio.primers.thresholds` judges it. `OrderedOligo` is that oligo, and its row carries no
verdict rather than a pass nothing measured.
"""

from collections.abc import Sequence
from dataclasses import dataclass

from liulab_mbio.bench.pcr import PRIMER_STOCK_UM
from liulab_mbio.primers.evaluation import PrimerReport
from liulab_mbio.primers.thresholds import Thresholds, reading
from liulab_mbio.protocol.model import Check, Oligo

#: The columns of the primer order sheet.
SHEET_COLUMNS = ("name", "sequence", "length", "tm_c")


@dataclass(frozen=True, slots=True)
class OrderedOligo:
    """One oligo to order that primes nothing, so no sourced threshold judges it.

    Parameters
    ----------
    name, sequence
        What to order it under, and its bases 5' to 3'.
    purpose
        What it is for.
    stock
        The concentration its source asks for it at.
    note
        Its purity, and why the row carries no verdict.
    """

    name: str
    sequence: str
    purpose: str = ""
    stock: str = ""
    note: str = ""


def primer_sheet(reports: Sequence[PrimerReport], *, oligos: Sequence[OrderedOligo] = ()) -> str:
    """Return these oligos as a tab-separated sheet, one row each, in the order given.

    The columns are `SHEET_COLUMNS`: the name to order it under, the sequence 5' to 3', its
    length, and the Tm of the part that anneals. `oligos` are the ones that prime nothing, and
    they follow the primers with their melting temperature column left empty.
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
    rows.extend(
        "\t".join((oligo.name, oligo.sequence, str(len(oligo.sequence)), "")) for oligo in oligos
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


def ordered_row(oligo: OrderedOligo) -> Oligo:
    """Return one oligo that primes nothing as a row of a protocol's order sheet.

    The row carries no verdict and says in its note why, so a silence is never read as a pass.
    """
    return Oligo(
        oligo.name,
        oligo.sequence,
        purpose=oligo.purpose,
        stock=oligo.stock,
        note=oligo.note,
        status=None,
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

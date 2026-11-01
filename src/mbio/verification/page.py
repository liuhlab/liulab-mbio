"""The page sequence verification writes for one clone: its verdict, the map, each disagreement.

One HTML file that loads nothing over the network: every drawing is molviz's still SVG, with
molviz's fonts written once. It reads top to bottom:

1. whether the clone is verified;
2. the table the command prints: a row per region, one per disagreement outside every region
   naming the features it falls in, with no verdict, and one per result;
3. the map: the record with each result a read over the bases it trusts, on its own strand,
   and each disagreement coloured by its region's verdict, each result and each disagreement
   lit, as `molviz_verification` hands them over;
4. a close-up of each disagreement: a sequence view a few dozen bases either side, each result's
   own disagreements marked on its read, and each Sanger result's trace under the bases.

A result that does not read as the record has its row and nothing drawn. When no result reads
as it, the page holds only those rows and the verdict.
"""

import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from html import escape
from importlib.resources import files
from pathlib import Path

import molviz

from mbio.checks import Status, counted
from mbio.plot.converter import molviz_record
from mbio.sequence import SequenceRecord
from mbio.verification.converter import molviz_verification
from mbio.verification.judge import CONSENSUS_CAVEAT, Disagreement, Placement, Verification
from mbio.verification.result import SequencingResult
from mbio.verification.trace import Channels

#: How many bases a close-up draws either side of its disagreement.
FLANK = 30

#: The most bases an insertion's caption spells out; a longer one is named by how many it adds.
_SPELLED = 6

#: What a verdict reads as on the page; a check no threshold judges says so.
_WORDS: Mapping[Status | None, str] = {
    "pass": "pass",
    "warn": "warn",
    "fail": "fail",
    None: "no verdict",
}


@dataclass(frozen=True, slots=True)
class TableRow:
    """One row of the page's table.

    Parameters
    ----------
    where
        A region's name, the features a disagreement outside every region falls in, or a
        result's name.
    verdict
        The row's verdict; ``None`` where nothing judges it.
    detail
        What was found there.
    """

    where: str
    verdict: Status | None
    detail: str


@dataclass(frozen=True, slots=True, eq=False)
class CloseUp:
    """One disagreement drawn close: its name, what it is, and the drawing."""

    disagreement: Disagreement
    name: str
    caption: str
    drawing: molviz.Drawing


@dataclass(frozen=True, slots=True, eq=False)
class Page:
    """One clone's verification laid out, to write as one page.

    Parameters
    ----------
    title
        What the page is called.
    verified
        Whether the clone is verified.
    rows
        The table, in the page's order.
    map
        molviz's drawing of the record with its reads and disagreements; ``None`` where no
        result reads as it.
    close_ups
        One per disagreement, in record order.
    """

    title: str
    verified: bool
    rows: tuple[TableRow, ...]
    map: molviz.Drawing | None
    close_ups: tuple[CloseUp, ...]

    @property
    def html(self) -> str:
        """The page, one self-contained HTML file."""
        verdict = "Verified" if self.verified else "Not verified"
        state = "pass" if self.verified else "fail"
        body = [
            f'<header><h1 class="verdict is-{state}">{verdict}</h1>'
            f"<p>{escape(self.title)}</p></header>",
            _table(self.rows),
        ]
        if self.map is not None:
            body.append(
                '<section class="map"><h2>Map</h2>'
                f"<figure>{self.map.still_svg()}</figure></section>"
            )
        if self.close_ups:
            figures = "".join(
                f'<figure class="close-up"><figcaption><strong>{escape(one.name)}</strong> '
                f"{escape(one.caption)}</figcaption>{one.drawing.still_svg()}</figure>"
                for one in self.close_ups
            )
            body.append(f'<section class="close-ups"><h2>Each disagreement</h2>{figures}</section>')
        fonts = f"{molviz.font_rules()}\n" if self.map is not None else ""
        return (
            '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
            '<meta name="color-scheme" content="light">\n'
            f"<title>{escape(self.title)}</title>\n<style>\n{fonts}{_asset('page.css')}</style>\n"
            f'</head>\n<body>\n<main class="verification">\n{"\n".join(body)}\n</main>\n'
            "</body>\n</html>\n"
        )

    def write(self, path: str | os.PathLike[str]) -> Path:
        """Write the page to `path`, making its directory if it is not there, and return it."""
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(self.html, encoding="utf-8")
        return out


def draw_page(
    expected: SequenceRecord,
    results: Sequence[SequencingResult],
    verification: Verification,
    *,
    channels: Sequence[Channels | None] = (),
) -> Page:
    """Lay out the page for one clone's `verification` of `results` against `expected`.

    Parameters
    ----------
    expected
        The record the results were held against.
    results
        The results `verify` was given.
    verification
        What `verify` made of them.
    channels
        Each result's channels and peaks, in the order of `results`, drawn under its bases in
        every close-up it reaches: ``None`` for a result that is no Sanger trace, and empty
        where none is drawn. They pair with the results by position, as the checks and the
        placements do, so two results of one name each draw their own.

    Raises
    ------
    ValueError
        If `channels` is given, but not one for each result.
    """
    added = molviz_verification(expected, results, verification, channels=channels)
    statuses: dict[str, Status | None] = {one.name: one.status for one in verification.checks}
    # Withheld from a result that does not read as the record, and from one that trusts nothing.
    placed = [one for one in verification.placements if one.strand is not None]
    results_rows = tuple(
        TableRow(check.name, check.status, _said(check.detail, result, placement))
        for check, result, placement in zip(
            verification.result_checks, results, verification.placements, strict=True
        )
    )
    title = f"{expected.name or 'record'}: {counted(len(results), 'result')}"
    if not placed:
        return Page(title, verification.verified, results_rows, None, ())
    n = len(expected)
    record = molviz_record(expected) | added
    named = [
        (given["name"], one)
        for given, one in zip(
            record.get("disagreements", []), verification.disagreements, strict=True
        )
    ]
    rows = (
        *(TableRow(check.name, check.status, check.detail) for check in verification.checks),
        *(
            TableRow(_outside(one), None, f"{one.said(n)}, read by {', '.join(one.results)}")
            for _, one in named
            if not one.regions
        ),
        *results_rows,
    )
    reads = [read["name"] for read in record.get("reads", [])]
    mapped = molviz.draw(
        record, kinds_off=["cuts"], highlight=[*(name for name, _ in named), *reads]
    )
    close_ups = []
    for name, one in named:
        start, end = _window(one, expected)
        close = molviz.draw(
            record,
            span={"start": start, "end": end},
            map="line",
            sequence_view=True,
            bases_per_row=end - start,
            kinds_off=["cuts"],
            highlight=[name, *reads],
        )
        close_ups.append(CloseUp(one, name, _caption(one, statuses), close))
    return Page(title, verification.verified, rows, mapped, tuple(close_ups))


def _said(detail: str, result: SequencingResult, placement: Placement) -> str:
    """Return a result's row as the command says it: a placed consensus says what it cannot show."""
    if placement.strand is not None and result.quality is None and result.depth is None:
        return f"{detail}; {CONSENSUS_CAVEAT}"
    return detail


def _outside(one: Disagreement) -> str:
    """Return what a disagreement outside every region falls in."""
    if one.features:
        return f"{', '.join(one.features)}, outside every region"
    return "outside every region and feature"


def _caption(one: Disagreement, statuses: Mapping[str, Status | None]) -> str:
    """Return what a close-up says of its disagreement: what it is, where, and who read it."""
    if one.regions:
        where = "in " + ", ".join(
            f"{region} ({_WORDS[statuses.get(region)]})" for region in one.regions
        )
    elif one.features:
        where = f"in {', '.join(one.features)}, outside every region"
    else:
        where = "outside every region and feature"
    kind = one.kind
    if kind == "insertion":
        added = one.bases if len(one.bases) <= _SPELLED else counted(len(one.bases), "base")
        kind = f"insertion of {added}"
    return f"{kind} {where}, read by {', '.join(one.results)}."


def _window(one: Disagreement, record: SequenceRecord) -> tuple[int, int]:
    """Return the stretch a close-up draws: `FLANK` bases either side, within one turn."""
    n = len(record)
    start, end = one.start - FLANK, one.end + FLANK
    if record.topology == "linear":
        return max(0, start), min(n, end)
    if end - start >= n:
        return 0, n
    if start < 0:
        start, end = start + n, end + n
    return start, end


def _table(rows: Sequence[TableRow]) -> str:
    body = "".join(
        f'<tr class="is-{row.verdict or "none"}"><th scope="row">{escape(row.where)}</th>'
        f'<td class="verdict">{_WORDS[row.verdict]}</td><td>{escape(row.detail)}</td></tr>'
        for row in rows
    )
    return (
        '<table class="checks"><thead><tr><th scope="col">Where</th>'
        '<th scope="col">Verdict</th><th scope="col">Detail</th></tr></thead>'
        f"<tbody>{body}</tbody></table>"
    )


@cache
def _asset(name: str) -> str:
    return files("mbio.verification").joinpath(name).read_text(encoding="utf-8")

"""The page sequence verification writes for one clone: its verdict, the map, each disagreement.

One HTML file that loads nothing over the network, written as `mbio.plot.page` writes a map's:
every drawing is the SVG `mbio.plot` lays out, with the fonts it was measured in embedded. It
reads top to bottom:

1. whether the clone is verified;
2. the table the command prints: a row per region, one per disagreement outside every region
   naming the features it falls in, with no verdict, and one per result;
3. the map: a copy of the record with each result a primer over the bases it trusts, on its own
   strand, and each disagreement a feature coloured by its region's verdict, every one lit;
4. a close-up of each disagreement: a sequence view a few dozen bases either side, each result's
   own disagreements marked on its arrow, and each Sanger result's trace under the bases.

A result that does not read as the record has its row and nothing drawn. When no result reads
as it, the page holds only those rows and the verdict.
"""

import dataclasses
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from html import escape
from importlib.resources import files
from pathlib import Path

from mbio.checks import Status, counted, worst
from mbio.plot import Drawing, draw_map, layers
from mbio.plot import sequence_view as view
from mbio.plot.fonts import BOLD, MONO, SANS
from mbio.plot.page import font_face
from mbio.sequence import (
    BindingSite,
    Feature,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    position_text,
    reverse_complement,
    span_text,
)
from mbio.verification.align import place
from mbio.verification.judge import CONSENSUS_CAVEAT, Disagreement, Placement, Verification
from mbio.verification.result import SequencingResult
from mbio.verification.trace import Channels

#: How many bases a close-up draws either side of its disagreement.
FLANK = 30

#: What a disagreement is coloured by its region's verdict; one outside every region by `None`.
_VERDICT_COLORS: Mapping[Status | None, str] = {
    "pass": "#009988",
    "warn": "#ee7733",
    "fail": "#cc3311",
    None: "#bbbbbb",
}

#: Each base's curve under a close-up, in the colours a Sanger trace is read in.
_BASE_COLORS: Mapping[str, str] = {
    "A": "#228833",
    "C": "#4477aa",
    "G": "#000000",
    "T": "#cc3311",
}

#: The most bases an insertion's name spells out; a longer one is named by how many it adds.
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
    drawing: Drawing


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
        The record with its reads and disagreements; ``None`` where no result reads as it.
    close_ups
        One per disagreement, in record order.
    """

    title: str
    verified: bool
    rows: tuple[TableRow, ...]
    map: Drawing | None
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
                f'<section class="map"><h2>Map</h2><figure>{self.map.element()}</figure></section>'
            )
        if self.close_ups:
            traced = any(
                row.strips
                for one in self.close_ups
                if one.drawing.sequence_view is not None
                for row in one.drawing.sequence_view.rows
            )
            legend = (
                '<p class="legend">Each trace: '
                + ", ".join(
                    f'<span style="color: {color}">{base}</span>'
                    for base, color in _BASE_COLORS.items()
                )
                + ", over the top strand's bases.</p>"
                if traced
                else ""
            )
            figures = "".join(
                f'<figure class="close-up"><figcaption><strong>{escape(one.name)}</strong> '
                f"{escape(one.caption)}</figcaption>{one.drawing.element()}</figure>"
                for one in self.close_ups
            )
            body.append(
                f'<section class="close-ups"><h2>Each disagreement</h2>{legend}{figures}</section>'
            )
        fonts = "".join(font_face(font) for font in (SANS, BOLD, MONO))
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
    if channels and len(channels) != len(results):
        raise ValueError(
            f"{counted(len(channels), 'set')} of channels were given for "
            f"{counted(len(results), 'result')}; give one for each, None for one with no trace"
        )
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
    named = [(_name(one, expected), one) for one in verification.disagreements]
    rows = (
        *(TableRow(check.name, check.status, check.detail) for check in verification.checks),
        *(
            TableRow(_outside(one), None, f"{one.said(n)}, read by {', '.join(one.results)}")
            for _, one in named
            if not one.regions
        ),
        *results_rows,
    )
    colored = {name: _VERDICT_COLORS[_verdict(one, statuses)] for name, one in named}
    features = tuple(
        Feature(name, "misc_difference", (Segment(one.start, one.end),), color=colored[name])
        for name, one in named
        if one.kind != "insertion"
    )
    insertions = tuple(
        layers.Insertion(name, one.start, len(one.bases), colored[name])
        for name, one in named
        if one.kind == "insertion"
    )
    primers = tuple(_primer(expected, one, verification.disagreements) for one in placed)
    copy = dataclasses.replace(
        expected,
        features=(*expected.features, *features),
        primers=(*expected.primers, *primers),
    )
    lit = [name for name, _ in named]
    drawn = draw_map(copy, cut_sites=False, insertions=insertions, highlight=lit)
    tracks = _tracks(expected, results, verification, channels or [None] * len(results))
    close_ups = []
    for name, one in named:
        start, end = _window(one, expected)
        close = draw_map(
            copy,
            region=(start, end),
            sequence_view=True,
            bases_per_row=end - start,
            cut_sites=False,
            insertions=insertions,
            tracks=tracks,
            highlight=[name, *(primer.name for primer in primers)],
        )
        close_ups.append(CloseUp(one, name, _caption(one, statuses), close))
    return Page(title, verification.verified, rows, drawn, tuple(close_ups))


def _said(detail: str, result: SequencingResult, placement: Placement) -> str:
    """Return a result's row as the command says it: a placed consensus says what it cannot show."""
    if placement.strand is not None and result.quality is None and result.depth is None:
        return f"{detail}; {CONSENSUS_CAVEAT}"
    return detail


def _verdict(one: Disagreement, statuses: Mapping[str, Status | None]) -> Status | None:
    """Return the verdict a disagreement takes its colour from: the worst of its regions'."""
    found: list[Status | None] = [statuses[name] for name in one.regions if name in statuses]
    if not found or all(status is None for status in found):
        return None
    return worst(found)


def _name(one: Disagreement, record: SequenceRecord) -> str:
    """Return a disagreement's name as a reader reads it, 1-based: ``1234 C>T``."""
    n = len(record)
    if one.kind == "insertion":
        added = one.bases if len(one.bases) <= _SPELLED else counted(len(one.bases), "base")
        point = layers.point_text(one.start, n, circular=record.topology == "circular")
        return f"{point} ins {added}"
    if one.kind == "deletion":
        if one.end - one.start == 1:
            return f"{position_text(one.start, n)} del"
        return f"{span_text(one.start, one.end, n)} del"
    expected = record.bases(one.start, one.end)
    return f"{position_text(one.start, n)} {expected}>{one.bases}"


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
    return f"{one.kind} {where}, read by {', '.join(one.results)}."


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


def _primer(
    record: SequenceRecord, placement: Placement, disagreements: Sequence[Disagreement]
) -> Primer:
    """Return a placed result as a primer over its trusted span, its bases those it reads.

    Each disagreement it shows is written into the record's bases there, a deleted base as
    ``N``, so the sequence view marks it on this result's arrow and on no other.
    """
    span, strand = placement.span, placement.strand
    assert span is not None
    assert strand is not None
    n = len(record)
    bases = list(record.bases(span.start, span.end))
    for one in disagreements:
        if placement.result not in one.results or one.kind == "insertion":
            continue
        for offset, position in enumerate(range(one.start, one.end)):
            # Counted from the read's start round a circular record, as its span is.
            at = (position - span.start) % n
            if at < len(bases):
                bases[at] = one.bases[offset] if one.bases else "N"
    sequence = "".join(bases)
    if strand is Strand.REVERSE:
        sequence = reverse_complement(sequence)
    site = BindingSite(span.start, span.end, strand)
    return Primer(placement.result, sequence, binding_sites=(site,))


def _tracks(
    record: SequenceRecord,
    results: Sequence[SequencingResult],
    verification: Verification,
    channels: Sequence[Channels | None],
) -> tuple[view.Track, ...]:
    """Return a track of each placed Sanger result's trace, each peak over the base it reads.

    Only a base laid over one base places a peak; a result read along the bottom strand draws
    each channel as the complement it reads on the top one.
    """
    tracks = []
    n = len(record)
    for result, placement, trace in zip(results, verification.placements, channels, strict=True):
        if trace is None or placement.strand is None:
            continue
        laid = place(record, result)
        if laid is None:
            continue
        peaks = sorted(
            (trace.peaks[column.behind[0]], column.start % n)
            for column in laid.placed
            if column.end - column.start == 1 and len(column.behind) == 1
        )
        reverse = laid.strand is Strand.REVERSE
        curves = tuple(
            view.Curve(
                tuple(trace.scans[reverse_complement(base) if reverse else base]),
                color,
            )
            for base, color in _BASE_COLORS.items()
        )
        tracks.append(view.Track(result.name, curves, tuple(peaks)))
    return tuple(tracks)


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

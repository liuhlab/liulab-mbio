"""A clone's verification in molviz's record format: its reads, disagreements and traces.

`molviz_verification` returns what the verification adds to its record, to merge into
`mbio.plot.converter.molviz_record` of the same record. Each placed result is a read over the
bases it trusts; one that does not read as the record is left out. molviz refuses two reads of
one name, so a later result of a name already taken is numbered, ``name (2)``, and a
disagreement names the read of each result that shows it. A disagreement is named as a reader
reads it, ``696 A>G``, and coloured by the worst verdict of the regions it falls in, grey outside
every region. `molviz_names` gives every name, for a page to light and caption what is drawn. A
Sanger result's trace is a track: each peak over the one base it was read as, and on a
bottom-strand read each channel is the base it pairs with on the top strand.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from mbio.checks import Status, counted, worst
from mbio.sequence import SequenceRecord, Strand, position_text, reverse_complement, span_text
from mbio.verification.align import place
from mbio.verification.judge import Disagreement, Placement, Verification
from mbio.verification.result import SequencingResult
from mbio.verification.trace import Channels

#: What a disagreement is coloured by its regions' worst verdict; one outside every region by
#: `None`.
_VERDICT_COLORS: Mapping[Status | None, str] = {
    "pass": "#009988",
    "warn": "#ee7733",
    "fail": "#cc3311",
    None: "#bbbbbb",
}


@dataclass(frozen=True, slots=True)
class Names:
    """What a verification's drawing calls each result and each disagreement.

    Parameters
    ----------
    results
        Each result's read, in the order the results were given.
    disagreements
        Each disagreement's name, in record order.
    """

    results: tuple[str, ...]
    disagreements: tuple[str, ...]

    def readers(self, one: Disagreement) -> list[str]:
        """Return the read of each result that shows `one`."""
        return [self.results[at] for at in one.shown_by]


def molviz_names(expected: SequenceRecord, verification: Verification) -> Names:
    """Return what `molviz_verification` calls each result and each disagreement.

    A result is its name, numbered from 2 where an earlier result has taken it; a disagreement
    is named as a reader reads it, 1-based: ``1234 C>T``, ``1234 del``, or, as the command
    names an insertion, ``insertion after 1234``.
    """
    return Names(
        tuple(_unique(verification.placements)),
        tuple(_name(one, expected) for one in verification.disagreements),
    )


def molviz_verification(
    expected: SequenceRecord,
    results: Sequence[SequencingResult],
    verification: Verification,
    *,
    channels: Sequence[Channels | None] = (),
) -> dict[str, Any]:
    """Return the reads, disagreements and tracks `verification` adds to `expected` for molviz.

    Each key is left out where its list would be empty.

    Parameters
    ----------
    expected
        The record the results were held against.
    results
        The results `verify` was given.
    verification
        What `verify` made of them.
    channels
        Each result's channels, in the order of `results`: ``None`` for a result that is no
        Sanger trace, and empty where no trace is drawn.

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
    named = molviz_names(expected, verification)
    reads: list[dict[str, Any]] = []
    tracks: list[dict[str, Any]] = []
    for name, placement, result, trace in zip(
        named.results,
        verification.placements,
        results,
        channels or [None] * len(results),
        strict=True,
    ):
        span, strand = placement.span, placement.strand
        if span is None or strand is None:
            continue
        reads.append({"name": name, "start": span.start, "end": span.end, "strand": int(strand)})
        if trace is not None and (track := _track(name, expected, result, trace)):
            tracks.append(track)
    verdicts = verification.verdicts
    disagreements = [
        _disagreement(one, name, named.readers(one), expected, verdicts)
        for name, one in zip(named.disagreements, verification.disagreements, strict=True)
    ]
    out: dict[str, Any] = {}
    for key, values in (("reads", reads), ("disagreements", disagreements), ("tracks", tracks)):
        if values:
            out[key] = values
    return out


def _name(one: Disagreement, record: SequenceRecord) -> str:
    n = len(record)
    if one.kind == "insertion":
        return one.said(n)
    if one.kind == "deletion":
        if one.end - one.start == 1:
            return f"{position_text(one.start, n)} del"
        return f"{span_text(one.start, one.end, n)} del"
    return f"{position_text(one.start, n)} {record.bases(one.start, one.end)}>{one.bases}"


def _unique(placements: Sequence[Placement]) -> list[str]:
    """Return each placement's result name, a name already taken numbered from 2."""
    taken: set[str] = set()
    out = []
    for one in placements:
        name, n = one.result, 1
        while name in taken:
            n += 1
            name = f"{one.result} ({n})"
        taken.add(name)
        out.append(name)
    return out


def _disagreement(
    one: Disagreement,
    name: str,
    reads: Sequence[str],
    record: SequenceRecord,
    verdicts: Mapping[str, Status | None],
) -> dict[str, Any]:
    where = (
        {"at": one.start % len(record) if record.topology == "circular" else one.start}
        if one.kind == "insertion"
        else {"start": one.start, "end": one.end}
    )
    return {
        "name": name,
        **where,
        "bases": one.bases,
        "reads": list(reads),
        "color": _VERDICT_COLORS[_verdict(one, verdicts)],
    }


def _verdict(one: Disagreement, verdicts: Mapping[str, Status | None]) -> Status | None:
    """Return the verdict a disagreement takes its colour from: the worst of its regions'."""
    found: list[Status | None] = [verdicts[name] for name in one.regions if name in verdicts]
    if not found or all(status is None for status in found):
        return None
    return worst(found)


def _track(
    name: str, record: SequenceRecord, result: SequencingResult, trace: Channels
) -> dict[str, Any] | None:
    """Return a trace as a track, each peak over the one base its called base lies over.

    ``None`` where the result does not lie on the record.
    """
    laid = place(record, result)
    if laid is None:
        return None
    n = len(record)
    peaks = sorted(
        (trace.peaks[column.behind[0]], column.start % n)
        for column in laid.placed
        if column.end - column.start == 1 and len(column.behind) == 1
    )
    reverse = laid.strand is Strand.REVERSE
    return {
        "name": name,
        "curves": {
            base: list(trace.scans[reverse_complement(base) if reverse else base])
            for base in "ACGT"
        },
        "peaks": [list(peak) for peak in peaks],
    }

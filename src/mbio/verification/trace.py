"""Read a Sanger trace, an ``.ab1`` file, into a sequencing result, and into its four channels.

Only what the basecaller wrote is read, never a person's edits: ``PBAS2`` the bases, ``PCON2``
their Phred qualities, ``PLOC2`` the scan at each base's peak, and the analysed channels
``DATA9`` to ``DATA12``, which hold the four bases in ``FWO_1``'s order. `read_trace` reads what
a verification judges, and `read_channels` what a page draws under the bases, so a sequencing
result carries nothing only a trace has.
"""

import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from mbio.sequence import IUPAC_BASES
from mbio.verification.result import SequencingResult

#: A second peak this share of the called base's peak, or more, makes the base mixed: the top of
#: Applied Biosystems' recommended 15 to 25%, so a clean trace reads the fewest
#: (``docs/research/sequencing-read-evidence.md`` section 2.5).
MIXED_SHARE = 0.25

#: Mott's cutoff and the shortest read worth trimming, Biopython's ``abi-trim`` defaults
#: (section 2.4).
_MOTT_CUTOFF = 0.05
_MOTT_SHORTEST = 20

_CHANNELS = ("DATA9", "DATA10", "DATA11", "DATA12")
_PAIR_CODE = {bases: code for code, bases in IUPAC_BASES.items() if len(bases) == 2}


def read_trace(path: str | os.PathLike[str]) -> SequencingResult:
    """Read an ``.ab1`` trace into a sequencing result.

    Parameters
    ----------
    path
        An ``.ab1`` file a basecaller wrote.

    Returns
    -------
    SequencingResult
        Named for the file. ``bases`` are the bases the basecaller called, ``quality`` their
        Phred values, and ``trusted`` the span Biopython's ``abi-trim`` keeps by Mott's rule. A
        trace whose qualities are all zero trusts nothing.

        A base called A, C, G or T is mixed, and written as the IUPAC code of its two bases,
        where another channel's signal at its peak reaches 25% of the called base's. That signal
        is counted from the higher of its lowest points between this peak and each neighbouring
        base's, so the tail of a neighbouring peak in the same channel makes no base mixed. Any
        other code the basecaller wrote is kept.

    Raises
    ------
    ValueError
        If the file lacks the bases, qualities, peak positions or analysed channels.
    """
    raw = _tags(path)
    called = raw["PBAS2"].decode("ascii").upper()
    quality = tuple(raw["PCON2"])
    channels = dict(zip(raw["FWO_1"].decode("ascii"), (raw[tag] for tag in _CHANNELS), strict=True))
    bases = "".join(
        _read_base(base, raw["PLOC2"], index, channels) if base in channels else base
        for index, base in enumerate(called)
    )
    trusted = _mott(quality) if any(quality) else (0, 0)
    return SequencingResult(Path(path).stem, bases, quality=quality, trusted=trusted)


@dataclass(frozen=True, slots=True)
class Channels:
    """What a trace's bases were called from, for a page to draw under them.

    Parameters
    ----------
    scans
        Each base's analysed channel, one value a scan, by the base: A, C, G and T.
    peaks
        The scan at each base's peak, one for each base `read_trace` reads from the same file.
    """

    scans: Mapping[str, tuple[int, ...]]
    peaks: tuple[int, ...]


def read_channels(path: str | os.PathLike[str]) -> Channels:
    """Read an ``.ab1`` trace's four channels and the scan at each base's peak.

    A page takes it keyed by the name of the result `read_trace` reads from the same file.

    Raises
    ------
    ValueError
        As `read_trace` raises.
    """
    raw = _tags(path)
    scans = {
        base: tuple(raw[tag])
        for base, tag in zip(raw["FWO_1"].decode("ascii"), _CHANNELS, strict=True)
    }
    return Channels(scans, tuple(raw["PLOC2"]))


def _tags(path: str | os.PathLike[str]) -> Mapping[str, Any]:
    """Return a trace's tags, refusing one that lacks a tag either reader takes."""
    from Bio import SeqIO

    raw = SeqIO.read(path, "abi").annotations["abif_raw"]
    if missing := [
        tag for tag in ("PBAS2", "PCON2", "PLOC2", "FWO_1", *_CHANNELS) if tag not in raw
    ]:
        raise ValueError(f"{path} is not a base-called trace: it holds no {', '.join(missing)}")
    return raw


def _read_base(
    base: str, peaks: Sequence[int], index: int, channels: dict[str, Sequence[int]]
) -> str:
    """Return the base at `index`, or the IUPAC code of it and a second peak reaching its share."""
    scan = peaks[index]
    before = peaks[index - 1] if index else 0
    after = peaks[index + 1] if index + 1 < len(peaks) else len(channels[base]) - 1
    height = channels[base][scan]
    if height <= 0:
        return base
    rising = {
        other: signal[scan] - max(min(signal[before : scan + 1]), min(signal[scan : after + 1]))
        for other, signal in channels.items()
        if other != base
    }
    rival = max(rising, key=rising.__getitem__)
    return _PAIR_CODE[frozenset(base + rival)] if rising[rival] >= MIXED_SHARE * height else base


def _mott(quality: Sequence[int]) -> tuple[int, int]:
    """Return the half-open span Mott's rule keeps, as Biopython's ``abi-trim`` computes it."""
    if len(quality) <= _MOTT_SHORTEST:
        return 0, len(quality)
    total, best, start, end = 0.0, 0.0, None, 0
    for index in range(1, len(quality)):
        total += _MOTT_CUTOFF - 10 ** (-quality[index] / 10)
        if total < 0:
            total = 0.0
        elif start is None:
            start = index
        if total > best:
            best, end = total, index
    return (start or 0), end

"""Read a Sanger trace, an ``.ab1`` file, into a sequencing result.

Only the basecaller's own tags are read, never a person's edits: ``PBAS2`` the calls, ``PCON2``
their Phred qualities, ``PLOC2`` the scan at each call's peak, and the analysed channels
``DATA9`` to ``DATA12``, which hold the bases in ``FWO_1``'s order.
"""

from collections.abc import Sequence
from pathlib import Path

from mbio.sequence import IUPAC_BASES
from mbio.verification.result import SequencingResult

#: A second peak this share of the called base's peak, or more, makes the call mixed: the top of
#: Applied Biosystems' recommended 15 to 25%, so a clean trace calls the fewest
#: (``docs/research/sequencing-read-evidence.md`` section 2.5).
MIXED_SHARE = 0.25

_PAIR_CODE = {bases: code for code, bases in IUPAC_BASES.items() if len(bases) == 2}


def read_trace(path: str | Path) -> SequencingResult:
    """Read an ``.ab1`` trace into the calls it trusts.

    Parameters
    ----------
    path
        An ``.ab1`` file a basecaller wrote.

    Returns
    -------
    SequencingResult
        Named for the file. ``bases`` are the basecaller's calls, ``quality`` their Phred values,
        and ``trusted`` the span Biopython's ``abi-trim`` keeps: each base scores 0.05 less its
        chance of being wrong, and the span is the run whose total is greatest. A trace whose
        qualities are all zero trusts nothing.

        A call of A, C, G or T is mixed, and written as the IUPAC code of its two bases, where
        another channel's signal at the call's peak reaches 25% of the called base's. That signal
        is counted from the higher of its lowest points between this call and each neighbouring
        call, so the tail of a neighbouring peak in the same channel makes no call mixed. Any
        other code the basecaller wrote is kept.

    Raises
    ------
    ValueError
        If the file lacks the calls, qualities, peak positions or analysed channels.
    """
    from Bio import SeqIO

    raw = SeqIO.read(path, "abi").annotations["abif_raw"]
    missing = [
        tag
        for tag in ("PBAS2", "PCON2", "PLOC2", "FWO_1", "DATA9", "DATA10", "DATA11", "DATA12")
        if tag not in raw
    ]
    if missing:
        raise ValueError(f"{path} is not a base-called trace: it holds no {', '.join(missing)}")
    called = raw["PBAS2"].decode("ascii").upper()
    quality = tuple(raw["PCON2"])
    peaks = raw["PLOC2"]
    channels = {
        base: raw[f"DATA{9 + index}"] for index, base in enumerate(raw["FWO_1"].decode("ascii"))
    }
    bases = "".join(
        _mixed(
            base,
            peaks[index],
            peaks[index - 1] if index else 0,
            peaks[index + 1] if index + 1 < len(peaks) else len(channels[base]) - 1,
            channels,
        )
        if base in channels
        else base
        for index, base in enumerate(called)
    )
    trusted = _mott(quality) if any(quality) else (0, 0)
    return SequencingResult(name=Path(path).stem, bases=bases, quality=quality, trusted=trusted)


def _mixed(
    base: str, scan: int, before: int, after: int, channels: dict[str, Sequence[int]]
) -> str:
    """Return `base`, or the IUPAC code of it and the second peak that reaches its share."""
    first = channels[base][scan]
    if first <= 0:
        return base
    heights = {
        other: signal[scan] - max(min(signal[before : scan + 1]), min(signal[scan : after + 1]))
        for other, signal in channels.items()
        if other != base
    }
    second = max(heights, key=heights.__getitem__)
    return _PAIR_CODE[frozenset(base + second)] if heights[second] >= MIXED_SHARE * first else base


def _mott(quality: Sequence[int]) -> tuple[int, int]:
    """Return the half-open span Mott's rule keeps, as Biopython's ``abi-trim`` computes it."""
    if len(quality) <= 20:
        return 0, len(quality)
    total, best, start, end = 0.0, 0.0, None, 0
    for index in range(1, len(quality)):
        total += 0.05 - 10 ** (-quality[index] / 10)
        if total < 0:
            total = 0.0
        elif start is None:
            start = index
        if total > best:
            best, end = total, index
    return (start or 0), end

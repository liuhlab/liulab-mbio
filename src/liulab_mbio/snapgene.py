"""Read and write SnapGene ``.dna`` files as `SequenceRecord`.

A ``.dna`` file is a stream of packets: a type byte, a big-endian 4-byte length, then the
payload. SnapGene's 1-based inclusive ranges are converted here and nowhere else.
"""

import os
import struct
import xml.etree.ElementTree as ET
from collections.abc import Iterator
from pathlib import Path

from liulab_mbio.sequence import BindingSite, Feature, Primer, Segment, SequenceRecord, Strand

_SEQUENCE = 0x00
_PRIMERS = 0x05
_NOTES = 0x06
_COOKIE = 0x09
_FEATURES = 0x0A
_CIRCULAR = 0x01
#: SnapGene's name for the record on a map; the model holds it as `SequenceRecord.name`.
_MAP_LABEL = "CustomMapLabel"

_STRAND_OF_DIRECTIONALITY = {"1": Strand.FORWARD, "2": Strand.REVERSE, "3": Strand.BOTH}


def _span(text: str, length: int, *, base: int) -> tuple[int, int]:
    """Convert an inclusive ``"first-last"`` range; ``last < first`` runs across the origin."""
    first, last = (int(part) for part in text.split("-"))
    start, end = first - base, last - base + 1
    return (start, end) if end > start else (start, end + length)


def _qualifier_value(node: ET.Element) -> str | int:
    if (number := node.get("int")) is not None:
        return int(number)
    text, predef = node.get("text"), node.get("predef")
    if text is not None and predef is not None:
        return f"{predef}:{text}"
    return text if text is not None else predef or ""


def _read_features(payload: bytes, length: int) -> tuple[Feature, ...]:
    features = []
    for node in ET.fromstring(payload).iter("Feature"):
        parts = [
            (_span(segment.get("range", ""), length, base=1), segment)
            for segment in node.iter("Segment")
            if segment.get("type") != "gap"
        ]
        color = parts[0][1].get("color") if parts else None
        segments = tuple(
            Segment(
                start,
                end,
                name=segment.get("name", ""),
                color=None if segment.get("color") == color else segment.get("color"),
            )
            for (start, end), segment in parts
        )
        qualifiers = {
            q.get("name", ""): tuple(_qualifier_value(v) for v in q.iter("V"))
            for q in node.iter("Q")
        }
        features.append(
            Feature(
                node.get("name", ""),
                node.get("type", "misc_feature"),
                segments,
                strand=_STRAND_OF_DIRECTIONALITY.get(node.get("directionality", ""), Strand.NONE),
                qualifiers=qualifiers,
                color=color,
            )
        )
    return tuple(features)


def _packets(data: bytes) -> Iterator[tuple[int, bytes]]:
    offset = 0
    while offset < len(data):
        if offset + 5 > len(data):
            raise ValueError("truncated packet header")
        kind, length = struct.unpack_from(">BI", data, offset)
        payload = data[offset + 5 : offset + 5 + length]
        if len(payload) < length:
            raise ValueError(f"truncated packet 0x{kind:02X}")
        yield kind, payload
        offset += 5 + length


def _read_primers(payload: bytes, length: int) -> tuple[Primer, ...]:
    primers = []
    for node in ET.fromstring(payload).iter("Primer"):
        sites = tuple(
            BindingSite(
                *_span(site.get("location", ""), length, base=0),
                Strand.REVERSE if site.get("boundStrand") == "1" else Strand.FORWARD,
            )
            for site in node.iter("BindingSite")
            if site.get("simplified") != "1"
        )
        primers.append(
            Primer(
                node.get("name", ""),
                node.get("sequence", ""),
                binding_sites=sites,
                description=node.get("description", ""),
            )
        )
    return tuple(primers)


def _read_notes(payload: bytes) -> tuple[str, dict[str, str]]:
    """Return the map label, which the model holds as a name, and the remaining notes."""
    name, notes = "", {}
    for child in ET.fromstring(payload):
        if len(child):  # A note holding elements of its own, such as References.
            continue
        if child.tag == _MAP_LABEL:
            name = child.text or ""
        else:
            notes[child.tag] = child.text or ""
    return name, notes


def read_dna(path: str | os.PathLike[str]) -> SequenceRecord:
    """Read a SnapGene ``.dna`` file.

    SnapGene's HTML markup in note and qualifier text is kept as it is written.

    Raises
    ------
    ValueError
        If the file is not a SnapGene DNA file.
    """
    packets = list(_packets(Path(path).read_bytes()))
    if not packets or packets[0][0] != _COOKIE or not packets[0][1].startswith(b"SnapGene"):
        raise ValueError(f"{os.fspath(path)} is not a SnapGene file")
    sequences = [payload for kind, payload in packets if kind == _SEQUENCE]
    if len(sequences) != 1:
        raise ValueError(f"{os.fspath(path)} holds {len(sequences)} DNA sequence packets, not 1")
    flags, bases = sequences[0][0], sequences[0][1:].decode("ascii")
    by_kind = dict(packets)
    name, notes = _read_notes(by_kind[_NOTES]) if _NOTES in by_kind else ("", {})
    return SequenceRecord(
        bases,
        topology="circular" if flags & _CIRCULAR else "linear",
        name=name,
        features=_read_features(by_kind[_FEATURES], len(bases)) if _FEATURES in by_kind else (),
        primers=_read_primers(by_kind[_PRIMERS], len(bases)) if _PRIMERS in by_kind else (),
        notes=notes,
    )

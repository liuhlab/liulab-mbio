"""Read and write SnapGene ``.dna`` files as `SequenceRecord`.

A ``.dna`` file is a stream of packets: a type byte, a big-endian 4-byte length, then the
payload. SnapGene's 1-based inclusive ranges are converted here and nowhere else.

What the model does not hold — the cut-site cache, history, display settings, and the extra
attributes SnapGene writes on a feature — is kept in ``record.extras`` and written back
verbatim. A packet that describes the bases is dropped once the sequence changes, so that
SnapGene recomputes it.
"""

import os
import struct
import xml.etree.ElementTree as ET
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path

from liulab_mbio.sequence import (
    BindingSite,
    Feature,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    Topology,
    reverse_complement,
)

_SEQUENCE = 0x00
_PRIMERS = 0x05
_NOTES = 0x06
_COOKIE = 0x09
_FEATURES = 0x0A

#: Bit 0 of the sequence packet's flag byte. The others are strandedness and methylation.
_CIRCULAR = 0x01
_DOUBLE_STRANDED = 0x02

#: SnapGene's name for the record on a map; the model holds it as `SequenceRecord.name`.
_MAP_LABEL = "CustomMapLabel"
_DEFAULT_COLOR = "#a6acb3"
_COOKIE_PAYLOAD = b"SnapGene\x00\x01\x00\x0f\x00\x14"

#: The packets this module builds from the model rather than keeping.
_MODELLED = (_SEQUENCE, _FEATURES, _PRIMERS, _NOTES)

#: Packets that survive a sequence change. Everything else either describes the bases (the
#: ``0x03`` cut-site cache, the ``0x07`` and ``0x0B`` history, aligned sequences) or is
#: unknown, and SnapGene rebuilds it.
_SEQUENCE_INDEPENDENT = (_COOKIE, 0x08, 0x0D, 0x0E, 0x1C)

_STRAND_OF_DIRECTIONALITY = {"1": Strand.FORWARD, "2": Strand.REVERSE, "3": Strand.BOTH}
_DIRECTIONALITY_OF_STRAND = {strand: text for text, strand in _STRAND_OF_DIRECTIONALITY.items()}
_HYBRIDIZATION = (
    '<HybridizationParams minContinuousMatchLen="10" allowMismatch="1"'
    ' minMeltingTemperature="40" showAdditionalFivePrimeMatches="1"'
    ' minimumFivePrimeAnnealing="15"/>'
)


@dataclass(frozen=True, slots=True)
class _Kept:
    """What a file held that the model does not, ready for the writer."""

    sequence: str
    topology: Topology
    flags: int
    packets: tuple[tuple[int, bytes | None], ...]
    features: Mapping[Feature, str] = field(default_factory=dict)
    primers: Mapping[Primer, str] = field(default_factory=dict)
    hybridization: str = _HYBRIDIZATION
    notes_packet: bytes | None = None
    name: str = ""
    notes: Mapping[str, str] = field(default_factory=dict)


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


def _span(text: str, length: int, *, base: int) -> tuple[int, int]:
    """Convert an inclusive ``"first-last"`` range; ``last < first`` runs across the origin."""
    first, last = (int(part) for part in text.split("-"))
    start, end = first - base, last - base + 1
    return (start, end) if end > start else (start, end + length)


def _range(start: int, end: int, length: int, *, base: int) -> str:
    """Write a span as SnapGene's inclusive range, across the origin where it runs past."""
    last = end - 1 + base
    return f"{start + base}-{last - length if last >= length + base else last}"


def _qualifier_value(node: ET.Element) -> str | int:
    if (number := node.get("int")) is not None:
        return int(number)
    text, predef = node.get("text"), node.get("predef")
    if text is not None and predef is not None:
        return f"{predef}:{text}"
    return text if text is not None else predef or ""


def _read_features(payload: bytes, length: int) -> tuple[tuple[Feature, ...], dict[Feature, str]]:
    features, kept = [], {}
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
        feature = Feature(
            node.get("name", ""),
            node.get("type", "misc_feature"),
            segments,
            strand=_STRAND_OF_DIRECTIONALITY.get(node.get("directionality", ""), Strand.NONE),
            qualifiers=qualifiers,
            color=color,
        )
        features.append(feature)
        kept[feature] = ET.tostring(node, encoding="unicode")
    return tuple(features), kept


def _read_primers(payload: bytes, length: int) -> tuple[tuple[Primer, ...], dict[Primer, str], str]:
    root = ET.fromstring(payload)
    params = root.find("HybridizationParams")
    primers, kept = [], {}
    for node in root.iter("Primer"):
        sites = tuple(
            BindingSite(
                *_span(site.get("location", ""), length, base=0),
                Strand.REVERSE if site.get("boundStrand") == "1" else Strand.FORWARD,
            )
            for site in node.iter("BindingSite")
            if site.get("simplified") != "1"
        )
        primer = Primer(
            node.get("name", ""),
            node.get("sequence", ""),
            binding_sites=sites,
            description=node.get("description", ""),
        )
        primers.append(primer)
        kept[primer] = ET.tostring(node, encoding="unicode")
    hybridization = (
        ET.tostring(params, encoding="unicode") if params is not None else _HYBRIDIZATION
    )
    return tuple(primers), kept, hybridization


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
    topology: Topology = "circular" if flags & _CIRCULAR else "linear"
    by_kind = dict(packets)
    features, feature_xml = (
        _read_features(by_kind[_FEATURES], len(bases)) if _FEATURES in by_kind else ((), {})
    )
    primers, primer_xml, hybridization = (
        _read_primers(by_kind[_PRIMERS], len(bases))
        if _PRIMERS in by_kind
        else ((), {}, _HYBRIDIZATION)
    )
    name, notes = _read_notes(by_kind[_NOTES]) if _NOTES in by_kind else ("", {})
    kept = _Kept(
        sequence=bases,
        topology=topology,
        flags=flags,
        packets=tuple((kind, None if kind in _MODELLED else payload) for kind, payload in packets),
        features=feature_xml,
        primers=primer_xml,
        hybridization=hybridization,
        notes_packet=by_kind.get(_NOTES),
        name=name,
        notes=notes,
    )
    return SequenceRecord(
        bases,
        topology=topology,
        name=name,
        features=features,
        primers=primers,
        notes=notes,
        extras={"snapgene": kept},
    )


def write_dna(record: SequenceRecord, path: str | os.PathLike[str]) -> None:
    """Write `record` as a SnapGene ``.dna`` file.

    A feature or primer unchanged since the file it was read from is written back verbatim, so
    that what the model does not hold survives. A feature with no colour takes SnapGene's
    default grey.
    """
    Path(path).write_bytes(_dump(record))


def _dump(record: SequenceRecord) -> bytes:
    kept = record.extras.get("snapgene")
    kept = kept if isinstance(kept, _Kept) else None
    fresh = (
        kept is not None and kept.sequence == record.sequence and kept.topology == record.topology
    )
    layout = kept.packets if kept else ((_COOKIE, _COOKIE_PAYLOAD), *((k, None) for k in _MODELLED))
    written, out = set(), bytearray()
    for kind, payload in layout:
        if payload is None:
            written.add(kind)
            out += _packet(kind, _modelled_packet(kind, record, kept))
        elif fresh or kind in _SEQUENCE_INDEPENDENT:
            out += _packet(kind, payload)
    for kind, wanted in (
        (_SEQUENCE, True),
        (_FEATURES, bool(record.features)),
        (_PRIMERS, bool(record.primers)),
        (_NOTES, bool(record.notes or record.name)),
    ):
        if wanted and kind not in written:
            out += _packet(kind, _modelled_packet(kind, record, kept))
    return bytes(out)


def _packet(kind: int, payload: bytes) -> bytes:
    return struct.pack(">BI", kind, len(payload)) + payload


def _modelled_packet(kind: int, record: SequenceRecord, kept: _Kept | None) -> bytes:
    if kind == _SEQUENCE:
        flags = (kept.flags if kept else _DOUBLE_STRANDED) & ~_CIRCULAR
        flags |= _CIRCULAR if record.topology == "circular" else 0
        return bytes([flags]) + record.sequence.encode("ascii")
    if kind == _FEATURES:
        return _xml(_features_node(record, kept), declaration=True)
    if kind == _PRIMERS:
        return _xml(_primers_node(record, kept), declaration=True)
    return _notes_packet(record, kept)


def _xml(node: ET.Element, *, declaration: bool) -> bytes:
    text = ET.tostring(node, encoding="unicode")
    return (('<?xml version="1.0"?>' if declaration else "") + text).encode("utf-8")


def _features_node(record: SequenceRecord, kept: _Kept | None) -> ET.Element:
    root = ET.Element("Features", {"nextValidID": str(len(record.features))})
    for index, feature in enumerate(record.features):
        stored = kept.features.get(feature) if kept else None
        node = ET.fromstring(stored) if stored else _feature_node(feature, len(record))
        node.set("recentID", str(index))
        root.append(node)
    return root


def _feature_node(feature: Feature, length: int) -> ET.Element:
    node = ET.Element("Feature", {"recentID": "0", "name": feature.name})
    if (directionality := _DIRECTIONALITY_OF_STRAND.get(feature.strand)) is not None:
        node.set("directionality", directionality)
    node.set("type", feature.type)
    node.set("allowSegmentOverlaps", "0")
    node.set("consecutiveTranslationNumbering", "1")
    previous: int | None = None
    for segment in feature.segments:
        start = (
            segment.start
            if previous is None or segment.start >= previous
            else segment.start + length
        )
        if previous is not None and start > previous:
            ET.SubElement(
                node,
                "Segment",
                {
                    "range": _range(previous, start, length, base=1),
                    "color": "noColor",
                    "type": "gap",
                },
            )
        attributes = {"range": _range(segment.start, segment.end, length, base=1)}
        if segment.name:
            attributes["name"] = segment.name
        attributes["color"] = segment.color or feature.color or _DEFAULT_COLOR
        attributes["type"] = "standard"
        if feature.type == "CDS":
            attributes["translated"] = "1"
        ET.SubElement(node, "Segment", attributes)
        previous = start + segment.end - segment.start
    for name, values in feature.qualifiers.items():
        qualifier = ET.SubElement(node, "Q", {"name": name})
        for value in values:
            key = "int" if isinstance(value, int) else "text"
            ET.SubElement(qualifier, "V", {key: str(value)})
    return node


def _primers_node(record: SequenceRecord, kept: _Kept | None) -> ET.Element:
    root = ET.Element("Primers", {"nextValidID": str(len(record.primers))})
    root.append(ET.fromstring(kept.hybridization if kept else _HYBRIDIZATION))
    for index, primer in enumerate(record.primers):
        stored = kept.primers.get(primer) if kept else None
        node = ET.fromstring(stored) if stored else _primer_node(primer, record)
        node.set("recentID", str(index))
        root.append(node)
    return root


def _primer_node(primer: Primer, record: SequenceRecord) -> ET.Element:
    node = ET.Element(
        "Primer",
        {
            "recentID": "0",
            "name": primer.name,
            "sequence": primer.sequence,
            "description": primer.description,
        },
    )
    for site in primer.binding_sites:
        annealed = record.extract(Segment(site.start, site.end))
        if site.strand == Strand.REVERSE:
            annealed = reverse_complement(annealed)
        location = _range(site.start, site.end, len(record), base=0)
        binding = ET.SubElement(
            node,
            "BindingSite",
            {
                "location": location,
                "boundStrand": "1" if site.strand == Strand.REVERSE else "0",
                "annealedBases": annealed,
            },
        )
        if tail := primer.sequence[: len(primer.sequence) - (site.end - site.start)]:
            ET.SubElement(binding, "Component", {"bases": tail})
        ET.SubElement(binding, "Component", {"hybridizedRange": location, "bases": annealed})
    return node


def _notes_packet(record: SequenceRecord, kept: _Kept | None) -> bytes:
    unchanged = (
        kept is not None
        and kept.notes_packet is not None
        and kept.name == record.name
        and dict(kept.notes) == dict(record.notes)
    )
    if unchanged and kept is not None and kept.notes_packet is not None:
        return kept.notes_packet
    root = ET.Element("Notes")
    for tag, text in record.notes.items():
        ET.SubElement(root, tag).text = text
    if record.name:
        ET.SubElement(root, _MAP_LABEL).text = record.name
    return _xml(root, declaration=False)

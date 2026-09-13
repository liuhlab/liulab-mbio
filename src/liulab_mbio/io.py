"""Read a sequence file into the shared model, choosing the reader by its suffix."""

import os
from pathlib import Path
from typing import TYPE_CHECKING, cast

from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand, Topology

if TYPE_CHECKING:
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqRecord import SeqRecord

#: Suffixes Biopython reads for us, and the format name it knows them by.
_BIOPYTHON = {
    ".gb": "genbank",
    ".gbk": "genbank",
    ".genbank": "genbank",
    ".fa": "fasta",
    ".fasta": "fasta",
    ".fna": "fasta",
}

_STRANDS = {1: Strand.FORWARD, -1: Strand.REVERSE}
#: Qualifiers that name a feature, in the order they are tried.
_NAMING = ("label", "gene", "product", "locus_tag")
#: What Biopython calls a record that names itself nowhere.
_UNKNOWN = ("<unknown name>", "<unknown id>", "")


def read_record(path: str | os.PathLike[str]) -> SequenceRecord:
    """Read a ``.dna``, GenBank or FASTA file holding one sequence.

    Raises
    ------
    ValueError
        If the suffix names no reader, or the file holds no single record.
    """
    suffix = Path(path).suffix.lower()
    if suffix == ".dna":
        from liulab_mbio.snapgene import read_dna

        return read_dna(path)
    if (fmt := _BIOPYTHON.get(suffix)) is None:
        raise ValueError(f"no reader for a '{suffix.lstrip('.')}' file: {os.fspath(path)}")
    from Bio import SeqIO

    return _converted(SeqIO.read(os.fspath(path), fmt))


def _converted(record: "SeqRecord") -> SequenceRecord:
    sequence = str(record.seq)
    topology: Topology = (
        "circular" if record.annotations.get("topology") == "circular" else "linear"
    )
    name = str(next((text for text in (record.name, record.id) if text not in _UNKNOWN), ""))
    description = str(record.description) if record.description not in _UNKNOWN else ""
    return SequenceRecord(
        sequence,
        topology=topology,
        name=name,
        features=tuple(_feature(feature, len(sequence), topology) for feature in record.features),
        notes={"Description": description} if description != name else {},
    )


def _feature(feature: "SeqFeature", length: int, topology: Topology) -> Feature:
    qualifiers = {key: tuple(values) for key, values in feature.qualifiers.items()}
    name = ""
    for key in _NAMING:
        if values := qualifiers.get(key):
            name = str(values[0])
            if key == "label":
                del qualifiers[key]
            break
    strand = feature.location.strand if feature.location is not None else None
    return Feature(
        name,
        feature.type,
        _segments(feature, length, topology),
        strand=_STRANDS.get(strand or 0, Strand.NONE),
        qualifiers={
            key: tuple(int(value) if str(value).isdigit() else str(value) for value in values)
            for key, values in qualifiers.items()
        },
    )


def _segments(feature: "SeqFeature", length: int, topology: Topology) -> tuple[Segment, ...]:
    """Take the parts in top-strand order, joining the two a feature across the origin holds."""
    # Biopython's positions are ints, but an inexact one is a class its stubs do not narrow.
    parts = sorted(
        (cast("int", part.start), cast("int", part.end))
        for part in (feature.location.parts if feature.location is not None else ())
    )
    if topology == "circular" and len(parts) > 1 and parts[0][0] == 0 and parts[-1][1] == length:
        first, last = parts.pop(0), parts.pop()
        parts.append((last[0], length + first[1]))
        parts.sort()
    return tuple(Segment(start, end) for start, end in parts)

"""Read a sequence file, or one region of an indexed FASTA, into the shared model."""

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


def read_region(
    path: str | os.PathLike[str], sequence_name: str, start: int, end: int
) -> SequenceRecord:
    """Read one span of one sequence of a FASTA, through the `.fai` index beside it.

    Only the bytes the index points at are read, so a region of a chromosome costs the same on
    a whole genome as on a small file. `start` and `end` are 0-based and half-open: a `start`
    below zero reads from the first base, and an `end` past the sequence reads to its last. The
    record comes back linear and upper-case, named as the FASTA spells it.

    Raises
    ------
    FileNotFoundError
        If no index lies beside the FASTA.
    ValueError
        If the index names no such sequence, the span holds no base of it, or the FASTA does
        not read as the index describes.

    Examples
    --------
    >>> read_region("hg38.fa", "chr1", 156710401, 156710421).sequence  # doctest: +SKIP
    'ACCTAGGAGAAGTGGCCAGC'
    """
    fasta = Path(path)
    index = fasta.with_name(fasta.name + ".fai")
    if not index.is_file():
        raise FileNotFoundError(
            f"{fasta} has no index beside it: `genome assembly register` writes {index.name}"
        )
    length, offset, bases_per_line, bytes_per_line = _indexed(index, sequence_name)
    first, last = max(start, 0), min(end, length)
    if first >= last:
        raise ValueError(f"{sequence_name} is {length} bases: none of it lies in {start}-{end}")
    here = offset + _byte(first, bases_per_line, bytes_per_line)
    with fasta.open("rb") as handle:
        handle.seek(here)
        read = handle.read(offset + _byte(last, bases_per_line, bytes_per_line) - here)
    bases = read.replace(b"\r", b"").replace(b"\n", b"").decode()
    if len(bases) != last - first:
        raise ValueError(f"{fasta} does not read as {index.name} describes it")
    return SequenceRecord(bases, name=sequence_name)


def _indexed(index: Path, sequence_name: str) -> tuple[int, int, int, int]:
    """Return a sequence's length, offset, bases per line and bytes per line, from a `.fai`."""
    for row in index.read_text().splitlines():
        fields = row.split("\t")
        if fields[0] == sequence_name and len(fields) >= 5:
            length, offset, bases_per_line, bytes_per_line = (int(field) for field in fields[1:5])
            return length, offset, bases_per_line, bytes_per_line
    raise ValueError(f"{index} names no sequence {sequence_name!r}")


def _byte(position: int, bases_per_line: int, bytes_per_line: int) -> int:
    """Return how far into a sequence the base at `position` lies, its line endings counted."""
    return position // bases_per_line * bytes_per_line + position % bases_per_line


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

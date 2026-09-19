"""Read a sequence file, or one region of an indexed FASTA, into the shared model."""

import os
import re
from io import StringIO
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
#: SnapGene's colour note on a feature of one segment, which some types give a direction.
_COLOR_NOTE = re.compile(r"color: (#[0-9a-fA-F]{6})(?:; direction: (?:LEFT|RIGHT))?")
#: SnapGene's note on a feature of several segments, as the file wraps it: a segment a line.
_SEGMENT_LIST = re.compile(r'/note="(This [^"\n]* segments:\n[^"]*)"')
_LIST_HEADER = re.compile(r"This (?:\w+ )*feature has (\d+) segments:")
#: One listed segment: its 1-based inclusive range, its colour, and its name if it has one.
_LISTED_SEGMENT = re.compile(r"\d+: (\d+) \.\. (\d+) / (#[0-9a-fA-F]{6})(?: / (.+))?")


def read_record(path: str | os.PathLike[str]) -> SequenceRecord:
    """Read a ``.dna``, GenBank or FASTA file holding one sequence.

    The colours SnapGene writes into a GenBank file's notes become `Feature.color` and
    `Segment.color`, and those notes leave the qualifiers.

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

    text = Path(path).read_text()
    return _converted(SeqIO.read(StringIO(text), fmt), _segment_lists(text))


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


def _segment_lists(text: str) -> dict[str, list[str]]:
    """Map each segment list, as Biopython joins its lines, to the lines themselves.

    Only the line breaks tell a last segment's name from text SnapGene adds after the list.
    """
    lists = {}
    for match in _SEGMENT_LIST.finditer(text):
        lines = [line.strip() for line in match[1].splitlines() if line.strip()]
        lists[" ".join(lines)] = lines
    return lists


def _converted(record: "SeqRecord", lists: dict[str, list[str]]) -> SequenceRecord:
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
        features=tuple(
            _feature(feature, len(sequence), topology, lists) for feature in record.features
        ),
        notes={"Description": description} if description != name else {},
    )


def _feature(
    feature: "SeqFeature", length: int, topology: Topology, lists: dict[str, list[str]]
) -> Feature:
    qualifiers = {key: tuple(values) for key, values in feature.qualifiers.items()}
    name = ""
    for key in _NAMING:
        if values := qualifiers.get(key):
            name = str(values[0])
            if key == "label":
                del qualifiers[key]
            break
    plain = _segments(feature, length, topology)
    segments, color, notes = plain, None, []
    for note in map(str, qualifiers.get("note", ())):
        if match := _COLOR_NOTE.fullmatch(note):
            segments, color = plain, match[1]
        elif (listed := _listed(lists.get(note, []), plain, length)) is not None:
            color, segments, rest = listed
            if rest:
                notes.append(rest)
        else:
            notes.append(note)
    if notes:
        qualifiers["note"] = tuple(notes)
    else:
        qualifiers.pop("note", None)
    strand = feature.location.strand if feature.location is not None else None
    return Feature(
        name,
        feature.type,
        segments,
        strand=_STRANDS.get(strand or 0, Strand.NONE),
        qualifiers={
            key: tuple(int(value) if str(value).isdigit() else str(value) for value in values)
            for key, values in qualifiers.items()
        },
        color=color,
    )


def _listed(
    lines: list[str], segments: tuple[Segment, ...], length: int
) -> tuple[str, tuple[Segment, ...], str] | None:
    """Read a segment list into colours and names, as the ``.dna`` reader holds them.

    Returns the feature's colour, its segments, and the text after the list; `None` unless the
    list names every segment, and nothing else, each with a ``#rrggbb`` colour.
    """
    if not lines or (header := _LIST_HEADER.fullmatch(lines[0])) is None:
        return None
    count = int(header[1])
    listed: dict[tuple[int, int], tuple[str, str]] = {}
    for line in lines[1 : count + 1]:
        if (entry := _LISTED_SEGMENT.fullmatch(line)) is None:
            return None
        start, end = int(entry[1]) - 1, int(entry[2])
        listed[start, end if end > start else end + length] = entry[3], entry[4] or ""
    if sorted(listed) != sorted((segment.start, segment.end) for segment in segments):
        return None
    first, coloured = next(iter(listed.values()))[0], []
    for segment in segments:
        color, name = listed[segment.start, segment.end]
        coloured.append(
            Segment(segment.start, segment.end, name=name, color=None if color == first else color)
        )
    return first, tuple(coloured), " ".join(lines[count + 1 :])


def _segments(feature: "SeqFeature", length: int, topology: Topology) -> tuple[Segment, ...]:
    """Take the parts in the order the top strand reads them, joining two that meet at the origin.

    Biopython lists a join's parts in the order the feature reads, so a reverse-strand one's are
    turned round.
    """
    if (location := feature.location) is None:
        return ()
    spans: list[tuple[int, int]] = []
    for part in location.parts[::-1] if location.strand == -1 else location.parts:
        # Biopython's positions are ints, but an inexact one is a class its stubs do not narrow.
        start, end = cast("int", part.start), cast("int", part.end)
        if topology == "circular" and spans and spans[-1][1] == length and start == 0:
            spans[-1] = (spans[-1][0], length + end)
        else:
            spans.append((start, end))
    return tuple(Segment(start, end) for start, end in spans)

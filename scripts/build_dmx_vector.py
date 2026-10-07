#!/usr/bin/env python3
"""Rebuild the iGGA destination vector from its Addgene parent.

    pixi run python scripts/build_dmx_vector.py [--from PATH] [--out PATH]

The method's rounds run in a DMX vector, which the lab builds from DMX0001 by the steps
`docs/synthesis-and-assembly.md` lists under *Lab resources*. Those steps are this script, so the
record ships as a measurement of the parent rather than as a sequence someone typed:
`docs/research/dmx-destination.md` records what each one does and why, and every number the
script prints is reproduced from the parent on each run.

The parent is a `.dna` file under `reference_docs/`, which is not in the repository. Without it
the script says so and stops; it never reaches the network.

A DNA sequence carries no licence. The parent is deposited and the rebuild is the lab's own
molecule, which is why it ships: #287 settled the same question for the ccdB cassette.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

from liulab_mbio.edits import insert, replace
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.io import read_record
from liulab_mbio.sequence import Feature, Segment, SequenceRecord
from liulab_mbio.sites import CutSite, domesticate, find_sites
from liulab_synbio.igga.method import IGGA
from liulab_synbio.igga.vector import released_cargo, round_cassette

REPO = Path(__file__).resolve().parents[1]

#: DMX0001, Addgene 247434, the depositor's own map. 5,839 bp circular.
PARENT = REPO / (
    "reference_docs/synthesis_and_assembly/dmx/addgene/addgene-plasmid-247434-sequence-494299.dna"
)

#: Where the rebuilt record goes: the destination the worked example is planned against.
OUT = REPO / "docs/examples/ap1-library/vector.gb"

NAME = "DMX-iGGA"
DESCRIPTION = (
    "iGGA destination rebuilt from DMX0001 (Addgene 247434): the method's internal stuffer "
    "where the parent's cassette was, PmeI outboard of each BsaI site, and no BbsI left in the "
    "backbone. Built by scripts/build_dmx_vector.py."
)

#: Every enzyme the method names. No step may spell a new site of one of these.
METHOD = tuple(get_enzyme(name) for name in ("BsaI", "BsmBI", "BbsI", "PaqCI", "PmeI", "SrfI"))

#: Where each blunt site goes, as an offset from the cassette the parent gives up. The method
#: asks for one PmeI site between each BsaI cut and the primer that reads a well, and no more.
#: `gate.check_dmx_vector` computes those windows off the record and judges the result, which is
#: what holds these two offsets honest; section 3 of the research note says why these two.
BLUNT_BEFORE = 12
BLUNT_AFTER = 33

#: What a blunt chopper's site is spelled as when one is put in.
PMEI_SITE = "GTTTAAAC"


@dataclass(frozen=True)
class Change:
    """One base this script changed outside a coding sequence, and what it took away."""

    position: int
    before: str
    after: str
    enzyme: str


def rebuild(parent: SequenceRecord) -> tuple[SequenceRecord, list[Change]]:
    """Return the rebuilt destination, and every base changed outside a coding sequence.

    Four steps, in order. The parent's cassette is read off a digest and replaced by the method's
    internal stuffer, so a round opens the result on the overhangs its parts carry. A blunt site
    then goes outboard of each releasing cut. What is left of the method's enzymes in the
    backbone goes last, synonymously where a coding sequence spells it and by one stated base
    where none does.
    """
    span = released_cargo(parent)
    if span is None:
        raise ValueError(
            f"{parent.name or 'the parent'} gives up no piece on {IGGA.entry_overhang} and "
            f"{IGGA.scar_overhang}: this is not the DMX parent the rebuild starts from"
        )
    record, _ = replace(parent, span.start, span.end, round_cassette().bases)
    stuffer = _cassette(record)
    record, _ = insert(record, stuffer.start - BLUNT_BEFORE, PMEI_SITE)
    stuffer = _cassette(record)
    record, _ = insert(record, stuffer.end + BLUNT_AFTER, PMEI_SITE)
    record, changes = _domesticated(record)
    return _annotated(record), changes


def _cassette(record: SequenceRecord) -> Segment:
    """Where the piece a round gives up lies now, read off a digest after each edit."""
    found = released_cargo(record)
    if found is None:  # pragma: no cover - the step before this one put it there
        raise ValueError("the rebuilt record gives up no cassette, so nothing locates the edits")
    return found


def _domesticated(record: SequenceRecord) -> tuple[SequenceRecord, list[Change]]:
    """Take the method's enzymes out of the backbone, leaving the cassette's own sites alone.

    A site inside a coding sequence goes by a synonymous codon, which is `liulab_mbio.sites`'
    own rule. A site in no coding sequence is one that rule leaves, because changing it changes
    what the record spells; here it is in backbone nothing annotates, so this script changes one
    base and says which.
    """
    changes: list[Change] = []
    for enzyme in (get_enzyme("BbsI"), get_enzyme("BsaI")):
        record, _ = domesticate(record, (enzyme,), avoid=METHOD)
        for site in _strayed(record, enzyme):
            record, change = _one_base(record, site)
            changes.append(change)
    return record, changes


def _strayed(record: SequenceRecord, enzyme: Enzyme) -> tuple[CutSite, ...]:
    """Return the sites of `enzyme` the cassette is not read from: the ones cutting the backbone.

    A site that excises the cassette reads from within the piece, as the internal enzyme's two
    do, or from just outboard of it, as the releasing enzyme's two do. Either way it lies no
    further off than the enzyme's own cut reaches, which is what tells it from a stray.
    """
    cassette = _cassette(record)
    reach = enzyme.bottom_cut + len(enzyme.site)
    return tuple(
        site
        for site in find_sites(record, enzyme)
        if not cassette.start - reach <= site.span.start <= cassette.end + reach
    )


def _one_base(record: SequenceRecord, site: CutSite) -> tuple[SequenceRecord, Change]:
    """Change the first base of `site` that takes it away and spells no new method site.

    Read left to right over the site and A, C, G, T in turn, so the same parent gives the same
    record. A change is refused rather than guessed at where none of the twenty-odd candidates
    leaves the record free.
    """
    before = len(find_sites(record, METHOD))
    at = site.span.start % len(record)
    for offset in range(len(site.enzyme.site)):
        position = (at + offset) % len(record)
        was = record.sequence[position]
        for base in "ACGT":
            if base == was:
                continue
            tried, _ = replace(record, position, position + 1, base)
            if len(find_sites(tried, METHOD)) == before - 1:
                return tried, Change(position, was, base, site.enzyme.name)
    raise ValueError(
        f"no single base takes the {site.enzyme.name} site at {site.start} away without "
        "spelling another of this method's sites: this one needs a decision, not a script"
    )


def _annotated(record: SequenceRecord) -> SequenceRecord:
    """Name the record and annotate what the rebuild put in: the cassette and the blunt sites."""
    cassette = _cassette(record)
    features = [
        *record.features,
        Feature("internal stuffer", "misc_feature", (cassette,)),
        *(
            Feature("PmeI", "protein_bind", (site.span,))
            for site in find_sites(record, get_enzyme("PmeI"))
        ),
    ]
    return SequenceRecord(
        record.sequence,
        topology=record.topology,
        name=NAME,
        features=tuple(features),
        primers=record.primers,
        notes={"Description": DESCRIPTION},
    )


def write_genbank(record: SequenceRecord, path: Path) -> None:
    """Write `record` as GenBank, carrying its features and its topology.

    Imported here and not at module scope so that importing this script costs nothing.
    """
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    length = len(record)
    made: list[SeqFeature] = []
    for feature in record.features:
        strand = int(feature.strand) or None
        parts = [
            SimpleLocation(start, end, strand=strand)
            for segment in feature.segments
            for start, end in _wrapped(segment, length)
        ]
        # Biopython lists a join's parts in the order the feature reads, which the reader undoes.
        if strand == -1:
            parts.reverse()
        location = parts[0] if len(parts) == 1 else CompoundLocation(parts)
        made.append(SeqFeature(location, type=feature.type, qualifiers={"label": [feature.name]}))
    written = SeqRecord(
        Seq(record.sequence),
        id=record.name,
        name=record.name,
        description=record.notes.get("Description", ""),
        features=made,
        annotations={"molecule_type": "DNA", "topology": record.topology},
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        from Bio import SeqIO

        SeqIO.write(written, handle, "genbank")


def _wrapped(segment: Segment, length: int) -> list[tuple[int, int]]:
    """Cut a span running past the origin into the one or two GenBank can spell."""
    start, end = segment.start % length, segment.end
    if end - start <= length - start:
        return [(start, start + (end - segment.start))]
    return [(start, length), (0, end - length)]


def _check_written(record: SequenceRecord, path: Path) -> None:
    """Refuse to leave behind a file that reads back as a different record.

    The rebuild is only worth what the file carries, and a feature's segments survive a GenBank
    round trip only where they were written in the order the reader expects. The strand is not
    compared: GenBank cannot spell a feature belonging to neither strand, so one reads back
    forward.

    Raises
    ------
    ValueError
        Naming the first thing that came back differently.
    """
    again = read_record(path)

    def spelled(one: SequenceRecord) -> list[tuple[str, str, tuple[tuple[int, int], ...]]]:
        return [
            (held.name, held.type, tuple((int(s.start), int(s.end)) for s in held.segments))
            for held in one.features
        ]

    was, now = spelled(record), spelled(again)
    for left, right in zip(was, now, strict=False):
        if left != right:
            raise ValueError(f"{path} reads back {right} where the rebuild made {left}")
    if (again.sequence, again.topology, len(was)) != (record.sequence, record.topology, len(now)):
        raise ValueError(f"{path} does not read back as the record the rebuild made")


def main(argv: list[str] | None = None) -> int:
    """Rebuild the destination and write it, printing what the result measures."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from", dest="source", type=Path, default=PARENT, help="the parent")
    parser.add_argument("--out", type=Path, default=OUT, help="where the record goes")
    args = parser.parse_args(argv)
    if not args.source.exists():
        print(
            f"{args.source} is not here. The parent is a reference document, which this "
            "repository does not carry; reference_docs/README.md says where it comes from.",
            file=sys.stderr,
        )
        return 1
    record, changes = rebuild(read_record(args.source))
    write_genbank(record, args.out)
    _check_written(record, args.out)
    print(f"{args.out}: {record.name}, {len(record)} bp {record.topology}")
    for enzyme in METHOD:
        found = find_sites(record, enzyme)
        where = ", ".join(f"{one.start} {one.strand.name[0]}" for one in found)
        print(f"  {enzyme.name:6} {len(found)}  {where}")
    print(f"  cassette {_cassette(record)}")
    for change in changes:
        print(f"  {change.enzyme} at {change.position}: {change.before} -> {change.after}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Rebuild the shipped codon usage tables by counting codons in annotated genomes.

    python scripts/build_codon_usage.py [--cds PATH] [--out PATH]

A table here is a count of codons over the complete coding sequences of one genome. Counting a
genome is what makes the table shippable: the sequence and its annotation carry no restriction
on use or redistribution, and the counts are then this package's own measurement rather than
someone else's compilation. `docs/research/codon-usage.md` records that position, where each
genome comes from, and why no published codon usage database is redistributed instead.

A bacterium is fetched from NCBI's E-utilities as spliced CDS nucleotides. A mammal is read
from liulab-genome, one canonical transcript per gene, so the script runs where that assembly
and annotation are prepared, under an interpreter that imports `genome`.

Without `--cds` the bacterial sequences are downloaded. Tests count excerpts and never reach
the network or a genome.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from collections.abc import Callable, Iterable, Iterator, Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "src/liulab_mbio/data/codon_usage.json"

#: NCBI asks scripted callers to name themselves. No address is sent with it.
EFETCH = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id={accession}&rettype=fasta_cds_na&retmode=text&tool=liulab-mbio"
)

BASES = frozenset("ACGT")

#: Tags GENCODE puts on a coding sequence it says runs off one end.
INCOMPLETE = ('tag "cds_start_NF"', 'tag "cds_end_NF"')


@dataclass(frozen=True)
class Host:
    """One organism whose coding sequences a table is counted from.

    A host naming an `assembly` is read from liulab-genome with that `annotation`, and
    `accession` is the assembly's. Otherwise `accession` is the GenBank record fetched.
    """

    name: str
    organism: str
    taxid: int
    accession: str
    note: str = ""
    assembly: str = ""
    annotation: str = ""


#: The hosts shipped. Adding one is an entry here and a rerun.
HOSTS: tuple[Host, ...] = (
    Host(
        name="e-coli-k12",
        organism="Escherichia coli str. K-12 substr. MG1655",
        taxid=511145,
        accession="U00096.3",
        note="Whole genome, every complete CDS. Not a highly expressed reference set.",
    ),
    Host(
        name="human",
        organism="Homo sapiens",
        taxid=9606,
        accession="GCF_000001405.40",
        note=(
            "GRCh38 (hg38), GENCODE v50: one Ensembl canonical CDS per protein-coding gene. "
            "Not a highly expressed reference set."
        ),
        assembly="hg38",
        annotation="gencode_v50",
    ),
    Host(
        name="mouse",
        organism="Mus musculus",
        taxid=10090,
        accession="GCF_000001635.27",
        note=(
            "GRCm39 (mm39), GENCODE vM39: one Ensembl canonical CDS per protein-coding gene. "
            "Not a highly expressed reference set."
        ),
        assembly="mm39",
        annotation="gencode_vM39",
    ),
)


def fetch_cds(accession: str) -> str:
    """Download one genome's spliced coding sequences as FASTA."""
    with urllib.request.urlopen(EFETCH.format(accession=accession), timeout=300) as response:
        return response.read().decode("ascii", errors="replace")


def canonical_cds(gtf: Iterable[str]) -> dict[str, tuple[str, str, list[tuple[int, int]]]]:
    """Read each canonical protein-coding transcript's chromosome, strand and coding spans.

    A GENCODE transcript counts when it is tagged `Ensembl_canonical`, codes for protein, is
    not tagged as running off either end, and lies off the mitochondrion, which reads another
    genetic code. GENCODE writes the stop codon on its own row, and it is kept, so stops are
    counted as they are for a bacterium. Spans come back 0-based and half-open, in file order.
    """
    transcripts: dict[str, tuple[str, str, list[tuple[int, int]]]] = {}
    for line in gtf:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] not in ("CDS", "stop_codon"):
            continue
        chrom, start, end, strand, attributes = (fields[i] for i in (0, 3, 4, 6, 8))
        if (
            chrom == "chrM"
            or 'tag "Ensembl_canonical"' not in attributes
            or 'transcript_type "protein_coding"' not in attributes
            or any(tag in attributes for tag in INCOMPLETE)
        ):
            continue
        transcript = attributes.split('transcript_id "', 1)[1].split('"', 1)[0]
        _, _, spans = transcripts.setdefault(transcript, (chrom, strand, []))
        spans.append((int(start) - 1, int(end)))
    return transcripts


def spliced(
    fetch: Callable[[str, int, int, str], str],
    chrom: str,
    strand: str,
    spans: Iterable[tuple[int, int]],
) -> str:
    """Join one transcript's coding spans 5' to 3', upper-cased.

    `fetch` reads a span on a strand, the minus strand reverse-complemented.
    """
    ordered = sorted(spans, reverse=strand == "-")
    return "".join(fetch(chrom, start, end, strand) for start, end in ordered).upper()


def genome_cds(host: Host) -> str:
    """Splice one host's canonical coding sequences out of liulab-genome, as FASTA."""
    from genome import Genome, Region  # pyright: ignore[reportMissingImports]

    reference = Genome(host.assembly, progressbar=False)
    gtf = reference.default_gtf_path
    if gtf is None or gtf.stem != host.annotation:
        raise SystemExit(f"{host.assembly} has {gtf} as its annotation, not {host.annotation}")

    def fetch(chrom: str, start: int, end: int, strand: str) -> str:
        return str(reference.fetch_sequence(Region(chrom, start, end, strand)))

    with gtf.open() as lines:
        transcripts = canonical_cds(lines)
    return "".join(
        f">{transcript}\n{spliced(fetch, chrom, strand, spans)}\n"
        for transcript, (chrom, strand, spans) in transcripts.items()
    )


def records(fasta: str) -> Iterator[str]:
    """Yield each sequence of a FASTA file, upper-cased, without its header."""
    bases: list[str] = []
    for line in fasta.splitlines():
        if line.startswith(">"):
            if bases:
                yield "".join(bases).upper()
            bases = []
        elif line.strip():
            bases.append(line.strip())
    if bases:
        yield "".join(bases).upper()


def count_codons(fasta: str) -> tuple[dict[str, int], int]:
    r"""Count every codon of every complete coding sequence, and say how many were counted.

    A sequence whose length is not a whole number of codons, or that carries a base outside
    ACGT, is skipped rather than counted in part: a partial or ambiguous record would bias the
    table by the codons it does spell.

    Examples
    --------
    >>> count_codons(">one\\nATGAAATAA\\n")
    ({'ATG': 1, 'AAA': 1, 'TAA': 1}, 1)
    """
    counts: dict[str, int] = {}
    used = 0
    for sequence in records(fasta):
        if len(sequence) % 3 or set(sequence) - BASES:
            continue
        used += 1
        for index in range(0, len(sequence), 3):
            codon = sequence[index : index + 3]
            counts[codon] = counts.get(codon, 0) + 1
    return counts, used


def build(sequences: Mapping[str, str]) -> dict[str, Any]:
    """Turn each host's coding sequences into the shipped structure."""
    tables: list[dict[str, Any]] = []
    for host in HOSTS:
        counts, cds_count = count_codons(sequences[host.name])
        if not counts:
            raise SystemExit(f"no complete coding sequence found for {host.accession}")
        tables.append(
            {
                "name": host.name,
                "organism": host.organism,
                "taxid": host.taxid,
                "accession": host.accession,
                "note": host.note,
                "cds_count": cds_count,
                "codon_count": sum(counts.values()),
                "counts": {codon: counts.get(codon, 0) for codon in sorted(_every_codon())},
            }
        )
    return {
        "built": datetime.now(UTC).date().isoformat(),
        "source": {
            "title": (
                "NCBI E-utilities spliced coding sequences for a bacterium; liulab-genome "
                "assemblies and GENCODE annotations for a mammal"
            ),
            "url": EFETCH.format(accession="<accession>"),
            "retrieved": datetime.now(UTC).date().isoformat(),
            "terms": (
                "NCBI itself places no restrictions on the use or distribution of the data "
                "contained therein. EMBL-EBI itself places no additional restrictions on the "
                "use or redistribution of the data available via its Data Resources and Tools "
                "other than those provided by the original data owners."
            ),
        },
        "tables": tables,
    }


def _every_codon() -> Iterator[str]:
    """All 64 codons, so a codon absent from a genome is written as zero rather than missing."""
    for first in "ACGT":
        for second in "ACGT":
            for third in "ACGT":
                yield first + second + third


def main() -> int:
    """Write the data file, and say where it went."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cds", type=Path, help="a downloaded fasta_cds_na file for the bacterium")
    parser.add_argument("--out", type=Path, default=DATA)
    arguments = parser.parse_args()
    fetched = [host for host in HOSTS if not host.assembly]
    if arguments.cds is not None and len(fetched) > 1:
        raise SystemExit("--cds names one file and there is more than one host to fetch")
    sequences = {
        host.name: genome_cds(host)
        if host.assembly
        else arguments.cds.read_text()
        if arguments.cds
        else fetch_cds(host.accession)
        for host in HOSTS
    }
    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(json.dumps(build(sequences), indent=2) + "\n")
    print(f"wrote {arguments.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

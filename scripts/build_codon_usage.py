#!/usr/bin/env python3
"""Rebuild the shipped codon usage tables by counting codons in an annotated genome.

    pixi run python scripts/build_codon_usage.py [--cds PATH] [--out PATH]

A table here is a count of every codon over every complete coding sequence of one genome,
fetched from NCBI's E-utilities as spliced CDS nucleotides. Counting the genome is what makes
the table shippable: NCBI places no restrictions on the use or distribution of the sequence
data, and the counts are then this package's own measurement rather than someone else's
compilation. `docs/research/codon-usage.md` records that position and why no published codon
usage database is redistributed instead.

Without `--cds` the sequences are downloaded. Tests count an excerpt and never reach the
network.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from collections.abc import Iterator, Mapping
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


@dataclass(frozen=True)
class Host:
    """One organism whose genome a table is counted from."""

    name: str
    organism: str
    taxid: int
    accession: str
    note: str = ""


#: The hosts shipped. Adding one is a line here and a rerun.
HOSTS: tuple[Host, ...] = (
    Host(
        name="e-coli-k12",
        organism="Escherichia coli str. K-12 substr. MG1655",
        taxid=511145,
        accession="U00096.3",
        note="Whole genome, every complete CDS. Not a highly expressed reference set.",
    ),
)


def fetch_cds(accession: str) -> str:
    """Download one genome's spliced coding sequences as FASTA."""
    with urllib.request.urlopen(EFETCH.format(accession=accession), timeout=300) as response:
        return response.read().decode("ascii", errors="replace")


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
            "title": "NCBI E-utilities, spliced coding sequences of one GenBank record",
            "url": EFETCH.format(accession="<accession>"),
            "retrieved": datetime.now(UTC).date().isoformat(),
            "terms": (
                "NCBI itself places no restrictions on the use or distribution of the data "
                "contained therein."
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
    parser.add_argument("--cds", type=Path, help="a downloaded fasta_cds_na file for one host")
    parser.add_argument("--out", type=Path, default=DATA)
    arguments = parser.parse_args()
    if arguments.cds is not None and len(HOSTS) > 1:
        raise SystemExit("--cds names one file and there is more than one host to build")
    sequences = {
        host.name: arguments.cds.read_text() if arguments.cds else fetch_cds(host.accession)
        for host in HOSTS
    }
    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(json.dumps(build(sequences), indent=2) + "\n")
    print(f"wrote {arguments.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

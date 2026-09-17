"""Crowded records: pUC19 with its 99 unique 6+ cutters, which each map's cap is set to show whole,
and one crowded past the cap with the enzymes the package ships."""

import csv
from functools import cache
from pathlib import Path

from liulab_mbio.io import read_record
from liulab_mbio.plot import layers
from liulab_mbio.sequence import BindingSite, Primer, SequenceRecord, Strand

DATA = Path(__file__).parents[1] / "data"


@cache
def puc19() -> SequenceRecord:
    return read_record(DATA / "pUC19.dna")


@cache
def unique_6_cutters() -> tuple[tuple[str, int], ...]:
    """Each commercial enzyme with a site of six or more bases that cuts pUC19 once, and where.

    The cut is the 0-based boundary where it cuts the top strand. The package ships few of these
    enzymes, so the list is pinned rather than searched for.
    """
    with (DATA / "pUC19-unique-6-cutters.tsv").open(newline="") as table:
        return tuple(
            (row["name"], int(row["cut"])) for row in csv.DictReader(table, delimiter="\t")
        )


def puc19_crowded() -> tuple[layers.Item, ...]:
    """pUC19's features and primers, and a cut site for each position its unique 6+ cutters cut."""
    record = puc19()
    return (
        *layers.items(record, cut_sites=False),
        *layers.merge_cuts(unique_6_cutters(), len(record)),
    )


def ecori_crowd() -> SequenceRecord:
    """A circular record with fifty EcoRI sites close together, a primer among them, and two
    HindIII sites well apart from them.

    Each EcoRI site cuts after base 101, 109, and so on up to 493.
    """
    bases = list("A" * 3000)
    for start in [*range(100, 500, 8), 1700, 2000]:
        bases[start : start + 6] = "GAATTC" if start < 500 else "AAGCTT"
    primer = Primer("among them", "ACGT", binding_sites=(BindingSite(200, 220, Strand.FORWARD),))
    return SequenceRecord("".join(bases), topology="circular", name="crowd", primers=(primer,))

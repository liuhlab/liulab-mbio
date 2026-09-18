"""The worked-example records, built here rather than downloaded from a vendor.

Both vectors are put together from the att sequences `liulab_mbio.cloning.gateway.att` ships
and the two attP arms `docs/research/gateway-cloning.md` §2 prints, by the formulas that same
section gives: attL is attP's first 75 bases with the junction's region, and attR is the
region with the first 100 bases of attP's other arm. The insert is the GFP record already
under `tests/data/`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from liulab_mbio.sequence import SequenceRecord

#: attP1 and attP2 in full, 233 bases each, from US 7,670,823 B1 FIG. 9 as the note restates it.
ATTP = {
    1: (
        "CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTTTTATAATG"
        "CCAACTTTGTACAAAAAAGCTGAAC"
        "GAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTA"
        "AAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA"
    ),
    2: (
        "CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTCTTATAATG"
        "CCAACTTTGTACAAGAAAGCTGAAC"
        "GAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTA"
        "AAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA"
    ),
}

#: Where attP's recombination region sits inside it, and how much of the far arm attR keeps.
_REGION_AT = 75
_ARM_BP = 100

#: Filler standing in for a backbone: no att overlap reads in it, on either strand.
_BACKBONE = "TTGACCGCTTAAGGCTAACGTCAGTTGCATGCCTTAAGGCTTGACCAATGCGTTCAGGCT"

#: A stand-in for the ccdB and chloramphenicol genes a destination vector carries between its
#: att sites, and which the reaction throws away.
_CASSETTE = "ATGCAGTTTAAGGTTTACACCTATAAAAGAGAGAGCCGTTATCGTCTGTTTGTGGATGTACAGAGTGATATTTAA"


def att_site(name: str) -> str:
    """Return one whole att site, arms and all, built by the note's own formulas.

    attB is the recombination region alone, attP the 233 bp the patent prints, attL its first
    75 bases with the region, and attR the region with 100 bases of the other arm.
    """
    from liulab_mbio.cloning.gateway.att import REGIONS

    kind, number = name[:-1], int(name[-1])
    region, arm = REGIONS[name], ATTP[number]
    if kind == "attB":
        return region
    if kind == "attP":
        return arm
    if kind == "attL":
        return arm[:_REGION_AT] + region
    if kind == "attR":
        return region + arm[_REGION_AT + len(region) : _REGION_AT + len(region) + _ARM_BP]
    raise ValueError(f"no att site called {name!r}")


def entry_clone(insert: str, *, name: str = "pENTR-GFP", marker: str = "KanR") -> SequenceRecord:
    """Return an entry clone carrying `insert` between attL1 and attL2, with a backbone marker."""
    from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand, reverse_complement

    left, right = att_site("attL1"), reverse_complement(att_site("attL2"))
    bases = _BACKBONE + left + insert + right + _BACKBONE
    at = len(_BACKBONE) + len(left)
    return SequenceRecord(
        bases,
        topology="circular",
        name=name,
        features=(
            Feature(marker, "CDS", (Segment(0, len(_BACKBONE)),), strand=Strand.FORWARD),
            Feature("GFP", "CDS", (Segment(at, at + len(insert)),), strand=Strand.FORWARD),
        ),
    )


def destination_vector(*, name: str = "pDEST-test", marker: str = "AmpR") -> SequenceRecord:
    """Return a destination vector: a ccdB cassette between attR1 and attR2, and a marker."""
    from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand, reverse_complement

    promoter, site = "TAATACGACTCACTATAGGG", "AAGGAGAT"
    left, right = att_site("attR1"), reverse_complement(att_site("attR2"))
    lead = _BACKBONE + promoter + site
    bases = lead + left + _CASSETTE + right + _BACKBONE
    return SequenceRecord(
        bases,
        topology="circular",
        name=name,
        features=(
            Feature(marker, "CDS", (Segment(0, len(_BACKBONE)),), strand=Strand.FORWARD),
            Feature(
                "T7 promoter",
                "promoter",
                (Segment(len(_BACKBONE), len(_BACKBONE) + len(promoter)),),
                strand=Strand.FORWARD,
            ),
            Feature(
                "RBS",
                "RBS",
                (Segment(len(lead) - len(site), len(lead)),),
                strand=Strand.FORWARD,
            ),
            Feature(
                "ccdB",
                "CDS",
                (Segment(len(lead) + len(left), len(lead) + len(left) + len(_CASSETTE)),),
                strand=Strand.FORWARD,
            ),
        ),
    )

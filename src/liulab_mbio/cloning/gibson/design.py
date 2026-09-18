"""The overlap a Gibson junction takes: how long it is, and which part it is taken from.

An **overlap** is the bases two fragments both spell where they meet. One of them already spells
them; the other carries them as a 5' tail on its primer, so the junction adds nothing and reads
exactly what the two parts read.

Two rules, both the supplier's, through ``docs/research/gibson-assembly.md``:

- **Length and melting temperature are one rule, not two** (note §3 and §4). The product
  documents a band, and NEB's worked example lengthens the overlap inside it until the Wallace
  melting temperature reaches the floor: "The length of the overlap sequence is determined by
  the number of nucleotides needed to reach a Tm >= 48 C". So `overlap_before` and
  `overlap_after` return the shortest overlap in the band that reaches it, and the longest the
  band allows where nothing in it does.
- **Where it is taken from** (note §3). Both NEB manuals let either fragment carry it, and both
  ask for the whole overlap to come from the vector where the vector is PCR-generated and
  reused: "the entire overlap sequence must originate from the vector sequence and must be
  added to primers that will be used to amplify the insert". So a junction with the vector on
  one side takes its bases from the vector, and the insert tails them.

Coordinates are the model's, 0-based and half-open, and an overlap read off a circular vector
may run across the origin.
"""

from collections.abc import Iterator

from liulab_mbio.cloning.gibson.bench import OverlapRule
from liulab_mbio.sequence import Segment, SequenceRecord

#: The Wallace rule's degrees per base pair, which is the only way either NEB manual names for
#: an overlap's melting temperature: "AT pair = 2 C and GC pair = 4 C" (note §4).
WALLACE_AT = 2.0
WALLACE_GC = 4.0


def wallace_tm(bases: str) -> float:
    """Return the melting temperature of `bases` by the Wallace rule, °C.

    Two degrees for each A or T and four for each G or C. A base that is neither adds nothing,
    because the rule has nothing to say about one.

    Examples
    --------
    >>> wallace_tm("GCGCGCGCGCGCGCG")
    60.0
    >>> wallace_tm("ATATATATATATATA")
    30.0
    """
    upper = bases.upper()
    return WALLACE_AT * sum(upper.count(base) for base in "AT") + WALLACE_GC * sum(
        upper.count(base) for base in "GC"
    )


def overlap_before(record: SequenceRecord, at: int, rule: OverlapRule) -> str:
    """Return the overlap ending at `at`, taken from `record`'s own bases before that point.

    Raises
    ------
    ValueError
        If `record` holds fewer bases than the rule's longest overlap.

    Examples
    --------
    >>> rule = OverlapRule(4, 6, 18.0)
    >>> overlap_before(SequenceRecord("AAAAGGCCTTTT"), 8, rule)
    'AGGCC'
    """
    return _chosen((_window(record, at - size, size) for size in _sizes(record, rule)), rule)


def overlap_after(record: SequenceRecord, at: int, rule: OverlapRule) -> str:
    """Return the overlap beginning at `at`, taken from `record`'s own bases from that point.

    Raises
    ------
    ValueError
        If `record` holds fewer bases than the rule's longest overlap.

    Examples
    --------
    >>> rule = OverlapRule(4, 6, 18.0)
    >>> overlap_after(SequenceRecord("AAAAGGCCTTTT"), 3, rule)
    'AGGCC'
    """
    return _chosen((_window(record, at, size) for size in _sizes(record, rule)), rule)


def _sizes(record: SequenceRecord, rule: OverlapRule) -> range:
    """Return the overlap lengths the rule allows, shortest first.

    Raises
    ------
    ValueError
        If the record cannot spell the longest of them.
    """
    if len(record) < rule.longest:
        raise ValueError(
            f"{record.name or 'this record'} is {len(record)} bases, too short for the "
            f"{rule.shortest}-{rule.longest} bp overlap this assembly product documents"
        )
    return range(rule.shortest, rule.longest + 1)


def _window(record: SequenceRecord, start: int, size: int) -> str:
    """Return `size` bases of `record` from `start`, reading across the origin where it must.

    Raises
    ------
    ValueError
        If a linear record runs out before the overlap does.
    """
    if record.topology == "circular":
        first = start % len(record)
    elif start >= 0 and start + size <= len(record):
        first = start
    else:
        raise ValueError(
            f"an overlap of {size} bases does not fit {record.name or 'this record'} at {start}"
        )
    return record.extract(Segment(first, first + size))


def _chosen(candidates: Iterator[str], rule: OverlapRule) -> str:
    """Return the shortest candidate reaching the rule's floor, or the longest one offered.

    A product stating no melting temperature takes the shortest its band allows, which is what
    its manual asks for.
    """
    taken = ""
    for bases in candidates:
        taken = bases
        if rule.tm_floor is None or wallace_tm(bases) >= rule.tm_floor:
            return bases
    return taken

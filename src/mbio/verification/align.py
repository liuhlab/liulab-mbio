"""Laying one sequencing result on the record it should be, column by column.

Biopython's local aligner places the trusted bases, scoring both strands and aligning only the
better one in full, which keeps memory linear for the losing strand. A circular record is
doubled, so a result read across its origin lies in one piece, and every position comes back
into the record's own coordinates by one modulo (`docs/adr/0001-coordinates.md`).
`docs/research/sequencing-read-packages.md` sections 3.2 to 3.4 measured each choice. Which bases
count, and what the columns come to, is `mbio.verification.judge`'s.
"""

from dataclasses import dataclass
from itertools import pairwise

from mbio.sequence import Segment, SequenceRecord, Strand, reverse_complement
from mbio.verification.result import SequencingResult

#: Biopython's named scoring the strand was chosen with in the packages note, section 3.2.
SCORING = "blastn"


@dataclass(frozen=True, slots=True)
class Column:
    """One step of a result laid on a record.

    Parameters
    ----------
    start, end
        The record's bases it lies over, as an ADR 0001 span; empty where the result holds bases
        the record lacks, at the point they go in.
    bases
        What the result reads there, on the record's top strand: one base over one, the bases it
        inserts, or nothing over a deletion.
    behind
        The result's own indices behind it: the base it reads, the bases it inserts, or the two
        either side of a deletion.
    """

    start: int
    end: int
    bases: str
    behind: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class Aligned:
    """Where one result's trusted bases landed on a record.

    Parameters
    ----------
    strand
        The record's strand the result reads along.
    span
        The record's bases the result lies over, as an ADR 0001 span, at most one turn.
    placed
        Every column, in record order. The trusted bases the local alignment left off either end
        are laid without gaps against the record bases they face.
    """

    strand: Strand
    span: Segment
    placed: tuple[Column, ...]


def place(expected: SequenceRecord, result: SequencingResult) -> Aligned | None:
    """Return where `result`'s trusted bases lie on `expected`, or ``None`` if nowhere.

    The local alignment's ends are extended without gaps, so a mismatch near either end of the
    trusted bases reads as a substitution rather than a stretch left over. Trusted bases running
    past either end of a linear record face no base and are left out, because the record says
    nothing of what lies beyond it. Nothing trusted lies nowhere, and so do trusted bases no
    stretch of either strand scores above zero against.

    Examples
    --------
    >>> laid = place(SequenceRecord("GATTACAGGC"), SequencingResult("r", "TTACAG"))
    >>> laid.span, laid.strand.name
    (Segment(start=2, end=8, name='', color=None), 'FORWARD')
    """
    from Bio import Align

    first, last = result.trusted_span
    trusted = result.bases[first:last]
    if not trusted:
        return None
    n = len(expected)
    circular = expected.topology == "circular"
    target = expected.sequence * 2 if circular else expected.sequence
    aligner = Align.PairwiseAligner(scoring=SCORING)
    aligner.mode = "local"
    forward, reverse = trusted, reverse_complement(trusted)
    forward_score, reverse_score = aligner.score(target, forward), aligner.score(target, reverse)
    if max(forward_score, reverse_score) <= 0:
        return None
    strand = Strand.REVERSE if reverse_score > forward_score else Strand.FORWARD
    query = reverse if strand is Strand.REVERSE else forward
    coordinates = aligner.align(target, query)[0].coordinates
    assert coordinates is not None
    targets, queries = ([int(one) for one in row] for row in coordinates)

    def at(position: int) -> int:
        return position % n if circular else position

    def behind(start: int, end: int) -> tuple[int, ...]:
        if strand is Strand.REVERSE:
            return tuple(last - 1 - one for one in range(start, end))
        return tuple(first + one for one in range(start, end))

    def bases(target_start: int, query_start: int, width: int) -> list[Column]:
        return [
            Column(
                at(target_start + k),
                at(target_start + k) + 1,
                query[query_start + k],
                behind(query_start + k, query_start + k + 1),
            )
            for k in range(width)
        ]

    lead, tail = queries[0], len(query) - queries[-1]
    before = lead if circular else min(lead, targets[0])
    after = tail if circular else min(tail, n - targets[-1])
    placed = bases(targets[0] - before, lead - before, before)
    for (t0, t1), (q0, q1) in zip(pairwise(targets), pairwise(queries), strict=True):
        if t1 == t0:
            placed.append(Column(at(t0), at(t0), query[q0:q1], behind(q0, q1)))
        elif q1 == q0:
            placed.append(Column(at(t0), at(t0) + t1 - t0, "", behind(q0 - 1, q0 + 1)))
        else:
            placed.extend(bases(t0, q0, t1 - t0))
    placed.extend(bases(targets[-1], queries[-1], after))
    start = at(targets[0] - before)
    width = targets[-1] + after - targets[0] + before
    return Aligned(strand, Segment(start, start + min(width, n)), tuple(placed))

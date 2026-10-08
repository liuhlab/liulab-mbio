"""The two figures only this method needs, each chosen from what a build designed.

A function here returns a `liulab_mbio.protocol.model.Figure` and never a drawing, as
`liulab_mbio.protocol.figures` does: the records are ones `liulab_synbio.igga.plan` already
writes, and the page draws them. What each is equivalent to, and what it has to vary, is in
``docs/research/figure-sources.md``.

The method's own figures are here and not in `liulab_mbio` because each encodes one method's
choices: which enzyme opens a destination and which releases a donor, and which of three primer
roles each PCR pairs.
"""

from collections.abc import Sequence
from typing import Literal

from liulab_mbio.cloning.plan import PRODUCT_FILE
from liulab_mbio.protocol.figures import SOURCE, SOURCE_KEY
from liulab_mbio.protocol.model import Citation, Figure, counted
from liulab_mbio.sequence import Segment, across_the_origin
from liulab_synbio.igga.cargo import PoolPlan
from liulab_synbio.igga.rounds import ROUND_FILE, Round

__all__ = [
    "OLIGO_FILE",
    "SOURCE",
    "SOURCE_KEY",
    "PoolStage",
    "assembly_rows",
    "pool_pcr_figure",
]

#: What `liulab_synbio.igga.plan` calls the representative oligo it writes for the pool figure.
OLIGO_FILE = "oligo.dna"

#: The three roles one oligo spells, in the order it spells them: the batch's forward primer,
#: then the block's inner primer and the batch's outer primer, both on the bottom strand.
ROLES = ("forward", "inner", "outer")

#: Which of the two nested PCRs a pool figure lights: the outer pair that pulls a batch out of
#: the pool, or the inner pair that pulls one block out of that batch.
type PoolStage = Literal["PCR1", "PCR2"]

#: How many times the stuffer's own length a map draws either side of what it shows, so the
#: flanks it sits between are in the picture: either side of the stuffer in a vector's map, and
#: either side of the cassette in a round's. A margin for the drawing, measured from nothing.
STUFFER_MARGIN = 3

#: Where in the note each figure's published equivalent stands.
ASSEMBLY_CITATION = Citation(SOURCE_KEY, "section 4, iGGA")
POOL_CITATION = Citation(SOURCE_KEY, "section 1, Baker's three-primer PCR design")


def assembly_rows(
    rounds: Sequence[Round],
    *,
    lit: int,
    vector: str = "",
    at: str = "",
    caption: str = "",
    citation: Citation | None = ASSEMBLY_CITATION,
) -> Figure:
    """Return what one round opens and what it makes, as two rows, the part it joins lit.

    The library the round opens is drawn above the library it makes, which is how a figure shows
    one molecule becoming the next, and the part the round joined is lit. Both rows are drawn to
    the cassette wherever one span falls on both records, not to the kilobases of backbone the
    round leaves alone.

    The rows are the records `liulab_synbio.igga.rounds.write_records` writes, named as it names
    them, so the figure and the files cannot disagree. Round 1 opens the vector instead of a
    round, and a build writing no record of that vector draws the row it makes alone.

    Parameters
    ----------
    rounds
        Every round of the build, in order.
    lit
        Which round the step is at, counting from one: the two rows are that round's, and the
        part it joins is what the figure lights.
    vector
        What the plan calls the vector round 1 opens, relative to `at`; empty where it writes no
        record of it.
    at
        Where the records sit, relative to the directory the protocol is read from: empty where
        they sit beside it, ``"../"`` where the protocol is one directory down.
    caption
        What the figure shows. One is written from the round where this is empty.
    citation
        Where the published figure this is equivalent to was read.

    Raises
    ------
    ValueError
        If no round is given, or `lit` names no round of this build.
    """
    if not rounds:
        raise ValueError("an assembly figure draws at least one round")
    if not 1 <= lit <= len(rounds):
        raise ValueError(f"this build runs {len(rounds)} round(s), so round {lit} is not one")
    here = rounds[lit - 1]
    scheme = here.scheme
    made = at + (PRODUCT_FILE if lit == len(rounds) else ROUND_FILE.format(number=here.number))
    opened = _opened(rounds, lit, vector, at)
    return Figure(
        (opened, made) if opened else (made,),
        caption or _caption(rounds, here, opened=bool(opened)),
        span=_cassette(here),
        linear=True,
        enzymes=(scheme.internal_enzyme, scheme.external_enzyme),
        highlight=(here.part.name,),
        citation=citation,
    )


def pool_pcr_figure(
    plan: PoolPlan,
    *,
    stage: PoolStage = "PCR1",
    path: str = OLIGO_FILE,
    caption: str = "",
    citation: Citation | None = POOL_CITATION,
) -> Figure:
    """Return one oligo of the pool, the pair the named PCR reads it by lit and the third dimmed.

    The three roles sit on every oligo: the batch's forward primer at the 5' end, then the
    block's inner primer and the batch's outer primer at the 3' end, both on the bottom strand.
    PCR1 pairs the forward with the outer and pulls a whole batch out of the pool; PCR2 pairs the
    same forward with the inner and pulls one block out of that batch. Which pair is lit is what
    a step showing one of the two reactions needs.

    Parameters
    ----------
    plan
        The pool as designed, read for the pair each stage uses.
    stage
        Which PCR the step is at.
    path
        What the protocol calls the oligo record, relative to the directory it is read from.
    caption
        What the figure shows. One is written from the stage where this is empty.
    citation
        Where the diagram this is equivalent to was read.

    Raises
    ------
    ValueError
        If the pool holds no oligo, so there is no representative to draw.
    """
    if not plan.pool.oligos:
        raise ValueError("a pool figure draws one oligo, and this pool holds none")
    # Read off the oligo the plan writes, so the lit pair and the record drawn cannot disagree.
    names = plan.pool.oligos[0].primers
    if len(names) != len(ROLES):
        listed = ", ".join(ROLES)
        raise ValueError(f"this figure draws the {listed} roles, and the oligo spells {names}")
    forward, inner, outer = names
    pair = (forward, outer) if stage == "PCR1" else (forward, inner)
    return Figure(
        (path,),
        caption or _pool_caption(plan, stage, pair),
        linear=True,
        enzymes=(plan.pool.layout.cutter.name,),
        highlight=pair,
        citation=citation,
    )


def _opened(rounds: Sequence[Round], lit: int, vector: str, at: str) -> str:
    """Return the record the round opens: the round before's, or the vector round 1 opens.

    Empty where round 1 is drawn and the build writes no record of its vector.
    """
    if lit > 1:
        return at + ROUND_FILE.format(number=rounds[lit - 2].number)
    return at + vector if vector else ""


def _cassette(one: Round) -> tuple[int, int] | None:
    """Return the stretch both rows are drawn to: the cassette the round changes, and a margin.

    The two records count their bases alike up to the stuffer, so one span starts on the same
    base in both. Nothing is spanned where the round's edit ran across the origin, which moves
    the origin of what it made, or where the span would end past the record the round opened,
    which on a circle is a span across that record's own origin. Either draws the records whole.
    """
    if across_the_origin(Segment(one.excised.start, one.excised.end), len(one.destination)):
        return None
    margin = (one.stuffer.end - one.stuffer.start) * STUFFER_MARGIN
    start = max(0, one.entry.start - margin)
    end = min(len(one.product), one.scar.end + margin)
    return (start, end) if end <= len(one.destination) else None


def _caption(rounds: Sequence[Round], here: Round, *, opened: bool) -> str:
    """Return what the assembly figure shows: the round a step is at, and which row is which."""
    rows = "what it opens above, what it makes below" if opened else "what it makes"
    return (
        f"Round {here.number} of {len(rounds)} joins {here.part.name} at position "
        f"{here.position}, entering on {here.entry_overhang} and leaving {here.scar_overhang}: "
        f"{rows}"
    )


def _pool_caption(plan: PoolPlan, stage: PoolStage, pair: tuple[str, str]) -> str:
    """Return what the pool figure shows: which pair reads the oligo, in how many wells."""
    wells = len(plan.batches) if stage == "PCR1" else len(plan.inner_pairs)
    pulls = "one batch out of the pool" if stage == "PCR1" else "one block out of its batch"
    return (
        f"{stage} pairs {pair[0]} with {pair[1]} and pulls {pulls}, in {counted(wells, 'reaction')}"
    )

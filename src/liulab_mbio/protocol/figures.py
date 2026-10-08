"""The figure library: named figures a step shows, each chosen from the design it is drawn from.

A function here returns a `Figure` and never a drawing. A step holds a spec, an agent edits it in
the JSON like any other field, and `liulab_mbio.protocol.render` draws every figure the same way,
so a figure cannot go stale against the record it names.

What each figure is equivalent to, and what it has to vary, is in
``docs/research/figure-sources.md``; `SOURCE_KEY` is what a figure's citation keys, and a protocol
showing one puts `SOURCE` in its `sources`.
"""

from collections.abc import Iterable

from liulab_mbio.plot.sequence_view import LIMIT
from liulab_mbio.protocol.model import Citation, Figure, Source
from liulab_mbio.sequence import SequenceRecord

#: What a figure of this library cites, and the note it cites.
SOURCE_KEY = "figure-sources"
SOURCE = Source(
    "liulab-mbio, Figure sources for four designs",
    note="docs/research/figure-sources.md",
)

#: Where in the note a ligation figure's equivalent stands.
LIGATION_CITATION = Citation(SOURCE_KEY, "section 2, Cargo GGA ligation")

#: How many bases a ligation figure draws either side of the junction, by default. Enough to
#: carry each cut's recognition site and the codons around it, and short enough to read.
LIGATION_CONTEXT = 24


def ligation_figure(
    record: SequenceRecord,
    *,
    path: str,
    junction: tuple[int, int],
    enzymes: Iterable[str],
    caption: str,
    context: int = LIGATION_CONTEXT,
    highlight: Iterable[str] = (),
    citation: Citation | None = LIGATION_CITATION,
) -> Figure:
    """Return a Golden Gate junction at base level: both strands, each cut, the frame above.

    What `liulab_mbio.plot.sequence_view` draws over the span around one junction, which is what a
    reader needs to see that the four bases match and that the join stays in frame. Any method's
    junction reads the same way, so the figure is the package's and not one method's.

    Parameters
    ----------
    record
        The molecule the junction lies on, read for its length and whether it is circular.
    path
        What the protocol calls that record, relative to the directory the protocol is read from.
    junction
        Where the junction lies, 0-based and half-open, ending past the record's length across
        the origin. The overhang, or the whole cassette a step is opening.
    enzymes
        The enzymes whose cuts are drawn through both strands.
    caption
        What the figure shows, in the words the step uses.
    context
        How many bases are drawn either side of `junction`. A span past the end of a linear
        record stops at the end; one past the end of a circular record runs across the origin.
    highlight
        What the figure points at: a feature, a primer or an enzyme. Every other item dims.
    citation
        Where the published figure this is equivalent to was read.

    Raises
    ------
    ValueError
        If `context` is negative, or the span drawn is longer than a sequence view holds.

    Examples
    --------
    >>> figure = ligation_figure(
    ...     SequenceRecord("ACGT" * 20, name="v"),
    ...     path="product.dna",
    ...     junction=(36, 40),
    ...     enzymes=["BsaI"],
    ...     caption="The entry junction, before the part goes in",
    ...     context=8,
    ... )
    >>> figure.span, figure.sequence_view, figure.linear
    ((28, 48), True, True)
    """
    if context < 0:
        raise ValueError(f"a ligation figure draws 0 or more bases of context, not {context}")
    span = _around(junction, context, len(record), circular=record.topology == "circular")
    if span[1] - span[0] > LIMIT:
        raise ValueError(
            f"a ligation figure draws at most {LIMIT:,} bases and {caption!r} would draw "
            f"{span[1] - span[0]:,}: give a shorter junction or less context"
        )
    return Figure(
        (path,),
        caption,
        span=span,
        linear=True,
        sequence_view=True,
        enzymes=tuple(enzymes),
        highlight=tuple(highlight),
        citation=citation,
    )


def _around(
    junction: tuple[int, int], context: int, length: int, *, circular: bool
) -> tuple[int, int]:
    """Return the span `context` bases either side of `junction`, as its record allows.

    A linear record stops at its ends. A circular one runs across the origin, which `start` is
    counted back round to and `end` past, as `docs/adr/0001-coordinates.md` has it.
    """
    start, end = junction[0] - context, junction[1] + context
    if not circular:
        return max(0, start), min(length, end)
    if end - start >= length:
        return 0, length
    return (start + length, end + length) if start < 0 else (start, end)

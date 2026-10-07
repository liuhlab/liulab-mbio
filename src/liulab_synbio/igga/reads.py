"""The reads of a finished library, designed against the records the plan already simulated.

Three pairs. `linkage` carries a member's parts and its barcode block in one molecule, which is
what says a barcode still names its part. `representation` spans the block alone, which is what
makes repeating it after every bottleneck cheap. `final_representation` is that same read taken
on the other side of the move into the working vector, where only the forward anchor survives:
the reverse anchor moved with the vector, so the pair is designed again rather than reused.

Each span is read off the record — the cargo is what the external enzyme frees, the stuffer is
what the terminal round keeps — so nothing here states a coordinate the design did not already
carry. The amplicon length follows from the two binding sites. The platform follows from what
the read has to carry in one molecule and not from a length anyone published: linkage has to
put coding bases and barcodes in one read, representation never does.
"""

from collections.abc import Sequence
from dataclasses import KW_ONLY, dataclass
from typing import Literal

from liulab_mbio.edits import replace
from liulab_mbio.primers.design import design_pair
from liulab_mbio.primers.evaluation import evaluate_pair
from liulab_mbio.primers.placement import Placement, amplicon_sizes
from liulab_mbio.sequence import Primer, Segment, SequenceRecord
from liulab_synbio.igga.method import Scheme
from liulab_synbio.igga.rounds import Round
from liulab_synbio.igga.vector import Working

#: Which sequencing a read's amplicon asks for. A read that has to carry a member's coding bases
#: and its barcode block in one molecule is a long read; one spanning the block alone is a short
#: one. ``docs/research/synthesis-and-assembly-materials.md`` reads the two apart the same way.
type Platform = Literal["long read", "short read"]

#: How far either side of the anchor a read primer may sit, bases. The same allowance a colony
#: PCR primer is given in `liulab_mbio.bench.validation`, for the same reason: a placement wider
#: than the window it anchors in stops being an anchor.
ALLOWANCE = 10

#: How far outside the cargo a linkage or representation primer anneals. Far enough that the
#: primer sits wholly in the vector rather than straddling the overhang the cargo enters on.
FLANK = 60


@dataclass(frozen=True, slots=True)
class ReadPair:
    """One designed read of the finished library: the pair, how long it runs, and on what.

    Parameters
    ----------
    name
        What the read is called, as the step names it.
    forward, reverse
        The pair, designed against `template`.
    amplicon_length
        Bases the pair amplifies, tails included, read off the two binding sites.
    platform
        What that amplicon asks of the sequencing.
    """

    name: str
    _: KW_ONLY
    forward: Primer
    reverse: Primer
    amplicon_length: int
    platform: Platform

    @property
    def primers(self) -> tuple[Primer, Primer]:
        """The pair, forward first."""
        return self.forward, self.reverse


@dataclass(frozen=True, slots=True)
class ReadPairs:
    """Every read the finished library is read by.

    Parameters
    ----------
    linkage
        Across the whole cargo, in the library backbone. Read once.
    representation
        Across the barcode block alone, in the library backbone. Read after every bottleneck.
    final_representation
        The same read after the move into the working vector, or ``None`` where the build names
        no working vector and the library stays where the rounds built it.
    """

    linkage: ReadPair
    representation: ReadPair
    final_representation: ReadPair | None = None

    @property
    def designed(self) -> tuple[ReadPair, ...]:
        """Every pair designed, in the order the bench runs them."""
        made = (self.linkage, self.representation, self.final_representation)
        return tuple(one for one in made if one is not None)


def read_pairs(
    scheme: Scheme, rounds: Sequence[Round], working: Working | None = None
) -> ReadPairs:
    """Design every read of the library these rounds built.

    Raises
    ------
    ValueError
        If no round is given, if the external enzyme frees no cargo from the last round's
        product, or if no annealing region fits one of the placements.

    Examples
    --------
    >>> read_pairs(scheme, rounds).linkage.platform  # doctest: +SKIP
    'long read'
    """
    if not rounds:
        raise ValueError("a library is read on the rounds that built it, and there are none")
    final = rounds[-1]
    product = final.product
    cargo = _cargo(rounds)
    return ReadPairs(
        _linkage(product, cargo),
        _representation(product, cargo, final.stuffer),
        None if working is None else _final_representation(product, cargo, final, working),
    )


def _cargo(rounds: Sequence[Round]) -> Segment:
    """Return where the whole cargo lies in the last round's product.

    It opens on the overhang the first round's part entered the vector on and closes on the
    cloning scar past the last barcode, so both ends are read off the design rather than off a
    length. Each end is checked against the bases there, because a later round that ran across
    the origin would have moved them.

    Raises
    ------
    ValueError
        If either end does not spell the overhang the method joins on.
    """
    final = rounds[-1]
    product = final.product
    entry, closing = rounds[0].entry_overhang, rounds[0].scar_overhang
    scar = len(closing)
    width = sum(len(entry) + one.coding.end - one.coding.start for one in rounds)
    width += final.block.end - final.stuffer.start + scar
    end = final.block.end + scar
    start = (end - width) % len(product)
    for at, overhang, which in (
        (start, entry, "entry"),
        (end - scar, closing, "cloning scar"),
    ):
        found = product.bases(at, at + len(overhang))
        if found != overhang:
            raise ValueError(
                f"the {which} overhang of this library should lie at {at} and the product spells "
                f"{found!r} there, not {overhang!r}: the cargo has no span to read across"
            )
    return Segment(start, end if end > start else end + len(product))


def _linkage(product: SequenceRecord, cargo: Segment) -> ReadPair:
    """Design the read spanning the whole cargo, anchored in the vector either side of it."""
    start, end = cargo.start - FLANK, cargo.end + FLANK
    forward, reverse = design_pair(
        product,
        start,
        end,
        forward_placement=Placement(five_prime=_near(start, product)),
        reverse_placement=Placement(five_prime=_near(end, product)),
        forward_name="Linkage read forward",
        reverse_name="Linkage read reverse",
    )
    return _pair("linkage", forward, reverse, product, "long read")


def _representation(product: SequenceRecord, cargo: Segment, stuffer: Segment) -> ReadPair:
    """Design the read spanning the barcode block alone, anchored in the retained stuffer."""
    end = cargo.end + FLANK
    forward, reverse = design_pair(
        product,
        stuffer.start,
        end,
        forward_placement=_inside(stuffer),
        reverse_placement=Placement(five_prime=_near(end, product)),
        forward_name="Representation read forward",
        reverse_name="Representation read reverse",
    )
    return _pair("representation", forward, reverse, product, "short read")


def _final_representation(
    product: SequenceRecord, cargo: Segment, final: Round, working: Working
) -> ReadPair:
    """Design that read again on the library in the working vector, where the backbone changed.

    The forward anchor is the same retained stuffer and travels with the cargo; the reverse
    anchor was in the library backbone and is gone, so it is placed in the working vector past
    where the cargo lands.
    """
    moved, at = _in_working(product, cargo, working)
    stuffer = Segment(at + final.stuffer.start - cargo.start, at + final.stuffer.end - cargo.start)
    end = at + (cargo.end - cargo.start) + FLANK
    forward, reverse = design_pair(
        moved,
        stuffer.start,
        end,
        forward_placement=_inside(stuffer),
        reverse_placement=Placement(five_prime=_near(end, moved)),
        forward_name="Final representation read forward",
        reverse_name="Final representation read reverse",
    )
    return _pair("final representation", forward, reverse, moved, "short read")


def _in_working(
    product: SequenceRecord, cargo: Segment, working: Working
) -> tuple[SequenceRecord, int]:
    """Return the library in the working vector, and where the cargo starts in it.

    The final assembly as the gate already composes it: the external enzyme frees the cargo, the
    cargo enzyme opens the working vector, and the ligase joins the two. A cassette running
    across the origin leaves the cargo at the end of the record, as `replace` moves the origin.
    """
    hole = working.destination.stuffer
    bases = product.bases(cargo.start, cargo.end)
    moved, _ = replace(working.record, hole.start, hole.end, bases)
    across = hole.end > len(working.record)
    return moved, len(moved) - len(bases) if across else hole.start


def _pair(
    name: str, forward: Primer, reverse: Primer, template: SequenceRecord, platform: Platform
) -> ReadPair:
    """Measure a designed pair's amplicon and carry it with the pair.

    Raises
    ------
    ValueError
        If the pair makes no product on the record it was designed against.
    """
    report = evaluate_pair(forward, reverse, template)
    length = report.amplicon_length
    if length is None:
        sizes = amplicon_sizes(forward, reverse, template)
        if not sizes:
            raise ValueError(f"the {name} pair amplifies nothing on the record it was designed on")
        length = max(sizes)
    return ReadPair(
        name,
        forward=forward,
        reverse=reverse,
        amplicon_length=length,
        platform=platform,
    )


def _near(position: int, record: SequenceRecord) -> Segment:
    """Return the positions `ALLOWANCE` either way of one, wrapped round a circular origin."""
    return _span(position - ALLOWANCE, 2 * ALLOWANCE + 1, record)


def _inside(span: Segment) -> Placement:
    """Return a placement holding a forward primer wholly inside one span.

    Both ends are bounded, which is what makes the primer the same for every member: a read
    anchored past the retained stuffer would anneal to one member's own coding bases.
    """
    return Placement(
        five_prime=Segment(span.start, span.start + ALLOWANCE),
        three_prime=Segment(span.start + ALLOWANCE, span.end + 1),
    )


def _span(start: int, width: int, record: SequenceRecord) -> Segment:
    """Return `width` positions from `start`, wrapped on a circular record and clipped on a line."""
    if record.topology == "circular":
        start %= len(record)
        return Segment(start, start + width)
    return Segment(max(start, 0), min(start + width, len(record)))


def read_sheet(pairs: ReadPairs) -> str:
    """Return the designed read pairs as a tab-separated sheet, one primer a row."""
    rows = ["\t".join(("read", "name", "sequence", "length", "amplicon_bp", "platform"))]
    for one in pairs.designed:
        for primer in one.primers:
            rows.append(
                "\t".join(
                    (
                        one.name,
                        primer.name,
                        primer.sequence,
                        str(len(primer.sequence)),
                        str(one.amplicon_length),
                        one.platform,
                    )
                )
            )
    return "\n".join(rows) + "\n"

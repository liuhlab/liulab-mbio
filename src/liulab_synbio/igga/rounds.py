"""Simulate each round, and annotate the records a lab member reads.

One round is two digests and a ligation: the internal enzyme excises the internal stuffer of the
library built so far, the external enzyme releases a part from its synthesised block, and the two
ligate on the same pair of overhangs. The product becomes the next round's destination, and it
keeps the internal enzyme's sites, which is what lets the next round open it.
``docs/adr/0004-library-rounds.md`` says why none of this is `liulab_mbio.cloning.goldengate`.

Every overhang is read off a digest rather than stated, so a round joins what the ordered DNA will
join, and ends that do not match are refused naming both. The product is the destination with its
stuffer replaced by the part, which is what `liulab_mbio.edits.replace` does: the vector keeps its
own coordinates through every round, and what it annotates travels with them.

The records are the deliverable. Each round's product is annotated where its part's coding bases,
internal stuffer and barcode lie, and at the two junctions the ligation made, so a round can be
read before the next is run. A stuffer is drawn over what the internal enzyme excises, so the
round that opens it leaves none of it behind.

Nothing here carries a verdict. A round refuses what it cannot build, and what it did build is
judged by `liulab_synbio.igga.gate`, which reads the finished records rather than these
coordinates, so a design an agent composed is judged the same way.

Coordinates are the model's, 0-based and half-open, and a span across the origin of a circular
record ends past the record's length.
"""

import dataclasses
from collections.abc import Sequence
from dataclasses import KW_ONLY, dataclass
from pathlib import Path

from liulab_mbio.cloning.plan import PRODUCT_FILE
from liulab_mbio.edits import EditReport, ordered, replace
from liulab_mbio.sequence import (
    Feature,
    Segment,
    SequenceRecord,
    Strand,
    across_the_origin,
    span_at,
)
from liulab_mbio.sites import Fragment, digest, released
from liulab_mbio.snapgene import write_dna
from liulab_synbio.igga.method import Scheme
from liulab_synbio.igga.parts import Part

#: What a junction is drawn in. A feature built in code has no colour of its own, and
#: `liulab_mbio.snapgene` writes SnapGene's default grey for one that has none.
JUNCTION_COLOR = "#ff9900"

#: What each earlier round's intermediate is called there.
ROUND_FILE = "round-{number}.dna"


@dataclass(frozen=True, slots=True)
class Round:
    """One round of a library assembly: what it opened, what it joined, and what it made.

    Parameters
    ----------
    number
        Which round it is, counting from one. The barcode block's length is counted from it.
    part
        The member of this round's part list that the record represents.
    scheme
        The method the build is given, which says which enzyme does which job.
    terminal
        Whether this is the last round, whose stuffer the product keeps.
    destination
        The library the round opened: the vector for the first round, the round before's product
        after that.
    product
        What the round made, circular: the destination with its internal stuffer replaced by the
        part, annotated where every region of it lies.
    released
        The digest fragment the external enzyme releases from the part's block, in the block's own
        coordinates.
    excised
        The digest fragment the internal enzyme excises from the destination, in the destination's
        own coordinates. It ends past the destination's length when it runs across the origin.
    edit
        What replacing the stuffer did to the features and primer binding sites the destination
        carried, beyond shifting them.
    coding, stuffer, barcode
        Where the part's own coding bases, the internal stuffer it brought and its barcode lie in
        `product`. `stuffer` is what a further round would excise, which is what the record draws.
    block
        Where the barcode block lies in `product`: this round's barcode, then every earlier
        round's, each separated by the cloning scar.
    entry, scar
        The two junctions the ligation made: the entry overhang the part came in on, and the
        cloning scar at its other end.
    """

    number: int
    part: Part
    _: KW_ONLY
    scheme: Scheme
    terminal: bool
    destination: SequenceRecord
    product: SequenceRecord
    released: Fragment
    excised: Fragment
    edit: EditReport
    coding: Segment
    stuffer: Segment
    barcode: Segment
    block: Segment
    entry: Segment
    scar: Segment

    @property
    def position(self) -> str:
        """The position this round fills."""
        return self.part.position

    @property
    def entry_overhang(self) -> str:
        """The overhang the part entered on, read off the digest that released it."""
        return self.released.left_overhang

    @property
    def scar_overhang(self) -> str:
        """The cloning scar the part's other end left, read off that same digest."""
        return self.released.right_overhang

    @property
    def retained(self) -> Segment:
        """What the product keeps past this round's part: its stuffer and the barcode block.

        After the terminal round this is the region the whole construct reads through, and its
        length is the method's
        `liulab_synbio.igga.method.Scheme.retained_length`.
        """
        carried = _barcode_at(self.part, self.scheme) - self.part.coding.end
        return span_at(self.product, self.coding.end, carried + _length(self.block))


def assemble_round(
    destination: SequenceRecord,
    part: Part,
    scheme: Scheme,
    *,
    number: int = 1,
    terminal: bool = False,
    name: str = "",
) -> Round:
    """Open `destination`, release `part` from its block, and ligate the two into one circle.

    Both digests are simulated with the enzyme that does the joining; the method's blunt enzymes
    cut only the pieces the round throws away, so the product is the same with or without them.

    Parameters
    ----------
    destination
        The library the round opens, circular: the vector for the first round, and the round
        before's product after that.
    part
        The part this round appends, as `liulab_synbio.igga.parts.design_parts` built it.
    scheme
        The method the build is given.
    number
        Which round this is, counting from one. The barcode block's length is counted from it.
    terminal
        Whether this is the last round, after which the product reads through what it keeps.
    name
        What to call the product; the round number is added to it.

    Returns
    -------
    Round
        The product, the two fragments it was made from, and where each region of it lies.

    Raises
    ------
    ValueError
        If `destination` is not circular, if the external enzyme does not release one part from
        the block, if the part does not carry its barcode where the method puts it, or if the two
        molecules do not present the same pair of overhangs — which names both.
    KeyError
        If the method names an enzyme this package does not ship.
    """
    if destination.topology != "circular":
        raise ValueError(
            f"a round opens a circular destination and closes it again, and this one is "
            f"{destination.topology}"
        )
    donor = SequenceRecord(part.sequence, name=part.name)
    released = _only_released(donor, part, scheme)
    excised = _excised(destination, released, part, scheme, number)
    barcode_at = _barcode_at(part, scheme)
    bases = donor.bases(released.start, released.end)
    product, edit = replace(destination, excised.start, excised.end, bases)
    total = len(product)
    crossing = across_the_origin(Segment(excised.start, excised.end), len(destination))
    at = total - len(bases) if crossing else excised.start
    inner = _inner(donor, released, part, scheme)

    def moved(index: int) -> int:
        """Where a base of the part's block lands in the product."""
        return at + index - released.start

    coding = span_at(product, moved(part.coding.start), part.coding.end - part.coding.start)
    stuffer = span_at(product, moved(inner.start), inner.end - inner.start)
    barcode = span_at(product, moved(barcode_at), len(part.barcode))
    block = span_at(
        product,
        barcode.start,
        number * len(part.barcode) + (number - 1) * len(scheme.cloning_scar),
    )
    entry = span_at(product, at, len(released.left_overhang))
    scar = span_at(product, at + len(bases), len(released.right_overhang))
    drawn = _drawn(part, scheme, number, coding=coding, stuffer=stuffer, barcode=barcode)
    joins = _joins(part, number, entry=entry, scar=scar)
    titled = f"{name} round {number}".strip() if name else f"round {number}"
    return Round(
        number,
        part,
        scheme=scheme,
        terminal=terminal,
        destination=destination,
        product=ordered(
            dataclasses.replace(product, name=titled, features=(*product.features, *drawn, *joins))
        ),
        released=released,
        excised=excised,
        edit=edit,
        coding=coding,
        stuffer=stuffer,
        barcode=barcode,
        block=block,
        entry=entry,
        scar=scar,
    )


def assemble_rounds(
    destination: SequenceRecord,
    parts: Sequence[Part],
    scheme: Scheme,
    positions: Sequence[str],
    *,
    name: str = "",
) -> tuple[Round, ...]:
    """Run every round in turn, each round's product opening the next.

    Parameters
    ----------
    destination
        The vector the first round opens, circular. `liulab_synbio.igga.vector` makes one.
    parts
        One part a position, in the order the rounds fill them: `representative` picks a set.
    scheme
        The method the build is given.
    positions
        The positions, in that same order.
    name
        What to call each product; the round number is added to it.

    Returns
    -------
    tuple[Round, ...]
        One round each, in the order they run. The last round's product is the fully assembled
        construct, and the rounds before it are the intermediates.

    Raises
    ------
    ValueError
        If the parts are not one a position in the build's own order, or for any reason
        `assemble_round` refuses.
    KeyError
        If the method names an enzyme this package does not ship.
    """
    _check_parts(parts, positions)
    made: list[Round] = []
    for number, part in enumerate(parts, start=1):
        made.append(
            assemble_round(
                destination,
                part,
                scheme,
                number=number,
                terminal=number == len(parts),
                name=name,
            )
        )
        destination = made[-1].product
    return tuple(made)


def representative(parts: Sequence[Part], positions: Sequence[str]) -> tuple[Part, ...]:
    """Pick one part a position, the first of each part list, in the order the rounds run.

    Every member of a part list carries the same stuffers and differs only in what it codes for
    and in its barcode, so any member represents the round; the first is the one the sheet is
    ordered from first.

    Raises
    ------
    ValueError
        If no part fills one of the positions.
    """
    chosen: list[Part] = []
    for index, position in enumerate(positions):
        found = next((part for part in parts if part.index == index), None)
        if found is None:
            raise ValueError(f"no part fills position {position!r}")
        chosen.append(found)
    return tuple(chosen)


def write_records(rounds: Sequence[Round], directory: Path) -> tuple[Path, ...]:
    """Write every round's record into `directory` as SnapGene ``.dna``, and return the paths.

    Each round before the last is written as `ROUND_FILE`, so a round can be checked before
    the next is run, and the last round's product as `liulab_mbio.cloning.plan.PRODUCT_FILE`,
    which is the representative construct. The same rounds write the same bytes.

    Raises
    ------
    ValueError
        If no round is given.
    """
    if not rounds:
        raise ValueError("a library has at least one round to write")
    directory.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for one in rounds[:-1]:
        written.append(directory / ROUND_FILE.format(number=one.number))
        write_dna(one.product, written[-1])
    written.append(directory / PRODUCT_FILE)
    write_dna(rounds[-1].product, written[-1])
    return tuple(written)


def _length(span: Segment) -> int:
    """How many bases a span holds."""
    return span.end - span.start


def _check_parts(parts: Sequence[Part], positions: Sequence[str]) -> None:
    """Refuse parts that are not one a position, in the order the rounds fill them."""
    if len(parts) != len(positions):
        names = ", ".join(positions)
        raise ValueError(
            f"this build has {len(positions)} position(s) — {names} — and "
            f"{len(parts)} part(s) were given, one a round"
        )
    for index, part in enumerate(parts):
        if part.index != index:
            raise ValueError(
                f"part {part.name!r} fills position {part.position!r}, which is not position "
                f"{positions[index]!r}: the rounds run in the build's own order"
            )


def _only_released(donor: SequenceRecord, part: Part, scheme: Scheme) -> Fragment:
    """Return the one piece the external enzyme releases from a part's block.

    It is the one piece with an overhang at each end: the remnants either side carry the ends of
    the block itself, which a blunt enzyme cuts again so that neither can ligate back.

    Raises
    ------
    ValueError
        Unless exactly one piece has two overhangs.
    """
    pieces = released(donor, scheme.external)
    if len(pieces) != 1:
        raise ValueError(
            f"{scheme.external.name} releases {len(pieces)} piece(s) with an overhang at each end "
            f"from part {part.name!r}, where one of them is the part: the block was not built for "
            "this method"
        )
    return pieces[0]


def _inner(donor: SequenceRecord, released: Fragment, part: Part, scheme: Scheme) -> Fragment:
    """Return the piece the internal enzyme excises from a part's block, which a round reopens.

    Raises
    ------
    ValueError
        Unless exactly one such piece lies inside the released part.
    """
    pieces = [
        piece
        for piece in digest(donor, scheme.internal)
        if piece.left_overhang
        and piece.right_overhang
        and released.start <= piece.start
        and piece.end <= released.end
    ]
    if len(pieces) != 1:
        raise ValueError(
            f"{scheme.internal.name} excises {len(pieces)} piece(s) from part {part.name!r}, where "
            "one of them is the stuffer the next round opens it on"
        )
    return pieces[0]


def _excised(
    destination: SequenceRecord, released: Fragment, part: Part, scheme: Scheme, number: int
) -> Fragment:
    """Return the stuffer the part replaces: the piece whose two overhangs are the part's own.

    Raises
    ------
    ValueError
        If no piece presents both of them, naming what each molecule offers, or if more than one
        does.
    """
    wanted = (released.left_overhang, released.right_overhang)
    pieces = digest(destination, scheme.internal)
    matching = [piece for piece in pieces if (piece.left_overhang, piece.right_overhang) == wanted]
    if len(matching) == 1:
        return matching[0]
    where = f"round {number} cannot ligate part {part.name!r}"
    if matching:
        raise ValueError(
            f"{where}: {scheme.internal.name} excises {len(matching)} pieces from the destination "
            f"on {wanted[0]} and {wanted[1]}, so which one the part replaces is ambiguous"
        )
    offered = (
        ", ".join(
            f"{piece.left_overhang or '-'} and {piece.right_overhang or '-'}" for piece in pieces
        )
        or "nothing at all"
    )
    raise ValueError(
        f"{where}: {scheme.external.name} releases the part on {wanted[0]} and {wanted[1]}, and "
        f"{scheme.internal.name} opens the destination on {offered}: a round ligates only where "
        "both overhangs match"
    )


def _barcode_at(part: Part, scheme: Scheme) -> int:
    """Where a part's barcode begins in its block, checked against the barcode itself.

    Raises
    ------
    ValueError
        If the part does not spell its barcode where the 3' external stuffer puts it.
    """
    at = part.length - len(scheme.external_stuffer_3) - len(part.barcode)
    if part.sequence[at : at + len(part.barcode)] != part.barcode:
        raise ValueError(
            f"part {part.name!r} does not carry its barcode where position "
            f"{part.position!r} puts one: the part and the method disagree about its block"
        )
    return at


def _drawn(
    part: Part,
    scheme: Scheme,
    number: int,
    *,
    coding: Segment,
    stuffer: Segment,
    barcode: Segment,
) -> tuple[Feature, ...]:
    """Draw what the part brought: its coding bases, the stuffer it carries and its barcode."""
    return (
        Feature(part.name, "CDS", (coding,), strand=Strand.FORWARD),
        Feature(
            f"{part.position} internal stuffer",
            "misc_feature",
            (stuffer,),
            qualifiers={
                "note": (f"round {number}: {scheme.internal.name} excises this to open it",)
            },
        ),
        Feature(
            f"{part.name} barcode",
            "misc_feature",
            (barcode,),
            qualifiers={
                "note": (f"round {number}: names {part.name} at position {part.position}",)
            },
        ),
    )


def _joins(part: Part, number: int, *, entry: Segment, scar: Segment) -> tuple[Feature, ...]:
    """Draw the two junctions the ligation made."""
    return (
        Feature(
            f"{part.position} entry junction",
            "misc_feature",
            (entry,),
            color=JUNCTION_COLOR,
            qualifiers={"note": (f"round {number}: {part.name} enters the library here",)},
        ),
        Feature(
            "cloning scar",
            "misc_feature",
            (scar,),
            color=JUNCTION_COLOR,
            qualifiers={
                "note": (f"round {number}: {part.name} joins what the rounds before it built",)
            },
        ),
    )

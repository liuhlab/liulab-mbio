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

Coordinates are the model's, 0-based and half-open, and a span across the origin of a circular
record ends past the record's length.
"""

import dataclasses
from collections.abc import Sequence
from dataclasses import KW_ONLY, dataclass
from pathlib import Path

from liulab_mbio.checks import Check, Status, worst_of
from liulab_mbio.cloning.plan import PRODUCT_FILE
from liulab_mbio.edits import EditReport, ordered, replace
from liulab_mbio.library.parts import Part
from liulab_mbio.library.scheme import Scheme
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand
from liulab_mbio.sites import Fragment, digest, find_sites
from liulab_mbio.snapgene import write_dna
from liulab_mbio.translate import translate

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
        The architecture the build is given, which says which enzyme does which job.
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
    def terminal(self) -> bool:
        """Whether this is the last round, whose stuffer the product keeps."""
        return self.part.index == self.scheme.position_count - 1

    @property
    def retained(self) -> Segment:
        """What the product keeps past this round's part: its stuffer and the barcode block.

        After the terminal round this is the region the whole construct reads through, and its
        length is the scheme's `liulab_mbio.library.scheme.Scheme.retained_length`.
        """
        carried = _barcode_at(self.part, self.scheme) - self.part.coding.end
        return _span(self.coding.end, carried + _length(self.block), len(self.product))

    @property
    def checks(self) -> tuple[Check, ...]:
        """Judge the product, as data a protocol can print.

        Whether the next round can open it comes first, then whether the enzyme that released the
        part still cuts it, and after the terminal round whether what the product keeps past its
        last part reads in frame without a stop.
        """
        made = [self._opens(), self._external()]
        if self.terminal:
            made.append(self._frame())
        return tuple(made)

    @property
    def status(self) -> Status:
        """The worst status of any check."""
        return worst_of(self.checks)

    def __getitem__(self, name: str) -> Check:
        """Return the check of that name.

        Raises
        ------
        KeyError
            If no check has it.
        """
        for check in self.checks:
            if check.name == name:
                return check
        raise KeyError(name)

    def _opens(self) -> Check:
        """Whether the internal enzyme still excises the stuffer this round's part brought.

        Every piece of a circular digest is bounded by two cuts, so what says the next round can
        open the product is that one of those pieces is that stuffer itself.
        """
        enzyme = self.scheme.internal
        pieces = digest(self.product, enzyme)
        found = [
            piece
            for piece in pieces
            if (piece.start, piece.end) == (self.stuffer.start, self.stuffer.end)
        ]
        where = f"{self.stuffer.start}-{self.stuffer.end}"
        if len(found) == 1:
            detail = (
                f"{enzyme.name} excises the stuffer at {where}, on {found[0].left_overhang} and "
                f"{found[0].right_overhang}"
            )
        else:
            detail = (
                f"{enzyme.name} cuts the product in {len(pieces)} place(s), none of which excises "
                f"the stuffer at {where}"
            )
        return Check("opens", "pass" if len(found) == 1 else "fail", len(found), detail)

    def _external(self) -> Check:
        """Whether the enzyme that released the part is gone from the product, as it should be."""
        enzyme = self.scheme.external
        left = len(find_sites(self.product, enzyme))
        return Check(
            "external sites",
            "pass" if left == 0 else "fail",
            left,
            f"{enzyme.name} no longer cuts the product"
            if left == 0
            else f"{enzyme.name} reads {left} site(s) the external stuffers should have taken away",
        )

    def _frame(self) -> Check:
        """Whether the retained region is whole codons and spells no stop in the product's frame."""
        bases = self.product.extract(self.retained)
        if len(bases) % 3:
            return Check(
                "reading frame",
                "fail",
                len(bases),
                f"the {len(bases)} bases the product keeps past its last part are not a whole "
                "number of codons",
            )
        spelled = translate(bases)
        stops = spelled.count("*")
        return Check(
            "reading frame",
            "pass" if stops == 0 else "fail",
            stops,
            f"{len(bases)} bases, {len(spelled)} codons, no stop"
            if stops == 0
            else f"{stops} stop codon(s) in the {len(spelled)} the product keeps past its last part",
        )


def assemble_round(
    destination: SequenceRecord,
    part: Part,
    scheme: Scheme,
    *,
    number: int = 1,
    name: str = "",
) -> Round:
    """Open `destination`, release `part` from its block, and ligate the two into one circle.

    Both digests are simulated with the enzyme that does the joining; the scheme's blunt enzymes
    cut only the pieces the round throws away, so the product is the same with or without them.

    Parameters
    ----------
    destination
        The library the round opens, circular: the vector for the first round, and the round
        before's product after that.
    part
        The part this round appends, as `liulab_mbio.library.parts.design_parts` built it.
    scheme
        The architecture the build is given.
    number
        Which round this is, counting from one. The barcode block's length is counted from it.
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
        the block, if the part does not carry its barcode where the scheme puts it, or if the two
        molecules do not present the same pair of overhangs — which names both.
    KeyError
        If the scheme names an enzyme this package does not ship.
    """
    if destination.topology != "circular":
        raise ValueError(
            f"a round opens a circular destination and closes it again, and this one is "
            f"{destination.topology}"
        )
    donor = SequenceRecord(part.sequence, name=part.name)
    released = _released(donor, part, scheme)
    excised = _excised(destination, released, part, scheme, number)
    barcode_at = _barcode_at(part, scheme)
    bases = donor.extract(Segment(released.start, released.end))
    product, edit = replace(destination, excised.start, excised.end, bases)
    total = len(product)
    at = excised.start if excised.end <= len(destination) else total - len(bases)
    inner = _inner(donor, released, part, scheme)

    def moved(index: int) -> int:
        """Where a base of the part's block lands in the product."""
        return at + index - released.start

    coding = _span(moved(part.coding.start), part.coding.end - part.coding.start, total)
    stuffer = _span(moved(inner.start), inner.end - inner.start, total)
    barcode = _span(moved(barcode_at), len(part.barcode), total)
    block = _span(
        barcode.start,
        number * scheme.barcode_length + (number - 1) * len(scheme.cloning_scar),
        total,
    )
    entry = _span(at, len(released.left_overhang), total)
    scar = _span(at + len(bases), len(released.right_overhang), total)
    drawn = _drawn(part, scheme, number, coding=coding, stuffer=stuffer, barcode=barcode)
    joins = _joins(part, number, entry=entry, scar=scar)
    titled = f"{name} round {number}".strip() if name else f"round {number}"
    return Round(
        number,
        part,
        scheme=scheme,
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
    *,
    name: str = "",
) -> tuple[Round, ...]:
    """Run every round in turn, each round's product opening the next.

    Parameters
    ----------
    destination
        The vector the first round opens, circular. `liulab_mbio.library.vector` makes one.
    parts
        One part a position, in the order the rounds fill them: `representative` picks a set.
    scheme
        The architecture the build is given.
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
        If the parts are not one a position in the scheme's own order, or for any reason
        `assemble_round` refuses.
    KeyError
        If the scheme names an enzyme this package does not ship.
    """
    _check_parts(parts, scheme)
    made: list[Round] = []
    for number, part in enumerate(parts, start=1):
        made.append(assemble_round(destination, part, scheme, number=number, name=name))
        destination = made[-1].product
    return tuple(made)


def representative(parts: Sequence[Part], scheme: Scheme) -> tuple[Part, ...]:
    """Pick one part a position, the first of each part list, in the order the rounds run.

    Every member of a part list carries the same stuffers and differs only in what it codes for
    and in its barcode, so any member represents the round; the first is the one the sheet is
    ordered from first.

    Raises
    ------
    ValueError
        If no part fills one of the scheme's positions.
    """
    chosen: list[Part] = []
    for index, position in enumerate(scheme.positions):
        found = next((part for part in parts if part.index == index), None)
        if found is None:
            raise ValueError(f"no part fills position {position.name!r}")
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


def _span(start: int, length: int, total: int) -> Segment:
    """Return a span of `length` bases from `start`, reduced onto a record of `total` bases."""
    first = start % total
    return Segment(first, first + length)


def _check_parts(parts: Sequence[Part], scheme: Scheme) -> None:
    """Refuse parts that are not one a position, in the order the rounds fill them."""
    if len(parts) != scheme.position_count:
        names = ", ".join(position.name for position in scheme.positions)
        raise ValueError(
            f"this scheme has {scheme.position_count} position(s) — {names} — and "
            f"{len(parts)} part(s) were given, one a round"
        )
    for index, part in enumerate(parts):
        if part.index != index:
            raise ValueError(
                f"part {part.name!r} fills position {part.position!r}, which is not position "
                f"{scheme.positions[index].name!r}: the rounds run in the scheme's own order"
            )


def _released(donor: SequenceRecord, part: Part, scheme: Scheme) -> Fragment:
    """Return the piece the external enzyme releases from a part's block.

    It is the one piece with an overhang at each end: the remnants either side carry the ends of
    the block itself, which a blunt enzyme cuts again so that neither can ligate back.

    Raises
    ------
    ValueError
        Unless exactly one piece has two overhangs.
    """
    pieces = [
        piece
        for piece in digest(donor, scheme.external)
        if piece.left_overhang and piece.right_overhang
    ]
    if len(pieces) != 1:
        raise ValueError(
            f"{scheme.external.name} releases {len(pieces)} piece(s) with an overhang at each end "
            f"from part {part.name!r}, where one of them is the part: the block was not built for "
            "this scheme"
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
        If the part does not spell its barcode where its position's 3' external stuffer puts it.
    """
    at = part.length - len(scheme.positions[part.index].external_stuffer_3) - len(part.barcode)
    if part.sequence[at : at + len(part.barcode)] != part.barcode:
        raise ValueError(
            f"part {part.name!r} does not carry its barcode where position "
            f"{part.position!r} puts one: the part and the scheme disagree about its block"
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

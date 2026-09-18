"""Join two digest pieces into one circular product, and say what each junction now spells.

Two rules run through the module:

- **An overhang is written on the top strand, so two pieces join by concatenating their spans.**
  `liulab_mbio.sites.digest` cuts at top-strand positions, so the spans of the pieces that go on
  lay end to end; `liulab_mbio.overhangs.compatible` is what says they may. Where the bases they
  paired on lie is a second question, answered by which strand overhangs: a 5' overhang is the
  following piece's first bases, a 3' overhang the last bases of the piece before.
- **This method's junction is not scarless.** The two ends came from a recognition site, and
  ligating them puts that site back, so the product gains those bases. A junction says which
  site it spells, read off the product rather than assumed.

Coordinates are the model's, 0-based and half-open. The product keeps the vector's origin, so
the vector's own coordinates still read true and no junction sits at base zero.
"""

import dataclasses
from dataclasses import dataclass

from liulab_mbio.checks import Check, Status, worst_of
from liulab_mbio.cloning.restriction.digest import Piece, closes
from liulab_mbio.edits import carried, ordered, rotate
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.sequence import Feature, Primer, Segment, SequenceRecord, reverse_complement
from liulab_mbio.sites import CutSite, find_sites

#: What a junction is drawn in. A feature built in code has no colour of its own, and
#: `liulab_mbio.snapgene` writes SnapGene's default grey for one that has none.
JUNCTION_COLOR = "#ff9900"


@dataclass(frozen=True, slots=True)
class Junction:
    """Where two pieces meet in the product, and the recognition site the join puts back.

    Parameters
    ----------
    start
        0-based index in the product of the first base the two overhangs paired on.
    overhang
        What they spell, written on the top strand, and ``""`` for a blunt join, which pairs on
        no base at all.
    before, after
        The pieces either side, named as they were cut.
    spells
        The bases the product reads across the join, or ``""`` where the join puts back no
        recognition site at all.
    enzyme
        Whose site that is, or ``""`` where it spells none.
    site
        Where those bases lie in the product, or ``None`` where the join spells no site.
    """

    start: int
    overhang: str
    before: str
    after: str
    spells: str = ""
    enzyme: str = ""
    site: Segment | None = None

    @property
    def end(self) -> int:
        """Where the paired bases end."""
        return self.start + len(self.overhang)

    @property
    def span(self) -> Segment | None:
        """Those bases as a segment of the product, and ``None`` where a blunt join paired on none."""
        return Segment(self.start, self.end) if self.overhang else None

    @property
    def marked(self) -> Segment:
        """The bases a map draws for it: the site it spells, or the bases the overhangs paired on.

        A blunt join spelling no site has neither, and is marked at the base it begins on.
        """
        if self.site is not None:
            return self.site
        return self.span or Segment(self.start, self.start + 1)

    @property
    def label(self) -> str:
        """What the junction is called: the site it spells, or the overhang where it spells none."""
        return self.spells or self.overhang or "blunt"


@dataclass(frozen=True, slots=True)
class Ligation:
    """What one ligation makes.

    Parameters
    ----------
    product
        The circular plasmid: every piece's features carried to their new coordinates, the
        primers annotated where they anneal, and each junction marked.
    pieces
        The pieces that went in, the backbone first.
    junctions
        Where they meet, in the product's own order.
    """

    product: SequenceRecord
    pieces: tuple[Piece, ...]
    junctions: tuple[Junction, ...]

    @property
    def junction_positions(self) -> tuple[int, ...]:
        """Where each junction begins, which is what a validation design reads across."""
        return tuple(one.start for one in self.junctions)

    @property
    def blunt(self) -> bool:
        """Whether any junction is blunt, which is what decides the ligation's own incubation."""
        return any(not one.overhang for one in self.junctions)

    @property
    def checks(self) -> tuple[Check, ...]:
        """Judge the product, as data a protocol can print.

        One check per piece counting the whole copies of it the product holds, then how many
        junctions read across the join what they claim. A piece is counted from what the digest
        released and not from the product's own construction, so the count answers for the
        ligation rather than repeating it.
        """
        checks = [
            Check(
                piece.name,
                "pass" if copies == 1 else "fail",
                copies,
                f"{piece.length} bases, whole and once"
                if copies == 1
                else f"{copies} whole copies of its {piece.length} bases",
            )
            for piece in self.pieces
            for copies in (_copies(self.product, piece.bases),)
        ]
        matched = sum(_spells(self.product, one) for one in self.junctions)
        checks.append(
            Check(
                "junctions",
                "pass" if matched == len(self.junctions) else "fail",
                matched,
                ", ".join(f"{one.label} at {one.start}" for one in self.junctions),
            )
        )
        return tuple(checks)

    @property
    def status(self) -> Status:
        """The worst status of any check."""
        return worst_of(self.checks)


def ligate(backbone: Piece, insert: Piece, *, name: str = "") -> Ligation:
    """Ligate `insert` into `backbone` and return the circular product it makes.

    The backbone's right end meets the insert's left, and the insert's right end closes the
    circle. The product is turned so the vector's own first base keeps the place it had, which
    leaves the vector's coordinates readable and keeps a junction off base zero.

    Raises
    ------
    ValueError
        If the two pieces' ends do not anneal into a circle.
    """
    if not closes(backbone, insert):
        raise ValueError(
            f"{backbone.name} and {insert.name} do not ligate into a circle; their ends were "
            "not checked before they were joined"
        )
    pieces = (backbone, insert)
    bases = ""
    features: list[Feature] = []
    primers: list[Primer] = []
    joins: list[tuple[int, Piece, Piece]] = []
    for place, piece in enumerate(pieces):
        at = len(bases)
        bases += piece.bases
        over, kept = carried(piece.source, piece.start, piece.end, offset=at - piece.start)
        features.extend(over)
        primers.extend(kept)
        joins.append((at, pieces[place - 1], piece))
    origin = _origin(backbone)
    product = SequenceRecord(
        bases, topology="circular", name=name, features=tuple(features), primers=tuple(primers)
    )
    turned = rotate(product, origin) if origin else product
    junctions = tuple(
        sorted(
            (
                _junction((at - origin) % len(bases), before, after, turned)
                for at, before, after in joins
            ),
            key=lambda junction: junction.start,
        )
    )
    marked = dataclasses.replace(
        turned, features=(*turned.features, *(_junction_feature(one) for one in junctions))
    )
    return Ligation(ordered(marked), pieces, junctions)


def _spells(product: SequenceRecord, junction: Junction) -> bool:
    """Whether the product reads across the join what the two overhangs paired on.

    A blunt join pairs on no base at all, so there is nothing to read and nothing to contradict.
    """
    return junction.span is None or product.extract(junction.span) == junction.overhang


def _junction(seam: int, before: Piece, after: Piece, product: SequenceRecord) -> Junction:
    """Build the junction whose pieces meet at `seam`, reading off the product what it put back."""
    overhang = after.fragment.left_overhang
    start = _paired(seam, after.left_enzyme, len(product))
    site = _restored(product, start, len(overhang), (before.right_enzyme, after.left_enzyme))
    if site is None:
        return Junction(start, overhang, before.name, after.name)
    return Junction(
        start,
        overhang,
        before.name,
        after.name,
        product.extract(site.span),
        site.enzyme.name,
        site.span,
    )


def _paired(seam: int, enzyme: Enzyme, total: int) -> int:
    """Return where the bases the two overhangs paired on begin, given the pieces meet at `seam`.

    A fragment is bounded by the top-strand cuts and writes both overhangs top-strand, so which
    strand overhangs decides where those bases lie. A 5' overhang is the following piece's own
    first bases and begins at the seam; a 3' overhang is the last bases of the piece before, so
    it begins that many earlier. A blunt end leaves none and begins at the seam either way.
    """
    if enzyme.end != "3'":
        return seam
    return (seam - enzyme.overhang_length) % total


def _restored(
    product: SequenceRecord, start: int, length: int, enzymes: tuple[Enzyme, Enzyme]
) -> CutSite | None:
    """Return the site the product reads across the join at `start`, or ``None`` where it reads none.

    Both ends' enzymes are looked for: a join between two ends of one enzyme puts that site
    back, and one between two compatible ends of different enzymes may put back either or
    neither. A blunt join pairs on no base at all, so the join has to fall inside the site
    rather than under it.
    """
    total = len(product)
    covered = {(start + step) % total for step in range(length)}
    for site in find_sites(product, enzymes):
        reach = len(site.enzyme.site)
        span = {(site.start + step) % total for step in range(reach)}
        if covered <= span and (covered or 0 < (start - site.start) % total < reach):
            return site
    return None


def _junction_feature(junction: Junction) -> Feature:
    """Draw the junction, named for the site it spells."""
    note = f"{junction.before} to {junction.after}"
    if junction.enzyme:
        note += f", {junction.enzyme} site restored"
    return Feature(
        f"{junction.label} junction",
        "misc_feature",
        (junction.marked,),
        color=JUNCTION_COLOR,
        qualifiers={"note": (note,)},
    )


def _origin(piece: Piece) -> int:
    """Where the piece's own source origin falls in the product, or 0 when it is not there."""
    turns = (0, len(piece.source)) if piece.source.topology == "circular" else (0,)
    for turn in turns:
        if piece.start <= turn < piece.end:
            return turn - piece.start
    return 0


def _copies(record: SequenceRecord, bases: str) -> int:
    """Count where `bases` reads whole in `record`, on either strand and across the origin."""
    if not bases or len(bases) > len(record):
        return 0
    haystack = record.sequence
    if record.topology == "circular":
        haystack += record.sequence[: len(bases) - 1]
    return _occurrences(haystack, bases) + _occurrences(haystack, reverse_complement(bases))


def _occurrences(haystack: str, needle: str) -> int:
    """Count where `needle` reads in `haystack`, overlapping copies counted separately."""
    found, at = 0, haystack.find(needle)
    while at >= 0:
        found, at = found + 1, haystack.find(needle, at + 1)
    return found

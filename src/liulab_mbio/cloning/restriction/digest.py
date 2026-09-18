"""Cut a record with the enzymes a restriction and ligation cloning names, and name the pieces.

`liulab_mbio.sites.digest` leaves fragments; what a plan needs besides their bases is which
enzyme cut each end. The overhang alone does not say whether two ends anneal -- a 5' overhang
and a 3' overhang spelling the same bases run the wrong way for each other -- and only the
enzyme record carries the end type. So a `Piece` is a fragment with its two enzymes, and
`liulab_mbio.overhangs.compatible` is the rule every join here is held to.

This method's digest makes at most two cuts in the record the insert comes out of, and the piece
between them is what goes in: one site of each enzyme, or two sites of one. A vector is held to
less, because a further site there only shortens what the gel drops -- until it drops more than
it keeps, which is when the longest piece is no longer a backbone. Anything else is refused,
naming the sites that refused it. `diagnostic` is the one digest held to no such rule: it cuts a
finished miniprep to say whether the insert is in it, so it reads whatever bands the record
gives.

Two ends this method leaves may anneal to each other. `self_closing` says the backbone can shut
with nothing in it, and `reversible` says the insert goes in either way round; neither is refused
here, because a dephosphorylation answers the first and a colony PCR reading out of the insert
answers the second.

Coordinates are the model's, 0-based and half-open, and a piece across the origin of a circular
record ends past the record's length.
"""

import dataclasses
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from liulab_mbio.bench.steps import listed
from liulab_mbio.edits import flipped
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.overhangs import End, compatible
from liulab_mbio.sequence import Segment, SequenceRecord, reverse_complement
from liulab_mbio.sites import CutSite, EnzymeLike, Fragment, digest, find_sites

#: How many enzymes one digest of this method uses. Two of them cut the insert out
#: directionally; one cuts both of its ends, and then the insert goes in either way round.
MAX_ENZYMES = 2

#: How many cuts one digest makes. A piece comes out between two of them, whether two enzymes
#: make one cut each or one enzyme makes both; a third cut falls inside that piece.
MAX_CUTS = 2


@dataclass(frozen=True, slots=True)
class Piece:
    """One fragment a digest released, and the enzyme that cut each of its ends.

    Parameters
    ----------
    name
        What its tube and its gel lane are labelled.
    source
        The record it was cut out of.
    fragment
        Where it lay in that record, and the overhang at each end.
    left_enzyme, right_enzyme
        The enzyme that made each cut, in the same order the overhangs are written.
    """

    name: str
    source: SequenceRecord
    fragment: Fragment
    left_enzyme: Enzyme
    right_enzyme: Enzyme

    @property
    def start(self) -> int:
        """Where it begins in `source`."""
        return self.fragment.start

    @property
    def end(self) -> int:
        """Where it ends, passing the source's length when it runs across the origin."""
        return self.fragment.end

    @property
    def length(self) -> int:
        """Bases of top strand, which is what a gel measures."""
        return self.fragment.length

    @property
    def bases(self) -> str:
        """The top strand it carries, 5' to 3'."""
        return self.source.extract(Segment(self.start, self.end))

    @property
    def left_end(self) -> End:
        """Its left end, as the rule joining two ends sees it."""
        return End.cut_by(self.left_enzyme, self.fragment.left_overhang)

    @property
    def right_end(self) -> End:
        """Its right end."""
        return End.cut_by(self.right_enzyme, self.fragment.right_overhang)


def cut(
    record: SequenceRecord, enzymes: Sequence[Enzyme], *, unique: bool = True
) -> tuple[Piece, ...]:
    """Digest `record` with every one of `enzymes` and return the pieces, in top-strand order.

    Every enzyme must cut, and where `unique` they cut no more than `MAX_CUTS` times between
    them: this method cuts a piece out between two cuts, so a further cut falls inside that piece
    and cuts it in two. One site of each enzyme and two sites of one are the two ways to make
    those cuts. `unique` set aside drops that to one site at least, which is what a vector is
    held to.

    A linear record also has the two ends it came with. No enzyme made them and nothing ligates
    to them, so the pieces carrying them are left out and what is returned is what the digest
    released.

    Parameters
    ----------
    record
        The plasmid or fragment to cut.
    enzymes
        One or two enzymes, which go in the same reaction.
    unique
        Refuse an enzyme reading more than one site. `opened` is what sets it aside.

    Raises
    ------
    ValueError
        If no enzyme is given or more than `MAX_ENZYMES`, if an enzyme reads no site or, where
        `unique`, the enzymes read more than `MAX_CUTS` between them, or if no piece has an
        enzyme at both of its ends.
    """
    if not 1 <= len(enzymes) <= MAX_ENZYMES:
        raise ValueError(
            f"this method digests with one or two enzymes, got {len(enzymes)}: "
            f"{', '.join(one.name for one in enzymes) or 'none'}"
        )
    what = record.name or "the record"
    _sites_each(record, enzymes, what, unique=unique)
    length = len(record)
    at: dict[int, Enzyme] = {}
    for site in find_sites(record, enzymes):
        if site.cuts:
            at.setdefault(site.top_cut % length, site.enzyme)
    released: list[tuple[Fragment, Enzyme, Enzyme]] = []
    for fragment in digest(record, enzymes):
        left, right = at.get(fragment.start % length), at.get(fragment.end % length)
        if left is not None and right is not None:
            released.append((fragment, left, right))
    if not released:
        raise ValueError(
            f"the digest of {what} releases nothing with an enzyme at both ends; a linear "
            "fragment is cut out between one site of each enzyme, so it needs both of them"
        )
    labels = _labels(f"{record.name} fragment".strip(), len(released))
    return tuple(
        Piece(label, record, fragment, left, right)
        for label, (fragment, left, right) in zip(labels, released, strict=True)
    )


def opened(vector: SequenceRecord, enzymes: Sequence[Enzyme]) -> tuple[Piece, ...]:
    """Return what `enzymes` leave of `vector`, its backbone first.

    The backbone is the longest piece the digest gives; the rest is what the cloning site
    releases, and a gel is how the backbone is separated from it. A further site of either enzyme
    is allowed here and only lengthens what the gel drops, until the drop outweighs the backbone
    -- past that the longest piece is a piece of the vector rather than its backbone. A backbone
    `self_closing` calls true is opened all the same: it closes on itself with no insert, and a
    dephosphorylation is what stops it.

    Raises
    ------
    ValueError
        If the vector is not circular, for any reason `cut` refuses, or if the digest drops as
        much as it keeps.
    """
    if vector.topology != "circular":
        raise ValueError("a vector is cut open and closed again, so it must be circular")
    pieces = cut(vector, enzymes, unique=False)
    backbone = max(pieces, key=lambda piece: piece.length)
    dropped = sum(piece.length for piece in pieces if piece is not backbone)
    if dropped >= backbone.length:
        raise ValueError(
            f"{listed([one.name for one in enzymes])} cut {vector.name or 'the vector'} into "
            f"{len(pieces)} pieces, the longest {backbone.length} bp against {dropped} bp "
            "dropped, so what the gel keeps is not a backbone; a further site of one of them "
            "lies outside the cloning site"
        )
    return _named(vector, backbone, pieces, "backbone")


def excised(source: SequenceRecord, enzymes: Sequence[Enzyme], *, into: Piece) -> tuple[Piece, ...]:
    """Return what `enzymes` cut out of `source`, the insert first and turned to go into `into`.

    The insert is the shortest piece the digest gives. It is written on whichever strand closes
    the circle with `into`: which way round it lay in `source` is how that plasmid was built and
    says nothing about which way it goes in here.

    Raises
    ------
    ValueError
        For any reason `cut` refuses, or if the insert's ends do not anneal to the backbone's on
        either strand, naming the two ends that did not.
    """
    tried: list[tuple[Piece, Piece]] = []
    for record in (source, flipped(source)):
        pieces = cut(record, enzymes)
        insert = min(pieces, key=lambda piece: piece.length)
        if closes(into, insert):
            return _named(source, insert, pieces, "insert")
        tried.append((into, insert))
    raise ValueError(_not_annealing(tried))


@dataclass(frozen=True, slots=True)
class Diagnostic:
    """The digest that says a miniprep carries the insert, and the bands it should give.

    The enzymes are the cloning pair itself: they cut their own sites back out of the product,
    so a correct clone drops the insert and a colony carrying the vector it went into gives what
    that plasmid gives. Nothing else has to be chosen, and no band is a number this package
    invented.

    Parameters
    ----------
    enzymes
        What the miniprep is cut with.
    clone, empty
        The bands a correct clone gives and the bands the vector gives, largest first as a gel
        reads them.
    names
        What each of those two lanes is labelled.
    """

    enzymes: tuple[Enzyme, ...]
    clone: tuple[int, ...]
    empty: tuple[int, ...]
    names: tuple[str, str]

    @property
    def tells_them_apart(self) -> bool:
        """Whether the two lanes differ at all, which is the whole point of running it."""
        return self.clone != self.empty

    @property
    def bands(self) -> tuple[int, ...]:
        """Every band on the gel, which is what a ladder and an agarose percentage are chosen for."""
        return (*self.clone, *self.empty)


def diagnostic(
    product: SequenceRecord, vector: SequenceRecord, enzymes: Sequence[Enzyme]
) -> Diagnostic:
    """Return the diagnostic digest of `product` against the `vector` it went into."""
    return Diagnostic(
        tuple(enzymes),
        _bands(product, enzymes),
        _bands(vector, enzymes),
        (product.name or "clone", vector.name or "vector"),
    )


def _bands(record: SequenceRecord, enzymes: Sequence[Enzyme]) -> tuple[int, ...]:
    """Return the bands one digest gives, largest first; a record nothing cuts runs as one."""
    pieces = digest(record, enzymes)
    return tuple(sorted((piece.length for piece in pieces), reverse=True)) or (len(record),)


def _named(
    record: SequenceRecord, kept: Piece, pieces: Sequence[Piece], what: str
) -> tuple[Piece, ...]:
    """Label the piece that goes on for what it is, and the rest for what a gel drops."""
    rest = [piece for piece in pieces if piece is not kept]
    labels = _labels(f"{record.name} offcut".strip(), len(rest))
    return (
        dataclasses.replace(kept, name=f"{record.name} {what}".strip()),
        *(
            dataclasses.replace(piece, name=label)
            for piece, label in zip(rest, labels, strict=True)
        ),
    )


def closes(backbone: Piece, insert: Piece) -> bool:
    """Whether these two pieces ligate into a circle, the backbone's right end meeting first."""
    return compatible(backbone.right_end, insert.left_end) and compatible(
        insert.right_end, backbone.left_end
    )


def self_closing(piece: Piece) -> bool:
    """Whether a piece's own two ends anneal, so it can close into a circle with nothing in it."""
    return compatible(piece.right_end, piece.left_end)


def reversible(backbone: Piece, insert: Piece) -> bool:
    """Whether `insert` also closes `backbone` turned the other way round.

    One enzyme at both ends, two enzymes leaving compatible ends, and blunt ends all leave an
    insert that goes in either way. Turning a piece over swaps its two ends and writes each
    overhang as the other strand reads it, which is its reverse complement.
    """
    left = End.cut_by(insert.right_enzyme, reverse_complement(insert.fragment.right_overhang))
    right = End.cut_by(insert.left_enzyme, reverse_complement(insert.fragment.left_overhang))
    return compatible(backbone.right_end, left) and compatible(right, backbone.left_end)


def turned(insert: Piece) -> Piece:
    """Return `insert` turned over: the same piece, cut from the other strand of its source.

    A cohesive end's two strands stop at different bases, so which bases turn over and where
    they land both move with the overhang. It is cut with the enzymes that cut its own two ends.
    """
    enzymes = tuple(dict.fromkeys((insert.left_enzyme, insert.right_enzyme)))
    pieces = cut(flipped(insert.source), enzymes)
    return dataclasses.replace(min(pieces, key=lambda piece: piece.length), name=insert.name)


def _not_annealing(tried: Sequence[tuple[Piece, Piece]]) -> str:
    """Say which join refused, on the strand that came closest."""
    backbone, insert = tried[0]
    if compatible(backbone.right_end, insert.left_end):
        one, other, where = insert.right_end, backbone.left_end, "closing the circle"
    else:
        one, other, where = backbone.right_end, insert.left_end, "where the insert goes in"
    return (
        f"the ends do not anneal {where}: {said_end(one)} meets {said_end(other)}, on either strand "
        "of the insert"
    )


def said_end(end: End) -> str:
    """Name one cut end the way a refusal or a verdict reads it."""
    return "a blunt end" if end.end_type == "blunt" else f"a {end.end_type} {end.overhang}"


def said_ends(piece: Piece) -> str:
    """Name a piece's two ends the same way, in one phrase where they read alike."""
    left, right = piece.left_end, piece.right_end
    if left != right:
        return f"{said_end(left)} and {said_end(right)}"
    if left.end_type == "blunt":
        return "two blunt ends"
    return f"two {left.end_type} {left.overhang} ends"


def _sites_each(
    record: SequenceRecord, enzymes: Iterable[Enzyme], what: str, *, unique: bool
) -> None:
    """Refuse an enzyme reading no site in `record`, or the enzymes too many, naming what it read.

    Every enzyme has to cut. Where `unique` they cut no more than `MAX_CUTS` times between them,
    which is one site of each enzyme or two sites of one; a vector is held only to the first
    rule, a further site there lengthening what the gel drops rather than cutting the piece.

    Raises
    ------
    ValueError
        If an enzyme reads no site, or, where `unique`, the enzymes read more than `MAX_CUTS`
        between them, which is reported against the one that read most.
    """
    found = [
        (enzyme.name, [site for site in find_sites(record, enzyme) if site.cuts])
        for enzyme in enzymes
    ]
    for name, sites in found:
        if not sites:
            why = (
                "cuts a piece out between two cuts, so every enzyme named has to make one"
                if unique
                else "opens a plasmid between one site of each enzyme, so it needs one of each"
            )
            raise ValueError(f"{name} cuts {what} 0 time(s); this method {why}")
    if unique and sum(len(sites) for _, sites in found) > MAX_CUTS:
        name, sites = max(found, key=lambda one: len(one[1]))
        raise ValueError(
            f"{name} cuts {what} {len(sites)} time(s){_where(sites)}; this method cuts a piece "
            f"out between {MAX_CUTS} cuts, so a further cut falls inside that piece and cuts it "
            "in two"
        )


def _where(found: Sequence[CutSite]) -> str:
    """Name where those sites are, or nothing at all when there are none."""
    return f", at {', '.join(str(site.start) for site in found)}" if found else ""


def _labels(name: str, count: int) -> tuple[str, ...]:
    """Label each piece, numbered only where there is more than one to tell apart."""
    return (name,) if count == 1 else tuple(f"{name} {number}" for number in range(1, count + 1))


def resolve(enzymes: Iterable[EnzymeLike]) -> tuple[Enzyme, ...]:
    """Read the enzymes by name or by record, keeping the order they were named in.

    Raises
    ------
    KeyError
        If a name is not one the package ships.
    """
    return tuple(get_enzyme(one) if isinstance(one, str) else one for one in enzymes)

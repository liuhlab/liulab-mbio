"""Cut a record with the enzymes a restriction and ligation cloning names, and name the pieces.

`liulab_mbio.sites.digest` leaves fragments; what a plan needs besides their bases is which
enzyme cut each end. The overhang alone does not say whether two ends anneal -- a 5' overhang
and a 3' overhang spelling the same bases run the wrong way for each other -- and only the
enzyme record carries the end type. So a `Piece` is a fragment with its two enzymes, and
`liulab_mbio.overhangs.compatible` is the rule every join here is held to.

This method needs one site of each enzyme in the record it cuts: that is what makes the pieces
the backbone and what it releases. Anything else is refused, naming the sites that refused it.

Coordinates are the model's, 0-based and half-open, and a piece across the origin of a circular
record ends past the record's length.
"""

import dataclasses
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from liulab_mbio.edits import flipped
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.overhangs import End, compatible
from liulab_mbio.sequence import Segment, SequenceRecord
from liulab_mbio.sites import CutSite, EnzymeLike, Fragment, digest, find_sites

#: How many enzymes one digest of this method uses. Two of them cut the insert out
#: directionally; `self_closing` is what refuses a pair that leaves the vector free to rejoin.
MAX_ENZYMES = 2


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


def cut(record: SequenceRecord, enzymes: Sequence[Enzyme]) -> tuple[Piece, ...]:
    """Digest `record` with every one of `enzymes` and return the pieces, in top-strand order.

    Each enzyme must read exactly one site: this method cuts a piece out between one site of
    each, so a second site of either cuts that piece in two.

    Parameters
    ----------
    record
        The circular plasmid to cut.
    enzymes
        One or two enzymes, which go in the same reaction.

    Raises
    ------
    ValueError
        If no enzyme is given or more than `MAX_ENZYMES`, if an enzyme reads anything but one
        site, or if a cut leaves an end no enzyme made.
    """
    if not 1 <= len(enzymes) <= MAX_ENZYMES:
        raise ValueError(
            f"this method digests with one or two enzymes, got {len(enzymes)}: "
            f"{', '.join(one.name for one in enzymes) or 'none'}"
        )
    _one_site_each(record, enzymes, record.name or "the record")
    length = len(record)
    at: dict[int, Enzyme] = {}
    for site in find_sites(record, enzymes):
        if site.cuts:
            at.setdefault(site.top_cut % length, site.enzyme)
    pieces = digest(record, enzymes)
    labels = _labels(f"{record.name} fragment".strip(), len(pieces))
    return tuple(
        Piece(
            label, record, fragment, _at(at, fragment.start, length), _at(at, fragment.end, length)
        )
        for label, fragment in zip(labels, pieces, strict=True)
    )


def opened(vector: SequenceRecord, enzymes: Sequence[Enzyme]) -> tuple[Piece, ...]:
    """Return what `enzymes` leave of `vector`, its backbone first.

    The backbone is the longest piece the digest gives; the rest is what the cloning site
    releases, and a gel is how the backbone is separated from it.

    Raises
    ------
    ValueError
        If the vector is not circular, for any reason `cut` refuses, or if the backbone's two
        ends anneal to each other, which leaves it free to close on itself with no insert.
    """
    if vector.topology != "circular":
        raise ValueError("a vector is cut open and closed again, so it must be circular")
    pieces = cut(vector, enzymes)
    backbone = max(pieces, key=lambda piece: piece.length)
    if self_closing(backbone):
        raise ValueError(
            f"the backbone's two ends anneal to each other ({_said(backbone.left_end)} and "
            f"{_said(backbone.right_end)}), so it closes on itself and the insert could go in "
            "either way round; name two enzymes leaving ends that do not"
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


def _not_annealing(tried: Sequence[tuple[Piece, Piece]]) -> str:
    """Say which join refused, on the strand that came closest."""
    backbone, insert = tried[0]
    if compatible(backbone.right_end, insert.left_end):
        one, other, where = insert.right_end, backbone.left_end, "closing the circle"
    else:
        one, other, where = backbone.right_end, insert.left_end, "where the insert goes in"
    return (
        f"the ends do not anneal {where}: {_said(one)} meets {_said(other)}, on either strand "
        "of the insert"
    )


def _said(end: End) -> str:
    """Name one end the way a refusal reads it."""
    return "a blunt end" if end.end_type == "blunt" else f"a {end.end_type} {end.overhang}"


def _one_site_each(record: SequenceRecord, enzymes: Iterable[Enzyme], what: str) -> None:
    """Refuse an enzyme reading anything but one site in `record`, naming the sites it reads.

    Raises
    ------
    ValueError
        If an enzyme reads no site, or more than one.
    """
    for enzyme in enzymes:
        found = [site for site in find_sites(record, enzyme) if site.cuts]
        if len(found) != 1:
            raise ValueError(
                f"{enzyme.name} cuts {what} {len(found)} time(s){_where(found)}; this method "
                "cuts a piece out between one site of each enzyme, so a second site of either "
                "cuts that piece in two"
            )


def _where(found: Sequence[CutSite]) -> str:
    """Name where those sites are, or nothing at all when there are none."""
    return f", at {', '.join(str(site.start) for site in found)}" if found else ""


def _at(found: dict[int, Enzyme], position: int, length: int) -> Enzyme:
    """Return the enzyme that cut at `position`.

    Raises
    ------
    ValueError
        If no enzyme cut there, which is the end of a linear record rather than a cut.
    """
    enzyme = found.get(position % length)
    if enzyme is None:
        raise ValueError(f"the end at {position} is not a cut, so nothing can be ligated to it")
    return enzyme


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

"""Take the user's own destination vector, and make it one where it is not.

A vector is a destination when the scheme's internal enzyme excises one piece from it, leaving
position one's entry overhang at one end and the cloning scar at the other. That is what a round
opens, so it is read off a digest rather than off a string: a vector already carrying an internal
stuffer is taken as it stands, and one that does not has the terminal position's stuffer put at a
site the user names. The stuffer a vector carries is not a special shape — it is a position's own,
the terminal one's, whose prefix ends with the overhang position one enters on.

An insertion is reported as an `liulab_mbio.edits.EditReport`, so the user reads what it changed
rather than trusting a new file. Coordinates are the model's, and a stuffer across the origin ends
past the record's length.
"""

from dataclasses import dataclass

from liulab_mbio.edits import EditReport, insert
from liulab_mbio.library.scheme import Scheme
from liulab_mbio.sequence import Feature, Segment, SequenceRecord
from liulab_mbio.sites import digest, find_sites

#: Where the internal stuffer goes: the name of a feature, or a ``(start, end)`` span.
type Site = str | tuple[int, int]


@dataclass(frozen=True, slots=True)
class Destination:
    """A vector a round can open, and what making it one changed.

    Parameters
    ----------
    record
        The vector. The record it was given, unchanged, where that already carried a stuffer.
    stuffer
        The piece the internal enzyme excises to open the vector, bounded by its two cuts. It
        ends past the record's length when it runs across the origin.
    edit
        What inserting the stuffer did to the features and primer binding sites it did not simply
        shift, or ``None`` where the vector already carried one.
    """

    record: SequenceRecord
    stuffer: Segment
    edit: EditReport | None = None


def destination_vector(
    vector: SequenceRecord, scheme: Scheme, *, site: Site | None = None
) -> Destination:
    """Return `vector` as a destination the scheme's first round can open.

    A vector whose internal enzyme already excises one piece leaving the scheme's two overhangs
    is returned as it stands, `site` unread. One that does not has the terminal position's
    internal stuffer put at the start of `site`, and the result is held to the same rule.

    Parameters
    ----------
    vector
        The user's own vector, circular.
    scheme
        The architecture the build is given, which says what opens the vector and on what.
    site
        Where to put a stuffer, as a feature name or a ``(start, end)`` span. Only read where one
        has to be put.

    Raises
    ------
    ValueError
        If any of the scheme's enzymes reads a site outside the stuffer, if a stuffer has to be
        put and none is named, if `site` names no feature or does not fit, if the stuffer is not
        a whole number of codons where it would land inside a coding sequence, or if the scheme's
        own stuffer does not carry both of the cuts that open a vector.
    KeyError
        If the scheme names an enzyme this package does not ship.
    """
    found = _stuffer(vector, scheme)
    if found is not None:
        _check_clean(vector, scheme, found)
        return Destination(vector, found)
    block = scheme.internal_stuffer(-1)
    _check_block(scheme, block)
    at = _position(vector, scheme, site)
    _check_frame(vector, at, block)
    edited, report = insert(vector, at, block)
    made = _stuffer(edited, scheme)
    if made is None:
        raise ValueError(
            f"the stuffer put at {at} left a vector {scheme.internal.name} does not open on "
            f"{scheme.entry_overhang(0)!r} and {scheme.scar_overhang!r}: the bases it now joins "
            "spell a further site. Name another site"
        )
    _check_clean(edited, scheme, made)
    return Destination(edited, made, report)


def _stuffer(record: SequenceRecord, scheme: Scheme) -> Segment | None:
    """Return the piece the internal enzyme excises to open `record`, or ``None`` where none is.

    The piece is the one whose two ends are the overhangs a first-round part enters and leaves on,
    which is what makes the vector a destination rather than a plasmid with a stuffer-shaped gap.
    """
    entry, scar = scheme.entry_overhang(0), scheme.scar_overhang
    for piece in digest(record, scheme.internal):
        if piece.left_overhang == entry and piece.right_overhang == scar:
            return Segment(piece.start, piece.end)
    return None


def _check_block(scheme: Scheme, block: str) -> None:
    """Refuse a scheme whose own stuffer does not carry both of the cuts that open a vector."""
    entry, scar = scheme.entry_overhang(0), scheme.scar_overhang
    left = {site.overhang for site in find_sites(SequenceRecord(block), scheme.internal)}
    if missing := [wanted for wanted in (entry, scar) if wanted not in left]:
        spells = ", ".join(sorted(one for one in left if one)) or "no overhang at all"
        raise ValueError(
            f"the internal stuffer of position {scheme.positions[-1].name!r} leaves {spells} "
            f"where {scheme.internal.name} cuts it, and not "
            f"{' and '.join(repr(one) for one in missing)}: the stuffer a vector carries holds "
            "both of the cuts that open it"
        )


def _position(vector: SequenceRecord, scheme: Scheme, site: Site | None) -> int:
    """Return where the stuffer goes, from the feature or the span the user names."""
    if site is None:
        raise ValueError(
            f"this vector carries no internal stuffer for {scheme.internal.name} to excise: name "
            "the site to put one at, as a feature name or a (start, end) span"
        )
    if isinstance(site, tuple):
        start, end = site
        if not 0 <= start < end <= len(vector):
            raise ValueError(
                f"the site {start}-{end} does not lie inside {len(vector)} bases of vector"
            )
        return start
    found = next(
        (feature for feature in vector.features if feature.name.lower() == site.lower()), None
    )
    if found is None:
        raise ValueError(f"this vector annotates no feature called {site!r}")
    return found.segments[0].start


def _check_frame(vector: SequenceRecord, at: int, block: str) -> None:
    """Refuse a stuffer that is not whole codons where it would land inside a coding sequence."""
    if len(block) % 3 == 0:
        return
    coding = _coding(vector, at)
    if coding is not None:
        raise ValueError(
            f"the internal stuffer is {len(block)} bases, which is not a whole number of codons, "
            f"and {at} lies inside the coding sequence {coding.name!r}: putting it there would "
            "shift the reading frame. Name a site outside that coding sequence"
        )


def _coding(record: SequenceRecord, at: int) -> Feature | None:
    """Return the first coding sequence `at` lies inside, or ``None`` where it lies in none.

    An insertion at a segment's first base lands ahead of it, which shifts the whole feature and
    leaves its frame alone, so that position counts as outside.
    """
    length = len(record)
    for feature in record.features:
        if feature.type == "CDS" and any(
            0 < (at - segment.start) % length < segment.end - segment.start
            for segment in feature.segments
        ):
            return feature
    return None


def _check_clean(record: SequenceRecord, scheme: Scheme, stuffer: Segment) -> None:
    """Refuse a vector reading any of the scheme's enzymes outside its internal stuffer.

    Inside it they are the design: the internal enzyme's two cuts and whichever blunt chopper
    shreds the excised piece. Outside, a round would cut the backbone.
    """
    length = len(record)
    inside = {index % length for index in range(stuffer.start, stuffer.end)}
    for site in find_sites(record, (scheme.internal, scheme.external, *scheme.blunt)):
        if not all((site.start + step) % length in inside for step in range(len(site.enzyme.site))):
            raise ValueError(
                f"{site.enzyme.name} reads a site at {site.start} on the "
                f"{site.strand.name.lower()} strand, outside the internal stuffer at "
                f"{stuffer.start}-{stuffer.end}: a round would cut the backbone there. Take that "
                "site out of the vector, or name a scheme whose enzymes it is free of"
            )

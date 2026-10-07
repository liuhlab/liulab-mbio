"""Take the user's own destination vector, and make it one where it is not.

A vector is a destination when the method's internal enzyme excises one piece from it, leaving
the entry overhang at one end and the cloning scar at the other. That is what a round opens, so
it is read off a digest rather than off a string: a vector already carrying an internal stuffer
is taken as it stands, and one that does not has one put at a site the user names. The stuffer a
vector carries is not a special shape — it is the method's own, whose prefix ends with the
overhang a part enters on.

An insertion is reported as an `liulab_mbio.edits.EditReport`, so the user reads what it changed
rather than trusting a new file. Coordinates are the model's, and a stuffer across the origin ends
past the record's length.
"""

from collections.abc import Iterable
from dataclasses import dataclass

from liulab_mbio.checks import Check
from liulab_mbio.codons import CodonUsage
from liulab_mbio.edits import EditReport, insert
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.sequence import Feature, Segment, SequenceRecord
from liulab_mbio.sites import (
    OLIGO_REACH,
    CutSite,
    DomesticationReport,
    EnzymeSearch,
    digest,
    domesticate,
    find_sites,
    free_enzyme_search,
    out_of_reach,
)
from liulab_synbio.library.method import IGGA, INTERFACE_OVERHANGS, Scheme, refuse

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
    """Return `vector` as a destination the method's first round can open.

    A vector whose internal enzyme already excises one piece leaving the method's two overhangs
    is returned as it stands, `site` unread. One that does not has the method's
    internal stuffer put at the start of `site`, and the result is held to the same rule.

    Parameters
    ----------
    vector
        The user's own vector, circular.
    scheme
        The method the build is given, which says what opens the vector and on what.
    site
        Where to put a stuffer, as a feature name or a ``(start, end)`` span. Only read where one
        has to be put.

    Raises
    ------
    ValueError
        If any of the method's enzymes reads a site outside the stuffer, if a stuffer has to be
        put and none is named, if `site` names no feature or does not fit, or if the stuffer is
        not a whole number of codons where it would land inside a coding sequence.
    KeyError
        If the method names an enzyme this package does not ship.
    """
    found = _stuffer(vector, scheme)
    if found is not None:
        _check_clean(vector, scheme, found)
        return Destination(vector, found)
    block = scheme.internal_stuffer
    at = _position(vector, scheme, site)
    _check_frame(vector, at, block)
    edited, report = insert(vector, at, block)
    made = _stuffer(edited, scheme)
    if made is None:
        raise ValueError(
            f"the stuffer put at {at} left a vector {scheme.internal.name} does not open on "
            f"{scheme.entry_overhang!r} and {scheme.scar_overhang!r}: the bases it now joins "
            "spell a further site. Name another site"
        )
    _check_clean(edited, scheme, made)
    return Destination(edited, made, report)


def _stuffer(record: SequenceRecord, scheme: Scheme) -> Segment | None:
    """Return the piece the internal enzyme excises to open `record`, or ``None`` where none is.

    The piece is the one whose two ends are the overhangs a first-round part enters and leaves on,
    which is what makes the vector a destination rather than a plasmid with a stuffer-shaped gap.
    """
    entry, scar = scheme.entry_overhang, scheme.scar_overhang
    for piece in digest(record, scheme.internal):
        if piece.left_overhang == entry and piece.right_overhang == scar:
            return Segment(piece.start, piece.end)
    return None


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
    leaves its frame alone, so each segment is asked about past that base.
    """
    for feature in record.features:
        if feature.type == "CDS" and any(
            segment.end - segment.start > 1
            and record.covers(Segment(segment.start + 1, segment.end), at)
            for segment in feature.segments
        ):
            return feature
    return None


def _check_clean(record: SequenceRecord, scheme: Scheme, stuffer: Segment) -> None:
    """Refuse a vector reading any of the method's enzymes outside its internal stuffer.

    Inside it they are the design: the internal enzyme's two cuts and whichever blunt chopper
    shreds the excised piece. Outside, a round would cut the backbone.
    """
    for site in find_sites(record, (scheme.internal, scheme.external, *scheme.blunt)):
        if not record.covers(stuffer, site.span):
            raise ValueError(
                f"{site.enzyme.name} reads a site at {site.start} on the "
                f"{site.strand.name.lower()} strand, outside the internal stuffer at "
                f"{stuffer.start}-{stuffer.end}: a round would cut the backbone there. Take that "
                "site out of the vector, or name a method whose enzymes it is free of"
            )


def cargo_candidates(scheme: Scheme = IGGA) -> tuple[str, ...]:
    """Return the enzymes a cargo enzyme is searched among, in the order the search reads them.

    The Type IIS enzymes `liulab_mbio`'s Golden Gate pipeline ships a reaction and cycling
    protocol for, less any reading a reserved enzyme's site and less its own last resort, which
    ships no protocol. Derived rather than listed, so an enzyme the package stops shipping a
    protocol for leaves this list with it.

    The list is read from the cloning pipeline that ships those protocols, inside this function
    so that importing a library plan does not load a cloning pipeline with it.
    """
    from liulab_mbio.cloning.goldengate.design import GOLDEN_GATE_ENZYMES, LAST_RESORT

    barred = {enzyme.site for enzyme in scheme.reserved_enzymes}
    return tuple(
        name
        for name in GOLDEN_GATE_ENZYMES
        if name not in LAST_RESORT and get_enzyme(name).site not in barred
    )


@dataclass(frozen=True, slots=True)
class Cargo:
    """Which enzyme may admit cargo to a working vector, and what ruled the others out.

    Parameters
    ----------
    enzyme
        The enzyme, or ``None`` where none of the candidates is free. That is a finding: `search`
        then names the sites blocking each one, which is what a staged digest, a bridging
        assembly or a re-tailoring PCR is reached for against.
    search
        What `liulab_mbio.sites.free_enzyme_search` returned.
    check
        The verdict, in the method's words.
    """

    enzyme: Enzyme | None
    search: EnzymeSearch
    check: Check

    @property
    def sites(self) -> tuple[CutSite, ...]:
        """Every site that blocked a candidate, so a failure can be drawn on the record."""
        return tuple(site for one in self.search.blocked for site in one.sites)


def cargo_enzyme(
    records: Iterable[SequenceRecord],
    *,
    scheme: Scheme = IGGA,
    candidates: Iterable[str] | None = None,
) -> Cargo:
    """Choose the enzyme that admits cargo to a working vector, over the molecules sharing its pot.

    The method fixes every other enzyme, because every other enzyme is carried by DNA already on
    the shelf. This one is searched per vector, against the rule
    `docs/synthesis-and-assembly.md` states: distinct from the enzyme that inserts the working
    cassette, a 4 nt 5' overhang so the method's four junction overhangs are unchanged, and no
    site in any molecule sharing its reaction.

    Parameters
    ----------
    records
        The working vector and every other molecule in the final-assembly pot: the cargo, the
        iGGA cargo and any part of the working cassette.
    scheme
        The method, whose reserved enzymes the candidates exclude.
    candidates
        The enzymes to search among, defaulting to `cargo_candidates`.
    """
    overhangs = {len(one) for pair in INTERFACE_OVERHANGS.values() for one in pair}
    if len(overhangs) != 1:
        refuse("interface-overhangs", f"the method's junctions are not all one length: {overhangs}")
    search = free_enzyme_search(
        records,
        cargo_candidates(scheme) if candidates is None else candidates,
        overhang_length=overhangs.pop(),
        end="5'",
    )
    chosen = search.best
    if chosen is None:
        blocked = "; ".join(f"{one.enzyme.name} {one.reason}" for one in search.blocked)
        detail = (
            "no candidate is free to admit cargo to this working vector, so the final assembly "
            f"cannot be one pot: {blocked}"
        )
    else:
        detail = f"{chosen.name} reads no site in any molecule sharing the final-assembly pot"
    return Cargo(
        chosen,
        search,
        Check("cargo enzyme", "pass" if chosen else "fail", len(search.free), detail),
    )


@dataclass(frozen=True, slots=True)
class Domesticated:
    """A vector with the method's sites taken out of it, and what could not be taken out.

    Parameters
    ----------
    record
        The vector as the changes left it.
    report
        What `liulab_mbio.sites.domesticate` changed, and the sites it left for someone to
        decide about.
    unreachable
        The sites no oligo reaches, refused before anything else was tried.
    check
        The verdict: a pass only where no site of these enzymes is left.
    """

    record: SequenceRecord
    report: DomesticationReport
    unreachable: tuple[CutSite, ...]
    check: Check

    @property
    def remaining(self) -> tuple[CutSite, ...]:
        """Every site still in the record, in top-strand order."""
        return tuple(
            sorted(
                (*self.report.outside_cds, *self.report.unchanged),
                key=lambda site: site.start,
            )
        )


def domesticate_vector(
    vector: SequenceRecord,
    enzymes: Iterable[Enzyme],
    *,
    usage: CodonUsage | None = None,
    reach: int = OLIGO_REACH,
) -> Domesticated:
    """Take the sites `enzymes` read out of `vector`, refusing the ones no oligo reaches.

    The order is the one `docs/research/domestication-methods.md` sets: count repeats first,
    then sites. A site carrying more than `reach` bases of identical sequence on both sides is
    named and left alone, because no oligo placed there — mutagenic primer, assembly overlap or
    Type IIS overhang — is unique to one copy of the repeat. Such a site is designed around by
    changing the enzyme, not edited at the bench.

    A site in no coding sequence is left alone too, by `liulab_mbio.sites.domesticate`'s own
    rule: changing it changes what the vector spells, which no method decides on its own.
    """
    held = tuple(enzymes)
    refused = out_of_reach(vector, find_sites(vector, held), reach=reach)
    record, report = domesticate(vector, held, usage=usage, avoid=held)
    left = tuple(sorted((*report.outside_cds, *report.unchanged), key=lambda one: one.start))
    if not left:
        detail = f"{len(report.changes)} site(s) changed; none of {_named(held)} is left"
    else:
        named = ", ".join(f"{one.enzyme.name} at {one.start}" for one in left)
        blocked = {one.start for one in refused}
        inside = ", ".join(str(one.start) for one in left if one.start in blocked)
        detail = f"{len(report.changes)} site(s) changed, {len(left)} left: {named}" + (
            f". No oligo reaches {inside}: each sits inside a repeat" if inside else ""
        )
    return Domesticated(
        record,
        report,
        refused,
        Check("vector domestication", "pass" if not left else "fail", len(left), detail),
    )


def _named(enzymes: Iterable[Enzyme]) -> str:
    """Return the enzymes, named in the order they were given."""
    return ", ".join(one.name for one in enzymes)

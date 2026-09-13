"""Find restriction sites in a record, put one there, and take one out of a coding sequence.

Two rules decide what counts as a site, and both matter to whoever asks whether an enzyme is
free to use:

- **A site is reported where the enzyme may cut.** Every position of the match needs a base the
  recognition site admits, so an IUPAC code in the template is read as the bases it stands for
  and a site it could spell is reported. `CutSite.certain` is ``False`` for such a hit. Erring
  towards reporting is the safe direction here: an enzyme called free of sites is one a design
  goes on to rely on.
- **A palindromic site is one site, not two.** Its recognition site reads the same on either
  strand, so only the forward match is reported.

Coordinates are the model's, 0-based and half-open, and a site across the origin of a circular
record ends past the record's length.
"""

import re
from collections.abc import Callable, Iterable, Mapping
from dataclasses import KW_ONLY, dataclass
from functools import cache
from itertools import islice, product

from liulab_mbio.edits import EditReport, insert
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.enzymes import enzymes as shipped
from liulab_mbio.sequence import Segment, SequenceRecord, Strand, reverse_complement

#: How many bases NEB recommends 5' of a recognition site for an enzyme to cut near an end.
FLANK_LENGTH = 6

#: Dcm methylates the inner cytosine of this site. Only an enzyme whose supplier says Dcm
#: affects it is held to the rule; `Enzyme.methylation` carries that answer.
DCM_SITE = "CCWGG"

#: How many candidate flanks or fillers to try before giving up on an overhang.
_TRIES = 4096

_RUN = re.compile(r"(.)\1{3}")

#: What each IUPAC code stands for. `liulab_mbio.sequence.IUPAC_DNA` names the same fifteen.
_BASES: Mapping[str, frozenset[str]] = {
    code: frozenset(bases)
    for code, bases in {
        "A": "A", "C": "C", "G": "G", "T": "T",
        "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
        "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
    }.items()
}  # fmt: skip

#: An enzyme, or the name of one the package ships.
type EnzymeLike = Enzyme | str


@dataclass(frozen=True, slots=True)
class CutSite:
    """One recognition site found in a record, and where its enzyme cuts there.

    Parameters
    ----------
    enzyme
        The enzyme that reads this site.
    start
        0-based index of the first base of the match on the TOP strand, whichever strand the
        site lies on.
    strand
        `Strand.FORWARD` when the top strand spells the recognition site, `Strand.REVERSE` when
        the bottom strand does. A palindromic site is always forward.
    top_cut, bottom_cut
        Where the two strands are severed, in the record's own coordinates and in that order
        whatever the strand. Reduced modulo the length on a circular record; on a linear one
        either may fall outside the record, which is what `overhang` being ``None`` says.
    overhang
        The top-strand bases the two cuts leave single-stranded, 5' to 3'; ``""`` for a blunt
        cutter and ``None`` when a cut falls off the end of a linear record.
    certain
        ``False`` when an IUPAC code in the template, rather than a definite base, is what
        allowed the match.
    """

    enzyme: Enzyme
    start: int
    strand: Strand
    _: KW_ONLY
    top_cut: int
    bottom_cut: int
    overhang: str | None
    certain: bool = True

    @property
    def end(self) -> int:
        """Where the match ends, passing the length when it runs across the origin."""
        return self.start + len(self.enzyme.site)

    @property
    def span(self) -> Segment:
        """The matched bases, as a segment of the top strand."""
        return Segment(self.start, self.end)

    @property
    def cuts(self) -> bool:
        """Whether both cuts land in the record, so the enzyme can actually cut here."""
        return self.overhang is not None


@dataclass(frozen=True, slots=True)
class Fragment:
    """One piece a digest leaves, bounded by the cuts in the top strand.

    Parameters
    ----------
    start, end
        The two top-strand cuts that bound it, 0-based and half-open. A fragment across the
        origin of a circular record keeps ``start < end`` and ends past the record's length.
    left_overhang, right_overhang
        The single-stranded bases at each end, written as the TOP strand reads them 5' to 3',
        and ``""`` for a blunt end. Written that way whether the overhang is 5' or 3', and
        whether it sits on the top strand of this fragment or the bottom, so two ends anneal
        exactly when the two strings are equal.
    """

    start: int
    end: int
    left_overhang: str
    right_overhang: str

    @property
    def length(self) -> int:
        """How many bases of top strand it carries, which is what a gel measures."""
        return self.end - self.start


def find_sites(
    record: SequenceRecord, enzymes: EnzymeLike | Iterable[EnzymeLike]
) -> tuple[CutSite, ...]:
    """Return every site one or more enzymes read in `record`, in top-strand order.

    Both strands are searched, and a circular record is searched across its origin.

    Parameters
    ----------
    record
        The record to search.
    enzymes
        One enzyme or several, each an `Enzyme` or a name `get_enzyme` answers to.

    Raises
    ------
    KeyError
        If a name is not one the package ships.

    Examples
    --------
    >>> record = SequenceRecord("AAAAGGTCTCGTTTTCCCC")
    >>> [(site.start, site.overhang) for site in find_sites(record, "BsaI")]
    [(4, 'TTTT')]
    """
    found: list[CutSite] = []
    for enzyme in _resolve(enzymes):
        reverse = reverse_complement(enzyme.site)
        needles = [(Strand.FORWARD, enzyme.site)]
        if reverse != enzyme.site:
            needles.append((Strand.REVERSE, reverse))
        for strand, needle in needles:
            found.extend(
                _hit(record, enzyme, start, strand, needle) for start in _starts(record, needle)
            )
    return tuple(sorted(found, key=lambda site: (site.start, site.enzyme.name, -int(site.strand))))


def has_site(record: SequenceRecord, enzymes: EnzymeLike | Iterable[EnzymeLike]) -> bool:
    """Whether `record` still holds a site any of these enzymes reads."""
    return bool(find_sites(record, enzymes))


def site_counts(
    records: Iterable[SequenceRecord], enzymes: Iterable[EnzymeLike] | None = None
) -> dict[str, int]:
    """Count the sites each enzyme reads across several records, keyed by enzyme name.

    Every enzyme asked about gets an entry, so a count of zero is stated rather than missing.
    `enzymes` defaults to every enzyme the package ships.
    """
    chosen = shipped() if enzymes is None else _resolve(enzymes)
    counts = {enzyme.name: 0 for enzyme in chosen}
    for record in records:
        for site in find_sites(record, chosen):
            counts[site.enzyme.name] += 1
    return counts


def free_enzymes(
    records: Iterable[SequenceRecord], enzymes: Iterable[EnzymeLike] | None = None
) -> tuple[Enzyme, ...]:
    """Return the enzymes with no site in any of `records`, in the order they were given.

    These are the enzymes a Golden Gate design can use without domesticating anything.
    """
    chosen = shipped() if enzymes is None else _resolve(enzymes)
    counts = site_counts(records, chosen)
    return tuple(enzyme for enzyme in chosen if counts[enzyme.name] == 0)


def digest(
    record: SequenceRecord, enzymes: EnzymeLike | Iterable[EnzymeLike]
) -> tuple[Fragment, ...]:
    """Cut `record` with one or more enzymes and return the fragments, in top-strand order.

    A site whose cut falls off the end of a linear record is read but not cut, so it splits
    nothing. An uncut linear record is one fragment, being a molecule already; an uncut
    circular record is none, nothing having been cut.

    Examples
    --------
    >>> [piece.length for piece in digest(SequenceRecord("AAAAGGTCTCGTTTTCCCC"), "BsaI")]
    [11, 8]
    """
    length = len(record)
    #: Top-strand cut -> the overhang it leaves. Two enzymes cutting one position share it.
    boundaries: dict[int, str] = {}
    for site in find_sites(record, enzymes):
        if site.cuts and (record.topology == "circular" or 0 < site.top_cut < length):
            boundaries.setdefault(site.top_cut, site.overhang or "")
    cuts = sorted(boundaries)
    if record.topology == "circular":
        if not cuts:
            return ()
        ends = [*cuts[1:], cuts[0] + length]
    else:
        cuts, ends = [0, *cuts], [*cuts, length]
    return tuple(
        Fragment(start, end, boundaries.get(start, ""), boundaries.get(end % length, ""))
        for start, end in zip(cuts, ends, strict=True)
    )


def insert_site(
    record: SequenceRecord,
    enzyme: EnzymeLike,
    position: int,
    *,
    strand: Strand = Strand.FORWARD,
) -> tuple[SequenceRecord, EditReport]:
    """Put one enzyme's recognition site into `record` before `position`.

    Everything after the site shifts, and the report says which features the insertion fell
    inside. A reverse-strand site is written as the reverse complement, so the enzyme reaches
    back towards lower coordinates to cut.

    Raises
    ------
    ValueError
        If the recognition site holds an IUPAC code, there being no one sequence to write.

    Examples
    --------
    >>> edited, _ = insert_site(SequenceRecord("AAAACCCC"), "BsaI", 4)
    >>> edited.sequence
    'AAAAGGTCTCCCCC'
    """
    one = _resolve(enzyme)[0]
    bases = _definite(one)
    return insert(record, position, reverse_complement(bases) if strand < 0 else bases)


def primer_tail(
    enzyme: EnzymeLike,
    overhang: str = "",
    *,
    flank: str | None = None,
    flank_length: int = FLANK_LENGTH,
    avoid: Iterable[EnzymeLike] = (),
) -> str:
    """Build the 5' tail of a cloning primer, 5' to 3', for the caller to put its own 3' end on.

    The tail is flanking bases, the recognition site, the bases the enzyme reaches over, and
    then `overhang` — so that cutting the amplicon leaves exactly `overhang` single-stranded.

    Parameters
    ----------
    enzyme
        The enzyme the tail is cut by.
    overhang
        The overhang the cut should leave, as long as the enzyme leaves. Empty, and only
        empty, for an enzyme that cuts inside its own site.
    flank
        The bases 5' of the site. Chosen when not given, and checked when it is.
    flank_length
        How many bases to choose when `flank` is not given. NEB recommends six.
    avoid
        Enzymes besides this one whose sites the tail must not spell.

    Raises
    ------
    ValueError
        If `overhang` is not one this enzyme leaves, if `flank` spells a further site or puts a
        Dcm site beside one this enzyme is impaired by, or if no bases could be found that
        avoid both.

    Examples
    --------
    >>> primer_tail("BsaI", "AATG")
    'AAACACGGTCTCAAATG'
    """
    one = _resolve(enzyme)[0]
    overhang = overhang.upper()
    _check_overhang(one, overhang)
    active = (one, *_resolve(avoid))
    site = _definite(one)
    filler = _search(
        max(0, one.top_cut - len(site)), lambda bases: site + bases + overhang, active, one
    )
    core = site + filler + overhang
    if flank is None:
        return _search(flank_length, lambda bases: bases + core, active, one) + core
    flank = flank.upper()
    if (reason := _problem(flank + core, active, one)) is not None:
        raise ValueError(f"flank {flank!r} cannot be used: {reason}")
    return flank + core


def _check_overhang(enzyme: Enzyme, overhang: str) -> None:
    """Refuse an overhang this enzyme would not leave."""
    if enzyme.type == "II":
        if overhang:
            raise ValueError(
                f"{enzyme.name} cuts inside its own site, so its overhang is not one to choose"
            )
    elif len(overhang) != enzyme.overhang_length:
        raise ValueError(
            f"{enzyme.name} leaves a {enzyme.overhang_length}-base overhang, "
            f"and {overhang!r} is {len(overhang)}"
        )


def _definite(enzyme: Enzyme) -> str:
    """Return the recognition site as bases to write, refusing one that holds an IUPAC code."""
    if ambiguous := sorted(set(enzyme.site) - set("ACGT")):
        raise ValueError(
            f"the {enzyme.name} site {enzyme.site} holds the IUPAC code(s) "
            f"{''.join(ambiguous)}, so write the bases you mean instead"
        )
    return enzyme.site


def _problem(tail: str, active: tuple[Enzyme, ...], enzyme: Enzyme) -> str | None:
    """Why this tail will not do, or ``None`` when it will."""
    found = find_sites(SequenceRecord(tail), active)
    if len(found) != 1:
        names = ", ".join(sorted({site.enzyme.name for site in found})) or "none"
        return f"it spells {len(found)} sites ({names}) where one {enzyme.name} site is wanted"
    if enzyme.methylation.get("dcm", "not sensitive") != "not sensitive" and _pattern(
        DCM_SITE
    ).search(tail):
        return f"it puts a Dcm site ({DCM_SITE}) across the {enzyme.name} site"
    return None


def _search(
    length: int, build: Callable[[str], str], active: tuple[Enzyme, ...], enzyme: Enzyme
) -> str:
    """Return the first candidate of `length` bases leaving one site and no Dcm site."""
    for candidate in islice(_candidates(length), _TRIES):
        if _problem(build(candidate), active, enzyme) is None:
            return candidate
    raise ValueError(
        f"no {length} bases leave a clean {enzyme.name} tail; choose the overhang again"
    )


def _candidates(length: int) -> Iterable[str]:
    """Bases to try, in a fixed order, skipping runs and lopsided GC that nobody would order."""
    if length == 0:
        return [""]
    return (
        candidate
        for choice in product("ACGT", repeat=length)
        if not _RUN.search(candidate := "".join(choice))
        and (length < 4 or 0.25 <= sum(base in "GC" for base in candidate) / length <= 0.75)
    )


def _resolve(enzymes: EnzymeLike | Iterable[EnzymeLike]) -> tuple[Enzyme, ...]:
    """Read one enzyme or several, by name or by record, into a tuple of records."""
    if isinstance(enzymes, Enzyme | str):
        enzymes = (enzymes,)
    return tuple(get_enzyme(one) if isinstance(one, str) else one for one in enzymes)


@cache
def _pattern(needle: str) -> re.Pattern[str]:
    """Return a regex matching every stretch of template the site `needle` may read.

    Wrapped in a lookahead so that sites overlapping one another are all found.
    """
    classes = "".join(
        "[" + "".join(sorted(code for code, bases in _BASES.items() if bases & _BASES[want])) + "]"
        for want in needle
    )
    return re.compile(f"(?=({classes}))")


def _starts(record: SequenceRecord, needle: str) -> list[int]:
    """Every top-strand index where `needle` may begin, reading across the origin if circular."""
    length = len(record)
    if len(needle) > length:
        return []
    haystack = record.sequence
    if record.topology == "circular":
        haystack += record.sequence[: len(needle) - 1]
    return [
        match.start() for match in _pattern(needle).finditer(haystack) if match.start() < length
    ]


def _hit(
    record: SequenceRecord, enzyme: Enzyme, start: int, strand: Strand, needle: str
) -> CutSite:
    """Build the hit for one match, reading its cuts and the overhang they leave."""
    length = len(record)
    top, bottom = enzyme.cut_positions(start, strand)
    if record.topology == "circular":
        top, bottom = top % length, bottom % length
    try:
        overhang = enzyme.overhang(record, start, strand)
    except ValueError:
        # A Type IIS site near the end of a linear record: read, but with nothing left to cut.
        overhang = None
    observed = record.extract(Segment(start, start + len(needle)))
    return CutSite(
        enzyme,
        start,
        strand,
        top_cut=top,
        bottom_cut=bottom,
        overhang=overhang,
        certain=all(
            _BASES[seen] <= _BASES[want] for seen, want in zip(observed, needle, strict=True)
        ),
    )

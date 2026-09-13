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
from collections.abc import Iterable, Mapping
from dataclasses import KW_ONLY, dataclass
from functools import cache

from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.enzymes import enzymes as shipped
from liulab_mbio.sequence import Segment, SequenceRecord, Strand, reverse_complement

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

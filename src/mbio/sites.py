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

from mbio.checks import counted
from mbio.codons import CodonUsage, amino_acid, codon_usage
from mbio.edits import EditReport, insert, replace
from mbio.enzymes import EndType, Enzyme, get_enzyme
from mbio.enzymes import enzymes as shipped
from mbio.sequence import Feature, Segment, SequenceRecord, Strand, reverse_complement

#: How many bases NEB recommends 5' of a recognition site for an enzyme to cut near an end.
SPACER_LENGTH = 6

#: Dcm methylates the inner cytosine of this site. Only an enzyme whose supplier says Dcm
#: affects it is held to the rule; `Enzyme.methylation` carries that answer.
DCM_SITE = "CCWGG"

#: Dam methylates the adenine of this site, on the same terms. `mbio.bench.steps` keeps
#: its own copy for DpnI, that being a question about a plasmid against a fresh amplicon and not
#: about whether an enzyme still cuts.
DAM_SITE = "GATC"

#: How many candidate spacers or fillers to try before giving up on an overhang.
_TRIES = 4096

#: How much identical sequence either side of a site puts it past every oligo placed there —
#: mutagenic primer, assembly overlap or Type IIS overhang — because none is unique to one copy.
#: `docs/research/domestication-methods.md`, the rule and section 3.
OLIGO_REACH = 100

_RUN = re.compile(r"(.)\1{3}")

#: What each IUPAC code stands for. `mbio.sequence.IUPAC_DNA` names the same fifteen.
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
        whether it sits on the top strand of this fragment or the bottom. The string is half of
        what says whether two ends anneal; the other half is the end type, which a caller reads
        off the enzyme that made the cut. `mbio.overhangs.compatible` is the rule.
    """

    start: int
    end: int
    left_overhang: str
    right_overhang: str

    @property
    def length(self) -> int:
        """How many bases of top strand it carries, which is what a gel measures."""
        return self.end - self.start


@dataclass(frozen=True, slots=True)
class Domestication:
    """One synonymous codon change, and the site it took away.

    Parameters
    ----------
    site
        The site that was there before the change.
    feature
        The coding sequence the codon belongs to, as it read before the change.
    position
        0-based index on the TOP strand where the codon's three bases begin. A codon of a
        reverse-strand coding sequence reads back from ``position + 2``.
    codon_index
        Which codon of the coding sequence this is, counting from zero.
    old_codon, new_codon
        The codon before and after, read 5' to 3' along the coding sequence.
    amino_acid
        The one-letter amino acid both codons spell.
    """

    site: CutSite
    feature: Feature
    position: int
    codon_index: int
    old_codon: str
    new_codon: str
    amino_acid: str


@dataclass(frozen=True, slots=True)
class DomesticationReport:
    """What domestication changed, and what it left for someone to decide.

    Attributes
    ----------
    changes
        Each site taken away, with the codon change that did it.
    outside_cds
        Sites lying in no coding sequence. Removing one of these changes what the record
        spells, so it is reported and the choice is left to the caller.
    unchanged
        Sites in a coding sequence that no synonymous change could take away.
    """

    changes: tuple[Domestication, ...] = ()
    outside_cds: tuple[CutSite, ...] = ()
    unchanged: tuple[CutSite, ...] = ()


@dataclass(frozen=True, slots=True)
class Blocked:
    """One candidate enzyme a search ruled out, and what ruled it out.

    Parameters
    ----------
    enzyme
        The candidate.
    reason
        Why it is not free, in the terms the search was given.
    sites
        The sites standing in its way, empty where the ends it leaves ruled it out instead.
    """

    enzyme: Enzyme
    reason: str
    sites: tuple[CutSite, ...] = ()


@dataclass(frozen=True, slots=True)
class EnzymeSearch:
    """Which candidates are free to cut alongside a set of records, and what blocked the rest.

    Attributes
    ----------
    free
        The free enzymes, the rarest recognition site first and ties in the order they were
        given. Empty is a finding, not an error.
    blocked
        Every candidate that is not free, in the order they were given.
    """

    free: tuple[Enzyme, ...] = ()
    blocked: tuple[Blocked, ...] = ()

    @property
    def best(self) -> Enzyme | None:
        """The first free enzyme, or ``None`` where none is free."""
        return self.free[0] if self.free else None


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


def free_enzyme_search(
    records: Iterable[SequenceRecord],
    enzymes: Iterable[EnzymeLike] | None = None,
    *,
    overhang_length: int | None = None,
    end: EndType | None = None,
) -> EnzymeSearch:
    """Search `enzymes` for the ones free to cut alongside `records`, and say what blocked the rest.

    A candidate is free when it leaves the end the caller asked for and reads no site in any of
    the records. An empty result is an answer: `EnzymeSearch.blocked` then names the sites that
    stand in each candidate's way, which is what a caller reads before reaching for a staged
    digest, a bridging assembly or a re-tailoring PCR.

    Parameters
    ----------
    records
        The molecules the enzyme would share a reaction with, outside what it is meant to cut.
    enzymes
        The candidates, defaulting to every enzyme the package ships.
    overhang_length
        How many bases the cut must leave single-stranded, where a design's overhangs are fixed.
    end
        The end the cut must leave: ``"5'"``, ``"3'"`` or ``"blunt"``.

    Examples
    --------
    >>> search = free_enzyme_search([SequenceRecord("AAAAGGTCTCGTTTT")], ("BsaI", "PaqCI"))
    >>> [enzyme.name for enzyme in search.free], search.blocked[0].sites[0].start
    (['PaqCI'], 4)
    """
    pool = tuple(records)
    free: list[Enzyme] = []
    blocked: list[Blocked] = []
    for enzyme in shipped() if enzymes is None else _resolve(enzymes):
        if overhang_length is not None and enzyme.overhang_length != overhang_length:
            blocked.append(
                Blocked(
                    enzyme,
                    f"leaves {enzyme.overhang_length} bases single-stranded, not {overhang_length}",
                )
            )
            continue
        if end is not None and enzyme.end != end:
            blocked.append(Blocked(enzyme, f"leaves a {enzyme.end} end, not a {end} one"))
            continue
        found = {
            record.name or f"record {index + 1}": find_sites(record, enzyme)
            for index, record in enumerate(pool)
        }
        if sites := tuple(site for hits in found.values() for site in hits):
            where = ", ".join(f"{len(hits)} in {name}" for name, hits in found.items() if hits)
            blocked.append(Blocked(enzyme, f"reads {counted(len(sites), 'site')}: {where}", sites))
        else:
            free.append(enzyme)
    return EnzymeSearch(tuple(sorted(free, key=lambda one: -len(one.site))), tuple(blocked))


def repeat_context(record: SequenceRecord, span: Segment) -> int:
    """Return how much identical sequence flanks the best other copy of `span` in `record`.

    Walking outwards from each further copy of the bases `span` spells, this is the largest
    flank the two copies share on **both** sides. Zero where the bases appear only once.

    No oligo shorter than this reaches one copy and not the other, which is what puts a site
    inside a long repeat past domestication: see `OLIGO_REACH`.

    Examples
    --------
    >>> repeat_context(SequenceRecord("TTGGCCAATTGGCCAA"), Segment(2, 6))
    2
    """
    bases = record.bases(span.start, span.end)
    length = len(record)
    circular = record.topology == "circular"
    here = span.start % length

    def base(index: int) -> str | None:
        """Return the base at `index`, or ``None`` off the end of a linear record."""
        if circular:
            return record.sequence[index % length]
        return record.sequence[index] if 0 <= index < length else None

    best = 0
    for start in _starts(record, bases):
        if start % length == here:
            continue
        shared = 0
        while shared < (length - len(bases)) // 2:
            pairs = (
                (base(here - shared - 1), base(start - shared - 1)),
                (base(here + len(bases) + shared), base(start + len(bases) + shared)),
            )
            if any(one is None or one != other for one, other in pairs):
                break
            shared += 1
        best = max(best, shared)
    return best


def out_of_reach(
    record: SequenceRecord, sites: Iterable[CutSite], *, reach: int = OLIGO_REACH
) -> tuple[CutSite, ...]:
    """Return the sites no oligo can reach, because a long repeat carries each more than once.

    Checked before anything else a domestication route decides: a site here is not a bench job
    at all, whatever its codon context says.
    """
    return tuple(site for site in sites if repeat_context(record, site.span) >= reach)


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


def released(
    record: SequenceRecord, enzymes: EnzymeLike | Iterable[EnzymeLike]
) -> tuple[Fragment, ...]:
    """Return the digest's fragments with an overhang at each end: what a ligation takes.

    A fragment blunt at either end carries an end of a linear record, or an end a blunt cutter
    made, so matched overhangs cannot put it back in.

    Examples
    --------
    >>> cut = released(SequenceRecord("AAAAGGTCTCGTTTTCCCCCCGAGACCAAAA"), "BsaI")
    >>> [fragment.length for fragment in cut]
    [5]
    """
    return tuple(one for one in digest(record, enzymes) if one.left_overhang and one.right_overhang)


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
    spacer: str | None = None,
    spacer_length: int = SPACER_LENGTH,
    avoid: Iterable[EnzymeLike] = (),
) -> str:
    """Build the 5' tail of a cloning primer, 5' to 3', for the caller to put its own 3' end on.

    The tail is a spacer, the recognition site, the bases the enzyme reaches over, and
    then `overhang` — so that cutting the amplicon leaves exactly `overhang` single-stranded.

    Parameters
    ----------
    enzyme
        The enzyme the tail is cut by.
    overhang
        The overhang the cut should leave, as long as the enzyme leaves. Empty, and only
        empty, for an enzyme that cuts inside its own site.
    spacer
        The bases 5' of the site. Chosen when not given, and checked when it is.
    spacer_length
        How many bases to choose when `spacer` is not given. NEB recommends six.
    avoid
        Enzymes besides this one whose sites the tail must not spell.

    Raises
    ------
    ValueError
        If `overhang` is not one this enzyme leaves, if `spacer` spells a further site or puts a
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
    if spacer is None:
        return _search(spacer_length, lambda bases: bases + core, active, one) + core
    spacer = spacer.upper()
    if (reason := _problem(spacer + core, active, one)) is not None:
        raise ValueError(f"spacer {spacer!r} cannot be used: {reason}")
    return spacer + core


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


def domesticate(
    record: SequenceRecord,
    enzymes: EnzymeLike | Iterable[EnzymeLike],
    *,
    usage: CodonUsage | None = None,
    avoid: Iterable[EnzymeLike] = (),
) -> tuple[SequenceRecord, DomesticationReport]:
    """Take away every site of `enzymes` that a synonymous codon change can reach.

    A site inside a coding sequence goes by changing one codon for another spelling the same
    amino acid, keeping the reading frame and the protein. The replacement is the one the host
    uses most often, among those that take the site away and spell no new site for `enzymes` or
    `avoid`. A site lying in no coding sequence is reported and **not** edited: removing it
    would change what the record spells, which is the caller's decision to make.

    A coding sequence is one the record types ``CDS``, **or** any feature whose own bases spell
    a whole protein: `ATG`, whole codons of definite bases, and one stop as the last codon. So a
    vector whose open reading frames arrive as `misc_feature` is still domesticated, while an
    annotation too loose to read in one frame is left outside. `docs/adr/0013-coding-by-bases.md`
    says why the bases decide and no frame is searched for.

    Parameters
    ----------
    record
        The record to domesticate.
    enzymes
        The enzyme or enzymes whose sites should go.
    usage
        The host's codon usage. The shipped *E. coli* K-12 table by default.
    avoid
        Further enzymes whose sites no change may create.

    Returns
    -------
    tuple[SequenceRecord, DomesticationReport]
        The edited record, and what changed and what did not.
    """
    table = usage if usage is not None else codon_usage()
    targets = _resolve(enzymes)
    active = targets + _resolve(avoid)
    changes: list[Domestication] = []
    outside: list[CutSite] = []
    unchanged: list[CutSite] = []
    handled: set[tuple[str, int, int]] = set()
    while True:
        pending = [
            site
            for site in find_sites(record, targets)
            if (site.enzyme.name, site.start, int(site.strand)) not in handled
        ]
        if not pending:
            return record, DomesticationReport(tuple(changes), tuple(outside), tuple(unchanged))
        site = pending[0]
        handled.add((site.enzyme.name, site.start, int(site.strand)))
        feature = _coding_feature(record, site)
        if feature is None:
            outside.append(site)
        elif (swap := _synonymous(record, site, feature, table, active)) is None:
            unchanged.append(site)
        else:
            record, change = swap
            changes.append(change)


def _coding_feature(record: SequenceRecord, site: CutSite) -> Feature | None:
    """Return the first coding sequence the site touches, or ``None`` when it touches none."""
    for feature in record.features:
        if any(record.covers(feature, at) for at in range(site.start, site.end)) and (
            feature.type == "CDS" or _reads_as_coding(record, feature)
        ):
            return feature
    return None


def _reads_as_coding(record: SequenceRecord, feature: Feature) -> bool:
    """Whether a feature's own bases spell a whole protein, whatever the record types it.

    `ATG`, then whole codons of definite bases, then one stop as the last codon and none before
    it. The feature gives the frame and the strand, so nothing is searched for.
    """
    bases = record.extract(feature)
    if len(bases) < 6 or len(bases) % 3 or set(bases) - set("ACGT") or bases[:3] != "ATG":
        return False
    codons = [bases[at : at + 3] for at in range(0, len(bases), 3)]
    return amino_acid(codons[-1]) == "*" and all(amino_acid(one) != "*" for one in codons[:-1])


def _coding_positions(feature: Feature) -> list[int]:
    """Return each base of a coding sequence as its segments index it, in codon order."""
    positions = [
        index for segment in feature.segments for index in range(segment.start, segment.end)
    ]
    if feature.strand == Strand.REVERSE:
        positions.reverse()
    return positions


def _codon_span(record: SequenceRecord, three: list[int]) -> tuple[int, int] | None:
    """One codon's span on the top strand, or ``None`` when its bases do not run together."""
    if len(three) != 3:
        return None
    top = three if _follows(record, three[0], three[1]) else three[::-1]
    if not (_follows(record, top[0], top[1]) and _follows(record, top[1], top[2])):
        return None
    # `replace` and the report both take a start inside the record.
    start = top[0] % len(record)
    return start, start + 3


def _follows(record: SequenceRecord, base: int, after: int) -> bool:
    """Whether the top strand reads `after` straight after `base`, across the origin too."""
    return record.covers(Segment(base + 1, base + 2), after)


def _synonymous(
    record: SequenceRecord,
    site: CutSite,
    feature: Feature,
    table: CodonUsage,
    active: tuple[Enzyme, ...],
) -> tuple[SequenceRecord, Domestication] | None:
    """Change one codon to take the site away, or ``None`` when no synonymous codon does."""
    positions = _coding_positions(feature)
    coding = record.extract(feature)
    if len(coding) != len(positions):
        return None
    where = {position: index for index, position in enumerate(positions)}
    before = len(find_sites(record, active))
    candidates: list[tuple[float, int, int, str, str, str, tuple[int, int]]] = []
    touched = {index // 3 for at, index in where.items() if record.covers(site.span, at)}
    for index in sorted(touched):
        old = coding[3 * index : 3 * index + 3]
        span = _codon_span(record, positions[3 * index : 3 * index + 3])
        if len(old) != 3 or span is None or set(old) - set("ACGT"):
            continue
        amino = table.amino_acid(old)
        candidates.extend(
            (
                -table.fraction(new),
                sum(a != b for a, b in zip(old, new, strict=True)),
                index,
                old,
                new,
                amino,
                span,
            )
            for new in table.synonymous(old)
            if new != old
        )
    for _, _, index, old, new, amino, span in sorted(candidates):
        bases = reverse_complement(new) if feature.strand == Strand.REVERSE else new
        edited, _ = replace(record, span[0], span[1], bases)
        after = find_sites(edited, active)
        if len(after) != before - 1 or any(_is(one, site) for one in after):
            continue
        return edited, Domestication(site, feature, span[0], index, old, new, amino)
    return None


def _is(one: CutSite, other: CutSite) -> bool:
    """Whether two hits are the same site of the same enzyme on the same strand."""
    return (one.enzyme.name, one.start, one.strand) == (
        other.enzyme.name,
        other.start,
        other.strand,
    )


def _resolve(enzymes: EnzymeLike | Iterable[EnzymeLike]) -> tuple[Enzyme, ...]:
    """Read one enzyme or several, by name or by record, into a tuple of records.

    An enzyme named twice is kept once, where it first appears: a site is one site however
    many times its enzyme was listed.
    """
    if isinstance(enzymes, Enzyme | str):
        enzymes = (enzymes,)
    read = (get_enzyme(one) if isinstance(one, str) else one for one in enzymes)
    return tuple(dict.fromkeys(read))


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
    observed = record.bases(start, start + len(needle))
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

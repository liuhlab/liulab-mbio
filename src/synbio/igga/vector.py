"""Take the user's own destination vector, and make it one where it is not.

A vector is a destination when one enzyme excises one piece from it, leaving the entry overhang
at one end and the cloning scar at the other. That is what a round opens, so it is read off a
digest rather than off a string: a vector already carrying such a piece is taken as it stands,
and one that does not has one put at a site the user names.

A `Cassette` says which piece, to which enzyme, and in which tube. A round's destination gives up
the method's internal stuffer to the internal enzyme; a donor backbone gives the same stuffer up
to the external enzyme; the working vector gives up its ccdB cassette to the cargo enzyme, and
that cassette is the method's own DNA — `CCDB_PAYLOAD` between two of that enzyme's sites. What
a record may not read outside its cassette is the enzymes of its own tube and no others, which
is what lets a destination carry outboard the sites that release its cargo:
`docs/adr/0014-one-record-one-tube.md`.

An insertion is reported as an `mbio.edits.EditReport`, so the user reads what it changed
rather than trusting a new file. Coordinates are the model's, and a cassette across the origin
ends past the record's length.
"""

from collections.abc import Iterable
from dataclasses import dataclass

from mbio.bench.goldengate import GOLDEN_GATE_ENZYMES, LAST_RESORT
from mbio.checks import Check, counted
from mbio.codons import CodonUsage
from mbio.edits import EditReport, insert, replace
from mbio.enzymes import Enzyme, get_enzyme
from mbio.sequence import (
    Feature,
    Segment,
    SequenceRecord,
    position_text,
    reverse_complement,
)
from mbio.sites import (
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
from synbio.igga.method import IGGA, INTERFACE_OVERHANGS, Scheme, refuse

#: Where a cassette goes: the name of a feature, or a ``(start, end)`` span.
type Site = str | tuple[int, int]

#: ``BBa_J23119``, the promoter the working vector's cassette reads ccdB from, and ``BBa_B0034``,
#: its ribosome binding site, as the iGEM Registry publishes them. The promoter's ``TTGACA`` and
#: ``TATAAT`` are the sigma-70 consensus, so the cassette kills in any ccdB-sensitive strain
#: rather than only in one carrying T7 polymerase. #287 decided both.
J23119 = "TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC"
B0034 = "AAAGAGGAGAAA"

#: What spaces the ribosome binding site from the start codon: the RFC10 scar, the de facto
#: BioBrick spacing. No scar-free ``J23119``-``B0034``-CDS spacer is published, so #287 records
#: this as a choice and not a citation.
BIOBRICK_SCAR = "TACTAG"

#: The ccdB coding sequence the cassette counter-selects with, read off SnapGene's own offline
#: ``pET-53-DEST`` map. It is Invitrogen's Gateway ccdB, which it agrees with over its first 45
#: bases, carrying one synonymous change that takes a BsaI site out of it. #287 decided it and
#: names the lab's own copy, DMX0001, Addgene 247434.
CCDB = (
    "ATGCAGTTTAAGGTTTACACCTATAAAAGAGAGAGCCGTTATCGTCTGTTTGTGGATGTACAGAGTGATATT"
    "ATTGACACGCCCGGGCGACGGATGGTGATCCCCCTGGCCAGTGCACGTCTGCTGTCAGATAAAGTCTCCCGT"
    "GAACTTTACCCGGTGGTGCATATCGGGGATGAAAGCTGGCGCATGATGACCACCGATATGGCCAGTGTGCCG"
    "GTTTCCGTTATCGGGGAAGAAGTGGCTGATCTCAGCCACCGCGAAAATGACATCAAAAACGCCATTAACCTG"
    "ATGTTCTGGGGAATATAA"
)

#: What the working vector's cassette carries between the two cargo-enzyme sites that release it:
#: promoter, ribosome binding site, scar and ccdB, 359 bases (#287). #287 measured it free of
#: BsaI, BsmBI, BbsI, PaqCI, SapI and PmeI on both strands, with one SrfI inside ccdB, which the
#: method permits inside a cassette.
CCDB_PAYLOAD = J23119 + B0034 + BIOBRICK_SCAR + CCDB


@dataclass(frozen=True, slots=True)
class Cassette:
    """The piece a vector gives up to admit a part, and what may not be read outside it.

    Parameters
    ----------
    bases
        The piece as the vector carries it: the entry overhang first, the cloning scar last.
    enzyme
        The enzyme that excises it, which is what opens the vector.
    free_of
        The enzymes no site of which may lie outside it: the ones acting in `tube` that would cut
        what this record keeps. An enzyme acting in another tube never meets this record, and a
        chopper aimed at the piece this tube throws away is meant to cut outside.
    tube
        The digest this record is cut in, named as a message can say it.
    """

    bases: str
    enzyme: Enzyme
    free_of: tuple[Enzyme, ...]
    tube: str


def round_cassette(scheme: Scheme = IGGA) -> Cassette:
    """Return the cassette a round's destination gives up: the method's internal stuffer.

    Barred outside it are the enzymes of that one digest: the internal enzyme, which opens the
    library, and the chopper that shreds the stuffer it gives up. The external enzyme and the
    donor's own chopper act in the donor's tube, and the destination carries their sites outboard
    of the cassette on purpose -- that is what releases the cargo once the rounds are done.
    """
    return Cassette(
        scheme.internal_stuffer,
        scheme.internal,
        (scheme.internal, *scheme.blunt_for_the_destination),
        "the digest that opens a round's destination",
    )


def donor_cassette(scheme: Scheme = IGGA) -> Cassette:
    """Return the cassette a donor backbone gives up: the same stuffer, to the external enzyme.

    One stuffer, two tubes. A part is held in a DMX backbone whose flanks carry the external
    enzyme's sites, so what that backbone gives up is the piece between them.

    Nothing is barred outside it, because this is the one tube that keeps nothing outside it: the
    backbone is what the digest throws away. The releasing cuts lie outboard of the cassette, and
    so does the blunt chopper aimed at the backbone they leave. A third releasing site is caught
    where it does harm, by the gate counting the cuts in the tube, and where a chopper may sit is
    `synbio.igga.gate.check_dmx_vector`'s question.

    The sites that excise this cassette lie in the backbone and not in the cassette, so a backbone
    carrying none cannot be made a donor by putting one in.
    """
    return Cassette(
        scheme.internal_stuffer,
        scheme.external,
        (),
        "the digest that releases a part from its donor",
    )


def ccdb_cassette(
    cargo: Enzyme, *, scheme: Scheme = IGGA, n_part: bool = True, c_part: bool = True
) -> Cassette:
    """Return the cassette a working vector gives up when `cargo` admits the library.

    `CCDB_PAYLOAD` between two of `cargo`'s sites, each facing out, so the digest leaves the
    entry overhang at one end and the cloning scar at the other -- the pair the library already
    presents. Only BsmBI and `cargo` are barred outside it: those are the two enzymes that ever
    share a tube with this vector, and holding it to the rest would turn away a backbone whose
    other sites never meet them.

    Parameters
    ----------
    cargo
        The enzyme `cargo_enzyme` chose for this vector.
    scheme
        The method, which says what the cassette's two ends are and which enzyme is reserved.
    n_part
        Whether the working cassette carries a part 5' of this one. Where it does not, the
        cassette carries that end itself, over the one glycine codon the method's layout spells
        there.
    c_part
        The same, 3' of it.

    Raises
    ------
    ValueError
        If `cargo` does not read exactly the two sites that release the payload.
    """
    entry, scar = scheme.entry_overhang, scheme.scar_overhang
    start, end = INTERFACE_OVERHANGS["working cassette"]
    reach = ("TA" * cargo.top_cut)[: cargo.top_cut - len(cargo.site)]
    bases = (
        entry + reach + reverse_complement(cargo.site) + CCDB_PAYLOAD + cargo.site + reach + scar
    )
    if not n_part:
        bases = start + "GG" + bases
    if not c_part:
        bases = bases + "GG" + end
    read = find_sites(SequenceRecord(bases), cargo)
    if len(read) != len(INTERFACE_OVERHANGS["cargo"]):
        raise ValueError(
            f"{cargo.name} reads {len(read)} site(s) in this ccdB cassette, where the two that "
            "release the payload are all it may read: choose another cargo enzyme"
        )
    return Cassette(
        bases, cargo, (cargo, *scheme.reserved_enzymes), "the digest that opens the working vector"
    )


@dataclass(frozen=True, slots=True)
class Destination:
    """A vector one enzyme can open, and what making it one changed.

    Parameters
    ----------
    record
        The vector. The record it was given, unchanged, where that already carried a cassette.
    stuffer
        The piece the cassette's enzyme excises to open the vector, bounded by its two cuts. It
        ends past the record's length when it runs across the origin.
    edit
        What inserting the cassette did to the features and primer binding sites it did not simply
        shift, or ``None`` where the vector already carried one.
    """

    record: SequenceRecord
    stuffer: Segment
    edit: EditReport | None = None


def destination_vector(
    vector: SequenceRecord,
    scheme: Scheme,
    *,
    site: Site | None = None,
    cassette: Cassette | None = None,
) -> Destination:
    """Return `vector` as a destination `cassette`'s enzyme can open.

    A vector that already gives up one piece leaving the method's two overhangs is returned as it
    stands, `site` unread. One that does not has `cassette` put at the start of `site`, and the
    result is held to the same rule. The default is `round_cassette`, which is what a round's own
    destination gives up; `donor_cassette` and `ccdb_cassette` are what a donor backbone and a
    working vector give up instead, each barring the enzymes of its own tube.

    Parameters
    ----------
    vector
        The user's own vector, circular.
    scheme
        The method the build is given, which says what the two overhangs are.
    site
        Where to put a cassette, as a feature name or a ``(start, end)`` span. Only read where one
        has to be put.
    cassette
        What the vector gives up, and the enzymes its backbone must be free of.

    Raises
    ------
    ValueError
        If one of the barred enzymes reads a site outside the cassette, if a cassette has to be
        put and no site is named, if `site` names no feature or does not fit, or if the cassette
        is not a whole number of codons where it would land inside a coding sequence.
    KeyError
        If the method names an enzyme this package does not ship.
    """
    held = round_cassette(scheme) if cassette is None else cassette
    found = _stuffer(vector, scheme, held.enzyme)
    if found is not None:
        _check_clean(vector, held, found)
        return Destination(vector, found)
    at = _position(vector, held, site)
    _check_frame(vector, at, held.bases)
    edited, report = insert(vector, at, held.bases)
    made = _stuffer(edited, scheme, held.enzyme)
    if made is None:
        raise ValueError(
            f"the cassette put at {at} left a vector {held.enzyme.name} does not open on "
            f"{scheme.entry_overhang!r} and {scheme.scar_overhang!r}: the bases it now joins "
            "spell a further site. Name another site"
        )
    _check_clean(edited, held, made)
    return Destination(edited, made, report)


def entry_destination(destination: Destination, scheme: Scheme) -> Destination:
    """Return `destination` respelt so the piece it gives up enters on `scheme`'s entry overhang.

    A part enters on its own position's overhang, so a build over several positions needs one
    destination each. They differ in the bases that overhang spells and in nothing else: those
    bases lie between the two cuts, not in the sites that make them, so the enzyme still opens
    the vector there. Four new bases can still spell a site the backbone did not read before,
    so the result goes back through `destination_vector` rather than being taken on trust.

    Raises
    ------
    ValueError
        For any reason `destination_vector` refuses the respelt vector.
    """
    at = destination.stuffer.start
    entry = scheme.entry_overhang
    if destination.record.bases(at, at + len(entry)) == entry:
        return destination
    record, _ = replace(destination.record, at, at + len(entry), entry)
    made = destination_vector(record, scheme)
    return Destination(made.record, made.stuffer, destination.edit)


def released_cargo(record: SequenceRecord, scheme: Scheme = IGGA) -> Segment | None:
    """Return the span `scheme`'s external enzyme frees from `record`, or ``None`` where none.

    What the final assembly moves: the piece bounded by the entry overhang and the cloning scar,
    which is the whole cargo a library's rounds built. It is read off a digest rather than off a
    length, so a record carrying no such site answers ``None`` and the step says so.

    Examples
    --------
    >>> released_cargo(SequenceRecord("ACGT")) is None
    True
    """
    return _stuffer(record, scheme, scheme.external)


def _stuffer(record: SequenceRecord, scheme: Scheme, enzyme: Enzyme) -> Segment | None:
    """Return the piece `enzyme` excises from `record`, or ``None`` where none is.

    The piece is the one whose two ends are the overhangs a part enters and leaves on, which is
    what makes the vector a destination rather than a plasmid with a cassette-shaped gap.
    """
    entry, scar = scheme.entry_overhang, scheme.scar_overhang
    for piece in digest(record, enzyme):
        if piece.left_overhang == entry and piece.right_overhang == scar:
            return Segment(piece.start, piece.end)
    return None


def _position(vector: SequenceRecord, cassette: Cassette, site: Site | None) -> int:
    """Return where the cassette goes, from the feature or the span the user names."""
    if site is None:
        raise ValueError(
            f"this vector carries no cassette for {cassette.enzyme.name} to excise: name the "
            "site to put one at, as a feature name or a (start, end) span"
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
    """Refuse a cassette that is not whole codons where it would land inside a coding sequence."""
    if len(block) % 3 == 0:
        return
    coding = _coding(vector, at)
    if coding is not None:
        raise ValueError(
            f"the cassette is {len(block)} bases, which is not a whole number of codons, "
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


def _check_clean(record: SequenceRecord, cassette: Cassette, stuffer: Segment) -> None:
    """Refuse a vector reading outside its cassette an enzyme of the tube it is digested in.

    Only that one tube is asked about. Inside the cassette its enzymes are the design: the cuts
    that excise the cassette, and whichever blunt chopper shreds the excised piece. Outside, that
    digest would cut the backbone. A method enzyme acting in some other tube is not barred, and a
    destination carries the ones that release its cargo outboard of the cassette on purpose.
    """
    for site in find_sites(record, cassette.free_of):
        if not record.covers(stuffer, site.span):
            raise ValueError(
                f"{site.enzyme.name} reads a site at {site.start} on the "
                f"{site.strand.name.lower()} strand, outside the cassette at "
                f"{stuffer.start}-{stuffer.end}: {site.enzyme.name} acts in {cassette.tube}, "
                "which would cut the backbone there. Take that site out of the vector, or name a "
                "method whose enzymes it is free of"
            )


def cargo_candidates(scheme: Scheme = IGGA) -> tuple[str, ...]:
    """Return the enzymes a cargo enzyme is searched among, in the order the search reads them.

    The Type IIS enzymes `mbio` ships a Golden Gate reaction and cycling protocol for,
    less any reading a reserved enzyme's site and less its own last resort, which ships no
    protocol. Derived rather than listed, so an enzyme the package stops shipping a protocol for
    leaves this list with it.
    """
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
        What `mbio.sites.free_enzyme_search` returned.
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
    the shelf. This one is searched per vector, against the method's own rule: distinct from the
    enzyme that inserts the working cassette, a 4 nt 5' overhang so the method's four junction
    overhangs are unchanged, and no site in any molecule sharing its reaction.

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
        # The cassette is built around this enzyme after it is chosen, so the pot does carry its
        # sites: what the search settles is that nothing else there does.
        detail = f"{chosen.name} cuts nothing in the final-assembly pot but the cassette it opens"
    return Cargo(
        chosen,
        search,
        Check("cargo enzyme", "pass" if chosen else "fail", len(search.free), detail),
    )


@dataclass(frozen=True, slots=True)
class Working:
    """The vector a finished library moves into, and the enzyme that admits it.

    Parameters
    ----------
    destination
        The working vector, and what putting its ccdB cassette in changed.
    cargo
        Which enzyme admits the library, and what ruled the others out.
    cassette
        The ccdB cassette that enzyme releases, which is what the library displaces.
    """

    destination: Destination
    cargo: Cargo
    cassette: Cassette

    @property
    def record(self) -> SequenceRecord:
        """The working vector itself."""
        return self.destination.record

    @property
    def enzyme(self) -> Enzyme:
        """The cargo enzyme. A `Working` is only built where one was found."""
        found = self.cargo.enzyme
        if found is None:  # pragma: no cover - working_vector refuses before building one
            raise ValueError("this working vector has no cargo enzyme")
        return found


def working_vector(
    backbone: SequenceRecord,
    cargo: Iterable[SequenceRecord],
    *,
    scheme: Scheme = IGGA,
    site: Site | None = None,
    enzyme: Enzyme | None = None,
) -> Working:
    """Return `backbone` as the working vector `cargo` is assembled into.

    The enzyme is chosen first, over every molecule sharing the final-assembly pot, and the ccdB
    cassette is then written with that enzyme's sites and put in. A backbone already carrying one
    is taken as it stands.

    Parameters
    ----------
    backbone
        The user's own vector, circular, before its ccdB cassette.
    cargo
        The finished library, and anything else sharing the final-assembly pot.
    scheme
        The method, which says what the cassette's two ends are.
    site
        Where the cassette goes, read only where the backbone carries none.
    enzyme
        The cargo enzyme, where a caller chose it from `backbone` alone before designing what is
        in `cargo`. It is the only candidate searched, so the pot still has to be free of it.

    Raises
    ------
    ValueError
        If no candidate enzyme is free of every molecule, or for any reason
        `destination_vector` refuses the backbone.
    """
    chosen = cargo_enzyme(
        (backbone, *cargo), scheme=scheme, candidates=None if enzyme is None else (enzyme.name,)
    )
    if chosen.enzyme is None:
        raise ValueError(chosen.check.detail)
    held = ccdb_cassette(chosen.enzyme, scheme=scheme)
    return Working(destination_vector(backbone, scheme, site=site, cassette=held), chosen, held)


@dataclass(frozen=True, slots=True)
class Domesticated:
    """A vector with the method's sites taken out of it, and what could not be taken out.

    Parameters
    ----------
    record
        The vector as the changes left it.
    report
        What `mbio.sites.domesticate` changed, and the sites it left for someone to
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

    A site in no coding sequence is left alone too, by `mbio.sites.domesticate`'s own
    rule: changing it changes what the vector spells, which no method decides on its own.
    """
    held = tuple(enzymes)
    refused = out_of_reach(vector, find_sites(vector, held), reach=reach)
    record, report = domesticate(vector, held, usage=usage, avoid=held)
    left = tuple(sorted((*report.outside_cds, *report.unchanged), key=lambda one: one.start))
    if not left:
        detail = f"{counted(len(report.changes), 'site')} changed; none of {_named(held)} is left"
    else:
        named = ", ".join(
            f"{one.enzyme.name} at {position_text(one.start, len(vector))}" for one in left
        )
        blocked = {one.start for one in refused}
        inside = ", ".join(
            position_text(one.start, len(vector)) for one in left if one.start in blocked
        )
        detail = f"{counted(len(report.changes), 'site')} changed, {len(left)} left: {named}" + (
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

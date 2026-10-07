"""The method's gate: which predicate, with which parameters, on which molecule.

Every rule here is one of `liulab_mbio`'s predicates called with this method's own parameters,
and every message is in the method's words. The gate reads a reaction -- the molecules in one
tube, each in its role, and the enzymes acting there -- so it is blind to how a design was
reached and says which tube a failure belongs to. A design an agent composed and one
`plan_igga` wrote are judged the same way.

A `Judgement` carries the `liulab_mbio.checks.Check` and the findings behind it, each a domain
object mbio already returns, so a failure can be drawn on the record it occurred in. A failing
check names what is wrong and need not name a remedy.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from liulab_mbio.barcodes import check_barcodes, separation
from liulab_mbio.checks import STATUSES, Check, Status, worst_of
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.overhangs import fidelity
from liulab_mbio.primers.placement import find_binding_sites
from liulab_mbio.reaction import Pool, Reaction, Role
from liulab_mbio.sequence import Segment, SequenceRecord, Strand, reverse_complement
from liulab_mbio.sites import CutSite, Fragment, digest, find_sites
from liulab_mbio.translate import stop_codons
from liulab_synbio.igga.parts import barcode_rules
from liulab_synbio.igga.project import Project

#: What stands behind a check: a domain object `liulab_mbio` already returns, so `plot` can draw
#: a failure on the record it occurred in. No finding type of its own.
type Finding = CutSite | Segment

#: How many pieces the method's enzymes cut one molecule of a round into: a destination gives up
#: its stuffer, and a block gives up its cargo, so each is cut in two places.
CUTS = 2

#: How many overhangs one round's ligation joins on: the one a part enters by, and the cloning
#: scar. Both ends of every molecule in that tube present one of these two.
LIGATION_OVERHANGS = 2

#: What the primers that read a well anneal to, 5' to 3'. LevSeq publishes them and they bind
#: the DMX vector unchanged, so a design holds them rather than choosing them:
#: ``docs/research/route-b-index-primers.md``.
WELL_PRIMERS: tuple[str, ...] = (
    "ATCTCGATCCCGCGAAATTAATACGACTCAC",
    "GCCCCAAGGGGTTATGCTAGTTATTGCTC",
)


@dataclass(frozen=True, slots=True)
class Judgement:
    """One check, the molecule it judged, and the findings behind it.

    Parameters
    ----------
    check
        The verdict and the value it judged.
    where
        What was judged, in the method's words: a cargo, a destination, one tube.
    findings
        What mbio returned where the gate looked, empty where a check found nothing to point at.
    """

    check: Check
    where: str
    findings: tuple[Finding, ...] = ()

    @property
    def name(self) -> str:
        """What was measured."""
        return self.check.name

    @property
    def status(self) -> Status | None:
        """The verdict, or ``None`` where no threshold judges it."""
        return self.check.status


@dataclass(frozen=True, slots=True)
class Verdict:
    """What the gate hands back: its judgements, and the worst verdict among them.

    Parameters
    ----------
    judgements
        One per check made, in the order the gate made them.
    """

    judgements: Sequence[Judgement]

    def __post_init__(self) -> None:
        """Fix the order the judgements were made in."""
        object.__setattr__(self, "judgements", tuple(self.judgements))

    @property
    def checks(self) -> tuple[Check, ...]:
        """Every check, for a protocol to print."""
        return tuple(one.check for one in self.judgements)

    @property
    def findings(self) -> tuple[Finding, ...]:
        """Every finding behind those checks, for a map to draw."""
        return tuple(finding for one in self.judgements for finding in one.findings)

    @property
    def status(self) -> Status:
        """The worst verdict any check carries."""
        return worst_of(self.checks)

    @property
    def failures(self) -> tuple[Judgement, ...]:
        """The judgements that failed, which is what a report leads with."""
        return tuple(one for one in self.judgements if one.status == "fail")

    @property
    def summary(self) -> tuple[Check, ...]:
        """One check a name, carrying the worst judgement of that name and how many there were.

        A design is judged once a molecule and read once, so a page that prints a badge a check
        prints this rather than `checks`.
        """
        grouped: dict[str, list[Check]] = {}
        for one in self.judgements:
            grouped.setdefault(one.name, []).append(one.check)
        return tuple(_worst(name, group) for name, group in grouped.items())

    def __getitem__(self, name: str) -> Judgement:
        """Return the first judgement of that check name.

        Raises
        ------
        KeyError
            If no check has it.
        """
        for one in self.judgements:
            if one.name == name:
                return one
        raise KeyError(name)

    def __len__(self) -> int:
        """How many checks were made."""
        return len(self.judgements)


def check_cargo(
    cargo: SequenceRecord,
    *,
    project: Project,
    where: str = "cargo",
    ends_chain: bool = False,
    annealing: Sequence[str] = WELL_PRIMERS,
) -> tuple[Judgement, ...]:
    """Judge one cargo: the sites it must not carry, the frame it keeps, and what it must not read.

    A cargo is what the external enzyme releases from a synthesised block, counted from the
    first base of the overhang a part enters on. Every enzyme the project reserves is one no
    block spells anywhere, its stuffers included, because the step that reserved it cuts the
    cargo in a tube of its own. A cargo another part follows is whole codons; the one that ends
    the chain is one base past them, which is the method's capping block.

    A cargo also gives the primers that read a well nowhere to bind. One that carries an
    annealing region hands its well a second priming site, and that well's read is no longer one
    anybody can call: ``docs/research/route-b-index-primers.md`` §2.7.

    Parameters
    ----------
    cargo
        The released molecule, its entry overhang first.
    project
        What the build chose, which carries the method and the reserved enzymes.
    where
        What to call this cargo in a message.
    ends_chain
        Whether no further part follows this one.
    annealing
        What each primer that reads a well anneals to, 5' to 3'.
    """
    held = project.reserved_enzymes
    found = find_sites(cargo, held)
    named = ", ".join(sorted({site.enzyme.name for site in found}))
    sites = Judgement(
        Check(
            "cargo sites",
            "pass" if not found else "fail",
            len(found),
            f"{where} carries no site of {', '.join(one.name for one in held)}"
            if not found
            else f"{where} carries {len(found)} site(s) it must not, of {named}",
        ),
        where,
        tuple(found),
    )
    wanted = 1 if ends_chain else 0
    over = len(cargo) % 3
    kind = "the block that ends the chain" if ends_chain else "a cargo another part follows"
    frame = Judgement(
        Check(
            "cargo frame",
            "pass" if over == wanted else "fail",
            over,
            f"{where} is {len(cargo)} bases from its entry overhang, {over} past a whole number "
            f"of codons, which is what {kind} is",
        ),
        where,
    )
    bound = tuple(
        Segment(site.start, site.end)
        for region in annealing
        for site in find_binding_sites(region, cargo)
    )
    primed = Judgement(
        Check(
            "well primers",
            "pass" if not bound else "fail",
            len(bound),
            f"no primer that reads a well binds {where}"
            if not bound
            else f"a primer that reads a well binds {where} in {len(bound)} place(s), which "
            "leaves that well a second priming site and no read anyone can call",
        ),
        where,
        bound,
    )
    return (sites, frame, primed)


def check_barcode_set(
    codes: Sequence[str], *, project: Project, where: str = "this part list"
) -> tuple[Judgement, ...]:
    """Judge one part list's barcodes: how far apart they stand, and what each one reads as.

    The rules are the project's length and distance with the method's cloning scar, the frame the
    finished block reads the barcode in, and the enzymes a block is kept clear of.
    `liulab_mbio.barcodes` holds all of them; this names them for this method.
    """
    rules = barcode_rules(
        project.scheme,
        project.barcode.length,
        distance=project.barcode.min_distance,
        reserved=project.reserved_extra,
    )
    apart = separation(codes, metric=rules.metric) if len(codes) > 1 else rules.distance
    spacing = Judgement(
        Check(
            "barcode spacing",
            "pass" if apart >= rules.distance else "fail",
            apart,
            f"the closest two barcodes of {where} stand {apart} apart, "
            f"{'at' if apart == rules.distance else 'over' if apart > rules.distance else 'under'}"
            f" the {rules.distance} two parts of one list need to be told apart",
        ),
        where,
    )
    broken = check_barcodes(codes, rules)
    reading = Judgement(
        Check(
            "barcode reading",
            "pass" if not broken else "fail",
            len(broken),
            f"every barcode of {where} reads as the method needs it to"
            if not broken
            else f"{len(broken)} barcode rule(s) broken in {where}: {'; '.join(broken[:3])}",
        ),
        where,
    )
    return (spacing, reading)


def check_product(
    product: SequenceRecord,
    *,
    project: Project,
    barcodes: Mapping[str, Sequence[str]],
    where: str = "the library product",
) -> tuple[Judgement, ...]:
    """Judge the finished product: whether it opens, its barcode block, and what it keeps.

    The last round's product is the one molecule no round opens, so the chain of reactions never
    judges it as a destination; it is judged here instead, on the same rule. The block is read off
    the record rather than stated: the terminal stuffer is found, then one barcode a round follows
    it joined by the cloning scar, in the reverse of the order the rounds ran. What the product
    keeps past its last part -- that stuffer and that block -- is then read as codons, because
    every library member translates through it.
    """
    opened = _opening(product, project, where)
    bases = str(product.sequence)
    circular = product.topology == "circular"
    haystack = bases * 2 if circular else bases
    core = project.scheme.internal_stuffer_core
    at = haystack.find(core)
    opens = at - len(project.scheme.internal_stuffer_prefix)
    if at < 0 or opens < 0:
        return (
            opened,
            Judgement(
                Check(
                    "barcode block",
                    "fail",
                    0,
                    f"{where} carries no terminal stuffer, so the next round could not open it "
                    "and nothing says where its barcode block begins",
                ),
                where,
            ),
        )
    read, missing, ends = _block(haystack, at + len(core), project, barcodes)
    block = Judgement(
        Check(
            "barcode block",
            "pass" if not missing else "fail",
            len(read),
            f"{where} keeps one barcode a round, {' then '.join(read)}, "
            "in the reverse of the order the rounds ran"
            if not missing
            else f"{where} keeps {len(read)} of the {project.position_count} barcodes a round "
            f"leaves: {missing}",
        ),
        where,
    )
    start = opens % len(product) if circular else opens
    span = Segment(start, start + ends - opens)
    over = (span.end - span.start) % 3
    terminal = Judgement(
        Check(
            "terminal block",
            "pass" if over == 0 else "fail",
            over,
            f"the {span.end - span.start} bases {where} keeps past its last part are "
            f"{'whole codons' if over == 0 else f'{over} past a whole number of codons'}",
        ),
        where,
    )
    stops = stop_codons(product, span) if over == 0 else ()
    reading = Judgement(
        Check(
            "terminal stop",
            "pass" if not stops else "fail",
            len(stops),
            f"{where} reads through what it keeps without a stop"
            if not stops
            else f"{len(stops)} stop codon(s) in what {where} keeps past its last part",
        ),
        where,
        tuple(stops),
    )
    return (opened, block, terminal, reading)


def check_dmx_vector(
    vector: SequenceRecord,
    *,
    project: Project,
    annealing: Sequence[str] = WELL_PRIMERS,
    where: str = "the DMX vector",
) -> tuple[Judgement, ...]:
    """Judge where the DMX vector's blunt sites sit: clear of the primers that read a well.

    The method asks for one blunt site outboard of each site the external enzyme reads, and says
    no more. Reading a well amplifies the cassette whole with one primer pair, so a blunt site in
    either primer's footprint leaves that well no product and no verdict either, which no other
    check would say. The bases such a site may take are therefore read off this record, never
    stated: each primer's footprint and the releasing cut it faces bound them.

    Only the blunt sites beside the cassette are judged. One inside the cargo is the part's own,
    and one elsewhere in the backbone reaches no primer.

    Parameters
    ----------
    vector
        The rebuilt DMX vector, as it is held.
    project
        What the build chose, which carries the method's enzymes.
    annealing
        What each primer that reads a well anneals to, 5' to 3'.
    where
        What to call this vector in a message.
    """
    external = project.scheme.external
    flanks, missing = _flanks(vector, external, annealing)
    if missing:
        return (
            Judgement(
                Check(
                    "blunt sites",
                    "fail",
                    0,
                    f"nothing in {where} says where a blunt site may sit: {missing}",
                ),
                where,
            ),
        )
    judged: list[CutSite] = []
    stray: list[tuple[CutSite, Segment]] = []
    for site in find_sites(vector, project.scheme.blunt):
        for flank, clear in flanks:
            if not _meets(vector, flank, site.span):
                continue
            judged.append(site)
            if not vector.covers(clear, site.span):
                stray.append((site, clear))
            break
    spans = " or ".join(_printed(vector, clear) for _, clear in flanks)
    named = ", ".join(
        f"{one.enzyme.name} at {_printed(vector, one.span)}, outside {_printed(vector, clear)}"
        for one, clear in stray
    )
    return (
        Judgement(
            Check(
                "blunt sites",
                "pass" if not stray else "fail",
                len(judged),
                f"the {len(judged)} blunt site(s) beside {where}'s cassette sit between a "
                f"{external.name} site and the primer that reads a well, in {spans}"
                if not stray
                else f"{len(stray)} blunt site(s) of {where} miss the bases a {external.name} "
                f"site and the primer that reads a well leave clear: {named}",
            ),
            where,
            tuple(one for one, _ in stray),
        ),
    )


def check_reaction(reaction: Reaction, *, project: Project) -> tuple[Judgement, ...]:
    """Judge one tube: whether every molecule in it is cut where this method cuts it.

    In a digest, each molecule carries the acting enzymes' sites exactly where it is meant to be
    cut and nowhere else, which for this method is two cuts: a destination gives up its stuffer
    and a block gives up its cargo. In a ligation, the ends meeting in the tube are mutually
    distinguishable, which for this method is the entry overhang and the cloning scar and nothing
    else.
    """
    if reaction.kind == "ligation":
        return _ligation(reaction, project=project)
    return tuple(_cutting(reaction, pool) for pool in reaction.pools)


def library_reactions(
    project: Project,
    *,
    destination: SequenceRecord,
    blocks: Mapping[str, Sequence[SequenceRecord]],
    products: Sequence[SequenceRecord],
) -> tuple[Reaction, ...]:
    """Compose this method's reactions for one build: two digests and a ligation a round.

    A round opens the library built so far in one tube and releases its parts from their blocks
    in another, then joins the two in a third. One round's product is the next one's destination,
    and the chain records that order and nothing else.

    Raises
    ------
    ValueError
        If there is not one product and one part list a position.
    """
    if len(products) != project.position_count:
        raise ValueError(
            f"this build fills {project.position_count} position(s) and {len(products)} "
            "product(s) were given, one a round"
        )
    made: list[Reaction] = []
    standing = destination
    for number, position in enumerate(project.positions, start=1):
        if position not in blocks:
            raise ValueError(f"no block was given for position {position!r}")
        donors = tuple(blocks[position])
        opened = Pool("destination", (standing,), enzyme=project.scheme.internal)
        donated = Pool("donor", donors, enzyme=project.scheme.external)
        made.extend(
            (
                Reaction(
                    f"round {number} opening",
                    "digest",
                    pools=(opened,),
                    enzymes=(project.scheme.internal,),
                ),
                Reaction(
                    f"round {number} release",
                    "digest",
                    pools=(donated,),
                    enzymes=(project.scheme.external,),
                ),
                Reaction(f"round {number} ligation", "ligation", pools=(opened, donated)),
            )
        )
        standing = products[number - 1]
    return tuple(made)


def final_assembly_reactions(
    project: Project, *, library: SequenceRecord, working: SequenceRecord, cargo: Enzyme
) -> tuple[Reaction, ...]:
    """Compose the final assembly: two digests and the ligation that joins what they leave.

    One tube, staged. The external enzyme frees the cargo from the backbone the rounds ran in and
    is killed; the working vector then gives its ccdB cassette up to `cargo`, and the ligase joins
    the two. The gate reads it as three reactions because the ends meeting at the end were made by
    two different enzymes, which is the one thing a round never does.
    """
    freed = Pool("donor", (library,), enzyme=project.scheme.external)
    opened = Pool("destination", (working,), enzyme=cargo)
    return (
        Reaction(
            "final assembly release", "digest", pools=(freed,), enzymes=(project.scheme.external,)
        ),
        Reaction("final assembly opening", "digest", pools=(opened,), enzymes=(cargo,)),
        Reaction("final assembly ligation", "ligation", pools=(opened, freed)),
    )


def check_library(
    project: Project,
    *,
    destination: SequenceRecord,
    blocks: Mapping[str, Sequence[SequenceRecord]],
    products: Sequence[SequenceRecord],
    barcodes: Mapping[str, Sequence[str]],
    working: SequenceRecord | None = None,
    cargo: Enzyme | None = None,
) -> Verdict:
    """Judge a whole build: its vector, every tube, every cargo, every barcode set and the product.

    Nothing here is read off a plan. The molecules are finished records, the barcodes are the
    table that decodes the sequencing, and the reactions are composed from the chain's own order,
    so a design an agent wrote is judged exactly as one `plan_igga` wrote is.

    The rounds run in the DMX vector, so the destination is where a build accepts one and
    `check_dmx_vector` judges it here. A destination no primer that reads a well binds carries no
    DMX backbone -- a minimal stand-in -- and has no blunt site any primer could reach, so
    nothing judges it.

    Parameters
    ----------
    project
        What the build chose.
    destination
        The vector the first round opens, judged as a rebuilt DMX vector where one is given.
    blocks
        The synthesised blocks of each position, keyed by position name.
    products
        One record a round, in the order the rounds ran. The last is the library product.
    barcodes
        Each position's barcodes, keyed by position name.
    working, cargo
        The working vector the library is moved into and the enzyme that admits it, where the
        build names one. Both or neither: the final assembly is judged only where there is a
        vector to judge it in.

    Raises
    ------
    ValueError
        If there is not one product, one part list and one barcode list a position, or if one of
        `working` and `cargo` is given without the other.
    """
    if (working is None) != (cargo is None):
        raise ValueError(
            "a final assembly is judged on a working vector and the enzyme that admits cargo to "
            "it: name both, or neither"
        )
    made: list[Judgement] = []
    if _read_by_a_well_primer(destination, WELL_PRIMERS):
        made.extend(check_dmx_vector(destination, project=project))
    reactions = library_reactions(
        project, destination=destination, blocks=blocks, products=products
    )
    if working is not None and cargo is not None:
        reactions += final_assembly_reactions(
            project, library=products[-1], working=working, cargo=cargo
        )
    for reaction in reactions:
        made.extend(check_reaction(reaction, project=project))
    last = project.positions[-1]
    for position in project.positions:
        if position not in barcodes:
            raise ValueError(f"no barcode list was given for position {position!r}")
        made.extend(check_barcode_set(barcodes[position], project=project, where=position))
        for index, block in enumerate(blocks[position]):
            made.extend(
                check_cargo(
                    _cargo(block, project),
                    project=project,
                    where=f"the {position} cargo {index + 1}",
                    ends_chain=position == last,
                )
            )
    made.extend(check_product(products[-1], project=project, barcodes=barcodes))
    return Verdict(made)


def _opening(product: SequenceRecord, project: Project, where: str) -> Judgement:
    """Whether the internal enzyme still cuts the finished product where a further round opens it.

    The method leaves the library openable after its last round, which is what a further round or
    a transfer into a working vector reads.
    """
    enzyme = project.scheme.internal
    found = find_sites(product, (enzyme,))
    opens = len(found) == CUTS
    return Judgement(
        Check(
            "product opens",
            "pass" if opens else "fail",
            len(found),
            f"{enzyme.name} cuts {where} in the {CUTS} places a further round opens it on"
            if opens
            else f"{enzyme.name} reads {len(found)} site(s) in {where}, where a further round "
            f"opens it on {CUTS}",
        ),
        where,
        () if opens else tuple(found),
    )


def _read_by_a_well_primer(record: SequenceRecord, annealing: Sequence[str]) -> bool:
    """Whether a primer that reads a well binds this record at all.

    Both bind a rebuilt DMX vector, and neither binds a stand-in carrying a cassette and nothing
    else. Where one binds and the other does not, the rebuild broke a primer site and
    `check_dmx_vector` says so rather than passing.
    """
    return any(find_binding_sites(region, record) for region in annealing)


def _worst(name: str, group: Sequence[Check]) -> Check:
    """Return the worst check of one name, saying how many carried it where more than one did."""
    kept = max(group, key=lambda one: STATUSES.index(one.status) if one.status else -1)
    detail = kept.detail if len(group) == 1 else f"{len(group)} judged, worst: {kept.detail}"
    return Check(name, kept.status, kept.value, detail)


def _cargo(block: SequenceRecord, project: Project) -> SequenceRecord:
    """Return what the external enzyme releases from one block, its entry overhang first.

    Raises
    ------
    ValueError
        If the block does not give up exactly one piece with an overhang at each end.
    """
    released = _released(block, project.scheme.external)
    if len(released) != 1:
        raise ValueError(
            f"{project.scheme.external.name} releases {len(released)} cargo from this block, "
            "and a block carries one"
        )
    return SequenceRecord(block.bases(released[0].start, released[0].end))


def _released(record: SequenceRecord, enzyme: Enzyme) -> tuple[Fragment, ...]:
    """Return the pieces of a digest with an overhang at each end: what a ligation takes."""
    return tuple(
        piece for piece in digest(record, enzyme) if piece.left_overhang and piece.right_overhang
    )


def _cutting(reaction: Reaction, pool: Pool) -> Judgement:
    """Whether every molecule of one pool is cut in the places this method cuts it."""
    acting = reaction.acting
    offenders: list[SequenceRecord] = []
    findings: list[Finding] = []
    for record in pool:
        found = find_sites(record, acting)
        if len(found) != CUTS:
            offenders.append(record)
            findings.extend(found)
    named = ", ".join(one.name for one in acting)
    verb = "opens" if pool.role == "destination" else "releases"
    return Judgement(
        Check(
            f"{pool.role} {verb}",
            "pass" if not offenders else "fail",
            len(offenders),
            f"{named} cuts each of the {len(pool)} {pool.role}(s) of {reaction.name} in "
            f"{CUTS} places"
            if not offenders
            else f"{len(offenders)} of the {len(pool)} {pool.role}(s) of {reaction.name} are not "
            f"cut by {named} in the {CUTS} places this method cuts them",
        ),
        reaction.name,
        tuple(findings),
    )


def _ligation(reaction: Reaction, *, project: Project) -> tuple[Judgement, ...]:
    """Whether the ends meeting in one tube can be told apart, and how well they ligate."""
    ends = sorted(_ends(reaction, project))
    distinguishable = len(ends) == LIGATION_OVERHANGS and not any(
        one == reverse_complement(other) for one in ends for other in ends
    )
    told = Judgement(
        Check(
            "ends distinguishable",
            "pass" if distinguishable else "fail",
            len(ends),
            f"the ends meeting in {reaction.name} read as {', '.join(ends)}, which a ligation "
            "pairs one way"
            if distinguishable
            else f"the ends meeting in {reaction.name} read as "
            f"{', '.join(ends) or 'nothing at all'}, where this method joins on "
            f"{LIGATION_OVERHANGS} overhangs no two of which anneal",
        ),
        reaction.name,
    )
    if not ends:
        return (told,)
    report = fidelity(ends, project.scheme.internal)
    return (
        told,
        Judgement(
            Check(
                "ligation fidelity",
                None,
                report.value,
                f"{reaction.name} scores {report.value:.3f} on {report.source}",
            ),
            reaction.name,
        ),
    )


def _ends(reaction: Reaction, project: Project) -> set[str]:
    """Every overhang the molecules of a ligation present, read off the digest that made it.

    A pool naming its own cutter is read off that one. The method's round is what the rest fall
    back to, where the internal enzyme opens the destination and the external releases everything
    else; a final assembly does not follow it, because the working vector gives its cassette up to
    the cargo enzyme instead.
    """
    cutters: dict[Role, Enzyme] = {
        "destination": project.scheme.internal,
        "donor": project.scheme.external,
        "insert": project.scheme.external,
        "carrier": project.scheme.external,
    }
    found: set[str] = set()
    for pool in reaction.pools:
        enzyme = pool.cutter or cutters[pool.role]
        for record in pool:
            for piece in _released(record, enzyme):
                found.update((piece.left_overhang, piece.right_overhang))
    return found


def _block(
    bases: str, at: int, project: Project, barcodes: Mapping[str, Sequence[str]]
) -> tuple[tuple[str, ...], str, int]:
    """Walk the barcode block from `at`, one barcode a round joined by the cloning scar.

    Returns what was read, a phrase saying where the walk stopped -- empty where it read one
    barcode for every position -- and where it stopped.
    """
    scar = project.scheme.cloning_scar
    length = project.barcode.length
    read: list[str] = []
    for position in reversed(project.positions):
        code = bases[at : at + length]
        if code not in barcodes.get(position, ()):
            return tuple(read), f"{code!r} names no part of position {position}", at
        read.append(code)
        at += length
        if len(read) < project.position_count:
            if bases[at : at + len(scar)] != scar:
                return tuple(read), f"no cloning scar follows the {position} barcode", at
            at += len(scar)
    return tuple(read), "", at


def _flanks(
    record: SequenceRecord, external: Enzyme, annealing: Sequence[str]
) -> tuple[tuple[tuple[Segment, Segment], ...], str]:
    """Where a blunt site may sit beside the cassette, read off the record: one pair a primer.

    The first span of a pair runs from the primer's footprint to the nearest releasing cut it
    faces, so a blunt site there is one placed against that cut; the second is what the first
    leaves clear of the footprint itself. A footprint is the whole annealing region, counted
    back from where its 3' end binds, so bases a site already broke are inside it. The phrase
    says what stopped the reading, empty where nothing did.
    """
    opened = find_sites(record, external)
    found: list[tuple[Segment, Segment]] = []
    for region in annealing:
        sites = find_binding_sites(region, record)
        if len(sites) != 1:
            return (), (
                f"a primer that reads a well binds it in {len(sites)} place(s), where one "
                "amplification binds it in 1"
            )
        site = sites[0]
        forward = site.strand is Strand.FORWARD
        best: tuple[Segment, Segment] | None = None
        for cut in opened:
            clear = _span(record, *((site.end, cut.start) if forward else (cut.end, site.start)))
            if clear is None or (best and clear.end - clear.start >= best[1].end - best[1].start):
                continue
            flank = (
                _span(record, site.end - len(region), clear.end)
                if forward
                else _span(record, clear.start, site.start + len(region))
            )
            if flank is not None:
                best = (flank, clear)
        if best is None:
            return (), f"no {external.name} site faces the primer binding at {site.start}"
        found.append(best)
    return tuple(found), ""


def _span(record: SequenceRecord, start: int, end: int) -> Segment | None:
    """Return the bases from `start` to `end` as the top strand reads them, or ``None``.

    There are none where the two meet. A span across the origin ends past the record's length.
    """
    length = len(record)
    if record.topology == "circular":
        start %= length
        reach = (end - start) % length
    else:
        start = max(start, 0)
        reach = end - start
    return Segment(start, start + reach) if reach > 0 else None


def _meets(record: SequenceRecord, span: Segment, other: Segment) -> bool:
    """Whether two spans of one record share a base, across the origin too."""
    return record.covers(span, other.start) or record.covers(other, span.start)


def _printed(record: SequenceRecord, span: Segment) -> str:
    """One span as a person reads it: 1-based, inclusive, and counted round the origin."""
    return f"{span.start + 1}..{(span.end - 1) % len(record) + 1}"

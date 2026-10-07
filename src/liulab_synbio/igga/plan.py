"""One combinatorial library, planned from a project file and nothing else.

`plan_igga` is the one way in. It reads the project, sorts the part lists, chooses the overhang
standard the proteins cost least, builds every part's synthesis sequence, makes the vector a
destination, simulates every round, works out what each round takes at the bench and how many
colonies it needs, and hands the finished records to `liulab_synbio.igga.gate` to be judged.
`LibraryPlan.write` puts one directory's worth of output in one place: the
synthesis order sheet, the barcode table, the amino-acid change table, a record for every round,
the protocol as JSON data, and the page rendered from that data.

The method is `liulab_synbio.igga.method.IGGA` and is not an argument. What a build chooses is
the project's, and `docs/adr/0010-method-in-code.md` draws the line between them.

**The vector and the standard have to agree about position one.** A vector already carrying an
internal stuffer spells an entry overhang in DNA that exists, so position one's overhang is
pinned to it and the standard designs around that. A vector that has to be retrofitted has its
stuffer synthesised now, so the standard chooses position one freely and the stuffer put in
carries what it chose. Nothing else in the method moves either way.
"""

import dataclasses
import os
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from liulab_mbio.barcodes import BarcodeRules
from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.pools import pool_sheet, primer_inventory
from liulab_mbio.bench.prices import PriceRecord, read_prices
from liulab_mbio.checks import Check, Status
from liulab_mbio.cloning.plan import as_record, status, write_protocol_files
from liulab_mbio.codons import codon_usage
from liulab_mbio.overhangs import MIN_DISTANCE
from liulab_mbio.protocol.model import Protocol
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.sites import digest
from liulab_mbio.translate import translate
from liulab_synbio import dmx
from liulab_synbio.igga.bench import digest_amount, ligation_amounts, transformation_amount
from liulab_synbio.igga.cargo import PoolPlan, design_pool, read_bands, read_primers
from liulab_synbio.igga.coverage import RoundCoverage, constructs, plan_coverage
from liulab_synbio.igga.gate import Verdict, check_library
from liulab_synbio.igga.method import Scheme
from liulab_synbio.igga.parts import (
    Part,
    barcode_rules,
    barcode_table,
    change_table,
    design_parts,
    synthesis_sheet,
)
from liulab_synbio.igga.project import Project, read_project
from liulab_synbio.igga.rounds import Round, assemble_rounds, representative, write_records
from liulab_synbio.igga.standard import PartList, Standard, design_standard
from liulab_synbio.igga.steps import RoundBench
from liulab_synbio.igga.steps import protocol as protocol_for
from liulab_synbio.igga.vector import Destination, Site, destination_vector

#: What a part list holds: the proteins each member codes for, or the DNA it is already coded in.
type Kind = Literal["protein", "dna"]

#: How a record's name says which part list it belongs to. ``{position}`` stands for a position's
#: name, and what is around it holds that name to a word boundary, so ``N_ATF2`` names position
#: ``N`` and ``C_NFAT`` does not. A caller whose names read another way passes another pattern.
NAME_PATTERN = r"(?<![A-Za-z0-9]){position}(?![A-Za-z0-9])"

#: What `LibraryPlan.write` calls the sheets it writes. The records are named by
#: `liulab_synbio.igga.rounds`, which writes one for each round and the product for the
#: last, and the protocol pair by `liulab_mbio.cloning.plan`.
PARTS_FILE = "parts.tsv"
BARCODE_FILE = "barcodes.tsv"
CHANGE_FILE = "changes.tsv"
POOL_FILE = "pool.tsv"
POOL_PRIMER_FILE = "pool-primers.tsv"


@dataclass(frozen=True, slots=True)
class Files:
    """The files a library plan writes.

    Parameters
    ----------
    parts
        The synthesis order sheet: every block to order, with its barcode on the same row.
    barcodes
        Which barcode names which part, and where it sits in the finished block.
    changes
        Every amino acid the overhang standard moved, wild type beside synthesised.
    records
        One annotated record a round, the last of them the representative construct.
    protocol_data
        The bench protocol as JSON, which ``protocol render`` turns back into a page.
    protocol
        The interactive bench protocol, as one self-contained HTML page rendered from
        `protocol_data`.
    pool, pool_primers
        The oligo pool to order and the primers that amplify it. Both are ``None`` where the
        project names no primer set, because the primer sites are templated on the oligo and
        nothing can be written without them.
    """

    parts: Path
    barcodes: Path
    changes: Path
    records: tuple[Path, ...]
    protocol_data: Path
    protocol: Path
    pool: Path | None = None
    pool_primers: Path | None = None

    @property
    def paths(self) -> tuple[Path, ...]:
        """Every file, in the order they were written: the sheets, the records, the protocol."""
        return (
            self.parts,
            self.barcodes,
            self.changes,
            *self.records,
            self.protocol_data,
            self.protocol,
        )


@dataclass(frozen=True, slots=True)
class LibraryPlan:
    """One planned library: what it costs, what it makes, and what the bench has to do.

    Parameters
    ----------
    project
        What this build chose: its positions, its inputs and its dials.
    scheme
        The method the build was planned by.
    vector
        The record the plan was made from, before any retrofit.
    destination
        The vector a round can open, and what making it one changed.
    part_lists
        One per position, in the project's order: each member's name and the protein it codes for,
        whether it was given as protein or read off DNA.
    coding
        One per position: the coding sequence a member came already coded in, empty for a member
        written from the host's codon usage.
    standard
        The overhang standard the build works to, and what it charges each part.
    parts
        Every part's synthesis sequence, in position order and then list order.
    rounds
        One a position, in the order they run. The last round's product is the representative
        construct and the ones before it are the intermediates.
    coverage
        What each round has to cover, and what the colonies asked for leave out.
    bench
        What each round takes at the bench, computed from that round's own lengths.
    verdict
        What `liulab_synbio.igga.gate` made of the finished design: every tube, every cargo,
        every barcode set and the product.
    host
        The codon usage table the coding bases were written for.
    name
        What each round's product is called.
    prices
        The price record the protocol's bill is costed against, if the caller holds one. Its
        quantities compute either way; with no record every money cell is a hole.
    pool
        The oligo pool every block is synthesised from, where the project names a primer set.
        `liulab_synbio.igga.cargo` designs it.
    """

    project: Project
    scheme: Scheme
    vector: SequenceRecord
    destination: Destination
    part_lists: tuple[PartList, ...]
    coding: tuple[Mapping[str, str], ...]
    standard: Standard
    parts: tuple[Part, ...]
    rounds: tuple[Round, ...]
    coverage: tuple[RoundCoverage, ...]
    bench: tuple[RoundBench, ...]
    verdict: Verdict
    host: str
    name: str = ""
    prices: PriceRecord | None = None
    pool: PoolPlan | None = None

    @property
    def product(self) -> SequenceRecord:
        """The representative construct: the last round's product."""
        return self.rounds[-1].product

    @property
    def constructs(self) -> int:
        """How many distinct constructs the part lists yield."""
        return constructs([len(one) for one in self.part_lists])

    @property
    def representative_parts(self) -> tuple[Part, ...]:
        """The one part a position the records were simulated from."""
        return tuple(one.part for one in self.rounds)

    @property
    def validation(self) -> dmx.Validation | None:
        """What reading these designs back takes, or `None` where the project reads none.

        The project's floor chooses the designs and its route reads them. The bench is sized
        from that set and not from the part list, so a design the floor leaves out costs no
        well, no plate and no reagent.
        """
        if self.project.route is None:
            return None
        return dmx.validation(
            dmx.ROUTES[self.project.route],
            designs(self.parts, self.project.oligo_length),
            self.project.validate_from,
        )

    @property
    def checks(self) -> tuple[Check, ...]:
        """Every check the gate made, in the order it made them.

        The plan judges nothing itself: these are `verdict`'s, so a design this pipeline wrote
        and one an agent composed are judged by the same code.
        """
        return self.verdict.checks

    @property
    def status(self) -> Status:
        """The worst status of any check."""
        return status(self.checks)

    def protocol(self) -> Protocol:
        """Return the bench protocol for this plan, covering every round as one experiment."""
        return protocol_for(
            scheme=self.scheme,
            positions=self.project.positions,
            barcode_length=self.project.barcode.length,
            vector=self.vector,
            destination=self.destination,
            part_lists=self.part_lists,
            standard=self.standard,
            parts=self.parts,
            rounds=self.rounds,
            bench=self.bench,
            constructs=self.constructs,
            checks=self.verdict.summary,
            host=self.host,
            sheet=PARTS_FILE,
            barcodes=BARCODE_FILE,
            validation=self.validation,
            prices=self.prices,
        )

    def write(self, directory: str | os.PathLike[str]) -> Files:
        """Write the sheets, the records, the protocol data and its page into `directory`.

        The directory is made when it is not there. The files are named by `PARTS_FILE`,
        `BARCODE_FILE` and `CHANGE_FILE`, by `liulab_synbio.igga.rounds` for the records
        and by `liulab_mbio.cloning.plan` for the protocol pair, and a second run over the
        same inputs writes the same bytes.
        """
        out = Path(directory)
        out.mkdir(parents=True, exist_ok=True)
        sheet = out / PARTS_FILE
        sheet.write_text(synthesis_sheet(self.parts), encoding="utf-8")
        barcodes = out / BARCODE_FILE
        barcodes.write_text(
            barcode_table(self.parts, self.project.position_count), encoding="utf-8"
        )
        changes = out / CHANGE_FILE
        changes.write_text(change_table(self.standard), encoding="utf-8")
        records = write_records(self.rounds, out)
        written = write_protocol_files(self.protocol(), out)
        pool = primers = None
        if self.pool is not None:
            pool = out / POOL_FILE
            pool.write_text(pool_sheet(self.pool.pool), encoding="utf-8")
            primers = out / POOL_PRIMER_FILE
            primers.write_text(primer_inventory(self.pool.pool), encoding="utf-8")
        return Files(sheet, barcodes, changes, records, written.data, written.page, pool, primers)


def designs(parts: Sequence[Part], oligo_length: int) -> tuple[dmx.Design, ...]:
    """Return one design a part, in as many pieces as the project's oligo length forces it into.

    A block no longer than one oligo of the pool is ordered whole; a longer one is ordered in
    pieces and joined, and that count is what a design's chance of a clean colony falls with.

    Raises
    ------
    ValueError
        If the oligo length is not positive.

    Examples
    --------
    >>> [one.fragments for one in designs(plan.parts, 350)]  # doctest: +SKIP
    [1, 2, 4]
    """
    if oligo_length < 1:
        raise ValueError(f"an oligo is at least one base, got {oligo_length}")
    return tuple(dmx.Design(one.name, -(-one.length // oligo_length)) for one in parts)


def plan_igga(
    project: Project | str | os.PathLike[str],
    *,
    parts: Sequence[Mapping[str, str]] | None = None,
    kind: Kind = "protein",
    site: Site | None = None,
    pattern: str = NAME_PATTERN,
    rules: BarcodeRules | None = None,
    min_distance: int = MIN_DISTANCE,
    allow_uniform: bool = False,
    prices: PriceRecord | str | os.PathLike[str] | None = None,
) -> LibraryPlan:
    """Plan the whole library `project` asks for, by the method `project` is built under.

    One round appends one part list to every member of the library at once, so the rounds run in
    the project's own order and the product of each opens the next. The same inputs return the
    same design: the barcodes are drawn from the project's seed, and nothing else here is random.

    Parameters
    ----------
    project
        What this build chooses, or a path to the JSON holding it;
        `liulab_synbio.igga.project.read_project` reads one. It names the parts FASTA and the
        vector by path.
    parts
        One already-sorted mapping per position, in the project's order. The project's own parts
        FASTA is read where this is not given.
    kind
        What the part lists hold. ``"dna"`` is checked and kept rather than written again; only a
        codon a junction or a forbidden site moves is this package's.
    site
        Where to put an internal stuffer, as a feature name or a ``(start, end)`` span. Read only
        where the vector carries none.
    pattern
        How a record's name says which part list it belongs to; see `NAME_PATTERN`.
    rules
        What every barcode holds to. Built from the project's own barcode length and distance by
        default; see `liulab_synbio.igga.parts.barcode_rules`.
    min_distance, allow_uniform
        How far apart the standard's overhangs must stand, and whether one base kind is allowed.
    prices
        A price record the user holds, or a path to one;
        `liulab_mbio.bench.prices.read_prices` reads one. The protocol's bill computes its
        quantities either way, and prices nothing without this.

    Returns
    -------
    LibraryPlan
        The design, the simulated rounds, what the bench has to do, and the gate's verdict on all
        of it. `liulab_synbio.igga.gate.check_library` judges the finished records, so this
        plan and a design an agent composed are judged by the same code.

    Raises
    ------
    ValueError
        If a record's name says no position of the project or says more than one, if two records
        share a name, if a part is not a protein or not a coding sequence, if no overhang
        standard fits the part lists, if a block spells a site the method does not expect or
        gives up no cargo, if the vector cannot be made a destination, if a round cannot ligate,
        or if `prices` names a file that is not a price record.
    KeyError
        If the project names a codon usage table or an enzyme this package does not ship.
    liulab_mbio.barcodes.SpaceExhaustedError
        If a part list is larger than the barcodes its rules allow.

    Examples
    --------
    >>> plan = plan_igga("project.json")  # doctest: +SKIP
    >>> plan.write("library/")  # doctest: +SKIP
    """
    chosen = project if isinstance(project, Project) else read_project(project)
    design = chosen.scheme
    positions = chosen.positions
    one = as_record(chosen.vector)
    given = (
        parts if parts is not None else read_part_lists(chosen.parts, positions, pattern=pattern)
    )
    lists, coded = _sequences(_checked(given, positions), kind)
    compatible = _compatible(one, design)
    pinned = {positions[0]: design.entry_overhang} if compatible else {}
    standard = design_standard(
        design,
        positions,
        lists,
        pinned=pinned,
        usage=codon_usage(chosen.host),
        min_distance=min_distance,
        allow_uniform=allow_uniform,
    )
    held = (
        rules
        if rules is not None
        else barcode_rules(
            design,
            chosen.barcode.length,
            distance=chosen.barcode.min_distance,
            reserved=chosen.reserved_extra,
        )
    )
    built = design_parts(
        design,
        positions,
        lists,
        standard,
        host=chosen.host,
        rules=held,
        coding=coded,
        seed=chosen.seed,
        reserved=chosen.reserved_extra,
    )
    destination = destination_vector(
        one, design if compatible else _restandardised(design, standard), site=site
    )
    named = chosen.name or one.name
    rounds = assemble_rounds(
        destination.record, representative(built, positions), design, positions, name=named
    )
    rows = plan_coverage([len(each) for each in lists], coverage=chosen.coverage)
    judged = check_library(
        chosen,
        destination=destination.record,
        blocks=_blocks(built, positions),
        products=[one.product for one in rounds],
        barcodes=_barcodes(built, positions),
    )
    return LibraryPlan(
        chosen,
        design,
        one,
        destination,
        tuple(lists),
        tuple(coded) if coded is not None else tuple({} for _ in lists),
        standard,
        built,
        rounds,
        rows,
        _bench(rounds, built, rows),
        judged,
        chosen.host,
        named,
        prices if prices is None or isinstance(prices, PriceRecord) else read_prices(prices),
        _pool(chosen, built),
    )


def _pool(project: Project, parts: Sequence[Part]) -> PoolPlan | None:
    """Design the oligo pool, or none where the project names no primer set.

    The sites that cut a fragment out are templated on the oligo rather than carried by a
    primer, so a pool cannot be written at all without the set. A project naming none still
    plans every other output.
    """
    if project.primers is None:
        return None
    return design_pool(
        parts,
        project,
        primers=read_primers(project.primers),
        bands=read_bands(project.bands),
        seed=project.seed,
    )


def read_part_lists(
    path: str | os.PathLike[str], positions: Sequence[str], *, pattern: str = NAME_PATTERN
) -> tuple[dict[str, str], ...]:
    """Read one FASTA and sort its records into one part list a position, in the project's order.

    A record is ordered under the name the FASTA gives it, and that name is what says which
    position it fills.

    Raises
    ------
    ValueError
        If the file holds no record, if two records share a name, if a name says no position or
        says more than one, or if a position ends up with nothing to fill it.

    Examples
    --------
    >>> read_part_lists("parts.fasta", ("N", "DBD", "C"))  # doctest: +SKIP
    ({'N_ATF2': 'MKT...'}, {'DBD_JUN': 'WQA...'}, {'C_VP64': 'MKT...'})
    """
    return _sorted(_fasta(path), positions, pattern)


def _fasta(path: str | os.PathLike[str]) -> dict[str, str]:
    """Read every record of a FASTA as its name and the letters it spells."""
    from Bio import SeqIO

    found: dict[str, str] = {}
    # Opened here rather than by name, so the handle is closed rather than left to the collector.
    with Path(path).open(encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            named = str(record.id)
            if named in found:
                raise ValueError(
                    f"{os.fspath(path)} names {named!r} twice: every part is ordered under its "
                    "own name"
                )
            found[named] = str(record.seq)
    if not found:
        raise ValueError(f"{os.fspath(path)} holds no FASTA record")
    return found


def _sorted(
    records: Mapping[str, str], positions: Sequence[str], pattern: str
) -> tuple[dict[str, str], ...]:
    """Put each record in the one part list its name names.

    Raises
    ------
    ValueError
        If a name says no position or says more than one, or a position is left empty.
    """
    names = list(positions)
    matchers = [
        re.compile(pattern.replace("{position}", re.escape(one)), re.IGNORECASE) for one in names
    ]
    lists: list[dict[str, str]] = [{} for _ in names]
    for named, sequence in records.items():
        hit = [index for index, matcher in enumerate(matchers) if matcher.search(named)]
        if not hit:
            raise ValueError(
                f"the name {named!r} says no position of this project, whose positions are "
                f"{', '.join(names)}. Rename the record, or pass a pattern that reads the "
                "names you already have"
            )
        if len(hit) > 1:
            said = ", ".join(names[index] for index in hit)
            raise ValueError(
                f"the name {named!r} says {len(hit)} positions of this project — {said} — so "
                "which part list it belongs to is ambiguous"
            )
        lists[hit[0]][named] = sequence
    if empty := [one for one, parts in zip(names, lists, strict=True) if not parts]:
        raise ValueError(
            f"no record names position(s) {', '.join(empty)}: every position needs a part list"
        )
    return tuple(lists)


def _checked(
    given: Sequence[Mapping[str, str]], positions: Sequence[str]
) -> tuple[Mapping[str, str], ...]:
    """Refuse part lists that are not one a position, or that name one part twice.

    Raises
    ------
    ValueError
        Naming the count, or the name two part lists share.
    """
    lists = tuple(given)
    if len(lists) != len(positions):
        said = ", ".join(positions)
        raise ValueError(
            f"this project has {len(positions)} position(s) — {said} — and "
            f"{len(lists)} part list(s) were given"
        )
    seen: set[str] = set()
    for parts in lists:
        for name in parts:
            if name in seen:
                raise ValueError(
                    f"{name!r} names a part in two part lists: every part is ordered under its "
                    "own name, so the order sheet says which is which"
                )
            seen.add(name)
    return lists


def _sequences(
    given: Sequence[Mapping[str, str]], kind: Kind
) -> tuple[tuple[PartList, ...], tuple[Mapping[str, str], ...] | None]:
    """Return the protein each part codes for, and the bases it came coded in where it did.

    Raises
    ------
    ValueError
        If `kind` is neither, or a part is not a coding sequence or spells a stop.
    """
    if kind == "protein":
        return tuple({name: one.upper() for name, one in parts.items()} for parts in given), None
    if kind != "dna":
        raise ValueError(f"kind is 'protein' or 'dna', got {kind!r}")
    lists: list[PartList] = []
    coded: list[Mapping[str, str]] = []
    for parts in given:
        lists.append({name: _protein(name, dna) for name, dna in parts.items()})
        coded.append({name: dna.upper() for name, dna in parts.items()})
    return tuple(lists), tuple(coded)


def _protein(name: str, dna: str) -> str:
    """Return what one coding sequence spells, checking it rather than writing it again.

    Raises
    ------
    ValueError
        If it is not a whole number of definite codons, or spells a stop anywhere. A part of a
        library construct is read through in frame, so a stop truncates everything after it.
    """
    try:
        spelled = translate(dna)
    except ValueError as error:
        raise ValueError(f"part {name!r} is not a coding sequence: {error}") from error
    if "*" in spelled:
        raise ValueError(
            f"part {name!r} spells a stop at amino acid {spelled.index('*')}, and every part of "
            "a construct is read through in frame"
        )
    return spelled


def _compatible(vector: SequenceRecord, scheme: Scheme) -> bool:
    """Whether the internal enzyme already excises one piece of `vector` on the method's overhangs.

    That piece is DNA that physically exists, so the overhang it spells cannot be re-chosen, and
    position one's entry overhang is pinned to it.
    """
    entry, scar = scheme.entry_overhang, scheme.scar_overhang
    return any(
        (piece.left_overhang, piece.right_overhang) == (entry, scar)
        for piece in digest(vector, scheme.internal)
    )


def _ending(bases: str, overhang: str) -> str:
    """Return `bases` with its last bases replaced by `overhang`."""
    return bases[: len(bases) - len(overhang)] + overhang


def _restandardised(scheme: Scheme, standard: Standard) -> Scheme:
    """Return `scheme` with the standard's first entry overhang written into its stuffers.

    A vector that has to be retrofitted has its internal stuffer synthesised now, so the stuffer
    put in carries the overhang the standard chose for position one rather than the one the
    method's own DNA carries. Only the bases an overhang occupies move; every other base of every
    stuffer stays.

    Raises
    ------
    ValueError
        If the stuffers that overhang makes are not ones the method allows, naming the invariant
        that refused them.
    """
    first = standard.entry_overhangs[0]
    try:
        return dataclasses.replace(
            scheme,
            internal_stuffer_prefix=_ending(scheme.internal_stuffer_prefix, first),
            external_stuffer_5=_ending(scheme.external_stuffer_5, first),
        )
    except ValueError as error:
        raise ValueError(
            "this vector has to be retrofitted, and the internal stuffer carrying the entry "
            f"overhang this design chose is not one method {scheme.name!r} allows: {error}"
        ) from error


def _blocks(parts: Sequence[Part], positions: Sequence[str]) -> dict[str, list[SequenceRecord]]:
    """Every part's synthesis block as a record, keyed by the position it fills."""
    return {
        position: [
            SequenceRecord(part.sequence, name=part.name) for part in parts if part.index == index
        ]
        for index, position in enumerate(positions)
    }


def _barcodes(parts: Sequence[Part], positions: Sequence[str]) -> dict[str, list[str]]:
    """Every part's barcode, keyed by the position it fills: the table that decodes a read."""
    return {
        position: [part.barcode for part in parts if part.index == index]
        for index, position in enumerate(positions)
    }


def _mean(lengths: Sequence[int]) -> int:
    """Return the mean of these lengths, to the nearest whole base."""
    return round(sum(lengths) / len(lengths))


def _bench(
    rounds: Sequence[Round], parts: Sequence[Part], coverage: Sequence[RoundCoverage]
) -> tuple[RoundBench, ...]:
    """Return what each round takes at the bench, weighed at that round's own lengths.

    A round cuts a whole part list in one tube, so the donor is weighed at the pool's mean block
    length; every member carries the same stuffers, so the released fragment is that mean less
    what the external stuffers keep. The destination is one molecule and is weighed exactly.
    """
    made: list[RoundBench] = []
    for one, row in zip(rounds, coverage, strict=True):
        pool = [part for part in parts if part.index == one.part.index]
        block = _mean([part.length for part in pool])
        released = block - (one.part.length - one.released.length)
        destination = one.destination.name or "library"
        donor = f"{one.position} part list"
        made.append(
            RoundBench(
                one.number,
                one.position,
                destination_digest=digest_amount((destination, len(one.destination))),
                donor_digest=digest_amount((donor, block)),
                ligation=_ligation(
                    destination, len(one.destination) - one.excised.length, donor, released
                ),
                transformation=transformation_amount(
                    (f"round {one.number} ligation", len(one.product))
                ),
                coverage=row,
            )
        )
    return tuple(made)


def _ligation(destination: str, opened: int, donor: str, released: int) -> tuple[Amount, Amount]:
    """Return what one round's ligation takes, the opened destination first."""
    return ligation_amounts((f"{destination}, opened", opened), (f"{donor}, released", released))

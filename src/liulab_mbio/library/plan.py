"""One combinatorial library, planned from its part lists, a scheme and a destination vector.

`plan_library` is the one way in. It reads the part lists, chooses the overhang standard the
proteins cost least, builds every part's synthesis sequence, makes the vector a destination,
simulates every round, and works out what each round takes at the bench and how many colonies it
needs. `LibraryPlan.write` puts one directory's worth of output in one place: the synthesis order
sheet, the barcode table, the amino-acid change table, a record for every round, the protocol as
JSON data, and the page rendered from that data.

**The vector and the standard have to agree about position one.** A vector already carrying an
internal stuffer spells an entry overhang in DNA that exists, so position one's overhang is
pinned to it and the standard designs around that. A vector that has to be retrofitted has its
stuffer synthesised now, so the standard chooses position one freely and the stuffer put in
carries what it chose. Nothing else in the scheme moves either way.
"""

import dataclasses
import os
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from liulab_mbio.barcodes import SEED, BarcodeRules
from liulab_mbio.bench.amounts import Amount
from liulab_mbio.checks import Check, Status, worst
from liulab_mbio.codons import codon_usage
from liulab_mbio.goldengate.design import MIN_DISTANCE
from liulab_mbio.io import read_record
from liulab_mbio.library.bench import digest_amount, ligation_amounts, transformation_amount
from liulab_mbio.library.coverage import RoundCoverage, constructs, plan_coverage
from liulab_mbio.library.parts import (
    Part,
    barcode_table,
    change_table,
    design_parts,
    synthesis_sheet,
)
from liulab_mbio.library.rounds import Round, assemble_rounds, representative, write_records
from liulab_mbio.library.scheme import Scheme, read_scheme
from liulab_mbio.library.standard import PartList, Standard, design_standard
from liulab_mbio.library.steps import RoundBench
from liulab_mbio.library.steps import protocol as protocol_for
from liulab_mbio.library.vector import Destination, Site, destination_vector
from liulab_mbio.protocol import Protocol, read_protocol, write_html, write_protocol
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.sites import digest
from liulab_mbio.translate import translate

#: What a part list holds: the proteins each member codes for, or the DNA it is already coded in.
type Kind = Literal["protein", "dna"]

#: How a record's name says which part list it belongs to. ``{position}`` stands for a position's
#: name, and what is around it holds that name to a word boundary, so ``N_ATF2`` names position
#: ``N`` and ``C_NFAT`` does not. A caller whose names read another way passes another pattern.
NAME_PATTERN = r"(?<![A-Za-z0-9]){position}(?![A-Za-z0-9])"

#: What `LibraryPlan.write` calls the files it writes. The records are named by
#: `liulab_mbio.library.rounds`, which writes one for each round and the product for the last.
PARTS_FILE = "parts.tsv"
BARCODE_FILE = "barcodes.tsv"
CHANGE_FILE = "changes.tsv"
PROTOCOL_DATA_FILE = "protocol.json"
PROTOCOL_FILE = "protocol.html"


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
    """

    parts: Path
    barcodes: Path
    changes: Path
    records: tuple[Path, ...]
    protocol_data: Path
    protocol: Path


@dataclass(frozen=True, slots=True)
class LibraryPlan:
    """One planned library: what it costs, what it makes, and what the bench has to do.

    Parameters
    ----------
    scheme
        The architecture the build was given.
    vector
        The record the plan was made from, before any retrofit.
    destination
        The vector a round can open, and what making it one changed.
    part_lists
        One per position, in the scheme's order: each member's name and the protein it codes for,
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
    host
        The codon usage table the coding bases were written for.
    name
        What each round's product is called.
    """

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
    host: str
    name: str = ""

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
    def checks(self) -> tuple[Check, ...]:
        """Every round's verdicts, each named for the round that made it.

        A round judges the product it made: whether the next round can open it, whether the
        enzyme that released the part is gone, and after the last round whether what the product
        keeps past its final part reads in frame without a stop.
        """
        return tuple(
            dataclasses.replace(check, name=f"round {one.number} {check.name}")
            for one in self.rounds
            for check in one.checks
        )

    @property
    def status(self) -> Status:
        """The worst status of any check."""
        return worst(check.status for check in self.checks)

    def protocol(self) -> Protocol:
        """Return the bench protocol for this plan, covering every round as one experiment."""
        return protocol_for(
            scheme=self.scheme,
            vector=self.vector,
            destination=self.destination,
            part_lists=self.part_lists,
            standard=self.standard,
            parts=self.parts,
            rounds=self.rounds,
            bench=self.bench,
            constructs=self.constructs,
            checks=self.checks,
            host=self.host,
            sheet=PARTS_FILE,
            barcodes=BARCODE_FILE,
        )

    def write(self, directory: str | os.PathLike[str]) -> Files:
        """Write the sheets, the records, the protocol data and its page into `directory`.

        The directory is made when it is not there. The files are named by `PARTS_FILE`,
        `BARCODE_FILE`, `CHANGE_FILE`, `liulab_mbio.library.rounds.ROUND_FILE` and
        `PRODUCT_FILE`, `PROTOCOL_DATA_FILE` and `PROTOCOL_FILE`, and a second run over the same
        inputs writes the same bytes. The page is rendered from the data as written, so the two
        cannot disagree.
        """
        out = Path(directory)
        out.mkdir(parents=True, exist_ok=True)
        sheet = out / PARTS_FILE
        sheet.write_text(synthesis_sheet(self.parts), encoding="utf-8")
        barcodes = out / BARCODE_FILE
        barcodes.write_text(barcode_table(self.parts, self.scheme), encoding="utf-8")
        changes = out / CHANGE_FILE
        changes.write_text(change_table(self.standard), encoding="utf-8")
        records = write_records(self.rounds, out)
        data = write_protocol(self.protocol(), out / PROTOCOL_DATA_FILE)
        page = write_html(read_protocol(data), out / PROTOCOL_FILE)
        return Files(sheet, barcodes, changes, records, data, page)


def plan_library(
    parts: str | os.PathLike[str] | Sequence[Mapping[str, str]],
    scheme: Scheme | str | os.PathLike[str],
    vector: SequenceRecord | str | os.PathLike[str],
    *,
    host: str,
    coverage: float,
    kind: Kind = "protein",
    site: Site | None = None,
    pattern: str = NAME_PATTERN,
    rules: BarcodeRules | None = None,
    seed: int = SEED,
    name: str = "",
    min_distance: int = MIN_DISTANCE,
    allow_uniform: bool = False,
) -> LibraryPlan:
    """Plan the whole library `parts` makes under `scheme`, built into `vector`.

    One round appends one part list to every member of the library at once, so the rounds run in
    the scheme's own order and the product of each opens the next. The same inputs return the
    same design: the barcodes are drawn from `seed`, and nothing else here is random.

    Parameters
    ----------
    parts
        A FASTA holding every part list, each record named so that `pattern` says which position
        it fills, or one already-sorted mapping per position in the scheme's order.
    scheme
        The architecture the build is given, or a path to the JSON holding it.
    vector
        The destination, circular: a record, or a path to a ``.dna``, GenBank or FASTA file.
    host
        The name of the codon usage table the coding bases are written for. Never assumed:
        `liulab_mbio.codons.codon_tables` lists the tables that ship.
    coverage
        How many times over each round's colonies have to hold every product it can make. No
        default: how much of a library a round may lose is the caller's call.
    kind
        What the part lists hold. ``"dna"`` is checked and kept rather than written again; only a
        codon a junction or a forbidden site moves is this package's.
    site
        Where to put an internal stuffer, as a feature name or a ``(start, end)`` span. Read only
        where the vector carries none.
    pattern
        How a record's name says which part list it belongs to; see `NAME_PATTERN`.
    rules
        What every barcode holds to. `liulab_mbio.library.parts.barcode_rules` by default.
    seed
        The seed the barcodes are drawn with.
    name
        What to call each round's product. The vector's own name by default.
    min_distance, allow_uniform
        How far apart the standard's overhangs must stand, and whether one base kind is allowed.

    Returns
    -------
    LibraryPlan
        The design, the simulated rounds, and what the bench has to do.

    Raises
    ------
    ValueError
        If a record's name says no position of the scheme or says more than one, if two records
        share a name, if a part is not a protein or not a coding sequence, if no overhang
        standard fits the part lists, if a block spells a site the scheme does not expect, if the
        vector cannot be made a destination, or if a round cannot ligate.
    KeyError
        If no shipped codon usage table is called `host`, or the scheme names an enzyme this
        package does not ship.
    liulab_mbio.barcodes.SpaceExhaustedError
        If a part list is larger than the barcodes its rules allow.

    Examples
    --------
    >>> plan = plan_library("parts.fasta", "scheme.json", "vector.dna",
    ...                     host="human", coverage=10)  # doctest: +SKIP
    >>> plan.write("library/")  # doctest: +SKIP
    """
    design = _scheme(scheme)
    one = _record(vector)
    if isinstance(parts, str | os.PathLike):
        given: Sequence[Mapping[str, str]] = read_part_lists(parts, design, pattern=pattern)
    else:
        given = parts
    lists, coded = _sequences(_checked(given, design), kind)
    compatible = _compatible(one, design)
    pinned = {design.positions[0].name: design.entry_overhang(0)} if compatible else {}
    standard = design_standard(
        design,
        lists,
        pinned=pinned,
        usage=codon_usage(host),
        min_distance=min_distance,
        allow_uniform=allow_uniform,
    )
    built = design_parts(design, lists, standard, host=host, coding=coded, rules=rules, seed=seed)
    destination = destination_vector(
        one, design if compatible else _restandardised(design, standard), site=site
    )
    rounds = assemble_rounds(
        destination.record, representative(built, design), design, name=name or one.name
    )
    rows = plan_coverage([len(each) for each in lists], coverage=coverage)
    return LibraryPlan(
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
        host,
        name or one.name,
    )


def read_part_lists(
    path: str | os.PathLike[str], scheme: Scheme, *, pattern: str = NAME_PATTERN
) -> tuple[dict[str, str], ...]:
    """Read one FASTA and sort its records into one part list a position, in the scheme's order.

    A record is ordered under the name the FASTA gives it, and that name is what says which
    position it fills.

    Raises
    ------
    ValueError
        If the file holds no record, if two records share a name, if a name says no position of
        the scheme or says more than one, or if a position ends up with nothing to fill it.

    Examples
    --------
    >>> read_part_lists("parts.fasta", scheme)  # doctest: +SKIP
    ({'N_ATF2': 'MKT...'}, {'bZIP_JUN': 'WQA...'}, {'C_VP64': 'MKT...'})
    """
    return _sorted(_fasta(path), scheme, pattern)


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


def _sorted(records: Mapping[str, str], scheme: Scheme, pattern: str) -> tuple[dict[str, str], ...]:
    """Put each record in the one part list its name names.

    Raises
    ------
    ValueError
        If a name says no position or says more than one, or a position is left empty.
    """
    names = [position.name for position in scheme.positions]
    matchers = [
        re.compile(pattern.replace("{position}", re.escape(one)), re.IGNORECASE) for one in names
    ]
    lists: list[dict[str, str]] = [{} for _ in names]
    for named, sequence in records.items():
        hit = [index for index, matcher in enumerate(matchers) if matcher.search(named)]
        if not hit:
            raise ValueError(
                f"the name {named!r} says no position of scheme {scheme.name!r}, whose positions "
                f"are {', '.join(names)}. Rename the record, or pass a pattern that reads the "
                "names you already have"
            )
        if len(hit) > 1:
            said = ", ".join(names[index] for index in hit)
            raise ValueError(
                f"the name {named!r} says {len(hit)} positions of scheme {scheme.name!r} — "
                f"{said} — so which part list it belongs to is ambiguous"
            )
        lists[hit[0]][named] = sequence
    if empty := [one for one, parts in zip(names, lists, strict=True) if not parts]:
        raise ValueError(
            f"no record names position(s) {', '.join(empty)}: every position needs a part list"
        )
    return tuple(lists)


def _checked(given: Sequence[Mapping[str, str]], scheme: Scheme) -> tuple[Mapping[str, str], ...]:
    """Refuse part lists that are not one a position, or that name one part twice.

    Raises
    ------
    ValueError
        Naming the count, or the name two part lists share.
    """
    lists = tuple(given)
    if len(lists) != scheme.position_count:
        said = ", ".join(position.name for position in scheme.positions)
        raise ValueError(
            f"this scheme has {scheme.position_count} position(s) — {said} — and "
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
    """Whether the internal enzyme already excises one piece of `vector` on the scheme's overhangs.

    That piece is DNA that physically exists, so the overhang it spells cannot be re-chosen, and
    position one's entry overhang is pinned to it.
    """
    entry, scar = scheme.entry_overhang(0), scheme.scar_overhang
    return any(
        (piece.left_overhang, piece.right_overhang) == (entry, scar)
        for piece in digest(vector, scheme.internal)
    )


def _ending(bases: str, overhang: str) -> str:
    """Return `bases` with its last bases replaced by `overhang`."""
    return bases[: len(bases) - len(overhang)] + overhang


def _restandardised(scheme: Scheme, standard: Standard) -> Scheme:
    """Return `scheme` with the standard's entry overhangs written into its stuffers.

    A vector that has to be retrofitted has its internal stuffer synthesised now, so the stuffer
    put in carries the overhang the standard chose rather than the one the scheme was written
    with. Only the bases an overhang occupies move; every other base of every stuffer stays.

    Raises
    ------
    ValueError
        If the stuffers those overhangs make are not ones this scheme allows, naming the
        invariant that refused them.
    """
    entry = standard.entry_overhangs
    positions = tuple(
        dataclasses.replace(
            position,
            internal_stuffer_prefix=_ending(
                position.internal_stuffer_prefix, entry[(index + 1) % len(entry)]
            ),
            external_stuffer_5=_ending(position.external_stuffer_5, entry[index]),
        )
        for index, position in enumerate(scheme.positions)
    )
    try:
        return dataclasses.replace(scheme, positions=positions)
    except ValueError as error:
        raise ValueError(
            "this vector has to be retrofitted, and the internal stuffer carrying the entry "
            f"overhangs this design chose is not one scheme {scheme.name!r} allows: {error}"
        ) from error


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


def _scheme(value: Scheme | str | os.PathLike[str]) -> Scheme:
    """Read a scheme, or take one already read."""
    return value if isinstance(value, Scheme) else read_scheme(value)


def _record(value: SequenceRecord | str | os.PathLike[str]) -> SequenceRecord:
    """Read a record, or take one already read."""
    return value if isinstance(value, SequenceRecord) else read_record(value)

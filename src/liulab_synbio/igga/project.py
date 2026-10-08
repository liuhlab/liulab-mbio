"""The `Build`: what one build chooses, read from `project.json` and checked as it is read.

A build holds what the method leaves open — which positions, which proteins, which vector, how
deep to sample — and nothing the method's own molecules already carry. A second build is a
second file and no change to this package. `docs/adr/0010-method-in-code.md` draws the line.

The file keeps the name `project.json`, which is what a user writes. A `Project` is something
else here: `liulab_mbio.protocol.model.Project`, the chain of protocols a build writes.

Nothing a build states replaces a method constant. Where the two meet they compose: the
enzymes a block is kept clear of are the method's unioned with `reserved_extra`, and the barcode
length is checked against the method's cloning scar rather than against a number stated here.
"""

import os
from collections.abc import Mapping
from dataclasses import KW_ONLY, dataclass, field
from pathlib import Path
from typing import Any

from liulab_mbio import jsonfile
from liulab_mbio.barcodes import MIN_DISTANCE, SEED
from liulab_mbio.bench.coverage import REPRESENTATION_MARKS, RepresentationMarks
from liulab_mbio.bench.pcr import PRIMER_STOCK_UM
from liulab_mbio.codons import codon_tables
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_synbio.dmx import ROUTE_INDEX_PCR, ROUTES, refuse_unclonal
from liulab_synbio.igga.method import IGGA, Scheme, refuse

#: The method's own barcode length, which a build takes unless it states another.
BARCODE_LENGTH = 11

#: The plate format a primer order comes in unless a build names another, and how many working
#: copies it splits. Neither is a measurement: one is the format every supplier quotes a plated
#: oligo order in, the other is one plate a run.
PRIMER_PLATE_WELLS = 96
PRIMER_PLATE_COPIES = 1


@dataclass(frozen=True, slots=True)
class Barcode:
    """What one part's barcode holds to.

    Parameters
    ----------
    length
        How many bases name one part.
    min_distance
        How far apart two barcodes of one part list must stand, by the metric
        `liulab_mbio.barcodes` measures with.
    """

    length: int = BARCODE_LENGTH
    min_distance: int = MIN_DISTANCE


@dataclass(frozen=True, slots=True)
class FinalAssembly:
    """What a lab that has run the one-pot assembly states about it.

    Both numbers travel together: a mass with no ratio sizes one side of the pot, and a ratio
    with no mass sizes neither. Nothing published sizes this reaction, so a build that has
    measured it states both and a build that has not leaves the hole standing.

    Parameters
    ----------
    vector_ng
        How much working vector goes into the one-pot assembly, ng.
    ratio
        How much cargo meets it, as a molar ratio of cargo to vector.
    """

    vector_ng: float
    ratio: float

    def __post_init__(self) -> None:
        """Refuse a mass or a ratio that is not positive."""
        for key, value in (("vector_ng", self.vector_ng), ("ratio", self.ratio)):
            if value <= 0:
                raise ValueError(
                    f"final_assembly.{key} is {value}, and a build states a positive one"
                )


@dataclass(frozen=True, slots=True)
class PrimerPlates:
    """How a run lays its routine primers out, which is the lab's own and not the method's.

    The amounts are the build author's: what the vendor delivers in a well and what this lab
    resuspends and dilutes to. None of them is published anywhere this package can cite, so a
    build that wants the plates states them and one that does not gets no such protocol.

    Parameters
    ----------
    nanomoles
        What the vendor delivers in one well, which sets the volume it is resuspended in.
    stock_um
        What the stock plate is resuspended to, µM.
    working_ul
        What one well of a working plate holds.
    working_um
        What a working plate is diluted to, µM. It defaults to what a PCR is pipetted from,
        `liulab_mbio.bench.pcr.PRIMER_STOCK_UM`.
    wells
        The format every plate comes in.
    copies
        How many working plates to split, one run each.
    """

    nanomoles: float
    stock_um: float
    working_ul: float
    working_um: float = PRIMER_STOCK_UM
    wells: int = PRIMER_PLATE_WELLS
    copies: int = PRIMER_PLATE_COPIES

    def __post_init__(self) -> None:
        """Refuse an amount that is not positive, or a working plate that dilutes nothing."""
        for key, value in (
            ("nanomoles", self.nanomoles),
            ("stock_um", self.stock_um),
            ("working_ul", self.working_ul),
            ("working_um", self.working_um),
            ("wells", self.wells),
            ("copies", self.copies),
        ):
            if value <= 0:
                raise ValueError(
                    f"primer_plates.{key} is {value}, and a build states a positive one"
                )
        if self.working_um >= self.stock_um:
            raise ValueError(
                f"primer_plates.working_um is {self.working_um} and stock_um is "
                f"{self.stock_um}, so a working plate dilutes nothing"
            )


@dataclass(frozen=True, slots=True)
class Build:
    """What one build chooses, checked on construction.

    Parameters
    ----------
    name
        What the run and each round's product is called.
    positions
        The positions, in the order the rounds fill them, so the last is the terminal one. A
        position is a name: every one of them carries the method's own stuffers.
    parts
        The FASTA of every part list, each record named for the position it fills.
    vector
        The destination vector, circular: a ``.dna``, GenBank or FASTA file.
    working_vector
        The vector the finished library is moved into, circular, before its ccdB cassette: a
        stock the lab holds and an application chooses. Omitted, the library stays in the
        destination vector and the final assembly is written as what it cannot say.
    host
        The codon usage table the coding bases are written for.
    oligo_length
        How long one synthesised oligo of the pool may be.
    primers
        The orthogonal primer set the pool is amplified by, as a two-column sheet the user
        holds. Without one no oligo can be written, because the sites are templated on it.
    bands
        The vendor's bands for each quantity the pool is banded by, written ``low-high`` with an
        empty top for a tier with no top. They are what headroom is measured against, and they
        are reported whether or not anyone holds a price.
    batch_size
        How many parts one batch of the bench work carries.
    completeness
        The chance each round is sized for that none of its products is missing.
    validate_from
        The fragment-count floor at or above which a design is read back one well at a time.
        Omitted, nothing is read and the library stays polyclonal; ``0`` reads every design.
        There is no default: `liulab_mbio.bench.readback.clean_colony_chance` gives a design's
        chance of a clean colony, not the chance worth paying to check.
    routes
        Which of `liulab_synbio.dmx.ROUTES` read those wells back, one or more. Named exactly
        when `validate_from` is, because an unread build needs no route. Name two and the run
        writes a page for each, as the two ways of one job the bench does one of. A build that
        reads anything back names `primers` too: without a pool each block arrives as the vendor
        ships it, and DMX has no colony to pick.
    index_plate
        What this lab calls its prepared plate of barcoded primer pairs. Only the index PCR
        route takes one, so only that route may name it.
    seed
        The seed the barcodes are drawn with.
    reserved_extra
        Further enzymes this build needs a block kept clear of, added to the method's own.
    representation_seen, representation_skew, reads_per_member
        What this build holds the representation read to: the share of combinations that must be
        read at all, the skew ratio the counts must stay under, and the depth both are judged at.
        Each takes `liulab_mbio.bench.coverage.REPRESENTATION_MARKS` where the build states
        none, and each may only be tightened.
    linkage_fidelity
        The share of reads whose barcode must still name its part. No default: nothing published
        sets a mark for it, so a build that states none is read against no mark at all.
    final_assembly
        What this lab measured the one-pot assembly at: the working vector's mass and the molar
        ratio the cargo meets it at. No default, as for the two cycle counts below: a build that
        has run the pilot states it and a build that has not leaves the hole standing.
    pcr1_cycles, pcr2_cycles
        The cycle counts this lab measured for PCR1, against the polymerase it runs, and for
        PCR2. No default: nobody published either, so a build stating neither prints the band
        Twist gives for another polymerase and no PCR2 count at all.
    barcode
        What one part's barcode holds to.
    primer_plates
        How this lab lays the pool's primers out as a stock plate and its working copies.
        Omitted, no such protocol is written: an empty page is worse than a step, and the
        amounts are nobody's to guess. It needs `primers`, which is what there is to plate.
    scheme
        The method the build is made by. There is one, and it is `IGGA`.

    Raises
    ------
    ValueError
        If a position is repeated or missing, a number is not positive, the completeness does not
        lie between 0 and 1, the floor is negative, a route is neither of the two or named twice,
        the floor and the routes are not both there or both absent, an index plate is named
        without the route that takes one, a route is named over cargo DMX cannot pick, a
        representation mark loosens the sourced one, or the barcode and the method's cloning
        scar are not whole codons together — which names ``barcode-frame``.
    KeyError
        If `reserved_extra` names an enzyme this package does not ship, or `host` no shipped
        codon usage table.
    """

    name: str
    _: KW_ONLY
    positions: tuple[str, ...]
    parts: Path
    vector: Path
    host: str
    oligo_length: int
    batch_size: int
    completeness: float
    primers: Path | None = None
    working_vector: Path | None = None
    bands: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    validate_from: int | None = None
    routes: tuple[str, ...] = ()
    index_plate: str = ""
    seed: int = SEED
    reserved_extra: tuple[str, ...] = ()
    representation_seen: float | None = None
    representation_skew: float | None = None
    reads_per_member: int | None = None
    linkage_fidelity: float | None = None
    final_assembly: FinalAssembly | None = None
    pcr1_cycles: int | None = None
    pcr2_cycles: int | None = None
    barcode: Barcode = field(default_factory=Barcode)
    primer_plates: PrimerPlates | None = None
    scheme: Scheme = IGGA

    def __post_init__(self) -> None:
        """Normalise the lists and paths, then check everything a designer will read."""
        object.__setattr__(self, "positions", tuple(self.positions))
        object.__setattr__(self, "reserved_extra", tuple(self.reserved_extra))
        object.__setattr__(self, "parts", Path(self.parts))
        object.__setattr__(self, "vector", Path(self.vector))
        if self.working_vector is not None:
            object.__setattr__(self, "working_vector", Path(self.working_vector))
        object.__setattr__(self, "bands", dict(self.bands))
        self._check_positions()
        self._check_numbers()
        self._check_validation()
        self._check_plates()
        self._check_marks()
        self._check_barcode_frame()
        if self.host not in codon_tables():
            raise KeyError(
                f"no shipped codon usage table is called {self.host!r}: {', '.join(codon_tables())}"
            )
        for named in self.reserved_extra:
            get_enzyme(named)

    @property
    def position_count(self) -> int:
        """How many positions one product joins."""
        return len(self.positions)

    @property
    def reserved(self) -> tuple[str, ...]:
        """Every enzyme a block is kept clear of: the method's, then this build's own.

        A build adds and never replaces, so an enzyme the method reserves stays reserved
        whatever a build says.
        """
        return tuple(dict.fromkeys((*self.scheme.reserved, *self.reserved_extra)))

    @property
    def reserved_enzymes(self) -> tuple[Enzyme, ...]:
        """Those enzymes, in that order."""
        return tuple(get_enzyme(name) for name in self.reserved)

    @property
    def barcode_block_length(self) -> int:
        """How long the finished barcode block is: one barcode a position, joined by the scar."""
        return self.scheme.barcode_block_length(self.barcode.length, self.position_count)

    @property
    def retained_length(self) -> int:
        """What the finished product keeps past its last part: stuffer and barcode block."""
        return self.scheme.retained_length(self.barcode.length, self.position_count)

    def _check_positions(self) -> None:
        """Refuse a build with no position, an unnamed one, or one named twice."""
        if not self.positions:
            raise ValueError("a build needs at least one position")
        if any(not one for one in self.positions):
            raise ValueError("every position of a build is named")
        if len(set(self.positions)) != len(self.positions):
            raise ValueError(f"the positions {', '.join(self.positions)} name one of them twice")

    def _check_numbers(self) -> None:
        """Refuse a build whose lengths or counts are not positive, or whose chance is not one."""
        for named, value in (
            ("oligo_length", self.oligo_length),
            ("batch_size", self.batch_size),
            ("barcode.length", self.barcode.length),
            ("barcode.min_distance", self.barcode.min_distance),
        ):
            if value <= 0:
                raise ValueError(f"{named} is {value}, and a build states a positive one")
        if not 0.0 < self.completeness < 1.0:
            raise ValueError(
                f"completeness is {self.completeness}, and a build states a chance between 0 and 1"
            )
        for named, cycles in (("pcr1_cycles", self.pcr1_cycles), ("pcr2_cycles", self.pcr2_cycles)):
            if cycles is not None and cycles <= 0:
                raise ValueError(f"{named} is {cycles}, and a build states a positive one")

    def _check_validation(self) -> None:
        """Refuse a bad floor or route, one of the two alone, a stray plate, or unclonal cargo.

        The floor and the routes travel together: a floor with no route says which designs are
        read and not how, and a route with no floor names a read nobody asked for. Whether
        there is anything to read back at all is DMX's own question, not a build's.

        Naming a route twice is refused rather than collapsed: it would write one page twice and
        ask the bench to choose between a page and itself.
        """
        if self.validate_from is not None and self.validate_from < 0:
            raise ValueError(
                f"validate_from is {self.validate_from}, and a fragment-count floor counts "
                "fragments; omit it to read nothing, or set 0 to read every design"
            )
        for one in self.routes:
            if one not in ROUTES:
                raise ValueError(
                    f"route is {one!r}, and a build reads its wells back on one of "
                    f"{', '.join(repr(each) for each in ROUTES)}"
                )
        if len(set(self.routes)) != len(self.routes):
            raise ValueError(
                "routes names the same route twice, and the ways of one job are different ways"
            )
        if (self.validate_from is None) != (not self.routes):
            raise ValueError(
                "validate_from and routes are stated together: a build that reads designs back "
                f"says which, and on which of {', '.join(repr(one) for one in ROUTES)}"
            )
        if self.index_plate and ROUTE_INDEX_PCR.name not in self.routes:
            raise ValueError(
                "index_plate names a plate of barcoded primer pairs, which only the "
                f"{ROUTE_INDEX_PCR.name!r} route takes"
            )
        if self.routes:
            refuse_unclonal(
                "a block as the vendor ships it, which is what a build naming no primer set orders",
                clonal=self.primers is not None,
            )

    def _check_plates(self) -> None:
        """Refuse primer plates where there is no primer set to plate.

        The plates seat the primers that amplify the pool, and a build naming no primer set
        writes no pool and no primer.
        """
        if self.primer_plates is not None and self.primers is None:
            raise ValueError(
                "primer_plates lays out the primers that amplify the pool, so a build naming "
                "it names primers too"
            )

    def _check_marks(self) -> None:
        """Refuse a representation mark that loosens the sourced one, or a share outside 0 to 1.

        A build tightens and never loosens, the same rule `liulab_synbio.dmx.depth_check`
        holds a read depth to: the sourced mark is the floor the method stands on.
        """
        for named, value in (
            ("representation_seen", self.representation_seen),
            ("linkage_fidelity", self.linkage_fidelity),
        ):
            if value is not None and not 0.0 < value <= 1.0:
                raise ValueError(f"{named} is {value}, and a build states a share of 1 or less")
        sourced = REPRESENTATION_MARKS
        if self.representation_seen is not None and self.representation_seen < sourced.seen:
            raise ValueError(
                f"representation_seen is {self.representation_seen}, and a build may only "
                f"raise the {sourced.seen} share Joung sets, not lower it"
            )
        if self.representation_skew is not None and self.representation_skew > sourced.skew:
            raise ValueError(
                f"representation_skew is {self.representation_skew}, and a build may only "
                f"lower the skew ratio of {sourced.skew} Joung allows, not raise it"
            )
        if self.representation_skew is not None and self.representation_skew < 1.0:
            raise ValueError(
                f"representation_skew is {self.representation_skew}, and the evenest library "
                "a count can describe has a ratio of 1"
            )
        if self.reads_per_member is not None and self.reads_per_member < sourced.reads_per_member:
            raise ValueError(
                f"reads_per_member is {self.reads_per_member}, and a build may only raise the "
                f"{sourced.reads_per_member} reads a member Joung judges at, not lower it"
            )

    @property
    def marks(self) -> RepresentationMarks:
        """What this build holds its representation reads to: the sourced marks, as tightened."""
        return RepresentationMarks(
            seen=(
                REPRESENTATION_MARKS.seen
                if self.representation_seen is None
                else self.representation_seen
            ),
            skew=(
                REPRESENTATION_MARKS.skew
                if self.representation_skew is None
                else self.representation_skew
            ),
            reads_per_member=(
                REPRESENTATION_MARKS.reads_per_member
                if self.reads_per_member is None
                else self.reads_per_member
            ),
        )

    def _check_barcode_frame(self) -> None:
        """Check a barcode and the scar joining it to the last make whole codons together."""
        scar = len(self.scheme.cloning_scar)
        unit = self.barcode.length + scar
        if unit % 3:
            refuse(
                "barcode-frame",
                f"a barcode of {self.barcode.length} bases and a cloning scar of {scar} make "
                f"{unit}, which is not a whole number of codons",
            )


def read_build(path: str | os.PathLike[str]) -> Build:
    """Read a build from JSON, resolving every file path against the file's own directory.

    Raises
    ------
    ValueError
        On a missing or unknown key, a value of another JSON type, a path naming no file, or any
        check `Build` makes.
    KeyError
        If it names an enzyme or a codon usage table this package does not ship.

    Examples
    --------
    >>> read_build("project.json").positions  # doctest: +SKIP
    ('N', 'DBD', 'C')
    """
    file = Path(path)
    data = jsonfile.read_object(path)
    jsonfile.refuse_keys(data, _BUILD_KEYS, _BUILD_OPTIONAL, "a build")
    given = dict(data)
    return Build(
        jsonfile.text(given, "name", "a build's"),
        positions=tuple(
            jsonfile.one_text(one, f"positions[{index}]")
            for index, one in enumerate(jsonfile.listing(given, "positions", "a build's"))
        ),
        parts=jsonfile.named_file(
            file, jsonfile.text(given, "parts", "a build's"), "parts", "a build's"
        ),
        vector=jsonfile.named_file(
            file, jsonfile.text(given, "vector", "a build's"), "vector", "a build's"
        ),
        host=jsonfile.text(given, "host", "a build's"),
        oligo_length=jsonfile.whole(given, "oligo_length", "a build's"),
        batch_size=jsonfile.whole(given, "batch_size", "a build's"),
        completeness=jsonfile.number(given, "completeness", "a build's"),
        primers=(
            jsonfile.named_file(
                file, jsonfile.text(given, "primers", "a build's"), "primers", "a build's"
            )
            if "primers" in given
            else None
        ),
        working_vector=(
            jsonfile.named_file(
                file,
                jsonfile.text(given, "working_vector", "a build's"),
                "working_vector",
                "a build's",
            )
            if "working_vector" in given
            else None
        ),
        bands=_bands(given.get("bands")),
        validate_from=(
            jsonfile.whole(given, "validate_from", "a build's")
            if "validate_from" in given
            else None
        ),
        routes=tuple(
            jsonfile.one_text(one, f"routes[{index}]")
            for index, one in enumerate(
                jsonfile.listing(given, "routes", "a build's") if "routes" in given else ()
            )
        ),
        index_plate=jsonfile.text(given, "index_plate", "a build's")
        if "index_plate" in given
        else "",
        representation_seen=(
            jsonfile.number(given, "representation_seen", "a build's")
            if "representation_seen" in given
            else None
        ),
        representation_skew=(
            jsonfile.number(given, "representation_skew", "a build's")
            if "representation_skew" in given
            else None
        ),
        reads_per_member=(
            jsonfile.whole(given, "reads_per_member", "a build's")
            if "reads_per_member" in given
            else None
        ),
        linkage_fidelity=(
            jsonfile.number(given, "linkage_fidelity", "a build's")
            if "linkage_fidelity" in given
            else None
        ),
        final_assembly=_final_assembly(given.get("final_assembly")),
        pcr1_cycles=jsonfile.whole(given, "pcr1_cycles", "a build's")
        if "pcr1_cycles" in given
        else None,
        pcr2_cycles=jsonfile.whole(given, "pcr2_cycles", "a build's")
        if "pcr2_cycles" in given
        else None,
        seed=jsonfile.whole(given, "seed", "a build's") if "seed" in given else SEED,
        reserved_extra=tuple(
            jsonfile.one_text(one, f"reserved_extra[{index}]")
            for index, one in enumerate(
                jsonfile.listing(given, "reserved_extra", "a build's")
                if "reserved_extra" in given
                else ()
            )
        ),
        barcode=_barcode(given.get("barcode")),
        primer_plates=_primer_plates(given.get("primer_plates")),
    )


#: The keys a build is written with, and the ones it may leave out.
_BUILD_KEYS = frozenset(
    {
        "name",
        "positions",
        "parts",
        "vector",
        "host",
        "oligo_length",
        "batch_size",
        "completeness",
    }
)
_BUILD_OPTIONAL = frozenset(
    {
        "seed",
        "reserved_extra",
        "barcode",
        "primers",
        "working_vector",
        "bands",
        "validate_from",
        "routes",
        "index_plate",
        "representation_seen",
        "representation_skew",
        "reads_per_member",
        "linkage_fidelity",
        "final_assembly",
        "pcr1_cycles",
        "pcr2_cycles",
        "primer_plates",
    }
)
_BARCODE_OPTIONAL = frozenset({"length", "min_distance"})
_PLATES_REQUIRED = frozenset({"nanomoles", "stock_um", "working_ul"})
_PLATES_OPTIONAL = frozenset({"working_um", "wells", "copies"})
_ASSEMBLY_REQUIRED = frozenset({"vector_ng", "ratio"})


def _bands(entry: Any) -> Mapping[str, tuple[str, ...]]:
    """Read the vendor's bands a build states, as a quantity naming its tiers.

    Raises
    ------
    ValueError
        If it is not an object of lists of strings.
    """
    if entry is None:
        return {}
    if not isinstance(entry, Mapping):
        raise ValueError(f"a build's bands are {type(entry).__name__}, not an object")
    return {
        quantity: tuple(
            jsonfile.one_text(one, f"bands {quantity}[{index}]")
            for index, one in enumerate(jsonfile.listing(entry, quantity, "a build's bands"))
        )
        for quantity in entry
    }


def _barcode(entry: Any) -> Barcode:
    """Build the barcode rules from parsed JSON, the method's own where it states none.

    Raises
    ------
    ValueError
        If it is not an object of the two keys a barcode is written with.
    """
    if entry is None:
        return Barcode()
    if not isinstance(entry, Mapping):
        raise ValueError(f"a build's barcode is {type(entry).__name__}, not an object")
    jsonfile.refuse_keys(entry, frozenset(), _BARCODE_OPTIONAL, "a build's barcode")
    return Barcode(
        jsonfile.whole(entry, "length", "a build's barcode")
        if "length" in entry
        else BARCODE_LENGTH,
        jsonfile.whole(entry, "min_distance", "a build's barcode")
        if "min_distance" in entry
        else MIN_DISTANCE,
    )


def _final_assembly(entry: Any) -> FinalAssembly | None:
    """Build the one-pot assembly's amounts from parsed JSON, or `None` where a build states none.

    Raises
    ------
    ValueError
        If it is not an object, or a key is missing, unknown or of another JSON type.
    """
    if entry is None:
        return None
    if not isinstance(entry, Mapping):
        raise ValueError(f"a build's final_assembly is {type(entry).__name__}, not an object")
    where = "a build's final_assembly"
    jsonfile.refuse_keys(entry, _ASSEMBLY_REQUIRED, frozenset(), where)
    return FinalAssembly(
        jsonfile.number(entry, "vector_ng", where), jsonfile.number(entry, "ratio", where)
    )


def _primer_plates(entry: Any) -> PrimerPlates | None:
    """Build the primer-plate amounts from parsed JSON, or `None` where a build states none.

    Three amounts are required because nothing publishes them: what the vendor delivers, what
    this lab resuspends to, and what a working well holds. The other three have a default.

    Raises
    ------
    ValueError
        If it is not an object, or a key is missing, unknown or of another JSON type.
    """
    if entry is None:
        return None
    if not isinstance(entry, Mapping):
        raise ValueError(f"a build's primer_plates is {type(entry).__name__}, not an object")
    jsonfile.refuse_keys(entry, _PLATES_REQUIRED, _PLATES_OPTIONAL, "a build's primer_plates")
    where = "a build's primer_plates"
    return PrimerPlates(
        jsonfile.number(entry, "nanomoles", where),
        jsonfile.number(entry, "stock_um", where),
        jsonfile.number(entry, "working_ul", where),
        **(
            {"working_um": jsonfile.number(entry, "working_um", where)}
            if "working_um" in entry
            else {}
        ),
        **({"wells": jsonfile.whole(entry, "wells", where)} if "wells" in entry else {}),
        **({"copies": jsonfile.whole(entry, "copies", where)} if "copies" in entry else {}),
    )

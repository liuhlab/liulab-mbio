"""One library build, as every protocol of its chain reads it.

`Run` carries what a build settled — the scheme, the design, the files it wrote and the lab's
own choices — so a protocol module takes one argument instead of the twenty-five a chain of
dispatches used to thread. What more than one protocol shares also sits here: the items one
hands the next by name, the reagents the rounds buy, and the hardware they need.
"""

from collections.abc import Sequence
from dataclasses import KW_ONLY, dataclass

from liulab_mbio import checks as judged
from liulab_mbio.bench.amounts import REFERENCES as AMOUNT_REFERENCES
from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.coverage import REFERENCES as COVERAGE_REFERENCES
from liulab_mbio.bench.coverage import (
    REPRESENTATION_MARKS,
    RepresentationMarks,
    reads_for_representation,
)
from liulab_mbio.bench.goldengate import REFERENCES as GOLDEN_GATE_REFERENCES
from liulab_mbio.bench.materials import material
from liulab_mbio.bench.phenotype import selection_marker
from liulab_mbio.bench.prices import PriceRecord
from liulab_mbio.bench.steps import enzyme_material, listed
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.protocol.model import Item as Handed
from liulab_mbio.protocol.model import Material, Reference
from liulab_mbio.sequence import SequenceRecord
from liulab_synbio import dmx
from liulab_synbio.igga import stages
from liulab_synbio.igga.bench import (
    CUTSMART,
    DIGEST_CELSIUS,
    ENZYME_UL,
    GROWTH_CELSIUS,
    SPRI_AFTER_DIGEST,
    SPRI_AFTER_LIGATION,
    SPRI_BEADS,
    STRAIN,
    STRAIN_CATALOG,
    STRAIN_SUPPLIER,
    RoundBench,
    choppers,
    pulse,
)
from liulab_synbio.igga.bench import REFERENCES as BENCH_REFERENCES
from liulab_synbio.igga.cargo import PoolPlan
from liulab_synbio.igga.figures import OLIGO_FILE
from liulab_synbio.igga.method import Scheme
from liulab_synbio.igga.parts import Part
from liulab_synbio.igga.project import FinalAssembly, PrimerPlates
from liulab_synbio.igga.reads import ReadPair, ReadPairs
from liulab_synbio.igga.rounds import Round
from liulab_synbio.igga.standard import PartList, Standard
from liulab_synbio.igga.vector import Destination, Working

#: Where the records a figure draws sit, relative to the pages that draw them:
#: `liulab_synbio.igga.plan.PROTOCOL_DIR` puts the pages one directory below the records.
RECORDS_AT = "../"

#: The hardware a round needs, which no reagent table covers.
ROUND_EQUIPMENT: tuple[str, ...] = (
    f"Incubator or heat block at {DIGEST_CELSIUS:g} °C",
    "Magnetic rack for the bead clean-ups",
    "Electroporator and cuvettes",
    f"Shaking incubator at {GROWTH_CELSIUS:g} °C",
    "Spectrophotometer or fluorometer",
    "Long-read sequencer, for the two reads of the finished library",
)

#: What one protocol hands the next, by name. The name is the contract, so the protocol that
#: produces it and the one that consumes it spell it alike.
POOL_ITEM = "oligo pool"
BLOCKS_ITEM = "synthesised blocks"
ARCHIVE_ITEM = "cargo archive plate"
PICKED_ITEM = "clonal picked plate"
CALLS_ITEM = "well calls"
PREP_ITEM = "library prep, round {number}"
LIBRARY_ITEM = "the library in its working vector"
BLOCK_VECTOR_ITEM = "block vector {number}"
WORKING_ITEM = "working vector"

#: What a synthesised block is resuspended in, and the floor both vendors publish for it. Each
#: gives 10 ng/µL as a least, not a value: at it a 1,000 ng pool fills 100 µL, more than twice
#: what the round's digest leaves, so the pool's own floor is computed and these set the buffer.
POOL_RESUSPENSION_REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Twist Bioscience, How should Multiplexed Gene Fragments be resuspended? — nuclease-free "
        "TE pH 8.0 or 10 mM Tris-HCl pH 8.0, at least 10 ng/µL for the stock dilution",
        url="https://www.twistbioscience.com/faq/multiplexed-gene-fragments/how-should-multiplexed-gene-fragments-be-resuspended",
    ),
    Reference(
        "Integrated DNA Technologies, gBlocks Gene Fragments resuspension — spin down, add IDTE "
        "or molecular-grade water to 10 ng/µL, vortex, 50 °C for 15-20 min, then verify",
        url="https://www.idtdna.com/page/?p=890",
    ),
)

#: What the two readout cautions of the confirming step are measured by.
READOUT_REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Van Nieuwerburgh, F. et al. (2011) Quantitative bias in Illumina TruSeq and a novel "
        "post amplification barcoding strategy for multiplexed DNA and small RNA deep "
        "sequencing. PLoS ONE 6, e26969, for a barcode 3 bp from the insert giving up to "
        "100-fold differences in read counts, where one introduced during the PCR 34 bp away "
        "gave R2 = 0.9977",
        url="https://doi.org/10.1371/journal.pone.0026969",
    ),
    Reference(
        "Alon, S. et al. (2011) Barcoding bias in high-throughput multiplex sequencing of "
        "miRNA. Genome Res. 21, 1506-1511, for ligated barcodes spreading read counts about "
        "twofold where the same barcodes introduced during the PCR did not",
        url="https://doi.org/10.1101/gr.121715.111",
    ),
)


@dataclass(frozen=True)
class Run:
    """What one library build settled, which every protocol of its chain is written against.

    Each field is the `liulab_synbio.igga.plan.LibraryPlan` field or property of that name;
    `sheet`, `barcodes`, `changes`, `pool_sheet`, `primer_sheet` and `read_sheet` are what the
    plan calls the sheets its pages point at, and `files` gathers them, with the block vectors'
    own records, as the pages link them.

    `validation` is `None` for a build that states no fragment-count floor, and no read-back
    protocol is written: the library stays polyclonal, which is the default.

    `working` is the vector the finished library is moved into, and `None` for a build naming
    none. Without one the final protocol's steps still run, because it is a stage of the method,
    and what the vector would have fixed is a hole instead.

    `pool` is the oligo pool the blocks are built from, and `None` for a build that writes
    none. With one, the blocks are not ordered: they are amplified out of the pool and
    assembled, so ordering and making the cargo are two protocols instead of one step.

    `block_vectors` is what a block closes into, one a position: what each is called and the
    file the plan wrote it to. `block_records` is the same vectors themselves, in that order,
    which a step's figure is drawn over. `records_at` says where those records sit relative to
    the pages, which is what every figure's path is written against.

    `primer_plates` is how this lab lays the pool's primers out, and `None` for a build
    stating none.

    `final_assembly`, `pcr1_cycles` and `pcr2_cycles` are what this build measured where the
    method leaves the number open, and `None` where it measured none. Stated, the step prints
    the number and says it is this run's own; left out, the hole stands as it does today.
    """

    _: KW_ONLY
    scheme: Scheme
    positions: Sequence[str]
    barcode_length: int
    vector: SequenceRecord
    destination: Destination
    part_lists: Sequence[PartList]
    standard: Standard
    parts: Sequence[Part]
    rounds: Sequence[Round]
    bench: Sequence[RoundBench]
    constructs: int
    checks: Sequence[judged.Check]
    host: str
    sheet: str
    barcodes: str
    changes: str = ""
    read_sheet: str = ""
    validation: dmx.Validation | None = None
    prices: PriceRecord | None = None
    pool: PoolPlan | None = None
    pool_sheet: str = ""
    primer_sheet: str = ""
    working: Working | None = None
    reads: ReadPairs | None = None
    marks: RepresentationMarks = REPRESENTATION_MARKS
    linkage_fidelity: float | None = None
    final_assembly: FinalAssembly | None = None
    pcr1_cycles: int | None = None
    pcr2_cycles: int | None = None
    block_vectors: Sequence[tuple[str, str]] = ()
    block_records: Sequence[Destination] = ()
    primer_plates: PrimerPlates | None = None
    records_at: str = RECORDS_AT

    @property
    def files(self) -> tuple[str, ...]:
        """Every sheet this run writes, as a path from a page, for the pages that name one.

        The run states them once and each page links whichever of them its own text says, so
        no page carries a filename its reader cannot open. A run without a pool writes neither
        pool sheet, so neither is here to be linked.
        """
        pooled = (self.pool_sheet, self.primer_sheet) if self.pool is not None else ()
        return tuple(
            self.records_at + name
            for name in (
                self.sheet,
                self.barcodes,
                self.changes,
                *pooled,
                self.read_sheet,
                *(file for _, file in self.block_vectors),
            )
            if name
        )

    @property
    def inside(self) -> tuple[Enzyme, ...]:
        """The blunt enzymes that cut the stuffer the internal digest excises."""
        return choppers(self.scheme)[0]

    @property
    def outside(self) -> tuple[Enzyme, ...]:
        """The blunt enzymes that cut the external stuffers the donor digest leaves behind."""
        return choppers(self.scheme)[1]

    @property
    def selection(self) -> str:
        """What a round's transformants are plated on, read off the destination's own marker."""
        return stages.selection_for(self.vector)

    @property
    def plated(self) -> bool:
        """Whether this lab lays the pool's primers out as plates of its own."""
        return self.pool is not None and self.primer_plates is not None

    @property
    def oligo(self) -> str:
        """Where the oligo record a pool figure is drawn over sits, from a page."""
        return self.records_at + OLIGO_FILE

    @property
    def destinations(self) -> Sequence[tuple[str, str]]:
        """What a block closes into, one a position, falling back to the round's destination."""
        return self.block_vectors or ((self.rounds[0].destination.name or "the destination", ""),)

    @property
    def block_vector(self) -> tuple[Destination, str, str] | None:
        """The first block vector, what it is called and where its record sits, or nothing.

        A build whose blocks are ordered whole writes no block vector, and nothing is drawn.
        """
        if not self.block_records or not self.block_vectors or not self.block_vectors[0][1]:
            return None
        name, file = self.block_vectors[0]
        return self.block_records[0], name, self.records_at + file

    @property
    def opened(self) -> Amount:
        """What the first round's ligation opens, which a cargo assembly also goes into."""
        return self.bench[0].ligation[0]

    @property
    def ordered(self) -> Handed:
        """What the ordering protocol leaves behind: a pool, or the blocks themselves."""
        if self.pool:
            return Handed(POOL_ITEM, "the pool resuspended", storage="-20 °C")
        return Handed(BLOCKS_ITEM, "every block as the vendor shipped it", storage="-20 °C")

    @property
    def blocks(self) -> tuple[Handed, ...]:
        """Each position's block vector, ready for the digest that opens it."""
        return tuple(
            Handed(BLOCK_VECTOR_ITEM.format(number=number), f"{name}, ready to open")
            for number, (name, _) in enumerate(self.block_vectors, 1)
        ) or (
            Handed(BLOCK_VECTOR_ITEM.format(number=1), f"{self.vector.name or 'the destination'}"),
        )

    @property
    def archive(self) -> Handed:
        """One well a design, which is what making the cargo leaves."""
        return Handed(
            ARCHIVE_ITEM,
            "one well a design, polyclonal, sealed and frozen",
            storage="-80 °C glycerol stock",
        )

    @property
    def picked(self) -> Handed:
        """One well a picked colony, which the read-back takes the archive to."""
        return Handed(PICKED_ITEM, "one well a picked colony, each well one clone")

    @property
    def calls(self) -> Handed:
        """A pass or a fail a well, which is the read-back's other product."""
        return Handed(CALLS_ITEM, "a pass or a fail a well, and which design each well holds")

    @property
    def prep(self) -> Handed:
        """The pooled library after the last round, as a plasmid prep."""
        return Handed(
            PREP_ITEM.format(number=len(self.rounds)),
            f"the pooled library after round {len(self.rounds)}, as a plasmid prep",
            storage="-20 °C",
        )

    @property
    def library(self) -> Handed:
        """The finished library in its working vector, which is what the run is for."""
        return Handed(
            LIBRARY_ITEM,
            f"{len(self.part_lists)} part lists joined, pooled and ready for the screen",
            storage="-20 °C",
        )

    @property
    def marking_stock(self) -> Handed | None:
        """The lab stock the chosen marking route takes, which no protocol of this run makes."""
        return dmx.marking_stock(self.validation)

    @property
    def round_references(self) -> tuple[Reference, ...]:
        """Where a round's numbers come from, and where the scheme itself came from."""
        items = [
            *BENCH_REFERENCES,
            *COVERAGE_REFERENCES,
            *AMOUNT_REFERENCES,
            *READOUT_REFERENCES,
            *POOL_RESUSPENSION_REFERENCES,
            *GOLDEN_GATE_REFERENCES,
        ]
        if self.scheme.source:
            items.append(Reference(f"The method this build was planned by: {self.scheme.source}"))
        return tuple(items)

    @property
    def round_materials(self) -> tuple[Material, ...]:
        """Every reagent and consumable the rounds ask for.

        Where there is a pool the blocks are not bought, so the part lists say what they are
        built from; what the pool itself buys is the ordering protocol's.
        """
        came_from = (
            f"assembled from the oligo pool; {self.pool_sheet} says which oligos"
            if self.pool
            else f"synthesised blocks, ordered from {self.sheet}"
        )
        made = [
            Material(
                f"{position} part list",
                storage="-20 °C",
                amount=f"{len(one)} members, pooled",
                note=came_from,
            )
            for position, one in zip(self.positions, self.part_lists, strict=True)
        ]
        made.append(
            Material(
                f"{self.vector.name or 'destination'} vector",
                storage="-20 °C",
                amount=(
                    ""
                    if len(self.block_vectors) < 2
                    else f"{len(self.block_vectors)}, one a position"
                ),
                note=_destination_note(self.block_vectors),
            )
        )
        if self.working is not None:
            made.append(
                Material(
                    f"{self.working.record.name or 'Working'} vector",
                    storage="-20 °C",
                    note=f"the library's final home; {self.working.enzyme.name} releases its ccdB "
                    "cassette to admit the cargo",
                )
            )
            made.append(
                enzyme_amount(
                    self.working.enzyme, "admits the finished cargo to the working vector"
                )
            )
        made.append(enzyme_amount(self.scheme.internal, "opens the library, excising its stuffer"))
        made.append(enzyme_amount(self.scheme.external, "releases a part from its own block"))
        inside, outside = self.inside, self.outside
        for one in self.scheme.blunt:
            jobs = [
                what
                for what, named in (
                    ("the excised stuffer", inside),
                    ("the donor backbone", outside),
                )
                if one in named
            ]
            made.append(enzyme_amount(one, f"cuts {listed(jobs)}, so it cannot ligate back"))
        made.append(Material(CUTSMART, storage="-20 °C", note="both digests run in it"))
        # Both carry their own rules, so neither can appear in a protocol that does not show them.
        made += list(stages.method_materials())
        made.append(
            Material(
                SPRI_BEADS,
                amount=f"{SPRI_AFTER_DIGEST:g} volumes after a digest, "
                f"{SPRI_AFTER_LIGATION:g} after a ligation",
                note="the two ratios differ; both elute in water",
            )
        )
        shot = pulse()
        made.append(
            material(
                STRAIN,
                supplier=STRAIN_SUPPLIER,
                catalog=STRAIN_CATALOG,
                storage="-80 °C",
                amount=f"one aliquot per round, {shot.cells_ul:g} µL a pulse",
                note=f"pulsed at {shot.volts:g} V, {shot.ohms:g} Ω and {shot.microfarads:g} µF "
                f"in a {shot.cuvette_mm:g} mm cuvette",
                citation=shot.citation,
            )
        )
        made.append(Material("Recovery medium", amount="one outgrowth per round"))
        made.append(
            Material("Selective broth and plates", note=_selection_note(self.vector, "destination"))
        )
        if self.working is not None:
            made.append(
                Material(
                    "Selective plates for the final transfer",
                    note=_selection_note(self.working.record, "working"),
                )
            )
        made.append(Material("Plasmid prep kit", amount="one prep per round"))
        made.append(
            Material("Electroporation cuvettes", amount=f"{len(self.part_lists)}, one per round")
        )
        return tuple(made)


def enzyme_amount(enzyme: Enzyme, note: str) -> Material:
    """Return one enzyme as a material, carrying what every digest takes of it."""
    return enzyme_material(enzyme, amount=f"{ENZYME_UL:g} µL per digest", note=note)


def vector_names(destinations: Sequence[tuple[str, str]]) -> list[str]:
    """Say each destination as the bench reads it: what it is called, and the file holding it."""
    return [f"{name} ({file})" if file else name for name, file in destinations]


def with_pair(pair: ReadPair | None, sheet: str) -> str:
    """Name the designed pair, the sheet holding it and its amplicon, or say nothing.

    A pair the bench has to order is a pair the step says where to find, so the sentence names
    the sheet beside the two primers rather than leaving the reader to look for them.
    """
    if pair is None:
        return ""
    from_sheet = f" from {sheet}" if sheet else ""
    return (
        f", on {pair.forward.name} and {pair.reverse.name}{from_sheet}, which give a "
        f"{pair.amplicon_length} bp amplicon"
    )


def as_platform(pair: ReadPair | None) -> str:
    """Say what the amplicon asks of the sequencing, which follows from what it has to carry."""
    return "" if pair is None else f" as a {pair.platform}"


def marks_sentence(constructs: int, marks: RepresentationMarks, what: str) -> str:
    """State the three marks the counts are judged on, and the depth they are judged at."""
    depth = reads_for_representation(constructs, marks.reads_per_member) if constructs else 0
    at = f", which is {depth:,} reads over {constructs:,} combinations" if depth else ""
    return (
        f"At least {marks.seen:.1%} {what} seen; a 90th/10th percentile skew ratio below "
        f"{marks.skew:g}, judged at {marks.reads_per_member} or more reads a member{at}."
    )


def _destination_note(block_vectors: Sequence[tuple[str, str]]) -> str:
    """Say what the destination is for, and that a later position has one of its own."""
    opened = "opened by the internal digest in round 1"
    if len(block_vectors) < 2:
        return opened
    return (
        f"{opened}, and to take each block's cargo. A part enters on its own position's entry "
        f"overhang, so there is one a position: {listed(vector_names(block_vectors))}"
    )


def _selection_note(record: SequenceRecord, which: str) -> str:
    """Return what this vector's transformants are selected on, read off its own marker.

    A record annotating no marker this method can name a drug for leaves the plate to the reader
    rather than naming one, and `stages.ROUND_SELECTION` then stands as the protocol's hole.
    """
    plate = stages.selection_for(record)
    marker = selection_marker(record)
    if not plate or marker is None:
        return f"the {which} vector's own antibiotic, which this plan does not name"
    return f"LB with {plate}, the {which} vector's own marker ({marker.name})"

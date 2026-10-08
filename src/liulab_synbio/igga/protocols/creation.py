"""Protocol 03: pull every block out of the ordered pool and close each cargo into its vector.

Two nested PCRs and one assembly: the first pulls a batch out of the whole pool, the second
pulls one block out of its batch, and the third closes that block's cargo into the vector for
the position it fills. One archived well a design comes out.

This protocol only exists where a pool does; a build ordering its blocks whole makes its
cargo in the vendor's tube.
"""

from collections.abc import Iterable, Mapping, Sequence

from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.gels import choose_ladder
from liulab_mbio.bench.goldengate import (
    KIT,
    assembly_amounts,
    assembly_program,
    assembly_reaction,
)
from liulab_mbio.bench.pcr import pcr_program, pcr_reaction
from liulab_mbio.bench.plates import plate, seat
from liulab_mbio.bench.steps import listed
from liulab_mbio.enzymes import get_enzyme
from liulab_mbio.primers.polymerase import melting_temperature
from liulab_mbio.protocol.figures import ligation_figure
from liulab_mbio.protocol.model import (
    Citation,
    Figure,
    Gel,
    Lane,
    Material,
    Plate,
    Reference,
    Source,
    Step,
    Troubleshooting,
)
from liulab_mbio.protocol.model import (
    Item as Handed,
)
from liulab_synbio.igga.cargo import Batch, PoolPlan
from liulab_synbio.igga.figures import PoolStage, pool_pcr_figure
from liulab_synbio.igga.method import SYNTHESIS_ENZYME, Scheme
from liulab_synbio.igga.protocols.ordering import (
    POOL_POLYMERASE,
    POOL_REFERENCES,
    pool_materials,
)
from liulab_synbio.igga.protocols.primer_plates import working_plate
from liulab_synbio.igga.protocols.protocol import Protocol, figured, labelled
from liulab_synbio.igga.protocols.run import Run, vector_names

#: What the page is headed and what the chain names it by.
CREATION = "Cargo creation"

#: The two nested PCRs that pull one block out of the pool, in the order the steps run them.
POOL_STAGES: tuple[PoolStage, ...] = ("PCR1", "PCR2")

#: How many wells the plate PCR2 runs in holds. One batch is one plate of PCR2, which is what
#: fixes `liulab_synbio.igga.method.ORTHOGONAL_SPLIT`'s 96 inner primers.
PCR2_WELLS = 96
PCR2_PLATE = "PCR2 plate"

#: Twist's own cycle counts for amplifying an oligo pool: each band's longest oligo in nt, then
#: the fewest and the most cycles it allows. FRM-001034 REV 8 p. 2 and DOC-4060 REV 1.0 give the
#: same three. ``docs/research/oligo-pool-pcr-cycles.md`` section 2.
POOL_CYCLE_BANDS: tuple[tuple[int, int, int], ...] = ((100, 6, 10), (150, 10, 12), (350, 12, 14))

#: The document PCR1's cycle count is cited to, and where in it the count stands.
POOL_CYCLE_SOURCE_KEY = "FRM-001034"
POOL_CYCLE_CITATION = Citation(POOL_CYCLE_SOURCE_KEY, "p. 2, cycle chart")
POOL_SOURCES: dict[str, Source] = {
    POOL_CYCLE_SOURCE_KEY: Source(
        "Twist Bioscience, Amplifying Twist Oligo Pools",
        edition="REV 8",
        url="https://www.twistbioscience.com/content/dam/twistbioscience/resources/2026-01/"
        "FRM-001034-AmplifyingOligoPools-REV8%20singles.pdf",
        read_as="plain curl",
        date="2026-10-07",
    )
}

#: The hardware the pool route needs on top of a round's, which no reagent table covers.
POOL_EQUIPMENT: tuple[str, ...] = (
    "Thermocycler taking a 96-well plate",
    "Gel tank and a transilluminator",
)


class Creation(Protocol):
    """Amplify the pool twice and close each block's cargo into its position's vector."""

    round_reagents = True

    def title(self, run: Run) -> str:
        """Return the page's own heading."""
        return CREATION

    def summary(self, run: Run) -> str:
        """Return what the two PCRs and the assembly come to."""
        return (
            "Pull every block out of the ordered pool and close each one's cargo into the "
            "vector for the position it fills, one archived well a design."
        )

    def steps(self, run: Run) -> tuple[Step, ...]:
        """Return the two amplifications and the assembly, each showing what it makes."""
        pool = _pool_of(run)
        sequences = {one.name: one.sequence for one in pool.pool.primers}
        layout = pool.pool.layout
        batches = pool.batches
        first = _annealing(((one.forward, one.outer) for one in batches), sequences)
        second = _annealing(pool.inner_pairs, sequences)
        made = (
            _pcr1_step(batches, first, layout.length),
            _pcr2_step(pool, batches, len(run.parts), second, layout.length - layout.primer_length),
        )
        return (
            *(
                labelled(
                    figured(one, pool_pcr_figure(pool, stage=stage, path=run.oligo)),
                    "Amplify the pool",
                )
                for one, stage in zip(made, POOL_STAGES, strict=True)
            ),
            labelled(
                figured(_assembly_step(run, pool), _cargo_ligation_figure(run)),
                "Close each cargo into its vector",
            ),
        )

    def consumes(self, run: Run) -> tuple[Handed, ...]:
        """Return the pool, the working plate of primers where one was poured, and the block vectors."""
        return (run.ordered, *working_plate(run), *run.blocks)

    def produces(self, run: Run) -> tuple[Handed, ...]:
        """One well a design, sealed and frozen."""
        return (run.archive,)

    def carried(self, run: Run) -> tuple[Material, ...]:
        """Return what the pool route buys, which the ordering protocol bought most of."""
        return pool_materials(_pool_of(run), run.pool_sheet, run.primer_sheet)

    def equipment(self, run: Run) -> tuple[str, ...]:
        """Return what the pool route needs on top of a round's."""
        return POOL_EQUIPMENT

    def plates(self, run: Run) -> tuple[Plate, ...]:
        """Return the plate PCR2 runs in, one block a well."""
        return (_pcr2_plate(_pool_of(run)),)

    def references(self, run: Run) -> tuple[Reference, ...]:
        """Return the pool's own references, which the protocol before it printed too."""
        return POOL_REFERENCES

    def sources(self, run: Run) -> dict[str, Source]:
        """Return the document PCR1's cycle count is cited to."""
        return POOL_SOURCES


def pool_cycles(length_nt: int) -> tuple[int, int]:
    """Return the fewest and the most cycles Twist allows a pool of this length.

    The count is banded by length, so it follows the build's own oligo length rather than
    sitting fixed in the source. A pool longer than the last band takes that band's count: the
    table stops where the product does.

    Examples
    --------
    >>> pool_cycles(350)
    (12, 14)
    >>> pool_cycles(60)
    (6, 10)
    >>> pool_cycles(400)
    (12, 14)
    """
    for longest, fewest, most in POOL_CYCLE_BANDS:
        if length_nt <= longest:
            return (fewest, most)
    _, fewest, most = POOL_CYCLE_BANDS[-1]
    return (fewest, most)


def _pool_of(run: Run) -> PoolPlan:
    """Return the run's pool, which this protocol is only written for a run that has one.

    Raises
    ------
    ValueError
        If the run orders its blocks whole, where this protocol does not belong in the chain.
    """
    if run.pool is None:
        raise ValueError("this run orders its blocks whole and makes no cargo from a pool")
    return run.pool


def _pcr2_plate(pool: PoolPlan) -> Plate:
    """Return the plate PCR2 runs in, one block a well, drawn as the first batch fills it.

    The wells are labelled rather than seated: a block is no material, oligo, vessel or plate,
    so its name describes the well instead of naming an occupant the protocol declares.
    """
    batches = pool.batches
    return plate(
        PCR2_PLATE,
        PCR2_WELLS,
        holds="one block a well, in the order its batch allotted the inner primers",
        labels=seat(batches[0].blocks, PCR2_WELLS) if batches else {},
        note=(
            f"{_counted(len(batches), 'plate')}, one a batch; "
            f"{', '.join(f'batch {one.number} holds {len(one.blocks)}' for one in batches)}"
        ),
    )


def _annealing(
    pairs: Iterable[tuple[str, str]], sequences: Mapping[str, str]
) -> tuple[float, float]:
    """Return the lowest and the highest annealing temperature these primer pairs ask for, °C."""
    found = [
        POOL_POLYMERASE.annealing_temperature(
            melting_temperature(sequences[one], POOL_POLYMERASE),
            melting_temperature(sequences[other], POOL_POLYMERASE),
        )
        for one, other in pairs
    ]
    return min(found), max(found)


def _counted(number: int, word: str) -> str:
    """Say a count and its noun, the noun plural only where the count is not one."""
    return f"{number} {word}" if number == 1 else f"{number} {word}s"


def _one_band(name: str, length_bp: int, *, title: str) -> Gel:
    """Return the gel one lane of a PCR makes, which is the same lane in every tube of it."""
    return Gel(choose_ladder((length_bp,)), (Lane(name, (length_bp,)),), title=title)


def _band_note(low: float, high: float) -> str:
    """Say what annealing temperature a block of tubes runs at, and the spread it covers."""
    if low == high:
        return f"Every pair anneals at {low:g} °C, so one block of tubes takes them all."
    return (
        f"The pairs anneal between {low:g} and {high:g} °C. The program runs at {low:g}, the "
        "lowest of them, so one block of tubes takes them all."
    )


def _pcr1_step(batches: Sequence[Batch], annealing: tuple[float, float], length_bp: int) -> Step:
    """Pull one batch of blocks out of the whole pool, which is what PCR1 is for."""
    from liulab_synbio.igga import stages

    low, high = annealing
    fewest, most = pool_cycles(length_bp)
    pairs = "; ".join(f"batch {one.number}: {one.forward} with {one.outer}" for one in batches)
    return Step(
        f"PCR1: pull {_counted(len(batches), 'batch')} out of the pool",
        key="pcr1",
        instructions=(
            f"Set up {_counted(len(batches), 'reaction')}, one a batch, with the pool as template.",
            f"Give each its own pair: {pairs}.",
            "Run the program below. On a real-time instrument, add an intercalating dye and stop "
            "before the curve plateaus; the printed count is what to run without one.",
        ),
        tables=(
            pcr_reaction(POOL_POLYMERASE, reactions=len(batches), title="PCR1, one tube a batch"),
        ),
        programs=(
            pcr_program(
                POOL_POLYMERASE,
                annealing_temperature=low,
                amplicon_length=length_bp,
                cycles=fewest,
                cycles_citation=POOL_CYCLE_CITATION,
                title="PCR1",
            ),
        ),
        gels=(_one_band("PCR1 product", length_bp, title="PCR1, any batch"),),
        expected=(
            f"One band at {length_bp} bp in every batch, which is the whole oligo.",
            f"Batch sizes: "
            f"{', '.join(f'{one.number} holds {len(one.blocks)}' for one in batches)}.",
        ),
        notes=(
            _band_note(low, high),
            f"Twist's band for a {length_bp} nt pool is {fewest} to {most} cycles; Twist Oligo "
            "Pools Amplification Protocol DOC-4060 REV 1.0 gives the same three bands. Its FAQ "
            f"answers that more cycles give worse uniformity, so {fewest} is what prints.",
            "The outer primer is what makes a batch a batch: it is dropped at PCR2, so a block "
            "cannot be pulled out of a batch it does not sit in.",
        ),
        holes=(stages.PCR1_POLYMERASE,),
        troubleshooting=(
            Troubleshooting(
                "No band",
                "Drop the annealing temperature by 3 °C and check the pool went in.",
            ),
            Troubleshooting(
                "A hump after the peak on a Bioanalyzer or TapeStation trace",
                "Heteroduplexes, which is what over-amplification leaves. Run it again with "
                "fewer cycles. What too many cycles cost is dropout, chimeras and polymerase "
                "error; pooling the blocks equimolar later restores the evenness but none of "
                "those.",
            ),
        ),
    )


def _pcr2_step(
    pool: PoolPlan,
    batches: Sequence[Batch],
    blocks: int,
    annealing: tuple[float, float],
    length_bp: int,
) -> Step:
    """Pull one block out of its batch, after which a block is named by the well it sits in.

    The program's cycle count is blank. Twist's band is for amplifying the pool as it arrives,
    and this reaction's template is PCR1's product, so the reader is given the rule that stops
    the reaction and `liulab_synbio.igga.stages.PCR2_CYCLES` stands where the number would be.
    """
    from liulab_synbio.igga import stages

    low, high = annealing
    return Step(
        f"PCR2: pull each of the {blocks} blocks out of its batch",
        key="pcr2",
        instructions=(
            f"Set up one reaction a block, {blocks} in all, in {PCR2_PLATE}.",
            "Give each its batch's PCR1 product as template, that batch's forward primer, and "
            "the block's own inner primer.",
            "Run the program below. Nobody published a cycle count for this reaction, so run it "
            "on a real-time instrument with an intercalating dye and stop before the curve "
            "plateaus.",
        ),
        tables=(pcr_reaction(POOL_POLYMERASE, reactions=blocks, title="PCR2, one well a block"),),
        programs=(
            pcr_program(
                POOL_POLYMERASE,
                annealing_temperature=low,
                amplicon_length=length_bp,
                # A blank count has nothing to cite; `stages.PCR2_CYCLES` stands where it would.
                cycles=None,
                cycles_citation=None,
                title="PCR2",
            ),
        ),
        gels=(_one_band("PCR2 product", length_bp, title="PCR2, any well"),),
        expected=(
            f"One band at {length_bp} bp in every well, the outer primer's "
            f"{pool.pool.layout.primer_length} nt shorter than PCR1's.",
            f"{blocks} wells filled across {_counted(len(batches), 'plate')}.",
        ),
        notes=(
            _band_note(low, high),
            "A block is named by its well from here on, not by anything in the tube.",
        ),
        holes=(stages.PCR2_CYCLES,),
        troubleshooting=(
            Troubleshooting(
                "A band at PCR1's length",
                "The outer primer region was not dropped. Check the inner primer went in, and "
                "that the template is PCR1's product and not the pool.",
            ),
            Troubleshooting(
                "A hump after the peak on a Bioanalyzer or TapeStation trace",
                "Heteroduplexes, which is what over-amplification leaves. Run it again with "
                "fewer cycles; the forward primer is the whole batch's, so over-cycling here "
                "also pulls a neighbour's product into the well.",
            ),
        ),
    )


def _assembly_step(run: Run, pool: PoolPlan) -> Step:
    """Clone each block's cargo into its position's opened destination, one well a block.

    The destination supplies the stuffers the cargo is no longer synthesised with, and its own
    marker is what selects a well that closed. A part enters on its own position's entry
    overhang, so there is one destination a position and a block goes into its own. The reaction
    and its cycling are NEB's kit table for the most pieces any block takes, so one master mix
    covers the plate.
    """
    scheme = run.scheme
    opened: Amount = run.opened
    enzyme = get_enzyme(SYNTHESIS_ENZYME)
    pieces = {split.pieces for split in pool.splits}
    most = max(pieces)
    piece_bp = pool.pool.layout.length - pool.pool.layout.primer_length
    cargo = [len(split.cargo.sequence) for split in pool.splits]
    amounts = assembly_amounts(
        (opened.name, opened.length_bp),
        tuple((f"PCR2 piece {number}", piece_bp) for number in range(1, most + 1)),
    )
    named = listed(vector_names(run.destinations))
    plated = run.selection or "the destinations' own marker, which this plan does not name"
    return Step(
        f"Assemble each cargo into its position's destination, from its {min(pieces)} to "
        f"{most} pieces",
        key="assemble-cargo",
        instructions=(
            f"Open {named} with {scheme.internal.name} and "
            f"{listed([one.name for one in run.inside])}, the digest a round opens the "
            "library with, and clean each up.",
            "Set up one assembly a well, holding that block's own PCR2 pieces, the opened "
            "destination of the position that block fills, and nothing from another well.",
            f"Transform, plate on {plated}, and pick one colony a block.",
        ),
        tables=(assembly_reaction(enzyme, amounts, system=KIT, reactions=len(run.parts)),),
        programs=(assembly_program(enzyme, fragments=most + 1, system=KIT),),
        expected=(
            f"{len(run.parts)} plasmids, one a block: its position's destination carrying that "
            f"block's cargo, {min(cargo):,} to {max(cargo):,} bp of it.",
            f"Each cargo reads: the overhang its part enters on, its coding bases, the internal "
            f"stuffer, its barcode, and the {scheme.cloning_scar} cloning scar. The external "
            f"stuffers {run.sheet} spells either side of it are the destination's own bases.",
        ),
        notes=(
            f"{SYNTHESIS_ENZYME} cuts each oligo back to its fragment, so the primer sites, the "
            "recognition sites and the padding all stay outside the cargo.",
            "The table is sized at the most pieces any block takes, so one master mix covers "
            "the plate; a well with fewer pieces fills fewer of its DNA rows.",
            f"A part enters on its own position's entry overhang, so there is one destination a "
            f"position. They are one vector differing in the {len(scheme.entry_overhang)} bases "
            "of that overhang; a block put in the wrong one cannot close.",
        ),
        troubleshooting=(
            Troubleshooting(
                "A block comes out short",
                "A piece was missing from the well. Check that well's PCR2 lane before "
                "assembling it again; the pieces of one block are not interchangeable.",
            ),
            Troubleshooting(
                "Colonies carrying the destination with no cargo",
                "The destination was not opened to completion. Run a little of its digest on a "
                "gel before setting the plate up again.",
            ),
        ),
    )


def _cargo_ligation_figure(run: Run) -> Figure | None:
    """Return the stuffer a cargo replaces, at base level: where the assembly joins it."""
    found = run.block_vector
    if found is None:
        return None
    vector, name, path = found
    scheme: Scheme = run.scheme
    return ligation_figure(
        vector.record,
        path=path,
        junction=(vector.stuffer.start, vector.stuffer.end),
        enzymes=(scheme.internal.name, SYNTHESIS_ENZYME),
        caption=(
            f"The stuffer {scheme.internal.name} cuts out of {name}, which a cargo replaces. "
            "The cargo enters on the four bases its position spells here."
        ),
    )

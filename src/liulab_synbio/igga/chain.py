"""A library run as the chain of protocols the bench works through, in order.

This module holds the order and nothing else about any one protocol: which protocols a run
includes, how the run's reagents spread across them, and what no single protocol owns — the
background a reader is told first, the bill for the whole run, and the checks that judge the
design. Each protocol is a module of `liulab_synbio.igga.protocols`, and adding a fact to one
is an edit there.

The chain is the order someone does the work in: plate the primers, order the blocks, make the
cargo, read back the designs the build asks for, join the part lists round by round, and move
the finished library into a working vector.
"""

from collections.abc import Iterable, Sequence

from liulab_mbio.bench.prices import Item
from liulab_mbio.bench.prices import bill as priced
from liulab_mbio.bench.steps import badges
from liulab_mbio.protocol.figures import SOURCE as FIGURE_SOURCE
from liulab_mbio.protocol.figures import SOURCE_KEY as FIGURE_SOURCE_KEY
from liulab_mbio.protocol.model import (
    Bill,
    Material,
    Project,
    Source,
    Step,
    Topic,
    by_place,
    counted,
    names,
)
from liulab_mbio.protocol.model import Item as Handed
from liulab_mbio.protocol.model import Protocol as Page
from liulab_synbio import dmx
from liulab_synbio.igga import stages
from liulab_synbio.igga.bench import ENZYME_UL, STRAIN, STRAIN_CATALOG, choppers
from liulab_synbio.igga.cargo import PoolPlan
from liulab_synbio.igga.protocols.assembly import Assembly
from liulab_synbio.igga.protocols.creation import Creation
from liulab_synbio.igga.protocols.final import FinalLigation
from liulab_synbio.igga.protocols.ordering import Ordering
from liulab_synbio.igga.protocols.primer_plates import PrimerPlating
from liulab_synbio.igga.protocols.protocol import Protocol
from liulab_synbio.igga.protocols.run import BLOCK_VECTOR_ITEM, WORKING_ITEM, Run
from liulab_synbio.igga.protocols.validation import ReadBack

#: What a price record prices the synthesis order and the bill's own source by. Neither has a
#: catalogue number, so each is keyed by the name the protocol's own row already carries. The
#: pool's own key is the build's, through `liulab_synbio.igga.cargo.design_pool`.
BLOCKS_KEY = "synthesised blocks"
POOL_PRIMER_KEY = "pool-primers"
PRICES_SOURCE = "prices"


def ordered(run: Run) -> tuple[Protocol, ...]:
    """Return the protocols this run writes, in the order the bench works through them.

    A run plating no primers opens at the ordering protocol; one ordering its blocks whole
    makes its cargo in the vendor's tube and writes no creation protocol; one stating no
    fragment-count floor reads nothing back. A run naming more than one route writes one
    read-back protocol per route, next to each other, because they are the ways of one job.
    """
    made: list[Protocol] = []
    if run.plated:
        made.append(PrimerPlating())
    made.append(Ordering())
    if run.pool:
        made.append(Creation())
    made += [ReadBack(one) for one in run.validations]
    made += [Assembly(), FinalLigation()]
    return tuple(made)


def project(run: Run) -> Project:
    """Return the run as the chain of protocols the bench works through, in order.

    Each protocol says what it is handed and what it leaves, and nothing else about the run.

    The bill and the design checks are the run's, not one protocol's: two protocols buying the
    same cells would otherwise be counted twice, and no one protocol judges the design. Money
    comes only from `run.prices`, and every row it does not price carries a hole. The record
    those rows cite is the run's own source, since the bill sits on a page no protocol owns.
    """
    made = ordered(run)
    staged = [one.steps(run) for one in made]
    spread = _spread(run, made, staged)
    sources = _sources(run, made)
    pages = tuple(
        one.page(run, steps=steps, materials=bought, sources=sources)
        for one, steps, bought in zip(made, staged, spread, strict=True)
    )
    return Project(
        f"{run.vector.name or 'Library'}: {len(run.part_lists)} part lists in "
        f"{len(run.rounds)} rounds",
        summary=(
            f"Join {len(run.parts)} synthesised parts into {run.constructs:,} distinct constructs "
            f"in {len(run.rounds)} rounds, over {_places(pages)} protocols. Each round opens "
            f"the library with {run.scheme.internal.name}, releases one part list with "
            f"{run.scheme.external.name}, ligates the two, and transforms, grows and preps the "
            "result for the round after it."
        ),
        background=(Topic("How this library is designed", _highlights(run)), *_choosing(run)),
        files=run.files,
        inputs=_inputs(run),
        protocols=pages,
        checks=badges(run.checks),
        sources={PRICES_SOURCE: run.prices.source} if run.prices else {},
        bill=_consumed(run),
    )


def _places(pages: Sequence[Page]) -> int:
    """Return how many protocols the bench works through, which is fewer than the pages.

    The ways of one job stand in one place and the bench does one of them.
    """
    return len(by_place(pages, lambda one: one.choice))


def _choosing(run: Run) -> tuple[Topic, ...]:
    """Return what to weigh where the bench picks a route, and nothing where it does not.

    The topic is `dmx`'s: comparing that method's own routes is the method's knowledge.
    """
    return (dmx.route_choice(),) if run.read_back_is_a_choice else ()


def _spread(
    run: Run, made: Sequence[Protocol], staged: Sequence[Sequence[Step]]
) -> list[tuple[Material, ...]]:
    """Return each protocol's reagent list: what its own steps use, and nothing else.

    One list a protocol, in the chain's own order. A reagent the run shares goes to the
    protocols whose steps name it and to the ones that claim it in `Protocol.shares`; a
    reagent one protocol buys alone goes to it, or to the protocols naming it where another
    does. A page can therefore list nothing no step of it reaches. The lists are keyed by
    where a protocol stands in the chain and not by its heading, so two protocols could share
    a heading without sharing a reagent.
    """
    found: list[list[Material]] = [[] for _ in made]

    def naming(material: Material) -> list[int]:
        return [
            place
            for place, steps in enumerate(staged)
            if any(names(material.name, step.named) for step in steps)
        ]

    for place, one in enumerate(made):
        for material in one.carried(run):
            _bought(found, material, naming(material) or [place])
    for material in run.round_materials:
        claiming = [
            place
            for place, one in enumerate(made)
            if any(names(claim, (material.name,)) for claim in one.shares(run) if claim)
        ]
        _bought(found, material, sorted({*naming(material), *claiming}))
    return [tuple(group) for group in found]


def _bought(found: list[list[Material]], material: Material, places: Iterable[int]) -> None:
    """Put `material` on each of those protocols' lists, once each."""
    for place in places:
        if material not in found[place]:
            found[place].append(material)


def _sources(run: Run, made: Sequence[Protocol]) -> dict[str, Source]:
    """Return every document this run could cite; `citing` drops the ones a protocol did not.

    Every protocol is handed the whole set, because one of them draws a figure another's
    module sources and a reference the reader follows back has to resolve either way.
    """
    found = dict(stages.SOURCES) | {FIGURE_SOURCE_KEY: FIGURE_SOURCE}
    for one in made:
        found |= one.sources(run)
    return found


def _inputs(run: Run) -> tuple[Handed, ...]:
    """Return what the bench holds before the first protocol: stock the run does not make."""
    made = [
        Handed(
            BLOCK_VECTOR_ITEM.format(number=number),
            f"{name}, which a block of position {number} closes into",
            storage="-20 °C",
        )
        for number, (name, _) in enumerate(run.block_vectors, 1)
    ] or [
        Handed(
            BLOCK_VECTOR_ITEM.format(number=1),
            f"{run.vector.name or 'the destination'}, which round 1 opens",
            storage="-20 °C",
        )
    ]
    made += run.marking_stocks
    if run.working is not None:
        made.append(
            Handed(
                WORKING_ITEM,
                f"{run.working.record.name or 'the backbone'} the finished library ends in",
                storage="-20 °C",
            )
        )
    return tuple(made)


def _highlights(run: Run) -> tuple[str, ...]:
    """Return what the facts mean, a sentence each."""
    scheme, standard, destination = run.scheme, run.standard, run.destination
    last = run.bench[-1]
    said = [
        f"One round appends one part list to every member of the library at once, so "
        f"{len(run.rounds)} rounds make {run.constructs:,} constructs out of "
        f"{sum(one.coverage.part_list_size for one in run.bench)} synthesised parts.",
        f"The product keeps {scheme.internal.name}'s sites, which is what lets the next round "
        f"open it, and loses {scheme.external.name}'s with the external stuffers.",
    ]
    if destination.edit is None:
        said.append(
            f"Your vector already carried an internal stuffer, so every part enters position "
            f"{run.positions[0]} on {standard.entry_overhangs[0]}, which is the overhang "
            "that stuffer spells: DNA that exists cannot be re-chosen."
        )
    else:
        said.append(
            f"Your vector carried no internal stuffer, so one was put in at "
            f"{destination.stuffer.start} carrying {standard.entry_overhangs[0]}, the overhang "
            "this design chose. Read the edit before you order the vector."
        )
    if standard.cost:
        said.append(
            f"The overhang standard changes {counted(standard.cost, 'amino acid')} across "
            f"{counted(len(standard.changes), 'part end')}; {run.changes} has wild type beside "
            "synthesised, and that is what you are about to pay for."
        )
    else:
        said.append(
            "The overhang standard changes no amino acid: every part list already spells its "
            "junctions."
        )
    said.append(
        f"The last round needs {last.coverage.colonies:,} colonies for a "
        f"{last.coverage.completeness:g} chance that none of its {last.coverage.products:,} "
        f"products is missing, which is {last.coverage.coverage:.0f}x its products. A round "
        "short of that floor loses members no later round can put back."
    )
    said.append(
        f"The finished barcode block is "
        f"{scheme.barcode_block_length(run.barcode_length, len(run.positions))} bp and reads in "
        f"the reverse of the order the rounds ran; {run.barcodes} says which barcode names which "
        "part."
    )
    return tuple(said)


def _consumed(run: Run) -> Bill:
    """Return what this build buys, in the quantities the design computes.

    Only what the design fixes is billed. The buffer, the beads and the medium scale with volumes
    the method leaves to the supplier, so a row for one would be a quantity nobody computed.

    **The blocks are bought once or not at all.** A project with a pool buys oligos and the
    primers that amplify them, and assembles its blocks from those; one without buys the blocks
    themselves. Billing both would charge the same DNA twice.
    """
    scheme, pool, record = run.scheme, run.pool, run.prices
    inside, outside = choppers(scheme)
    rounds = len(run.rounds)
    digests = [(scheme.internal, rounds), (scheme.external, rounds)]
    digests += [(one, rounds * ((one in inside) + (one in outside))) for one in scheme.blunt]
    items = [
        *(
            (
                Item(
                    "Synthesised blocks",
                    len(run.parts),
                    unit="blocks",
                    key=BLOCKS_KEY,
                    quantities={
                        "count": len(run.parts),
                        "length_nt": max(one.length for one in run.parts),
                    },
                ),
            )
            if pool is None
            else (pool.item, _primer_item(pool))
        ),
        *(
            Item(
                one.commercial_name or one.name,
                ENZYME_UL * uses,
                unit="µL",
                key=one.catalog_number or "",
                quantities={"volume_ul": ENZYME_UL * uses},
            )
            for one, uses in digests
        ),
        Item(STRAIN, rounds, unit="aliquots", key=STRAIN_CATALOG, quantities={"count": rounds}),
        Item(
            "Electroporation cuvettes",
            rounds,
            unit="cuvettes",
            key="cuvettes",
            quantities={"count": rounds},
        ),
        Item(
            "Plasmid preps", rounds, unit="preps", key="plasmid prep", quantities={"count": rounds}
        ),
    ]
    return priced(items, record, source_key=PRICES_SOURCE)


def _primer_item(pool: PoolPlan) -> Item:
    """Return the amplification primers as one line of the bill, beside the pool's own."""
    primers = pool.pool.primers
    return Item(
        "Pool amplification primers",
        len(primers),
        unit="primers",
        key=POOL_PRIMER_KEY,
        quantities={"count": len(primers), "length": max(len(one) for one in primers)},
    )

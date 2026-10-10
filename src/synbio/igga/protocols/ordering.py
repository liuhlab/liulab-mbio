"""Protocol 03: order the synthesised material a library is built from, and put it away.

A pool splits ordering from making the cargo: the pool is ordered and stored in one sitting,
and pulled apart into blocks in another. Without one the blocks are bought whole and this is
the only protocol before the rounds.

What the pool route buys sits here, because this is the protocol that buys it: the polymerase
the two amplifications run in, the primers that pull a block out, and the documents both
protocols cite.
"""

from mbio.bench.gels import SOURCES as GEL_SOURCES
from mbio.bench.pcr import SOURCES as PCR_SOURCES
from mbio.bench.pcr import polymerase_name
from mbio.bench.steps import catalogued, listed
from mbio.primers.polymerase import Q5, Polymerase
from mbio.protocol.model import Citation, Material, Note, Source, Step, Troubleshooting, sectioned
from mbio.protocol.model import Item as Handed
from synbio.igga.cargo import PoolPlan
from synbio.igga.method import SYNTHESIS_ENZYME
from synbio.igga.protocols.protocol import Protocol
from synbio.igga.protocols.run import Run

#: What the page is headed and what the chain names it by.
ORDERING = "Cargo ordering and pool preparation"

#: What `pool_materials` calls the primers that amplify the pool. Where a run plates them,
#: they are ordered in that protocol instead and this material is only used in the next one.
PRIMER_MATERIAL = "Pool amplification primers"

#: What `pool_materials` calls the pool. The bill names it by the pool's own name instead, so
#: the chain matches the two by this rather than by spelling.
POOL_MATERIAL = "Oligo pool"

#: What amplifies the pool. The package's high-fidelity default, named here so the two PCRs
#: and the material agree about which buffer the annealing temperatures were computed in.
POOL_POLYMERASE: Polymerase = Q5
#: The tube under the one name the reaction pipettes it by, NEB's M0491.
POOL_POLYMERASE_PRODUCT = f"{polymerase_name(POOL_POLYMERASE)} (M0491)"

#: What the dried pool is dissolved in, and the least it may be left at, ng/µL. Both are the
#: vendor's, from the document the cycle count is read in: DOC-4060 REV 1.0, "Before You Begin".
POOL_STOCK_BUFFER = "10 mM Tris buffer, pH 8.0"
POOL_STOCK_NG_PER_UL = 20.0

#: The document the dried pool's own buffer and floor are read from. It carries the cycle
#: bands too, so the protocol that amplifies the pool cites it as well.
POOL_STOCK_SOURCE_KEY = "DOC-4060"

#: Where the pool's buffer and floor stand in it, and where its answer on cycle count does.
POOL_STOCK_CITATION = Citation(POOL_STOCK_SOURCE_KEY, "Before You Begin")
POOL_UNIFORMITY_CITATION = Citation(POOL_STOCK_SOURCE_KEY, "Appendix B")

#: What the pool route's own rows and sentences are cited to, beside the round's own.
POOL_SOURCES: dict[str, Source] = {
    **PCR_SOURCES,
    **GEL_SOURCES,
    POOL_STOCK_SOURCE_KEY: Source(
        "Twist Bioscience, Twist Oligo Pools Amplification Protocol",
        edition="DOC-4060 REV 1.0",
        read_as="held under reference_docs/",
        note="docs/research/oligo-pool-pcr-cycles.md",
    ),
}


class Ordering(Protocol):
    """Order the pool, or every block, and store it ready for the bench."""

    def title(self, run: Run) -> str:
        """Return the page's own heading."""
        return ORDERING

    def summary(self, run: Run) -> str:
        """Return what ordering comes to, which is more than one line only where there is a pool."""
        if run.pool:
            return (
                "Order the synthesised material this library is built from, and put it away "
                "ready for the bench."
            )
        return "Order every block this library is built from."

    def steps(self, run: Run) -> tuple[Step, ...]:
        """Return the order itself, and for a pool the resuspension that readies it."""
        if run.pool:
            made = (
                _pool_order_step(run.pool, run.pool_sheet, run.primer_sheet, run.plated),
                _pool_resuspend_step(),
            )
        else:
            made = (_order_step(run),)
        return sectioned("Order and store", *made)

    def produces(self, run: Run) -> tuple[Handed, ...]:
        """Return the pool, or the blocks the vendor shipped."""
        return (run.ordered,)

    def carried(self, run: Run) -> tuple[Material, ...]:
        """Return what the pool route buys here; a run that plates its primers buys them there."""
        if run.pool is None:
            return ()
        bought = pool_materials(run.pool, run.pool_sheet, run.primer_sheet)
        if run.plated:
            return tuple(one for one in bought if one.name != PRIMER_MATERIAL)
        return bought

    def sources(self, run: Run) -> dict[str, Source]:
        """Return the documents this page could cite.

        Ordering a pool and resuspending it runs no PCR and pours no gel, so the cycle count,
        the ladders and the gel resolution are cited by the protocol after it and not here.
        """
        return dict(POOL_SOURCES)


def pool_materials(pool: PoolPlan, pool_sheet: str, primer_sheet: str) -> tuple[Material, ...]:
    """Return what the pool route buys that a build ordering its blocks does not."""
    layout = pool.pool.layout
    roles = ", ".join(f"{count} {role}" for role, count in primer_roles(pool).items())
    return (
        Material(
            POOL_MATERIAL,
            storage="-20 °C",
            amount=f"{pool.pool.count} oligos, every one {layout.length} nt, resuspended to "
            f"{POOL_STOCK_NG_PER_UL:g} ng/µL",
            note=f"ordered from {pool_sheet}, which says which block each oligo is a piece of",
        ),
        Material(
            PRIMER_MATERIAL,
            storage="-20 °C",
            amount=f"{len(pool.pool.primers)} primers: {roles}",
            note=f"ordered from {primer_sheet}; a pair a batch, an inner primer a block",
        ),
        catalogued(
            POOL_POLYMERASE_PRODUCT,
            supplier="NEB",
            storage="-20 °C",
            amount="one reaction per batch, then one per block",
            note="the buffer the annealing temperatures below were computed in",
        ),
        catalogued(
            "NEBridge Golden Gate Assembly Kit, BsmBI-v2 (E1602)",
            supplier="NEB",
            storage="-20 °C",
            amount="one assembly per block",
            note=f"its mix carries the {SYNTHESIS_ENZYME} that cuts every oligo back to its own "
            "fragment, which this method reserves, and the T4 DNA Ligase that joins them",
        ),
    )


def primer_roles(pool: PoolPlan) -> dict[str, int]:
    """How many primers of the pool take each role, in the order the set allots them."""
    counted: dict[str, int] = {}
    for one in pool.pool.primers:
        counted[one.role] = counted.get(one.role, 0) + 1
    return counted


def _order_step(run: Run) -> Step:
    """Order every block, which is what the whole design comes down to."""
    lengths = [part.length for part in run.parts]
    sheet = run.sheet
    return Step(
        "Order the synthesised parts",
        key="order-blocks",
        instructions=(
            f"Order every row of {sheet} as a double-stranded synthesised block.",
            "Ask for sequence-verified material: a block carries the scheme's sites at fixed "
            "offsets, and a base out of place stops it being cut where the design says.",
        ),
        expected=(
            f"{len(run.parts)} blocks, {min(lengths)} to {max(lengths)} bp.",
            "Each block reads: 5' external stuffer, coding bases, internal stuffer, barcode, "
            "3' external stuffer.",
        ),
        notes=(
            f"{sheet} carries each block's barcode on its own row, so the sheet you order from "
            "is also what decodes the sequencing afterwards.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The vendor cannot synthesise a block",
                "It is usually a repeat or a GC run in the coding bases. Plan again with another "
                "codon usage table; the stuffers and the overhangs are fixed by the scheme.",
            ),
        ),
    )


def _pool_order_step(pool: PoolPlan, pool_sheet: str, primer_sheet: str, plated: bool) -> Step:
    """Order the pool, and its primers where no other protocol has ordered them already."""
    layout = pool.pool.layout
    roles = ", ".join(f"{count} {role}" for role, count in primer_roles(pool).items())
    over = pool.over_floor
    spent = (
        f"{pool.floor} oligos is what the length budget allows for; this design spends "
        f"{pool.pool.count}, one more on {listed(list(over))}, whose forced cuts spell no legal "
        "overhang."
        if over
        else f"{pool.pool.count} oligos is the fewest the length budget allows for."
    )
    return Step(
        "Order the oligo pool"
        if plated
        else "Order the oligo pool and the primers that amplify it",
        key="order-pool",
        instructions=(
            f"Order every row of {pool_sheet} as one synthesised oligo pool, {pool.pool.count} "
            f"members at {layout.length} nt.",
            *(
                ()
                if plated
                else (f"Order every row of {primer_sheet} as an ordinary oligo: {roles}.",)
            ),
        ),
        expected=(
            f"One pool, {pool.pool.count} oligos, every one {layout.length} nt and no length "
            "spread, which is what the padding is for.",
            *(
                ()
                if plated
                else (
                    f"{len(pool.pool.primers)} primers, each "
                    f"{min(len(one) for one in pool.pool.primers)} nt.",
                )
            ),
        ),
        notes=(
            f"{pool_sheet} names each oligo's block, which piece of it that is, and the primers "
            "that pull it out, so an oligo is traceable to its protein.",
            spent,
        ),
        troubleshooting=(
            Troubleshooting(
                "The vendor bins the pool by length",
                "Every oligo is padded to one length, so the spread is zero. A vendor asking "
                "for a different length wants this run's oligo length changed and the "
                "design run again.",
            ),
        ),
    )


def _pool_resuspend_step() -> Step:
    """Put the dried pool into buffer, which is what both amplifications pipette from."""
    floor = f"{POOL_STOCK_NG_PER_UL:g}"
    return Step(
        "Resuspend the oligo pool",
        key="resuspend-pool",
        instructions=(
            f"Divide the total yield in ng printed on the shipping tube label by {floor}, "
            "rounding down, to get the resuspension volume in µL.",
            f"Add that volume of {POOL_STOCK_BUFFER} to the tube.",
            "Vortex the tube until nothing is left undissolved.",
        ),
        cautions=("Spin the tube down before taking the cap off.",),
        expected=(
            f"One tube of pool in solution at {floor} ng/µL or above, with nothing left "
            "undissolved on the wall of the tube.",
        ),
        notes=(
            Note(
                f"Dividing by {floor} is what the floor of at least {floor} ng/µL comes to, "
                f"and the amplification then pipettes 1 µL of {floor} ng/µL.",
                citation=POOL_STOCK_CITATION,
            ),
        ),
    )

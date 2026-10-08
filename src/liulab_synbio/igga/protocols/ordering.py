"""Protocol 02: order the synthesised material a library is built from, and put it away.

A pool splits ordering from making the cargo: the pool is ordered and stored in one sitting,
and pulled apart into blocks in another. Without one the blocks are bought whole and this is
the only protocol before the rounds.

What the pool route buys sits here, because this is the protocol that buys it: the polymerase
the two amplifications run in, the primers that pull a block out, and the references both
protocols print.
"""

from liulab_mbio.bench.gels import REFERENCES as GEL_REFERENCES
from liulab_mbio.bench.pcr import REFERENCES as PCR_REFERENCES
from liulab_mbio.bench.steps import catalogued, listed
from liulab_mbio.primers.polymerase import Q5, Polymerase
from liulab_mbio.protocol.model import Item as Handed
from liulab_mbio.protocol.model import Material, Reference, Step, Troubleshooting
from liulab_synbio.igga.cargo import PoolPlan
from liulab_synbio.igga.method import SYNTHESIS_ENZYME
from liulab_synbio.igga.protocols.protocol import Protocol, labelled
from liulab_synbio.igga.protocols.run import Run

#: What the page is headed and what the chain names it by.
ORDERING = "Cargo ordering and pool preparation"

#: What `pool_materials` calls the primers that amplify the pool. Where a run plates them,
#: they are ordered in that protocol instead and this material is only used in the next one.
PRIMER_MATERIAL = "Pool amplification primers"

#: What amplifies the pool. The package's high-fidelity default, named here so the two PCRs
#: and the material agree about which buffer the annealing temperatures were computed in.
POOL_POLYMERASE: Polymerase = Q5
POOL_POLYMERASE_PRODUCT = "Q5 High-Fidelity DNA Polymerase (M0491)"

#: Where PCR1's cycle count is read. Two independently revised Twist documents give the same
#: three length bands, and the second's appendix answers what more cycles cost.
POOL_CYCLE_REFERENCES: tuple[Reference, ...] = (
    Reference(
        "Twist Bioscience, Amplifying Twist Oligo Pools, FRM-001034 REV 8, p. 2, for the cycle "
        "count banded by the pool's length"
    ),
    Reference(
        "Twist Bioscience, Twist Oligo Pools Amplification Protocol, DOC-4060 REV 1.0, for the "
        "same three bands, and for the FAQ answering that more cycles give worse uniformity"
    ),
)

#: What the two PCRs and their gel cite, beside the round's own references.
POOL_REFERENCES: tuple[Reference, ...] = (
    *PCR_REFERENCES,
    *GEL_REFERENCES,
    *POOL_CYCLE_REFERENCES,
)


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
        """Return the one step this protocol is: the order itself."""
        if run.pool:
            made = _pool_order_step(run.pool, run.pool_sheet, run.primer_sheet, run.plated)
        else:
            made = _order_step(run)
        return (labelled(made, "Order and store"),)

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

    def references(self, run: Run) -> tuple[Reference, ...]:
        """Return the pool's own references, which the protocol after it prints too."""
        return POOL_REFERENCES if run.pool else ()


def pool_materials(pool: PoolPlan, pool_sheet: str, primer_sheet: str) -> tuple[Material, ...]:
    """Return what the pool route buys that a project ordering its blocks does not."""
    layout = pool.pool.layout
    roles = ", ".join(f"{count} {role}" for role, count in primer_roles(pool).items())
    return (
        Material(
            "Oligo pool",
            storage="-20 °C",
            amount=f"{pool.pool.count} oligos, every one {layout.length} nt",
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
                "for a different length wants the project's oligo length changed and the "
                "design run again.",
            ),
        ),
    )

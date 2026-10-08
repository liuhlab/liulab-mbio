"""Protocol 01: put every part in its own carrier plasmid, one part a well.

A part is a lab resource built once, so it goes in before any round runs and is kept until a
round wants it. The reaction, the plate and the plasmid are `liulab_synbio.dmx.carrier`'s: any
method keeping a part that way keeps it the same way, and this protocol only places that step
in the chain and says what the page around it shows.

This protocol only exists where the build names a carrier; a build naming none holds its parts
already and the run opens at the one after it.
"""

from liulab_mbio.checks import counted
from liulab_mbio.protocol.model import Item as Handed
from liulab_mbio.protocol.model import Material, Plate, Reference, Step
from liulab_synbio.dmx.carrier import (
    CARRIER,
    CARRIER_KIT,
    CARRIER_MARKER,
    ENZYME,
    REFERENCES,
    SeatedParts,
    carrier_step,
    materials,
)
from liulab_synbio.igga.protocols.protocol import Protocol
from liulab_synbio.igga.protocols.run import Run

#: What the page is headed and what the chain names it by.
SEATING = "Part carrier"


class Seating(Protocol):
    """Seat every part in its own carrier, which ends in one plasmid a well and no library."""

    def title(self, run: Run) -> str:
        """Return the page's own heading."""
        return SEATING

    def summary(self, run: Run) -> str:
        """Return what the one reaction comes to, which is a part collection and not a library."""
        seated = _seated(run)
        return (
            f"Put each of {counted(seated.products, 'part')} in its own {CARRIER}, one part a "
            f"well, so the collection keeps and {ENZYME} can release any of them again."
        )

    def steps(self, run: Run) -> tuple[Step, ...]:
        """Return the one step, which the carrier module writes."""
        return (carrier_step(_seated(run)),)

    def consumes(self, run: Run) -> tuple[Handed, ...]:
        """Return the parts themselves, which the bench holds before the run opens."""
        return (run.to_seat,)

    def produces(self, run: Run) -> tuple[Handed, ...]:
        """Return the plate of carrier plasmids the rounds pool their donors from."""
        return (run.carrier_plate,)

    def carried(self, run: Run) -> tuple[Material, ...]:
        """Return the carrier and its kit, which no other protocol of the run buys."""
        return materials()

    def references(self, run: Run) -> tuple[Reference, ...]:
        """Return where the kit's own reaction is read from."""
        return REFERENCES

    def plates(self, run: Run) -> tuple[Plate, ...]:
        """Return the carrier plate, so every well the step names has one."""
        return (_seated(run).plate,)

    def overview(self, run: Run) -> dict[str, str]:
        """Return the facts to check before starting, each short enough to be a card."""
        seated = _seated(run)
        return {
            "Parts": f"{seated.products}, one a well",
            "Plate": f"{seated.plate.wells} wells",
            "Carrier": CARRIER,
            "Kit": CARRIER_KIT,
            "Selection": CARRIER_MARKER,
            "Released by": ENZYME,
        }


def _seated(run: Run) -> SeatedParts:
    """Return the run's seated parts, which this protocol is only written for a run that has.

    Raises
    ------
    ValueError
        If the build names no carrier, where this protocol does not belong in the chain.
    """
    if run.seated is None:
        raise ValueError("this run names no carrier and seats no part")
    return run.seated

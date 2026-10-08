"""Protocol 01: lay the pool's amplification primers out as the plates a run pipettes from.

The primers are the pool's own, so there is nothing to plate without one. Every amount is the
build's; this invents none of them. The page is the shared primer-plate protocol, so its
steps, materials and plates are that builder's and the chain spreads nothing into it.
"""

from collections.abc import Mapping, Sequence
from dataclasses import replace

from liulab_mbio.bench.steps import primer_plate_protocol
from liulab_mbio.protocol import model
from liulab_mbio.protocol.model import Item as Handed
from liulab_mbio.protocol.model import Material, Oligo, Source, Step
from liulab_synbio.igga.protocols.protocol import Protocol
from liulab_synbio.igga.protocols.run import Run

#: What the page is headed and what the chain names it by.
PRIMER_PLATES = "Primer plates"


class PrimerPlating(Protocol):
    """Order the pool's primers and split them into the working plates a run pipettes from."""

    def title(self, run: Run) -> str:
        """Return the page's own heading."""
        return PRIMER_PLATES

    def summary(self, run: Run) -> str:
        """Return the shared builder writes this protocol's summary with the rest of its page."""
        return _plates(run).summary

    def produces(self, run: Run) -> tuple[Handed, ...]:
        """Return the stock plate and the working plate the shared builder's own steps leave."""
        return _plates(run).produces

    def page(
        self,
        run: Run,
        *,
        steps: Sequence[Step],
        materials: Sequence[Material],
        sources: Mapping[str, Source],
    ) -> model.Protocol:
        """Return the shared primer-plate protocol, which holds its own steps and materials.

        The chain spreads nothing into this one, so the arguments every other protocol is
        handed are not read here.
        """
        return replace(_plates(run), files=run.files)


def _plates(run: Run) -> model.Protocol:
    """Return the shared primer-plate protocol for this run, built once.

    Raises
    ------
    ValueError
        If the run lays out no primer plates, where this protocol is not in the chain.
    """
    if run.pool is None or run.primer_plates is None:
        raise ValueError("this run lays out no primer plates")
    asked = run.primer_plates
    return primer_plate_protocol(
        [Oligo(one.name, one.sequence, purpose=one.role) for one in run.pool.pool.primers],
        wells=asked.wells,
        copies=asked.copies,
        nanomoles=asked.nanomoles,
        stock_um=asked.stock_um,
        working_um=asked.working_um,
        working_ul=asked.working_ul,
        title=PRIMER_PLATES,
        # The pages sit in their own folder, so the sheet beside them is a step up.
        order_sheet=f"../{run.primer_sheet}" if run.primer_sheet else "",
    )


def working_plate(run: Run) -> tuple[Handed, ...]:
    """Return the plate of diluted primers this protocol leaves for the one that pipettes it.

    Empty where the run pours no plates, so the protocol after it takes the pool's primers
    from the tubes they were ordered in.
    """
    if run.pool is None or run.primer_plates is None:
        return ()
    # The stock plate stays here; only the working plate is pipetted from downstream.
    _stock, working = _plates(run).produces[:2]
    return (working,)

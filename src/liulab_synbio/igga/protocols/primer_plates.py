"""Protocol 01: lay the pool's amplification primers out as the plates a run pipettes from.

The primers are the pool's own, so there is nothing to plate without one. Every amount is the
project's; this invents none of them. The page is the shared primer-plate protocol, so its
steps, materials and plates are that builder's and the chain spreads nothing into it.
"""

from collections.abc import Mapping, Sequence

from liulab_mbio.bench.steps import primer_plate_protocol
from liulab_mbio.protocol.model import Item as Handed
from liulab_mbio.protocol.model import Material, Oligo, Protocol, Source, Step
from liulab_synbio.igga.protocols.run import Run
from liulab_synbio.igga.protocols.sitting import Sitting

#: What the page is headed and what the chain names it by.
PRIMER_PLATES = "Primer plates"


class PrimerPlating(Sitting):
    """Order the pool's primers and split them into the working plates a run pipettes from."""

    def title(self, run: Run) -> str:
        """Return the page's own heading."""
        return PRIMER_PLATES

    def summary(self, run: Run) -> str:
        """Return the shared builder writes this protocol's summary with the rest of its page."""
        return self.page(run, steps=(), materials=(), sources={}).summary

    def produces(self, run: Run) -> tuple[Handed, ...]:
        """Return the stock plate and the working plate the shared builder's own steps leave."""
        return self.page(run, steps=(), materials=(), sources={}).produces

    def page(
        self,
        run: Run,
        *,
        steps: Sequence[Step],
        materials: Sequence[Material],
        sources: Mapping[str, Source],
    ) -> Protocol:
        """Return the shared primer-plate protocol, which holds its own materials and sources.

        Raises
        ------
        ValueError
            If the run plates no primers; `included` is what says whether it does.
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
        )

"""One protocol of a library run, and what every one of them answers.

A protocol declares what it consumes and what it produces and nothing else about the run, as
`CONTEXT.md` has it. `Protocol` is what a protocol module implements, one method a thing its
page prints, so adding a fact to a protocol is one edit in one module and the chain never
dispatches on its title. `liulab_mbio.protocol.model.Protocol` is the page it writes, which is
why that one is reached through `model` here.
"""

from collections.abc import Mapping, Sequence
from dataclasses import replace

from liulab_mbio.protocol import model
from liulab_mbio.protocol.model import (
    Figure,
    Hole,
    Material,
    Plate,
    Reference,
    Source,
    Step,
    citing,
)
from liulab_mbio.protocol.model import (
    Item as Handed,
)
from liulab_synbio.igga.protocols.run import Run


class Protocol:
    """One protocol of a library run, which owns everything its own page prints.

    Every method takes the run and nothing else, so no protocol reads another's arguments. The
    defaults are the empty answers: a protocol that needs no plates, no equipment or no
    references of its own says nothing rather than being dispatched past.
    """

    #: Whether this protocol buys from the reagents the rounds share. A reagent two protocols
    #: name is bought once, so the run holds that list and the chain spreads it over these.
    round_reagents: bool = False

    def title(self, run: Run) -> str:
        """Return what the page is headed and what the chain names it by, written for the bench."""
        raise NotImplementedError

    def summary(self, run: Run) -> str:
        """One paragraph saying what this protocol of the run does."""
        raise NotImplementedError

    def steps(self, run: Run) -> tuple[Step, ...]:
        """Return the steps, in the order the bench works through them."""
        return ()

    def consumes(self, run: Run) -> tuple[Handed, ...]:
        """Return what this protocol is handed, by the names the one before it produced them under."""
        return ()

    def produces(self, run: Run) -> tuple[Handed, ...]:
        """Return what this protocol leaves for the one after it."""
        return ()

    def carried(self, run: Run) -> tuple[Material, ...]:
        """Every reagent bought for this protocol alone; the rounds' own are the run's."""
        return ()

    def equipment(self, run: Run) -> tuple[str, ...]:
        """Return the hardware this protocol needs that no reagent table covers."""
        return ()

    def plates(self, run: Run) -> tuple[Plate, ...]:
        """Return the plates this protocol fills, so every well a transfer names has one."""
        return ()

    def overview(self, run: Run) -> dict[str, str]:
        """Return the facts to check before starting, each short enough to be a card."""
        return {}

    def references(self, run: Run) -> tuple[Reference, ...]:
        """Where this protocol's own numbers are read from."""
        return ()

    def sources(self, run: Run) -> dict[str, Source]:
        """Every document this protocol could cite; `citing` drops the ones it did not."""
        return {}

    def holes(self, run: Run) -> tuple[Hole, ...]:
        """Return what this protocol leaves to the reader because no source settles it."""
        return ()

    def page(
        self,
        run: Run,
        *,
        steps: Sequence[Step],
        materials: Sequence[Material],
        sources: Mapping[str, Source],
    ) -> model.Protocol:
        """Return the protocol as the page prints it.

        `steps` is what `steps` returned, built once because the chain reads them to spread the
        reagents. `materials` is this protocol's share of those reagents, and `sources` is every
        document the run could cite, so a protocol citing another's document still resolves it.
        """
        return citing(
            model.Protocol(
                self.title(run),
                summary=self.summary(run),
                overview=self.overview(run),
                consumes=self.consumes(run),
                produces=self.produces(run),
                materials=tuple(materials),
                equipment=self.equipment(run),
                plates=self.plates(run),
                steps=tuple(steps),
                references=self.references(run),
                sources=dict(sources),
                holes=self.holes(run),
            )
        )


def labelled(step: Step, section: str) -> Step:
    """Return `step` under the stage of its protocol it belongs to."""
    return replace(step, section=section)


def figured(step: Step, figure: Figure | None) -> Step:
    """Return `step` showing `figure`, or unchanged where there is no record to draw."""
    return step if figure is None else replace(step, figures=(figure,))

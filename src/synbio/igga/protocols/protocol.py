"""One protocol of a library run, and what every one of them answers.

A protocol declares what it consumes and what it produces and nothing else about the run, as
`CONTEXT.md` has it. `Protocol` is what a protocol module implements, one method a thing its
page prints, so adding a fact to a protocol is one edit in one module and the chain never
dispatches on its title. `mbio.protocol.model.Protocol` is the page it writes, which is
why that one is reached through `model` here.
"""

from collections.abc import Mapping, Sequence

from mbio.protocol import model
from mbio.protocol.model import (
    Hole,
    Material,
    Plate,
    Source,
    Step,
    citing,
)
from mbio.protocol.model import (
    Item as Handed,
)
from synbio.igga.protocols.run import Run

#: The reads that close the run, which two protocols both label.
READ_BACK_SECTION = "Read the library back"


class Protocol:
    """One protocol of a library run, which owns everything its own page prints.

    Every method takes the run and nothing else, so no protocol reads another's arguments. The
    defaults are the empty answers: a protocol that needs no plates, no equipment or no
    sources of its own says nothing rather than being dispatched past.
    """

    def title(self, run: Run) -> str:
        """Return what the page is headed and what the chain names it by, written for the bench."""
        raise NotImplementedError

    def summary(self, run: Run) -> str:
        """One paragraph saying what this protocol of the run does."""
        raise NotImplementedError

    def choice(self, run: Run) -> str:
        """Return the job this protocol is one way of doing, where the run offers several ways.

        Empty is what a protocol that is simply a step of the chain says, which is most of them.
        """
        return ""

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

    def shares(self, run: Run) -> tuple[str, ...]:
        """Name the run's shared reagents this protocol's steps use but never spell.

        A reagent two protocols need is bought once, so the run holds it and each page claims
        its own. A reagent a step spells reaches the page without being claimed; this is for
        the rest -- a buffer, a plate, a consumable no sentence mentions. Each name is matched
        inside the reagent's own, so ``"part list"`` claims every one of them.
        """
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
                choice=self.choice(run),
                consumes=self.consumes(run),
                produces=self.produces(run),
                materials=tuple(materials),
                files=run.files,
                equipment=self.equipment(run),
                plates=self.plates(run),
                steps=tuple(steps),
                sources=dict(sources),
                holes=self.holes(run),
            )
        )

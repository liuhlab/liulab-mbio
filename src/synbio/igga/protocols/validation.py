"""Protocol 05: take the archived designs to clonal wells and read each one back.

`synbio.dmx` is a method of its own and takes any cargo; what a library run builds is one
cargo among others. So nothing of the read-back is written here: this page is only its title,
the figure of what a well holds, and its place in a library run's chain. The steps, the
materials and the plates are `dmx`'s, for the route this page is written for.

A build stating no fragment-count floor reads nothing back and this protocol is not in its
chain. A build naming more than one route writes one of these pages per route, and they are the
ways of one job: the bench does one of them.
"""

from mbio.protocol.figures import SOURCE_KEY as FIGURE_SOURCE_KEY
from mbio.protocol.model import Citation, Figure, Material, Plate, Reference, Source, Step
from mbio.protocol.model import Item as Handed
from mbio.sequence import Segment, SequenceRecord
from synbio import dmx
from synbio.igga.figures import STUFFER_MARGIN
from synbio.igga.protocols.protocol import Protocol
from synbio.igga.protocols.run import Run

#: What the page is headed, before the route that marks its wells is named after it.
VALIDATION = "Cargo validation"

#: What a vector's own records call the piece its cassette's enzyme cuts out.
STUFFER_FEATURE = "internal stuffer"

#: Where in ``docs/research/figure-sources.md`` the read-back figure's equivalent stands.
VALIDATION_CITATION = Citation(FIGURE_SOURCE_KEY, "section 3, DMX")


def validation_title(one: dmx.Validation) -> str:
    """Name the read-back protocol for the route that marks its wells."""
    return f"{VALIDATION}: {one.route.name}"


class ReadBack(Protocol):
    """Array, pick, mark and call every design the build's floor asks to read back.

    One of these is written for each route the build named, so a run offering two routes writes
    two pages and the bench does one of them.
    """

    def __init__(self, validation: dmx.Validation) -> None:
        """Hold the read-back this page is written for, which is one route's."""
        self.validation = validation

    def title(self, run: Run) -> str:
        """Return the page's heading, which names the route that marks its wells."""
        return validation_title(self.validation)

    def choice(self, run: Run) -> str:
        """Name the job this page is one way of doing, where the run offers more than one way."""
        return dmx.READ_BACK if run.read_back_is_a_choice else ""

    def summary(self, run: Run) -> str:
        """Return what reading the designs back comes to."""
        return (
            "Take the archived designs to clonal wells, mark each well so sequencing says "
            "which well it came from, and call a pass or a fail per well."
        )

    def steps(self, run: Run) -> tuple[Step, ...]:
        """Return the route's own steps, each under the section `dmx` puts it in."""
        return dmx.validation_steps(self.validation, marking=_marking_figure(run, self.validation))

    def consumes(self, run: Run) -> tuple[Handed, ...]:
        """Return the archive plate, and the lab stock this page's own route marks with."""
        stock = dmx.marking_stock(self.validation)
        return (run.archive, *((stock,) if stock else ()))

    def produces(self, run: Run) -> tuple[Handed, ...]:
        """One clone a well, and a call for each of them."""
        return (run.picked, run.calls)

    def carried(self, run: Run) -> tuple[Material, ...]:
        """Return what the route buys, which no other protocol of the run does."""
        return dmx.validation_materials(self.validation)

    def equipment(self, run: Run) -> tuple[str, ...]:
        """Return the hardware the route needs that no reagent table covers."""
        return dmx.validation_equipment(self.validation)

    def plates(self, run: Run) -> tuple[Plate, ...]:
        """Return the plates the read-back fills, so every well a transfer names has one."""
        return self.validation.plates

    def references(self, run: Run) -> tuple[Reference, ...]:
        """Where the route's own numbers are read from."""
        return dmx.REFERENCES

    def sources(self, run: Run) -> dict[str, Source]:
        """Every document the route could cite."""
        return dmx.SOURCES


def _marking_figure(run: Run, validation: dmx.Validation) -> Figure | None:
    """Return what the route's own marking step works on: what a picked well holds.

    The cassette its design sits in, with the stuffer lit. Drawn as a map and not at base
    level on either route: the kit's barcodes chain on through its own chemistry and the index
    pair's sites are the lab's, so neither route's cut is a cut this record spells.
    """
    found = run.block_vector
    if found is None:
        return None
    vector, name, path = found
    marks = (
        f"The kit chains its {dmx.GROUPS} barcodes on, {dmx.CHAIN[0]} through "
        f"{dmx.CHAIN[-1]}, reading on the strand the cargo reads on."
        if validation.route is dmx.ROUTE_LIGATION
        else "The pair reads across it, and the band is that stretch plus the two marks the "
        "well's address names."
    )
    return Figure(
        (path,),
        f"The cassette {name} carries, drawn empty: a picked well holds one with its design in "
        f"place of the stuffer. {marks}",
        span=_cassette(vector.stuffer, len(vector.record)),
        linear=True,
        highlight=_drawn_as(vector.record, STUFFER_FEATURE),
        citation=VALIDATION_CITATION,
    )


def _cassette(stuffer: Segment, length: int) -> tuple[int, int]:
    """Return the stretch a map of the cassette draws: the stuffer and its neighbourhood.

    The margin is the stuffer's own length, so the cassette fills a known share of the picture
    whatever a vector spells there. A linear record stops at its ends.
    """
    margin = (stuffer.end - stuffer.start) * STUFFER_MARGIN
    return max(0, stuffer.start - margin), min(length, stuffer.end + margin)


def _drawn_as(record: SequenceRecord, name: str) -> tuple[str, ...]:
    """Return `name` where the record draws a feature under it, and nothing where it does not.

    A vector that came in already carrying its cassette names its own features, and a highlight
    lights what a record holds rather than what this module would like it to.
    """
    return (name,) if any(one.name == name for one in record.features) else ()

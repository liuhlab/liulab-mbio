"""Protocol 04: take the archived designs to clonal wells and read each one back.

The route the project chose is what marks a well, and `liulab_synbio.dmx` owns both routes:
this protocol is where a library run calls for one. A project stating no fragment-count floor
reads nothing back and this protocol is not in its chain.
"""

from liulab_mbio.protocol.figures import SOURCE_KEY as FIGURE_SOURCE_KEY
from liulab_mbio.protocol.model import Citation, Figure, Material, Plate, Reference, Source, Step
from liulab_mbio.protocol.model import Item as Handed
from liulab_mbio.sequence import Segment, SequenceRecord
from liulab_synbio import dmx
from liulab_synbio.igga.protocols.run import Run
from liulab_synbio.igga.protocols.sitting import Sitting, labelled

#: What the page is headed, before the route that marks its wells is named after it.
VALIDATION = "Cargo validation"

#: What a vector's own records call the piece its cassette's enzyme cuts out.
STUFFER_FEATURE = "internal stuffer"

#: How many times its own length a map of the cassette draws either side of the stuffer, so the
#: flanks it sits between are in the picture. A margin for the drawing, measured from nothing.
STUFFER_MARGIN = 3

#: Where in ``docs/research/figure-sources.md`` the read-back figure's equivalent stands.
VALIDATION_CITATION = Citation(FIGURE_SOURCE_KEY, "section 3, DMX")


def validation_title(one: dmx.Validation) -> str:
    """Name the read-back protocol for the route that marks its wells."""
    return f"{VALIDATION}: {one.route.name}"


class ReadBack(Sitting):
    """Array, pick, mark and call every design the project's floor asks to read back."""

    def title(self, run: Run) -> str:
        """Return the page's heading, which names the route that marks its wells."""
        return validation_title(_validation_of(run))

    def summary(self, run: Run) -> str:
        """Return what reading the designs back comes to."""
        return (
            "Take the archived designs to clonal wells, mark each well so sequencing says "
            "which well it came from, and call a pass or a fail per well."
        )

    def steps(self, run: Run) -> tuple[Step, ...]:
        """Return the route's own steps, under the four stages a reader works through them in."""
        made = dmx.validation_steps(_validation_of(run), marking=_marking_figure(run))
        return (
            *(labelled(one, "Array and pick") for one in made[:2]),
            *(labelled(one, "Mark every well") for one in made[2:-1]),
            labelled(made[-1], "Call the wells"),
        )

    def consumes(self, run: Run) -> tuple[Handed, ...]:
        """Return the archive plate, and the lab stock the chosen route marks with."""
        stock = run.marking_stock
        return (run.archive, *((stock,) if stock else ()))

    def produces(self, run: Run) -> tuple[Handed, ...]:
        """One clone a well, and a call for each of them."""
        return (run.picked, run.calls)

    def carried(self, run: Run) -> tuple[Material, ...]:
        """Return what the route buys, which no other protocol of the run does."""
        return dmx.validation_materials(_validation_of(run))

    def equipment(self, run: Run) -> tuple[str, ...]:
        """Return the hardware the route needs that no reagent table covers."""
        return dmx.validation_equipment(_validation_of(run))

    def plates(self, run: Run) -> tuple[Plate, ...]:
        """Return the plates the read-back fills, so every well a transfer names has one."""
        one = _validation_of(run)
        return (*one.picked, *(one.index if one.route is dmx.ROUTE_INDEX_PCR else ()))

    def references(self, run: Run) -> tuple[Reference, ...]:
        """Where the route's own numbers are read from."""
        return dmx.REFERENCES

    def sources(self, run: Run) -> dict[str, Source]:
        """Every document the route could cite."""
        return dmx.SOURCES


def _validation_of(run: Run) -> dmx.Validation:
    """Return the run's read-back, which this protocol is only written for a run that has one.

    Raises
    ------
    ValueError
        If the project states no fragment-count floor, where this protocol is not in the chain.
    """
    if run.validation is None:
        raise ValueError("this run states no validation floor and reads no design back")
    return run.validation


def _marking_figure(run: Run) -> Figure | None:
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
        if _validation_of(run).route is dmx.ROUTE_LIGATION
        else "The pair reads across it, and the band is that stretch plus the two marks the "
        "well's address names."
    )
    return Figure(
        (path,),
        f"What a picked well holds: its design in the cassette {name} carries. {marks}",
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

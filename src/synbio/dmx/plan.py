"""`plan_dmx`: the way in. One build file in, one protocol page out.

A run reads designs the lab already holds, so there is no DNA to design and no record to write:
what a plan is here is the protocol someone works through, as a run of one through
`mbio.cloning.plan.as_project`. A run of one is written as one page, so what the build says
about the shape of the run rides on that page rather than on an index of one.

The run is one protocol in three sections, as `synbio.igga` has rendered these same steps
since it first chained them.
"""

import os
from dataclasses import dataclass, replace
from pathlib import Path

from mbio.checks import counted
from mbio.cloning.plan import as_project
from mbio.protocol.model import Item, Project, Protocol, Topic, citing
from mbio.protocol.render import write_run_files
from synbio.dmx.build import Build, read_build, read_designs
from synbio.dmx.method import (
    PICKED_PLATE,
    REFERENCES,
    SOURCES,
    WELL_CALLS,
    Design,
    Validation,
    marking_stock,
    selected_on,
    validation,
    validation_equipment,
    validation_materials,
)
from synbio.dmx.steps import pooled_plates, validation_steps


@dataclass(frozen=True, slots=True)
class Files:
    """The files a read-back plan writes.

    Parameters
    ----------
    protocol_data
        The run as one JSON file, which ``protocol render`` turns back into the page.
    protocol
        The interactive bench protocol, as one self-contained HTML page rendered from
        `protocol_data`.
    """

    protocol_data: Path
    protocol: Path

    @property
    def paths(self) -> tuple[Path, ...]:
        """Both, in the order they were written."""
        return (self.protocol_data, self.protocol)


@dataclass(frozen=True, slots=True)
class ReadBackPlan:
    """One planned read-back: which designs are read, on which route, and what that takes.

    Parameters
    ----------
    build
        What this run chose: its designs, its archive, its route and its floor.
    designs
        Every design the sheet named, before the floor chose among them.
    validation
        What reading the chosen designs back takes: the wells, the plates and the depth.
    """

    build: Build
    designs: tuple[Design, ...]
    validation: Validation

    @property
    def archive(self) -> Item:
        """The plate of clonal stock each design is spotted from, which this run only reads."""
        return Item(
            self.build.archive,
            f"one clonal stock a design, {len(self.designs)} of them, sealed and frozen",
            storage="-80 °C glycerol stock",
        )

    def protocol(self) -> Protocol:
        """Return the page someone works through, its three sections labelled."""
        one = self.validation
        stock = marking_stock(one)
        return citing(
            Protocol(
                f"Design read-back: {one.route.name}",
                summary=(
                    f"Take {counted(len(one.designs), 'archived design')} to {one.wells} clonal "
                    "wells, "
                    "mark each well so sequencing says which well it came from, and call a pass "
                    "or a fail per well."
                ),
                overview=self._overview(),
                highlights=self._highlights(),
                consumes=(self.archive, *((stock,) if stock else ())),
                produces=(PICKED_PLATE, WELL_CALLS),
                materials=validation_materials(one),
                equipment=validation_equipment(one),
                plates=one.plates,
                steps=validation_steps(one),
                references=REFERENCES,
                sources=dict(SOURCES),
            )
        )

    def chain(self) -> Project:
        """Return the run as a project of one, carrying why it is shaped as it is."""
        return replace(
            as_project(self.protocol()),
            title=self.build.name,
            background=self._background(),
        )

    def write(self, directory: str | os.PathLike[str]) -> Files:
        """Write the run's data and its page into `directory`.

        The directory is made when it is not there, the files are named by
        `mbio.protocol.render`, and a second run over the same inputs writes the same
        bytes. A run of one protocol is one page, so nothing here is a folder.
        """
        written = write_run_files(self.chain(), Path(directory))
        return Files(written.data, written.page)

    def _overview(self) -> dict[str, str]:
        """Return the facts a reader scans before starting."""
        one = self.validation
        return {
            "Designs": f"{len(one.designs)} of {len(self.designs)} on the sheet",
            "Route": one.route.name,
            "Wells": f"{one.wells} over {counted(len(one.picked), 'picked plate')}",
            "Selection": selected_on(one.selection),
        }

    def _highlights(self) -> tuple[str, ...]:
        """Return what those facts mean, a sentence each."""
        return (
            "A design is read back as a name and a count of the pieces it was built from, never "
            "as bases, so any clonal stock this lab holds can be read back here.",
            f"The bench is sized from the {counted(len(self.validation.designs), 'design')} read "
            "and not "
            "from the sheet, so a design the floor leaves out costs no well, no plate and no "
            "reagent.",
        )

    def _background(self) -> tuple[Topic, ...]:
        """Return why the run is shaped as it is, a topic at a time."""
        one = self.validation
        floor = (
            "every design on the sheet, which a floor of zero does"
            if one.floor == 0
            else f"every design built from {counted(one.floor, 'fragment')} or more"
        )
        pooled = counted(len(pooled_plates(one)), "plate")
        apart = counted(one.route.plate_axis, "plate")
        return (
            Topic(
                "Which designs are read back",
                (
                    f"This run reads {floor}: {len(one.designs)} of {len(self.designs)}. The "
                    "rest stay polyclonal and are never read one design at a time.",
                    f"Every design is spotted from {self.build.archive}, which this run reads a "
                    "copy of and never changes. A design whose spot grows nothing has no clone "
                    "to read, and is re-transformed before anyone re-synthesises it.",
                ),
            ),
            Topic(
                "How a well is told from the rest",
                (
                    f"The {one.route.name} route marks each well: {one.route.marking}. The other "
                    "route marks the same wells a different way, and the picking, the pass rule "
                    "and the reformat are the same either way.",
                    "A well's marks are worked out from where the well is rather than looked up, "
                    "so a demultiplexer can check an address instead of trusting a file. One "
                    f"axis of the address is the plate: this run pools {pooled}, and the marks "
                    f"tell {apart} apart on one flow cell.",
                    f"A well is called above {one.route.wanted_reads} reads. The other route's "
                    "floor is not a stricter or a looser setting of the same scale, so neither "
                    "number carries over to the other route.",
                ),
            ),
        )


def plan_dmx(build: Build | str | os.PathLike[str]) -> ReadBackPlan:
    """Plan the read-back `build` asks for: which designs, on which route, and what that takes.

    Parameters
    ----------
    build
        What this run chooses, or a path to the JSON holding it;
        `synbio.dmx.build.read_build` reads one. It names the designs sheet by path.

    Returns
    -------
    ReadBackPlan
        The plan, whose `write` puts the run's data and its page in one directory.

    Raises
    ------
    ValueError
        If the build or its designs sheet cannot be read, or its floor reads no design back.

    Examples
    --------
    >>> plan_dmx("build.json").validation.wells  # doctest: +SKIP
    8
    """
    one = build if isinstance(build, Build) else read_build(build)
    designs = read_designs(one.designs)
    sized = validation(
        one.marking_route,
        designs,
        one.validate_from,
        selection=one.selection,
        index_plate=one.index_plate,
    )
    if sized is None:
        raise ValueError(
            f"validate_from is {one.validate_from}, and no design on {os.fspath(one.designs)} is "
            f"built from that many fragments: {counted(len(designs), 'design')} read, "
            "none of them back"
        )
    return ReadBackPlan(one, designs, sized)

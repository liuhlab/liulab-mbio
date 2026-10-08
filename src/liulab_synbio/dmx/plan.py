"""`plan_dmx`: the way in. One build file in, one protocol folder out.

A run reads designs the lab already holds, so there is no DNA to design and no record to write:
what a plan is here is the protocol someone works through, as a run of one through
`liulab_mbio.cloning.plan.as_project`.

The run is one protocol in three sections, as `liulab_synbio.igga` has rendered these same steps
since it first chained them.
"""

import os
from dataclasses import dataclass, replace
from pathlib import Path

from liulab_mbio.cloning.plan import as_project
from liulab_mbio.protocol.model import Item, Project, Protocol, Topic, citing
from liulab_mbio.protocol.render import write_project_files
from liulab_synbio.dmx.build import Build, read_build, read_designs
from liulab_synbio.dmx.method import (
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
from liulab_synbio.dmx.steps import validation_steps


@dataclass(frozen=True, slots=True)
class Files:
    """The files a read-back plan writes.

    Parameters
    ----------
    protocol_data
        The run as one JSON file, which ``protocol render`` turns back into pages.
    protocol
        Every page rendered from `protocol_data`: the index first, then the protocol, then the
        two pages the run shares.
    """

    protocol_data: Path
    protocol: tuple[Path, ...]

    @property
    def paths(self) -> tuple[Path, ...]:
        """Every file written, the first written first."""
        return (self.protocol_data, *self.protocol)


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
                    f"Take {len(one.designs)} archived design(s) to {one.wells} clonal wells, "
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
        """Write the run's data and its pages into `directory`, which becomes the folder.

        The directory is made when it is not there, the files are named by
        `liulab_mbio.protocol.render`, and a second run over the same inputs writes the same
        bytes.
        """
        written = write_project_files(self.chain(), Path(directory))
        return Files(
            written.data,
            (written.index, *written.protocols, written.reagents, written.references),
        )

    def _overview(self) -> dict[str, str]:
        """Return the facts a reader scans before starting."""
        one = self.validation
        return {
            "Designs": f"{len(one.designs)} of {len(self.designs)} on the sheet",
            "Route": one.route.name,
            "Wells": f"{one.wells} over {len(one.picked)} picked plate(s)",
            "Selection": selected_on(one.selection),
        }

    def _highlights(self) -> tuple[str, ...]:
        """Return what those facts mean, a sentence each."""
        return (
            "A design is read back as a name and a count of the pieces it was built from, never "
            "as bases, so any clonal stock this lab holds can be read back here.",
            f"The bench is sized from the {len(self.validation.designs)} design(s) read and not "
            "from the sheet, so a design the floor leaves out costs no well, no plate and no "
            "reagent.",
        )

    def _background(self) -> tuple[Topic, ...]:
        """Return why the run is shaped as it is, a topic at a time."""
        one = self.validation
        floor = (
            "every design on the sheet, which a floor of zero does"
            if one.floor == 0
            else f"every design built from {one.floor} fragment(s) or more"
        )
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
                    "so two plates on one flow cell are told apart by construction and a "
                    "demultiplexer can check an address instead of trusting a file.",
                    f"A well is called above {one.route.wanted_reads} reads, which is the depth "
                    f"this route was measured at and no other.",
                ),
            ),
        )


def plan_dmx(build: Build | str | os.PathLike[str]) -> ReadBackPlan:
    """Plan the read-back `build` asks for: which designs, on which route, and what that takes.

    Parameters
    ----------
    build
        What this run chooses, or a path to the JSON holding it;
        `liulab_synbio.dmx.build.read_build` reads one. It names the designs sheet by path.

    Returns
    -------
    ReadBackPlan
        The plan, whose `write` puts the run's data and its pages in one directory.

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
            f"built from that many fragments: {len(designs)} design(s) read, none of them back"
        )
    return ReadBackPlan(one, designs, sized)

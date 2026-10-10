"""What every cloning plan is: the files it writes, its status, and the records it was made from.

A method's own plan module holds its design, its `Plan` and the file set it writes. What every
one of them repeats is here: the shape a plan and its written files take, the names of the
records and sheets any plan writes, a protocol read as the run of one it is, the worst-of rule a
plan's status follows, where a method puts its inserts and which way round they go, which oligos
the product it writes is drawn with, how its oligos are judged as one, and taking a record the
caller already read. What a run writes to disk is `mbio.protocol.render.write_run_files`, which
every pipeline shares.

The library pipeline is not a cloning method -- `docs/adr/0004-library-rounds.md` says why --
but its plan writes the same files, so it uses this module too.
"""

import os
import typing  # Spelled out: `Protocol` here is the bench protocol imported below.
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import replace
from pathlib import Path
from typing import Literal

from mbio.bench.validation import ColonyCheck, SangerRead
from mbio.checks import Check, Status, worst, worst_of
from mbio.io import read_record
from mbio.primers.evaluation import PrimerReport
from mbio.primers.thresholds import reading
from mbio.protocol.model import Project, Protocol, slug
from mbio.sequence import Feature, SequenceRecord, across_the_origin

#: What a plan calls the product it writes: the annotated plasmid the design makes.
PRODUCT_FILE = "product.dna"

#: What a plan calls the sheet it orders its oligos from.
PRIMER_FILE = "primers.tsv"

#: What a plan calls the record of one amplicon a step makes, after the part it makes.
AMPLICON_FILE = "{name}-amplicon.dna"

#: The stretch where the pieces to be joined are made, named by every method that amplifies one.
MAKE_SECTION = "Make the fragments"

#: The stretch where the pieces become one molecule.
JOIN_SECTION = "Assemble"

#: The stretch where the joined molecule goes into cells and the right colony is found.
SCREEN_SECTION = "Transform and screen"

#: Where a method puts its inserts: a feature name, a span, or `None` to look for `MCS_FEATURE`.
type Site = str | tuple[int, int] | None

#: The feature a vector names its cloning site with, looked for when the caller names none.
MCS_FEATURE = "MCS"

#: Which way round an insert goes into the vector.
type Orientation = Literal["forward", "reverse"]


class Written(typing.Protocol):
    """The files one plan wrote, as anything reporting them reads that record.

    A pipeline's own `Files` names each file it writes as a field. `paths` is the order those
    files are reported in, which is the order they were written.
    """

    @property
    def paths(self) -> tuple[Path, ...]:
        """Every file the plan wrote, the first written first."""
        ...


class Planned(typing.Protocol):
    """What every pipeline's plan promises: it writes its files into one directory."""

    def write(self, directory: str | os.PathLike[str], /) -> Written:
        """Write every file into `directory`, and say where each one went."""
        ...


def amplicon_files(names: Iterable[str]) -> tuple[str, ...]:
    """Return what a plan calls the amplicon of each part named, in order and no two alike.

    Examples
    --------
    >>> amplicon_files(["pUC19 backbone", "GFP", "GFP"])
    ('puc19-backbone-amplicon.dna', 'gfp-amplicon.dna', 'gfp-2-amplicon.dna')
    """
    taken: set[str] = set()
    files = []
    for name in names:
        stem = slug(name) or "part"
        chosen, n = stem, 1
        while chosen in taken:
            n += 1
            chosen = f"{stem}-{n}"
        taken.add(chosen)
        files.append(AMPLICON_FILE.format(name=chosen))
    return tuple(files)


def ordered_from_sheet(protocol: Protocol) -> Protocol:
    """Return `protocol` naming the primer sheet and the product its plan writes beside its page.

    Every cloning plan writes `PRIMER_FILE` and `PRODUCT_FILE` into the directory it writes the
    page into, so the page links the sheet rather than being the thing a supplier is sent, and
    any text naming either file links it.
    """
    return replace(protocol, order_sheet=PRIMER_FILE, files=(PRIMER_FILE, PRODUCT_FILE))


def as_project(protocol: Protocol) -> Project:
    """Return `protocol` as a run of one, which is what every one-protocol method plans.

    What it consumes is what the bench already holds, so those are the project's inputs and the
    chain resolves. Whatever reads a run then reads one shape, and a run of several protocols is
    no special case of it. What that run writes still follows its length:
    `mbio.protocol.render.write_run_files` turns a run of one back into one page.
    """
    return Project(
        protocol.title,
        summary=protocol.summary,
        inputs=protocol.consumes,
        protocols=(protocol,),
    )


def status(checks: Iterable[Check]) -> Status:
    """Return a plan's status: the worst verdict any of its checks carries."""
    return worst_of(checks)


def as_record(value: SequenceRecord | str | os.PathLike[str]) -> SequenceRecord:
    """Read a record from a path, or take one already read."""
    return value if isinstance(value, SequenceRecord) else read_record(value)


def insertion_span(vector: SequenceRecord, site: Site) -> tuple[int, int]:
    """Return the vector bases a method's inserts replace, named as `Site` allows.

    Every method reads its insertion site the same way, so a user moving between them learns
    nothing new. The span is the model's: 0-based, half-open, and inside the vector.

    Raises
    ------
    ValueError
        If no site is named and the vector annotates no `MCS_FEATURE` feature, if it annotates
        no feature of the name given, if that feature runs across the origin of a circular
        vector or lies at both ends of a linear one, or if the span does not lie inside the
        vector.
    """
    if isinstance(site, tuple):
        start, end = site
    else:
        if site is None:
            found = _named(vector, MCS_FEATURE)
            if found is None:
                raise ValueError(
                    f"name the insertion site: this vector annotates no {MCS_FEATURE!r} "
                    "feature. Pass a feature name or a (start, end) span"
                )
        else:
            found = _named(vector, site)
            if found is None:
                raise ValueError(f"this vector annotates no feature called {site!r}")
        if across_the_origin(found, len(vector)):
            what = vector.name or "the vector"
            if vector.topology == "linear":
                raise ValueError(
                    f"{what} is linear, and the insertion site {found.name!r} lies in two pieces "
                    "at its two ends; pass a (start, end) span that lies between them"
                )
            raise ValueError(
                f"the insertion site {found.name!r} runs across the origin of {what}, and a "
                "plan cannot replace bases there; move the vector's origin away from the site, "
                "or pass a (start, end) span that stays on one side of the origin"
            )
        start = min(one.start for one in found.segments)
        end = max(one.end for one in found.segments)
    if not 0 <= start < end <= len(vector):
        raise ValueError(
            f"the insertion site {start}-{end} does not lie inside {len(vector)} bases"
        )
    return start, end


def orientations(
    orientation: Orientation | Sequence[Orientation], count: int
) -> tuple[Orientation, ...]:
    """Spread one orientation over every insert, or take the one given for each.

    Raises
    ------
    ValueError
        If a value is neither ``forward`` nor ``reverse``, or there is not one per insert.

    Examples
    --------
    >>> orientations("reverse", 2)
    ('reverse', 'reverse')
    """
    given = (orientation,) * count if isinstance(orientation, str) else tuple(orientation)
    if len(given) != count:
        raise ValueError(f"orientation has {len(given)} values for {count} insert(s)")
    for one in given:
        if one not in ("forward", "reverse"):
            raise ValueError(f"orientation is 'forward' or 'reverse', got {one!r}")
    return given


def _named(record: SequenceRecord, name: str) -> Feature | None:
    """Return the first feature of that name, whatever its case."""
    return next(
        (feature for feature in record.features if feature.name.lower() == name.lower()), None
    )


def annotated(
    product: SequenceRecord, colony: ColonyCheck, reads: Iterable[SangerRead]
) -> SequenceRecord:
    """Return `product` drawing the oligos a plan designed on it, each where it anneals.

    Every cloning plan annotates every primer it designed, so a map of one method's product
    reads like a map of another's. A part's amplification primers reach the product with the
    part; the colony PCR and sequencing primers are designed on the product itself, so each
    already carries the site it binds and what it was missing is a place on the record the user
    opens.
    """
    designed = (*colony.primers, *(read.primer for read in reads))
    return replace(product, primers=(*product.primers, *designed))


def primer_check(reports: Sequence[PrimerReport]) -> Check:
    """Return one verdict over every oligo a plan designed, the worst of theirs.

    A count alone cannot be acted on, so the detail names each kind of check that fired with
    the rows it covers, and names apart what no sourced threshold judged.
    """
    warned = sum(1 for report in reports if report.status == "warn")
    failed = sum(1 for report in reports if report.status == "fail")
    counted = Counter(
        reading(check).label
        for report in reports
        for check in report.checks
        if check.status not in (None, "pass")
    )
    unjudged = dict.fromkeys(
        reading(check).label
        for report in reports
        for check in report.checks
        if check.status is None
    )
    said = f"{len(reports)} designed, {warned} with a warning, {failed} failing"
    if counted:
        kinds = ", ".join(
            f"{label} on {rows}"
            for label, rows in sorted(counted.items(), key=lambda one: (-one[1], one[0]))
        )
        said += f": {kinds}"
    detail = f"{said}; not judged: {', '.join(unjudged)}" if unjudged else said
    return Check("primers", worst(report.status for report in reports), len(reports), detail)

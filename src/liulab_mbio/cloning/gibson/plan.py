"""One Gibson experiment, planned from a vector and the inserts that go round it.

`plan_gibson` runs the whole design: open the vector by PCR across the span the inserts
replace, choose the overlap at each junction, design the primers that carry it, simulate the
product, work out what the assembly reaction takes, and design the colony PCR and the
sequencing that say whether the clone is the one the design asked for. `Plan.write` puts four
files in one directory -- the annotated product, an oligo order sheet, the protocol as JSON
data, and the interactive HTML page rendered from that data.

Every number the protocol prints is computed here or by the modules this one calls, and every
supplier's number behind them is `liulab_mbio.cloning.gibson.bench`, through
``docs/research/gibson-assembly.md``. What the protocol says about the phenotype is
`liulab_mbio.bench.phenotype`, read off the product's own features.
"""

import dataclasses
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.oligos import OrderedOligo, primer_sheet
from liulab_mbio.bench.phenotype import Phenotype, read_phenotype
from liulab_mbio.bench.validation import (
    ColonyCheck,
    SangerRead,
    colony_pcr_check,
    sanger_primers,
)
from liulab_mbio.checks import Check, Status
from liulab_mbio.cloning.gibson.assembly import (
    Assembly,
    Part,
    amplify,
    assemble,
    given_vector,
    open_vector,
    stitch,
)
from liulab_mbio.cloning.gibson.bench import (
    NEBUILDER_HIFI,
    AssemblyProduct,
    OverlapRule,
    assembly_amounts,
    assembly_dna_check,
    fragment_check,
)
from liulab_mbio.cloning.gibson.design import (
    BRIDGE_OLIGO_PMOL,
    NO_VERDICT,
    OLIGO_PURITY,
    STITCH_OLIGO_NM,
    bridging_oligo,
    overlap_after,
    overlap_before,
    overlap_checks,
    requires_oligos,
    stitch_checks,
    stitch_oligos,
)
from liulab_mbio.cloning.gibson.oligos import DesignedOligo
from liulab_mbio.cloning.gibson.steps import DEFAULT_HOST
from liulab_mbio.cloning.gibson.steps import protocol as protocol_for
from liulab_mbio.cloning.plan import (
    PRIMER_FILE,
    PRODUCT_FILE,
    Orientation,
    Site,
    as_record,
    insertion_span,
    orientations,
    primer_check,
    status,
    write_protocol_files,
)
from liulab_mbio.edits import flipped
from liulab_mbio.primers.evaluation import PrimerReport, evaluate_primer
from liulab_mbio.primers.polymerase import ONETAQ, Q5, Polymerase
from liulab_mbio.primers.thresholds import THRESHOLDS_FOR, PrimerRole, Thresholds
from liulab_mbio.protocol.model import Protocol
from liulab_mbio.sequence import Primer, SequenceRecord
from liulab_mbio.snapgene import write_dna

#: How a part reaches the reaction: amplified by PCR, or stitched out of overlapping oligos.
type Route = Literal["amplify", "stitch"]


@dataclass(frozen=True, slots=True)
class Files:
    """The four files a plan writes.

    Parameters
    ----------
    product
        The annotated product, as a SnapGene ``.dna`` file.
    primers
        Every designed oligo, as a tab-separated sheet to order from.
    protocol_data
        The bench protocol as JSON, which ``protocol render`` turns back into a page.
    protocol
        The interactive bench protocol, as one self-contained HTML page rendered from
        `protocol_data`.
    """

    product: Path
    primers: Path
    protocol_data: Path
    protocol: Path

    @property
    def paths(self) -> tuple[Path, ...]:
        """The four, in the order they were written."""
        return (self.product, self.primers, self.protocol_data, self.protocol)


@dataclass(frozen=True, slots=True)
class Plan:
    """One planned Gibson experiment.

    Parameters
    ----------
    vector, inserts
        The records the plan was made from, each insert on the strand that goes in.
    span
        The vector bases the inserts replace. Empty, at the vector's end, where the vector was
        handed in linear and nothing was taken out of it.
    product
        The assembly product on the bench, which sets the overlap rule, the reaction and the
        incubation.
    linearised_vector
        The part the vector is opened into.
    insert_parts
        The part each insert is amplified into, in insert order.
    assembly
        The simulated plasmid, the parts and the junctions.
    colony
        The colony PCR that reads every junction and tells a correct clone from an empty vector.
        No insert can go in the other way round, so there is no lane for one.
    reads
        A sequencing primer reading in from outside the first junction and the last.
    amounts
        What to put in the assembly reaction, vector first.
    phenotype
        What the plasmid says about itself.
    designed_oligos
        Every designed oligo and what it is for, in order: each part's forward and reverse
        primer, the colony PCR, the reads.
    ordered_oligos
        Every oligo the plan orders that primes nothing: a stitched part's, then any bridging
        oligo. Nothing measures one, so its row carries no verdict.
    host, polymerase
        The choices the protocol names.
    thresholds
        What those oligos were designed and judged by, for each role, so a page prints the band
        beside the value.
    """

    vector: SequenceRecord
    inserts: tuple[SequenceRecord, ...]
    span: tuple[int, int]
    product: AssemblyProduct
    linearised_vector: Part
    insert_parts: tuple[Part, ...]
    assembly: Assembly
    colony: ColonyCheck
    reads: tuple[SangerRead, SangerRead]
    amounts: tuple[Amount, ...]
    phenotype: Phenotype
    designed_oligos: tuple[DesignedOligo, ...]
    ordered_oligos: tuple[OrderedOligo, ...]
    host: str
    polymerase: Polymerase
    thresholds: Mapping[PrimerRole, Thresholds] = THRESHOLDS_FOR

    @property
    def plasmid(self) -> SequenceRecord:
        """The circular plasmid the assembly makes.

        Named apart from `product`, which is the kit the reaction is run with.
        """
        return self.assembly.product

    @property
    def parts(self) -> tuple[Part, ...]:
        """The parts that go into the reaction, the linearised vector first."""
        return (self.linearised_vector, *self.insert_parts)

    @property
    def overlaps(self) -> tuple[str, ...]:
        """The bases carried at each junction, in the product's own order.

        `liulab_mbio.cloning.gibson.assembly.Junction` says which molecule carries them.
        """
        return tuple(one.overlap for one in self.assembly.junctions)

    @property
    def junction_names(self) -> tuple[tuple[str, str], ...]:
        """Each junction's name and the bases carried there, in the product's own order."""
        return tuple((f"{one.before}-{one.after}", one.overlap) for one in self.assembly.junctions)

    @property
    def oligos(self) -> tuple[Primer, ...]:
        """Every oligo the plan designs, in the order the sheet lists them."""
        return tuple(report.primer for report in self.reports)

    @property
    def reports(self) -> tuple[PrimerReport, ...]:
        """Every designed oligo's evaluation, in the order the sheet lists them."""
        return tuple(oligo.report for oligo in self.designed_oligos)

    @property
    def checks(self) -> tuple[Check, ...]:
        """What the plan is judged on: the plasmid, every overlap, the reaction and the oligos.

        The ones no sourced threshold judges carry no verdict, and say in their detail that
        nothing measured them.
        """
        return (
            *self.assembly.checks,
            *overlap_checks(self.junction_names, self.product),
            *stitch_checks(self.parts),
            fragment_check(self.product, len(self.insert_parts)),
            assembly_dna_check(self.product, self.amounts),
            primer_check(self.reports),
        )

    @property
    def status(self) -> Status:
        """The worst status of any check."""
        return status(self.checks)

    def protocol(self) -> Protocol:
        """Return the bench protocol for this plan."""
        return protocol_for(
            vector=self.vector,
            span=self.span,
            product=self.product,
            linearised_vector=self.linearised_vector,
            insert_parts=self.insert_parts,
            assembly=self.assembly,
            colony=self.colony,
            reads=self.reads,
            amounts=self.amounts,
            phenotype=self.phenotype,
            oligos=self.designed_oligos,
            ordered=self.ordered_oligos,
            checks=self.checks,
            host=self.host,
            polymerase=self.polymerase,
            thresholds=self.thresholds,
        )

    def write(self, directory: str | os.PathLike[str]) -> Files:
        """Write the plasmid, the oligo sheet, the protocol data and its page into `directory`.

        The directory is made when it is not there. The four files are named by
        `liulab_mbio.cloning.plan`, and a second run over the same inputs writes the same bytes.
        """
        out = Path(directory)
        out.mkdir(parents=True, exist_ok=True)
        plasmid = out / PRODUCT_FILE
        write_dna(self.plasmid, plasmid)
        sheet = out / PRIMER_FILE
        sheet.write_text(primer_sheet(self.reports, oligos=self.ordered_oligos), encoding="utf-8")
        written = write_protocol_files(self.protocol(), out)
        return Files(plasmid, sheet, written.data, written.page)


def plan_gibson(
    vector: SequenceRecord | str | os.PathLike[str],
    *inserts: SequenceRecord | str | os.PathLike[str],
    site: Site = None,
    orientation: Orientation | Sequence[Orientation] = "forward",
    route: Route | Sequence[Route] = "amplify",
    bridge: Sequence[tuple[str, str]] = (),
    product: AssemblyProduct = NEBUILDER_HIFI,
    polymerase: Polymerase = Q5,
    host: str = DEFAULT_HOST,
    name: str = "",
    thresholds: Mapping[PrimerRole, Thresholds] = THRESHOLDS_FOR,
) -> Plan:
    """Plan one Gibson experiment putting `inserts` into `vector`.

    The inserts go round the product in the order they are given. A circular vector is opened
    by PCR across the span they replace; one that is already linear is taken as the opened part
    as it is. Every junction takes its overlap at the length and melting temperature `product`
    documents, by the rules `liulab_mbio.cloning.gibson.design` states, and the part on the
    other side carries those bases as a 5' tail, so the junction adds nothing.

    **Which side the bases come from** is the note's vector rule: a junction with the vector on
    either side takes the vector's own bases, so the opened backbone needs no tail and can be
    made once and reused. Every other junction takes the bases of the insert before it.

    A colony PCR and two sequencing primers are then designed over
    `liulab_mbio.cloning.gibson.assembly.Assembly.boundaries`, where one part gives way to the
    next: the overlap is one part's own bases, so the boundary is not where the junction begins.

    Parameters
    ----------
    vector, inserts
        A record, or a path to a ``.dna``, GenBank or FASTA file holding one. At least one
        insert, in the order they go round the product.
    site
        Where the inserts go: a feature name, a ``(start, end)`` span of the vector, or
        ``None`` to use the vector's own `liulab_mbio.cloning.plan.MCS_FEATURE` feature. A
        vector handed in linear replaces nothing, so it takes no site.
    orientation
        ``"reverse"`` puts the other strand of an insert into the product. One value covers
        every insert; a sequence gives one for each.
    route
        ``"stitch"`` makes an insert out of overlapping oligos tiling both strands instead of
        out of a PCR, which suits a linker, a tag or a short promoter. One value covers every
        insert; a sequence gives one for each.
    bridge
        Junctions joined by one oligo carrying homology to both ends, each named by the two
        parts it joins in the order they go round the product. Neither of the two then carries
        a tail, so a fragment amplified for something else goes in as it is.
    product
        The assembly product on the bench, which sets the overlap rule, the reaction, the
        incubation, the molar ratio and the fragment count it is documented for.
    polymerase
        For the PCRs.
    host, name
        The strain the protocol names, and what to call the plasmid.
    thresholds
        For each role, what its oligos are designed and judged by in `liulab_mbio.primers`.

    Returns
    -------
    Plan
        The design, the simulated plasmid and what the bench takes.

    Raises
    ------
    ValueError
        If no insert is given, if an orientation is neither ``forward`` nor ``reverse`` or a
        route neither ``amplify`` nor ``stitch``, or there is not one per insert, if a bridge
        names a junction this assembly does not have, if a site is named for a vector that is
        already linear, if no insertion site is named and a circular vector annotates none, if
        `product` documents no single-stranded oligo and either oligo route is asked for, if a
        part is too long to stitch or too short to take an overlap from, or if no annealing
        region fits a part's end.
    """
    if not inserts:
        raise ValueError("a Gibson plan needs a vector and at least one insert")
    one = as_record(vector)
    linear = one.topology == "linear"
    if linear and site is not None:
        raise ValueError(
            f"{one.name or 'the vector'} is already linear, so it is used as it is and there "
            "is no span to replace; drop the site"
        )
    ways = orientations(orientation, len(inserts))
    going = tuple(
        flipped(read) if way == "reverse" else read
        for read, way in (
            (as_record(record), way) for record, way in zip(inserts, ways, strict=True)
        )
    )
    labels = [record.name or f"insert {number}" for number, record in enumerate(going, start=1)]
    routes = _routes(route, len(going))
    span = (len(one), len(one)) if linear else insertion_span(one, site)
    opened = (0, len(one)) if linear else (span[1], span[0] + len(one))
    backbone = f"{one.name} backbone".strip()
    joins = tuple(zip((backbone, *labels), (*labels, backbone), strict=True))
    bridged = _bridged(bridge, joins)
    if "stitch" in routes:
        requires_oligos(product, "stitch a part out of overlapping oligos")
    if bridged:
        requires_oligos(product, "join two fragments with a bridging oligo")
    rule = product.tier(len(going) + 1).overlap
    tails = _tails(one, going, opened, rule, bridged)
    linearised_vector = (
        given_vector(one, name=backbone)
        if linear
        else open_vector(
            one,
            *span,
            name=backbone,
            polymerase=polymerase,
            thresholds=thresholds["amplification"],
        )
    )
    insert_parts = tuple(
        _insert_part(
            record,
            label,
            way,
            (left, right),
            product=product,
            polymerase=polymerase,
            thresholds=thresholds["amplification"],
        )
        for label, record, way, (left, right) in zip(labels, going, routes, tails, strict=True)
    )
    parts, bridges = _bridge_parts(
        [linearised_vector, *insert_parts], joins, bridged, product=product
    )
    linearised_vector, insert_parts = parts[0], tuple(parts[1:])
    built = assemble(parts, name=name or "-".join([one.name, *labels]).strip("-"))
    boundaries = built.boundaries
    colony = colony_pcr_check(
        built.product,
        boundaries,
        vector=one,
        insert_primer=True,
        polymerase=ONETAQ,
        thresholds=thresholds["colony PCR"],
    )
    reads = sanger_primers(built.product, boundaries, thresholds=thresholds["sequencing"])
    return Plan(
        one,
        going,
        span,
        product,
        linearised_vector,
        insert_parts,
        built,
        colony,
        reads,
        assembly_amounts(
            (linearised_vector.name, linearised_vector.fragment_length),
            tuple((part.name, part.fragment_length) for part in insert_parts),
            product=product,
        ),
        read_phenotype(built.product, built.insert_span, vector=one, span=span),
        (
            *(
                DesignedOligo(report, "amplification", part)
                for part in parts
                if part.report is not None
                for report in (part.report.forward, part.report.reverse)
            ),
            *(DesignedOligo(report, "colony PCR") for report in colony.reports),
            *(
                DesignedOligo(
                    evaluate_primer(
                        read.primer, built.product, thresholds=thresholds["sequencing"]
                    ),
                    "sequencing",
                )
                for read in reads
            ),
        ),
        _ordered_oligos(parts, bridges),
        host,
        polymerase,
        thresholds,
    )


def _routes(route: Route | Sequence[Route], count: int) -> tuple[Route, ...]:
    """Spread one route over every insert, or take the one given for each.

    Raises
    ------
    ValueError
        If a value is neither ``amplify`` nor ``stitch``, or there is not one per insert.
    """
    given = (route,) * count if isinstance(route, str) else tuple(route)
    if len(given) != count:
        raise ValueError(f"route has {len(given)} values for {count} insert(s)")
    for one in given:
        if one not in ("amplify", "stitch"):
            raise ValueError(f"route is 'amplify' or 'stitch', got {one!r}")
    return given


def _bridged(
    bridge: Sequence[tuple[str, str]], joins: Sequence[tuple[str, str]]
) -> tuple[int, ...]:
    """Return which junction each named pair of parts is, in the order the product runs.

    Raises
    ------
    ValueError
        If a pair does not name two parts that meet.
    """
    found = []
    for pair in bridge:
        named = tuple(pair)
        if named not in joins:
            meeting = ", ".join(f"{before}-{after}" for before, after in joins)
            raise ValueError(
                f"no junction joins {' to '.join(str(one) for one in named)}; this assembly "
                f"joins {meeting}"
            )
        found.append(joins.index(named))
    return tuple(sorted(set(found)))


def _tails(
    vector: SequenceRecord,
    inserts: Sequence[SequenceRecord],
    opened: tuple[int, int],
    rule: OverlapRule,
    bridged: Sequence[int],
) -> tuple[tuple[str, str], ...]:
    """Return the two tails each insert carries, in insert order.

    `opened` is the vector span that reaches the product, so the outer junctions read the
    vector's bases at its two ends whether it was opened by PCR or handed in linear. An inner
    junction reads the last bases of the insert before it, and the insert after it tails them,
    so exactly one side of every junction carries the bases and neither the vector nor the last
    insert is asked for a tail it cannot spell. A bridged junction leaves both sides bare: its
    oligo carries the bases instead.
    """
    first, last = opened
    lefts = ["" if 0 in bridged else overlap_before(vector, last, rule)]
    lefts += [
        "" if index in bridged else overlap_before(record, len(record), rule)
        for index, record in enumerate(inserts[:-1], start=1)
    ]
    rights = [""] * (len(inserts) - 1)
    rights.append("" if len(inserts) in bridged else overlap_after(vector, first, rule))
    return tuple(zip(lefts, rights, strict=True))


def _insert_part(
    record: SequenceRecord,
    label: str,
    route: Route,
    tails: tuple[str, str],
    *,
    product: AssemblyProduct,
    polymerase: Polymerase,
    thresholds: Thresholds,
) -> Part:
    """Make one insert the way it was asked for: out of a PCR, or out of overlapping oligos."""
    left, right = tails
    if route == "stitch":
        molecule = left + record.sequence + right
        return stitch(
            record,
            stitch_oligos(molecule, name=label, product=product),
            left_tail=left,
            right_tail=right,
            name=label,
        )
    return amplify(
        record,
        0,
        len(record),
        left_tail=left,
        right_tail=right,
        name=label,
        polymerase=polymerase,
        thresholds=thresholds,
    )


def _bridge_parts(
    parts: list[Part],
    joins: Sequence[tuple[str, str]],
    bridged: Sequence[int],
    *,
    product: AssemblyProduct,
) -> tuple[list[Part], tuple[tuple[str, str, str], ...]]:
    """Design an oligo for every bridged junction and hang it on the part that follows it.

    Returns the parts and, for each bridge, the two parts it joins and its sequence.
    """
    designed = []
    for index in bridged:
        before, after = joins[index]
        following = (index + 1) % len(parts)
        oligo = bridging_oligo(
            parts[index].bases,
            parts[following].bases,
            name=f"{before}-{after} bridge",
            product=product,
        )
        parts[following] = dataclasses.replace(parts[following], bridge=oligo)
        designed.append((before, after, oligo.sequence))
    return parts, tuple(designed)


def _ordered_oligos(
    parts: Sequence[Part], bridges: Sequence[tuple[str, str, str]]
) -> tuple[OrderedOligo, ...]:
    """Return every oligo the plan orders that primes nothing, with what it is for.

    A stitched part's oligos come first, then any bridging oligo. None of them anneals to a
    template, so each row carries its concentration and its purity and no verdict at all.
    """
    return (
        *(
            OrderedOligo(
                oligo.name,
                oligo.sequence,
                purpose=f"Stitch {part.name}",
                stock=f"{STITCH_OLIGO_NM:g} nM of each in the assembly",
                note=f"{OLIGO_PURITY} {NO_VERDICT}",
            )
            for part in parts
            for oligo in part.oligos
        ),
        *(
            OrderedOligo(
                f"{before}-{after} bridge",
                sequence,
                purpose=f"Bridge {before} to {after}",
                stock=f"{BRIDGE_OLIGO_PMOL:g} pmol in the assembly",
                note=f"{OLIGO_PURITY} {NO_VERDICT}",
            )
            for before, after, sequence in bridges
        ),
    )

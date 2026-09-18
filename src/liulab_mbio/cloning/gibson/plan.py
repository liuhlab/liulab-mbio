"""One Gibson experiment, planned from a vector and the insert that goes into it.

`plan_gibson` runs the whole design: open the vector by PCR across the span the insert
replaces, choose the overlap at each junction, design the primers that carry it, simulate the
product and work out what the assembly reaction takes. `Plan.write` puts four files in one
directory -- the annotated product, an oligo order sheet, the protocol as JSON data, and the
interactive HTML page rendered from that data.

Every number the protocol prints is computed here or by the modules this one calls, and every
supplier's number behind them is `liulab_mbio.cloning.gibson.bench`, through
``docs/research/gibson-assembly.md``. What the protocol says about the phenotype is
`liulab_mbio.bench.phenotype`, read off the product's own features.
"""

import os
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.oligos import primer_sheet
from liulab_mbio.bench.phenotype import Phenotype, read_phenotype
from liulab_mbio.checks import Check, Status
from liulab_mbio.cloning.gibson.assembly import Assembly, Part, amplify, assemble, open_vector
from liulab_mbio.cloning.gibson.bench import (
    NEBUILDER_HIFI,
    AssemblyProduct,
    assembly_amounts,
)
from liulab_mbio.cloning.gibson.design import overlap_after, overlap_before
from liulab_mbio.cloning.gibson.oligos import DesignedOligo
from liulab_mbio.cloning.gibson.steps import DEFAULT_HOST
from liulab_mbio.cloning.gibson.steps import protocol as protocol_for
from liulab_mbio.cloning.plan import (
    PRIMER_FILE,
    PRODUCT_FILE,
    Site,
    as_record,
    insertion_span,
    primer_check,
    status,
    write_protocol_files,
)
from liulab_mbio.primers.evaluation import PrimerReport
from liulab_mbio.primers.polymerase import Q5, Polymerase
from liulab_mbio.primers.thresholds import THRESHOLDS_FOR, PrimerRole, Thresholds
from liulab_mbio.protocol.model import Protocol
from liulab_mbio.sequence import Primer, SequenceRecord
from liulab_mbio.snapgene import write_dna


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
        The records the plan was made from.
    span
        The vector bases the insert replaces.
    product
        The assembly product on the bench, which sets the overlap rule, the reaction and the
        incubation.
    linearised_vector
        The part the vector is opened into.
    insert_parts
        The part each insert is amplified into, in insert order.
    assembly
        The simulated plasmid, the parts and the junctions.
    amounts
        What to put in the assembly reaction, vector first.
    phenotype
        What the plasmid says about itself.
    designed_oligos
        Every designed oligo and what it is for, in order: each part's forward and reverse
        primer.
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
    amounts: tuple[Amount, ...]
    phenotype: Phenotype
    designed_oligos: tuple[DesignedOligo, ...]
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
        """The bases shared at each junction, in the product's own order."""
        return tuple(one.overlap for one in self.assembly.junctions)

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
        """The plasmid's checks, with one more for the oligos."""
        return (*self.assembly.checks, primer_check(self.reports))

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
            amounts=self.amounts,
            phenotype=self.phenotype,
            oligos=self.designed_oligos,
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
        sheet.write_text(primer_sheet(self.reports), encoding="utf-8")
        written = write_protocol_files(self.protocol(), out)
        return Files(plasmid, sheet, written.data, written.page)


def plan_gibson(
    vector: SequenceRecord | str | os.PathLike[str],
    *inserts: SequenceRecord | str | os.PathLike[str],
    site: Site = None,
    product: AssemblyProduct = NEBUILDER_HIFI,
    polymerase: Polymerase = Q5,
    host: str = DEFAULT_HOST,
    name: str = "",
    thresholds: Mapping[PrimerRole, Thresholds] = THRESHOLDS_FOR,
) -> Plan:
    """Plan one Gibson experiment putting `inserts` into `vector`.

    The vector is opened by PCR across the span the insert replaces, and each of the two
    junctions takes its overlap from the vector's own bases there -- the length and melting
    temperature `product` documents, by the rules
    `liulab_mbio.cloning.gibson.design` states. The insert's primers carry those bases as 5'
    tails, so the junctions add nothing and the vector's PCR needs no tail at all.

    Parameters
    ----------
    vector, inserts
        A record, or a path to a ``.dna``, GenBank or FASTA file holding one. One insert; a
        circular vector.
    site
        Where the insert goes: a feature name, a ``(start, end)`` span of the vector, or
        ``None`` to use the vector's own `liulab_mbio.cloning.plan.MCS_FEATURE` feature.
    product
        The assembly product on the bench, which sets the overlap rule, the reaction, the
        incubation and the molar ratio.
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
        If no insert is given or more than one, if the vector is not circular, if no insertion
        site is named and the vector annotates none, if the vector is too short to take an
        overlap from, or if no annealing region fits a part's end.
    """
    if len(inserts) != 1:
        raise ValueError(
            f"a Gibson plan takes a vector and one insert, got {len(inserts)}: several inserts "
            "are not designed yet"
        )
    one = as_record(vector)
    if one.topology != "circular":
        raise ValueError(
            f"{one.name or 'the vector'} is linear, and a Gibson plan opens a circular vector "
            "by PCR; a backbone that is already linear is not taken yet"
        )
    going = tuple(as_record(record) for record in inserts)
    labels = [record.name or f"insert {number}" for number, record in enumerate(going, start=1)]
    start, end = insertion_span(one, site)
    rule = product.tier(len(going) + 1).overlap
    left, right = overlap_before(one, start, rule), overlap_after(one, end, rule)
    linearised_vector = open_vector(
        one,
        start,
        end,
        name=f"{one.name} backbone".strip(),
        polymerase=polymerase,
        thresholds=thresholds["amplification"],
    )
    insert_parts = tuple(
        amplify(
            record,
            0,
            len(record),
            left_tail=left,
            right_tail=right,
            name=label,
            polymerase=polymerase,
            thresholds=thresholds["amplification"],
        )
        for label, record in zip(labels, going, strict=True)
    )
    parts = (linearised_vector, *insert_parts)
    built = assemble(parts, name=name or "-".join([one.name, *labels]).strip("-"))
    return Plan(
        one,
        going,
        (start, end),
        product,
        linearised_vector,
        insert_parts,
        built,
        assembly_amounts(
            (linearised_vector.name, linearised_vector.fragment_length),
            tuple((part.name, part.fragment_length) for part in insert_parts),
            tier=product.tier(len(parts)),
        ),
        read_phenotype(built.product, built.insert_span, vector=one, span=(start, end)),
        tuple(
            DesignedOligo(report, "amplification", part)
            for part in parts
            for report in (part.report.forward, part.report.reverse)
        ),
        host,
        polymerase,
        thresholds,
    )

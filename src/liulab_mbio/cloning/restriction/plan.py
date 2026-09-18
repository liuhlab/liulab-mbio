"""One restriction and ligation cloning, planned from a vector and the plasmid the insert is in.

`plan_restriction` runs the whole design: cut both plasmids with the enzymes named, take the
backbone out of one and the insert out of the other, check the two ends anneal, simulate the
ligation, and design the colony PCR and the sequencing that confirm the clone. `Plan.write` puts
four files in one directory -- the annotated product, an oligo order sheet, the protocol as JSON
data, and the interactive HTML page rendered from that data.

Every number the protocol prints is computed here or by the modules this one calls, and
`liulab_mbio.cloning.restriction.bench` is where each bench number's source is written down.
What the protocol says about the phenotype is `liulab_mbio.bench.phenotype`, read off the
product's own features.
"""

import dataclasses
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.oligos import primer_sheet
from liulab_mbio.bench.phenotype import Phenotype, read_phenotype
from liulab_mbio.bench.validation import ColonyCheck, SangerRead, colony_pcr_check, sanger_primers
from liulab_mbio.checks import Check, Status
from liulab_mbio.cloning.plan import (
    PRIMER_FILE,
    PRODUCT_FILE,
    as_record,
    primer_check,
    status,
    write_protocol_files,
)
from liulab_mbio.cloning.restriction.bench import digest_amount, ligation_amounts
from liulab_mbio.cloning.restriction.digest import (
    Diagnostic,
    Piece,
    diagnostic,
    excised,
    opened,
    resolve,
)
from liulab_mbio.cloning.restriction.ligation import Junction, Ligation, ligate
from liulab_mbio.cloning.restriction.steps import DEFAULT_HOST
from liulab_mbio.cloning.restriction.steps import protocol as protocol_for
from liulab_mbio.cloning.restriction.verdicts import (
    buffer_check,
    cleanup_check,
    diagnostic_check,
    frame_check,
    methylation_check,
    ratio_check,
    temperature_check,
)
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.primers.evaluation import PrimerReport, evaluate_primer
from liulab_mbio.primers.polymerase import ONETAQ
from liulab_mbio.primers.thresholds import THRESHOLDS_FOR, PrimerRole, Thresholds
from liulab_mbio.protocol.model import Protocol
from liulab_mbio.sequence import Primer, SequenceRecord
from liulab_mbio.sites import EnzymeLike
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
    """One planned restriction and ligation cloning.

    Parameters
    ----------
    vector, source
        The records the plan was made from: the plasmid the backbone comes out of, and the one
        the insert is cut out of.
    enzymes
        The enzymes both digests use, in the order they were named.
    vector_pieces, source_pieces
        Every piece each digest leaves, the one that goes on first. The rest is what the gel
        separates it from.
    ligation
        The simulated product, the pieces and the junctions.
    span
        The vector bases the insert replaces.
    digests
        What each digest takes, the vector's first.
    amounts
        What to put in the ligation, the backbone first.
    diagnostic
        The digest that says a miniprep carries the insert, and the bands it should give.
    colony
        The colony PCR that reads both junctions and tells a correct clone from an empty vector.
    reads
        A sequencing primer reading in from outside the first junction and the last.
    read_reports
        What each of those two scored on the product.
    phenotype
        What the product says about itself.
    host
        The competent strain the protocol names.
    thresholds
        What the oligos were designed and judged by, for each role, so a page prints the band
        beside the value.
    """

    vector: SequenceRecord
    source: SequenceRecord
    enzymes: tuple[Enzyme, ...]
    vector_pieces: tuple[Piece, ...]
    source_pieces: tuple[Piece, ...]
    ligation: Ligation
    span: tuple[int, int]
    digests: tuple[Amount, Amount]
    amounts: tuple[Amount, Amount]
    diagnostic: Diagnostic
    colony: ColonyCheck
    reads: tuple[SangerRead, SangerRead]
    read_reports: tuple[PrimerReport, ...]
    phenotype: Phenotype
    host: str
    thresholds: Mapping[PrimerRole, Thresholds] = THRESHOLDS_FOR

    @property
    def backbone(self) -> Piece:
        """The vector piece the insert goes into."""
        return self.vector_pieces[0]

    @property
    def insert(self) -> Piece:
        """The source piece that goes in."""
        return self.source_pieces[0]

    @property
    def product(self) -> SequenceRecord:
        """The circular plasmid the ligation makes."""
        return self.ligation.product

    @property
    def junctions(self) -> tuple[Junction, ...]:
        """Where the two pieces meet, and what each join now spells."""
        return self.ligation.junctions

    @property
    def reports(self) -> tuple[PrimerReport, ...]:
        """Every designed oligo's evaluation, in the order the sheet lists them."""
        return (*self.colony.reports, *self.read_reports)

    @property
    def oligos(self) -> tuple[Primer, ...]:
        """Every oligo the plan designs, in the order the sheet lists them."""
        return tuple(report.primer for report in self.reports)

    @property
    def checks(self) -> tuple[Check, ...]:
        """Every verdict the plan carries, in the order the bench meets them.

        The digest's four first, then the product's own, then what the ligation takes and what
        confirms it. `liulab_mbio.cloning.restriction.verdicts` is where each one's threshold and
        its source are written down, and the buffer is the one nothing sourced can judge.
        """
        return (
            buffer_check(self.enzymes),
            temperature_check(self.enzymes),
            cleanup_check(self.enzymes),
            methylation_check(self.enzymes, (self.vector, self.source)),
            *self.ligation.checks,
            frame_check(self.product, self.junctions),
            ratio_check(*self.amounts),
            diagnostic_check(self.diagnostic),
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
            source=self.source,
            enzymes=self.enzymes,
            vector_pieces=self.vector_pieces,
            source_pieces=self.source_pieces,
            ligation=self.ligation,
            digests=self.digests,
            amounts=self.amounts,
            diagnostic=self.diagnostic,
            colony=self.colony,
            reads=self.reads,
            read_reports=self.read_reports,
            phenotype=self.phenotype,
            checks=self.checks,
            host=self.host,
            thresholds=self.thresholds,
        )

    def write(self, directory: str | os.PathLike[str]) -> Files:
        """Write the product, the oligo sheet, the protocol data and its page into `directory`.

        The directory is made when it is not there. All four names are
        `liulab_mbio.cloning.plan`'s, and a second run over the same inputs writes the same bytes.
        """
        out = Path(directory)
        out.mkdir(parents=True, exist_ok=True)
        product = out / PRODUCT_FILE
        write_dna(self.product, product)
        sheet = out / PRIMER_FILE
        sheet.write_text(primer_sheet(self.reports), encoding="utf-8")
        written = write_protocol_files(self.protocol(), out)
        return Files(product, sheet, written.data, written.page)


def plan_restriction(
    vector: SequenceRecord | str | os.PathLike[str],
    insert: SequenceRecord | str | os.PathLike[str],
    *,
    enzymes: Sequence[EnzymeLike],
    host: str = DEFAULT_HOST,
    name: str = "",
    thresholds: Mapping[PrimerRole, Thresholds] = THRESHOLDS_FOR,
) -> Plan:
    """Plan cutting an insert out of one plasmid and ligating it into another.

    Both plasmids are cut with the same enzymes in one reaction. The longest piece of the vector
    is its backbone and the shortest piece of `insert` is what goes in; each is gel-purified from
    what the same digest leaves beside it. The insert is written on whichever strand closes the
    circle, so which way round it lay in its own plasmid does not matter.

    The junction is not scarless. The two ends came from a recognition site and ligating them
    puts that site back, so the product gains those bases and `Plan.junctions` says what each one
    spells.

    Parameters
    ----------
    vector, insert
        A record, or a path to a ``.dna``, GenBank or FASTA file holding one. Both are circular
        plasmids: the vector the insert goes into, and the plasmid it is cut out of.
    enzymes
        One or two enzymes, each an `liulab_mbio.enzymes.Enzyme` or a name the package ships.
        Each must read exactly one site in each plasmid.
    host, name
        The competent strain the protocol names, and what to call the product. The colony PCR
        uses OneTaq, which is what NEB's protocol asks for.
    thresholds
        For each role, what its oligos are designed and judged by in `liulab_mbio.primers`.

    Returns
    -------
    Plan
        The digests, the simulated product and the validation.

    Raises
    ------
    KeyError
        If a name is not one the package ships.
    ValueError
        If either plasmid is not circular, if an enzyme reads anything but one site in either,
        if the backbone's two ends anneal to each other, or if the insert's ends do not anneal
        to the backbone's.
    """
    into = as_record(vector)
    holder = as_record(insert)
    chosen = resolve(enzymes)
    vector_pieces = opened(into, chosen)
    source_pieces = excised(holder, chosen, into=vector_pieces[0])
    backbone, released = vector_pieces[0], source_pieces[0]
    built = ligate(
        backbone,
        released,
        name=name or "-".join(part for part in (into.name, holder.name) if part) or "product",
    )
    span = _replaced(backbone, len(into))
    junctions = built.junction_positions
    colony = colony_pcr_check(
        built.product,
        junctions,
        vector=into,
        polymerase=ONETAQ,
        thresholds=thresholds["colony PCR"],
    )
    reads = sanger_primers(built.product, junctions, thresholds=thresholds["sequencing"])
    built = _annotated(built, (*colony.primers, *(read.primer for read in reads)))
    return Plan(
        into,
        holder,
        chosen,
        vector_pieces,
        source_pieces,
        built,
        span,
        (digest_amount((into.name, len(into))), digest_amount((holder.name, len(holder)))),
        ligation_amounts((backbone.name, backbone.length), (released.name, released.length)),
        diagnostic(built.product, into, chosen),
        colony,
        reads,
        tuple(
            evaluate_primer(read.primer, built.product, thresholds=thresholds["sequencing"])
            for read in reads
        ),
        read_phenotype(built.product, (junctions[0], junctions[-1]), vector=into, span=span),
        host,
        thresholds,
    )


def _annotated(built: Ligation, designed: Sequence[Primer]) -> Ligation:
    """Return the ligation with its designed oligos drawn on the product where they anneal.

    They were designed on the product, so each already carries the site it binds; what they were
    missing is a place on the record the user opens.
    """
    product = built.product
    return dataclasses.replace(
        built,
        product=dataclasses.replace(product, primers=(*product.primers, *designed)),
    )


def _replaced(backbone: Piece, length: int) -> tuple[int, int]:
    """Return the vector bases the insert replaces: everything the backbone is not."""
    return backbone.end % length, backbone.start

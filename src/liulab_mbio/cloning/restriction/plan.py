"""One restriction and ligation cloning, planned from a vector and where the insert comes from.

`plan_restriction` runs the whole design: open the vector with the enzymes named, get the insert
to the same two ends, check they anneal, simulate the ligation, and design the colony PCR and
the sequencing that confirm the clone. `Plan.write` puts four files in one directory -- the
annotated product, an oligo order sheet, the protocol as JSON data, and the interactive HTML page
rendered from that data.

There are two routes to the insert and the record handed in picks one. A record already carrying
the enzymes' sites is cut and the piece between them goes in. A record carrying none is amplified
first, with a spacer and a recognition site on each primer's 5' tail, and the amplicon is cut
instead -- `liulab_mbio.cloning.restriction.amplify` is that route.

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
from liulab_mbio.cloning.restriction.amplify import Amplicon, amplified
from liulab_mbio.cloning.restriction.bench import digest_amount, ligation_amounts
from liulab_mbio.cloning.restriction.design import Refusal, choose_pair, refusal
from liulab_mbio.cloning.restriction.digest import (
    Diagnostic,
    Piece,
    diagnostic,
    excised,
    opened,
    resolve,
)
from liulab_mbio.cloning.restriction.ligation import Junction, Ligation, ligate
from liulab_mbio.cloning.restriction.oligos import DesignedOligo
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
from liulab_mbio.primers.polymerase import ONETAQ, Q5, Polymerase
from liulab_mbio.primers.thresholds import THRESHOLDS_FOR, PrimerRole, Thresholds
from liulab_mbio.protocol.model import Protocol
from liulab_mbio.sequence import Primer, SequenceRecord
from liulab_mbio.sites import EnzymeLike, find_sites
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
        the insert comes from -- the plasmid it is cut out of, or the insert itself.
    enzymes
        The enzymes both digests use, in the order they were named.
    vector_pieces, source_pieces
        Every piece each digest leaves, the one that goes on first. The rest is what the gel
        separates it from.
    amplicon
        The PCR that put a site on each end of the insert, or ``None`` where `source` already
        carried them and the insert was cut out of it.
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
    designed_oligos
        Every oligo the plan designs and what it is for, in the order the sheet lists them: the
        amplification pair where there is one, then the colony PCR, then the sequencing.
    phenotype
        What the product says about itself.
    host, polymerase
        The competent strain the protocol names, and what the insert is amplified with.
    thresholds
        What the oligos were designed and judged by, for each role, so a page prints the band
        beside the value.
    refusals
        Every pair the chooser weighed and refused, with the rule that refused each. Empty where
        the enzymes were named, because then nothing was weighed against them.
    """

    vector: SequenceRecord
    source: SequenceRecord
    enzymes: tuple[Enzyme, ...]
    vector_pieces: tuple[Piece, ...]
    source_pieces: tuple[Piece, ...]
    amplicon: Amplicon | None
    ligation: Ligation
    span: tuple[int, int]
    digests: tuple[Amount, Amount]
    amounts: tuple[Amount, Amount]
    diagnostic: Diagnostic
    colony: ColonyCheck
    reads: tuple[SangerRead, SangerRead]
    designed_oligos: tuple[DesignedOligo, ...]
    phenotype: Phenotype
    host: str
    polymerase: Polymerase
    thresholds: Mapping[PrimerRole, Thresholds] = THRESHOLDS_FOR
    refusals: tuple[Refusal, ...] = ()

    @property
    def backbone(self) -> Piece:
        """The vector piece the insert goes into."""
        return self.vector_pieces[0]

    @property
    def insert(self) -> Piece:
        """The source piece that goes in."""
        return self.source_pieces[0]

    @property
    def digested(self) -> SequenceRecord:
        """The record the insert's own digest cuts: the amplicon, or `source` where there is none.

        An amplicon carries no methylation, whatever the plasmid that templated it carried.
        """
        return self.source if self.amplicon is None else self.amplicon.record

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
        return tuple(oligo.report for oligo in self.designed_oligos)

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
            methylation_check(self.enzymes, (self.vector, self.digested)),
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
            amplicon=self.amplicon,
            ligation=self.ligation,
            digests=self.digests,
            amounts=self.amounts,
            diagnostic=self.diagnostic,
            colony=self.colony,
            reads=self.reads,
            oligos=self.designed_oligos,
            phenotype=self.phenotype,
            checks=self.checks,
            refusals=self.refusals,
            host=self.host,
            polymerase=self.polymerase,
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
    enzymes: Sequence[EnzymeLike] = (),
    polymerase: Polymerase = Q5,
    host: str = DEFAULT_HOST,
    name: str = "",
    thresholds: Mapping[PrimerRole, Thresholds] = THRESHOLDS_FOR,
) -> Plan:
    """Plan getting an insert to two cut ends and ligating it into `vector`.

    The vector is cut with the enzymes named and its longest piece is the backbone, gel-purified
    from what the same digest leaves beside it.

    `insert` reaches the same two ends by one of two routes, and what it carries picks one. A
    record already reading a site of either enzyme is cut, and the piece between the two sites is
    what goes in -- a plasmid the insert already sits in, or a fragment ordered with the sites on
    it. A record reading no site of either is amplified first: each primer carries a spacer and a
    recognition site as a 5' tail, and the amplicon is cut instead. Either way the insert is
    written on whichever strand closes the circle, so which way round it lay does not matter.

    Naming no enzyme hands the choice to `liulab_mbio.cloning.restriction.design.choose_pair`,
    which weighs every pair the package ships and leaves `Plan.refusals` saying what refused the
    rest. Naming enzymes wins over it, and an unusable pair is refused in the chooser's own words.

    The junction is not scarless. The two ends came from a recognition site and ligating them
    puts that site back, so the product gains those bases and `Plan.junctions` says what each one
    spells.

    Parameters
    ----------
    vector, insert
        A record, or a path to a ``.dna``, GenBank or FASTA file holding one. The vector is a
        circular plasmid; `insert` is the insert, or the plasmid it is cut out of.
    enzymes
        One or two enzymes, each an `liulab_mbio.enzymes.Enzyme` or a name the package ships.
        Each must read a site in the vector, and one in whatever is cut for the insert. Chosen
        for the caller where none is named.
    polymerase
        For the insert's PCR, where one is run.
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
        If the vector is not circular, if an enzyme reads anything but one site in what is cut,
        if the backbone's two ends anneal to each other, if no tail spells one site, if the
        insert's ends do not anneal to the backbone's, or if no pair can be chosen at all.
    """
    into = as_record(vector)
    holder = as_record(insert)
    chosen, refusals = _pair(into, holder, enzymes)
    vector_pieces = opened(into, chosen)
    backbone = vector_pieces[0]
    amplicon = (
        None
        if _carries_a_site(holder, chosen)
        else amplified(
            holder,
            into=backbone,
            polymerase=polymerase,
            thresholds=thresholds["amplification"],
        )
    )
    digested = holder if amplicon is None else amplicon.record
    source_pieces = excised(digested, chosen, into=backbone)
    released = source_pieces[0]
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
        amplicon,
        built,
        span,
        (digest_amount((into.name, len(into))), digest_amount((digested.name, len(digested)))),
        ligation_amounts((backbone.name, backbone.length), (released.name, released.length)),
        diagnostic(built.product, into, chosen),
        colony,
        reads,
        _designed(amplicon, colony, reads, built.product, thresholds),
        read_phenotype(built.product, (junctions[0], junctions[-1]), vector=into, span=span),
        host,
        polymerase,
        thresholds,
        refusals,
    )


def _pair(
    vector: SequenceRecord, source: SequenceRecord, enzymes: Sequence[EnzymeLike]
) -> tuple[tuple[Enzyme, ...], tuple[Refusal, ...]]:
    """Return the pair to cut with, and what the chooser refused where it was the one choosing.

    The freezer beats the chooser: enzymes the caller names are used, and refused in the
    chooser's own words where they will not do.

    Raises
    ------
    ValueError
        If the named pair will not do, or if no pair can be chosen.
    """
    if not enzymes:
        picked = choose_pair(vector, source)
        return picked.enzymes, picked.refusals
    named = resolve(enzymes)
    if (found := refusal(vector, source, named)) is not None:
        raise ValueError(found.detail)
    return named, ()


def _carries_a_site(record: SequenceRecord, enzymes: Sequence[Enzyme]) -> bool:
    """Whether any of these enzymes cuts `record`, which is what decides the route.

    A record reading a site is cut as it is, and refused where the count is wrong. One reading
    none of them has the sites put on its two ends by PCR instead.
    """
    return any(site.cuts for site in find_sites(record, enzymes))


def _designed(
    amplicon: Amplicon | None,
    colony: ColonyCheck,
    reads: Sequence[SangerRead],
    product: SequenceRecord,
    thresholds: Mapping[PrimerRole, Thresholds],
) -> tuple[DesignedOligo, ...]:
    """Return every oligo the plan orders, in the order the bench uses them."""
    pair = () if amplicon is None else (amplicon.report.forward, amplicon.report.reverse)
    return (
        *(DesignedOligo(report, "amplification") for report in pair),
        *(DesignedOligo(report, "colony PCR") for report in colony.reports),
        *(
            DesignedOligo(
                evaluate_primer(read.primer, product, thresholds=thresholds["sequencing"]),
                "sequencing",
            )
            for read in reads
        ),
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

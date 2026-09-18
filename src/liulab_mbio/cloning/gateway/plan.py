"""One Gateway experiment, planned from an entry clone and a destination vector.

`plan_gateway` runs the whole design: find the att sites in both records, simulate the LR
reaction, work out what the reaction takes, and read what the product says about itself.
`Plan.write` puts three files in one directory -- the annotated expression clone, the protocol
as JSON data, and the interactive HTML page rendered from that data. No oligo is designed yet,
so no order sheet is written; the file set grows with the tickets that design one.

Every number the protocol prints is computed here or is one `liulab_mbio.cloning.gateway.bench`
cites from `docs/research/gateway-cloning.md`.
"""

import os
from dataclasses import dataclass
from pathlib import Path

from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.phenotype import Phenotype, read_phenotype
from liulab_mbio.checks import Check, Status
from liulab_mbio.cloning.gateway.bench import DEFAULT_HOST, lr_amounts
from liulab_mbio.cloning.gateway.recombination import Junction, Recombination, recombine
from liulab_mbio.cloning.gateway.steps import protocol as protocol_for
from liulab_mbio.cloning.plan import PRODUCT_FILE, as_record, status, write_protocol_files
from liulab_mbio.protocol.model import Protocol
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.snapgene import write_dna


@dataclass(frozen=True, slots=True)
class Files:
    """The files a plan writes.

    Parameters
    ----------
    product
        The annotated expression clone, as a SnapGene ``.dna`` file.
    protocol_data
        The bench protocol as JSON, which ``protocol render`` turns back into a page.
    protocol
        The interactive bench protocol, as one self-contained HTML page rendered from
        `protocol_data`.
    """

    product: Path
    protocol_data: Path
    protocol: Path

    @property
    def paths(self) -> tuple[Path, ...]:
        """Every file written, the first written first."""
        return (self.product, self.protocol_data, self.protocol)


@dataclass(frozen=True, slots=True)
class Plan:
    """One planned Gateway experiment.

    Parameters
    ----------
    entry, destination
        The records the plan was made from, each read on the strand the reaction ran it on.
    recombination
        The simulated LR reaction: the expression clone, the two pieces that made it, and the
        att junctions between them.
    amounts
        What to put in the reaction, the entry clone first.
    phenotype
        What the product says about itself.
    host
        The strain the protocol names.
    """

    entry: SequenceRecord
    destination: SequenceRecord
    recombination: Recombination
    amounts: tuple[Amount, ...]
    phenotype: Phenotype
    host: str

    @property
    def product(self) -> SequenceRecord:
        """The expression clone the reaction makes."""
        return self.recombination.product

    @property
    def junctions(self) -> tuple[Junction, Junction]:
        """The two att sites the reaction wrote, site 1 first."""
        return self.recombination.junctions

    @property
    def checks(self) -> tuple[Check, ...]:
        """Every verdict the plan carries."""
        return self.recombination.checks

    @property
    def status(self) -> Status:
        """The worst status of any check."""
        return status(self.checks)

    def protocol(self) -> Protocol:
        """Return the bench protocol for this plan."""
        return protocol_for(
            entry=self.entry,
            destination=self.destination,
            recombination=self.recombination,
            amounts=self.amounts,
            phenotype=self.phenotype,
            checks=self.checks,
            host=self.host,
        )

    def write(self, directory: str | os.PathLike[str]) -> Files:
        """Write the product, the protocol data and its page into `directory`.

        The directory is made when it is not there. The files are named by `PRODUCT_FILE` and
        by `liulab_mbio.cloning.plan` for the protocol pair, and a second run over the same
        inputs writes the same bytes.
        """
        out = Path(directory)
        out.mkdir(parents=True, exist_ok=True)
        product = out / PRODUCT_FILE
        write_dna(self.product, product)
        written = write_protocol_files(self.protocol(), out)
        return Files(product, written.data, written.page)


def plan_gateway(
    entry: SequenceRecord | str | os.PathLike[str],
    destination: SequenceRecord | str | os.PathLike[str],
    *,
    host: str = DEFAULT_HOST,
    name: str = "",
) -> Plan:
    """Plan one LR reaction putting an entry clone's insert into a destination vector.

    The att sites are found in both records rather than matched against a constant, because the
    vendor states that the sequences drift between products. Each record is then read on the
    strand its site 1 lies on, the segment between the entry clone's attL sites is recombined
    into the destination vector's backbone, and the product is simulated base for base.

    Parameters
    ----------
    entry
        The entry clone, or a path to a ``.dna``, GenBank or FASTA file holding it.
    destination
        The destination vector, the same way. It must be circular.
    host
        The strain the protocol names for selecting the clone.
    name
        What to call the product.

    Returns
    -------
    Plan
        The simulated expression clone, what the reaction takes, and the checks.

    Raises
    ------
    ValueError
        If either record does not carry the pair of att sites its side of the reaction needs,
        pointing at each other; if the destination vector is not circular; or if the DNA does
        not fit the reaction.
    """
    one, other = as_record(entry), as_record(destination)
    made = recombine(
        one,
        other,
        reaction="LR",
        name=name or "-".join(part for part in (other.name, one.name) if part),
    )
    first, second = made.junctions
    return Plan(
        made.moved.record,
        made.backbone.record,
        made,
        lr_amounts((one.name, len(one)), (other.name, len(other))),
        read_phenotype(
            made.product,
            (first.start, second.end),
            vector=made.backbone.record,
            span=made.cassette,
        ),
        host,
    )

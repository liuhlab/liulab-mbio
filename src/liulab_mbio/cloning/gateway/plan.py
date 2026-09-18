"""One Gateway experiment, planned over both reactions from the records a bench holds.

`plan_gateway` runs the whole design: find the att sites in every record, simulate the
recombinations, work out what each reaction takes, and read what each product says about
itself. Given an entry clone it plans the LR reaction alone; given an attB-flanked insert and a
donor vector it plans BP first and feeds the entry clone that makes into LR.

`Plan.write` puts three files in one directory -- the annotated expression clone, the protocol
as JSON data, and the interactive HTML page rendered from that data -- and a fourth, the entry
clone, where BP was planned. No oligo is designed yet, so no order sheet is written; the file
set grows with the tickets that design one.

Every number the protocol prints is computed here or is one `liulab_mbio.cloning.gateway.bench`
cites from `docs/research/gateway-cloning.md`.
"""

import os
from dataclasses import dataclass
from pathlib import Path

from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.phenotype import read_phenotype
from liulab_mbio.checks import Check, Status
from liulab_mbio.cloning.gateway.bench import DEFAULT_HOST, bp_amounts, lr_amounts
from liulab_mbio.cloning.gateway.recombination import (
    Junction,
    PlannedReaction,
    Recombination,
    recombine,
)
from liulab_mbio.cloning.gateway.steps import protocol as protocol_for
from liulab_mbio.cloning.plan import PRODUCT_FILE, as_record, status, write_protocol_files
from liulab_mbio.protocol.model import Protocol
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.snapgene import write_dna

#: What a Gateway plan calls the entry clone BP makes, which is the one file no other method
#: writes. The rest are named by `liulab_mbio.cloning.plan`.
ENTRY_FILE = "entry-clone.dna"


@dataclass(frozen=True, slots=True)
class Files:
    """The files a plan writes.

    Parameters
    ----------
    entry
        The entry clone BP makes, as a SnapGene ``.dna`` file, or ``None`` where one was given
        and no BP reaction was planned.
    product
        The annotated expression clone, the same way.
    protocol_data
        The bench protocol as JSON, which ``protocol render`` turns back into a page.
    protocol
        The interactive bench protocol, as one self-contained HTML page rendered from
        `protocol_data`.
    """

    entry: Path | None
    product: Path
    protocol_data: Path
    protocol: Path

    @property
    def paths(self) -> tuple[Path, ...]:
        """Every file written, the first written first, whichever size the set is."""
        written = (self.entry, self.product, self.protocol_data, self.protocol)
        return tuple(path for path in written if path is not None)


@dataclass(frozen=True, slots=True)
class Plan:
    """One planned Gateway experiment: the LR reaction, and the BP that fed it or not.

    Parameters
    ----------
    lr
        The LR reaction, which every plan runs.
    bp
        The BP reaction before it, or ``None`` where an entry clone was given.
    host
        The strain the protocol names for selecting each clone.
    """

    lr: PlannedReaction
    bp: PlannedReaction | None
    host: str

    @property
    def product(self) -> SequenceRecord:
        """The expression clone the LR reaction makes."""
        return self.lr.product

    @property
    def entry(self) -> SequenceRecord:
        """The entry clone LR ran on: the one given, or the one BP made."""
        return self.lr.recombination.moved.record

    @property
    def destination(self) -> SequenceRecord:
        """The destination vector, read on the strand the reaction ran it on."""
        return self.lr.recombination.backbone.record

    @property
    def insert(self) -> SequenceRecord | None:
        """The attB-flanked DNA BP ran on, or ``None`` where no BP reaction was planned."""
        return None if self.bp is None else self.bp.recombination.moved.record

    @property
    def donor(self) -> SequenceRecord | None:
        """The donor vector, or ``None`` where no BP reaction was planned."""
        return None if self.bp is None else self.bp.recombination.backbone.record

    @property
    def route(self) -> str:
        """Which reactions were planned, for a summary line to name."""
        return "LR" if self.bp is None else "BP then LR"

    @property
    def junctions(self) -> tuple[Junction, Junction]:
        """The two att sites the LR reaction wrote, site 1 first."""
        return self.lr.junctions

    @property
    def checks(self) -> tuple[Check, ...]:
        """Every verdict the plan carries, the earlier reaction's first."""
        earlier = () if self.bp is None else self.bp.recombination.checks
        return (*earlier, *self.lr.recombination.checks)

    @property
    def status(self) -> Status:
        """The worst status of any check."""
        return status(self.checks)

    def protocol(self) -> Protocol:
        """Return the bench protocol for this plan."""
        return protocol_for(lr=self.lr, bp=self.bp, checks=self.checks, host=self.host)

    def write(self, directory: str | os.PathLike[str]) -> Files:
        """Write the entry clone, the product, the protocol data and its page into `directory`.

        The directory is made when it is not there. The entry clone is written only where BP
        was planned, under `ENTRY_FILE`; the rest are named by `PRODUCT_FILE` and by
        `liulab_mbio.cloning.plan` for the protocol pair. A second run over the same inputs
        writes the same bytes.
        """
        out = Path(directory)
        out.mkdir(parents=True, exist_ok=True)
        entry = None
        if self.bp is not None:
            entry = out / ENTRY_FILE
            write_dna(self.bp.product, entry)
        product = out / PRODUCT_FILE
        write_dna(self.product, product)
        written = write_protocol_files(self.protocol(), out)
        return Files(entry, product, written.data, written.page)


def plan_gateway(
    carrier: SequenceRecord | str | os.PathLike[str],
    destination: SequenceRecord | str | os.PathLike[str],
    *,
    donor: SequenceRecord | str | os.PathLike[str] | None = None,
    host: str = DEFAULT_HOST,
    name: str = "",
) -> Plan:
    """Plan the reactions that put one insert into a destination vector.

    With no donor vector the carrier is an entry clone and one LR reaction is planned. With a
    donor vector the carrier is an attB-flanked insert -- a synthesised fragment, or an
    amplicon made elsewhere -- and BP is planned first, its entry clone feeding LR.

    The att sites are found in every record rather than matched against a constant, because the
    vendor states that the sequences drift between products. Each record is then read on the
    strand its site 1 lies on, the segment between the carrier's sites is recombined into the
    acceptor's backbone, and each product is simulated base for base.

    Parameters
    ----------
    carrier
        The entry clone, or the attB-flanked insert when `donor` is given, or a path to a
        ``.dna``, GenBank or FASTA file holding it. An attB insert may be linear.
    destination
        The destination vector, the same way. It must be circular.
    donor
        The donor vector, the same way, to plan the BP reaction that makes the entry clone. It
        must be circular.
    host
        The strain the protocol names for selecting each clone.
    name
        What to call the expression clone.

    Returns
    -------
    Plan
        The reactions planned, what each takes, and the checks.

    Raises
    ------
    ValueError
        If any record does not carry the pair of att sites its side of a reaction needs,
        pointing at each other; if a vector that must be circular is not; or if the DNA does
        not fit a reaction.
    """
    start, other = as_record(carrier), as_record(destination)
    one, bp = start, None
    if donor is not None:
        vector = as_record(donor)
        bp = _planned(
            recombine(one, vector, reaction="BP", name=_named(vector.name, one.name)),
            bp_amounts((one.name, len(one)), (vector.name, len(vector))),
        )
        one = bp.product
    lr = _planned(
        recombine(one, other, reaction="LR", name=name or _named(other.name, start.name)),
        lr_amounts((one.name, len(one)), (other.name, len(other))),
    )
    return Plan(lr, bp, host)


def _planned(made: Recombination, amounts: tuple[Amount, Amount]) -> PlannedReaction:
    """Return one simulated reaction with what it takes and what its product says about itself."""
    first, second = made.junctions
    return PlannedReaction(
        made,
        amounts,
        read_phenotype(
            made.product,
            (first.start, second.end),
            vector=made.backbone.record,
            span=made.cassette,
        ),
    )


def _named(vector: str, carrier: str) -> str:
    """Return what to call a product, from the names of the two records that made it."""
    return "-".join(part for part in (vector, carrier) if part)

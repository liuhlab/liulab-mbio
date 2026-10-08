"""The part carrier: one part into its own plasmid, one well a part, nothing pooled at the end.

A part is a reusable in-frame element — a tag, a linker, a signal peptide, a localization
signal, a degron — flanked by inward-facing BsmBI sites. Each one goes into the part carrier so
the collection can be kept and re-cut later. It is not a round: a round joins one part list to a
library in one tube and hands back a pool, and this hands back one plasmid a well. The method
reserves BsmBI for releasing a part from its carrier and for the last transfer into a working
vector, which is why neither has a role inside the rounds to be named by.

The carrier is supplied already linearised, with topoisomerase I bound to each 3' end, so the
reaction adds no enzyme and the backbone must stay free of BsmBI. This module simulates the
plasmid that comes out: it says where each part sits, what cuts it back out and what comes out,
so a step that seats parts can be written without the round model pretending a library came of
it.
"""

import dataclasses
from collections.abc import Collection, Sequence
from dataclasses import dataclass

from liulab_mbio.bench import plates
from liulab_mbio.edits import EditReport, insert, ordered
from liulab_mbio.protocol.model import FORMATS, Material, Plate, Step, Troubleshooting
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand
from liulab_mbio.sites import find_sites

#: The enzyme that releases a part from its carrier, which the method reserves rather than
#: naming inside a round. It does not seat the part: nothing cuts the carrier open.
ENZYME = "BsmBI"

#: The carrier, whose backbone carries no BsmBI site, and the marker it selects on.
CARRIER = "pCR-Blunt II-TOPO"
CARRIER_MARKER = "KanR"

#: The kit the carrier is supplied in, linearised with topoisomerase I covalently bound.
CARRIER_KIT = "Zero Blunt TOPO PCR Cloning Kit"

#: The head-to-head pair topoisomerase I sits either side of. The enzyme cleaves after 5'-CCCTT
#: on each strand, so a carrier is supplied opened halfway through this run.
TOPO_SITE = "GCCCTTAAGGGC"


@dataclass(frozen=True, slots=True)
class SeatedParts:
    """Every part seated in its own carrier, and where each one sits.

    Parameters
    ----------
    plate
        One well a part, seated in reading order.
    parts
        The part names, in that order.
    records
        The seated plasmid of each well, in that order.

    Notes
    -----
    There is no library field, and that is the point: this ends in one plasmid a well, each
    still identified by where it sits. Neither this class nor this module is called `Seating`:
    `liulab_mbio` already uses that word for where things sit in a plate.
    """

    plate: Plate
    parts: tuple[str, ...]
    records: tuple[SequenceRecord, ...] = ()

    @property
    def products(self) -> int:
        """How many plasmids the carrier reaction makes: one a part, never a pool."""
        return len(self.parts)


def topo_site(carrier: SequenceRecord) -> int:
    """Return the blunt point `carrier` is supplied linearised at.

    Raises
    ------
    ValueError
        If the carrier carries no head-to-head `TOPO_SITE`, or more than one, since either
        leaves the point a part is seated at undetermined.

    Examples
    --------
    >>> topo_site(SequenceRecord("AAAGCCCTTAAGGGCTTT"))
    9
    """
    length = len(carrier)
    reach = carrier.sequence
    if carrier.topology == "circular":
        reach += carrier.sequence[: len(TOPO_SITE) - 1]
    starts = []
    at = reach.find(TOPO_SITE)
    while at != -1 and at < length:
        starts.append(at)
        at = reach.find(TOPO_SITE, at + 1)
    named = carrier.name or "the carrier"
    if not starts:
        raise ValueError(
            f"{named} carries no {TOPO_SITE} run, so there is no point topoisomerase I opens it "
            "at and nothing says where a part would be seated"
        )
    if len(starts) > 1:
        raise ValueError(
            f"{named} carries {len(starts)} {TOPO_SITE} runs, at "
            f"{', '.join(str(start) for start in starts)}: a part would be seated at whichever "
            "one, so the product is undetermined"
        )
    return (starts[0] + len(TOPO_SITE) // 2) % length


def released(part: SequenceRecord) -> tuple[str, str]:
    """Return the overhang pair `ENZYME` releases `part` on, left end first.

    Raises
    ------
    ValueError
        If the part is not flanked by exactly one inward-facing `ENZYME` pair, since a part
        that cannot be cut back out of its carrier is not a part.

    Examples
    --------
    >>> released(SequenceRecord("CGTCTCATATGGGCAGGATGAGACG"))
    ('TATG', 'AGGA')
    """
    sites = find_sites(part, ENZYME)
    named = part.name or "the part"
    if len(sites) != 2:
        raise ValueError(
            f"{named} carries {len(sites)} {ENZYME} site(s): a part is flanked by two, facing "
            "inward, so it can be cut back out of its carrier"
        )
    left, right = sites
    inward = left.strand is Strand.FORWARD and right.strand is Strand.REVERSE
    if not (inward and left.cuts and right.cuts and left.top_cut < right.top_cut):
        raise ValueError(
            f"{named} carries two {ENZYME} sites that do not face inward, so cutting it would "
            "not release the part between them"
        )
    return left.overhang or "", right.overhang or ""


def seat(
    carrier: SequenceRecord,
    part: SequenceRecord,
    *,
    overhangs: Collection[tuple[str, str]] | None = None,
) -> tuple[SequenceRecord, EditReport]:
    """Seat one part at the carrier's blunt point, and say what the insertion moved.

    Parameters
    ----------
    carrier
        The part carrier, free of `ENZYME` so that releasing the part leaves it whole.
    part
        The part, flanked by an inward-facing `ENZYME` pair.
    overhangs
        The pairs a released part may carry, such as
        `liulab_synbio.igga.method.INTERFACE_OVERHANGS` values. ``None`` accepts any pair.

    Raises
    ------
    ValueError
        If the carrier carries an `ENZYME` site or no single blunt point, if the part carries
        no inward `ENZYME` pair, or if it is released on a pair `overhangs` does not hold.
    """
    sites = find_sites(carrier, ENZYME)
    if sites:
        raise ValueError(
            f"{carrier.name or 'the carrier'} carries {len(sites)} {ENZYME} site(s): releasing "
            "a part from it would cut the backbone as well"
        )
    pair = released(part)
    if overhangs is not None and pair not in overhangs:
        wanted = ", ".join(f"{one}/{other}" for one, other in sorted(overhangs))
        raise ValueError(
            f"{part.name or 'the part'} is released on {pair[0]}/{pair[1]}, which is no pair the "
            f"method declares: {wanted}"
        )
    at = topo_site(carrier)
    seated, report = insert(carrier, at, part.sequence)
    feature = Feature(
        part.name or "part", "misc_feature", (Segment(at, at + len(part)),), strand=Strand.FORWARD
    )
    seated = dataclasses.replace(
        seated,
        name=f"{part.name} in {carrier.name}" if part.name and carrier.name else seated.name,
        features=(*seated.features, feature),
    )
    return ordered(seated), report


def seat_parts(
    parts: Sequence[SequenceRecord],
    *,
    carrier: SequenceRecord,
    name: str = "parts",
    overhangs: Collection[tuple[str, str]] | None = None,
) -> SeatedParts:
    """Seat each part in its own well of the smallest plate format that holds them all.

    A well is named by the part sitting in it, so a part is named by its record.

    Raises
    ------
    ValueError
        If no part is given, if one has no name, if two share a name, if more parts are given
        than the largest format holds, or if any part will not seat in the carrier.

    Examples
    --------
    >>> carrier = SequenceRecord("AAAAGCCCTTAAGGGCTTTT", topology="circular", name="carrier")
    >>> part = SequenceRecord("CGTCTCATATGGGCAGGATGAGACG", name="FLAG")
    >>> seated = seat_parts([part], carrier=carrier)
    >>> seated.plate.seating["A1"], seated.products, len(seated.records[0])
    ('FLAG', 1, 45)
    """
    named = tuple(part.name for part in parts)
    if not named:
        raise ValueError("no part was given, and a carrier plate holds at least one")
    if not all(named):
        raise ValueError("a part has no name, and a well is named by the part sitting in it")
    if len(set(named)) != len(named):
        raise ValueError("two parts share a name, and a well is named by the part sitting in it")
    fits = [wells for wells in sorted(FORMATS) if wells >= len(named)]
    if not fits:
        raise ValueError(
            f"{len(named)} parts do not fit the largest plate format, {max(FORMATS)} wells"
        )
    records = tuple(seat(carrier, part, overhangs=overhangs)[0] for part in parts)
    return SeatedParts(
        plates.plate(
            name,
            fits[0],
            seating=plates.seat(named, fits[0]),
            holds=f"one seated part each, selected on {CARRIER_MARKER}",
            note=f"{len(named)} of {fits[0]} wells used",
        ),
        named,
        records,
    )


def materials() -> tuple[Material, ...]:
    """Return what the carrier reaction consumes besides the parts themselves."""
    return (
        Material(
            CARRIER,
            storage="-20 °C",
            note=(
                f"the part carrier; its backbone carries no {ENZYME} site, which is what lets "
                f"{ENZYME} release the part again later and leave the backbone whole"
            ),
        ),
        Material(
            CARRIER_KIT,
            storage="-20 °C",
            amount="one reaction per well",
            note=(
                "supplies the carrier already linearised with topoisomerase I bound to each 3' "
                "end, so the reaction adds no ligase and no other enzyme"
            ),
        ),
    )


def carrier_step(seated: SeatedParts) -> Step:
    """Return the step that seats every part, which ends in plasmids and not in a library."""
    return Step(
        f"Seat {seated.products} part(s) in {CARRIER}",
        key="seat-parts",
        instructions=(
            f"Set up one {CARRIER_KIT} reaction per well of {seated.plate.name}, each holding "
            "one blunt part, the linearised carrier and the kit's salt solution.",
            "Leave 5 min at room temperature. Topoisomerase I comes bound to the carrier, so "
            "nothing is added to join the two.",
            f"Transform each well on its own and select on {CARRIER_MARKER}.",
        ),
        expected=(
            f"{seated.products} carrier plasmids, one a part, each still named by the well it "
            "sits in.",
            "No pool and no library: nothing here is mixed, so nothing has to be told apart "
            "afterwards.",
        ),
        notes=(
            f"{ENZYME} releases a part from its carrier afterwards, and makes the last transfer "
            "into a working vector. The method reserves it for those two, so a part block "
            "spells its site nowhere else and the carrier backbone carries none.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Two parts end up in one well",
                "The wells were pooled. One part a well is on purpose: a part that has met "
                "another cannot be re-cut back out on its own.",
            ),
        ),
    )

"""Seating: one part into its own carrier, one well a part, with nothing pooled at the end.

A part is a reusable in-frame element — a tag, a linker, a signal peptide, a localization
signal, a degron — flanked by inward-facing BsmBI sites. Seating puts each one in the part
carrier so the collection can be kept and re-cut later. It is not a round: a round joins one
part list to a library in one tube and hands back a pool, and this hands back one plasmid a
well. The method reserves BsmBI for exactly this and for the last transfer into a working
vector, which is why neither has a role inside the rounds to be named by.

Nothing here designs DNA. It says where each part sits, what cuts it and what comes out, so a
step that seats parts can be written without the round model pretending a library came of it.
"""

from collections.abc import Sequence
from dataclasses import dataclass

from liulab_mbio.bench import plates
from liulab_mbio.protocol.model import FORMATS, Material, Plate, Step, Troubleshooting

#: The enzyme that releases a part from its carrier and seats it there, which the method
#: reserves rather than naming inside a round.
ENZYME = "BsmBI"

#: The carrier, whose backbone carries no BsmBI site, and the marker it selects on.
CARRIER = "pCR-Blunt II-TOPO"
CARRIER_MARKER = "KanR"


@dataclass(frozen=True, slots=True)
class Seating:
    """Every part seated in its own carrier, and where each one sits.

    Parameters
    ----------
    plate
        One well a part, seated in reading order.
    parts
        The part names, in that order.

    Notes
    -----
    There is no library field, and that is the point: seating ends in one plasmid a well, each
    still identified by where it sits.
    """

    plate: Plate
    parts: tuple[str, ...]

    @property
    def products(self) -> int:
        """How many plasmids seating makes: one a part, never a pool."""
        return len(self.parts)


def seat_parts(parts: Sequence[str], *, name: str = "parts") -> Seating:
    """Seat each part in its own well of the smallest plate format that holds them all.

    Raises
    ------
    ValueError
        If no part is given, if two share a name, or if more parts are given than the largest
        format holds.

    Examples
    --------
    >>> seated = seat_parts(["FLAG", "linker", "degron"])
    >>> seated.plate.wells, seated.plate.seating["A3"], seated.products
    (12, 'degron', 3)
    """
    named = tuple(parts)
    if not named:
        raise ValueError("seating takes at least one part")
    if len(set(named)) != len(named):
        raise ValueError("two parts share a name, and a well is named by the part sitting in it")
    fits = [wells for wells in sorted(FORMATS) if wells >= len(named)]
    if not fits:
        raise ValueError(
            f"{len(named)} parts do not fit the largest plate format, {max(FORMATS)} wells"
        )
    return Seating(
        plates.plate(
            name,
            fits[0],
            seating=plates.seat(named, fits[0]),
            holds=f"one seated part each, selected on {CARRIER_MARKER}",
            note=f"{len(named)} of {fits[0]} wells used",
        ),
        named,
    )


def materials() -> tuple[Material, ...]:
    """Return what seating consumes besides the parts themselves."""
    return (
        Material(
            CARRIER,
            storage="-20 °C",
            note=f"the part carrier; its backbone carries no {ENZYME} site",
        ),
        Material(
            f"{ENZYME} and ligase",
            storage="-20 °C",
            amount="one reaction per well",
            note="one part a well, so nothing is pooled and no part meets another",
        ),
    )


def seating_step(seated: Seating) -> Step:
    """Return the step that seats every part, which ends in plasmids and not in a library."""
    return Step(
        f"Seat {seated.products} part(s) in {CARRIER}",
        instructions=(
            f"Set up one {ENZYME} and ligase reaction per well of {seated.plate.name}, each "
            "holding one part and the carrier.",
            f"Transform each well on its own and select on {CARRIER_MARKER}.",
        ),
        expected=(
            f"{seated.products} carrier plasmids, one a part, each still named by the well it "
            "sits in.",
            "No pool and no library: nothing here is mixed, so nothing has to be told apart "
            "afterwards.",
        ),
        notes=(
            f"{ENZYME} is reserved for this step and for the last transfer into a working "
            "vector, so a part block spells its site nowhere else.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Two parts end up in one well",
                "The wells were pooled. Seating keeps one part a well on purpose: a part that "
                "has met another cannot be re-cut back out on its own.",
            ),
        ),
    )

"""An enzyme's heat inactivation, as its supplier gives it."""

from collections.abc import Sequence

from liulab_mbio.bench.steps import listed
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.protocol.model import Incubation, Stage, ThermocyclerProgram


def heat_inactivation(enzyme: Enzyme) -> ThermocyclerProgram | None:
    """Return the supplier's heat inactivation for this enzyme, or ``None`` where it gives none.

    This is not the 60 °C end soak a Golden Gate program ends with, which is a digest.
    """
    return next(iter(heat_inactivations((enzyme,))), None)


def heat_inactivations(enzymes: Sequence[Enzyme]) -> tuple[ThermocyclerProgram, ...]:
    """Return the suppliers' heat inactivations, enzymes agreeing on one sharing a program.

    One tube holding several enzymes is held once for each temperature and time its suppliers
    give, so enzymes killed alike are killed together. An enzyme whose supplier states none is
    left out, and a tube of nothing but those gets no program at all.
    """
    shared: dict[tuple[int, int], list[Enzyme]] = {}
    for one in enzymes:
        celsius, minutes = one.heat_inactivation_celsius, one.heat_inactivation_minutes
        if celsius is not None and minutes is not None:
            shared.setdefault((celsius, minutes), []).append(one)
    return tuple(
        ThermocyclerProgram(
            (Stage((Incubation("Heat inactivation", float(celsius), minutes * 60),)),),
            title=f"{listed([one.name for one in named])} heat inactivation",
        )
        for (celsius, minutes), named in shared.items()
    )

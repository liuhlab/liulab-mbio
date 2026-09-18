"""An enzyme's heat inactivation, as its supplier gives it."""

from liulab_mbio.enzymes import Enzyme
from liulab_mbio.protocol.model import Incubation, Stage, ThermocyclerProgram


def heat_inactivation(enzyme: Enzyme) -> ThermocyclerProgram | None:
    """Return the supplier's heat inactivation for this enzyme, or ``None`` where it gives none.

    This is not the 60 °C end soak a Golden Gate program ends with, which is a digest.
    """
    celsius, minutes = enzyme.heat_inactivation_celsius, enzyme.heat_inactivation_minutes
    if celsius is None or minutes is None:
        return None
    return ThermocyclerProgram(
        (Stage((Incubation("Heat inactivation", float(celsius), minutes * 60),)),),
        title=f"{enzyme.name} heat inactivation",
    )

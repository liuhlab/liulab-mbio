"""One reaction table, filled to volume, however a pipeline's method mixes it.

The rule is the same wherever it is used: the DNA goes in each tube and everything else into the
master mix, what the components take is summed, DNA that cannot fit the reaction is refused with
the fix named, and one last line fills the rest.
"""

from collections.abc import Sequence

from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.steps import listed
from liulab_mbio.protocol.model import Component, ReactionTable

#: What fills a reaction to volume, where the method names nothing of its own to go in with it.
WATER = "Nuclease-free water"


def reaction_table(
    amounts: Sequence[Amount] = (),
    components: Sequence[Component] = (),
    *,
    volume_ul: float,
    title: str = "",
    filler: str = WATER,
    reactions: int = 1,
    what: str = "",
) -> ReactionTable:
    """Return one reaction: the DNA per tube, the rest of it, and the line that fills the volume.

    Parameters
    ----------
    amounts
        The DNA, one row each, carrying its picomoles and its weight. Each goes in its own tube
        rather than into the master mix.
    components
        Everything else, in pipetting order.
    volume_ul
        What the reaction is made up to.
    title
        The table's.
    filler
        What the last line is called, where a method's own buffer goes in with the water.
    reactions
        How many reactions the table is scaled to.
    what
        What a refusal calls this reaction; the title by default.

    Raises
    ------
    ValueError
        If the DNA and the components do not fit `volume_ul`.

    Examples
    --------
    >>> from liulab_mbio.bench.amounts import dna_amount
    >>> table = reaction_table((dna_amount("pUC19", 2686, pmol=0.05),), volume_ul=15.0)
    >>> [(one.name, one.volume_ul) for one in table.components]
    [('pUC19', 1.0), ('Nuclease-free water', 14.0)]
    """
    taken_ul = sum(component.volume_ul for component in components)
    fits(amounts, volume_ul=volume_ul, taken_ul=taken_ul, what=what or title)
    used = sum(amount.volume_ul for amount in amounts) + taken_ul
    rows = (
        *dna_components(amounts),
        *components,
        Component(filler, round(volume_ul - used, 2), final=f"to {volume_ul:g} µL"),
    )
    return ReactionTable(rows, title=title, reactions=reactions)


def dna_components(amounts: Sequence[Amount]) -> tuple[Component, ...]:
    """Return one row a DNA, each added to its own tube rather than to the master mix."""
    return tuple(
        Component(
            amount.name,
            amount.volume_ul,
            final=f"{amount.pmol:g} pmol ({amount.nanograms:g} ng)",
            master_mix=False,
        )
        for amount in amounts
    )


def fits(
    amounts: Sequence[Amount], *, volume_ul: float, taken_ul: float = 0.0, what: str = ""
) -> None:
    """Refuse DNA that does not fit its reaction, naming what to concentrate or to scale.

    `taken_ul` is what everything else takes, which is room the DNA does not have. A caller
    building the table hands the rule its components; one measuring the DNA before there is a
    table, as `liulab_mbio.library.bench` does, hands it the volume they will take.

    Raises
    ------
    ValueError
        If the DNA and `taken_ul` reach `volume_ul`, leaving nothing to fill the rest with.
    """
    used = sum(amount.volume_ul for amount in amounts) + taken_ul
    if used < volume_ul:
        return
    fix = "scale the reaction up"
    room = volume_ul - taken_ul
    if room > 0 and amounts:
        needed = sum(amount.nanograms for amount in amounts) / room
        names = listed([amount.name for amount in amounts])
        fix = f"concentrate {names} to {needed:.3g} ng/µL or more, or {fix}"
    raise ValueError(
        f"{what}: {used:g} µL of components exceeds the {volume_ul:g} µL reaction; {fix}"
    )

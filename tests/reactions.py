"""A reaction table read the way the bench tests check one."""

from liulab_mbio.protocol import ReactionTable


def volumes(table: ReactionTable) -> dict[str, float]:
    """Each component's volume in one reaction, by name."""
    return {component.name: component.volume_ul for component in table.components}


def total(table: ReactionTable) -> float:
    """The volume of one reaction."""
    return sum(component.volume_ul for component in table.components)

"""A reaction: one tube, the molecules in it, and the enzymes acting on them.

A reaction is a digest, a ligation, or a one-pot that does both. It binds each role to the
molecules that fill it, names the enzymes acting there, and carries how many vessels differ only
in which species fills a role — one reaction over a 96-well array is one reaction, not 96.

It holds no verdict and no intent. What must be true inside one is the method's to judge, which
is why a gate reads reactions rather than bare records: the composition says which tube a
failure belongs to without saying what the designer was thinking. One reaction's product is the
next one's input, and a sequence of them records that order and no constraint of its own.
"""

from collections.abc import Iterable, Iterator, Sequence
from dataclasses import KW_ONLY, dataclass
from typing import Literal

from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.sites import EnzymeLike

#: What one tube does.
type Kind = Literal["digest", "ligation", "one-pot"]

#: What a molecule is in the tube for: the thing opened, the thing cut out of its block, the
#: thing going in, or the thing holding a part until it is wanted.
type Role = Literal["destination", "donor", "insert", "carrier"]


@dataclass(frozen=True, slots=True)
class Pool:
    """The molecules filling one role of a reaction: one species, or many.

    Parameters
    ----------
    role
        What these molecules are in the tube for.
    records
        The species, in the order the design names them. A pool holds at least one.

    Raises
    ------
    ValueError
        If no molecule fills the role.
    """

    role: Role
    records: Sequence[SequenceRecord]

    def __post_init__(self) -> None:
        """Fix the order, and refuse a role nothing fills."""
        object.__setattr__(self, "records", tuple(self.records))
        if not self.records:
            raise ValueError(f"the {self.role} role of this reaction holds no molecule")

    @property
    def one(self) -> SequenceRecord:
        """The single species filling this role.

        Raises
        ------
        ValueError
            If the pool holds more than one, where the caller expected a single molecule.
        """
        if len(self.records) != 1:
            raise ValueError(
                f"the {self.role} role holds {len(self.records)} species, and one was asked for"
            )
        return self.records[0]

    def __len__(self) -> int:
        """How many species fill the role."""
        return len(self.records)

    def __iter__(self) -> Iterator[SequenceRecord]:
        """Read the species in the order the design names them."""
        return iter(self.records)


@dataclass(frozen=True, slots=True)
class Reaction:
    """One tube: what is in it, what acts on it, and how many vessels it runs over.

    Parameters
    ----------
    name
        What the tube is called, so a report says which one a failure belongs to.
    kind
        ``"digest"``, ``"ligation"``, or ``"one-pot"`` for a tube that does both.
    pools
        One per role, each holding one species or many. At least one, and no role twice.
    enzymes
        The enzymes acting in this tube, each an `Enzyme` or a name
        `liulab_mbio.enzymes.get_enzyme` answers to. A ligation a ligase alone drives names
        none.
    vessels
        How many tubes differ only in which species fills a role, one and up.

    Raises
    ------
    ValueError
        If no pool is given, a role is bound twice, or `vessels` is under one.
    KeyError
        If a named enzyme is not one this package ships.
    """

    name: str
    kind: Kind
    _: KW_ONLY
    pools: Sequence[Pool]
    enzymes: Iterable[EnzymeLike] = ()
    vessels: int = 1

    def __post_init__(self) -> None:
        """Fix the order of the pools and the enzymes, then refuse a tube nothing could be."""
        object.__setattr__(self, "pools", tuple(self.pools))
        object.__setattr__(self, "enzymes", tuple(self.enzymes))
        if not self.pools:
            raise ValueError(f"the reaction {self.name!r} holds no molecule")
        if len(set(self.roles)) != len(self.pools):
            named = ", ".join(sorted(self.roles))
            raise ValueError(
                f"the reaction {self.name!r} binds a role twice, over {named}: one pool a role"
            )
        if self.vessels < 1:
            raise ValueError(
                f"the reaction {self.name!r} runs over {self.vessels} vessel(s), and a reaction "
                "is at least one tube"
            )
        for named in self.enzymes:
            if not isinstance(named, Enzyme):
                get_enzyme(named)

    @property
    def acting(self) -> tuple[Enzyme, ...]:
        """The acting enzymes, each read into a record, in the order this reaction names them."""
        return tuple(one if isinstance(one, Enzyme) else get_enzyme(one) for one in self.enzymes)

    @property
    def roles(self) -> tuple[Role, ...]:
        """The roles this reaction binds, in the order its pools are given."""
        return tuple(pool.role for pool in self.pools)

    @property
    def records(self) -> tuple[SequenceRecord, ...]:
        """Every molecule in the tube, pool by pool."""
        return tuple(record for pool in self.pools for record in pool)

    def __getitem__(self, role: Role) -> Pool:
        """Return the pool filling that role.

        Raises
        ------
        KeyError
            If this reaction binds no such role.
        """
        for pool in self.pools:
            if pool.role == role:
                return pool
        raise KeyError(role)

    def __contains__(self, role: object) -> bool:
        """Whether this reaction binds that role."""
        return role in self.roles

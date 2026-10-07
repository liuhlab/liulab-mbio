"""The oligo pool a library's blocks are synthesised from, and the primers that amplify it.

A block longer than one oligo is ordered as several and assembled in one pot before it is a
part. This is where that happens for this method: each block's cargo is split by
`liulab_mbio.split`, each fragment is dressed as an oligo by `liulab_mbio.bench.pools`, and the
method's own choices are what this module supplies -- which enzyme cuts the oligo, which two
overhangs are held out by name, and which primer of the orthogonal set serves which role.

**The cargo is split, not the whole block.** A block's external stuffers are the flanks its
destination already carries, and they are each other's reverse complement but for the overhang:
synthesising them would give every block two ends that anneal to each other and nothing to
assemble into. `cargo_record` cuts them off, so the outermost overhangs are the ones the part
enters and leaves its destination on.

The two reserved overhangs are the method's own interface, read off its stuffers rather than
stated: a part enters the vector on one and leaves its cloning scar as the other, so an internal
junction spelling either would ligate a fragment straight into the vector.

**The pool's count is reported twice on purpose.** `PoolPlan.floor` is the arithmetic: the
fewest oligos the length budget could ever need. `Pool.count` is the design: what the split
actually spends, which is one more wherever a block's forced cuts spell no legal overhang. The
two differing is a fact about a block, not a defect, and the report names which block.
"""

import os
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import KW_ONLY, dataclass
from decimal import Decimal
from math import ceil
from pathlib import Path

from liulab_mbio.bench.pools import (
    Oligo,
    OligoLayout,
    Pool,
    PrimerSite,
    build_oligo,
    pool_item,
)
from liulab_mbio.bench.prices import Band, Item
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.split import CargoSplit, Fragment, fewest_pieces, split_cargo
from liulab_synbio.igga.method import (
    LUND_SUCCESS,
    ORTHOGONAL_SPLIT,
    SYNTHESIS_ENZYME,
    Scheme,
)
from liulab_synbio.igga.parts import Part
from liulab_synbio.igga.project import Project

#: How long one primer of the orthogonal set is.
PRIMER_LENGTH = 20

#: What the pool is priced and banded by: how many oligos, and how long each one is.
QUANTITIES = ("count", "length")


@dataclass(frozen=True, slots=True)
class Batch:
    """One PCR1 tube, and the plate of PCR2 reactions it feeds.

    Parameters
    ----------
    number
        Which batch it is, counting from one.
    forward, outer
        The pair PCR1 pulls the whole batch out of the pool by, named as the inventory names
        them.
    blocks
        The blocks in it, in the order their inner primers were allotted, which is the order
        they sit in the PCR2 plate.
    """

    number: int
    _: KW_ONLY
    forward: str
    outer: str
    blocks: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PoolPlan:
    """One library's oligo pool, what it was split from, and what it is banded by.

    Parameters
    ----------
    pool
        Every oligo to order, with the primers that pull each one out.
    splits
        One per block, in the order the parts were given.
    floor
        The fewest oligos the length budget could need, summed over the blocks. The design
        spends this or more.
    item
        The pool as one line of a bill: its count, its bands and its headroom. The money is
        `liulab_mbio.bench.prices.bill`'s, and a hole wherever nobody holds a tariff.
    """

    pool: Pool
    _: KW_ONLY
    splits: tuple[CargoSplit, ...]
    floor: int
    item: Item

    @property
    def over_floor(self) -> tuple[str, ...]:
        """Each block the design spends an oligo more on than the arithmetic allows for."""
        return tuple(
            str(split.cargo.name)
            for split in self.splits
            if split.pieces > fewest_pieces(len(split.cargo.sequence), split.budget)
        )

    @property
    def batches(self) -> tuple[Batch, ...]:
        """Which blocks share one PCR1, and the pair that pulls them out of the pool together.

        Read back off the oligos rather than worked out again: the forward and outer pair every
        one of a block's oligos carries is the pair PCR1 amplifies that block's batch by, so a
        batch here and an oligo on the sheet cannot disagree.
        """
        found: dict[tuple[str, str], list[str]] = {}
        for oligo in self.pool.oligos:
            forward, _, outer = oligo.primers
            blocks = found.setdefault((forward, outer), [])
            if oligo.source not in blocks:
                blocks.append(oligo.source)
        return tuple(
            Batch(number, forward=forward, outer=outer, blocks=tuple(blocks))
            for number, ((forward, outer), blocks) in enumerate(found.items(), 1)
        )

    @property
    def inner_pairs(self) -> tuple[tuple[str, str], ...]:
        """Each pair PCR2 pulls one block out of its batch by: the batch forward, the block inner."""
        return tuple(
            dict.fromkeys((oligo.primers[0], oligo.primers[1]) for oligo in self.pool.oligos)
        )

    def against_lund(self) -> tuple[tuple[int, int, float | None], ...]:
        """How many blocks take each fragment count, beside the share Lund saw clone perfectly.

        A count Lund did not measure carries ``None``: nothing here interpolates one.
        """
        counted = self.pool.fragment_counts()
        return tuple(
            (pieces, blocks, LUND_SUCCESS.get(pieces)) for pieces, blocks in counted.items()
        )


def read_primers(path: str | os.PathLike[str]) -> tuple[PrimerSite, ...]:
    """Read the orthogonal primer set from a two-column sheet the user holds.

    The columns are ``name`` and ``sequence``, tab separated, one primer a row, in the order the
    roles are allotted. `ORTHOGONAL_SPLIT` says how many each role takes.

    Raises
    ------
    ValueError
        If the file is not that sheet, or holds too few primers for the split.
    """
    rows = Path(path).read_text(encoding="utf-8-sig").splitlines()
    if not rows or rows[0].split("\t")[:2] != ["name", "sequence"]:
        raise ValueError(
            f"{os.fspath(path)} is not a primer sheet; expected tab-separated columns "
            "name, sequence"
        )
    found: list[PrimerSite] = []
    allotted = _roles()
    for line, row in enumerate(rows[1:], 2):
        if not row.strip():
            continue
        cells = row.split("\t")
        if len(cells) < 2 or not cells[1].strip():
            raise ValueError(f"line {line} of {os.fspath(path)} names a primer with no sequence")
        place = len(found)
        if place >= len(allotted):
            break
        found.append(PrimerSite(allotted[place], cells[0].strip(), cells[1].strip()))
    if len(found) < len(allotted):
        raise ValueError(
            f"{os.fspath(path)} holds {len(found)} primers and this method allots "
            f"{len(allotted)} across its {len(ORTHOGONAL_SPLIT)} roles"
        )
    return tuple(found)


def design_pool(
    parts: Sequence[Part],
    project: Project,
    *,
    primers: Sequence[PrimerSite],
    bands: Mapping[str, Sequence[Band]] | None = None,
    key: str = "oligo-pool",
    seed: int = 0,
) -> PoolPlan:
    """Split the cargo of every block of `parts`, each fragment one oligo of one pool.

    A batch is one PCR1 and one plate of PCR2, so a gene takes its own inner primer from the
    batch it sits in and the batch takes a forward and an outer primer of its own. The batches
    are equal, which is what the evenness measurement behind `Project.batch_size` argues for.

    Parameters
    ----------
    parts
        Every block to synthesise, as `liulab_synbio.igga.parts.design_parts` wrote them.
    project
        What this build chose: the oligo length, the batch size and the reserved enzymes.
    primers
        The orthogonal set, already in role order; `read_primers` reads one.
    bands
        The vendor's bands for each of `QUANTITIES`, which the headroom is measured against.
        A quantity with no band carries none, and nothing is estimated.
    key
        What a price record prices the pool by.
    seed
        What the filler is drawn from, so a second run writes the same oligos.

    Raises
    ------
    ValueError
        If the batch is wider than the plate of inner primers, or the primer set is short of
        what the method allots. Also if a block is empty, the split refuses it, or its oligo
        spells a reserved site -- each of those naming the block.
    """
    _check_batch(project.batch_size)
    scheme = project.scheme
    layout = OligoLayout(
        project.oligo_length,
        enzyme=SYNTHESIS_ENZYME,
        primers=len(ORTHOGONAL_SPLIT),
        primer_length=PRIMER_LENGTH,
    )
    held = (scheme.entry_overhang, scheme.scar_overhang)
    avoid = tuple(one.name for one in project.reserved_enzymes)
    batches = _batches(len(parts), project.batch_size)
    inner, forward, outer = _allot(primers)
    oligos: list[Oligo] = []
    splits: list[CargoSplit] = []
    floor = 0
    used: list[PrimerSite] = []
    for index, part in enumerate(parts):
        batch, place = batches[index]
        if not part.sequence:
            raise ValueError(f"block {part.name!r} is empty, so there is nothing to synthesise")
        record = cargo_record(part, scheme)
        try:
            split = split_cargo(
                record,
                layout.cutter,
                budget=layout.budget,
                reserved=held,
                avoid=avoid,
            )
        except ValueError as refusal:
            raise ValueError(f"block {part.name!r}: {refusal}") from refusal
        splits.append(split)
        floor += fewest_pieces(len(record.sequence), layout.budget)
        head = (forward[batch % len(forward)],)
        tail = (inner[place], outer[batch % len(outer)])
        used.extend((*head, *tail))
        for fragment in split.fragments:
            oligos.append(
                _oligo(
                    split,
                    fragment,
                    source=part.name,
                    layout=layout,
                    forward=head,
                    reverse=tail,
                    avoid=avoid,
                    seed=seed + index,
                )
            )
    pool = Pool(
        project.name,
        layout=layout,
        oligos=tuple(oligos),
        primers=tuple(dict.fromkeys(used)),
    )
    return PoolPlan(
        pool, splits=tuple(splits), floor=floor, item=pool_item(pool, key=key, bands=bands)
    )


def cargo_record(part: Part, scheme: Scheme) -> SequenceRecord:
    """Return the cargo of one block as a record named for its part: what the split carries.

    The cargo is what the external enzyme releases, the overhang the part enters on through the
    cloning scar. The stuffers either side of it are the destination's own bases and are not
    synthesised.
    """
    overhang = scheme.external.overhang_length
    start = len(scheme.external_stuffer_5) - overhang
    end = len(part.sequence) - len(scheme.external_stuffer_3) + overhang
    return SequenceRecord(part.sequence[start:end], name=part.name)


def _check_batch(size: int) -> None:
    """Refuse a batch holding more blocks than the method has inner primers to name them by.

    Raises
    ------
    ValueError
        Naming the batch size, the role it outruns and the rule that fixes it.
    """
    role, plate = ORTHOGONAL_SPLIT[0]
    if size > plate:
        raise ValueError(
            f"batch_size is {size} and this method allots {plate} {role} primers, one to a "
            "block, because one batch is exactly one plate of PCR2"
        )


def _roles() -> tuple[str, ...]:
    """Which role each primer of the set takes, in the order the set lists them."""
    return tuple(role for role, count in ORTHOGONAL_SPLIT for _ in range(count))


def _allot(
    primers: Sequence[PrimerSite],
) -> tuple[Sequence[PrimerSite], Sequence[PrimerSite], Sequence[PrimerSite]]:
    """Cut the set into its three roles, in the order it is allotted in.

    Raises
    ------
    ValueError
        If the set is short of the allotment, which `read_primers` holds a sheet to as well.
    """
    inner, forward, outer = (count for _, count in ORTHOGONAL_SPLIT)
    allotted = inner + forward + outer
    if len(primers) < allotted:
        raise ValueError(
            f"this method allots {allotted} primers across its {len(ORTHOGONAL_SPLIT)} roles "
            f"and {len(primers)} were given"
        )
    return (
        primers[:inner],
        primers[inner : inner + forward],
        primers[inner + forward : inner + forward + outer],
    )


def _oligo(
    split: CargoSplit,
    fragment: Fragment,
    *,
    source: str,
    layout: OligoLayout,
    forward: Sequence[PrimerSite],
    reverse: Sequence[PrimerSite],
    avoid: Sequence[str],
    seed: int,
) -> Oligo:
    """Dress one fragment as an oligo, naming the block and the draw where it is refused.

    A refusal and not a fresh draw: `liulab_mbio.bench.pools.pad_bases` already screens the
    filler against the reserved sites in the context it sits in, so a site reaching the oligo's
    own check is one the block carries on a strand, which no redraw clears.

    Raises
    ------
    ValueError
        Naming the block, which fragment of it, and the seed the filler was drawn from.
    """
    try:
        return build_oligo(
            split,
            fragment,
            name=f"{source}_f{fragment.index + 1}",
            source=source,
            layout=layout,
            forward=forward,
            reverse=reverse,
            avoid=avoid,
            seed=seed,
        )
    except ValueError as refusal:
        raise ValueError(
            f"block {source!r} fragment {fragment.index + 1} of {split.pieces}, filler drawn "
            f"from seed {seed}: {refusal}"
        ) from refusal


def _batches(parts: int, size: int) -> tuple[tuple[int, int], ...]:
    """Which batch each part sits in and where in it, over batches of equal size.

    A larger library divides into equal batches rather than full ones and a remainder, because
    the evenness a pool is judged on compares subpools against each other.
    """
    count = max(1, ceil(parts / size))
    each = ceil(parts / count)
    return tuple((index // each, index % each) for index in range(parts))


def read_bands(stated: Mapping[str, Iterable[str]]) -> Mapping[str, tuple[Band, ...]]:
    """Read the vendor's bands for each quantity, written ``low-high`` with no top left empty.

    Raises
    ------
    ValueError
        If a quantity is not one the pool is banded by, or a band is not written that way.

    Examples
    --------
    >>> read_bands({"length": ["301-350"]})["length"][0].holds(Decimal(350))
    True
    """
    out: dict[str, tuple[Band, ...]] = {}
    for quantity, written in stated.items():
        if quantity not in QUANTITIES:
            raise ValueError(
                f"a pool is banded by {', '.join(QUANTITIES)}, and not by {quantity!r}"
            )
        bands: list[Band] = []
        for one in written:
            low, dash, high = str(one).partition("-")
            if not dash or not low.strip():
                raise ValueError(f"band {one!r} of {quantity} is not written 'low-high'")
            top = high.strip()
            bands.append(Band(quantity, Decimal(low.strip()), Decimal(top) if top else None))
        out[quantity] = tuple(bands)
    return out

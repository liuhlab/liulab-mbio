"""The DMX barcode kit: 96 plasmids in four groups of 24, used as supplied.

The sequences are not shipped. They are read from a copy the user holds, the way
`liulab_mbio.bench.prices` and `liulab_mbio.ligase` read theirs; `read_kit` finds one.
"""

import csv
import os
from collections.abc import Sequence
from dataclasses import KW_ONLY, dataclass
from pathlib import Path

#: How the 96 are grouped. ``docs/research/synthesis-and-assembly-barcode-kit.md``.
GROUPS = 4
GROUP_SIZE = 24

#: The overhangs the four groups chain on, read on the strand the cargo reads on, so the four
#: barcodes assemble in one order and none can be left out.
CHAIN: tuple[str, ...] = ("AGGA", "GTTC", "CCTT", "TCAG", "TTCC")

#: The two universal primers flanking the design, which exist only after barcoding because they
#: come from the barcodes. Group 4's 3' constant and the reverse complement of group 1's 5'.
DMX0 = "TTTGATACGCAGGAAGATGGCCCAC"
DMX7 = "ATCGGTGACGGCGATTCTCACATTT"

#: Where the 96 barcode sequences are read from when a caller names no file. The package ships
#: none; this points at a copy the user holds.
KIT_ENV = "LIULAB_SYNBIO_DMX_BARCODES"

#: The columns one of those files holds, in any order.
KIT_COLUMNS = ("name", "group", "index", "overhang5", "umi", "overhang3", "final_seq")

#: What a file has to hold to be one of these, said in the refusal so a caller need not guess.
KIT_EXPECTED = (
    f"a tab-separated file with the columns {', '.join(KIT_COLUMNS)}, holding "
    f"{GROUPS * GROUP_SIZE} rows: {GROUPS} groups of {GROUP_SIZE}, indexed 1 to {GROUP_SIZE} "
    "within each group, every UMI distinct"
)


@dataclass(frozen=True, slots=True)
class KitBarcode:
    """One member of the kit: which group it belongs to, what it spells, and how it chains.

    Parameters
    ----------
    name
        What the kit calls it, such as ``"DMX_1_1"``.
    group, index
        Which of the `GROUPS` groups it is in, counting from one, and which of `GROUP_SIZE`.
    overhang5, overhang3
        The two overhangs it chains on, which are its group's and not its own.
    umi
        The bases that tell it from the rest of its group.
    sequence
        The whole barcode as it is ordered.
    """

    name: str
    _: KW_ONLY
    group: int
    index: int
    overhang5: str
    overhang3: str
    umi: str
    sequence: str


@dataclass(frozen=True, slots=True)
class Kit:
    """The 96 barcodes a user holds, grouped as the kit groups them.

    Parameters
    ----------
    path
        The file they were read from, so a report can say which copy was used.
    barcodes
        Every barcode, in group and then index order.
    """

    path: Path
    barcodes: tuple[KitBarcode, ...]

    def group(self, group: int) -> tuple[KitBarcode, ...]:
        """Return one group's barcodes, in index order.

        Raises
        ------
        ValueError
            If the kit has no such group.
        """
        if not 1 <= group <= GROUPS:
            raise ValueError(f"the kit has groups 1 to {GROUPS}, and {group} is not one of them")
        return tuple(one for one in self.barcodes if one.group == group)

    def at(self, group: int, index: int) -> KitBarcode:
        """Return the barcode at one group and index, both counting from one.

        Raises
        ------
        ValueError
            If the kit holds no barcode there.
        """
        found = [one for one in self.group(group) if one.index == index]
        if not found:
            raise ValueError(f"the kit holds no barcode {index} of group {group}")
        return found[0]


def read_kit(path: str | os.PathLike[str] | None = None) -> Kit:
    """Read the 96 barcode sequences from a copy the user holds.

    The package ships none. `path` names the file, or `KIT_ENV` does.

    Raises
    ------
    ValueError
        If no file is named, or the file is not one of these, saying what one is.

    Examples
    --------
    >>> read_kit("dmx-barcodes.tsv").barcodes[0].group  # doctest: +SKIP
    1
    """
    named = path if path is not None else os.environ.get(KIT_ENV)
    if not named:
        raise ValueError(
            f"the DMX barcode sequences are not shipped: name the file, or set {KIT_ENV}. "
            f"It is {KIT_EXPECTED}"
        )
    one = Path(named)
    try:
        barcodes = _kit_rows(one.read_text(encoding="utf-8-sig").splitlines())
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            f"{one.name} is not the DMX barcode kit ({error}); expected {KIT_EXPECTED}"
        ) from error
    return Kit(one, barcodes)


def _kit_rows(lines: Sequence[str]) -> tuple[KitBarcode, ...]:
    """Parse the barcode rows, refusing a file that is not the kit.

    Raises
    ------
    ValueError
        Naming what is wrong with it.
    """
    reader = csv.DictReader(lines, delimiter="\t")
    if missing := sorted(set(KIT_COLUMNS) - set(reader.fieldnames or ())):
        raise ValueError(f"missing column(s) {', '.join(missing)}")
    found = [
        KitBarcode(
            row["name"],
            group=int(row["group"]),
            index=int(row["index"]),
            overhang5=row["overhang5"].upper(),
            overhang3=row["overhang3"].upper(),
            umi=row["umi"].upper(),
            sequence=row["final_seq"].upper(),
        )
        for row in reader
    ]
    _check_kit(found)
    return tuple(sorted(found, key=lambda one: (one.group, one.index)))


def _check_kit(found: Sequence[KitBarcode]) -> None:
    """Refuse a kit that is not four groups of 24 with distinct UMIs chaining in one order.

    Raises
    ------
    ValueError
        Naming what is wrong with it.
    """
    if len(found) != GROUPS * GROUP_SIZE:
        raise ValueError(f"{len(found)} rows, not {GROUPS * GROUP_SIZE}")
    if len({one.umi for one in found}) != len(found):
        raise ValueError("two rows share a UMI")
    for group in range(1, GROUPS + 1):
        members = [one for one in found if one.group == group]
        if sorted(one.index for one in members) != list(range(1, GROUP_SIZE + 1)):
            raise ValueError(f"group {group} is not indexed 1 to {GROUP_SIZE}")
        ends = {(one.overhang5, one.overhang3) for one in members}
        if len(ends) != 1:
            raise ValueError(f"group {group} chains on {len(ends)} different overhang pairs")
    # Head to tail, so a barcode cannot land in the wrong position and none can be left out.
    # Checked strand-agnostically: a copy held on the map strand spells the chain reversed and
    # complemented, and the adjacency is the same either way.
    for group in range(1, GROUPS):
        left = next(one for one in found if one.group == group)
        right = next(one for one in found if one.group == group + 1)
        if left.overhang3 != right.overhang5:
            raise ValueError(
                f"group {group} ends on {left.overhang3} and group {group + 1} begins on "
                f"{right.overhang5}, so the four do not chain head to tail"
            )

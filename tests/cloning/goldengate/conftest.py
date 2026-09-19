"""What the Golden Gate tests share: the synthesised inserts, and one many-part plan.

The top-level fixtures make only a two-fragment assembly, so the many-part case brings its own
inserts, from `tests/cloning/goldengate/inserts.py`. `four` puts them in and the command-line
test plans the same case from `insert_files`, so both designs read one product.

The package is imported inside each fixture, as `tests/conftest.py` does.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from .inserts import LINKER, ORIENTATIONS, TAG, fasta

if TYPE_CHECKING:
    from liulab_mbio.cloning.goldengate import Plan
    from liulab_mbio.sequence import SequenceRecord


@pytest.fixture(scope="package")
def linker() -> SequenceRecord:
    """The linker as a record, its one feature drawn in a colour of its own."""
    return _synthesised("Linker", LINKER, "#3366cc")


@pytest.fixture(scope="package")
def tag() -> SequenceRecord:
    """The tag as a record, in another colour."""
    return _synthesised("Tag", TAG, "#cc6633")


@pytest.fixture(scope="package")
def insert_files(tmp_path_factory: pytest.TempPathFactory, gfp_file: Path) -> tuple[Path, ...]:
    """`four`'s three inserts as files, in the order they go round the product."""
    directory = tmp_path_factory.mktemp("inserts")
    written = []
    for name, sequence in (("Linker", LINKER), ("Tag", TAG)):
        path = directory / f"{name.lower()}.fasta"
        path.write_text(fasta(name, sequence), encoding="utf-8")
        written.append(path)
    return (gfp_file, *written)


@pytest.fixture(scope="package")
def four(
    puc19: SequenceRecord, gfp: SequenceRecord, linker: SequenceRecord, tag: SequenceRecord
) -> Plan:
    """A backbone and three inserts, in the order they go round the product.

    The middle one goes in the other way round, so one plan takes both the many-part route and
    the reversed-insert route.
    """
    from liulab_mbio.cloning.goldengate import plan_assembly

    return plan_assembly(puc19, gfp, linker, tag, orientation=ORIENTATIONS)


def _synthesised(name: str, sequence: str, color: str) -> SequenceRecord:
    """An insert written in code, its one feature drawn in a colour of its own."""
    from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand

    return SequenceRecord(
        sequence,
        name=name,
        features=(
            Feature(
                name,
                "misc_feature",
                (Segment(0, len(sequence)),),
                strand=Strand.FORWARD,
                color=color,
            ),
        ),
    )

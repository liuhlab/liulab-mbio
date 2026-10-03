"""What every test may share: the data, the pUC19 and GFP records, one across the origin, a plan.

The package is imported inside each fixture, so a module that fails to import fails the tests
that ask for it rather than every test in the run.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from liulab_mbio.cloning.goldengate import Plan
    from liulab_mbio.sequence import SequenceRecord


@pytest.fixture(scope="session")
def data_dir() -> Path:
    """The sequence files and the example protocol the tests read."""
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def puc19_file(data_dir: Path) -> Path:
    """The pUC19 vector as a SnapGene file."""
    return data_dir / "pUC19.dna"


@pytest.fixture(scope="session")
def gfp_file(data_dir: Path) -> Path:
    """The GFP coding sequence as a SnapGene file."""
    return data_dir / "GFP.dna"


@pytest.fixture(scope="session")
def puc19(puc19_file: Path) -> SequenceRecord:
    """The pUC19 vector, circular, with its multiple cloning site annotated."""
    from liulab_mbio.io import read_record

    return read_record(puc19_file)


@pytest.fixture(scope="session")
def gfp(gfp_file: Path) -> SequenceRecord:
    """The GFP coding sequence, linear."""
    from liulab_mbio.io import read_record

    return read_record(gfp_file)


@pytest.fixture(scope="session")
def across_origin() -> SequenceRecord:
    """100 bases, circular: C at 2..8 and G at 91..98 as a person counts them, A elsewhere.

    Feature ``g`` reads 91..98 then 2..8 across the origin, and ``site`` is one reverse segment
    across it, 97..4, where primer ``p`` binds.
    """
    from liulab_mbio.sequence import (
        BindingSite,
        Feature,
        Primer,
        Segment,
        SequenceRecord,
        Strand,
    )

    return SequenceRecord(
        "A" + "C" * 7 + "A" * 82 + "G" * 8 + "AA",
        topology="circular",
        name="across",
        features=(
            Feature(
                "g",
                "misc_feature",
                (Segment(90, 98), Segment(101, 108)),
                strand=Strand.FORWARD,
                color="#ff0000",
            ),
            Feature(
                "site", "misc_feature", (Segment(96, 104),), strand=Strand.REVERSE, color="#00ff00"
            ),
        ),
        primers=(Primer("p", "GGGTTTCC", binding_sites=(BindingSite(96, 104, Strand.REVERSE),)),),
    )


@pytest.fixture(scope="session")
def plan(puc19: SequenceRecord, gfp: SequenceRecord) -> Plan:
    """GFP into the pUC19 multiple cloning site, every option left at its default.

    A plan and its records are frozen, so a test wanting another builds it, with
    `dataclasses.replace` or `plan_assembly`.
    """
    from liulab_mbio.cloning.goldengate import plan_assembly

    return plan_assembly(puc19, gfp)

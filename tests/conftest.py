"""What every test may share: the data directory, the pUC19 and GFP files and records, and one plan.

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
def plan(puc19: SequenceRecord, gfp: SequenceRecord) -> Plan:
    """GFP into the pUC19 multiple cloning site, every option left at its default.

    A plan and its records are frozen, so a test wanting another builds it, with
    `dataclasses.replace` or `plan_assembly`.
    """
    from liulab_mbio.cloning.goldengate import plan_assembly

    return plan_assembly(puc19, gfp)

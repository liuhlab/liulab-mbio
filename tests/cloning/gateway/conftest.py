"""The fixtures the Gateway tests share: the two records, and one plan over them."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from .records import destination_vector, entry_clone

if TYPE_CHECKING:
    from liulab_mbio.cloning.gateway import Plan
    from liulab_mbio.sequence import SequenceRecord


@pytest.fixture(scope="session")
def entry(gfp: SequenceRecord) -> SequenceRecord:
    """An entry clone carrying the GFP coding sequence between its attL sites."""
    return entry_clone(gfp.sequence)


@pytest.fixture(scope="session")
def destination() -> SequenceRecord:
    """A destination vector: ampicillin resistance, a T7 promoter, and a ccdB cassette."""
    return destination_vector()


@pytest.fixture(scope="session")
def gateway_plan(entry: SequenceRecord, destination: SequenceRecord) -> Plan:
    """That entry clone into that destination vector, every option left at its default."""
    from liulab_mbio.cloning.gateway import plan_gateway

    return plan_gateway(entry, destination)

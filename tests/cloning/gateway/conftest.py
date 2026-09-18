"""The fixtures the Gateway tests share: the two records, and one plan over them."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from .records import attb_insert, destination_vector, donor_vector, entry_clone

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
def insert(gfp: SequenceRecord) -> SequenceRecord:
    """A linear attB-flanked fragment carrying the GFP coding sequence."""
    return attb_insert(gfp.sequence)


@pytest.fixture(scope="session")
def donor() -> SequenceRecord:
    """A donor vector: kanamycin resistance and a ccdB cassette between its attP sites."""
    return donor_vector()


@pytest.fixture(scope="session")
def gateway_plan(entry: SequenceRecord, destination: SequenceRecord) -> Plan:
    """That entry clone into that destination vector, every option left at its default."""
    from liulab_mbio.cloning.gateway import plan_gateway

    return plan_gateway(entry, destination)


@pytest.fixture(scope="session")
def staged_plan(insert: SequenceRecord, destination: SequenceRecord, donor: SequenceRecord) -> Plan:
    """That attB insert through BP into the donor, then LR into that destination vector."""
    from liulab_mbio.cloning.gateway import plan_gateway

    return plan_gateway(insert, destination, donor=donor)

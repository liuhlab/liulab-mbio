"""The record every restriction and ligation test reads: the plasmid the insert is cut out of."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from .records import carrying

if TYPE_CHECKING:
    from liulab_mbio.sequence import SequenceRecord


@pytest.fixture(scope="session")
def source(puc19: SequenceRecord, gfp: SequenceRecord) -> SequenceRecord:
    """The plasmid GFP is cut out of: pUC19 carrying it between its own EcoRI and BamHI sites.

    What a lab member subcloning out of one plasmid into another holds.
    """
    return carrying(puc19, gfp, "pTrc-GFP")

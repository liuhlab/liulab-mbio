"""What the destination's rebuild reads off its two parents.

The rebuild itself runs under `pixi run examples-check`, which is a step of the `docs` CI job
rather than of the gate. These are the parents' side of that promise: the two tracked GenBank
files carry the bases, the topology and the one coding sequence each that
`scripts/build_dmx_vector.py` takes, so a re-conversion that loses any of it fails here rather
than in the example.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from liulab_mbio.sequence import Feature, SequenceRecord


def _cds(record: SequenceRecord, name: str) -> Feature:
    found = [one for one in record.features if one.name == name and one.type == "CDS"]
    assert len(found) == 1, f"{name}: {len(found)} coding sequences"
    return found[0]


def test_dmx0001_carries_the_ampr_the_rebuild_replaces(dmx0001: SequenceRecord) -> None:
    """The parent, and the marker the rebuild takes out: a joined CDS on the reverse strand."""
    from liulab_mbio.sequence import Strand

    assert (len(dmx0001), dmx0001.topology) == (5839, "circular")
    ampr = _cds(dmx0001, "AmpR")
    assert ampr.strand is Strand.REVERSE
    assert [(one.start, one.end) for one in ampr.segments] == [(1524, 2316), (2316, 2385)]


def test_pcr_blunt_ii_topo_carries_the_marker_that_goes_in(
    pcr_blunt_ii_topo: SequenceRecord,
) -> None:
    """The carrier, and the coding sequence the rebuild lifts out of it."""
    from liulab_mbio.sequence import Strand

    assert (len(pcr_blunt_ii_topo), pcr_blunt_ii_topo.topology) == (3519, "circular")
    marker = _cds(pcr_blunt_ii_topo, "NeoR/KanR")
    assert marker.strand is Strand.FORWARD
    assert [(one.start, one.end) for one in marker.segments] == [(1236, 2031)]

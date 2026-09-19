"""The records the restriction and ligation tests are built from.

`tests/data/pUC19.dna` is the vector and `tests/data/GFP.dna` the insert, so every plan here
designs on sequences the rest of the suite designs on too. Nothing about them is hard-coded in
the package; the numbers below are read off the records and pinned here.
"""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from liulab_mbio.sequence import SequenceRecord

#: Where pUC19's own EcoRI and BamHI sites begin, and the bases the two cuts leave between them.
ECORI, BAMHI, STUFFER = 395, 416, 21


def plasmid(sequence: str, name: str) -> SequenceRecord:
    """A circular record written in code, for a case that needs no real plasmid."""
    from liulab_mbio.sequence import SequenceRecord

    return SequenceRecord(sequence, topology="circular", name=name)


def carrying(vector: SequenceRecord, insert: SequenceRecord, name: str) -> SequenceRecord:
    """The vector with `insert` put between its own EcoRI and BamHI sites, both sites kept.

    What the insert annotates comes with it, as it would in a plasmid someone built.
    """
    from liulab_mbio.edits import carried, replace
    from liulab_mbio.sites import find_sites

    first, second = sorted(find_sites(vector, ["EcoRI", "BamHI"]), key=lambda site: site.start)
    built, _ = replace(vector, first.end, second.start, insert.sequence)
    features, primers = carried(insert, 0, len(insert), offset=first.end)
    return dataclasses.replace(
        built,
        name=name,
        features=(*built.features, *features),
        primers=(*built.primers, *primers),
    )


def padded(record: SequenceRecord, bases: str) -> SequenceRecord:
    """The record with `bases` added at each end, what it annotates carried along."""
    from liulab_mbio.edits import insert as added

    one, _ = added(record, len(record), bases)
    one, _ = added(one, 0, bases)
    return one

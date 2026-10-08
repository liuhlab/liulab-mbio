"""DMX: reading a construct back one well at a time, many wells at once.

A method of its own, and it takes any cargo. A design sits one per well, the well is marked,
sequenced and called on its own, and identity stays with well position throughout. What
`liulab_synbio.igga` builds is one cargo among others: it chains this protocol as any caller
would, and nothing here is a gate on a pool nothing re-identifies per member.

Two routes mark a well and one judgement reads them. The barcode ligation route ligates four
DMX barcodes into the construct in lysate; the index PCR route amplifies each well with one
barcoded primer pair. The picking, the pass rule and the reformat are shared; the marking step,
the plate and the depth floor are the route's own. The two floors are not a strict and a lenient
setting of one scale: one is where consensus calling starts, the other where a reader stops
trusting a well, over different amplification and different read filters.
`docs/research/route-choice.md` section 4 compares them.

**A well's marks are arithmetic, not a recorded draw.** The address factorises across the
barcode axes and one axis is the plate, so two plates on one flow cell are told apart by
construction and a demultiplexer can check an address instead of trusting a file.

`method` holds the routes, the addresses, the plates and the protocol's own pieces; `steps`
writes them up as the steps a reader works through; `kit` holds the 96 barcodes, which are not
shipped and are read from a copy the user holds; `carrier` puts one part a well in its own
plasmid, which is where this method's cargo starts. What a call is worth once it is back is
`liulab_mbio.bench.readback`'s, because any method reading a construct back per well judges it
the same way.
"""

from liulab_synbio.dmx.kit import CHAIN, GROUPS, KIT_ENV, Kit, KitBarcode, read_kit
from liulab_synbio.dmx.method import (
    REFERENCES,
    ROUTE_INDEX_PCR,
    ROUTE_LIGATION,
    ROUTES,
    SOURCES,
    Design,
    Route,
    Validation,
    depth_check,
    judge_well,
    refuse_unclonal,
    validation,
    validation_equipment,
    validation_materials,
)
from liulab_synbio.dmx.steps import validation_steps

__all__ = [
    "CHAIN",
    "GROUPS",
    "KIT_ENV",
    "REFERENCES",
    "ROUTES",
    "ROUTE_INDEX_PCR",
    "ROUTE_LIGATION",
    "SOURCES",
    "Design",
    "Kit",
    "KitBarcode",
    "Route",
    "Validation",
    "depth_check",
    "judge_well",
    "read_kit",
    "refuse_unclonal",
    "validation",
    "validation_equipment",
    "validation_materials",
    "validation_steps",
]

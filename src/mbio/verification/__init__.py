"""Holding a clone's sequencing results against the record it should be.

Every public name in its modules imports from here as well:

- `result`: `SequencingResult`, one sample's called bases and what the route gave to trust
  them by. Each reader converts its file into one at its own boundary.
- `align`: laying one result on a record, on either strand and across a circular record's
  origin.
- `trace`: `read_trace`, a Sanger trace read into a sequencing result; `read_channels`, its
  four channels and the scan at each base's peak, for a page to draw.
- `judge`: `verify`, which gives each region the record must read true a verdict, and lists
  every disagreement and where each result landed; `regions`, which reads those regions off
  the junctions a plan tagged with `JUNCTION_TAG`.
- `page`: `draw_page`, which lays one clone's verification out as one page to write.
"""

from mbio.verification.align import Aligned, Column, place
from mbio.verification.judge import (
    CONSENSUS_CAVEAT,
    DEPTH_FLOOR,
    DEPTH_WANTED,
    JUNCTION_TAG,
    MISREAD_SHARE,
    QUALITY_FLOOR,
    Disagreement,
    Kind,
    Placement,
    Verification,
    regions,
    verify,
)
from mbio.verification.page import FLANK, CloseUp, Page, TableRow, draw_page
from mbio.verification.result import SequencingResult
from mbio.verification.trace import MIXED_SHARE, Channels, read_channels, read_trace

__all__ = [
    "CONSENSUS_CAVEAT",
    "DEPTH_FLOOR",
    "DEPTH_WANTED",
    "FLANK",
    "JUNCTION_TAG",
    "MISREAD_SHARE",
    "MIXED_SHARE",
    "QUALITY_FLOOR",
    "Aligned",
    "Channels",
    "CloseUp",
    "Column",
    "Disagreement",
    "Kind",
    "Page",
    "Placement",
    "SequencingResult",
    "TableRow",
    "Verification",
    "draw_page",
    "place",
    "read_channels",
    "read_trace",
    "regions",
    "verify",
]

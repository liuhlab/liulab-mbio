"""Assembling a barcoded library in rounds: the scheme, what it covers, and what it takes.

`scheme` holds the constant architecture of one build, validated as it is read. `coverage` counts
the distinct constructs a design yields and the colonies each round needs for the coverage asked
for. `bench` turns the method's own volumes and masses into amounts at the design's real lengths.
Import a module by path; the scheme's own types are re-exported here.
"""

from liulab_mbio.library.scheme import Position, Scheme, read_scheme

__all__ = ["Position", "Scheme", "read_scheme"]

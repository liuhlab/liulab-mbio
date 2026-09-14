"""A barcoded combinatorial library, assembled by iterative Golden Gate one round at a time.

`scheme` holds the architecture a build is given — the positions, the enzymes, the stuffers, the
cloning scar and the barcode length — read from JSON the user supplies and checked as it is read.
A round opens the library built so far with the internal enzyme, releases that round's part list
with the external enzyme, and ligates the two; ``docs/adr/0004-library-rounds.md`` says why that
is a second pipeline rather than a flag on `liulab_mbio.goldengate`.
"""

from liulab_mbio.library.scheme import Position, Scheme, read_scheme

__all__ = [
    "Position",
    "Scheme",
    "read_scheme",
]

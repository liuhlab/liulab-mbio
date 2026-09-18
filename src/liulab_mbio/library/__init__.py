"""A barcoded combinatorial library, assembled one round at a time.

`plan_library` is the way in: it takes lists of proteins or coding DNA, a scheme and a
destination vector, and runs the whole design. `LibraryPlan.write` puts the synthesis order
sheet, the barcode table, the amino-acid change table, a record for each round, the assembled
product, the protocol as data and the page rendered from it in one directory.

The submodules are the steps it is made of -- `scheme`, `standard`, `parts`, `vector`, `rounds`,
`coverage`, `bench` and `steps`. `scheme` holds the architecture a build is given, validated as
it is read; `standard` chooses the overhang set by what it costs the proteins; `parts` writes
each part's synthesis sequence; `vector` accepts or retrofits the destination; `rounds` simulates
each round; `coverage` counts what a round must sample; `bench` turns the method's volumes and
masses into amounts at the design's real lengths, and `steps` the protocol's own steps.

A library is assembled in rounds and not in one pot, because the product keeps the internal
enzyme's sites and that is what lets the next round open it -- `docs/adr/0004-library-rounds.md`
says why `liulab_mbio.cloning.goldengate` does not generalise to one. The scheme is data the user
supplies, per `docs/adr/0005-scheme-data.md`.

Only the way in and the scheme's own types are re-exported. Everything else is imported by
module, as `liulab_mbio.cloning.goldengate.design` is.
"""

from liulab_mbio.library.plan import Files, LibraryPlan, plan_library
from liulab_mbio.library.scheme import Position, Scheme, read_scheme

__all__ = [
    "Files",
    "LibraryPlan",
    "Position",
    "Scheme",
    "plan_library",
    "read_scheme",
]

"""One module per protocol of a library run, each owning everything its own page prints.

`Protocol` is what each module implements, one method a thing the page shows, and `Run` is the
build every one of them is written against.
`synbio.igga.chain` holds the order and the build assembly, and nothing else names a
protocol by its title.
"""

from synbio.igga.protocols.assembly import ASSEMBLY, Assembly
from synbio.igga.protocols.creation import CREATION, Creation
from synbio.igga.protocols.final import FINAL, FinalLigation
from synbio.igga.protocols.ordering import ORDERING, Ordering
from synbio.igga.protocols.primer_plates import PRIMER_PLATES, PrimerPlating
from synbio.igga.protocols.protocol import Protocol
from synbio.igga.protocols.run import Run
from synbio.igga.protocols.seating import SEATING, Seating
from synbio.igga.protocols.validation import VALIDATION, ReadBack

__all__ = [
    "ASSEMBLY",
    "CREATION",
    "FINAL",
    "ORDERING",
    "PRIMER_PLATES",
    "SEATING",
    "VALIDATION",
    "Assembly",
    "Creation",
    "FinalLigation",
    "Ordering",
    "PrimerPlating",
    "Protocol",
    "ReadBack",
    "Run",
    "Seating",
]

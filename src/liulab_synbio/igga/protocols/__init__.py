"""One module per protocol of a library run, each owning everything its own page prints.

A protocol of the chain is a sitting of its own: `Sitting` is what each module implements, one
method a thing the page shows, and `Run` is the build every one of them is written against.
`liulab_synbio.igga.chain` holds the order and the project assembly, and nothing else names a
protocol by its title.
"""

from liulab_synbio.igga.protocols.assembly import ASSEMBLY, Assembly
from liulab_synbio.igga.protocols.creation import CREATION, Creation
from liulab_synbio.igga.protocols.final import FINAL, FinalLigation
from liulab_synbio.igga.protocols.ordering import ORDERING, Ordering
from liulab_synbio.igga.protocols.primer_plates import PRIMER_PLATES, PrimerPlating
from liulab_synbio.igga.protocols.run import Run
from liulab_synbio.igga.protocols.sitting import Sitting
from liulab_synbio.igga.protocols.validation import VALIDATION, ReadBack

__all__ = [
    "ASSEMBLY",
    "CREATION",
    "FINAL",
    "ORDERING",
    "PRIMER_PLATES",
    "VALIDATION",
    "Assembly",
    "Creation",
    "FinalLigation",
    "Ordering",
    "PrimerPlating",
    "ReadBack",
    "Run",
    "Sitting",
]

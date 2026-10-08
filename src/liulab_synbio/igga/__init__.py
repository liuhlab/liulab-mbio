"""A barcoded combinatorial library, assembled one round at a time.

`plan_igga` is the way in: it takes one build, a file naming lists of proteins or coding DNA
and a destination vector, and runs the whole design. `LibraryPlan.write` puts the synthesis order
sheet, the barcode table, the amino-acid change table, a record for each round, the assembled
product, the block vector each position's cargo closes into, the protocol as data and the page
rendered from it in one directory.

The submodules are the steps it is made of -- `method`, `project`, `standard`, `parts`, `cargo`,
`vector`, `rounds`, `reads`, `bench` and `stages` -- and `gate`, which judges a finished design
whoever composed it. `method` holds `IGGA`, the one method, checked when it is imported;
`project` holds what one build chooses, checked as it is read, and keeps that name because
the file it reads is the user's own `project.json`; `gate` judges each reaction's
molecules by the method's rules; `standard` chooses the overhang set by what it costs the
proteins; `parts` writes each part's synthesis sequence and `cargo` the oligo pool it is ordered
as; `vector` accepts or retrofits the destination and holds the working vector's own ccdB
cassette; `rounds` simulates each round and `reads` designs the reads that check what they
built; `bench` turns the method's volumes and masses into amounts at the design's real lengths,
and `stages` says what each stage asks of the bench model. How many colonies a round takes is
`liulab_mbio.bench.coverage`'s.

The run's own pages are `protocols`, one module a protocol, with `chain` holding their order and
`figures` the two drawings only this method needs. `plan` is the way in, and `cli` is the verb
that runs it.

A library is assembled in rounds and not in one pot, because the product keeps the internal
enzyme's sites and that is what lets the next round open it -- `docs/adr/0004-library-rounds.md`
says why `liulab_mbio.cloning.goldengate` does not generalise to one. The method is code and a
build is a file, per `docs/adr/0010-method-in-code.md`.

Only the way in, the method, the build and the gate are re-exported. Everything else is
imported by module, as `liulab_mbio.cloning.goldengate.design` is.
"""

from liulab_synbio.igga.gate import Judgement, Verdict, check_library
from liulab_synbio.igga.method import IGGA, Scheme
from liulab_synbio.igga.plan import Files, LibraryPlan, plan_igga
from liulab_synbio.igga.project import Build, read_build

__all__ = [
    "IGGA",
    "Build",
    "Files",
    "Judgement",
    "LibraryPlan",
    "Scheme",
    "Verdict",
    "check_library",
    "plan_igga",
    "read_build",
]

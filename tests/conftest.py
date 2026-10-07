"""What every test may share: the data, the pUC19 and GFP records, one across the origin, a plan.

The package is imported inside each fixture, so a module that fails to import fails the tests
that ask for it rather than every test in the run.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from liulab_mbio.cloning.goldengate import Plan
    from liulab_mbio.sequence import SequenceRecord


@pytest.fixture(scope="session")
def data_dir() -> Path:
    """The sequence files and the example protocol the tests read."""
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def puc19_file(data_dir: Path) -> Path:
    """The pUC19 vector as a SnapGene file."""
    return data_dir / "pUC19.dna"


@pytest.fixture(scope="session")
def gfp_file(data_dir: Path) -> Path:
    """The GFP coding sequence as a SnapGene file."""
    return data_dir / "GFP.dna"


@pytest.fixture(scope="session")
def puc19(puc19_file: Path) -> SequenceRecord:
    """The pUC19 vector, circular, with its multiple cloning site annotated."""
    from liulab_mbio.io import read_record

    return read_record(puc19_file)


@pytest.fixture(scope="session")
def gfp(gfp_file: Path) -> SequenceRecord:
    """The GFP coding sequence, linear."""
    from liulab_mbio.io import read_record

    return read_record(gfp_file)


@pytest.fixture(scope="session")
def plvx(data_dir: Path) -> SequenceRecord:
    """pLVX-TetOne-Puro-GFP, the Tet-on lentiviral parent, circular.

    Addgene 171123, sequence 336792, with the 32 Addgene annotations.
    `docs/research/working-vector-plvx-tetone.md` holds its provenance and every site counted
    on it, and `reference_docs/synthesis_and_assembly/working-vector/` the files it was built
    from.
    """
    from liulab_mbio.io import read_record

    return read_record(data_dir / "pLVX-TetOne-Puro-GFP.gb")


@pytest.fixture(scope="session")
def dmx0001(data_dir: Path) -> SequenceRecord:
    """DMX0001, the iGGA destination's parent, 5,839 bp circular.

    Addgene 247434, the depositor's own map, converted from SnapGene to GenBank.
    `scripts/build_dmx_vector.py` rebuilds `docs/examples/ap1-library/vector.gb` from it, and
    `docs/research/synthesis-and-assembly.md` holds its provenance.
    """
    from liulab_mbio.io import read_record

    return read_record(data_dir / "dmx0001.gb")


@pytest.fixture(scope="session")
def pcr_blunt_ii_topo(data_dir: Path) -> SequenceRecord:
    """pCR-Blunt II-TOPO, the Zero Blunt TOPO carrier, 3,519 bp circular.

    The Thermo Fisher (Invitrogen) catalogue map, converted from SnapGene to GenBank. The same
    rebuild takes its `NeoR/KanR` coding sequence as the marker that replaces the parent's
    `AmpR`, and the same note holds its provenance.
    """
    from liulab_mbio.io import read_record

    return read_record(data_dir / "pcr-blunt-ii-topo.gb")


@pytest.fixture(scope="session")
def across_origin() -> SequenceRecord:
    """100 bases, circular: C at 2..8 and G at 91..98 as a person counts them, A elsewhere.

    Feature ``g`` reads 91..98 then 2..8 across the origin, and ``site`` is one reverse segment
    across it, 97..4, where primer ``p`` binds.
    """
    from liulab_mbio.sequence import (
        BindingSite,
        Feature,
        Primer,
        Segment,
        SequenceRecord,
        Strand,
    )

    return SequenceRecord(
        "A" + "C" * 7 + "A" * 82 + "G" * 8 + "AA",
        topology="circular",
        name="across",
        features=(
            Feature(
                "g",
                "misc_feature",
                (Segment(90, 98), Segment(101, 108)),
                strand=Strand.FORWARD,
                color="#ff0000",
            ),
            Feature(
                "site", "misc_feature", (Segment(96, 104),), strand=Strand.REVERSE, color="#00ff00"
            ),
        ),
        primers=(Primer("p", "GGGTTTCC", binding_sites=(BindingSite(96, 104, Strand.REVERSE),)),),
    )


@pytest.fixture(scope="session")
def plan(puc19: SequenceRecord, gfp: SequenceRecord) -> Plan:
    """GFP into the pUC19 multiple cloning site, every option left at its default.

    A plan and its records are frozen, so a test wanting another builds it, with
    `dataclasses.replace` or `plan_assembly`.
    """
    from liulab_mbio.cloning.goldengate import plan_assembly

    return plan_assembly(puc19, gfp)

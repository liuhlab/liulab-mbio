"""No protocol a shipped pipeline writes holds two steps under one key.

A key is a step's handle: the page anchors it there and the bench's check mark is kept under it,
so two steps sharing one would share a tick. A page must still render, so the renderer numbers a
repeat rather than refusing it; this is where a duplicate is caught instead, in the pipeline that
wrote it, the day it lands.

Every pipeline is planned on inputs small enough to cost the gate little, since what is asked of
each plan is only the titles and keys of its steps.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from liulab_mbio.cloning.goldengate import Plan
    from liulab_mbio.protocol import Protocol
    from liulab_mbio.sequence import SequenceRecord

#: Three part lists of two members each, which is a whole iGGA build the gate can afford.
PARTS = (
    {"N_a": "MKTAEK", "N_b": "MKTCEK"},
    {"bZIP_a": "WQAFAK", "bZIP_b": "WQAYAK"},
    {"C_a": "MKTHGK", "C_b": "MKTWGK"},
)

#: The orthogonal primer set this repo ships, which is what an iGGA pool is amplified by.
PRIMERS = Path(__file__).parents[1] / "docs" / "examples" / "ap1-library" / "primers.tsv"


def igga_protocols(directory: Path) -> tuple[Protocol, ...]:
    """Every protocol of an iGGA run, from a build written into `directory`.

    The build asks for the primer plates and the read-back, so the run is the whole chain and
    not the shorter one a plainer build writes.
    """
    from liulab_mbio.sequence import SequenceRecord
    from liulab_mbio.snapgene import write_dna
    from liulab_synbio.igga.method import IGGA
    from liulab_synbio.igga.plan import plan_igga
    from liulab_synbio.igga.project import Build, PrimerPlates

    pad = ("TA" * 80)[:80]
    (directory / "parts.fasta").write_text(
        "".join(f">{name} a part\n{seq}\n" for one in PARTS for name, seq in one.items()),
        encoding="utf-8",
    )
    write_dna(
        SequenceRecord(pad + IGGA.internal_stuffer + pad, topology="circular", name="carrier"),
        directory / "carrier.dna",
    )
    build = Build(
        "library",
        positions=("N", "bZIP", "C"),
        parts=directory / "parts.fasta",
        vector=directory / "carrier.dna",
        host="e-coli-k12",
        oligo_length=350,
        batch_size=96,
        completeness=0.99,
        primers=PRIMERS,
        primer_plates=PrimerPlates(nanomoles=25, stock_um=100, working_ul=100),
        validate_from=0,
        route="index PCR",
    )
    return plan_igga(build, parts=PARTS).chain().protocols


@pytest.fixture(scope="module")
def shipped(
    plan: Plan, puc19: SequenceRecord, gfp: SequenceRecord, tmp_path_factory: pytest.TempPathFactory
) -> tuple[tuple[str, Protocol], ...]:
    """Every protocol every shipped pipeline writes, each named by the pipeline that wrote it."""
    from liulab_mbio.cloning.gateway import plan_gateway
    from liulab_mbio.cloning.gibson import plan_gibson
    from liulab_mbio.cloning.restriction import plan_restriction

    from .cloning.gateway.records import destination_vector, entry_clone

    gateway = plan_gateway(entry_clone(gfp.sequence), destination_vector())
    one_each = {
        "goldengate": plan.protocol(),
        "gibson": plan_gibson(puc19, gfp).protocol(),
        "restriction": plan_restriction(puc19, gfp).protocol(),
        "gateway": gateway.protocol(),
    }
    igga = igga_protocols(tmp_path_factory.mktemp("igga"))
    return tuple(one_each.items()) + tuple(("igga", one) for one in igga)


def test_no_shipped_pipeline_writes_two_steps_under_one_key(
    shipped: tuple[tuple[str, Protocol], ...],
) -> None:
    counted = {
        f"{pipeline}: {protocol.title}": Counter(step.key for step in protocol.steps)
        for pipeline, protocol in shipped
    }
    # Every protocol holds steps, so nothing here can pass by having found none to count.
    assert all(counted.values())
    assert not {
        where: sorted(key for key, seen in keys.items() if seen > 1)
        for where, keys in counted.items()
        if any(seen > 1 for seen in keys.values())
    }

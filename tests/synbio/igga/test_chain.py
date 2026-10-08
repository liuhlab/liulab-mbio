"""Which protocols a run writes, in what order, and that each one answers for itself.

A protocol module reads only the part of the run it needs, so these build a `Run` carrying
nothing but the facts the chain's order turns on. What each protocol puts on its page is
covered where a whole plan is built, in `test_plan.py` and `test_ap1_demo.py`.
"""

from types import SimpleNamespace
from typing import Any, cast

from liulab_mbio.sequence import SequenceRecord
from liulab_synbio import dmx
from liulab_synbio.igga.chain import sittings
from liulab_synbio.igga.method import IGGA
from liulab_synbio.igga.protocols import (
    ASSEMBLY,
    CREATION,
    FINAL,
    ORDERING,
    PRIMER_PLATES,
    VALIDATION,
    Assembly,
    Creation,
    FinalLigation,
    Ordering,
    PrimerPlating,
    ReadBack,
    Run,
)


def run(**asked: Any) -> Run:
    """A run carrying the facts the chain's order turns on, and stand-ins for the rest."""
    nothing = cast(Any, None)
    return Run(
        scheme=IGGA,
        positions=("N",),
        barcode_length=12,
        vector=SequenceRecord("TA" * 200, topology="circular", name="pDest"),
        destination=nothing,
        part_lists=(nothing,),
        standard=nothing,
        parts=(),
        rounds=(nothing,),
        bench=(),
        constructs=4,
        checks=(),
        host="e-coli-k12",
        sheet="parts.tsv",
        barcodes="barcodes.tsv",
        **asked,
    )


def titles(one: Run) -> list[str]:
    """What each protocol of this run's chain is headed, in order."""
    return [sitting.title(one) for sitting in sittings(one)]


def test_a_plain_run_orders_its_blocks_assembles_and_moves_the_library() -> None:
    assert titles(run()) == [ORDERING, ASSEMBLY, FINAL]


def test_a_pool_splits_ordering_from_making_the_cargo() -> None:
    assert titles(run(pool=cast(Any, object()))) == [ORDERING, CREATION, ASSEMBLY, FINAL]


def test_plating_the_primers_opens_the_chain_with_its_own_protocol() -> None:
    one = run(pool=cast(Any, object()), primer_plates=cast(Any, object()))

    assert titles(one)[0] == PRIMER_PLATES
    assert isinstance(sittings(one)[0], PrimerPlating)


def test_a_run_reading_designs_back_names_the_route_on_its_own_page() -> None:
    one = run(validation=cast(Any, SimpleNamespace(route=dmx.ROUTE_INDEX_PCR)))

    assert titles(one) == [ORDERING, f"{VALIDATION}: index PCR", ASSEMBLY, FINAL]


def test_every_protocol_of_a_chain_is_its_own_module_and_knows_its_own_title() -> None:
    assert Ordering().title(run()) == ORDERING
    assert Creation().title(run()) == CREATION
    assert Assembly().title(run()) == ASSEMBLY
    assert FinalLigation().title(run()) == FINAL
    assert {type(one) for one in sittings(run())} == {Ordering, Assembly, FinalLigation}


def test_ordering_says_less_where_there_is_no_pool_to_store() -> None:
    assert Ordering().summary(run()) == "Order every block this library is built from."
    assert "put it away" in Ordering().summary(run(pool=cast(Any, object())))


def test_only_the_round_protocols_buy_the_reagents_the_rounds_share() -> None:
    assert [Creation.round_reagents, Assembly.round_reagents, FinalLigation.round_reagents] == [
        True,
        True,
        True,
    ]
    assert [Ordering.round_reagents, PrimerPlating.round_reagents, ReadBack.round_reagents] == [
        False,
        False,
        False,
    ]

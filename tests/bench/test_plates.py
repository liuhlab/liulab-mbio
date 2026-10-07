"""The plate: its formats, where a thing sits, and the moves between wells."""

import pytest

from liulab_mbio.bench import plates
from liulab_mbio.protocol.model import FORMATS, Well


@pytest.mark.parametrize("wells", sorted(FORMATS))
def test_every_format_names_each_of_its_wells_once(wells: int) -> None:
    one = plates.plate("plate", wells)
    assert len(one.well_names) == wells == one.rows * one.columns
    assert len(set(one.well_names)) == wells


def test_the_densest_format_reaches_1536_because_four_384s_compress_into_one() -> None:
    assert plates.plate("compressed", 1536).well_names[-1] == "AF48"


def test_a_plate_refuses_a_format_no_plate_comes_in() -> None:
    with pytest.raises(ValueError, match="100 wells is no format"):
        plates.plate("plate", 100)


def test_a_plate_refuses_a_well_off_its_array() -> None:
    with pytest.raises(ValueError, match="has no well 'Z1'"):
        plates.plate("plate", 96, seating={"Z1": "water"})


def test_seating_fills_the_wells_in_reading_order() -> None:
    assert plates.seat(["a", "b"], 96, start=11) == {"A12": "a", "B1": "b"}


def test_seating_refuses_more_names_than_the_format_holds() -> None:
    with pytest.raises(ValueError, match="run past"):
        plates.seat(["a"] * 13, 12)


def test_compaction_is_dense_and_leaves_out_the_wells_that_failed() -> None:
    picked = plates.plate("picked", 384)
    good = [Well("picked", name) for name in ("A1", "A3", "B2")]
    moved = plates.compact(good, plates.plate("lysate", 1536), 2.0, title="Compress")
    assert [move.destination.well for move in moved.moves] == ["A1", "A2", "A3"]
    assert [move.source.well for move in moved.moves] == ["A1", "A3", "B2"]
    assert picked.wells == 384


def test_compaction_refuses_more_wells_than_the_destination_holds() -> None:
    sources = plates.wells_of(plates.plate("big", 384))
    with pytest.raises(ValueError, match="do not fit"):
        plates.compact(sources, plates.plate("small", 96), 1.0, title="Compress")


def test_pooling_runs_many_wells_into_one() -> None:
    one = plates.plate("lysate", 1536)
    pooled = plates.pool(plates.wells_of(one), Well("reservoir", "1"), 3.5, title="Pool")
    assert len(pooled.moves) == 1536
    assert {move.destination for move in pooled.moves} == {Well("reservoir", "1")}
    assert pooled.plates == ("lysate", "reservoir")

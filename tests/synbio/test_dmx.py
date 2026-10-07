"""The DMX stage: how a well is addressed, how deep its read must be, and what passes."""

import pytest

from liulab_mbio.protocol.model import Well
from liulab_synbio import dmx

KIT = """name\tgroup\tindex\toverhang5\tumi\toverhang3\tfinal_seq
"""


def kit_text(groups: int = dmx.GROUPS, size: int = dmx.GROUP_SIZE) -> str:
    """Return a kit file of `groups` groups of `size`, chaining the way the real one does."""
    chain = ["AGGA", "GTTC", "CCTT", "TCAG", "TTCC"]
    lines = [KIT.strip()]
    for group in range(1, groups + 1):
        for index in range(1, size + 1):
            umi = f"{group:02d}{index:02d}" + "A" * 21
            lines.append(
                f"DMX_{group}_{index}\t{group}\t{index}\t{chain[group - 1]}\t{umi}\t"
                f"{chain[group]}\tGGTCTCA{umi}CGAGACC"
            )
    return "\n".join(lines) + "\n"


def write_kit(tmp_path, text: str):
    """Write `text` as a kit file and return the path."""
    path = tmp_path / "dmx-barcodes.tsv"
    path.write_text(text, encoding="utf-8")
    return path


def test_address_factorises_across_the_axes():
    """A well's marks are its index in mixed radix, and the plate is its own axis."""
    first = dmx.address(dmx.ROUTE_A, plate=0, well=0)
    assert first.marks == (1, 1, 1, 1)
    last = dmx.address(dmx.ROUTE_A, plate=0, well=dmx.ROUTE_A.wells_per_plate - 1)
    assert last.marks == (24, 24, 24, 1)
    assert len({dmx.address(dmx.ROUTE_A, plate=0, well=n).marks for n in range(500)}) == 500


def test_the_plate_axis_tells_two_plates_apart():
    """Two wells at the same position on different plates differ only in the plate's mark."""
    here = dmx.address(dmx.ROUTE_B, plate=0, well=7)
    there = dmx.address(dmx.ROUTE_B, plate=5, well=7)
    assert here.well_marks == there.well_marks
    assert (here.plate_mark, there.plate_mark) == (1, 6)


def test_route_a_addresses_a_compressed_plate_and_route_b_the_levseq_set():
    """Three groups reach past 1,536 wells, and 96 by 96 reaches 9,216."""
    assert dmx.ROUTE_A.wells_per_plate >= dmx.COMPRESSED_WELLS
    assert dmx.ROUTE_B.capacity == 9216


def test_an_address_past_the_route_is_refused():
    """A well or a plate the route cannot tell apart is named rather than wrapped."""
    with pytest.raises(ValueError, match="past the last"):
        dmx.address(dmx.ROUTE_B, plate=0, well=dmx.INDEX_WELLS)
    with pytest.raises(ValueError, match="past the last"):
        dmx.address(dmx.ROUTE_B, plate=dmx.INDEX_WELLS, well=0)


def test_below_the_floor_is_no_verdict_rather_than_a_fail():
    """A well nobody could call carries None, which is what keeps it out of the compaction."""
    assert dmx.depth_check(dmx.ROUTE_A, 151).status == "pass"
    assert dmx.depth_check(dmx.ROUTE_A, 150).status is None
    assert dmx.depth_check(dmx.ROUTE_B, 20).status == "pass"
    assert dmx.depth_check(dmx.ROUTE_B, 10).status == "warn"
    assert dmx.depth_check(dmx.ROUTE_B, 9).status is None


def test_the_floor_travels_with_the_route_and_a_project_may_only_raise_it():
    """Each floor was measured on its own library prep, so neither is the other's."""
    assert dmx.ROUTE_A.wanted_reads != dmx.ROUTE_B.wanted_reads
    assert dmx.depth_check(dmx.ROUTE_B, 25, wanted=30).status == "warn"
    with pytest.raises(ValueError, match="raise"):
        dmx.depth_check(dmx.ROUTE_B, 25, wanted=5)


def test_a_pass_is_an_exact_match_and_a_silent_change_is_not_one():
    """A barcode that no longer names its member cannot be put right by the linkage read."""
    designed = "AGGAATGAAACCGTTCC"
    assert dmx.identity_check((designed,), designed).status == "pass"
    assert dmx.identity_check((designed[:-1] + "A",), designed).status == "fail"
    assert dmx.identity_check((designed, designed[:-1] + "A"), designed).status == "fail"
    assert dmx.identity_check((), designed).status is None


def test_reformatting_compacts_out_failures_and_keeps_the_uncalled():
    """Compacting out a well nobody read throws away a design that may be clean."""
    designed = "AGGAATGTTCC"
    kept = dmx.judge_well(
        Well("picked", "A1"), route=dmx.ROUTE_A, reads=4, called=(), designed=designed
    )
    gone = dmx.judge_well(
        Well("picked", "A2"), route=dmx.ROUTE_A, reads=400, called=("AGGAA",), designed=designed
    )
    passed = dmx.judge_well(
        Well("picked", "A3"), route=dmx.ROUTE_A, reads=400, called=(designed,), designed=designed
    )
    assert (kept.status, gone.status, passed.status) == ("pass", "fail", "pass")
    assert kept.called is False
    assert [one.well.well for one in dmx.reformat((kept, gone, passed))] == ["A1", "A3"]


def test_the_clean_colony_curve_returns_its_measured_anchors():
    """Every anchor is Lund's own; between them the curve is interpolated and says so."""
    for fragments, chance in dmx.CLEAN_COLONY_CURVE:
        assert dmx.clean_colony_chance(fragments) == pytest.approx(chance)
    assert dmx.clean_colony_chance(1) == 1.0
    assert dmx.clean_colony_chance(20) == 0.0
    assert 0.40 < dmx.clean_colony_chance(8) < 0.846


def test_a_kit_the_user_holds_is_read_and_its_chain_checked(tmp_path):
    """The sequences are not package data, and a file that is not the kit is refused by name."""
    kit = dmx.read_kit(write_kit(tmp_path, kit_text()))
    assert len(kit.barcodes) == dmx.GROUPS * dmx.GROUP_SIZE
    assert kit.at(2, 3).name == "DMX_2_3"
    names = [one.name for one in dmx.barcodes_for(kit, dmx.address(dmx.ROUTE_A, plate=0, well=0))]
    assert names == ["DMX_1_1", "DMX_2_1", "DMX_3_1", "DMX_4_1"]


def test_a_file_that_is_not_the_kit_says_what_one_is(tmp_path):
    """A short file, a missing column and a broken chain are each refused with the shape."""
    with pytest.raises(ValueError, match="rows, not 96"):
        dmx.read_kit(write_kit(tmp_path, kit_text(size=2)))
    with pytest.raises(ValueError, match="missing column"):
        dmx.read_kit(write_kit(tmp_path, "name\tgroup\n"))


def test_the_kit_is_not_shipped_and_the_refusal_says_where_to_put_one(monkeypatch):
    """Finding it the way the ligase matrix is found means an env var and no default file."""
    monkeypatch.delenv(dmx.KIT_ENV, raising=False)
    with pytest.raises(ValueError, match=dmx.KIT_ENV):
        dmx.read_kit()


def test_route_b_seats_each_sample_under_the_pair_its_address_names(tmp_path):
    """The plate is derived from the address, so there is no plate map to carry."""
    plate = dmx.index_plate("index", 96, plate=3)
    assert plate.wells == dmx.INDEX_WELLS
    assert plate.seating["A1"] == "forward 1, reverse 4"
    assert plate.seating["H12"] == "forward 96, reverse 4"
    with pytest.raises(ValueError, match="do not fit"):
        dmx.index_plate("index", 97)

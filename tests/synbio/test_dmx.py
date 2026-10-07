"""The DMX stage: how a well is addressed, how deep its read must be, and what passes."""

import pytest

from liulab_mbio.protocol.model import Citation, Well
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
    assert dmx.depth_check(dmx.ROUTE_B, 21).status == "pass"
    assert dmx.depth_check(dmx.ROUTE_B, 11).status == "warn"
    assert dmx.depth_check(dmx.ROUTE_B, 9).status is None


def test_a_well_landing_exactly_on_route_bs_marks_has_not_met_them():
    """LevSeq's SI checklist wants an alignment count above 20, so 20 itself is only tolerated.

    The article reads the same number as a minimum; the SI governs, and the stricter reading
    stands where the two still contest the boundary.
    """
    assert dmx.depth_check(dmx.ROUTE_B, 20).status == "warn"
    assert dmx.depth_check(dmx.ROUTE_B, 10).status is None


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
    assert 0.846 < dmx.clean_colony_chance(4) < 0.938


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


def test_picking_fills_one_quarter_of_the_plate_at_a_time():
    """288 wells give three full index plates, not four part-filled ones."""
    picked = dmx.picked_plate("picked 1", 288)
    assert len(picked.seating) == 288
    assert picked.catalog == dmx.PICKED_CATALOG
    assert set(picked.seating.values()) == {"quarter 1", "quarter 2", "quarter 3"}
    assert picked.seating["A1"] == "quarter 1"
    assert picked.seating["B1"] == "quarter 3"
    assert "B2" not in picked.seating
    with pytest.raises(ValueError, match="do not fit"):
        dmx.picked_plate("picked 1", 385)


def test_route_b_samples_one_quarter_of_a_picked_plate_into_each_index_plate():
    """One well in four lines up under the head, so each pass is one full index plate."""
    picked = dmx.picked_plate("picked 1", 288)
    index = [dmx.index_plate(f"index {n}", 96, plate=n - 1) for n in (1, 2, 3)]
    moves = dmx.sampling(picked, index)
    assert [len(one.moves) for one in moves] == [96, 96, 96]
    assert moves[0].instrument == dmx.MULTICHANNEL
    assert moves[0].moves[0].volume_ul == dmx.SAMPLE_UL
    assert moves[0].moves[1].source.well == "A3"
    assert moves[0].moves[1].destination.well == "A2"
    assert index[0].catalog == dmx.INDEX_CATALOG
    with pytest.raises(ValueError, match="one index plate covers one quarter"):
        dmx.sampling(picked, index[:2])


def sized(route, designs, floor):
    """Return what these designs take to read back, refusing the None a floor reading none gives."""
    one = dmx.validation(route, designs, floor)
    assert one is not None
    return one


def test_a_floor_reads_back_every_design_in_that_many_fragments_or_more():
    """Omitted reads nothing, zero reads every design, and the rest is a comparison."""
    some = (dmx.Design("two", 2), dmx.Design("five", 5), dmx.Design("twelve", 12))
    assert dmx.validated(some, None) == ()
    assert dmx.validated(some, 0) == some
    assert [one.name for one in dmx.validated(some, 5)] == ["five", "twelve"]
    assert dmx.validated(some, 13) == ()
    with pytest.raises(ValueError, match="below none"):
        dmx.validated(some, -1)
    with pytest.raises(ValueError, match="at least one fragment"):
        dmx.Design("none", 0)


def test_the_bench_is_sized_from_the_designs_read_and_not_from_the_design_list():
    """A design the floor leaves out costs no well, so the plates shrink by exactly those wells."""
    some = tuple(dmx.Design(f"d{n}", 1 + n % 4) for n in range(72))
    whole = sized(dmx.ROUTE_B, some, 0)
    assert (whole.wells, len(whole.picked), len(whole.index)) == (288, 1, 3)
    assert [len(one.seating) for one in whole.picked] == [288]
    fewer = sized(dmx.ROUTE_B, some, 4)
    assert len(fewer.designs) == 18
    assert fewer.wells == 72
    assert [len(one.seating) for one in fewer.picked] == [72]
    assert len(fewer.index) == 1
    assert dmx.validation(dmx.ROUTE_B, some, 5) is None
    assert dmx.validation(dmx.ROUTE_B, some, None) is None


def test_route_a_compresses_four_picked_plates_into_one_and_route_b_neither():
    """Each route's plates follow from the shared wells, and neither pours the other's."""
    many = tuple(dmx.Design(f"d{n}", 2) for n in range(400))
    route_a = sized(dmx.ROUTE_A, many, 0)
    assert route_a.wells == 1600
    assert len(route_a.picked) == 5
    assert len(route_a.compressed) == 2
    with pytest.raises(ValueError, match="only route B"):
        assert route_a.index
    route_b = sized(dmx.ROUTE_B, many, 0)
    assert len(route_b.index) == 17
    assert [one.seating["A1"] for one in route_b.index[:2]] == [
        "forward 1, reverse 1",
        "forward 1, reverse 2",
    ]
    with pytest.raises(ValueError, match="only route A"):
        assert route_b.compressed


def test_the_steps_print_each_design_chance_beside_the_floor():
    """The number reads as a choice: the floor is stated and the curve is printed beside it."""
    some = (dmx.Design("two", 2), dmx.Design("eight", 8))
    one = sized(dmx.ROUTE_B, some, 2)
    steps = dmx.validation_steps(one)
    assert [step.title for step in steps][:2] == [
        "Array 2 design(s) and grow",
        "Pick 4 colonies of each design",
    ]
    assert "2 fragment(s) or more" in " ".join(steps[0].notes)
    assert "2 fragment(s): 1 design(s), 100.0% of picks clean" in steps[1].notes
    assert "8 fragment(s): 1 design(s), 66.7% of picks clean" in steps[1].notes


def test_route_b_carries_a_hole_at_the_marks_and_route_a_carries_none():
    """The 192 index sequences are lab stock, and no source gives the Taq stock they amplify on."""
    some = (dmx.Design("one", 2),)
    route_b = dmx.validation_steps(sized(dmx.ROUTE_B, some, 0))
    assert [hole.id for step in route_b for hole in step.holes] == ["B1", "B2"]
    route_a = dmx.validation_steps(sized(dmx.ROUTE_A, some, 0))
    assert [hole.id for step in route_a for hole in step.holes] == []


def test_the_plates_a_pick_fills_name_where_their_numbers_were_read():
    """Both routes pick into these two, so a reader of either page can follow the numbers back."""
    some = (dmx.Design("one", 2),)
    for route in (dmx.ROUTE_A, dmx.ROUTE_B):
        cited = {
            one.name: one.citation
            for one in dmx.validation_materials(sized(route, some, 0))
            if one.citation
        }
        assert cited["25 cm BioAssay plate"] == Citation("Qian SI", "Day 2")
        assert cited[f"{dmx.PICKED_WELLS}-well culture plate"] == Citation("Qian SI", "Day 3")
        assert all(one.source in dmx.SOURCES for one in cited.values())


def test_the_index_pcr_is_one_wells_share_of_levseqs_published_mix():
    """Scaled back to a full plate the table is the SI's own, and the Taq carries no unit count."""
    table = dmx.index_pcr_reaction()
    assert round(sum(one.volume_ul for one in table.components), 2) == dmx.INDEX_PCR_UL
    assert table.mix_volumes(dmx.INDEX_WELLS)[:4] == (144.0, 28.8, 7.2, 57.6)
    taq = next(one for one in table.components if one.name.startswith("Taq"))
    assert (taq.volume_ul, taq.stock, taq.final) == (0.05, "", "")


def test_the_index_pcr_touches_down_before_it_plateaus():
    """Ten cycles half a degree apart, then 25 more: 35 in all, as the SI's two loops spell out."""
    stages = dmx.index_pcr_program().stages
    assert [stage.cycles for stage in stages[1:-2]] == [1] * 10 + [25]
    annealing = [stage.incubations[1].temperature_c for stage in stages[1 : 1 + 10]]
    assert annealing == [68.0, 67.5, 67.0, 66.5, 66.0, 65.5, 65.0, 64.5, 64.0, 63.5]

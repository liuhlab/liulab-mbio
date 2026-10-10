"""DMX: how a well is addressed, how deep its read must be, and what the kit has to be."""

import pytest

from mbio.protocol.model import Citation, Step, Well
from synbio.dmx import kit, method, steps

KIT = """name\tgroup\tindex\toverhang5\tumi\toverhang3\tfinal_seq
"""


def kit_text(groups: int = kit.GROUPS, size: int = kit.GROUP_SIZE) -> str:
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
    first = method.address(method.ROUTE_LIGATION, plate=0, well=0)
    assert first.marks == (1, 1, 1, 1)
    last = method.address(
        method.ROUTE_LIGATION, plate=0, well=method.ROUTE_LIGATION.wells_per_plate - 1
    )
    assert last.marks == (24, 24, 24, 1)
    assert (
        len({method.address(method.ROUTE_LIGATION, plate=0, well=n).marks for n in range(500)})
        == 500
    )


def test_the_plate_axis_tells_two_plates_apart():
    """Two wells at the same position on different plates differ only in the plate's mark."""
    here = method.address(method.ROUTE_INDEX_PCR, plate=0, well=7)
    there = method.address(method.ROUTE_INDEX_PCR, plate=5, well=7)
    assert here.well_marks == there.well_marks
    assert (here.plate_mark, there.plate_mark) == (1, 6)


def test_ligation_addresses_a_compressed_plate_and_index_pcr_the_levseq_set():
    """Three groups reach past 1,536 wells, and 96 by 96 reaches 9,216."""
    assert method.ROUTE_LIGATION.wells_per_plate >= method.COMPRESSED_WELLS
    assert method.ROUTE_INDEX_PCR.capacity == 9216


def test_an_address_past_the_route_is_refused():
    """A well or a plate the route cannot tell apart is named rather than wrapped."""
    with pytest.raises(ValueError, match="past the last"):
        method.address(method.ROUTE_INDEX_PCR, plate=0, well=method.INDEX_WELLS)
    with pytest.raises(ValueError, match="past the last"):
        method.address(method.ROUTE_INDEX_PCR, plate=method.INDEX_WELLS, well=0)


def test_below_the_floor_is_no_verdict_rather_than_a_fail():
    """A well nobody could call carries None, which is what keeps it out of the compaction."""
    assert method.depth_check(method.ROUTE_LIGATION, 151).status == "pass"
    assert method.depth_check(method.ROUTE_LIGATION, 150).status is None
    assert method.depth_check(method.ROUTE_LIGATION, 149).status is None
    assert method.depth_check(method.ROUTE_INDEX_PCR, 21).status == "pass"
    assert method.depth_check(method.ROUTE_INDEX_PCR, 11).status == "warn"
    assert method.depth_check(method.ROUTE_INDEX_PCR, 9).status is None


def test_the_wanted_mark_is_exceeded_and_the_tolerable_one_is_reached():
    """LevSeq's SI wants an alignment count above 20 and puts detection's floor at 10 reads.

    So 20 itself is only tolerated, while 10 itself still carries a verdict.
    """
    assert method.depth_check(method.ROUTE_INDEX_PCR, 20).status == "warn"
    assert method.depth_check(method.ROUTE_INDEX_PCR, 10).status == "warn"


def test_the_floor_travels_with_the_route_and_a_build_may_only_raise_it():
    """Not a strict and a lenient setting of one scale, so neither floor is the other's."""
    assert method.ROUTE_LIGATION.wanted_reads != method.ROUTE_INDEX_PCR.wanted_reads
    assert method.depth_check(method.ROUTE_INDEX_PCR, 25, wanted=30).status == "warn"
    with pytest.raises(ValueError, match="raise"):
        method.depth_check(method.ROUTE_INDEX_PCR, 25, wanted=5)


def test_a_well_too_thin_to_call_is_judged_on_its_depth_alone():
    """The two questions run in order, which is why a shallow well is never a failure."""
    designed = "AGGAATGTTCC"
    thin = method.judge_well(
        Well("picked", "A1"), route=method.ROUTE_LIGATION, reads=4, called=(), designed=designed
    )
    deep = method.judge_well(
        Well("picked", "A2"),
        route=method.ROUTE_LIGATION,
        reads=400,
        called=("AGGAA",),
        designed=designed,
    )
    assert [one.name for one in thin.checks] == ["reads_per_well"]
    assert thin.called is False
    assert [one.name for one in deep.checks] == ["reads_per_well", "well_identity"]
    assert deep.status == "fail"


def test_a_kit_the_user_holds_is_read_and_its_chain_checked(tmp_path):
    """The sequences are not package data, and a file that is not the kit is refused by name."""
    held = kit.read_kit(write_kit(tmp_path, kit_text()))
    assert len(held.barcodes) == kit.GROUPS * kit.GROUP_SIZE
    assert held.at(2, 3).name == "DMX_2_3"
    names = [
        one.name
        for one in method.barcodes_for(held, method.address(method.ROUTE_LIGATION, plate=0, well=0))
    ]
    assert names == ["DMX_1_1", "DMX_2_1", "DMX_3_1", "DMX_4_1"]


def test_a_file_that_is_not_the_kit_says_what_one_is(tmp_path):
    """A short file, a missing column and a broken chain are each refused with the shape."""
    with pytest.raises(ValueError, match="rows, not 96"):
        kit.read_kit(write_kit(tmp_path, kit_text(size=2)))
    with pytest.raises(ValueError, match="missing column"):
        kit.read_kit(write_kit(tmp_path, "name\tgroup\n"))


def test_the_kit_is_not_shipped_and_the_refusal_says_where_to_put_one(monkeypatch):
    """Finding it the way the ligase matrix is found means an env var and no default file."""
    monkeypatch.delenv(kit.KIT_ENV, raising=False)
    with pytest.raises(ValueError, match=kit.KIT_ENV):
        kit.read_kit()


def test_index_pcr_seats_each_sample_under_the_pair_its_address_names(tmp_path):
    """The plate is derived from the address, so there is no plate map to carry."""
    plate = method.index_plate("index", 96, plate=3)
    assert plate.wells == method.INDEX_WELLS
    assert plate.labels["A1"] == "forward 1, reverse 4"
    assert plate.labels["H12"] == "forward 96, reverse 4"
    with pytest.raises(ValueError, match="do not fit"):
        method.index_plate("index", 97)


def test_picking_fills_one_quarter_of_the_plate_at_a_time():
    """288 wells give three full index plates, not four part-filled ones."""
    picked = method.picked_plate("picked 1", 288)
    assert len(picked.labels) == 288
    assert picked.catalog == method.PICKED_CATALOG
    assert set(picked.labels.values()) == {"quarter 1", "quarter 2", "quarter 3"}
    assert picked.labels["A1"] == "quarter 1"
    assert picked.labels["B1"] == "quarter 3"
    assert "B2" not in picked.labels
    with pytest.raises(ValueError, match="do not fit"):
        method.picked_plate("picked 1", 385)


def test_index_pcr_samples_one_quarter_of_a_picked_plate_into_each_index_plate():
    """One well in four lines up under the head, so each pass is one full index plate."""
    picked = method.picked_plate("picked 1", 288)
    index = [method.index_plate(f"index {n}", 96, plate=n - 1) for n in (1, 2, 3)]
    moves = method.sampling(picked, index)
    assert [len(one.moves) for one in moves] == [96, 96, 96]
    assert moves[0].instrument == method.MULTICHANNEL
    assert moves[0].moves[0].volume_ul == method.SAMPLE_UL
    assert moves[0].moves[1].source.well == "A3"
    assert moves[0].moves[1].destination.well == "A2"
    assert index[0].catalog == method.INDEX_CATALOG
    with pytest.raises(ValueError, match="one index plate covers one quarter"):
        method.sampling(picked, index[:2])


def sized(route, designs, floor):
    """Return what these designs take to read back, refusing the None a floor reading none gives."""
    one = method.validation(route, designs, floor)
    assert one is not None
    return one


def test_a_floor_reads_back_every_design_in_that_many_fragments_or_more():
    """Omitted reads nothing, zero reads every design, and the rest is a comparison."""
    some = (method.Design("two", 2), method.Design("five", 5), method.Design("twelve", 12))
    assert method.validated(some, None) == ()
    assert method.validated(some, 0) == some
    assert [one.name for one in method.validated(some, 5)] == ["five", "twelve"]
    assert method.validated(some, 13) == ()
    with pytest.raises(ValueError, match="below none"):
        method.validated(some, -1)
    with pytest.raises(ValueError, match="at least one fragment"):
        method.Design("none", 0)


def test_the_bench_is_sized_from_the_designs_read_and_not_from_the_design_list():
    """A design the floor leaves out costs no well, so the plates shrink by exactly those wells."""
    some = tuple(method.Design(f"d{n}", 1 + n % 4) for n in range(72))
    whole = sized(method.ROUTE_INDEX_PCR, some, 0)
    assert (whole.wells, len(whole.picked), len(whole.index)) == (288, 1, 3)
    assert [len(one.labels) for one in whole.picked] == [288]
    fewer = sized(method.ROUTE_INDEX_PCR, some, 4)
    assert len(fewer.designs) == 18
    assert fewer.wells == 72
    assert [len(one.labels) for one in fewer.picked] == [72]
    assert len(fewer.index) == 1
    assert method.validation(method.ROUTE_INDEX_PCR, some, 5) is None
    assert method.validation(method.ROUTE_INDEX_PCR, some, None) is None


def test_ligation_compresses_four_picked_plates_into_one_and_index_pcr_neither():
    """Each route's plates follow from the shared wells, and neither pours the other's."""
    many = tuple(method.Design(f"d{n}", 2) for n in range(400))
    ligation = sized(method.ROUTE_LIGATION, many, 0)
    assert ligation.wells == 1600
    assert len(ligation.picked) == 5
    assert len(ligation.compressed) == 2
    with pytest.raises(ValueError, match="only the index PCR route"):
        assert ligation.index
    index_pcr = sized(method.ROUTE_INDEX_PCR, many, 0)
    assert len(index_pcr.index) == 17
    assert [one.labels["A1"] for one in index_pcr.index[:2]] == [
        "forward 1, reverse 1",
        "forward 1, reverse 2",
    ]
    with pytest.raises(ValueError, match="only the barcode ligation route"):
        assert index_pcr.compressed


def test_the_steps_print_each_design_chance_beside_the_floor():
    """The number reads as a choice: the floor is stated and the curve is printed beside it."""
    some = (method.Design("two", 2), method.Design("eight", 8))
    one = sized(method.ROUTE_INDEX_PCR, some, 2)
    made = steps.validation_steps(one)
    assert [step.title for step in made][:2] == [
        "Array 2 designs and grow",
        "Pick 4 colonies of each design",
    ]
    assert "2 fragments or more" in " ".join(n.text for n in made[0].noted)
    assert "2 fragments: 1 design, 100.0% of picks clean" in [n.text for n in made[1].noted]
    assert "8 fragments: 1 design, 66.7% of picks clean" in [n.text for n in made[1].noted]


def test_no_step_spells_emphasis_the_page_renders_as_asterisks():
    """A page renders what a step says and nothing more, so Markdown reads as punctuation."""
    some = (method.Design("one", 2),)
    said = [
        text
        for route in (method.ROUTE_INDEX_PCR, method.ROUTE_LIGATION)
        for step in steps.validation_steps(sized(route, some, 2))
        for text in (
            *step.instructions,
            *(n.text for n in step.noted),
            *(one.text for one in step.expectations),
        )
    ]
    assert [text for text in said if "*" in text] == []


def pooling(one) -> Step:
    """Return the step that pools the marked plates and sequences them."""
    return next(step for step in steps.validation_steps(one) if step.key == "pool-and-sequence")


def test_the_pool_is_cleaned_up_before_it_is_amplified_on_the_kits_own_pairs():
    """Nothing counts the flanking pairs, so the page names them and prints no number."""
    said = pooling(sized(method.ROUTE_LIGATION, (method.Design("one", 2),), 0))
    assert [line.split(" ", 1)[0] for line in said.instructions] == [
        "Pool",
        "Clean",
        "Amplify",
        "Sequence",
    ]
    written = " ".join(said.instructions)
    assert f"{method.POOL_COLUMNS} miniprep columns" in written
    assert "the barcode kit's own flanking primer pairs" in written
    assert "three" not in written.lower()


def test_the_flow_cell_note_counts_the_plates_this_read_pools():
    """Both numbers follow from the run: what it pools, and what its address could tell apart."""
    one = sized(method.ROUTE_LIGATION, (method.Design("one", 2),), 0)
    assert len(steps.pooled_plates(one)) == 1
    note = " ".join(n.text for n in pooling(one).noted)
    assert "This read pools 1 plate." in note
    assert f"tell {method.ROUTE_LIGATION.plate_axis} plates apart" in note
    many = sized(method.ROUTE_INDEX_PCR, tuple(method.Design(f"d{n}", 2) for n in range(50)), 0)
    assert len(steps.pooled_plates(many)) == 3
    assert "This read pools 3 plates." in " ".join(n.text for n in pooling(many).noted)


def test_index_pcr_carries_a_hole_at_the_marks_and_ligation_carries_none():
    """The 192 index sequences are lab stock, which only the lab's own plate fills."""
    some = (method.Design("one", 2),)
    index_pcr = steps.validation_steps(sized(method.ROUTE_INDEX_PCR, some, 0))
    assert [hole.id for step in index_pcr for hole in step.holes] == ["IDX1"]
    ligation = steps.validation_steps(sized(method.ROUTE_LIGATION, some, 0))
    assert [hole.id for step in ligation for hole in step.holes] == []


def test_the_plates_a_pick_fills_name_where_their_numbers_were_read():
    """Both routes pick into these two, so a reader of either page can follow the numbers back."""
    some = (method.Design("one", 2),)
    for route in (method.ROUTE_LIGATION, method.ROUTE_INDEX_PCR):
        cited = {
            one.name: one.citation
            for one in method.validation_materials(sized(route, some, 0))
            if one.citation
        }
        assert cited["25 cm BioAssay plate"] == Citation("Qian SI", "Day 2")
        assert cited[f"{method.PICKED_WELLS}-well culture plate"] == Citation("Qian SI", "Day 3")
        assert all(one.source in method.SOURCES for one in cited.values())


def test_the_index_pcr_is_one_wells_share_of_levseqs_published_mix():
    """Scaled back to a full plate the table is the SI's own, and the Taq carries its units."""
    table = method.index_pcr_reaction()
    assert round(sum(one.volume_ul for one in table.components), 2) == method.INDEX_PCR_UL
    assert table.mix_volumes(method.INDEX_WELLS)[:4] == (144.0, 28.8, 7.2, 57.6)
    taq = next(one for one in table.components if one.name.startswith("Taq"))
    assert (taq.volume_ul, taq.stock, taq.final) == (0.05, "5 U/µL", "0.25 units")


def test_the_index_pcr_touches_down_before_it_plateaus():
    """One stepping stage of ten cycles, then 25 more: 35 in all, the SI's two loops as two."""
    stages = method.index_pcr_program().stages
    assert [stage.cycles for stage in stages] == [1, 10, 25, 1, 1]
    anneal = stages[1].incubations[1]
    assert (anneal.temperature_c, anneal.delta_c, anneal.last_c(10)) == (68.0, -0.5, 63.5)


def said_by(one) -> str:
    """Every line this read-back prints that could name a drug."""
    return " ".join(
        (
            *(plate.holds for plate in one.picked),
            *(material.note or "" for material in method.validation_materials(one)),
            *(line for step in steps.validation_steps(one) for line in step.instructions),
        )
    )


def test_every_plate_is_selected_on_the_drug_the_caller_read_off_the_vector():
    """This method rebuilt its DMX vector KanR, so nothing the read-back pours is carbenicillin."""
    one = method.validation(
        method.ROUTE_LIGATION, (method.Design("one", 2),), 0, selection="50 µg/mL kanamycin"
    )
    assert one is not None

    said = said_by(one)

    assert "carbenicillin" not in said
    assert said.count("50 µg/mL kanamycin") == 5


def test_a_read_that_cannot_name_the_drug_leaves_it_to_the_record():
    """A caller naming none prints the vector's own antibiotic rather than the paper's."""
    said = said_by(sized(method.ROUTE_LIGATION, (method.Design("one", 2),), 0))

    assert "carbenicillin" not in said
    assert said.count("the vector's own antibiotic") == 5

"""DMX seating: one part a well, simulated in its carrier, and no library at the end."""

import pytest

from liulab_mbio.sequence import Feature, Segment, SequenceRecord
from liulab_mbio.sites import digest, find_sites
from liulab_synbio.dmx import seating

#: A carrier the size of a test: one blunt point, no site the releasing enzyme reads, and a
#: feature the insertion falls inside so the report has something to say.
BACKBONE = "ACGTACGTACGT" + seating.TOPO_SITE + "TTGCATGCATGCAT"
CARRIER = SequenceRecord(
    BACKBONE,
    topology="circular",
    name="carrier",
    features=(Feature("marker", "CDS", (Segment(6, 30),)),),
)

#: A part: an inward-facing BsmBI pair around a core, released on TATG and AGGA.
PART = "CGTCTCATATGGGCAGCAGGATGAGACG"


def part(name: str) -> SequenceRecord:
    """Return the test part under a name, since a well is named by the part in it."""
    return SequenceRecord(PART, name=name)


def test_a_part_is_seated_at_the_blunt_point_and_can_be_cut_back_out():
    """The point of taking records: the carrier plasmid is simulated, not described."""
    seated, report = seating.seat(CARRIER, part("FLAG"))
    assert len(seated) == len(BACKBONE) + len(PART)
    assert seated.topology == "circular"
    assert [feature.name for feature in report.changed] == ["marker"]
    assert len(find_sites(seated, seating.ENZYME)) == 2
    release, backbone = sorted(digest(seated, seating.ENZYME), key=lambda piece: piece.length)
    assert (release.left_overhang, release.right_overhang) == seating.released(part("FLAG"))
    assert seated.bases(release.start, release.end) in PART
    assert backbone.length == len(seated) - release.length


def test_a_carrier_the_releasing_enzyme_reads_is_refused():
    """Releasing the part from such a carrier would cut the backbone with it."""
    leaky = SequenceRecord("CGTCTCAAAA" + BACKBONE, topology="circular", name="leaky")
    with pytest.raises(ValueError, match="cut the backbone"):
        seating.seat(leaky, part("FLAG"))


def test_a_carrier_without_exactly_one_blunt_point_is_refused():
    """Topoisomerase opens the carrier at one run, so none or two leaves the product unsaid."""
    with pytest.raises(ValueError, match="no GCCCTTAAGGGC run"):
        seating.seat(SequenceRecord("ACGTACGTACGT", topology="circular"), part("FLAG"))
    twice = SequenceRecord(BACKBONE + seating.TOPO_SITE, topology="circular", name="twice")
    with pytest.raises(ValueError, match="2 GCCCTTAAGGGC runs"):
        seating.topo_site(twice)


def test_a_part_the_enzyme_cannot_release_is_refused():
    """A part that cannot be cut back out of its carrier is not a part."""
    with pytest.raises(ValueError, match="carries 0 BsmBI site"):
        seating.seat(CARRIER, SequenceRecord("ACGTACGTACGT", name="dull"))
    with pytest.raises(ValueError, match="no pair the method declares"):
        seating.seat(CARRIER, part("FLAG"), overhangs=[("AGGA", "TTCC")])


def test_seating_ends_in_one_plasmid_a_well_and_no_pool():
    """A round hands back a library; this hands back parts that are still told apart."""
    seated = seating.seat_parts(
        [part(name) for name in ("FLAG", "GS linker", "NLS", "degron")],
        carrier=CARRIER,
        overhangs=[("TATG", "AGGA")],
    )
    assert seated.products == 4
    assert sorted(seated.plate.seating) == ["A1", "A2", "A3", "A4"]
    assert [record.name for record in seated.records] == [
        f"{name} in carrier" for name in seated.parts
    ]
    assert not hasattr(seated, "library")


def test_the_seating_reaction_is_topoisomerase_and_names_no_enzyme_the_user_adds():
    """The carrier backbone has no BsmBI site, so a BsmBI reaction on it cannot run."""
    step = seating.seating_step(seating.seat_parts([part("FLAG"), part("NLS")], carrier=CARRIER))
    assert seating.ENZYME not in " ".join(step.instructions)
    assert not any(seating.ENZYME in material.name for material in seating.materials())
    assert "topoisomerase" in " ".join(step.instructions).lower()
    assert "2 carrier plasmids" in " ".join(step.expected)
    assert "no pool and no library" in " ".join(step.expected).lower()


def test_the_seating_step_carries_the_handle_its_builder_assigned():
    """No pipeline writes this step, so `tests/test_step_keys.py` never sees its key."""
    seated = seating.seat_parts([part("FLAG")], carrier=CARRIER)
    assert seating.seating_step(seated).key == "seat-parts"


def test_the_plate_is_the_smallest_format_that_holds_the_parts():
    """Format is one parameter, so seating picks rather than hard-coding a plate."""
    eleven = [part(f"p{n}") for n in range(11)]
    assert seating.seat_parts(eleven, carrier=CARRIER).plate.wells == 12
    assert seating.seat_parts([*eleven, part("x"), part("y")], carrier=CARRIER).plate.wells == 24


def test_every_part_needs_its_own_name():
    """A well is named by the part sitting in it, so a repeat would lose one of them."""
    with pytest.raises(ValueError, match="share a name"):
        seating.seat_parts([part("NLS"), part("NLS")], carrier=CARRIER)
    with pytest.raises(ValueError, match="has no name"):
        seating.seat_parts([part("")], carrier=CARRIER)
    with pytest.raises(ValueError, match="at least one"):
        seating.seat_parts([], carrier=CARRIER)

"""The DMX part carrier: one part a well, simulated in its plasmid, and no library at the end."""

import pytest

from mbio.sequence import Feature, Segment, SequenceRecord
from mbio.sites import digest, find_sites
from synbio.dmx.carrier import (
    ENZYME,
    TOPO_SITE,
    carrier_step,
    materials,
    released,
    seat,
    seat_parts,
    topo_site,
)

#: A carrier the size of a test: one blunt point, no site the releasing enzyme reads, and a
#: feature the insertion falls inside so the report has something to say.
BACKBONE = "ACGTACGTACGT" + TOPO_SITE + "TTGCATGCATGCAT"
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
    seated, report = seat(CARRIER, part("FLAG"))
    assert len(seated) == len(BACKBONE) + len(PART)
    assert seated.topology == "circular"
    assert [feature.name for feature in report.changed] == ["marker"]
    assert len(find_sites(seated, ENZYME)) == 2
    release, backbone = sorted(digest(seated, ENZYME), key=lambda piece: piece.length)
    assert (release.left_overhang, release.right_overhang) == released(part("FLAG"))
    assert seated.bases(release.start, release.end) in PART
    assert backbone.length == len(seated) - release.length


def test_a_carrier_the_releasing_enzyme_reads_is_refused():
    """Releasing the part from such a carrier would cut the backbone with it."""
    leaky = SequenceRecord("CGTCTCAAAA" + BACKBONE, topology="circular", name="leaky")
    with pytest.raises(ValueError, match="cut the backbone"):
        seat(leaky, part("FLAG"))


def test_a_carrier_without_exactly_one_blunt_point_is_refused():
    """Topoisomerase opens the carrier at one run, so none or two leaves the product unsaid."""
    with pytest.raises(ValueError, match="no GCCCTTAAGGGC run"):
        seat(SequenceRecord("ACGTACGTACGT", topology="circular"), part("FLAG"))
    twice = SequenceRecord(BACKBONE + TOPO_SITE, topology="circular", name="twice")
    with pytest.raises(ValueError, match="2 GCCCTTAAGGGC runs"):
        topo_site(twice)


def test_a_part_the_enzyme_cannot_release_is_refused():
    """A part that cannot be cut back out of its carrier is not a part."""
    with pytest.raises(ValueError, match="carries 0 BsmBI site"):
        seat(CARRIER, SequenceRecord("ACGTACGTACGT", name="dull"))
    with pytest.raises(ValueError, match="no pair the method declares"):
        seat(CARRIER, part("FLAG"), overhangs=[("AGGA", "TTCC")])


def test_the_carrier_ends_in_one_plasmid_a_well_and_no_pool():
    """A round hands back a library; this hands back parts that are still told apart."""
    seated = seat_parts(
        [part(name) for name in ("FLAG", "GS linker", "NLS", "degron")],
        carrier=CARRIER,
        overhangs=[("TATG", "AGGA")],
    )
    assert seated.products == 4
    assert sorted(seated.plate.labels) == ["A1", "A2", "A3", "A4"]
    assert [record.name for record in seated.records] == [
        f"{name} in carrier" for name in seated.parts
    ]
    assert not hasattr(seated, "library")


def test_the_carrier_reaction_is_topoisomerase_and_names_no_enzyme_the_user_adds():
    """The carrier backbone has no BsmBI site, so a BsmBI reaction on it cannot run."""
    step = carrier_step(seat_parts([part("FLAG"), part("NLS")], carrier=CARRIER))
    assert ENZYME not in " ".join(step.instructions)
    assert not any(ENZYME in material.name for material in materials())
    assert "topoisomerase" in " ".join(step.instructions).lower()
    said = " ".join(one.text for one in step.expectations)
    assert "2 carrier plasmids" in said
    assert "no pool and no library" in said.lower()


def test_the_carrier_step_carries_the_handle_its_builder_assigned():
    """No pipeline writes this step, so `tests/test_step_keys.py` never sees its key."""
    seated = seat_parts([part("FLAG")], carrier=CARRIER)
    assert carrier_step(seated).key == "seat-parts"


def test_the_plate_is_the_smallest_format_that_holds_the_parts():
    """Format is one parameter, so the carrier step picks rather than hard-coding a plate."""
    eleven = [part(f"p{n}") for n in range(11)]
    assert seat_parts(eleven, carrier=CARRIER).plate.wells == 12
    assert seat_parts([*eleven, part("x"), part("y")], carrier=CARRIER).plate.wells == 24


def test_every_part_needs_its_own_name():
    """A well is named by the part sitting in it, so a repeat would lose one of them."""
    with pytest.raises(ValueError, match="share a name"):
        seat_parts([part("NLS"), part("NLS")], carrier=CARRIER)
    with pytest.raises(ValueError, match="has no name"):
        seat_parts([part("")], carrier=CARRIER)
    with pytest.raises(ValueError, match="at least one"):
        seat_parts([], carrier=CARRIER)

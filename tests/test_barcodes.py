import dataclasses
import json
from itertools import combinations
from pathlib import Path

import pytest

from liulab_mbio.barcodes import (
    GC_BAND,
    MAX_HOMOPOLYMER,
    METRIC,
    MIN_DISTANCE,
    BarcodeRules,
    SpaceExhaustedError,
    check_barcodes,
    deletion_ambiguity,
    design_barcodes,
)
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.sites import has_site

# The published scheme: an 11-base barcode joined to the one before it by a four-base cloning
# scar, so a barcode and its scar are five codons and the block stays in frame.
LENGTH = 11
SCAR = "AGCG"
FORBIDDEN = ["BsaI", "BbsI"]

# Measured below, not assumed: the published barcode begins one base into a codon.
PHASE = 1

RULES = BarcodeRules(LENGTH, scar=SCAR, forbidden=FORBIDDEN)

# Spells no BsmBI site on its own; the scar in front of it completes CGTCTC.
JUNCTION_SITE = "TCTCAAGT"

# CAC CAT TA|A: the stop is finished by the first base of the scar, not by the barcode.
JUNCTION_STOP = "CACCATTA"

# The reverse complement of the BsaI site, GGTCTC, inside the barcode itself.
REVERSE_SITE = "GAGACCAAGTT"


def distance(one: str, other: str) -> int:
    """Count the mismatches, independently of the module under test."""
    return sum(a != b for a, b in zip(one, other, strict=True))


@pytest.fixture(scope="module")
def published(data_dir: Path) -> dict[str, list[str]]:
    """The three published part lists, keyed by the position each fills."""
    data = json.loads((data_dir / "ap1-barcodes.json").read_text(encoding="utf-8"))
    return data["part_lists"]


def test_a_set_is_the_size_and_the_length_that_was_asked_for() -> None:
    found = design_barcodes(24, RULES)
    assert len(found) == 24
    assert len(set(found)) == 24
    assert {len(one) for one in found} == {LENGTH}


def test_every_pair_within_a_part_list_stands_the_distance_apart() -> None:
    found = design_barcodes(12, RULES)
    pairs = list(combinations(found, 2))
    assert len(pairs) == 66
    assert min(distance(one, other) for one, other in pairs) >= MIN_DISTANCE


def test_no_barcode_carries_a_forbidden_site_in_either_orientation() -> None:
    for barcode in design_barcodes(12, RULES):
        # With the scar on either side, so a site spanning the junction would be found too.
        assert not has_site(SequenceRecord(SCAR + barcode + SCAR), FORBIDDEN)


def test_a_site_the_barcode_spells_backwards_is_found() -> None:
    assert has_site(SequenceRecord(REVERSE_SITE), "BsaI")
    (problem,) = check_barcodes([REVERSE_SITE], RULES)
    assert problem == f"{REVERSE_SITE} spells a BsaI site"


def test_a_site_only_the_cloning_scar_completes_is_found() -> None:
    rules = BarcodeRules(len(JUNCTION_SITE), scar=SCAR, forbidden=["BsmBI"])
    assert not has_site(SequenceRecord(JUNCTION_SITE), "BsmBI")
    (problem,) = check_barcodes([JUNCTION_SITE], rules)
    assert problem == (
        f"{JUNCTION_SITE} spells a BsmBI site across the junction with the cloning scar"
    )


def test_no_barcode_spells_a_stop_in_the_frame_the_construct_reads() -> None:
    for barcode in design_barcodes(12, RULES):
        codons = [(barcode + SCAR)[at : at + 3] for at in range(0, LENGTH + len(SCAR), 3)]
        assert not {"TAA", "TAG", "TGA"} & set(codons)


def test_a_stop_only_the_cloning_scar_completes_is_found() -> None:
    rules = BarcodeRules(len(JUNCTION_STOP), scar=SCAR)
    (problem,) = check_barcodes([JUNCTION_STOP], rules)
    assert problem == f"{JUNCTION_STOP} spells TAA, a stop, in the frame the construct reads"
    # Nothing translates it, so neither the frame rule nor the stop rule applies.
    assert check_barcodes([JUNCTION_STOP], BarcodeRules(len(JUNCTION_STOP), phase=None)) == ()


def test_the_indel_aware_metric_reads_its_paper_s_worked_example_as_two_edits() -> None:
    # Buschmann & Bystrykh 2013: CAGG and CGTC differ at three positions, but deleting the A and
    # substituting one base leaves CGT, which a read running on into a C spells as CGTC.
    rules = BarcodeRules(4, phase=None, max_homopolymer=None, metric="sequence-levenshtein")
    assert check_barcodes(["CAGG", "CGTC"], rules) == (
        "CAGG and CGTC stand 2 edit(s) apart, under the 3 one part list needs",
    )
    assert check_barcodes(["CAGG", "CGTC"], dataclasses.replace(rules, metric="hamming")) == ()


def test_one_deletion_can_leave_two_barcodes_reading_alike_and_the_share_is_counted() -> None:
    # Worked by hand: dropping the first base of ACGTT and the last of CGTTG both leave CGTT, so
    # one deletion of each is ambiguous, of the ten two five-base barcodes hold between them.
    assert deletion_ambiguity(["ACGTT", "CGTTG"]) == 0.2
    assert deletion_ambiguity(["ACGTT"]) == 0.0


def test_the_cost_of_each_metric_is_the_one_the_other_does_not_pay() -> None:
    # A Hamming set leaves deletions that read as another barcode; the indel-aware set does not.
    # It pays in space instead: at five bases it runs out where the Hamming rule still fills.
    short = BarcodeRules(5, phase=None, metric="hamming")
    indel = dataclasses.replace(short, metric="sequence-levenshtein")
    assert deletion_ambiguity(design_barcodes(24, short)) > 0
    assert deletion_ambiguity(design_barcodes(8, indel)) == 0
    with pytest.raises(SpaceExhaustedError, match="24 barcodes of 5 bases"):
        design_barcodes(24, indel)


def test_a_metric_that_is_neither_is_refused_naming_both() -> None:
    with pytest.raises(ValueError, match="'hamming' or 'sequence-levenshtein'"):
        BarcodeRules(LENGTH, scar=SCAR, metric="levenshtein")  # type: ignore[arg-type]


def test_a_barcode_and_scar_that_are_not_whole_codons_are_refused() -> None:
    with pytest.raises(ValueError, match="not a whole number of codons"):
        BarcodeRules(10, scar=SCAR)


def test_the_dials_default_to_a_homopolymer_cap_and_to_no_gc_band() -> None:
    assert (MIN_DISTANCE, MAX_HOMOPOLYMER, GC_BAND, METRIC) == (3, 5, None, "sequence-levenshtein")
    # A run over the cap is rejected by default, and turning the cap off accepts it.
    run = "AAAAAACGTCG"
    assert check_barcodes([run], RULES)[0].endswith("run of 6 A, over the cap of 5")
    assert check_barcodes([run], BarcodeRules(LENGTH, scar=SCAR, max_homopolymer=None)) == ()
    # No GC band by default, and a band is there to be turned on.
    lean = "AACAAATTAGT"
    assert check_barcodes([lean], RULES) == ()
    banded = BarcodeRules(LENGTH, scar=SCAR, gc_band=(0.4, 0.6))
    assert check_barcodes([lean], banded)[0].endswith("outside the band 40% to 60%")


def test_the_same_request_returns_the_same_set_and_another_seed_another() -> None:
    assert design_barcodes(8, RULES) == design_barcodes(8, RULES)
    assert design_barcodes(8, RULES) != design_barcodes(8, RULES, seed=1)


def test_the_draw_is_shuffled_and_not_a_walk_of_the_space_in_order() -> None:
    # A walk in order opens on a run of one base and stays poor in G and C; a seeded shuffle
    # spans the range. `docs/research/barcode-design.md` measures both.
    content = sorted(sum(base in "GC" for base in one) for one in design_barcodes(24, RULES))
    assert content[0] <= 4
    assert content[-1] >= 7


def test_a_set_bigger_than_the_space_is_refused_naming_the_rule_that_took_it() -> None:
    with pytest.raises(SpaceExhaustedError) as refusal:
        design_barcodes(40, BarcodeRules(3))
    said = str(refusal.value)
    assert "40 barcodes of 3 bases were asked for" in said
    assert "the distance rule (at least 3 edits within a part list) rejected" in said
    assert "the stop-codon rule rejected" in said


def test_the_published_set_holds_every_rule_this_module_designs_to(
    published: dict[str, list[str]],
) -> None:
    # The cap of 5 is the one default the published set does not meet, so it is turned off here
    # and measured on its own below.
    rules = BarcodeRules(LENGTH, scar=SCAR, phase=PHASE, forbidden=FORBIDDEN, max_homopolymer=None)
    for name, barcodes in published.items():
        assert len(barcodes) == 24, name
        # Under both metrics. The set was designed on mismatches alone, and it stands three edits
        # apart as well, so the indel-aware default rejects nothing that was published — and no
        # deletion of one of its barcodes reads as another.
        assert check_barcodes(barcodes, rules) == (), name
        assert check_barcodes(barcodes, dataclasses.replace(rules, metric="hamming")) == (), name
        assert deletion_ambiguity(barcodes) == 0.0, name


def test_the_published_part_lists_stand_further_apart_within_than_across(
    published: dict[str, list[str]],
) -> None:
    within = {
        name: min(distance(one, other) for one, other in combinations(barcodes, 2))
        for name, barcodes in published.items()
    }
    assert within == {"N": 3, "bZIP": 4, "C": 4}
    across = min(
        distance(one, other)
        for first, second in combinations(published, 2)
        for one in published[first]
        for other in published[second]
    )
    # Two, which is why the distance rule is held within a part list and not across the library.
    assert across == 2


def test_the_default_cap_rejects_one_published_barcode_and_the_inherited_band_would_reject_half(
    published: dict[str, list[str]],
) -> None:
    every = [one for barcodes in published.values() for one in barcodes]
    assert len(every) == 72
    # Distance is a part-list rule, and this measures the composition dials over all three, so
    # it is turned down to the one mismatch that only says the barcodes are distinct.
    capped = BarcodeRules(LENGTH, scar=SCAR, phase=PHASE, distance=1)
    assert [problem.split()[0] for problem in check_barcodes(every, capped)] == ["GCTTTTTTGGC"]
    # The 40-60% band the literature inherited would throw out more than half of a set that
    # worked, which is the measurement that keeps the band off by default.
    banded = BarcodeRules(
        LENGTH, scar=SCAR, phase=PHASE, distance=1, max_homopolymer=None, gc_band=(0.4, 0.6)
    )
    assert len(check_barcodes(every, banded)) == 37


def test_the_published_block_is_read_one_base_into_a_codon(
    published: dict[str, list[str]],
) -> None:
    # Nothing states which frame the block is read in, so it is measured here rather than
    # assumed: at phase 1 no published barcode spells a stop, and at either other phase many do.
    every = [one for barcodes in published.values() for one in barcodes]
    stops: dict[int, int] = {}
    for phase in (0, 1, 2):
        rules = BarcodeRules(LENGTH, scar=SCAR, phase=phase, max_homopolymer=None, distance=1)
        stops[phase] = sum("a stop" in problem for problem in check_barcodes(every, rules))
    assert stops == {0: 16, 1: 0, 2: 30}
    assert PHASE == 1

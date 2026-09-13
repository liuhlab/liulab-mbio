import dataclasses

from liulab_mbio.checks import Check
from liulab_mbio.primers import THRESHOLDS, THRESHOLDS_FOR, Band, evaluate_primer, reading
from liulab_mbio.sequence import Primer

from .sequences import M13_FWD

#: Genewiz's free M13F universal sequencing primer, and a 16-mer sequencing primer annealing
#: upstream of the first pUC19-GFP junction.
GENEWIZ_M13F = "GTAAAACGACGGCCAG"
PLAN_SEQUENCING_FORWARD = "AACTGTTGGGAAGGGC"


def test_a_sixteen_mer_is_short_for_pcr_not_for_sequencing_and_a_poor_primer_still_warns() -> None:
    for sequence in (GENEWIZ_M13F, PLAN_SEQUENCING_FORWARD):
        lengths = {
            role: evaluate_primer(Primer(role, sequence), thresholds=THRESHOLDS_FOR[role])["length"]
            for role in ("amplification", "colony PCR", "sequencing")
        }
        assert {role: check.status for role, check in lengths.items()} == {
            "amplification": "warn",
            "colony PCR": "warn",
            "sequencing": "pass",
        }
    # Shorter than any universal primer Genewiz offers, and faultless otherwise.
    short = evaluate_primer(
        Primer("short", "ACGGTCACAGCTTGT"), thresholds=THRESHOLDS_FOR["sequencing"]
    )
    assert [check.name for check in short.checks if check.status not in (None, "pass")] == [
        "length"
    ]
    assert short.status == "warn"


def test_a_check_reads_as_a_label_a_value_and_the_band_it_was_held_to() -> None:
    report = evaluate_primer(Primer("M13 fwd", M13_FWD))
    length = reading(report["length"])
    assert (length.label, length.value, length.limit) == ("length", "17", "band 18-30")
    assert length.detail == "17 (band 18-30)"
    assert reading(report["gc_percent"]).detail == "53% (band 40-60)"
    assert reading(report["tm"]).detail == "62.3 °C (band 60-64)"
    # A band no published rule sets says so.
    assert reading(report["gc_clamp"]).detail == "3 (proposed band 1-3)"
    assert reading(Check("dinucleotide_repeat", "warn", 4)).limit == "proposed band 0-3"
    # A check nothing judged has no band to print, and says only what it measured.
    unjudged = reading(report["end_stability"])
    assert (unjudged.label, unjudged.limit) == ("3' end stability", "")
    assert unjudged.detail == "-4.0 kcal/mol"


def test_a_reading_is_held_to_the_thresholds_it_is_given() -> None:
    strict = dataclasses.replace(THRESHOLDS, gc_percent=Band(45.0, 55.0))
    assert reading(Check("gc_percent", "warn", 39.1), strict).limit == "band 45-55"
    assert reading(Check("hairpin", "warn", 52.0)).limit == "max 47"
    assert reading(Check("binding_sites", "fail", 0)).limit == "exactly 1"


def test_every_threshold_comes_from_one_place() -> None:
    strict = dataclasses.replace(THRESHOLDS, tm=Band(63.0, 64.0, 62.5, 65.0))
    report = evaluate_primer(Primer("M13 fwd", M13_FWD), thresholds=strict)
    assert report["tm"].status == "fail"

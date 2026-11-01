"""A record handed to molviz: what mbio converts, held to molviz's own check of its format."""

from pathlib import Path

import molviz

from mbio.io import read_record
from mbio.plot.converter import molviz_options, molviz_record
from mbio.sequence import Feature, Segment, SequenceRecord, Strand


def test_a_record_converts_to_one_molviz_accepts(
    puc19: SequenceRecord, gfp: SequenceRecord, data_dir: Path
) -> None:
    # Circular with a CDS across a join, linear, and with primers that mismatch.
    for record in (puc19, gfp, read_record(data_dir / "primer-test.dna")):
        converted = molviz_record(record)
        molviz.check(converted)
        assert converted["bases"] == record.sequence, record.name


def test_a_regulatory_feature_of_a_terminating_class_takes_a_terminators_colour() -> None:
    def regulatory(name: str, kind: str) -> Feature:
        return Feature(
            name,
            "regulatory",
            (Segment(0, 10),),
            strand=Strand.FORWARD,
            qualifiers={"regulatory_class": (kind,)},
        )

    record = SequenceRecord(
        "ACGT" * 10,
        features=(
            regulatory("rrnB T1", "terminator"),
            regulatory("SV40 pA", "polyA_signal_sequence"),
            regulatory("lac", "promoter"),
        ),
    )
    # molviz sees no qualifiers, so it would colour all three as promoters.
    colors = [feature.get("color") for feature in molviz_record(record)["features"]]
    assert colors == ["#ffaabb", "#ffaabb", None]


def test_every_shipped_enzyme_cutting_the_record_is_handed_over_however_often_it_cuts(
    puc19: SequenceRecord,
) -> None:
    cuts = molviz_record(puc19)["cuts"]
    enzymes = [cut["enzyme"] for cut in cuts]
    # BsaI cuts pUC19 once and BsmBI twice; molviz picks the unique cutters out itself.
    assert enzymes.count("BsaI") == 1
    assert enzymes.count("BsmBI") == 2
    [ecori] = [cut for cut in cuts if cut["enzyme"] == "EcoRI"]
    # G^AATTC: the top strand cut after the G, the bottom four bases on.
    assert ecori["at"] == ecori["start"] + 1
    assert ecori["stagger"] == 4


def test_a_cds_carries_the_protein_snapgene_reads_from_it(puc19: SequenceRecord) -> None:
    [feature] = [one for one in molviz_record(puc19)["features"] if one["name"] == "AmpR"]
    # SnapGene marks where AmpR's segments join with a comma.
    assert feature["translation"].startswith("MSIQHFRVALIPFFAAFCLPVFAHPETLV")
    assert feature["translation"].endswith("IKHW*")


def test_a_primer_carries_its_whole_sequence_and_where_it_mismatches(data_dir: Path) -> None:
    record = read_record(data_dir / "primer-test.dna")
    primers = {one["name"]: one for one in molviz_record(record)["primers"]}
    assert primers["rev-mismatch"]["bases"] == "CCGAATTCGAGTTAGCTAACTCATTAGGCAC"
    assert primers["rev-mismatch"]["sites"] == [
        {"start": 250, "end": 272, "strand": -1, "mismatches": [263]}
    ]
    assert primers["plain-fwd"]["sites"] == [{"start": 120, "end": 140, "strand": 1}]


def test_options_name_each_shipped_enzyme_and_leave_out_each_of_molvizs_defaults() -> None:
    every = molviz_options(
        span=(2679, 2696),
        linear=True,
        sequence_view=True,
        bases_per_row=30,
        both_strands=False,
        features=False,
        primers=False,
        cut_sites=False,
        hide_types=("CDS",),
        source=True,
        enzymes=("EcoRI-HF", "EcoRI", "BamHI"),
        highlight="lacZ",
    )
    assert every == {
        "span": {"start": 2679, "end": 2696},
        "map": "line",
        "sequence_view": True,
        "bases_per_row": 30,
        "bottom_strand": False,
        "kinds_off": ["features", "primers", "cuts"],
        "types_off": ["CDS"],
        # molviz ignores a name that cuts nothing, so a commercial one goes as the shipped enzyme.
        "enzymes": ["EcoRI", "BamHI"],
        "highlight": ["lacZ"],
    }
    # Left at molviz's defaults, nothing goes over: a source feature hidden, all else shown.
    assert molviz_options() == {}
    assert molviz_options(hide_types=("source",)) == {}
    assert molviz_options(source=True) == {"types_off": []}

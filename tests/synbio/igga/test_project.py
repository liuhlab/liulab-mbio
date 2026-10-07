"""What a project file holds, what it composes with the method, and what it refuses when read."""

import json
from pathlib import Path

import pytest

from liulab_synbio.igga.method import IGGA
from liulab_synbio.igga.project import Barcode, Project, read_project

DEMO = Path(__file__).parents[3] / "docs" / "examples" / "ap1-library" / "project.json"

#: A whole project, as a user writes one.
WRITTEN = {
    "name": "demo",
    "positions": ["N", "DBD", "C"],
    "parts": "parts.fasta",
    "vector": "vector.dna",
    "host": "human",
    "oligo_length": 350,
    "batch_size": 96,
    "coverage": 300,
}


def write(directory: Path, **changes) -> Path:
    """Write a project and the two files it names, with any key replaced."""
    (directory / "parts.fasta").write_text(">N_a\nMKTAEK\n")
    (directory / "vector.dna").write_text("not read here")
    path = directory / "project.json"
    path.write_text(json.dumps({**WRITTEN, **changes}), encoding="utf-8")
    return path


def test_the_ap1_project_reads_as_what_the_demo_plans():
    made = read_project(DEMO)

    assert made.name == "AP-1 DESynR"
    assert made.positions == ("N", "DBD", "C")
    assert made.parts.name == "parts.fasta"
    assert made.vector.is_file()
    assert (made.host, made.coverage, made.seed) == ("human", 300.0, 0)
    assert (made.oligo_length, made.batch_size) == (350, 96)
    assert (made.barcode.length, made.barcode.min_distance) == (11, 3)
    assert made.scheme is IGGA
    assert made.retained_length == 75


def test_a_path_is_resolved_against_the_project_file(tmp_path):
    made = read_project(write(tmp_path))

    assert made.parts == tmp_path / "parts.fasta"
    assert made.vector == tmp_path / "vector.dna"


def test_the_method_s_defaults_stand_where_a_project_states_nothing(tmp_path):
    made = read_project(write(tmp_path))

    assert (made.barcode.length, made.barcode.min_distance) == (11, 3)
    assert made.reserved_extra == ()


def test_reserved_enzymes_compose_rather_than_replace(tmp_path):
    made = read_project(write(tmp_path, reserved_extra=["EcoRI"]))

    assert made.reserved == ("BsmBI", "EcoRI")
    assert [one.name for one in made.reserved_enzymes] == ["BsmBI", "EcoRI"]
    # Whatever a project says, what the method reserves stays reserved.
    assert IGGA.reserved[0] in made.reserved


def test_the_barcode_frame_rule_refuses_12_and_admits_14(tmp_path):
    with pytest.raises(ValueError, match="barcode-frame"):
        read_project(write(tmp_path, barcode={"length": 12}))

    assert read_project(write(tmp_path, barcode={"length": 14})).barcode.length == 14


def test_a_position_named_twice_is_refused(tmp_path):
    with pytest.raises(ValueError, match="twice"):
        read_project(write(tmp_path, positions=["N", "C", "N"]))


def test_a_project_with_no_position_is_refused(tmp_path):
    with pytest.raises(ValueError, match="at least one position"):
        read_project(write(tmp_path, positions=[]))


def test_a_coverage_that_is_not_positive_is_refused(tmp_path):
    with pytest.raises(ValueError, match="coverage"):
        read_project(write(tmp_path, coverage=0))


def test_a_path_naming_no_file_is_refused(tmp_path):
    with pytest.raises(ValueError, match="is no file"):
        read_project(write(tmp_path, vector="nowhere.dna"))


def test_an_unshipped_host_is_refused(tmp_path):
    with pytest.raises(KeyError, match="codon usage table"):
        read_project(write(tmp_path, host="nowhere"))


def test_an_unshipped_reserved_enzyme_is_refused(tmp_path):
    with pytest.raises(KeyError):
        read_project(write(tmp_path, reserved_extra=["NotAnEnzyme"]))


def test_a_missing_key_is_refused_naming_it(tmp_path):
    path = tmp_path / "thin.json"
    path.write_text(json.dumps({"name": "thin"}), encoding="utf-8")

    with pytest.raises(ValueError, match="missing batch_size"):
        read_project(path)


def test_an_unknown_key_is_refused_naming_it(tmp_path):
    with pytest.raises(ValueError, match="unknown key"):
        read_project(write(tmp_path, scheme="iGGA"))


def test_a_value_of_another_json_type_is_refused(tmp_path):
    with pytest.raises(ValueError, match="not a whole number"):
        read_project(write(tmp_path, oligo_length="350"))


def test_a_project_built_in_code_is_checked_the_same_way(tmp_path):
    with pytest.raises(ValueError, match="barcode-frame"):
        Project(
            "demo",
            positions=("N",),
            parts=tmp_path / "parts.fasta",
            vector=tmp_path / "vector.dna",
            host="human",
            oligo_length=350,
            batch_size=96,
            coverage=10.0,
            barcode=Barcode(12),
        )

"""What a project file holds, what it composes with the method, and what it refuses when read."""

import json
from pathlib import Path

import pytest

from liulab_synbio.igga.coverage import REPRESENTATION_MARKS
from liulab_synbio.igga.method import IGGA
from liulab_synbio.igga.project import Barcode, Build, read_build

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
    "completeness": 0.99,
}


def write(directory: Path, **changes) -> Path:
    """Write a project and the two files it names, with any key replaced."""
    (directory / "parts.fasta").write_text(">N_a\nMKTAEK\n")
    (directory / "vector.dna").write_text("not read here")
    path = directory / "project.json"
    path.write_text(json.dumps({**WRITTEN, **changes}), encoding="utf-8")
    return path


def test_the_ap1_project_reads_as_what_the_demo_plans():
    made = read_build(DEMO)

    assert made.name == "AP-1 DESynR"
    assert made.positions == ("N", "DBD", "C")
    assert made.parts.name == "parts.fasta"
    assert made.vector.is_file()
    assert (made.host, made.completeness, made.seed) == ("human", 0.99, 0)
    assert (made.oligo_length, made.batch_size) == (350, 96)
    assert (made.barcode.length, made.barcode.min_distance) == (11, 3)
    assert made.scheme is IGGA
    assert made.retained_length == 75


def test_a_path_is_resolved_against_the_project_file(tmp_path):
    made = read_build(write(tmp_path))

    assert made.parts == tmp_path / "parts.fasta"
    assert made.vector == tmp_path / "vector.dna"


def test_the_method_s_defaults_stand_where_a_project_states_nothing(tmp_path):
    made = read_build(write(tmp_path))

    assert (made.barcode.length, made.barcode.min_distance) == (11, 3)
    assert made.reserved_extra == ()


def test_reserved_enzymes_compose_rather_than_replace(tmp_path):
    made = read_build(write(tmp_path, reserved_extra=["EcoRI"]))

    assert made.reserved == ("BsmBI", "EcoRI")
    assert [one.name for one in made.reserved_enzymes] == ["BsmBI", "EcoRI"]
    # Whatever a project says, what the method reserves stays reserved.
    assert IGGA.reserved[0] in made.reserved


def test_the_barcode_frame_rule_refuses_12_and_admits_14(tmp_path):
    with pytest.raises(ValueError, match="barcode-frame"):
        read_build(write(tmp_path, barcode={"length": 12}))

    assert read_build(write(tmp_path, barcode={"length": 14})).barcode.length == 14


def test_a_position_named_twice_is_refused(tmp_path):
    with pytest.raises(ValueError, match="twice"):
        read_build(write(tmp_path, positions=["N", "C", "N"]))


def test_a_project_with_no_position_is_refused(tmp_path):
    with pytest.raises(ValueError, match="at least one position"):
        read_build(write(tmp_path, positions=[]))


@pytest.mark.parametrize("given", [0, 1, 1.5])
def test_a_completeness_that_is_not_a_chance_is_refused(tmp_path, given):
    with pytest.raises(ValueError, match="completeness"):
        read_build(write(tmp_path, completeness=given))


def test_a_path_naming_no_file_is_refused(tmp_path):
    with pytest.raises(ValueError, match="is no file"):
        read_build(write(tmp_path, vector="nowhere.dna"))


def test_an_unshipped_host_is_refused(tmp_path):
    with pytest.raises(KeyError, match="codon usage table"):
        read_build(write(tmp_path, host="nowhere"))


def test_an_unshipped_reserved_enzyme_is_refused(tmp_path):
    with pytest.raises(KeyError):
        read_build(write(tmp_path, reserved_extra=["NotAnEnzyme"]))


def test_a_missing_key_is_refused_naming_it(tmp_path):
    path = tmp_path / "thin.json"
    path.write_text(json.dumps({"name": "thin"}), encoding="utf-8")

    with pytest.raises(ValueError, match="missing batch_size"):
        read_build(path)


def test_an_unknown_key_is_refused_naming_it(tmp_path):
    with pytest.raises(ValueError, match="unknown key"):
        read_build(write(tmp_path, scheme="iGGA"))


def test_a_value_of_another_json_type_is_refused(tmp_path):
    with pytest.raises(ValueError, match="not a whole number"):
        read_build(write(tmp_path, oligo_length="350"))


def test_a_project_built_in_code_is_checked_the_same_way(tmp_path):
    with pytest.raises(ValueError, match="barcode-frame"):
        Build(
            "demo",
            positions=("N",),
            parts=tmp_path / "parts.fasta",
            vector=tmp_path / "vector.dna",
            host="human",
            oligo_length=350,
            batch_size=96,
            completeness=0.99,
            barcode=Barcode(12),
        )


def test_a_project_states_no_floor_and_no_route_by_default(tmp_path):
    """Validation is optional and polyclonal by default, so the package ships no floor."""
    made = read_build(write(tmp_path))

    assert (made.validate_from, made.route) == (None, None)


def test_a_floor_and_a_route_are_stated_together(tmp_path):
    """A floor with no route says which designs are read and not how; a route alone reads none."""
    made = read_build(write(tmp_path, validate_from=0, route="barcode ligation"))
    assert (made.validate_from, made.route) == (0, "barcode ligation")

    with pytest.raises(ValueError, match="'barcode ligation', 'index PCR'"):
        read_build(write(tmp_path, validate_from=3))
    with pytest.raises(ValueError, match="stated together"):
        read_build(write(tmp_path, route="index PCR"))


def test_a_project_reads_no_route_but_the_two(tmp_path):
    with pytest.raises(ValueError, match="route is 'DMX'"):
        read_build(write(tmp_path, validate_from=0, route="DMX"))


def test_a_floor_counts_fragments(tmp_path):
    with pytest.raises(ValueError, match="omit it to read nothing"):
        read_build(write(tmp_path, validate_from=-1, route="barcode ligation"))


def test_the_ap1_project_reads_every_design_back_by_index_pcr():
    made = read_build(DEMO)

    assert (made.validate_from, made.route) == (0, "index PCR")


def test_a_project_names_the_working_vector_it_moves_into_or_none(tmp_path):
    """The vector the library ends in is an application's choice, so a project may name one."""
    (tmp_path / "pWORK.fasta").write_text(">pWORK\nACGT\n", encoding="utf-8")

    assert read_build(write(tmp_path)).working_vector is None
    named = read_build(write(tmp_path, working_vector="pWORK.fasta"))
    assert named.working_vector == tmp_path / "pWORK.fasta"


def test_a_project_may_tighten_each_representation_mark(tmp_path):
    """The sourced mark is a floor the method stands on, as a read depth is in `dmx`."""
    made = read_build(
        write(
            tmp_path,
            representation_seen=0.999,
            representation_skew=4.0,
            reads_per_member=200,
        )
    )

    assert made.marks.seen == 0.999
    assert made.marks.skew == 4.0
    assert made.marks.reads_per_member == 200


def test_a_project_stating_no_mark_takes_joungs(tmp_path):
    assert read_build(write(tmp_path)).marks == REPRESENTATION_MARKS


@pytest.mark.parametrize(
    ("field", "value", "says"),
    [
        ("representation_seen", 0.9, "raise"),
        ("representation_skew", 20.0, "lower the skew ratio"),
        ("reads_per_member", 50, "raise"),
    ],
)
def test_a_project_may_not_loosen_a_representation_mark(tmp_path, field, value, says):
    """Each one is Joung's, and a project that loosened it would be judged by nothing."""
    with pytest.raises(ValueError, match=says):
        read_build(write(tmp_path, **{field: value}))


def test_linkage_fidelity_has_no_default(tmp_path):
    """Nothing published sets a mark for it, so a project stating none is held to none."""
    assert read_build(write(tmp_path)).linkage_fidelity is None
    assert read_build(write(tmp_path, linkage_fidelity=0.9)).linkage_fidelity == 0.9

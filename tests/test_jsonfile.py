"""What a hand-written JSON file is refused for, in the words a reader of the refusal sees."""

import json
import re
from dataclasses import dataclass
from pathlib import Path

import pytest

from mbio import jsonfile


@dataclass(frozen=True)
class Toy:
    """Two fields a reader is built from: one a file must name, one it may leave out."""

    name: str
    seed: int = 0


def write(directory: Path, data: object) -> Path:
    """Write one JSON file and return it."""
    path = directory / "build.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def test_an_object_is_read_and_any_other_json_value_is_refused(tmp_path):
    assert jsonfile.read_object(write(tmp_path, {"name": "demo"})) == {"name": "demo"}

    path = write(tmp_path, ["demo"])
    with pytest.raises(ValueError, match=f"^{re.escape(str(path))} holds list, not an object$"):
        jsonfile.read_object(path)


def test_a_missing_key_and_a_key_nothing_reads_are_both_named(tmp_path):
    required, optional = frozenset({"name", "route"}), frozenset({"selection"})

    with pytest.raises(ValueError, match=r"^a build is missing name, route$"):
        jsonfile.refuse_keys({"selection": "carbenicillin"}, required, optional, "a build")
    with pytest.raises(ValueError, match=r"^a build carries unknown key\(s\) host, seed$"):
        jsonfile.refuse_keys(
            {"name": "demo", "route": "index PCR", "seed": 0, "host": "human"},
            required,
            optional,
            "a build",
        )


@pytest.mark.parametrize(
    ("read", "given", "said"),
    [
        (jsonfile.text, 350, "a build's name is int, not a string"),
        (jsonfile.whole, "350", "a build's name is str, not a whole number"),
        (jsonfile.whole, True, "a build's name is bool, not a whole number"),
        (jsonfile.number, True, "a build's name is bool, not a number"),
        (jsonfile.listing, "N", "a build's name is str, not a list"),
    ],
)
def test_a_value_of_another_json_type_is_refused_naming_the_key(read, given, said):
    """True is not a whole number here, whatever Python says of it."""
    with pytest.raises(ValueError, match=f"^{said}$"):
        read({"name": given}, "name", "a build's")


def test_a_value_read_on_its_own_is_named_by_where_it_sits():
    with pytest.raises(ValueError, match=r"^positions\[2\] is int, not a string$"):
        jsonfile.one_text(3, "positions[2]")


def test_a_named_file_resolves_against_the_file_that_names_it(tmp_path):
    (tmp_path / "parts.fasta").write_text(">N_a\nMKTAEK\n")
    file = write(tmp_path, {"parts": "parts.fasta"})

    assert (
        jsonfile.named_file(file, "parts.fasta", "parts", "a build's") == tmp_path / "parts.fasta"
    )
    said = f"a build's parts is 'nowhere.fasta', and {tmp_path / 'nowhere.fasta'} is no file"
    with pytest.raises(ValueError, match=f"^{re.escape(said)}$"):
        jsonfile.named_file(file, "nowhere.fasta", "parts", "a build's")


def test_a_dataclass_reader_says_its_key_faults_in_the_words_every_caller_shares():
    """`reader` refuses a key the way `refuse_keys` does, so one fault reads one way."""
    read = jsonfile.reader(Toy)

    assert read({"name": "demo", "seed": 3}, "a toy") == Toy("demo", 3)
    with pytest.raises(ValueError, match=r"^a toy is missing name$"):
        read({"seed": 3}, "a toy")
    with pytest.raises(ValueError, match=r"^a toy carries unknown key\(s\) host$"):
        read({"name": "demo", "host": "human"}, "a toy")
    with pytest.raises(ValueError, match=r"^a toy\.seed: expected a whole number, got a string$"):
        read({"name": "demo", "seed": "x"}, "a toy")

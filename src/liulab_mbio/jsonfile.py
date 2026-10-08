"""Read a JSON file of what one run chooses, and say what is wrong with it.

A method is code; what one run of it chooses is a file someone writes by hand. These read that
file's values and refuse a missing key, a key nothing reads, and a value of another JSON type,
in words naming the key and what belongs there.

Nothing here knows which method it reads for. Every refusal takes the subject it names as
`where`, so the reader supplies its own word for the file it is reading — `"a build"`, or the
object inside it a key sits in.
"""

import json
import os
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


def read_object(path: str | os.PathLike[str]) -> Mapping[str, Any]:
    """Parse a JSON file holding one object.

    Raises
    ------
    ValueError
        If the file holds any other JSON value.
    """
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, Mapping):
        raise ValueError(f"{os.fspath(path)} holds {type(data).__name__}, not an object")
    return data


def refuse_keys(
    data: Mapping[str, Any], required: frozenset[str], optional: frozenset[str], where: str
) -> None:
    """Refuse a mapping that is missing a key or carries one its reader does not read.

    Raises
    ------
    ValueError
        Naming the keys and `where` they are.
    """
    if missing := sorted(required - set(data)):
        raise ValueError(f"{where} is missing {', '.join(missing)}")
    if unknown := sorted(set(data) - required - optional):
        raise ValueError(f"{where} carries unknown key(s) {', '.join(unknown)}")


def text(data: Mapping[str, Any], key: str, where: str) -> str:
    """Return one string value.

    Raises
    ------
    ValueError
        If the value is of another JSON type.
    """
    return one_text(data[key], f"{where} {key}")


def one_text(value: Any, where: str) -> str:
    """Return one string.

    Raises
    ------
    ValueError
        If the value is of another JSON type.
    """
    if not isinstance(value, str):
        raise ValueError(f"{where} is {type(value).__name__}, not a string")
    return value


def whole(data: Mapping[str, Any], key: str, where: str) -> int:
    """Return one whole number.

    Raises
    ------
    ValueError
        If the value is of another JSON type, true and false among them.
    """
    value = data[key]
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{where} {key} is {type(value).__name__}, not a whole number")
    return value


def number(data: Mapping[str, Any], key: str, where: str) -> float:
    """Return one number, whole or not.

    Raises
    ------
    ValueError
        If the value is of another JSON type, true and false among them.
    """
    value = data[key]
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ValueError(f"{where} {key} is {type(value).__name__}, not a number")
    return float(value)


def sequence(data: Mapping[str, Any], key: str, where: str) -> Sequence[Any]:
    """Return one list of values.

    Raises
    ------
    ValueError
        If the value is of another JSON type, a string among them.
    """
    value = data[key]
    if isinstance(value, str) or not isinstance(value, Sequence):
        raise ValueError(f"{where} {key} is {type(value).__name__}, not a list")
    return value


def named_file(file: Path, named: str, key: str, where: str) -> Path:
    """Resolve a path the file names against the file's own directory.

    Raises
    ------
    ValueError
        If nothing is there to read.
    """
    found = Path(named)
    resolved = found if found.is_absolute() else file.parent / found
    if not resolved.is_file():
        raise ValueError(f"{where} {key} is {named!r}, and {os.fspath(resolved)} is no file")
    return resolved

"""Read a JSON file someone writes by hand, and say what is wrong with it.

A missing key, a key nothing reads, and a value of another JSON type are each refused in words
naming the key and what belongs there. A caller reads one key at a time, or hands `reader` a
dataclass whose fields and annotations say what its object holds. The two word a key
fault alike, and a type fault differently: one names the type given, the other the value.

Nothing here knows what it reads for. Every refusal takes as `where` the subject it names, so a
caller supplies its own word for the file, or for the object inside it a key sits in.
"""

import json
import os
from collections.abc import Callable, Mapping, Sequence
from dataclasses import MISSING, fields, is_dataclass
from pathlib import Path
from types import NoneType, UnionType
from typing import Any, Literal, TypeAliasType, Union, get_args, get_origin, get_type_hints


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

    A missing key is named first, and a key nothing reads only once none is missing.

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

    Examples
    --------
    A refusal reads as `where` then the key, so `where` is the subject that owns it.

    >>> text({"host": 1}, "host", "a build")
    Traceback (most recent call last):
        ...
    ValueError: a build host is int, not a string
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


def listing(data: Mapping[str, Any], key: str, where: str) -> Sequence[Any]:
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


type _Convert = Callable[[Any, str], Any]


def reader[T](cls: type[T]) -> Callable[[Any, str], T]:
    """Return a reader from a JSON object to the dataclass `cls`, checking its keys.

    Each field takes its annotation's JSON type, and a field with a default may be left out,
    its keys refused in `refuse_keys`'s words. The reader takes the parsed JSON and the `where`
    its refusals name, each key named under it.

    Raises
    ------
    TypeError
        If a field's annotation has no JSON form.
    """
    spec = fields(cls)  # pyright: ignore[reportArgumentType]
    hints = get_type_hints(cls)
    nested = {f.name: _converter(hints[f.name]) for f in spec}
    required = frozenset(
        f.name for f in spec if f.default is MISSING and f.default_factory is MISSING
    )

    def convert(data: Any, where: str) -> T:
        if not isinstance(data, Mapping):
            raise _refused(where, "an object", data)
        refuse_keys(data, required, frozenset(nested) - required, where)
        return cls(**{key: nested[key](value, f"{where}.{key}") for key, value in data.items()})

    return convert


def _refused(where: str, expected: str, value: Any) -> ValueError:
    match value:
        case str():
            got = "a string"
        case list():
            got = "a list"
        case Mapping():
            got = "an object"
        case None | bool():
            got = json.dumps(value)
        case int() | float():
            got = repr(value)
        case _:
            got = type(value).__name__
    return ValueError(f"{where}: expected {expected}, got {got}")


def _converter(hint: Any) -> _Convert:
    """Return a converter from parsed JSON to the annotation `hint`, refusing another JSON type."""
    if isinstance(hint, TypeAliasType):
        return _converter(hint.__value__)
    origin, args = get_origin(hint), get_args(hint)
    if origin is Literal:
        return _converter(type(args[0]))
    if origin in (Union, UnionType) and len(args) == 2 and NoneType in args:
        (inner,) = (arg for arg in args if arg is not NoneType)
        return _or_null(_converter(inner))
    if origin is tuple and args[1:] == (...,):
        return _list(_converter(args[0]))
    if origin is tuple:
        return _fixed(tuple(_converter(arg) for arg in args))
    if origin is Mapping and args[0] is str:
        return _mapping(_converter(args[1]))
    if isinstance(hint, type) and is_dataclass(hint):
        return reader(hint)
    if hint in _SCALARS:
        return _scalar(*_SCALARS[hint])
    raise TypeError(f"a field has no JSON form: {hint!r}")


def _list(item: _Convert) -> _Convert:
    def convert(data: Any, where: str) -> tuple[Any, ...]:
        if not isinstance(data, list):
            raise _refused(where, "a list", data)
        return tuple(item(value, f"{where}[{i}]") for i, value in enumerate(data))

    return convert


def _fixed(items: tuple[_Convert, ...]) -> _Convert:
    """Return a converter to a tuple of a fixed length, such as a span's two numbers."""
    expected = f"a list of {len(items)}"

    def convert(data: Any, where: str) -> tuple[Any, ...]:
        if not isinstance(data, list):
            raise _refused(where, expected, data)
        if len(data) != len(items):
            raise ValueError(f"{where}: expected {expected}, got {len(data)}")
        return tuple(
            item(value, f"{where}[{i}]")
            for i, (item, value) in enumerate(zip(items, data, strict=True))
        )

    return convert


def _mapping(value: _Convert) -> _Convert:
    def convert(data: Any, where: str) -> dict[str, Any]:
        if not isinstance(data, Mapping):
            raise _refused(where, "an object", data)
        return {
            key: value(item, f"{where}[{json.dumps(key, ensure_ascii=False)}]")
            for key, item in data.items()
        }

    return convert


def _or_null(convert: _Convert) -> _Convert:
    return lambda data, where: None if data is None else convert(data, where)


def _scalar(expected: str, accepts: Callable[[Any], bool]) -> _Convert:
    def convert(data: Any, where: str) -> Any:
        if not accepts(data):
            raise _refused(where, expected, data)
        return data

    return convert


# `bool` is a subclass of `int` in Python, and true is not a number in JSON.
_SCALARS: dict[type, tuple[str, Callable[[Any], bool]]] = {
    str: ("a string", lambda value: isinstance(value, str)),
    bool: ("true or false", lambda value: isinstance(value, bool)),
    int: ("a whole number", lambda value: type(value) is int),
    float: ("a number", lambda value: type(value) in (int, float)),
}

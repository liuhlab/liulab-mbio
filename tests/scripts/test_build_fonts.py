"""The font build script's rules, on plain values rather than on a font or the network."""

import importlib.util
import io
import sys
import zipfile
from pathlib import Path
from types import ModuleType

import pytest

REPO = Path(__file__).resolve().parents[2]


def _script() -> ModuleType:
    spec = importlib.util.spec_from_file_location("build_fonts", REPO / "scripts/build_fonts.py")
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # Registered before it runs: `dataclasses` resolves annotations through `sys.modules`.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


build = _script()
SANS, BOLD, MONO = build.FACES


def test_an_archive_other_than_the_pinned_release_is_refused() -> None:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        for face in build.FACES:
            archive.writestr(f"{build.FOLDER}/ttf/{face.source}", b"not a font")
        archive.writestr(f"{build.FOLDER}/LICENSE", b"")

    with pytest.raises(ValueError, match="SHA-256"):
        build.unpack(buffer.getvalue())


def test_a_label_face_ships_what_the_font_has_of_its_blocks_less_what_is_never_drawn() -> None:
    # A soft hyphen and a byte-order mark are invisible, and a CJK character lies in no block.
    cmap = {0x41, 0xE9, 0xAD, 0x3B2, 0x2026, 0x2122, 0xFEFF, 0xFFFD, 0x4E2D}
    drawn = [chr(point) for point in (0x41, 0xE9, 0x3B2, 0x2026, 0x2122, 0xFFFD)]

    assert build.characters(SANS, cmap) == drawn
    assert build.characters(BOLD, cmap) == drawn


def test_the_mono_face_ships_printable_ascii_and_the_replacement_character() -> None:
    every = {point for block in build.LABELS for point in block}

    shipped = build.characters(MONO, every)

    assert shipped == [chr(point) for point in range(0x20, 0x7F)] + [chr(0xFFFD)]


def test_each_subset_is_renamed_away_from_every_name_the_licence_reserves() -> None:
    full_names = set()
    for face in build.FACES:
        records = build.names(face)
        for value in records.values():
            for reserved in (*build.RESERVED_NAMES, "DejaVu"):
                assert reserved not in value
        full_names.add(records[4])

    assert len(full_names) == len(build.FACES)
    assert build.names(BOLD)[6] == "LiulabMbioSans-Bold"


def _kern(covered: str, first: dict, second: dict, values: list[list[int]]):
    return build.ClassKerning(frozenset(covered), first, second, values)


GLYPHS = {"A": "A", "V": "V", "T": "T", ".": "period"}


def test_a_pair_is_kerned_by_the_classes_of_its_two_glyphs() -> None:
    lookup = [_kern("AT", {"T": 1}, {"V": 1, "period": 2}, [[0, -80, 0], [0, 0, -120]])]

    assert build.kerning(GLYPHS, [lookup]) == {"AV": -80, "T.": -120}


def test_the_first_subtable_covering_a_glyph_decides_even_at_zero() -> None:
    first = _kern("A", {}, {"V": 1}, [[0, 0]])
    second = _kern("A", {}, {"V": 1}, [[0, -80]])

    assert build.kerning(GLYPHS, [[first, second]]) == {}


def test_lookups_add_up() -> None:
    one = [_kern("A", {}, {"V": 1}, [[0, -80]])]
    other = [_kern("A", {}, {"V": 1, "T": 1}, [[0, 30]])]

    assert build.kerning(GLYPHS, [one, other]) == {"AV": -50, "AT": 30}

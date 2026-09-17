#!/usr/bin/env python3
"""Rebuild the shipped font tables and subsets from the DejaVu 2.37 release.

    pixi run build-fonts [--release PATH] [--out DIR]

For each face `liulab_mbio.plot.fonts` draws with, DejaVu Sans, Sans Bold and Sans Mono, this
writes a table of advance widths, kerning pairs and glyph outlines, which the package measures
and draws from, and a WOFF2 subset a page embeds. Each subset is renamed, as DejaVu's licence
asks of a modified font, and keeps only the `kern` feature, so a browser shapes what the table
measured. The licence is written beside them. `docs/research/fonts.md` sources the data.

fontTools and brotli come from the `fonts` environment `build-fonts` runs in, and are imported
only where a font is read or written, so the package gains neither.

Without `--release` the archive is downloaded. Either way it is refused unless its SHA-256 is
the pinned one. Tests read no font and never reach the network.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
import unicodedata
import urllib.request
import zipfile
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from liulab_mbio.plot.fonts import INVISIBLE

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "src/liulab_mbio/data/fonts"

RELEASE = (
    "https://github.com/dejavu-fonts/dejavu-fonts/releases/download/version_2_37/"
    "dejavu-fonts-ttf-2.37.zip"
)
SHA256 = "7576310b219e04159d35ff61dd4a4ec4cdba4f35c00e002a136f00e96a908b0a"

#: The folder the archive holds everything in.
FOLDER = "dejavu-fonts-ttf-2.37"

#: What a label may spell: printable ASCII, Latin-1, Greek, General Punctuation, the trade mark
#: sign and the replacement character.
LABELS = (
    range(0x20, 0x7F),
    range(0xA0, 0x100),
    range(0x370, 0x400),
    range(0x2000, 0x2070),
    range(0x2122, 0x2123),
    range(0xFFFD, 0xFFFE),
)

#: What bases and translations spell: printable ASCII, and the replacement character.
BASES = (range(0x20, 0x7F), range(0xFFFD, 0xFFFE))

#: Whose kerning a label is measured with. Every pair DejaVu kerns among the characters shipped
#: holds a Latin letter, and text around a Latin letter is shaped as Latin.
SCRIPT = "latn"

#: The name records a subset keeps: the copyright, names and version, and the licence.
NAME_IDS = (0, 1, 2, 3, 4, 5, 6, 13, 14)

#: Words DejaVu's licence forbids in the name of a modified font.
RESERVED_NAMES = ("Bitstream", "Vera", "Arev", "Tavmjong Bah")


@dataclass(frozen=True)
class Face:
    """One face shipped: the font it is cut from, the name it is renamed to, what it covers."""

    name: str
    source: str
    family: str
    style: str
    blocks: tuple[range, ...]


FACES = (
    Face("sans", "DejaVuSans.ttf", "Liulab Mbio Sans", "Regular", LABELS),
    Face("sans-bold", "DejaVuSans-Bold.ttf", "Liulab Mbio Sans", "Bold", LABELS),
    Face("sans-mono", "DejaVuSansMono.ttf", "Liulab Mbio Sans Mono", "Regular", BASES),
)


@dataclass(frozen=True)
class ClassKerning:
    """One subtable kerning glyphs by class, as OpenType's pair adjustment format 2 holds it.

    For a first glyph in `covered`, ``values[first class][second class]`` moves its advance. A
    glyph absent from `first` or `second` is in class 0.
    """

    covered: frozenset[str]
    first: Mapping[str, int]
    second: Mapping[str, int]
    values: Sequence[Sequence[int]]


def unpack(archive: bytes) -> dict[str, bytes]:
    """Read each face's font file, and the licence as ``"LICENSE"``, out of the release.

    Raises
    ------
    ValueError
        If the archive's SHA-256 is not the pinned release's.
    """
    digest = hashlib.sha256(archive).hexdigest()
    if digest != SHA256:
        raise ValueError(f"the archive's SHA-256 is {digest}, not DejaVu 2.37's {SHA256}")
    with zipfile.ZipFile(io.BytesIO(archive)) as release:
        files = {face.source: release.read(f"{FOLDER}/ttf/{face.source}") for face in FACES}
        files["LICENSE"] = release.read(f"{FOLDER}/LICENSE")
    return files


def characters(face: Face, cmap: Collection[int]) -> list[str]:
    """Return the characters a face ships: those of its blocks the font has, in code point order.

    A character `liulab_mbio.plot.fonts` drops as invisible is left out.
    """
    return [
        chr(point)
        for block in face.blocks
        for point in block
        if point in cmap and unicodedata.category(chr(point)) not in INVISIBLE
    ]


def names(face: Face) -> dict[int, str]:
    """Return the name records a subset is renamed with, by name ID.

    Examples
    --------
    >>> names(FACES[1])[4]
    'Liulab Mbio Sans Bold'
    """
    full = face.family if face.style == "Regular" else f"{face.family} {face.style}"
    postscript = face.family.replace(" ", "")
    if face.style != "Regular":
        postscript += f"-{face.style}"
    return {1: face.family, 2: face.style, 3: full, 4: full, 6: postscript}


def kerning(glyphs: Mapping[str, str], lookups: Sequence[Sequence[ClassKerning]]) -> dict[str, int]:
    """Return how far each pair of characters is kerned, keyed by the pair, in font units.

    `glyphs` maps each character to its glyph. The lookups add up; within one, the first
    subtable covering a pair's first glyph decides the pair, even at zero, as HarfBuzz applies a
    class subtable. A pair left at zero is left out.
    """
    pairs: dict[str, int] = {}
    for first, first_glyph in glyphs.items():
        deciding = [
            subtable
            for lookup in lookups
            if (subtable := next((s for s in lookup if first_glyph in s.covered), None))
        ]
        for second, second_glyph in glyphs.items():
            total = sum(
                subtable.values[subtable.first.get(first_glyph, 0)][
                    subtable.second.get(second_glyph, 0)
                ]
                for subtable in deciding
            )
            if total:
                pairs[first + second] = total
    return pairs


def read(face: Face, font_file: bytes) -> tuple[dict[str, Any], list[str]]:
    """Read one face's table from its font file, and say which characters it ships."""
    from fontTools.pens.svgPathPen import SVGPathPen  # pyright: ignore[reportMissingImports]
    from fontTools.ttLib import TTFont  # pyright: ignore[reportMissingImports]

    font = TTFont(io.BytesIO(font_file))
    cmap = font.getBestCmap()
    shipped = characters(face, cmap)
    glyphs = {char: cmap[ord(char)] for char in shipped}
    glyph_set = font.getGlyphSet()
    outlines: dict[str, str] = {}
    for char, glyph in glyphs.items():
        pen = SVGPathPen(glyph_set, ntos=lambda number: f"{number:g}")
        glyph_set[glyph].draw(pen)
        outlines[char] = pen.getCommands()
    table = {
        "family": face.family,
        "style": face.style,
        "source": {
            "release": "DejaVu 2.37",
            "url": RELEASE,
            "sha256": SHA256,
            "file": f"{FOLDER}/ttf/{face.source}",
            "licence": "LICENSE, beside this file",
        },
        "units_per_em": font["head"].unitsPerEm,
        "ascender": font["hhea"].ascent,
        "descender": font["hhea"].descent,
        "advances": {char: font["hmtx"][glyph][0] for char, glyph in glyphs.items()},
        "kerning": kerning(glyphs, _kern_lookups(font)),
        "outlines": outlines,
    }
    return table, shipped


def subset(face: Face, font_file: bytes, shipped: Sequence[str]) -> bytes:
    """Cut a font down to the characters shipped, rename it, and write it as WOFF2.

    Every layout feature but `kern` is dropped, and so is hinting, which moves no advance width.
    """
    from fontTools import subset as subsetting  # pyright: ignore[reportMissingImports]
    from fontTools.ttLib import TTFont  # pyright: ignore[reportMissingImports]

    font = TTFont(io.BytesIO(font_file), recalcTimestamp=False)
    options = subsetting.Options()
    options.layout_features = ["kern"]
    options.hinting = False
    options.name_IDs = list(NAME_IDS)
    options.drop_tables += ["FFTM", "GSUB", "MATH"]
    subsetter = subsetting.Subsetter(options)
    subsetter.populate(unicodes=[ord(char) for char in shipped])
    subsetter.subset(font)
    for name_id, value in names(face).items():
        font["name"].removeNames(nameID=name_id)
        font["name"].setName(value, name_id, 3, 1, 0x409)
    font.flavor = "woff2"
    buffer = io.BytesIO()
    font.save(buffer)
    return buffer.getvalue()


def _kern_lookups(font: Any) -> list[list[ClassKerning]]:
    """Read the kerning lookups `SCRIPT` applies, refusing any this script does not read."""
    if "GPOS" not in font:
        return []
    gpos = font["GPOS"].table
    scripts = [record for record in gpos.ScriptList.ScriptRecord if record.ScriptTag == SCRIPT]
    if not scripts:
        return []
    features = [
        gpos.FeatureList.FeatureRecord[index].Feature
        for index in scripts[0].Script.DefaultLangSys.FeatureIndex
        if gpos.FeatureList.FeatureRecord[index].FeatureTag == "kern"
    ]
    lookups: list[list[ClassKerning]] = []
    for feature in features:
        for index in feature.LookupListIndex:
            lookup = gpos.LookupList.Lookup[index]
            subtables = []
            for pairs in lookup.SubTable:
                shape = (lookup.LookupType, pairs.Format, pairs.ValueFormat1, pairs.ValueFormat2)
                if shape != (2, 2, 4, 0):
                    raise SystemExit(f"kerning lookup {index} is {shape}, which is not read here")
                subtables.append(
                    ClassKerning(
                        covered=frozenset(pairs.Coverage.glyphs),
                        first=dict(pairs.ClassDef1.classDefs),
                        second=dict(pairs.ClassDef2.classDefs),
                        values=[
                            [getattr(cell.Value1, "XAdvance", 0) for cell in row.Class2Record]
                            for row in pairs.Class1Record
                        ],
                    )
                )
            lookups.append(subtables)
    return lookups


def fetch() -> bytes:
    """Download the release archive."""
    with urllib.request.urlopen(RELEASE, timeout=300) as response:
        return response.read()


def main(argv: list[str] | None = None) -> int:
    """Read the release, and write each face's table and subset and the licence."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", type=Path, help="the release archive, already downloaded")
    parser.add_argument("--out", type=Path, default=DATA, help="the directory to write into")
    arguments = parser.parse_args(argv)

    archive = arguments.release.read_bytes() if arguments.release else fetch()
    try:
        files = unpack(archive)
    except ValueError as error:
        raise SystemExit(str(error)) from None
    arguments.out.mkdir(parents=True, exist_ok=True)
    for face in FACES:
        table, shipped = read(face, files[face.source])
        (arguments.out / f"{face.name}.json").write_text(json.dumps(table, indent=1) + "\n")
        (arguments.out / f"{face.name}.woff2").write_bytes(
            subset(face, files[face.source], shipped)
        )
    (arguments.out / "LICENSE").write_bytes(files["LICENSE"])
    print(f"wrote {arguments.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

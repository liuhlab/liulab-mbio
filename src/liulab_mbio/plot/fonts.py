"""Measure and draw text in the three pinned faces, from tables shipped as package data.

Labels are set in DejaVu Sans, regular or bold, and bases and translations in DejaVu Sans Mono.
`scripts/build_fonts.py` builds each face's table of advance widths, kerning pairs and glyph
outlines, and its WOFF2 subset, from the DejaVu 2.37 release; the licence ships beside them, and
`docs/research/fonts.md` sources them. Nothing here reads a font file.

Measuring and drawing share one rule, so a label drawn from `Font.letters` fills the box
`Font.width` laid out: characters in an `INVISIBLE` category are dropped, and any other
character the face lacks is drawn as `REPLACEMENT`.
"""

import json
import unicodedata
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from importlib.resources import files
from typing import NamedTuple

#: The Unicode categories dropped rather than drawn: control and format characters, such as a
#: byte-order mark, and line and paragraph separators.
INVISIBLE = frozenset({"Cc", "Cf", "Zl", "Zp"})

#: What a character the face lacks is drawn as.
REPLACEMENT = "\ufffd"


class Letter(NamedTuple):
    """One letter of a text as a face draws it.

    `x` is where the letter starts from the start of the text and `advance` how far its own
    glyph moves on, both in the unit the size was given in. Kerning moves `x` and never
    `advance`, so the text is as wide as its last letter's `x` plus `advance`. `outline` is SVG
    path data in font units, `Font.units_per_em` to the em, with y pointing up from the baseline;
    a space has none.
    """

    char: str
    x: float
    advance: float
    outline: str


@dataclass(frozen=True, slots=True)
class _Table:
    family: str
    style: str
    units_per_em: int
    ascender: int
    descender: int
    advances: Mapping[str, int]
    kerning: Mapping[str, int]
    outlines: Mapping[str, str]


@dataclass(frozen=True, slots=True)
class Font:
    r"""One pinned face, read from package data the first time it is used.

    Parameters
    ----------
    name
        The name its table and subset ship under: ``"sans"``, ``"sans-bold"`` or
        ``"sans-mono"``.

    Examples
    --------
    Kerning draws ``V`` in under ``A``:

    >>> [round(letter.x, 2) for letter in SANS.letters("AV", 10)]
    [0.0, 6.2]
    >>> round(SANS.width("AV", 10), 2), round(2 * SANS.width("A", 10), 2)
    (13.04, 13.68)

    A byte-order mark is dropped, and a character the face lacks is drawn as ``REPLACEMENT``:

    >>> SANS.drawn("\ufeffΔlacZ 中")
    'ΔlacZ �'
    >>> MONO.drawn("ΔlacZ")
    '�lacZ'
    """

    name: str

    @property
    def family(self) -> str:
        """The family name the subset carries, which is not DejaVu's."""
        return _table(self.name).family

    @property
    def style(self) -> str:
        """``"Regular"`` or ``"Bold"``: the style the subset carries within its family."""
        return _table(self.name).style

    @property
    def units_per_em(self) -> int:
        """Font units to the em, the unit of `ascender`, `descender` and each outline."""
        return _table(self.name).units_per_em

    @property
    def ascender(self) -> int:
        """How far the face reaches above the baseline, in font units."""
        return _table(self.name).ascender

    @property
    def descender(self) -> int:
        """How far the face reaches below the baseline, in font units, as a negative number."""
        return _table(self.name).descender

    def drawn(self, text: str) -> str:
        """Return `text` as this face draws it, which is the text a page should hold.

        Invisible characters are dropped, and any other the face lacks becomes `REPLACEMENT`.
        """
        advances = _table(self.name).advances
        return "".join(
            char if char in advances else REPLACEMENT
            for char in text
            if unicodedata.category(char) not in INVISIBLE
        )

    def width(self, text: str, size: float) -> float:
        """Return how wide `text` is at `size`, kerning applied, in the unit of `size`."""
        letters = self.letters(text, size)
        return letters[-1].x + letters[-1].advance if letters else 0.0

    def letters(self, text: str, size: float) -> tuple[Letter, ...]:
        """Return each letter `text` is drawn with at `size`, placed as `width` measures it."""
        table = _table(self.name)
        scale = size / table.units_per_em
        drawn = self.drawn(text)
        letters: list[Letter] = []
        x = 0
        for index, char in enumerate(drawn):
            if index:
                x += table.kerning.get(drawn[index - 1] + char, 0)
            advance = table.advances[char]
            letters.append(Letter(char, x * scale, advance * scale, table.outlines[char]))
            x += advance
        return tuple(letters)

    def woff2(self) -> bytes:
        """Return the face's WOFF2 subset, for a page to embed."""
        return (files("liulab_mbio") / "data" / "fonts" / f"{self.name}.woff2").read_bytes()


#: DejaVu Sans, for labels.
SANS = Font("sans")

#: DejaVu Sans Bold, for the record's name and an enzyme that cuts once.
BOLD = Font("sans-bold")

#: DejaVu Sans Mono, for bases and translations.
MONO = Font("sans-mono")


@cache
def _table(name: str) -> _Table:
    data = json.loads((files("liulab_mbio") / "data" / "fonts" / f"{name}.json").read_text())
    return _Table(
        family=data["family"],
        style=data["style"],
        units_per_em=data["units_per_em"],
        ascender=data["ascender"],
        descender=data["descender"],
        advances=data["advances"],
        kerning=data["kerning"],
        outlines=data["outlines"],
    )

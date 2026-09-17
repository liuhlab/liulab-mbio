"""The shapes a view is drawn in, and the SVG they are written as.

A view computes every shape in points, x to the right and y down; this module knows no view. Text
is written as text in the face it was measured in, pinned to its measured width with
`textLength`, or letter by letter at the places it was laid out, so a browser draws it where it
was measured and it stays text to search and copy. Written as outlines instead, each letter is its
glyph's outline from the font tables, drawn at its place, so the SVG draws the same with no font
to look up.
"""

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from functools import lru_cache
from html import escape
from typing import NamedTuple

from liulab_mbio.plot.fonts import Font, Letter
from liulab_mbio.plot.labels import Box


@dataclass(frozen=True, slots=True)
class Path:
    """An outline in SVG path data, filled and stroked."""

    d: str
    fill: str
    stroke: str
    width: float = 1.0


@dataclass(frozen=True, slots=True)
class Line:
    """A straight line."""

    x1: float
    y1: float
    x2: float
    y2: float
    stroke: str
    width: float = 1.0


@dataclass(frozen=True, slots=True)
class Circle:
    """A circle, stroked and never filled."""

    cx: float
    cy: float
    r: float
    stroke: str
    width: float = 1.0


@dataclass(frozen=True, slots=True)
class Rect:
    """A box with rounded corners, filled and stroked."""

    box: Box
    fill: str
    stroke: str
    width: float = 1.0
    corner: float = 0.0


@dataclass(frozen=True, slots=True)
class Text:
    """A line of text, `x` where it starts and `y` its baseline."""

    x: float
    y: float
    text: str
    font: Font
    size: float
    fill: str = "#000000"


class Place(NamedTuple):
    """Where a letter starts on its baseline, and how many degrees clockwise it turns about it."""

    x: float
    y: float
    rotate: float


@dataclass(frozen=True, slots=True)
class Letters:
    """A line of text set letter by letter, each letter at its own place, as a name on an arc is.

    Parameters
    ----------
    places
        One for each letter ``font.letters(text, size)`` draws, in order.
    text, font, size, fill
        As a `Text` has them.
    """

    places: tuple[Place, ...]
    text: str
    font: Font
    size: float
    fill: str = "#000000"


@dataclass(frozen=True, slots=True)
class Group:
    """Shapes drawn together.

    Parameters
    ----------
    shapes
        In the order they are drawn.
    classes
        Names a page styles or switches the group by.
    data
        What a page reads from the group, such as its hover details, each written as a
        ``data-`` attribute.
    rotate
        Degrees clockwise about the centre (0, 0), if the group is turned.
    x, y
        Where the group's (0, 0) is drawn, if the group is moved; it turns before it moves.
    """

    shapes: tuple["Shape", ...]
    classes: tuple[str, ...] = ()
    data: Mapping[str, str] = field(default_factory=dict, hash=False)
    rotate: float = 0.0
    x: float = 0.0
    y: float = 0.0


type Shape = Path | Line | Circle | Rect | Text | Letters | Group


def document(shapes: Iterable[Shape], extent: Box, *, outlines: bool = False) -> str:
    """Return an SVG element drawing `shapes`, its view box `extent`.

    Text is written as text, or with `outlines` as the outlines of its letters, each glyph defined
    once and used wherever it is drawn.
    """
    view = " ".join(number(value) for value in (extent.x, extent.y, extent.width, extent.height))
    glyphs: dict[str, str] | None = {} if outlines else None
    body = "".join(_shape(shape, glyphs) for shape in shapes)
    if glyphs:
        defs = "".join(f'<path id="{key}" d="{outline}"/>' for key, outline in glyphs.items())
        body = f"<defs>{defs}</defs>{body}"
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{view}" '
        f'width="{number(extent.width)}" height="{number(extent.height)}">{body}</svg>'
    )


def _shape(shape: Shape, glyphs: dict[str, str] | None) -> str:
    """Return a shape as SVG, its text as outlines when `glyphs` collects the glyphs they use."""
    match shape:
        case Path(d, fill, stroke, width):
            return f'<path d="{d}"{_paint(fill, stroke, width)}/>'
        case Line(x1, y1, x2, y2, stroke, width):
            points = _attributes(x1=x1, y1=y1, x2=x2, y2=y2)
            return f"<line{points}{_paint('none', stroke, width)}/>"
        case Circle(cx, cy, r, stroke, width):
            return f"<circle{_attributes(cx=cx, cy=cy, r=r)}{_paint('none', stroke, width)}/>"
        case Rect(box, fill, stroke, width, corner):
            place = _attributes(x=box.x, y=box.y, width=box.width, height=box.height, rx=corner)
            return f"<rect{place}{_paint(fill, stroke, width)}/>"
        case Text():
            return _text(shape) if glyphs is None else _outlined(shape, glyphs)
        case Letters():
            return _letters(shape) if glyphs is None else _outlined(shape, glyphs)
        case Group(shapes, classes, data, rotate, x, y):
            attributes = f' class="{escape(" ".join(classes))}"' if classes else ""
            attributes += "".join(
                f' data-{escape(key)}="{escape(value)}"' for key, value in data.items()
            )
            moves = [f"translate({number(x)} {number(y)})"] if x or y else []
            moves += [f"rotate({number(rotate)})"] if rotate else []
            if moves:
                attributes += f' transform="{" ".join(moves)}"'
            return f"<g{attributes}>{''.join(_shape(one, glyphs) for one in shapes)}</g>"


def _text(text: Text) -> str:
    weight = "700" if text.font.style == "Bold" else "400"
    attributes = _attributes(
        x=text.x, y=text.y, font_size=text.size, textLength=text.font.width(text.text, text.size)
    )
    return (
        f'<text{attributes} font-family="\'{escape(text.font.family)}\'" font-weight="{weight}"'
        f' fill="{escape(text.fill)}">{escape(text.font.drawn(text.text))}</text>'
    )


def _letters(letters: Letters) -> str:
    """Write one text element placing each letter, so the whole line stays text in a page.

    A baseline the letters share is written once, and their turns only where one turns: a letter
    given no baseline keeps the one before it, and one given no turn is not turned.
    """
    xs = [place.x for place in letters.places]
    ys = {place.y for place in letters.places}
    places = f' x="{" ".join(map(number, xs))}"'
    if len(ys) == 1:
        places += f' y="{number(ys.pop())}"'
    else:
        places += f' y="{" ".join(number(place.y) for place in letters.places)}"'
    if any(place.rotate for place in letters.places):
        places += f' rotate="{" ".join(number(place.rotate) for place in letters.places)}"'
    weight = "700" if letters.font.style == "Bold" else "400"
    return (
        f"<text{places}{_attributes(font_size=letters.size)}"
        f' font-family="\'{escape(letters.font.family)}\'" font-weight="{weight}"'
        f' fill="{escape(letters.fill)}">{escape(letters.font.drawn(letters.text))}</text>'
    )


def _outlined(line: Text | Letters, glyphs: dict[str, str]) -> str:
    """Return each letter of `line` as a use of its glyph's outline, adding the glyph to `glyphs`.

    A glyph keeps the tables' path data, in font units with y up. Each use scales it to the line's
    size and flips it, then puts it where it was measured: along the baseline of a `Text`, or at a
    `Letters` place, turned about it.
    """
    font = line.font
    scale = line.size / font.units_per_em
    uses: list[str] = []
    match line:
        case Text(x, y):
            for letter in font.letters(line.text, line.size):
                if letter.outline:
                    at = f'x="{number(letter.x / scale)}"'
                    uses.append(f'<use href="#{_glyph(font, letter, glyphs)}" {at}/>')
            place = f' transform="{_matrix(scale, 0.0, x, y)}"'
        case Letters(places):
            for letter, (x, y, rotate) in zip(
                font.letters(line.text, line.size), places, strict=True
            ):
                if letter.outline:
                    at = f'transform="{_matrix(scale, rotate, x, y)}"'
                    uses.append(f'<use href="#{_glyph(font, letter, glyphs)}" {at}/>')
            place = ""
    if not uses:
        return ""
    return f'<g{place} fill="{escape(line.fill)}">{"".join(uses)}</g>'


def _glyph(font: Font, letter: Letter, glyphs: dict[str, str]) -> str:
    """Return the id a letter's glyph is defined under, adding it to `glyphs`."""
    key = f"{font.name}-{ord(letter.char):x}"
    glyphs[key] = letter.outline
    return key


def _matrix(scale: float, rotate: float, x: float, y: float) -> str:
    """Return a transform from font units to a place: scaled, y flipped, turned `rotate` degrees."""
    cos, sin = math.cos(math.radians(rotate)), math.sin(math.radians(rotate))
    a, b = scale * cos, scale * sin
    return f"matrix({a:.9g} {b:.9g} {b:.9g} {-a:.9g} {number(x)} {number(y)})"


def _paint(fill: str, stroke: str, width: float) -> str:
    return f' fill="{escape(fill)}" stroke="{escape(stroke)}" stroke-width="{number(width)}"'


def _attributes(**values: float) -> str:
    return "".join(f' {name.replace("_", "-")}="{number(value)}"' for name, value in values.items())


# A drawing repeats its coordinates row after row.
@lru_cache(maxsize=1 << 16)
def number(value: float) -> str:
    """Return a coordinate as the SVG writes it: to a hundredth of a point, no trailing zeros."""
    text = f"{value:.2f}".rstrip("0").rstrip(".")
    return "0" if text == "-0" else text

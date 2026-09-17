"""The shapes a view is drawn in, and the SVG they are written as.

A view computes every shape in points, x to the right and y down; this module knows no view. Text
is written as text in the face it was measured in, pinned to its measured width with
`textLength`, so a browser draws it the width it was laid out.
"""

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from html import escape

from liulab_mbio.plot.fonts import Font
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
    """

    shapes: tuple["Shape", ...]
    classes: tuple[str, ...] = ()
    data: Mapping[str, str] = field(default_factory=dict, hash=False)
    rotate: float = 0.0


type Shape = Path | Line | Circle | Rect | Text | Group


def document(shapes: Iterable[Shape], extent: Box) -> str:
    """Return an SVG element drawing `shapes`, its view box `extent`."""
    view = " ".join(number(value) for value in (extent.x, extent.y, extent.width, extent.height))
    body = "".join(_shape(shape) for shape in shapes)
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{view}" '
        f'width="{number(extent.width)}" height="{number(extent.height)}">{body}</svg>'
    )


def _shape(shape: Shape) -> str:
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
            return _text(shape)
        case Group(shapes, classes, data, rotate):
            attributes = f' class="{escape(" ".join(classes))}"' if classes else ""
            attributes += "".join(
                f' data-{escape(key)}="{escape(value)}"' for key, value in data.items()
            )
            if rotate:
                attributes += f' transform="rotate({number(rotate)})"'
            return f"<g{attributes}>{''.join(_shape(one) for one in shapes)}</g>"


def _text(text: Text) -> str:
    weight = "700" if text.font.style == "Bold" else "400"
    attributes = _attributes(
        x=text.x, y=text.y, font_size=text.size, textLength=text.font.width(text.text, text.size)
    )
    return (
        f'<text{attributes} font-family="\'{escape(text.font.family)}\'" font-weight="{weight}"'
        f' fill="{escape(text.fill)}">{escape(text.font.drawn(text.text))}</text>'
    )


def _paint(fill: str, stroke: str, width: float) -> str:
    return f' fill="{escape(fill)}" stroke="{escape(stroke)}" stroke-width="{number(width)}"'


def _attributes(**values: float) -> str:
    return "".join(f' {name.replace("_", "-")}="{number(value)}"' for name, value in values.items())


def number(value: float) -> str:
    """Return a coordinate as the SVG writes it: to a hundredth of a point, no trailing zeros."""
    text = f"{value:.2f}".rstrip("0").rstrip(".")
    return "0" if text == "-0" else text

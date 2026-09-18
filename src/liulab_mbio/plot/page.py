"""The HTML page a drawing is written as: one file that loads nothing over the network.

The page holds the drawing's SVG inline, its styles and script beside this module, and the three
font subsets the drawing was measured in, so text draws at its laid-out width and stays text to
search and copy. It keeps a white background whatever the browser's colour scheme.

Switches above the drawing show or hide each kind of item and each feature type in place, flip
the map between its shapes at top right, and show the sequence view with it: beside the map while
the page is wide enough for both, under it otherwise. A caption on the map zooms a shape that
zooms, and a drag pans it. A click on an item highlights it in both views, and scrolls the other
view to it. In the sequence view, a toggle shows one strand or both, hovering over a base shows
its position, and a drag selects bases to copy, scrolling the view while it passes the view's top
or bottom.
"""

import base64
import math
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from html import escape
from importlib.resources import files
from typing import Literal

from liulab_mbio.plot.fonts import BOLD, MONO, SANS, Font

#: How many times its size a shape that zooms is drawn at most.
_ZOOM = 8


@dataclass(frozen=True, slots=True)
class Switch:
    """A switch showing or hiding in place every item of one kind, or every feature of one type.

    Parameters
    ----------
    by
        What it matches: an item's ``data-kind``, or a feature's ``data-type``.
    value
        The kind or type it switches.
    text
        What it says.
    on
        Whether it starts on.
    """

    by: Literal["kind", "type"]
    value: str
    text: str
    on: bool


def render(
    maps: Mapping[str, str],
    *,
    title: str,
    shown: str,
    zooms: Collection[str] = (),
    switches: Sequence[Switch] = (),
    sequence_view: str | None = None,
    sequence_shown: bool = False,
    both_strands: bool = True,
) -> str:
    """Return a page showing a map under `title`, with `sequence_view`.

    `maps` holds the map as an SVG element in each shape it is drawn in, such as ``"circle"``
    and ``"line"``. The shape `shown` names shows first, and a switch flips between them when
    there are two. Each shape `zooms` names zooms in up to eight times by scaling its drawing,
    from a caption on the map that shows while it does. `sequence_view` is an SVG element too, as
    `sequence_view.layout` draws one with both strands, or ``None`` for the map alone. A switch
    shows it, first when `sequence_shown`, and a toggle in it shows its bottom strand, first when
    `both_strands`.
    """
    fonts = "".join(_font_face(font) for font in (SANS, BOLD, MONO))
    shapes = "".join(
        _shape(shape, image, shape == shown, shape in zooms) for shape, image in maps.items()
    )
    caption = _zoom(shown in zooms) if zooms else ""
    figures = f'<figure class="map">{caption}\n{shapes}</figure>\n'
    if sequence_view is not None:
        figures += _sequence_view(sequence_view, sequence_shown, both_strands)
    return (
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<meta name="color-scheme" content="light">\n'
        f"<title>{escape(title)}</title>\n<style>\n{fonts}{_asset('plot.css')}</style>\n"
        f"</head>\n<body>\n{_switches(switches, list(maps), shown, sequence_view, sequence_shown)}"
        f'<main class="plot">\n{figures}</main>\n'
        '<div class="hover" role="tooltip" hidden></div>\n'
        f"<script>\n{_asset('plot.js')}</script>\n</body>\n</html>\n"
    )


def _switches(
    switches: Sequence[Switch],
    shapes: Sequence[str],
    shown: str,
    sequence_view: str | None,
    sequence_shown: bool,
) -> str:
    """Return the bar of switches: kinds, then feature types, then the views at the right."""
    sets = [
        _fieldset(
            f"{by}s",
            words,
            [
                _input("checkbox", by, one.value, one.text, one.on)
                for one in switches
                if one.by == by
            ],
        )
        for by, words in (("kind", "Layers"), ("type", "Feature types"))
    ]
    inputs = []
    if len(shapes) > 1:
        inputs = [
            _input("radio", "shape", shape, shape.capitalize(), shape == shown) for shape in shapes
        ]
    if sequence_view is not None:
        inputs.append(_input("checkbox", "view", "sequence", "Sequence", sequence_shown))
    sets.append(_fieldset("shapes", "Views", inputs))
    body = "".join(sets)
    return f'<form class="switches" autocomplete="off">\n{body}</form>\n' if body else ""


def _sequence_view(image: str, shown: bool, both_strands: bool) -> str:
    """Return the sequence view's figure: its toggle, what is selected and a button to copy it."""
    classes = "sequence-view" if both_strands else "sequence-view one-strand"
    return (
        f'<figure class="{classes}"{"" if shown else " hidden"}>'
        f"<figcaption>{_input('checkbox', 'strands', 'both', 'Both strands', both_strands)}"
        '<output class="selection"></output><button type="button" class="copy" hidden>Copy'
        f"</button></figcaption>{image}</figure>\n"
    )


def _shape(shape: str, image: str, shown: bool, zooms: bool) -> str:
    zoom = ' data-zoom="free"' if zooms else ""
    hidden = "" if shown else " hidden"
    return f'<div class="shape" data-shape="{escape(shape)}"{zoom}{hidden}>{image}</div>\n'


def _zoom(shown: bool) -> str:
    """Return the map's caption: a zoom slider counting doublings, and a button to reset it."""
    slider = (
        f'<input type="range" name="zoom" min="0" max="{math.log2(_ZOOM):g}" step="any" value="0">'
    )
    return (
        f"<figcaption{'' if shown else ' hidden'}><label>Zoom {slider}</label>"
        '<button type="button" class="reset">Reset</button></figcaption>'
    )


def _fieldset(name: str, words: str, inputs: Sequence[str]) -> str:
    if not inputs:
        return ""
    return f'<fieldset class="{name}" aria-label="{words}">{"".join(inputs)}</fieldset>\n'


def _input(kind: str, name: str, value: str, text: str, on: bool) -> str:
    checked = " checked" if on else ""
    return (
        f'<label><input type="{kind}" name="{name}" value="{escape(value)}"{checked}>'
        f"{escape(text)}</label>"
    )


@cache
def _font_face(font: Font) -> str:
    data = base64.b64encode(font.woff2()).decode("ascii")
    weight = 700 if font.style == "Bold" else 400
    return (
        f"@font-face {{ font-family: '{font.family}'; font-weight: {weight}; "
        f"src: url(data:font/woff2;base64,{data}) format('woff2'); }}\n"
    )


@cache
def _asset(name: str) -> str:
    return files("liulab_mbio.plot").joinpath(name).read_text(encoding="utf-8")

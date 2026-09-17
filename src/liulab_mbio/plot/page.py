"""The HTML page a drawing is written as: one file that loads nothing over the network.

The page holds the drawing's SVG inline, its styles and script beside this module, and the three
font subsets the drawing was measured in, so text draws at its laid-out width and stays text to
search and copy. It keeps a white background whatever the browser's colour scheme.

Switches above the drawing show or hide each kind of item and each feature type in place, and flip
the map between its shapes at top right.
"""

import base64
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from html import escape
from importlib.resources import files
from typing import Literal

from liulab_mbio.plot.fonts import BOLD, MONO, SANS, Font


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
    switches: Sequence[Switch] = (),
    sequence_view: str | None = None,
) -> str:
    """Return a page showing a map under `title`, with `sequence_view` beside it.

    `maps` holds the map as an SVG element in each shape it is drawn in, such as ``"circle"``
    and ``"line"``. The shape `shown` names shows first, and a switch flips between them when
    there are two. `sequence_view` is an SVG element too, or ``None`` for the map alone.
    """
    fonts = "".join(_font_face(font) for font in (SANS, BOLD, MONO))
    shapes = "".join(
        f'<div class="shape" data-shape="{escape(shape)}"{"" if shape == shown else " hidden"}>'
        f"{image}</div>\n"
        for shape, image in maps.items()
    )
    figures = f'<figure class="map">\n{shapes}</figure>\n'
    if sequence_view is not None:
        figures += f'<figure class="sequence-view">{sequence_view}</figure>\n'
    return (
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<meta name="color-scheme" content="light">\n'
        f"<title>{escape(title)}</title>\n<style>\n{fonts}{_asset('plot.css')}</style>\n"
        f"</head>\n<body>\n{_switches(switches, list(maps), shown)}"
        f'<main class="plot">\n{figures}</main>\n'
        '<div class="hover" role="tooltip" hidden></div>\n'
        f"<script>\n{_asset('plot.js')}</script>\n</body>\n</html>\n"
    )


def _switches(switches: Sequence[Switch], shapes: Sequence[str], shown: str) -> str:
    """Return the bar of switches: kinds, then feature types, then the shapes at the right."""
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
    if len(shapes) > 1:
        inputs = [
            _input("radio", "shape", shape, shape.capitalize(), shape == shown) for shape in shapes
        ]
        sets.append(_fieldset("shapes", "Shape", inputs))
    body = "".join(sets)
    return f'<form class="switches" autocomplete="off">\n{body}</form>\n' if body else ""


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


def _font_face(font: Font) -> str:
    data = base64.b64encode(font.woff2()).decode("ascii")
    weight = 700 if font.style == "Bold" else 400
    return (
        f"@font-face {{ font-family: '{font.family}'; font-weight: {weight}; "
        f"src: url(data:font/woff2;base64,{data}) format('woff2'); }}\n"
    )


def _asset(name: str) -> str:
    return files("liulab_mbio.plot").joinpath(name).read_text(encoding="utf-8")

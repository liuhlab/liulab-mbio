"""The HTML page a drawing is written as: one file that loads nothing over the network.

The page holds the drawing's SVG inline, its styles and script beside this module, and the three
font subsets the drawing was measured in, so text draws at its laid-out width and stays text to
search and copy. It keeps a white background whatever the browser's colour scheme.
"""

import base64
from html import escape
from importlib.resources import files

from liulab_mbio.plot.fonts import BOLD, MONO, SANS, Font


def render(drawing: str, *, title: str, sequence_view: str | None = None) -> str:
    """Return a page showing `drawing`, an SVG element, under `title`, with `sequence_view` beside.

    `sequence_view` is an SVG element too, or ``None`` for the map alone.
    """
    fonts = "".join(_font_face(font) for font in (SANS, BOLD, MONO))
    figures = f'<figure class="map">{drawing}</figure>\n'
    if sequence_view is not None:
        figures += f'<figure class="sequence-view">{sequence_view}</figure>\n'
    return (
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<meta name="color-scheme" content="light">\n'
        f"<title>{escape(title)}</title>\n<style>\n{fonts}{_asset('plot.css')}</style>\n"
        f'</head>\n<body>\n<main class="plot">\n{figures}</main>\n'
        '<div class="hover" role="tooltip" hidden></div>\n'
        f"<script>\n{_asset('plot.js')}</script>\n</body>\n</html>\n"
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

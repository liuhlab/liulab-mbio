"""Turn a drawing's SVG into a PNG or a PDF.

vl-convert draws the SVG, and pypdf joins a PDF's pages; each is imported inside the function that
uses it. Hand these functions an SVG with its letters as outlines, `svg.document` with `outlines`:
vl-convert draws nothing for text in a font it has not been given, and raises no error.
"""

import io
from collections.abc import Sequence

#: The longest side, in pixels, a PNG is drawn with. vl-convert's renderer places shapes in 32-bit
#: floating point, which holds a position to the quarter pixel it anti-aliases at only below this
#: size; past it, shapes drift from where they were laid out.
PNG_SIDE = 2**22

#: The most pixels a PNG is drawn with. The renderer holds every pixel in memory, and more while it
#: writes the file, so past this a PNG could take more memory than a laptop has to spare.
PNG_PIXELS = 2**28

#: Points to the inch. A drawing is laid out in points.
POINTS_PER_INCH = 72


def png(document: str, *, dpi: float) -> bytes:
    """Return the SVG `document`, measured in points, drawn as a PNG at `dpi`.

    The PNG records `dpi`, so it prints at the size it was laid out. The caller keeps each side
    within `PNG_SIDE` pixels.
    """
    import vl_convert

    # vl-convert's linux-64 conda build ships no type stub, so pyright there sees no functions.
    return vl_convert.svg_to_png(document, ppi=dpi)  # pyright: ignore[reportAttributeAccessIssue]


def pdf(pages: Sequence[str]) -> bytes:
    """Return a PDF with each SVG document of `pages`, measured in points, on a page of its own."""
    import pypdf
    import vl_convert

    joined = pypdf.PdfWriter()
    for page in pages:
        # No type stub on linux-64, as in `png`.
        joined.append(pypdf.PdfReader(io.BytesIO(vl_convert.svg_to_pdf(page))))  # pyright: ignore[reportAttributeAccessIssue]
    out = io.BytesIO()
    joined.write(out)
    return out.getvalue()

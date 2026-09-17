"""PROTOTYPE — throwaway. Ticket #80: which converter turns our SVG into PNG and PDF.

Asks two things of CairoSVG and vl-convert:
1. Font: does the converter draw with the one pinned font file the layout measured?
2. Scale: how long does a 100 kb sequence view take as a paged PDF, and how tall a PNG works?

Run: pixi exec -s python=3.13 -s "vl-convert-python=1.9.*" -s cairosvg -s fonttools \
       -s pypdf -s pillow -s numpy -s uharfbuzz -s fontconfig -- python measure.py
"""

import base64
import io
import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
FONTDIR = HERE / "fontdir"
PINNED = FONTDIR / "PinnedSans.ttf"
OUT = HERE / "out"
FAMILY = "PinnedSans"
TEXT = "KpnI - Acc65I (1234) lacZα M13 fwd AVATAR"
SIZE = 40.0


# ---------------------------------------------------------------- prep


def prep():
    from fontTools.ttLib import TTFont

    FONTDIR.mkdir(exist_ok=True)
    OUT.mkdir(exist_ok=True)
    source = Path(os.environ["CONDA_PREFIX"]) / "fonts" / "DejaVuSans.ttf"
    font = TTFont(source)
    for rec in font["name"].names:
        if rec.nameID in (1, 4, 6, 16):
            rec.string = FAMILY
    font.save(PINNED)
    conf = HERE / "fonts.conf"
    conf.write_text(
        '<?xml version="1.0"?><!DOCTYPE fontconfig SYSTEM "fonts.dtd"><fontconfig>'
        f'<include ignore_missing="yes">{os.environ["CONDA_PREFIX"]}/etc/fonts/fonts.conf</include>'
        f"<dir>{FONTDIR}</dir><cachedir>{HERE / 'fccache'}</cachedir></fontconfig>"
    )


# ---------------------------------------------------------------- font test


def outline_svg(kern: bool) -> str:
    """The text as glyph outlines from the pinned file, shaped by HarfBuzz or not."""
    import uharfbuzz as hb
    from fontTools.pens.svgPathPen import SVGPathPen
    from fontTools.pens.transformPen import TransformPen
    from fontTools.ttLib import TTFont

    font = TTFont(PINNED)
    upem = font["head"].unitsPerEm
    glyphs = font.getGlyphSet()
    scale = SIZE / upem
    blob = hb.Blob.from_file_path(str(PINNED))
    hbfont = hb.Font(hb.Face(blob))
    buf = hb.Buffer()
    buf.add_str(TEXT)
    buf.guess_segment_properties()
    hb.shape(hbfont, buf, {"kern": kern})
    order = font.getGlyphOrder()
    x = 10.0
    paths = []
    for info, pos in zip(buf.glyph_infos, buf.glyph_positions, strict=True):
        name = order[info.codepoint]
        pen = SVGPathPen(glyphs)
        glyphs[name].draw(TransformPen(pen, (scale, 0, 0, -scale, x + pos.x_offset * scale, 55)))
        paths.append(f'<path d="{pen.getCommands()}"/>')
        x += pos.x_advance * scale
    return frame("".join(paths)), x - 10


def frame(body: str, style: str = "") -> str:
    return (
        '<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="80" viewBox="0 0 1400 80">'
        f'{style}<rect width="1400" height="80" fill="white"/>{body}</svg>'
    )


def text_svg(family: str, style: str = "") -> str:
    return frame(
        f'<text x="10" y="55" font-family="{family}" font-size="{SIZE}">{TEXT}</text>', style
    )


def embedded_svg() -> str:
    data = base64.b64encode(PINNED.read_bytes()).decode()
    style = (
        "<style>@font-face{font-family:Embedded;"
        f"src:url(data:font/ttf;base64,{data}) format('truetype');}}</style>"
    )
    return text_svg("Embedded", style)


def ink(png: bytes):
    import numpy as np
    from PIL import Image

    return np.asarray(Image.open(io.BytesIO(png)).convert("L")) < 128


def iou(a, b) -> float:
    return float((a & b).sum() / max((a | b).sum(), 1))


def width(mask) -> int:
    cols = mask.any(axis=0).nonzero()[0]
    return int(cols[-1] - cols[0]) if len(cols) else 0


def to_png(converter: str, svg: str) -> bytes:
    if converter == "vl-convert":
        import vl_convert

        return vl_convert.svg_to_png(svg)
    import cairosvg

    return cairosvg.svg2png(bytestring=svg.encode(), output_width=1400, output_height=80)


def font_test(converter: str, register: str):
    if converter == "vl-convert" and register == "yes":
        import vl_convert

        vl_convert.register_font_directory(str(FONTDIR))
    kerned, _ = outline_svg(True)
    plain, _ = outline_svg(False)
    ref_k, ref_p = ink(to_png(converter, kerned)), ink(to_png(converter, plain))
    rows = {}
    for label, svg in [
        (f"family {FAMILY}", text_svg(FAMILY)),
        ("family DejaVu Sans", text_svg("DejaVu Sans")),
        ("@font-face embedded", embedded_svg()),
    ]:
        try:
            m = ink(to_png(converter, svg))
            rows[label] = {
                "iou_kerned": round(iou(m, ref_k), 3),
                "iou_unkerned": round(iou(m, ref_p), 3),
                "ink_width": width(m),
            }
        except Exception as exc:  # noqa: BLE001
            rows[label] = {"error": repr(exc)[:120]}
    rows["reference ink width kerned / unkerned"] = [width(ref_k), width(ref_p)]
    print(json.dumps(rows))


# ---------------------------------------------------------------- sequence view pages

W, H = 595, 842
ROW = 100
ROWS_PER_PAGE = 7
BASE = 7.2
AA = ["Met", "Ser", "Lys", "Gly", "Glu", "Glu", "Leu", "Phe", "Thr", "Gly", "Val", "Val", "Pro"]


def row_svg(y: float, start: int, seq: str, rng: random.Random, family: str) -> str:
    x0 = 40
    out = []
    comp = seq.translate(str.maketrans("ACGT", "TGCA"))
    n = len(seq)
    for k in range(2):
        cut = x0 + rng.randrange(n) * BASE
        out.append(
            f'<g class="site" data-info="EcoRI"><polyline points="{cut},{y + 30} {cut},{y + 8 + 10 * k}" '
            f'stroke="#999" fill="none"/><text x="{cut + 2}" y="{y + 8 + 10 * k}" '
            f'font-family="{family}" font-size="8">EnzymeI{k}</text></g>'
        )
    ticks = "".join(f"M{x0 + i * BASE + BASE / 2},{y + 34}v{4 if i % 10 == 9 else 2}" for i in range(n))
    out.append(f'<path d="{ticks}" stroke="#666" stroke-width="0.5"/>')
    for i in range(9, n, 10):
        out.append(
            f'<text x="{x0 + i * BASE + BASE / 2}" y="{y + 32}" font-family="{family}" '
            f'font-size="6" text-anchor="middle">{start + i + 1}</text>'
        )
    out.append(
        f'<text x="{x0}" y="{y + 48}" font-family="{family}" font-size="10" '
        f'textLength="{n * BASE}" lengthAdjust="spacingAndGlyphs">{seq}</text>'
    )
    out.append(f'<line x1="{x0}" y1="{y + 51}" x2="{x0 + n * BASE}" y2="{y + 51}" stroke="#bbb"/>')
    out.append(
        f'<text x="{x0}" y="{y + 61}" font-family="{family}" font-size="10" '
        f'textLength="{n * BASE}" lengthAdjust="spacingAndGlyphs">{comp}</text>'
    )
    out.append(
        f'<text x="{x0 + n * BASE + 6}" y="{y + 48}" font-family="{family}" font-size="8">{start + n}</text>'
    )
    for c in range(n // 3):
        out.append(
            f'<text x="{x0 + c * 3 * BASE}" y="{y + 72}" font-family="{family}" font-size="8" '
            f'textLength="{3 * BASE}" lengthAdjust="spacingAndGlyphs">{rng.choice(AA)}</text>'
        )
    for k in range(2):
        a = rng.randrange(n // 2)
        b = min(n, a + rng.randrange(10, max(11, n // 2)))
        out.append(
            f'<g class="feature" data-info="CDS"><rect x="{x0 + a * BASE}" y="{y + 76 + 11 * k}" '
            f'width="{(b - a) * BASE}" height="9" fill="#77AADD" stroke="#333" stroke-width="0.5"/>'
            f'<text x="{x0 + a * BASE + 2}" y="{y + 83 + 11 * k}" font-family="{family}" '
            f'font-size="7">feature {k}</text></g>'
        )
    return "".join(out)


def page_svg(first_row: int, rows: int, bases: str, rng, family: str, height: float = H) -> str:
    body = []
    for r in range(rows):
        start = (first_row + r) * 60
        body.append(row_svg(30 + r * ROW, start, bases[start : start + 60], rng, family))
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{height}" '
        f'viewBox="0 0 {W} {height}"><rect width="{W}" height="{height}" fill="white"/>'
        + "".join(body)
        + "</svg>"
    )


def pages_test(converter: str, kb: int):
    import pypdf

    rng = random.Random(1)
    bases = "".join(rng.choice("ACGT") for _ in range(kb * 1000))
    total_rows = -(-len(bases) // 60)
    family = FAMILY if converter == "vl-convert" else "DejaVu Sans"
    if converter == "vl-convert":
        import vl_convert

        vl_convert.register_font_directory(str(FONTDIR))
        convert = vl_convert.svg_to_pdf
    else:
        import cairosvg

        def convert(svg):
            return cairosvg.svg2pdf(bytestring=svg.encode())

    t0 = time.perf_counter()
    svgs = [
        page_svg(p, min(ROWS_PER_PAGE, total_rows - p), bases, rng, family)
        for p in range(0, total_rows, ROWS_PER_PAGE)
    ]
    t_svg = time.perf_counter() - t0
    t0 = time.perf_counter()
    pdfs = [convert(s) for s in svgs]
    t_conv = time.perf_counter() - t0
    t0 = time.perf_counter()
    writer = pypdf.PdfWriter()
    for pdf in pdfs:
        writer.append(pypdf.PdfReader(io.BytesIO(pdf)))
    merged = io.BytesIO()
    writer.write(merged)
    t_merge = time.perf_counter() - t0
    (OUT / f"seqview-{converter}-{kb}kb.pdf").write_bytes(merged.getvalue())
    text = pypdf.PdfReader(io.BytesIO(pdfs[0])).pages[0].extract_text() or ""
    print(
        json.dumps(
            {
                "pages": len(svgs),
                "svg_MB": round(sum(map(len, svgs)) / 1e6, 2),
                "svg_s": round(t_svg, 2),
                "convert_s": round(t_conv, 2),
                "per_page_ms": round(1000 * t_conv / len(svgs)),
                "merge_s": round(t_merge, 2),
                "pdf_MB": round(len(merged.getvalue()) / 1e6, 2),
                "pdf_text_selectable": bases[:60] in text.replace("\n", "").replace(" ", ""),
            }
        )
    )


def png_test(converter: str, kb: int, scale: float):
    rng = random.Random(1)
    bases = "".join(rng.choice("ACGT") for _ in range(kb * 1000))
    rows = -(-len(bases) // 60)
    height = 60 + rows * ROW
    family = FAMILY if converter == "vl-convert" else "DejaVu Sans"
    svg = page_svg(0, rows, bases, rng, family, height)
    t0 = time.perf_counter()
    if converter == "vl-convert":
        import vl_convert

        vl_convert.register_font_directory(str(FONTDIR))
        png = vl_convert.svg_to_png(svg, scale=scale)
    else:
        import cairosvg

        png = cairosvg.svg2png(bytestring=svg.encode(), scale=scale)
    from PIL import Image

    Image.MAX_IMAGE_PIXELS = None
    size = Image.open(io.BytesIO(png)).size
    print(
        json.dumps(
            {"px": size, "s": round(time.perf_counter() - t0, 2), "png_MB": round(len(png) / 1e6, 2)}
        )
    )


# ---------------------------------------------------------------- driver


def child(*args: str, fontconfig: bool = False, timeout: int = 900) -> str:
    env = dict(os.environ)
    if fontconfig:
        env["FONTCONFIG_FILE"] = str(HERE / "fonts.conf")
    try:
        done = subprocess.run(
            [sys.executable, __file__, *args], env=env, capture_output=True, text=True, timeout=timeout
        )
    except subprocess.TimeoutExpired:
        return f"timeout after {timeout} s"
    return done.stdout.strip() or ("FAILED: " + done.stderr.strip().splitlines()[-1][:200])


def main():
    prep()
    print("== 1. Does the converter draw with the pinned font file?")
    print("   iou = ink overlap with the pinned file's own outlines; ~0.9+ means that font drew")
    for converter, register, fc in [
        ("vl-convert", "no", False),
        ("vl-convert", "yes", False),
        ("cairosvg", "no", False),
        ("cairosvg", "fontconfig-dir", True),
    ]:
        print(f"-- {converter}, font dir registered: {register}")
        print("  ", child("font", converter, register, fontconfig=fc))
    print("== 2. Sequence view as a paged PDF (7 rows of 60 bp a page), merged with pypdf")
    for kb in (10, 100):
        for converter in ("cairosvg", "vl-convert"):
            print(f"-- {converter}, {kb} kb")
            print("  ", child("pages", converter, str(kb)))
    print("== 3. Sequence view as one tall PNG")
    for kb, scale in ((5, 2), (10, 2), (20, 2), (50, 1)):
        for converter in ("cairosvg", "vl-convert"):
            print(f"-- {converter}, {kb} kb, scale {scale}")
            print("  ", child("png", converter, str(kb), str(scale)))


if __name__ == "__main__":
    if len(sys.argv) == 1:
        main()
    elif sys.argv[1] == "font":
        font_test(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "pages":
        pages_test(sys.argv[2], int(sys.argv[3]))
    elif sys.argv[1] == "png":
        png_test(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]))

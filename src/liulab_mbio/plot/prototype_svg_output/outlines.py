"""PROTOTYPE — throwaway. The sequence view page with every glyph drawn as an outline from the pinned font."""
import random, re, sys, time, io, json
import measure as m
from fontTools.ttLib import TTFont
from fontTools.pens.svgPathPen import SVGPathPen

font = TTFont(m.PINNED)
upem = font["head"].unitsPerEm
cmap = font.getBestCmap()
hmtx = font["hmtx"]
glyphs = font.getGlyphSet()


def defs(chars):
    out = []
    for ch in sorted(chars):
        name = cmap[ord(ch)]
        pen = SVGPathPen(glyphs)
        glyphs[name].draw(pen)
        out.append(f'<path id="g{ord(ch)}" d="{pen.getCommands()}"/>')
    return "<defs>" + "".join(out) + "</defs>"


TEXT = re.compile(r'<text x="([^"]+)" y="([^"]+)" font-family="[^"]+" font-size="([^"]+)"([^>]*)>([^<]*)</text>')


def to_outlines(svg):
    used = set()

    def repl(mt):
        x, y, size, rest, s = float(mt[1]), float(mt[2]), float(mt[3]), mt[4], mt[5]
        used.update(s)
        adv = [hmtx[cmap[ord(c)]][0] for c in s]
        natural = sum(adv) * size / upem
        tl = re.search(r'textLength="([^"]+)"', rest)
        stretch = float(tl[1]) / natural if tl else 1.0
        if 'text-anchor="middle"' in rest:
            x -= natural / 2
        k = size / upem
        uses, pos = [], 0
        for c, a in zip(s, adv):
            if c != " ":
                uses.append(f'<use href="#g{ord(c)}" x="{pos:.0f}"/>')
            pos += a * stretch
        return f'<g transform="matrix({k:.5f},0,0,{-k:.5f},{x:.2f},{y:.2f})">{"".join(uses)}</g>'

    body = TEXT.sub(repl, svg)
    return body.replace("<rect", defs(used) + "<rect", 1)


def main():
    import vl_convert, cairosvg
    rng = random.Random(1)
    bases = "".join(rng.choice("ACGT") for _ in range(420))
    page = m.page_svg(0, 7, bases, rng, "PinnedSans")
    outl = to_outlines(page)
    assert "<text" not in outl
    print("page bytes, text / outlines:", len(page), len(outl))
    for label, f in [
        ("vl-convert pdf, outlines", lambda: vl_convert.svg_to_pdf(outl)),
        ("vl-convert png x2, outlines", lambda: vl_convert.svg_to_png(outl, scale=2)),
        ("cairosvg pdf, outlines", lambda: cairosvg.svg2pdf(bytestring=outl.encode())),
        ("cairosvg png x2, outlines", lambda: cairosvg.svg2png(bytestring=outl.encode(), scale=2)),
    ]:
        f(); t0 = time.perf_counter()
        for _ in range(3): f()
        print(f"{label}: {1000 * (time.perf_counter() - t0) / 3:.0f} ms")
    vl_convert.register_font_directory(str(m.FONTDIR))
    open("out/page-text-vl.png", "wb").write(vl_convert.svg_to_png(page, scale=2))
    open("out/page-outlines-vl.png", "wb").write(vl_convert.svg_to_png(outl, scale=2))
    open("out/page-outlines-cairo.png", "wb").write(cairosvg.svg2png(bytestring=outl.encode(), scale=2))
    open("out/page-outlines.svg", "w").write(outl)


if __name__ == "__main__":
    main()

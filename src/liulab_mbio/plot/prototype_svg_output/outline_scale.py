"""PROTOTYPE — throwaway. The outline route at scale: a 100 kb sequence view as a paged PDF, and tall PNGs."""
import io, json, random, subprocess, sys, time
import measure as m


def pdf(kb):
    import pypdf, vl_convert
    from outlines import to_outlines
    rng = random.Random(1)
    bases = "".join(rng.choice("ACGT") for _ in range(kb * 1000))
    rows = -(-len(bases) // 60)
    t0 = time.perf_counter()
    svgs = [to_outlines(m.page_svg(p, min(m.ROWS_PER_PAGE, rows - p), bases, rng, "PinnedSans"))
            for p in range(0, rows, m.ROWS_PER_PAGE)]
    t_svg = time.perf_counter() - t0
    t0 = time.perf_counter()
    pdfs = [vl_convert.svg_to_pdf(s) for s in svgs]
    t_conv = time.perf_counter() - t0
    t0 = time.perf_counter()
    w = pypdf.PdfWriter()
    for b in pdfs:
        w.append(pypdf.PdfReader(io.BytesIO(b)))
    out = io.BytesIO(); w.write(out)
    t_merge = time.perf_counter() - t0
    (m.OUT / f"seqview-outlines-{kb}kb.pdf").write_bytes(out.getvalue())
    print(json.dumps({"pages": len(svgs), "svg_MB": round(sum(map(len, svgs)) / 1e6, 1),
                      "svg_s": round(t_svg, 1), "convert_s": round(t_conv, 1),
                      "per_page_ms": round(1000 * t_conv / len(svgs)), "merge_s": round(t_merge, 1),
                      "pdf_MB": round(len(out.getvalue()) / 1e6, 1)}))


def png(kb, scale):
    import vl_convert
    from outlines import to_outlines
    from PIL import Image
    rng = random.Random(1)
    bases = "".join(rng.choice("ACGT") for _ in range(kb * 1000))
    rows = -(-len(bases) // 60)
    svg = to_outlines(m.page_svg(0, rows, bases, rng, "PinnedSans", 60 + rows * m.ROW))
    t0 = time.perf_counter()
    data = vl_convert.svg_to_png(svg, scale=scale)
    Image.MAX_IMAGE_PIXELS = None
    print(json.dumps({"px": Image.open(io.BytesIO(data)).size, "s": round(time.perf_counter() - t0, 1),
                      "png_MB": round(len(data) / 1e6, 1)}))


def child(*args, timeout=600):
    try:
        d = subprocess.run([sys.executable, __file__, *args], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return f"timeout after {timeout} s"
    lines = (d.stderr or "").strip().splitlines()
    return d.stdout.strip() or f"FAILED (exit {d.returncode}): " + (lines[-1][:200] if lines else "no message")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        for kb in (10, 100):
            print("pdf", kb, "kb", child("pdf", str(kb)), flush=True)
        for kb, scale in ((5, 2), (10, 2), (20, 2), (50, 2), (100, 1)):
            print("png", kb, "kb scale", scale, child("png", str(kb), str(scale)), flush=True)
    elif sys.argv[1] == "pdf":
        pdf(int(sys.argv[2]))
    else:
        png(int(sys.argv[2]), float(sys.argv[3]))

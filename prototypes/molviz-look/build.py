"""PROTOTYPE, throwaway: one offline page comparing SeqViz, mbio today and SnapGene on pUC19-GFP."""

import base64
import html
import json
import sys
from pathlib import Path

build, out = Path(sys.argv[1]), Path(sys.argv[2])
seqviz = Path(sys.argv[3]).read_text()
mbio_plain = Path(sys.argv[4]).read_text()
props = json.loads((build / "props.json").read_text())
mbio_crowded = (build / "mbio-crowded.html").read_text()
snaps = {k: base64.b64encode((build / f).read_bytes()).decode() for k, f in (("snapgene-puc", "snapgene-map.jpg"), ("snapgene-ap1", "snapgene-ap1.jpg"))}

VARIANTS = [
    ("A", "SeqViz: pUC19-GFP", "seqviz-plain",
     "SeqViz's own look, drawn live from mbio's record: features, primers and the GFP translation. "
     "Look at the arrows, colours, label placement on the circle, and the sequence rows on the right. "
     "Primers show in the sequence but not on the circle. Rows you scroll past leave the page, so the browser's find misses them."),
    ("B", "SeqViz: crowded (AP-1 product, 28 enzymes)", "seqviz-crowded",
     "The AP-1 library product with every site of mbio's 28 enzymes, cut by SeqViz itself. "
     "Labels that don't fit fold into their neighbour, as 'SphI,+6' or a bare '+10'; hover a folded label to list what it holds."),
    ("C", "mbio today: pUC19-GFP", "mbio-plain",
     "mbio's current page for the same record: circle, line and sequence view, switches along the top."),
    ("D", "mbio today: crowded (AP-1 product, 28 enzymes)", "mbio-crowded",
     "The same AP-1 record and enzymes. mbio hides 6 enzyme-site labels it can't place and counts them in the notice; hover it to list them."),
    ("E", "SnapGene Viewer: pUC19-GFP", "snapgene-puc",
     "SnapGene Viewer's Map tab, with its own default enzymes (unique 6+ cutters). Every label sits outside the circle on a leader line."),
    ("F", "SnapGene Viewer: crowded (AP-1 product)", "snapgene-ap1",
     "SnapGene on the AP-1 product, its own enzyme set. Labels stack in columns on both sides with long leader lines; nothing is hidden or folded."),
]

def iframe(doc):
    return f'<iframe srcdoc="{html.escape(doc, quote=True)}"></iframe>'

panes = []
for key, title, kind, note in VARIANTS:
    if kind.startswith("seqviz"):
        body = f'<div class="sv" id="{kind}"></div>'
    elif kind == "mbio-plain":
        body = iframe(mbio_plain)
    elif kind == "mbio-crowded":
        body = iframe(mbio_crowded)
    else:
        body = f'<div class="img"><img alt="{html.escape(title)}" src="data:image/jpeg;base64,{snaps[kind]}"></div>'
    panes.append(f'<section class="pane" data-key="{key}" data-title="{html.escape(title)}" hidden>'
                 f'<p class="note"><b>{key}. {html.escape(title)}.</b> {html.escape(note)}</p>{body}</section>')

page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>PROTOTYPE: molviz look</title>
<style>
:root {{ color-scheme: light; --ink:#1d1d1f; --mute:#5b5b66; --line:#d9d9e0; --bar:#1d1d1f; }}
* {{ box-sizing: border-box; }}
html, body {{ margin:0; height:100%; background:#fff; color:var(--ink); font:15px/1.45 system-ui, sans-serif; }}
header {{ padding:12px 16px 8px; border-bottom:1px solid var(--line); }}
header h1 {{ font-size:16px; margin:0 0 4px; }}
header ol {{ margin:4px 0 0; padding-left:20px; color:var(--mute); }}
.pane {{ position:absolute; top:var(--head,120px); bottom:56px; left:0; right:0; display:flex; flex-direction:column; }}
.pane[hidden] {{ display:none; }}
.note {{ margin:0; padding:8px 16px; background:#f6f6f8; border-bottom:1px solid var(--line); font-size:14px; }}
.sv, iframe, .img {{ flex:1; min-height:0; width:100%; border:0; }}
.img {{ overflow:auto; text-align:center; }} .img img {{ max-width:100%; height:auto; }}
#bar {{ position:fixed; bottom:10px; left:50%; transform:translateX(-50%); display:flex; gap:8px; align-items:center;
  background:var(--bar); color:#fff; border-radius:999px; padding:6px 10px; box-shadow:0 4px 16px rgba(0,0,0,.25); font-size:14px; z-index:9; }}
#bar button {{ background:#3a3a40; color:#fff; border:0; border-radius:999px; width:32px; height:32px; font-size:16px; cursor:pointer; }}
#label {{ min-width:220px; text-align:center; }}
</style></head><body>
<header><h1>PROTOTYPE, throw away: what molviz copies from SeqViz's look</h1>
<ol><li>Which parts of SeqViz's look should molviz copy (arrows, colours, label placement, sequence rows, fonts)?</li>
<li>When labels crowd: group them as SeqViz does (B), or hide and count them as mbio does (D)?</li></ol></header>
{''.join(panes)}
<div id="bar"><button id="prev" aria-label="Previous">&larr;</button><span id="label"></span><button id="next" aria-label="Next">&rarr;</button></div>
<script>{seqviz}</script>
<script>
const PROPS = {json.dumps(props)};
const panes = [...document.querySelectorAll('.pane')];
const keys = panes.map(p => p.dataset.key);
const drawn = new Set();
document.documentElement.style.setProperty('--head', document.querySelector('header').offsetHeight + 'px');
function show(key) {{
  if (!keys.includes(key)) key = keys[0];
  for (const p of panes) p.hidden = p.dataset.key !== key;
  const pane = panes.find(p => p.dataset.key === key);
  document.getElementById('label').textContent = key + ' (' + pane.dataset.title + ')';
  const sv = pane.querySelector('.sv');
  if (sv && !drawn.has(sv.id)) {{
    drawn.add(sv.id);
    const p = PROPS[sv.id === 'seqviz-plain' ? 'plain' : 'crowded'];
    window.seqviz.Viewer(sv.id, {{...p, viewer: 'both', showComplement: true,
      style: {{ height: '100%', width: '100%' }} }}).render();
  }}
  const url = new URL(location.href); url.searchParams.set('variant', key); history.replaceState(null, '', url);
}}
function step(d) {{
  const i = keys.indexOf(new URL(location.href).searchParams.get('variant') || keys[0]);
  show(keys[(i + d + keys.length) % keys.length]);
}}
document.getElementById('prev').onclick = () => step(-1);
document.getElementById('next').onclick = () => step(1);
addEventListener('keydown', e => {{
  if (e.target.closest && e.target.closest('input,textarea,[contenteditable]')) return;
  if (e.key === 'ArrowLeft') step(-1); if (e.key === 'ArrowRight') step(1);
}});
show(new URL(location.href).searchParams.get('variant') || 'A');
</script></body></html>"""
out.write_text(page)
print(out, len(page.encode()))

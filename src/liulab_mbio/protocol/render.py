"""Render a protocol to one self-contained HTML page."""

import hashlib
import os
from collections.abc import Callable, Iterable
from html import escape
from importlib.resources import files
from pathlib import Path

from liulab_mbio.protocol.model import (
    Check,
    Gel,
    Material,
    Oligo,
    Protocol,
    ReactionTable,
    Reference,
    Step,
    ThermocyclerProgram,
    Timer,
)


def render_html(protocol: Protocol) -> str:
    """Return `protocol` as one HTML page with its styles and script inline.

    The page loads nothing over the network, remembers check marks and reaction counts in the
    browser's local storage when it can, and prints without its controls.
    """
    key = hashlib.sha256(repr(protocol).encode()).hexdigest()[:16]
    body = "".join(
        [
            _header(protocol),
            _materials(protocol.materials, protocol.equipment),
            _oligos(protocol.oligos),
            *(_step(n, step) for n, step in enumerate(protocol.steps, 1)),
            _references(protocol.references),
        ]
    )
    return (
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{escape(protocol.title)}</title>\n<style>\n{_asset('protocol.css')}</style>\n"
        f'</head>\n<body data-protocol="{key}">\n<main class="page">\n{body}</main>\n'
        f"<script>\n{_asset('protocol.js')}</script>\n</body>\n</html>\n"
    )


def write_html(protocol: Protocol, path: str | os.PathLike[str]) -> Path:
    """Write `render_html(protocol)` to `path` as UTF-8 and return the path."""
    out = Path(path)
    out.write_text(render_html(protocol), encoding="utf-8")
    return out


def _asset(name: str) -> str:
    return files("liulab_mbio.protocol").joinpath(name).read_text(encoding="utf-8")


def _num(value: float) -> str:
    return f"{value:.3f}".rstrip("0").rstrip(".")


def _duration(seconds: float) -> str:
    hours, rest = divmod(round(seconds), 3600)
    minutes, secs = divmod(rest, 60)
    parts = [f"{hours} h"] if hours else []
    if minutes:
        parts.append(f"{minutes} min")
    if secs or not parts:
        parts.append(f"{secs} s")
    return " ".join(parts)


def _clock(seconds: float) -> str:
    hours, rest = divmod(round(seconds), 3600)
    minutes, secs = divmod(rest, 60)
    return f"{hours}:{minutes:02}:{secs:02}" if hours else f"{minutes}:{secs:02}"


def _bullets(items: Iterable[str]) -> str:
    lines = "".join(f"<li>{escape(item)}</li>" for item in items)
    return f"<ul>{lines}</ul>" if lines else ""


def _copy(text: str, label: str = "Copy") -> str:
    return f'<button type="button" class="copy" data-copy="{escape(text)}">{escape(label)}</button>'


def _header(protocol: Protocol) -> str:
    parts = [f'<header class="intro">\n<h1>{escape(protocol.title)}</h1>\n']
    if protocol.summary:
        parts.append(f'<p class="summary">{escape(protocol.summary)}</p>\n')
    if protocol.overview:
        facts = "".join(
            f"<div><dt>{escape(k)}</dt><dd>{escape(v)}</dd></div>"
            for k, v in protocol.overview.items()
        )
        parts.append(f'<dl class="overview">{facts}</dl>\n')
    if protocol.highlights:
        lines = "".join(f"<p>{escape(one)}</p>" for one in protocol.highlights)
        parts.append(f'<div class="highlights">{lines}</div>\n')
    parts.append(_checks(protocol.checks))
    if protocol.steps:
        count = len(protocol.steps)
        parts.append(
            f'<div class="toolbar"><span class="progress" aria-live="polite">0 of {count} steps'
            ' done</span><button type="button" class="print">Print</button>'
            '<button type="button" class="clear">Clear checks</button></div>\n'
        )
        links = "".join(
            f'<li><a href="#step-{n}">{escape(step.title)}</a></li>'
            for n, step in enumerate(protocol.steps, 1)
        )
        parts.append(f'<nav class="toc" aria-label="Steps"><ol>{links}</ol></nav>\n')
    parts.append("</header>\n")
    return "".join(parts)


def _checks(checks: tuple[Check, ...]) -> str:
    """One badge per verdict, and the detail of every verdict that is not a pass."""
    if not checks:
        return ""
    badges = "".join(
        f'<li class="check is-{check.status}"><span class="check-name">{escape(check.name)}</span>'
        f'<span class="verdict">{escape(check.status)}</span></li>'
        for check in checks
    )
    details = "".join(
        f'<p class="check-detail"><strong>{escape(check.name)} {escape(check.status)}:</strong> '
        f"{escape(check.detail)}</p>"
        for check in checks
        if check.status != "pass" and check.detail
    )
    return f'<ul class="checks" aria-label="Checks">{badges}</ul>\n{details}\n'


def _cell(tag: str, css: str, inner: str) -> str:
    return f'<{tag} class="{css}">{inner}</{tag}>' if css else f"<{tag}>{inner}</{tag}>"


def _materials(materials: tuple[Material, ...], equipment: tuple[str, ...]) -> str:
    """Everything that is not an oligo, and the hardware as one light line under it."""
    if not materials and not equipment:
        return ""
    columns: list[tuple[str, Callable[[Material], str]]] = [
        ("Supplier", lambda m: m.supplier),
        ("Catalogue", lambda m: m.catalog),
        ("Storage", lambda m: m.storage),
        ("Per run", lambda m: m.amount),
        ("Note", lambda m: m.note),
    ]
    shown = [(label, get) for label, get in columns if any(get(m) for m in materials)]
    table = ""
    if materials:
        head = "<th>Name</th>" + "".join(f"<th>{escape(label)}</th>" for label, _ in shown)
        rows = "".join(
            f"<tr><td>{escape(material.name)}</td>"
            + "".join(f"<td>{escape(get(material))}</td>" for _, get in shown)
            + "</tr>"
            for material in materials
        )
        table = (
            f'<div class="scroll"><table><thead><tr>{head}</tr></thead>'
            f"<tbody>{rows}</tbody></table></div>"
        )
    line = ""
    if equipment:
        line = (
            f'<p class="equipment"><strong>Equipment:</strong> {escape(", ".join(equipment))}</p>'
        )
    return f'<section class="block materials">\n<h2>Materials</h2>\n{table}{line}\n</section>\n'


def _oligos(oligos: tuple[Oligo, ...]) -> str:
    """Render the order sheet: one row each, every sequence with a copy button."""
    if not oligos:
        return ""
    columns: list[tuple[str, str, Callable[[Oligo], str]]] = [
        # One decimal, so the column reads as one: `_num` prints 63 beside 63.1.
        ("Tm (°C)", "num", lambda o: "" if o.tm_c is None else f"{o.tm_c:.1f}"),
        ("For", "", lambda o: o.purpose),
        ("Working stock", "", lambda o: o.stock),
        ("Note", "", lambda o: o.note),
    ]
    shown = [(label, css, get) for label, css, get in columns if any(get(o) for o in oligos)]
    written = escape("Sequence (5'→3')")
    head = f'<th>Name</th><th>{written}</th><th class="num">Length</th>' + "".join(
        _cell("th", css, escape(label)) for label, css, _ in shown
    )
    rows = []
    for oligo in oligos:
        cells = [
            f"<td>{escape(oligo.name)}</td>",
            f'<td class="seq-cell"><code class="seq">{escape(oligo.sequence)}</code>'
            f"{_copy(oligo.sequence)}</td>",
            f'<td class="num">{len(oligo.sequence)}</td>',
            *(_cell("td", css, escape(get(oligo))) for _, css, get in shown),
        ]
        rows.append(f"<tr>{''.join(cells)}</tr>")
    copy_all = ""
    if len(oligos) > 1:
        sheet = "\n".join(f"{oligo.name}\t{oligo.sequence}" for oligo in oligos)
        copy_all = f"<p>{_copy(sheet, 'Copy all sequences')}</p>"
    return (
        '<section class="block oligos">\n<h2>Oligos</h2>\n'
        f'<div class="scroll"><table><thead><tr>{head}</tr></thead>'
        f"<tbody>{''.join(rows)}</tbody></table></div>{copy_all}\n</section>\n"
    )


def _step(n: int, step: Step) -> str:
    key = f"step-{n}"
    parts = [
        f'<section class="step" id="{key}">\n<h2 class="step-title"><label>'
        f'<input type="checkbox" class="done" data-key="{key}">'
        f'<span class="step-n">{n}</span><span>{escape(step.title)}</span></label></h2>\n'
    ]
    parts += [
        f'<p class="caution"><strong>Caution:</strong> {escape(c)}</p>\n' for c in step.cautions
    ]
    if step.instructions:
        items = "".join(
            f'<li><label><input type="checkbox" data-key="{key}-{i}">'
            f"<span>{escape(text)}</span></label></li>"
            for i, text in enumerate(step.instructions, 1)
        )
        parts.append(f'<ol class="instructions">{items}</ol>\n')
    parts += [_table(f"{key}-table-{i}", t) for i, t in enumerate(step.tables, 1)]
    parts += [_program(p) for p in step.programs]
    if step.timers:
        parts.append(f'<div class="timers">{"".join(_timer(t) for t in step.timers)}</div>\n')
    if step.expected or step.gels:
        gels = "".join(_gel(g) for g in step.gels)
        parts.append(
            f'<div class="expected"><h3>Expected result</h3>{_bullets(step.expected)}{gels}</div>\n'
        )
    if step.troubleshooting:
        entries = "".join(
            f"<dt>{escape(t.problem)}</dt><dd>{escape(t.solution)}</dd>"
            for t in step.troubleshooting
        )
        parts.append(f'<div class="trouble"><h3>Troubleshooting</h3><dl>{entries}</dl></div>\n')
    if step.notes:
        parts.append(f'<div class="notes"><h3>Notes</h3>{_bullets(step.notes)}</div>\n')
    parts.append("</section>\n")
    return "".join(parts)


def _table(key: str, table: ReactionTable) -> str:
    stock = any(c.stock for c in table.components)
    final = any(c.final for c in table.components)
    scale = table.reactions * (1 + table.overage)
    rows = []
    for component, mix in zip(table.components, table.mix_volumes(table.reactions), strict=True):
        cells = [f"<td>{escape(component.name)}</td>"]
        cells += [f"<td>{escape(component.stock)}</td>"] if stock else []
        cells += [f"<td>{escape(component.final)}</td>"] if final else []
        cells.append(f'<td class="num">{_num(component.volume_ul)}</td>')
        if mix is None:
            cells.append('<td class="num per-tube">each tube</td>')
        else:
            cells.append(f'<td class="num mix" data-ul="{component.volume_ul!r}">{_num(mix)}</td>')
        rows.append(f"<tr>{''.join(cells)}</tr>")
    blanks = "<td></td>" * (stock + final)
    in_mix = sum(c.volume_ul for c in table.components if c.master_mix)
    total = sum(c.volume_ul for c in table.components)
    head = (
        "<th>Component</th>"
        + ("<th>Stock</th>" if stock else "")
        + ("<th>Final</th>" if final else "")
        + '<th class="num">1 rxn (µL)</th>'
        + f'<th class="num">Mix for <span class="rxn-n">{table.reactions}</span> (µL)</th>'
    )
    foot = (
        f'<tr><th>Total</th>{blanks}<td class="num">{_num(total)}</td>'
        f'<td class="num mix" data-ul="{in_mix!r}">{_num(round(in_mix * scale, 2))}</td></tr>'
    )
    per_tube = [c for c in table.components if not c.master_mix]
    dispense = ""
    if in_mix:
        then = ", ".join(f"{_num(c.volume_ul)} µL {c.name}" for c in per_tube)
        dispense = f"Put {_num(in_mix)} µL of mix in each tube" + (
            f", then add {then}" if then else ""
        )
        dispense = f'<p class="dispense">{escape(dispense)}.</p>'
    caption = f"<figcaption>{escape(table.title)}</figcaption>" if table.title else ""
    return (
        f'<figure class="reaction" data-overage="{table.overage!r}">{caption}'
        f'<label class="count">Reactions <input type="number" class="rxn-count" name="reactions"'
        f' min="1" step="1" inputmode="numeric" value="{table.reactions}" data-key="{key}"></label>'
        f'<span class="muted">mix includes {_num(table.overage * 100)}% extra</span>'
        f'<div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody>'
        f"<tfoot>{foot}</tfoot></table></div>{dispense}</figure>\n"
    )


def _program(program: ThermocyclerProgram) -> str:
    meta = []
    if program.lid_temperature_c is not None:
        meta.append(f"lid {_num(program.lid_temperature_c)} °C")
    meta.append(f"{_duration(program.duration_seconds)} plus ramps")
    title = escape(program.title or "Thermocycler program")
    bodies = []
    for stage in program.stages:
        rows = []
        for i, step in enumerate(stage.incubations):
            time = "∞" if step.seconds is None else _duration(step.seconds)
            cycles = (
                f'<td class="num" rowspan="{len(stage.incubations)}">{stage.cycles}</td>'
                if i == 0
                else ""
            )
            rows.append(
                f"<tr><td>{escape(step.label)}</td>"
                f'<td class="num">{_num(step.temperature_c)} °C</td>'
                f'<td class="num">{time}</td>{cycles}</tr>'
            )
        bodies.append(f'<tbody class="stage">{"".join(rows)}</tbody>')
    return (
        f'<figure class="program"><figcaption>{title} '
        f'<span class="muted">· {escape(" · ".join(meta))}</span></figcaption>'
        '<div class="scroll"><table><thead><tr><th>Step</th><th class="num">Temperature</th>'
        f'<th class="num">Time</th><th class="num">Cycles</th></tr></thead>{"".join(bodies)}'
        "</table></div></figure>\n"
    )


def _timer(timer: Timer) -> str:
    return (
        f'<button type="button" class="timer" data-seconds="{_num(timer.seconds)}">'
        f'<span class="timer-label">{escape(timer.label)}</span>'
        f'<span class="timer-time">{_clock(timer.seconds)}</span>'
        '<span class="timer-action">Start</span></button>'
    )


# Gel drawing geometry, in SVG user units.
_LABEL_W, _LANE_W, _GAP = 46, 38, 10
_GEL_TOP, _RUN_TOP, _RUN_H = 22, 48, 210


def _gel(gel: Gel) -> str:
    lanes = [
        ("L", gel.ladder.name, gel.ladder.bands_bp),
        *((str(i), lane.label, lane.bands_bp) for i, lane in enumerate(gel.lanes, 1)),
    ]
    gel_w = len(lanes) * (_LANE_W + _GAP) + _GAP
    gel_h = _RUN_TOP - _GEL_TOP + _RUN_H + 16
    width, height = _LABEL_W + gel_w + 2, _GEL_TOP + gel_h + 2

    def y(bp: int) -> float:
        return _RUN_TOP + gel.migration(bp) * _RUN_H

    label = escape(gel.title or "Simulated agarose gel")
    svg = [
        f'<svg class="gel-image" viewBox="0 0 {width} {height}" role="img" aria-label="{label}">',
        f'<rect class="gel-bg" x="{_LABEL_W}" y="{_GEL_TOP}" width="{gel_w}" height="{gel_h}" rx="6"/>',
    ]
    last = float("-inf")
    for bp in sorted(set(gel.ladder.bands_bp), reverse=True):
        if y(bp) - last >= 10:
            svg.append(f'<text class="size" x="{_LABEL_W - 5}" y="{y(bp) + 3.5:.1f}">{bp}</text>')
            last = y(bp)
    for k, (mark, name, bands) in enumerate(lanes):
        x = _LABEL_W + _GAP + k * (_LANE_W + _GAP)
        kind = "lane ladder" if k == 0 else "lane"
        svg.append(f'<g class="{kind}" data-lane="{escape(name)}">')
        svg.append(f'<text class="lane-mark" x="{x + _LANE_W / 2}" y="14">{escape(mark)}</text>')
        svg.append(f'<rect class="well" x="{x}" y="{_GEL_TOP + 8}" width="{_LANE_W}" height="6"/>')
        svg += [
            f'<rect class="band" data-bp="{bp}" x="{x + 3}" y="{y(bp) - 2:.1f}" width="{_LANE_W - 6}"'
            f' height="4" rx="1.5"><title>{bp} bp</title></rect>'
            for bp in sorted(bands, reverse=True)
        ]
        svg.append("</g>")
    svg.append("</svg>")
    legend = "".join(
        f"<li><strong>{mark}</strong> {escape(name)}"
        + ("" if k == 0 else f": {', '.join(map(str, bands))} bp" if bands else ": no band")
        + "</li>"
        for k, (mark, name, bands) in enumerate(lanes)
    )
    caption = f"<figcaption>{escape(gel.title)}</figcaption>" if gel.title else ""
    return (
        f'<figure class="gel">{caption}{"".join(svg)}<ol class="gel-legend">{legend}</ol></figure>'
    )


def _references(references: tuple[Reference, ...]) -> str:
    if not references:
        return ""
    items = "".join(
        f"<li>{escape(r.text)}"
        + (f' <a href="{escape(r.url)}" rel="noreferrer">{escape(r.url)}</a>' if r.url else "")
        + "</li>"
        for r in references
    )
    return (
        f'<section class="block references">\n<h2>References</h2>\n<ol>{items}</ol>\n</section>\n'
    )

"""Render a protocol to one self-contained HTML page."""

import hashlib
import os
from collections.abc import Callable, Iterable
from html import escape
from importlib.resources import files
from pathlib import Path

from liulab_mbio.checks import Status
from liulab_mbio.plot.drawing import draw_plate
from liulab_mbio.protocol.model import (
    Bill,
    Check,
    Citation,
    Gel,
    Hole,
    Material,
    Oligo,
    Plate,
    Protocol,
    ReactionTable,
    Reference,
    Rule,
    Step,
    ThermocyclerProgram,
    Timer,
    Transfer,
    number,
)

#: What the page reads where a check carries no verdict, so it is never taken for a pass.
NO_VERDICT = "not judged"

#: What stands where a number would, so a hole can never be read as a figure.
NO_NUMBER = "no sourced number"

#: What each kind of hole says it is waiting on.
HOLE_KINDS = {
    "undecided": "the method has not decided",
    "unpublished": "nobody published it",
    "lab": "the lab's own stock",
    "unread": "a source was not read",
    "price": "no price record prices it",
}


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
            _plates(protocol),
            _bill(protocol.bill),
            *(_step(n, step, protocol) for n, step in enumerate(protocol.steps, 1)),
            _holes(protocol),
            _sources(protocol),
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
    parts.append(_hole_count(protocol.all_holes))
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
    """One badge per verdict, and the detail of every verdict that is not a pass.

    A check no sourced threshold judges carries no verdict, and its badge says so rather than
    reading as a pass.
    """
    if not checks:
        return ""
    badges = "".join(
        f'<li class="check is-{check.status or "none"}">'
        f'<span class="check-name">{escape(check.name)}</span>'
        f'<span class="verdict">{escape(_word(check.status))}</span></li>'
        for check in checks
    )
    details = "".join(
        f'<p class="check-detail"><strong>{escape(check.name)} '
        f"{escape(_word(check.status))}:</strong> {escape(check.detail)}</p>"
        for check in checks
        if check.status != "pass" and check.detail
    )
    return f'<ul class="checks" aria-label="Checks">{badges}</ul>\n{details}\n'


def _word(status: Status | None) -> str:
    """Return what a badge reads: the verdict, or that there is none."""
    return status if status is not None else NO_VERDICT


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
    cited = any(m.citation for m in materials)
    shown = [(label, get) for label, get in columns if any(get(m) for m in materials)]
    table = ""
    if materials:
        head = (
            "<th>Name</th>"
            + "".join(f"<th>{escape(label)}</th>" for label, _ in shown)
            + ("<th>Source</th>" if cited else "")
        )
        rows = "".join(
            f"<tr><td>{escape(material.name)}</td>"
            + "".join(f"<td>{escape(get(material))}</td>" for _, get in shown)
            + (f"<td>{_cite(material.citation)}</td>" if cited else "")
            + "</tr>"
            for material in materials
        )
        table = (
            f'<div class="scroll"><table><thead><tr>{head}</tr></thead>'
            f"<tbody>{rows}</tbody></table></div>"
        )
    carried = [(m, rule) for m in materials for rule in m.rules]
    line = ""
    if equipment:
        line = (
            f'<p class="equipment"><strong>Equipment:</strong> {escape(", ".join(equipment))}</p>'
        )
    return (
        '<section class="block materials">\n<h2>Materials</h2>\n'
        f"{table}{line}{_rules(carried)}\n</section>\n"
    )


def _oligos(oligos: tuple[Oligo, ...]) -> str:
    """Render the order sheet: one row each, every sequence with a copy button."""
    if not oligos:
        return ""
    columns: list[tuple[str, str, Callable[[Oligo], str]]] = [
        # One decimal, so the column reads as one: `number` prints 63 beside 63.1.
        ("Tm (°C)", "num", lambda o: "" if o.tm_c is None else f"{o.tm_c:.1f}"),
        ("For", "", lambda o: o.purpose),
        ("Working stock", "", lambda o: o.stock),
        ("Note", "", lambda o: o.note),
    ]
    shown = [(label, css, get) for label, css, get in columns if any(get(o) for o in oligos)]
    judged = any(oligo.status is not None for oligo in oligos)
    written = escape("Sequence (5'→3')")
    head = (
        f'<th>Name</th><th>{written}</th><th class="num">Length</th>'
        + "".join(_cell("th", css, escape(label)) for label, css, _ in shown)
        + ("<th>Checks</th>" if judged else "")
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
        if judged:
            cells.append(f'<td class="verdict-cell">{_verdict(oligo.status)}</td>')
        rows.append(f"<tr>{''.join(cells)}</tr>")
    copy_all = ""
    if len(oligos) > 1:
        sheet = "\n".join(f"{oligo.name}\t{oligo.sequence}" for oligo in oligos)
        copy_all = f"<p>{_copy(sheet, 'Copy all sequences')}</p>"
    return (
        '<section class="block oligos">\n<h2>Oligos</h2>\n'
        f'<div class="scroll"><table><thead><tr>{head}</tr></thead>'
        f"<tbody>{''.join(rows)}</tbody></table></div>"
        f"{_oligo_checks(oligos)}{copy_all}\n</section>\n"
    )


def _verdict(status: Status | None) -> str:
    """One row's verdict, as a word: the fill is never what carries it."""
    if status is None:
        return f'<span class="muted">{NO_VERDICT}</span>'
    return f'<span class="check is-{status}"><span class="verdict">{escape(status)}</span></span>'


def _oligo_checks(oligos: tuple[Oligo, ...]) -> str:
    """One closed toggle under the sheet: which rows are not a plain pass, and why.

    Closed and out of the table, so the order sheet stays a sheet. What to do about a warning
    is the reader's call, not this page's.
    """
    flagged = [oligo for oligo in oligos if oligo.checks]
    if not flagged:
        return ""
    items = "".join(
        f"<li><strong>{escape(oligo.name)}</strong> {escape(oligo.status or '')}: "
        + ", ".join(f"{escape(check.name)} {escape(check.detail)}" for check in oligo.checks)
        + "</li>"
        for oligo in flagged
    )
    return (
        f'<details class="oligo-checks"><summary>Checks on {len(flagged)} of {len(oligos)} '
        f"oligos</summary><ul>{items}</ul></details>"
    )


def _cite(citation: Citation | None) -> str:
    """One citation, as the page shows it beside the number it carries."""
    if citation is None:
        return ""
    where = f" {citation.locator}" if citation.locator else ""
    return (
        f'<a class="cite" href="#source-{escape(_slug(citation.source))}">'
        f"{escape(citation.source + where)}</a>"
    )


def _after(citation: Citation | None) -> str:
    """`_cite`, set off from the text before it; nothing when uncited."""
    return f" {_cite(citation)}" if citation else ""


def _slug(text: str) -> str:
    return "".join(char if char.isalnum() else "-" for char in text.casefold())


def _rules(rules: Iterable[tuple[Material, Rule]]) -> str:
    """Every rule the materials in this step carry, computed from the material, never stored.

    A rule hangs on the material, so it shows wherever the material is and no edit to a step's
    prose can drop it.
    """
    items = "".join(
        f'<li class="rule is-{rule.kind}"><strong>{escape(material.name)}: '
        f"{escape('never' if rule.kind == 'forbids' else 'always')} "
        f"{escape(rule.subject)}</strong> {escape(rule.detail)}{_after(rule.citation)}</li>"
        for material, rule in rules
    )
    return f'<ul class="rules" aria-label="Rules">{items}</ul>\n' if items else ""


def _count(n: int, noun: str) -> str:
    """`n` with `noun`, given an s where there is more than one."""
    return f"{n} {noun}" if n == 1 else f"{n} {noun}s"


def _hole_count(holes: tuple[Hole, ...]) -> str:
    """Return the banner: how many numbers have no source, and how many of those are prices.

    A price is missing because no price record prices it, so it is counted apart from the
    numbers the bench needs. What a hole means is said once, in the section this links to.
    """
    if not holes:
        return ""
    prices = sum(1 for hole in holes if hole.kind == "price")
    split = (
        f": {_count(len(holes) - prices, 'bench number')} and {_count(prices, 'price')}"
        if prices
        else ""
    )
    one = len(holes) == 1
    return (
        f'<p class="hole-count"><a href="#holes"><strong>{len(holes)}</strong> '
        f"{'number' if one else 'numbers'} in this protocol {'has' if one else 'have'} "
        f"no source</a>{split}.</p>\n"
    )


def _hole(hole: Hole) -> str:
    """One hole, which reads as a hole and never as a value."""
    where = f"{escape(hole.where)}: " if hole.where else ""
    filled = f" <em>Filled by {escape(hole.filled_by)}.</em>" if hole.filled_by else ""
    issue = f' <span class="hole-issue">{escape(hole.issue)}</span>' if hole.issue else ""
    return (
        f'<li class="hole" id="hole-{escape(hole.id)}"><span class="hole-id">{escape(hole.id)}'
        f'</span> <span class="hole-none">{NO_NUMBER}</span> — {where}{escape(hole.missing)} '
        f'<span class="hole-kind">{escape(HOLE_KINDS[hole.kind])}</span>{filled}{issue}</li>'
    )


def _holes(protocol: Protocol) -> str:
    """Every hole the protocol carries, collected under its stable id."""
    holes = protocol.all_holes
    if not holes:
        return ""
    items = "".join(_hole(hole) for hole in holes)
    return (
        '<section class="block holes" id="holes">\n<h2>Holes</h2>\n'
        f"<p>{_count(len(holes), 'number')} this protocol would otherwise have to invent. A "
        "hole is a defect in what the package knows, not a failure of the run, and it is "
        f"never filled with a guess.</p>\n<ul>{items}</ul>\n</section>\n"
    )


def _plates(protocol: Protocol) -> str:
    """Each plate drawn as its wells, and each vessel named beside them."""
    if not protocol.plates and not protocol.vessels:
        return ""
    figures = "".join(_plate(one) for one in protocol.plates)
    vessels = ""
    if protocol.vessels:
        items = "".join(
            f"<li><strong>{escape(v.name)}</strong>"
            + "".join(f" · {escape(text)}" for text in (v.kind, v.catalog, v.holds, v.note) if text)
            + "</li>"
            for v in protocol.vessels
        )
        vessels = f'<ul class="vessels">{items}</ul>'
    return f'<section class="block plates">\n<h2>Plates</h2>\n{figures}{vessels}\n</section>\n'


def _plate(one: Plate) -> str:
    # Labels are drawn beside the seating and resolve against nothing, so a well a step names
    # keeps its own entry where a plate carries both.
    drawn = draw_plate(
        one.name,
        one.rows,
        one.columns,
        one.row_labels,
        seating={**one.labels, **one.seating},
    )
    legend = "".join(
        f'<li><span class="swatch" style="background:{fill}"></span>{escape(kind)}</li>'
        for kind, fill in drawn.layout.legend
    )
    facts = " · ".join(
        text for text in (f"{one.wells} wells", one.catalog, one.holds, one.note) if text
    )
    return (
        f'<figure class="plate" data-plate="{escape(one.name)}">{drawn.element()}'
        f'<figcaption>{escape(one.name)} <span class="muted">{escape(facts)}</span></figcaption>'
        + (f'<ul class="plate-legend">{legend}</ul>' if legend else "")
        + "</figure>\n"
    )


def _transfer(transfer: Transfer) -> str:
    """Return a transfer as a table: where each thing goes, so no step describes it."""
    rows = "".join(
        f"<tr><td>{escape(move.source.plate)} {escape(move.source.well)}</td>"
        f"<td>{escape(move.destination.plate)} {escape(move.destination.well)}</td>"
        f'<td class="num">{number(move.volume_ul)}</td></tr>'
        for move in transfer.moves
    )
    meta = " · ".join(
        text
        for text in (
            transfer.instrument,
            f"{len(transfer.moves)} wells",
            " → ".join(transfer.plates),
            transfer.note,
        )
        if text
    )
    return (
        f'<figure class="transfer"><figcaption>{escape(transfer.title)} '
        f'<span class="muted">{escape(meta)}</span>{_after(transfer.citation)}</figcaption>'
        '<div class="scroll"><table><thead><tr><th>From</th><th>To</th>'
        f'<th class="num">µL</th></tr></thead><tbody>{rows}</tbody></table></div></figure>\n'
    )


def _bill(bill: Bill | None) -> str:
    """Return the bill: what the run consumes, and a hole wherever no row priced it."""
    if bill is None:
        return ""
    money = f"Charge ({escape(bill.currency)})" if bill.currency else "Charge"
    headroom = any(row.headroom for row in bill.rows)
    rows = []
    for row in bill.rows:
        charge = (
            f'<span class="hole-none">{NO_NUMBER}</span>'
            if row.hole
            else escape(row.charge) + _after(row.citation)
        )
        cells = [
            f"<td>{escape(row.item)}</td>",
            f"<td>{escape(row.key)}</td>",
            f'<td class="num">{number(row.quantity)} {escape(row.unit)}</td>',
            f'<td class="num">{charge}</td>',
        ]
        if headroom:
            cells.append(f"<td>{escape(row.headroom)}</td>")
        rows.append(f"<tr{' class="is-holed"' if row.hole else ''}>{''.join(cells)}</tr>")
    head = f'<th>Item</th><th>Key</th><th class="num">Quantity</th><th class="num">{money}</th>' + (
        "<th>Headroom</th>" if headroom else ""
    )
    total = (
        f'<tfoot><tr><th>Total</th><td></td><td></td><td class="num">{escape(bill.total)}</td>'
        + ("<td></td>" if headroom else "")
        + "</tr></tfoot>"
        if bill.total
        else ""
    )
    return (
        f'<section class="block bill">\n<h2>{escape(bill.title or "Bill")}</h2>\n'
        "<p>Quantities come from the design and are here whatever is loaded. A charge comes "
        "only from a price record; where none prices a row, the money is a hole and no figure "
        "is estimated. Price steers no part of the design.</p>\n"
        f'<div class="scroll"><table><thead><tr>{head}</tr></thead>'
        f"<tbody>{''.join(rows)}</tbody>{total}</table></div>\n</section>\n"
    )


def _sources(protocol: Protocol) -> str:
    """Every document a number was read from, so a citation resolves on the page itself."""
    if not protocol.sources:
        return ""
    items = "".join(
        f'<li id="source-{escape(_slug(key))}"><strong>{escape(key)}</strong> '
        f"{escape(source.document)}"
        + "".join(
            f" · {escape(text)}" for text in (source.edition, source.read_as, source.date) if text
        )
        + (
            f' <a href="{escape(source.url)}" rel="noreferrer">{escape(source.url)}</a>'
            if source.url
            else ""
        )
        + "</li>"
        for key, source in protocol.sources.items()
    )
    return f'<section class="block sources">\n<h2>Sources</h2>\n<ul>{items}</ul>\n</section>\n'


def _step(n: int, step: Step, protocol: Protocol) -> str:
    key = f"step-{n}"
    parts = [
        f'<section class="step" id="{key}">\n<h2 class="step-title"><label>'
        f'<input type="checkbox" class="done" data-key="{key}">'
        f'<span class="step-n">{n}</span><span>{escape(step.title)}</span></label></h2>\n'
    ]
    parts.append(_rules(protocol.rules_for(step)))
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
    parts += [_transfer(t) for t in step.transfers]
    if step.holes:
        items = "".join(_hole(hole) for hole in step.holes)
        parts.append(f'<ul class="holes-here" aria-label="Holes">{items}</ul>\n')
    if step.timers:
        parts.append(f'<div class="timers">{"".join(_timer(t) for t in step.timers)}</div>\n')
    if step.expected or step.gels:
        gels = "".join(_gel(g) for g in step.gels)
        parts.append(
            f'<div class="expected"><h3>Expected result</h3>{_bullets(step.expected)}{gels}</div>\n'
        )
    if step.troubleshooting:
        entries = "".join(
            f"<dt>{escape(t.problem)}</dt><dd>{escape(t.solution)}{_after(t.citation)}</dd>"
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
    for component in table.components:
        cells = [f"<td>{escape(component.name)}{_after(component.citation)}</td>"]
        cells += [f"<td>{escape(component.stock)}</td>"] if stock else []
        cells += [f"<td>{escape(component.final)}</td>"] if final else []
        cells.append(f'<td class="num">{number(component.volume_ul)}</td>')
        if component.master_mix:
            # `protocol.js` writes this cell again from `data-ul`, by the same arithmetic.
            mix = number(component.volume_ul * scale)
            cells.append(f'<td class="num mix" data-ul="{component.volume_ul!r}">{mix}</td>')
        else:
            cells.append('<td class="num per-tube">each tube</td>')
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
        f'<tr><th>Total</th>{blanks}<td class="num">{number(total)}</td>'
        f'<td class="num mix" data-ul="{in_mix!r}">{number(in_mix * scale)}</td></tr>'
    )
    per_tube = [c for c in table.components if not c.master_mix]
    dispense = ""
    if in_mix:
        then = ", ".join(f"{number(c.volume_ul)} µL {c.name}" for c in per_tube)
        dispense = f"Put {number(in_mix)} µL of mix in each tube" + (
            f", then add {then}" if then else ""
        )
        dispense = f'<p class="dispense">{escape(dispense)}.</p>'
    caption = f"<figcaption>{escape(table.title)}</figcaption>" if table.title else ""
    return (
        f'<figure class="reaction" data-overage="{table.overage!r}">{caption}'
        f'<label class="count">Reactions <input type="number" class="rxn-count" name="reactions"'
        f' min="1" step="1" inputmode="numeric" value="{table.reactions}" data-key="{key}"></label>'
        f'<span class="muted">mix includes {number(table.overage * 100)}% extra</span>'
        f'<div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody>'
        f"<tfoot>{foot}</tfoot></table></div>{dispense}</figure>\n"
    )


def _program(program: ThermocyclerProgram) -> str:
    meta = []
    if program.lid_temperature_c is not None:
        meta.append(f"lid {number(program.lid_temperature_c)} °C")
    if program.duration_seconds is not None:
        meta.append(f"{_duration(program.duration_seconds)} plus ramps")
    title = escape(program.title or "Thermocycler program")
    caption = f' <span class="muted">· {escape(" · ".join(meta))}</span>' if meta else ""
    bodies = []
    for stage in program.stages:
        count = (
            f'<span class="hole-none">{NO_NUMBER}</span>'
            if stage.cycles is None
            else str(stage.cycles)
        ) + _after(stage.citation)
        rows = []
        for i, step in enumerate(stage.incubations):
            time = "∞" if step.seconds is None else _duration(step.seconds)
            cycles = (
                f'<td class="num" rowspan="{len(stage.incubations)}">{count}</td>' if i == 0 else ""
            )
            rows.append(
                f"<tr><td>{escape(step.label)}{_after(step.citation)}</td>"
                f'<td class="num">{number(step.temperature_c)} °C</td>'
                f'<td class="num">{time}</td>{cycles}</tr>'
            )
        bodies.append(f'<tbody class="stage">{"".join(rows)}</tbody>')
    return (
        f'<figure class="program"><figcaption>{title}{caption}</figcaption>'
        '<div class="scroll"><table><thead><tr><th>Step</th><th class="num">Temperature</th>'
        f'<th class="num">Time</th><th class="num">Cycles</th></tr></thead>{"".join(bodies)}'
        "</table></div></figure>\n"
    )


def _timer(timer: Timer) -> str:
    return (
        f'<button type="button" class="timer" data-seconds="{timer.seconds!r}">'
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

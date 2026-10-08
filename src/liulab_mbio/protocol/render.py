"""Render a protocol to one self-contained HTML page, and a project to a folder of them.

A protocol on its own is one page that loads nothing. A run of several is a folder: one page
each, an index, and two pages the run shares. Every link between them is relative and every
page is still self-contained, so the folder opens from disk and survives being zipped.
"""

import hashlib
import os
import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from html import escape
from importlib.resources import files
from pathlib import Path

from liulab_mbio.checks import Status
from liulab_mbio.plot.drawing import draw_map, draw_plate
from liulab_mbio.plot.fonts import BOLD, MONO, SANS
from liulab_mbio.plot.page import font_face
from liulab_mbio.protocol.model import (
    Bill,
    Check,
    Citation,
    Figure,
    Gel,
    Hole,
    Item,
    Material,
    Oligo,
    Plate,
    Project,
    Protocol,
    ReactionTable,
    Reference,
    Rule,
    Source,
    Step,
    ThermocyclerProgram,
    Timer,
    Transfer,
    Wait,
    number,
    read_project,
    write_project,
)

#: What the page reads where a check carries no verdict, so it is never taken for a pass.
NO_VERDICT = "not judged"

#: What stands where a number would, so a hole can never be read as a figure.
NO_NUMBER = "no sourced number"

#: What a project folder calls the data every one of its pages is rendered from.
PROJECT_DATA_FILE = "project.json"

#: What a project folder calls the three pages that are not one protocol's.
INDEX_FILE = "index.html"
REAGENTS_FILE = "reagents.html"
REFERENCES_FILE = "references.html"

#: The longest a protocol's title may run in the file its page is written to.
NAME_CHARS = 48

#: What each kind of hole says it is waiting on.
HOLE_KINDS = {
    "undecided": "the method has not decided",
    "unpublished": "nobody published it",
    "lab": "the lab's own stock",
    "unread": "a source was not read",
    "price": "no price record prices it",
}


@dataclass(frozen=True, slots=True)
class Page:
    """One protocol's page, as every other page of its folder sees it.

    Parameters
    ----------
    title
        What the protocol is called, which is how every other page names it.
    href
        The file it was written to, relative to the folder.
    steps
        How many numbered steps it holds.
    key
        What the page remembers its check marks under, so another page of the folder can read
        how far the bench got.
    """

    title: str
    href: str
    steps: int = 0
    key: str = ""


@dataclass(frozen=True, slots=True)
class Folder:
    """A project folder as one of its pages sees the rest: what the nav bar links to, in order.

    Parameters
    ----------
    pages
        One per protocol, in the order they are run.
    index
        The way in to the folder.
    reagents, references
        The two pages the whole run shares.
    """

    pages: tuple[Page, ...]
    index: str = INDEX_FILE
    reagents: str = REAGENTS_FILE
    references: str = REFERENCES_FILE


@dataclass(frozen=True, slots=True)
class ProjectFiles:
    """The folder one project wrote: the data it was rendered from, and the pages.

    Parameters
    ----------
    data
        The whole run as JSON, which ``protocol render`` turns back into these pages.
    index
        The way in.
    protocols
        One page per protocol, in the order they are run.
    reagents, references
        The two pages the whole run shares.
    """

    data: Path
    index: Path
    protocols: tuple[Path, ...]
    reagents: Path
    references: Path

    @property
    def paths(self) -> tuple[Path, ...]:
        """Every file written, the first written first."""
        return (self.data, *self.protocols, self.index, self.reagents, self.references)


def page_key(value: Protocol | Project) -> str:
    """Return what a page remembers its check marks under: a digest of its own content.

    Two pages of one folder never share it, and every `file://` page in a browser shares one
    store, so the key is what keeps one page's marks off another.

    Examples
    --------
    >>> page_key(Protocol("Demo")) == page_key(Protocol("Demo"))
    True
    """
    return hashlib.sha256(repr(value).encode()).hexdigest()[:16]


def page_name(place: int, title: str) -> str:
    """Return the file one protocol's page is written to: its place in the chain, then its title.

    Examples
    --------
    >>> page_name(2, "LR reaction")
    '02-lr-reaction.html'
    """
    slug = re.sub("-+", "-", _slug(title)).strip("-")[:NAME_CHARS].strip("-")
    return f"{place:02d}-{slug}.html" if slug else f"{place:02d}.html"


def render_html(
    protocol: Protocol,
    *,
    folder: Folder | None = None,
    here: str = "",
    base: str | os.PathLike[str] | None = None,
) -> str:
    """Return `protocol` as one HTML page with its styles and script inline.

    The page loads nothing over the network, remembers check marks and reaction counts in the
    browser's local storage when it can, and prints without its controls. Given the `folder` it
    stands in and its own address in it, the page gains the run's nav bar, a column either side
    and a line naming its neighbours, and the line prints with it.

    A `Figure` names its record by a path relative to `base`, the directory the protocol was read
    from; the current directory where nobody says.

    Raises
    ------
    FileNotFoundError
        If a figure's record is not there, naming the step and the path.
    """
    place = _place(folder, here) if folder else ""
    beside = Path() if base is None else Path(base)
    # A section labels the steps under it, so it is written once, where it changes.
    sections = [""] + [step.section for step in protocol.steps]
    body = "".join(
        [
            _header(protocol, toc=folder is None, place=place),
            _materials(protocol.materials, protocol.equipment),
            _oligos(protocol.oligos),
            _plates(protocol),
            _bill(protocol.bill),
            *(
                _step(
                    n,
                    step,
                    protocol,
                    beside,
                    step.section if step.section != sections[n - 1] else "",
                )
                for n, step in enumerate(protocol.steps, 1)
            ),
            _holes(protocol),
            _sources(protocol.sources),
            _references(protocol.references),
        ]
    )
    return _page(protocol.title, page_key(protocol), body, folder, here, _within(protocol))


def write_html(
    protocol: Protocol,
    path: str | os.PathLike[str],
    *,
    folder: Folder | None = None,
    base: str | os.PathLike[str] | None = None,
) -> Path:
    """Write `render_html(protocol)` to `path` as UTF-8 and return the path.

    Inside a `folder`, the file's own name is the address the other pages link it by. A figure's
    record is read from `base`, and from beside the page where nobody says: a pipeline writes the
    records, the data and the page into one directory.
    """
    out = Path(path)
    page = render_html(
        protocol, folder=folder, here=out.name, base=out.parent if base is None else base
    )
    out.write_text(page, encoding="utf-8")
    return out


def render_index(project: Project, folder: Folder) -> str:
    """Return the way in to a run: what it achieves, how it fits together and how long it takes.

    The explanation of why the run is shaped as it is stands here and on no protocol page, so a
    step never stops to explain a decision.
    """
    holes = _run_holes(project, folder)
    summary = f'<p class="summary">{escape(project.summary)}</p>\n' if project.summary else ""
    jumps = [(f"topic-{_slug(topic.title)}", topic.title) for topic in project.background]
    parts = [
        f'<header class="intro">\n<h1>{escape(project.title)}</h1>\n{summary}'
        f"{_checks(project.checks)}{_hole_count(holes, 'this run')}</header>\n",
        _background(project),
    ]
    for anchor, label, block in (
        ("flow", "How the run fits together", _flow(project, folder)),
        ("protocols", "Protocols", _protocols(folder)),
        ("schedule", "Schedule", _schedule(project, folder)),
        ("holes", "Holes", _run_hole_list(project, folder, holes)),
    ):
        if block:
            jumps.append((anchor, label))
            parts.append(block)
    return _page(
        project.title, page_key(project), "".join(parts), folder, folder.index, _jumps(jumps)
    )


def render_reagents(project: Project, folder: Folder) -> str:
    """Return the shopping list: everything the run buys and every instrument it uses, once.

    Merged across the protocols, each row naming which of them take it. A protocol page keeps
    its own list, which is what the bench reads while it works; this is what is ordered.
    """
    blocks = (
        ("materials", "Materials", _merged_materials(project)),
        ("kit", "Equipment and plasticware", _kit(project)),
        ("bill", "Bill", _bill(project.bill)),
    )
    jumps = [(anchor, label) for anchor, label, block in blocks if block]
    lead = (
        "Everything the whole run takes, each row naming the protocols that take it. Each "
        "protocol page carries its own list, which is what to lay out before that sitting."
        if jumps
        else "No protocol of this run lists a reagent or an instrument."
    )
    body = f"<h1>Reagents and equipment</h1>\n<p>{lead}</p>\n" + "".join(
        block for _, _, block in blocks
    )
    return _shared(project, folder, folder.reagents, "Reagents and equipment", body, jumps)


def render_references(project: Project, folder: Folder) -> str:
    """Return every document the run was built from, each naming the protocols that cite it."""
    blocks = (
        ("references", "References", _references(*_merged_references(project))),
        ("sources", "Sources", _sources(*_merged_sources(project))),
    )
    jumps = [(anchor, label) for anchor, label, block in blocks if block]
    empty = "" if jumps else "<p>No protocol of this run cites a document.</p>\n"
    return _shared(
        project,
        folder,
        folder.references,
        "References",
        f"<h1>References</h1>\n{empty}" + "".join(block for _, _, block in blocks),
        jumps,
    )


def _shared(
    project: Project,
    folder: Folder,
    here: str,
    heading: str,
    body: str = "",
    jumps: Iterable[tuple[str, str]] = (),
) -> str:
    """One of the two pages the whole run shares, which carry no protocol of their own."""
    return _page(
        f"{heading} — {project.title}",
        page_key(project) + _slug(here),
        body or f"<h1>{heading}</h1>\n",
        folder,
        here,
        _jumps(jumps),
    )


def _jumps(jumps: Iterable[tuple[str, str]]) -> str:
    """Return the right column of a page that holds no steps: where to jump inside it."""
    links = "".join(f'<li><a href="#{anchor}">{escape(label)}</a></li>' for anchor, label in jumps)
    if not links:
        return ""
    return (
        '<nav class="column within" aria-label="This page">'
        f'<h2 class="column-title">This page</h2><ol>{links}</ol></nav>\n'
    )


def _background(project: Project) -> str:
    """Return what the reader is told before the first protocol, one topic to a block."""
    return "".join(
        f'<section class="block topic" id="topic-{escape(_slug(topic.title))}">\n'
        f"<h2>{escape(topic.title)}</h2>\n"
        + "".join(f"<p>{escape(line)}</p>" for line in topic.body)
        + "\n</section>\n"
        for topic in project.background
    )


def _protocols(folder: Folder) -> str:
    """Return the run in order, each page carrying the key its own check marks are kept under.

    `protocol.js` reads that key on this page and writes how far the bench got, so the count
    comes from the pages themselves and is never stored twice.
    """
    if not folder.pages:
        return ""
    rows = "".join(
        f'<li data-page-key="{escape(page.key)}" data-steps="{page.steps}">'
        f'<a href="{escape(page.href)}">{escape(page.title)}</a> '
        f'<span class="page-progress muted">{_count(page.steps, "step")}</span></li>'
        for page in folder.pages
    )
    return (
        '<section class="block chain-list" id="protocols">\n<h2>Protocols</h2>\n'
        f'<ol class="chain-pages">{rows}</ol>\n</section>\n'
    )


def _flow(project: Project, folder: Folder) -> str:
    """Return the chain as boxes and the names handed between them, in HTML and never a drawing.

    Generated from what each protocol declares it consumes and produces: a box is a protocol,
    and above it stand the names it is handed, each saying where it came from. Drawn in text, so
    a browser's find reaches it, it is selectable and it prints with no font embedded.
    """
    if not folder.pages:
        return ""
    came_from = {item.name: "the bench already holds it" for item in project.inputs}
    left: dict[str, Item] = {}
    boxes = []
    for page, protocol in zip(folder.pages, project.protocols, strict=True):
        needs = "".join(
            _flow_item(item, came_from.get(item.name, "")) for item in protocol.consumes
        )
        for item in protocol.consumes:
            left.pop(item.name, None)
        hand = f'<ul class="flow-hand">{needs}</ul>' if needs else ""
        boxes.append(
            # The step count stands in the Protocols list, where `protocol.js` keeps it up to
            # date; printed here too it would be the same number twice, one of them stale.
            f'<li>{hand}<a class="flow-box" href="{escape(page.href)}">'
            f'<span class="flow-title">{escape(page.title)}</span></a></li>'
        )
        for item in protocol.produces:
            came_from[item.name] = f"from {page.title}"
            left[item.name] = item
    ends = "".join(_flow_item(item, "the run ends holding it") for item in left.values())
    return (
        '<section class="block flow" id="flow">\n<h2>How the run fits together</h2>\n'
        f'<ol class="flow-chain">{"".join(boxes)}</ol>'
        + (f'<ul class="flow-hand flow-end">{ends}</ul>' if ends else "")
        + "\n</section>\n"
    )


def _flow_item(item: Item, came_from: str) -> str:
    """One name handed between two protocols, and where it came from, or that nothing hands it."""
    said = came_from or "nothing in the run hands this over"
    dangling = "" if came_from else " is-dangling"
    wanted = escape(", ".join(item.spec))
    spec = f' <span class="muted">· {wanted}</span>' if item.spec else ""
    return (
        f'<li class="flow-hand-item{dangling}">'
        f'<span class="flow-item">{escape(item.name)}</span> '
        f'<span class="muted">{escape(item.what)}</span>{spec} '
        f'<span class="flow-from">{escape(said)}</span></li>'
    )


#: The schedule's columns after the protocol's own name. The waiting goes in a row of its own,
#: because how long a vendor takes is whatever they say and no column is wide enough for it.
SCHEDULE_COLUMNS = ("Steps", "Held", "Hands-on", "Unattended", "Holding nothing")


def _schedule(project: Project, folder: Folder) -> str:
    """Return the run's time plan as a table, with a hole wherever nothing states a number.

    A table and not a bar chart: most of a run's calendar is time nobody attends and nobody can
    date, and a bar draws an unknown wait as a length, which is a claim. Held is summed from the
    timers and thermocycler programs the steps already hold; hands-on is stated or it is a hole,
    and the unattended share is the difference wherever a step states both.

    A column no protocol states a number in is a hole in every row, which says the same thing
    once per row that one sentence under the table says once.
    """
    if not folder.pages:
        return ""
    times = [_time_of(protocol) for protocol in project.protocols]
    totals = _Time()
    for time in times:
        totals = totals.and_(time)
    kept = tuple(i for i, column in enumerate(totals.columns()) if column is not None)
    rows = "".join(
        f'<tr><td><a href="{escape(page.href)}">{escape(page.title)}</a></td>'
        f'{_cells(time, kept)}</tr><tr class="wait-row"><td colspan="{len(kept) + 1}">'
        f"{_waiting(time.waits)}</td></tr>"
        for page, time in zip(folder.pages, times, strict=True)
    )
    head = "<th>Protocol</th>" + "".join(
        f'<th class="num">{escape(SCHEDULE_COLUMNS[i])}</th>' for i in kept
    )
    missing = [SCHEDULE_COLUMNS[i].lower() for i in range(len(SCHEDULE_COLUMNS)) if i not in kept]
    left_out = (
        f" Nothing in this run states {_and(missing)}, so "
        f"{'that column is' if len(missing) == 1 else 'those columns are'} left out."
        if missing
        else ""
    )
    return (
        '<section class="block schedule" id="schedule">\n<h2>Schedule</h2>\n'
        '<div class="scroll"><table class="schedule"><thead><tr>'
        f"{head}</tr></thead><tbody>{rows}</tbody>"
        f"<tfoot><tr><th>Total</th>{_cells(totals, kept)}</tr></tfoot></table></div>\n"
        "<p>A blank here is a number nobody has stated, never a zero. Most of this run is time "
        "nobody attends and nobody can date, so it is written down rather than drawn as a "
        f"length.{left_out}</p>\n</section>\n"
    )


def _cells(time: "_Time", kept: tuple[int, ...]) -> str:
    """One row's kept columns, each a duration or the mark that stands for a missing one."""
    hole = f'<span class="hole-none">{NO_NUMBER}</span>'
    columns = time.columns()
    return "".join(f'<td class="num">{columns[i] or hole}</td>' for i in kept)


def _and(words: Sequence[str]) -> str:
    """Return a list of words as a sentence reads one: ``a``, ``a and b``, ``a, b and c``."""
    if len(words) < 2:
        return "".join(words)
    return f"{', '.join(words[:-1])} and {words[-1]}"


@dataclass(frozen=True, slots=True)
class _Time:
    """What one protocol, or a whole run, is known to take, and what nothing states."""

    held: float = 0.0
    hands_on: float = 0.0
    hands_said: int = 0
    unattended: float = 0.0
    unattended_said: int = 0
    steps: int = 0
    blank: int = 0
    waits: tuple[Wait, ...] = ()

    def and_(self, other: "_Time") -> "_Time":
        """Return the two added, which is how a run's row is the sum of its protocols'."""
        return _Time(
            self.held + other.held,
            self.hands_on + other.hands_on,
            self.hands_said + other.hands_said,
            self.unattended + other.unattended,
            self.unattended_said + other.unattended_said,
            self.steps + other.steps,
            self.blank + other.blank,
            self.waits + other.waits,
        )

    def columns(self) -> tuple[str | None, ...]:
        """Return this row under `SCHEDULE_COLUMNS`, each column a value or ``None`` for none."""
        return (
            str(self.steps),
            # The steps holding nothing have a column of their own, so held is not qualified.
            _taken(self.held, self.steps - self.blank),
            _taken(self.hands_on, self.hands_said, self.steps),
            _taken(self.unattended, self.unattended_said, self.steps),
            str(self.blank),
        )


def _taken(seconds: float, said: int, steps: int = 0) -> str | None:
    """One duration: the time, over how many steps it was stated, or ``None`` where none did.

    `steps` is left out where the row already says elsewhere how many steps stated nothing.
    """
    if not said:
        return None
    over = (
        f' <span class="muted">over {said} of {_count(steps, "step")}</span>'
        if 0 < said < steps
        else ""
    )
    return f"{_duration(seconds)}{over}"


def _waiting(waits: tuple[Wait, ...]) -> str:
    """Return the wait row: what is waited on and for how long, or that nothing recorded any."""
    if not waits:
        return '<span class="muted">No waiting recorded.</span>'
    return _waits(waits)


def _time_of(protocol: Protocol) -> _Time:
    """Return what one protocol is known to take, read from the steps and nowhere else."""
    held, blank = protocol.held_seconds
    hands = [step.hands_on_seconds for step in protocol.steps if step.hands_on_seconds is not None]
    spare = [
        step.held_seconds - step.hands_on_seconds
        for step in protocol.steps
        if step.held_seconds is not None and step.hands_on_seconds is not None
    ]
    return _Time(
        held,
        sum(hands),
        len(hands),
        sum(spare),
        len(spare),
        len(protocol.steps),
        blank,
        tuple(wait for step in protocol.steps for wait in step.waits),
    )


def _run_holes(project: Project, folder: Folder) -> tuple[Hole, ...]:
    """Every hole the run carries, the bill's first, then each protocol's, one id counted once."""
    return tuple(found for found, _ in _holes_found(project, folder).values())


def _holes_found(project: Project, folder: Folder) -> dict[str, tuple[Hole, tuple[str, str]]]:
    """Each hole by id, with the page holding it and what that page is called."""
    found: dict[str, tuple[Hole, tuple[str, str]]] = {}
    for row in project.bill.rows if project.bill else ():
        if row.hole:
            found.setdefault(row.hole.id, (row.hole, (folder.reagents, "Reagents and equipment")))
    for page, protocol in zip(folder.pages, project.protocols, strict=True):
        for hole in protocol.all_holes:
            found.setdefault(hole.id, (hole, (page.href, page.title)))
    return found


def _run_hole_list(project: Project, folder: Folder, holes: tuple[Hole, ...]) -> str:
    """Every hole the run carries, each linking to the page it stands on."""
    if not holes:
        return ""
    found = _holes_found(project, folder)
    items = "".join(
        _hole(
            hole,
            f' <span class="hole-page">on <a href="{escape(found[hole.id][1][0])}'
            f'#hole-{escape(hole.id)}">{escape(found[hole.id][1][1])}</a></span>',
        )
        for hole in holes
    )
    return (
        '<section class="block holes" id="holes">\n<h2>Holes</h2>\n'
        f"<p>{_count(len(holes), 'number')} this run would otherwise have to invent. A hole is "
        "a defect in what the package knows, not a failure of the run, and it is never filled "
        f"with a guess.</p>\n<ul>{items}</ul>\n</section>\n"
    )


def _merged_materials(project: Project) -> str:
    """Every material the run takes, merged by name and catalogue number.

    Two protocols buying one thing is one row on a shopping list. Where they state different
    amounts, both stand, because an amount is prose and nothing can add two of them.
    """
    found: dict[tuple[str, str], Material] = {}
    takers: dict[tuple[str, str], list[str]] = {}
    for protocol in project.protocols:
        for material in protocol.materials:
            key = (material.name, material.catalog)
            kept = found.get(key)
            found[key] = material if kept is None else replace(kept, amount=_join(kept, material))
            takers.setdefault(key, []).append(protocol.title)
    used = tuple(", ".join(dict.fromkeys(names)) for names in takers.values())
    return _materials(tuple(found.values()), (), used, note=False)


def _join(kept: Material, found: Material) -> str:
    """Both amounts where two protocols take different ones, and the one where they agree."""
    return "; ".join(dict.fromkeys(one for one in (kept.amount, found.amount) if one))


def _kit(project: Project) -> str:
    """Return the hardware, vessels and plates the run needs, each naming who needs it."""
    takers: dict[tuple[str, str], list[str]] = {}
    for protocol in project.protocols:
        named = [
            *((one, "") for one in protocol.equipment),
            *(
                (v.name, " · ".join(one for one in (v.kind, v.catalog, v.holds) if one))
                for v in protocol.vessels
            ),
            *(
                (p.name, " · ".join(one for one in (f"{p.wells} wells", p.catalog) if one))
                for p in protocol.plates
            ),
        ]
        for one in named:
            takers.setdefault(one, []).append(protocol.title)
    if not takers:
        return ""
    items = "".join(
        f'<li class="kit"><strong>{escape(name)}</strong>'
        + (f' <span class="muted">{escape(what)}</span>' if what else "")
        + f' <span class="used-in">used in {escape(", ".join(dict.fromkeys(names)))}</span></li>'
        for (name, what), names in takers.items()
    )
    return (
        '<section class="block kit-list" id="kit">\n<h2>Equipment and plasticware</h2>\n'
        f'<ul class="vessels">{items}</ul>\n</section>\n'
    )


def _merged_references(project: Project) -> tuple[tuple[Reference, ...], tuple[str, ...]]:
    """Every reference the run holds, once each, with the protocols holding it."""
    holders: dict[Reference, list[str]] = {}
    for protocol in project.protocols:
        for reference in protocol.references:
            holders.setdefault(reference, []).append(protocol.title)
    return tuple(holders), tuple(", ".join(dict.fromkeys(n)) for n in holders.values())


def _merged_sources(project: Project) -> tuple[dict[str, Source], dict[str, str]]:
    """Every document a number was read from, once each, with the protocols citing it."""
    found: dict[str, Source] = {}
    citers: dict[str, list[str]] = {}
    for protocol in project.protocols:
        for key, source in protocol.sources.items():
            found.setdefault(key, source)
            citers.setdefault(key, []).append(protocol.title)
    return found, {key: ", ".join(dict.fromkeys(names)) for key, names in citers.items()}


def write_project_files(project: Project, directory: str | os.PathLike[str]) -> ProjectFiles:
    """Write `project` into `directory` as one data file and the pages rendered from it.

    The directory is made when it is not there, and every page is rendered from the data as
    written, so the two cannot disagree. Every page is computed before any is written, so each
    knows every other's title, address, step count and key. The same project writes the same
    bytes.
    """
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    data = write_project(project, out / PROJECT_DATA_FILE)
    written = read_project(data)
    folder = Folder(
        tuple(
            Page(one.title, page_name(n, one.title), len(one.steps), page_key(one))
            for n, one in enumerate(written.protocols, 1)
        )
    )
    protocols = tuple(
        write_html(one, out / page.href, folder=folder)
        for one, page in zip(written.protocols, folder.pages, strict=True)
    )
    return ProjectFiles(
        data,
        _write(render_index(written, folder), out / folder.index),
        protocols,
        _write(render_reagents(written, folder), out / folder.reagents),
        _write(render_references(written, folder), out / folder.references),
    )


def _write(html: str, path: Path) -> Path:
    path.write_text(html, encoding="utf-8")
    return path


def _page(
    title: str, key: str, body: str, folder: Folder | None, here: str, within: str = ""
) -> str:
    """Wrap one page's body in the document, and in the run's frame where it is part of one."""
    frame = main = f'<main class="page">\n{body}</main>\n'
    if folder is not None:
        right = within or '<div class="column within"></div>'
        frame = (
            f'{_bar(folder, here)}<div class="frame">{_chain(folder, here)}{main}{right}</div>\n'
        )
    return (
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{escape(title)}</title>\n"
        f"<style>\n{_fonts(frame)}{_asset('protocol.css')}</style>\n"
        f'</head>\n<body data-protocol="{key}">\n{frame}'
        f"<script>\n{_asset('protocol.js')}</script>\n</body>\n</html>\n"
    )


def _bar(folder: Folder, here: str) -> str:
    """Return the run's nav bar: the four places a page can go, the one it is on marked."""
    addresses = [page.href for page in folder.pages]
    protocols = here if here in addresses else (addresses[0] if addresses else folder.index)
    items = (
        (folder.index, "Overview"),
        (protocols, "Protocols"),
        (folder.reagents, "Reagents and equipment"),
        (folder.references, "References"),
    )
    links = "".join(
        f'<a href="{escape(href)}"{' aria-current="page"' if href == here else ""}>'
        f"{escape(label)}</a>"
        for href, label in items
    )
    return f'<nav class="site" aria-label="This run">{links}</nav>\n'


def _chain(folder: Folder, here: str) -> str:
    """Return the left column: every protocol of the run, so one page switches to another."""
    if not folder.pages:
        return '<div class="column chain"></div>'
    items = "".join(
        f"<li{' class="is-here"' if page.href == here else ''}>"
        f'<a href="{escape(page.href)}">{escape(page.title)}</a></li>'
        for page in folder.pages
    )
    return (
        '<nav class="column chain" aria-label="Protocols">'
        f'<h2 class="column-title">Protocols</h2><ol>{items}</ol></nav>\n'
    )


def _within(protocol: Protocol) -> str:
    """Return the right column: every step of this page, so a reader jumps within it."""
    if not protocol.steps:
        return ""
    links = "".join(
        f'<li><a href="#step-{n}">{escape(step.title)}</a></li>'
        for n, step in enumerate(protocol.steps, 1)
    )
    return (
        '<nav class="column within" aria-label="Steps">'
        f'<h2 class="column-title">This protocol</h2><ol>{links}</ol></nav>\n'
    )


def _place(folder: Folder, here: str) -> str:
    """Where this page stands in the run, named in words, because it prints with the page."""
    addresses = [page.href for page in folder.pages]
    if here not in addresses:
        return ""
    at, total = addresses.index(here), len(folder.pages)
    before = f"Comes after {escape(folder.pages[at - 1].title)}." if at else "The run starts here."
    after = (
        f"Next is {escape(folder.pages[at + 1].title)}." if at + 1 < total else "The run ends here."
    )
    return f'<p class="neighbours">Protocol {at + 1} of {total}. {before} {after}</p>\n'


def _asset(name: str) -> str:
    return files("liulab_mbio.protocol").joinpath(name).read_text(encoding="utf-8")


def _fonts(body: str) -> str:
    """Each face a drawing in `body` was measured in, written once, or nothing where it draws none.

    A drawing pins every line of text to the width the package measured it at, so a page drawing
    one as text carries those faces or the browser stretches its own to fit. Every figure holding
    a drawing is marked `drawing`, whatever it draws, so the faces follow it onto any page.
    """
    drawn = 'class="drawing' in body
    return "".join(font_face(font) for font in (SANS, BOLD, MONO)) if drawn else ""


def _duration(seconds: float) -> str:
    """One duration in words, to the second under an hour and to the minute from an hour up.

    Nothing a bench reads in hours is planned to the second, and the seconds of a run's total
    are a precision the numbers behind it do not have.
    """
    if round(seconds) >= 3600:
        hours, minutes = divmod(round(seconds / 60), 60)
        return f"{hours} h {minutes} min" if minutes else f"{hours} h"
    minutes, secs = divmod(round(seconds), 60)
    parts = [f"{minutes} min"] if minutes else []
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


def _header(protocol: Protocol, *, toc: bool = True, place: str = "") -> str:
    parts = [f'<header class="intro">\n<h1>{escape(protocol.title)}</h1>\n{place}']
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
            f'<div class="toolbar"><span class="progress" aria-live="polite">0 of '
            f"{_count(count, 'step')} done</span>"
            '<button type="button" class="print">Print</button>'
            '<button type="button" class="reset">Reset page</button></div>\n'
        )
        if toc:
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


def _materials(
    materials: tuple[Material, ...],
    equipment: tuple[str, ...],
    used: tuple[str, ...] = (),
    *,
    note: bool = True,
) -> str:
    """Everything that is not an oligo, and the hardware as one light line under it.

    `used` names, row by row, which protocols of a run take each material, for the page a whole
    run shares. A protocol's own list leaves it empty. That page is what is ordered rather than
    what is laid out, so it drops the bench note and `note` is how.
    """
    if not materials and not equipment:
        return ""
    columns: list[tuple[str, Callable[[Material], str]]] = [
        ("Supplier", lambda m: m.supplier),
        ("Catalogue", lambda m: m.catalog),
        ("Storage", lambda m: m.storage),
        ("Per run", lambda m: m.amount),
        *([("Note", lambda m: m.note)] if note else []),
    ]
    cited = any(m.citation for m in materials)
    shown = [(label, get) for label, get in columns if any(get(m) for m in materials)]
    table = ""
    if materials:
        head = (
            "<th>Name</th>"
            + "".join(f"<th>{escape(label)}</th>" for label, _ in shown)
            + ("<th>Source</th>" if cited else "")
            + ("<th>Used in</th>" if used else "")
        )
        rows = "".join(
            f"<tr><td>{escape(material.name)}</td>"
            + "".join(f"<td>{escape(get(material))}</td>" for _, get in shown)
            + (f"<td>{_cite(material.citation)}</td>" if cited else "")
            + (f"<td>{escape(used[i])}</td>" if used else "")
            + "</tr>"
            for i, material in enumerate(materials)
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
        '<section class="block materials" id="materials">\n<h2>Materials</h2>\n'
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


def _hole_count(holes: tuple[Hole, ...], subject: str = "this protocol") -> str:
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
        f"{'number' if one else 'numbers'} in {subject} {'has' if one else 'have'} "
        f"no source</a>{split}.</p>\n"
    )


def _hole(hole: Hole, found: str = "") -> str:
    """One hole, which reads as a hole and never as a value.

    `Hole.issue` is not printed: the bench page is read by someone who cannot open a tracker.
    `found` is markup naming where the hole stands, for a page that is not the one holding it.
    """
    where = f"{escape(hole.where)}: " if hole.where else ""
    filled = f" <em>Filled by {escape(hole.filled_by)}.</em>" if hole.filled_by else ""
    return (
        f'<li class="hole" id="hole-{escape(hole.id)}"><span class="hole-id">{escape(hole.id)}'
        f'</span> <span class="hole-none">{NO_NUMBER}</span> — {where}{escape(hole.missing)} '
        f'<span class="hole-kind">{escape(HOLE_KINDS[hole.kind])}</span>{filled}{found}</li>'
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
    # The count stands for the legend where the plate holds more kinds than a legend can list.
    kinds = _count(drawn.layout.kinds, "kind") if drawn.layout.kinds else ""
    facts = " · ".join(
        text for text in (f"{one.wells} wells", kinds, one.catalog, one.holds, one.note) if text
    )
    return (
        f'<figure class="drawing plate" data-plate="{escape(one.name)}">{drawn.element()}'
        f'<figcaption>{escape(one.name)} <span class="muted">{escape(facts)}</span></figcaption>'
        + (f'<ul class="plate-legend">{legend}</ul>' if legend else "")
        + "</figure>\n"
    )


def _figure(figure: Figure, base: Path, where: str) -> str:
    """Return the figure's records drawn as maps and inlined, as `_plate` inlines a plate.

    Laid out here and never stored, so the figure follows the design it is drawn from. Several
    records stack as rows in the order they are named, each labelled by the record's own name,
    which is how a figure shows one molecule becoming the next. One record draws as it did.

    A highlight lights each row that answers to it and dims every other row whole, so a name only
    the lit round's record carries lights that round. A name no record answers to is a mistake.

    Raises
    ------
    FileNotFoundError
        If a record is not there, naming the step and the path.
    ValueError
        If no record of the figure answers to a name the highlight lights.
    """
    rows = [_row(figure, base / named, where) for named in figure.records]
    lit = frozenset().union(*(answering for _, answering in rows))
    unlit = [name for name in figure.highlight if name.casefold() not in lit]
    if unlit:
        listed = ", ".join(repr(name) for name in unlit)
        raise ValueError(f"{where}: no record of this figure draws {listed} to highlight")
    if len(rows) == 1:
        ((element, _),) = rows
        drawn, stacked = element, ""
    else:
        drawn = "".join(element for element, _ in rows)
        stacked = " rows"
    return (
        f'<figure class="drawing map{stacked}">{drawn}'
        f"<figcaption>{escape(figure.caption)}{_after(figure.citation)}</figcaption></figure>\n"
    )


def _row(figure: Figure, path: Path, where: str) -> tuple[str, frozenset[str]]:
    """Return one record of a figure as its SVG element, and every name it answers to.

    A stacked row carries the record's own name beside it, so a reader knows which molecule it is.
    """
    if not path.is_file():
        raise FileNotFoundError(f"{where}: no record at {path} to draw")
    drawn = draw_map(
        path,
        region=figure.span,
        linear=figure.linear,
        sequence_view=figure.sequence_view,
        enzymes=figure.enzymes,
    )
    if figure.highlight:
        # Lit here rather than by `draw_map`, which refuses a name its one record does not draw:
        # across rows a name belongs to the row it names, and dims every other row whole.
        drawn = replace(drawn, highlight=figure.highlight)
    answering = drawn.answering
    element = drawn.element()
    if len(figure.records) == 1:
        return element, answering
    label = f'<p class="row-name">{escape(drawn.record.name)}</p>'
    return f'<div class="row">{label}{element}</div>', answering


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
        f'<section class="block bill" id="bill">\n<h2>{escape(bill.title or "Bill")}</h2>\n'
        "<p>Quantities come from the design and are here whatever is loaded. A charge comes "
        "only from a price record; where none prices a row, the money is a hole and no figure "
        "is estimated. Price steers no part of the design.</p>\n"
        f'<div class="scroll"><table><thead><tr>{head}</tr></thead>'
        f"<tbody>{''.join(rows)}</tbody>{total}</table></div>\n</section>\n"
    )


def _sources(sources: Mapping[str, Source], cited: Mapping[str, str] | None = None) -> str:
    """Every document a number was read from, so a citation resolves on the page itself.

    `cited` names, key by key, which protocols of a run cite each, for the page a run shares.
    """
    if not sources:
        return ""
    items = "".join(
        f'<li id="source-{escape(_slug(key))}"><strong>{escape(key)}</strong> '
        f"{escape(source.document)}"
        + "".join(
            f" · {escape(text)}"
            for text in (source.edition, source.read_as, source.date, source.note)
            if text
        )
        + (
            f' <a href="{escape(source.url)}" rel="noreferrer">{escape(source.url)}</a>'
            if source.url
            else ""
        )
        + (f' <span class="cited-by">cited by {escape(cited[key])}</span>' if cited else "")
        + "</li>"
        for key, source in sources.items()
    )
    return (
        '<section class="block sources" id="sources">\n<h2>Sources</h2>\n'
        f"<ul>{items}</ul>\n</section>\n"
    )


def _step(n: int, step: Step, protocol: Protocol, base: Path, section: str = "") -> str:
    key = f"step-{n}"
    parts = [f'<p class="step-section">{escape(section)}</p>\n'] if section else []
    parts += [
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
    parts += [_figure(f, base, f"step {n} {step.title!r}") for f in step.figures]
    parts += [_table(f"{key}-table-{i}", t) for i, t in enumerate(step.tables, 1)]
    parts += [_program(p) for p in step.programs]
    parts += [_transfer(t) for t in step.transfers]
    if step.holes:
        items = "".join(_hole(hole) for hole in step.holes)
        parts.append(f'<ul class="holes-here" aria-label="Holes">{items}</ul>\n')
    if step.timers:
        timers = "".join(_timer(f"{key}-timer-{i}", t) for i, t in enumerate(step.timers, 1))
        parts.append(f'<div class="timers">{timers}</div>\n')
    parts.append(_waits(step.waits))
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


def _waits(waits: tuple[Wait, ...]) -> str:
    """Return what the step waits on and for how long, where the waiting falls.

    How long is whatever the vendor states, in their words; an unstated turnaround reads as a
    hole, because a wait nobody has timed is not a wait of no time.
    """
    if not waits:
        return ""
    items = "".join(
        f'<li class="wait">Waiting on {escape(wait.what)}: '
        + (
            escape(wait.duration)
            if wait.duration
            else f'<span class="hole-none">{NO_NUMBER}</span>'
        )
        + f"{_after(wait.citation)}</li>"
        for wait in waits
    )
    return f'<ul class="waits" aria-label="Waiting">{items}</ul>\n'


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
    per_tube = [c for c in table.components if not c.master_mix]
    # One tube takes every component; the mix holds only those the column adds up, so where the
    # two totals count different things the mix total says which it is.
    only = '<br><span class="muted">mix only</span>' if per_tube and in_mix else ""
    head = (
        "<th>Component</th>"
        + ("<th>Stock</th>" if stock else "")
        + ("<th>Final</th>" if final else "")
        + '<th class="num">1 rxn (µL)</th>'
        + f'<th class="num">Mix for <span class="rxn-n">{table.reactions}</span> (µL)</th>'
    )
    foot = (
        f'<tr><th>Total</th>{blanks}<td class="num">{number(total)}</td>'
        f'<td class="num mix"><span data-ul="{in_mix!r}">{number(in_mix * scale)}</span>'
        f"{only}</td></tr>"
    )
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


def _timer(key: str, timer: Timer) -> str:
    """One timer, keyed so `protocol.js` can give it back its deadline after a page turn."""
    return (
        f'<button type="button" class="timer" data-key="{key}" data-seconds="{timer.seconds!r}">'
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


def _references(references: tuple[Reference, ...], cited: tuple[str, ...] = ()) -> str:
    """Return the reading behind the run; `cited` names which protocols hold each of them."""
    if not references:
        return ""
    items = "".join(
        f"<li>{escape(r.text)}"
        + (f' <a href="{escape(r.url)}" rel="noreferrer">{escape(r.url)}</a>' if r.url else "")
        + (f' <span class="cited-by">cited by {escape(cited[i])}</span>' if cited else "")
        + "</li>"
        for i, r in enumerate(references)
    )
    return (
        '<section class="block references" id="references">\n<h2>References</h2>\n'
        f"<ol>{items}</ol>\n</section>\n"
    )

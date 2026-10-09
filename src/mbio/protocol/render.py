"""Render a protocol to one self-contained HTML page, and a project to a folder of them.

A protocol on its own is one page that loads nothing. A run of several is a folder: one page
each, an index, and two pages the run shares. Every link between them is relative and every
page is still self-contained, so the folder opens from disk and survives being zipped.

`write_run_files` is the way in for a pipeline: the shape follows the chain's length, so no
pipeline chooses whether it writes a page or a folder.
"""

import hashlib
import os
import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from html import escape
from importlib.resources import files
from itertools import groupby
from pathlib import Path
from typing import NamedTuple

from mbio.checks import Status
from mbio.plot.drawing import draw_map, draw_plate
from mbio.plot.fonts import BOLD, MONO, SANS
from mbio.plot.page import font_face
from mbio.protocol.model import (
    Bill,
    Caution,
    Check,
    Citation,
    Figure,
    Gel,
    Hole,
    Incubation,
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
    Stamp,
    Step,
    ThermocyclerProgram,
    Timer,
    Topic,
    Transfer,
    Wait,
    by_place,
    number,
    read_project,
    read_protocol,
    slug,
    well_at,
    write_project,
    write_protocol,
)

#: What the page reads where a check carries no verdict, so it is never taken for a pass.
NO_VERDICT = "not judged"

#: What stands where a number would, so a hole can never be read as a figure.
NO_NUMBER = "no sourced number"

#: What the references page prints as citing a source the run itself names.
CITED_BY_BILL = "the bill"

#: What one protocol written alone is called, as the data and as the page rendered from it.
PROTOCOL_DATA_FILE = "protocol.json"
PROTOCOL_FILE = "protocol.html"

#: What a project folder calls the data every one of its pages is rendered from.
PROJECT_DATA_FILE = "project.json"

#: What a project folder calls the three pages that are not one protocol's.
INDEX_FILE = "index.html"
REAGENTS_FILE = "reagents.html"
REFERENCES_FILE = "references.html"

#: The directory a run writes its folder into, below the directory it was given. A folder and
#: not a flat set: it names its own `project.json`, and a build file may be called that already.
PROTOCOL_DIR = "protocol"

#: From how many oligos a sheet is summarised. Fewer than this reads as a sheet; past it the
#: rows repeat one pattern and the page states the pattern instead.
OLIGO_SUMMARY = 20

#: What each kind of hole says it is waiting on, as a sentence of its own. Each reads the same
#: for one hole and for a group of them, so a group says it once instead of once a hole.
HOLE_KINDS = {
    "undecided": "Waiting on a bench to settle it.",
    "unpublished": "Waiting on a number nobody has published.",
    "lab": "Waiting on the lab's own stock.",
    "unread": "Waiting on a source nobody has read.",
    "price": "Waiting on a price record.",
}

#: What the holes block says above its list, on a protocol's page and on a run's index alike.
HOLES_INTRO = (
    "None is filled with a guess, and each says below what it waits on: a bench, a shelf, a "
    "price, or a number nobody has published."
)


@dataclass(frozen=True, slots=True)
class Page:
    """One protocol's page, as every other page of its folder sees it.

    Parameters
    ----------
    title
        What the protocol is called, which is how every other page names it.
    href
        The file it was written to, relative to the folder.
    step_keys
        One key per numbered step, as that page addresses it, so another page of the folder
        counts the marks rather than trusting a number written twice.
    key
        What the page remembers its check marks under, so another page of the folder can read
        how far the bench got.
    choice
        The job the protocol is one way of doing, carried from it, so every page of the folder
        can see which pages share one place in the run.
    """

    title: str
    href: str
    step_keys: tuple[str, ...] = ()
    key: str = ""
    choice: str = ""

    @classmethod
    def of(cls, place: int, protocol: Protocol) -> "Page":
        """Return the page `protocol` is written to, standing `place` in its run, from one.

        Everything the other pages see is derived here, so an index addresses exactly what the
        page it names wrote. Two ways of one job stand in one place, so one `place` writes two
        file names and the title is what tells them apart.
        """
        return cls(
            protocol.title,
            page_name(place, protocol.title),
            _step_keys(protocol.steps),
            page_key(protocol),
            protocol.choice,
        )


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
    explained
        What the overview's own topics are titled. A way's page links the guidance on its job
        only where a topic titles that job, so no page offers a reader a link to nothing.
    """

    pages: tuple[Page, ...]
    index: str = INDEX_FILE
    reagents: str = REAGENTS_FILE
    references: str = REFERENCES_FILE
    explained: tuple[str, ...] = ()

    @classmethod
    def of(cls, project: Project) -> "Folder":
        """Return the folder `project` is written into: one page per protocol, in its place.

        The ways of one job share a place, so this is the one spot the run's numbering is
        derived, and every page of the folder is addressed by what it derived.
        """
        return cls(
            tuple(
                Page.of(at, one)
                for at, group in enumerate(by_place(project.protocols, lambda one: one.choice), 1)
                for one in group
            ),
            explained=tuple(topic.title for topic in project.background),
        )

    @property
    def places(self) -> tuple[tuple[Page, ...], ...]:
        """The pages one place of the run at a time: the ways of one job together."""
        return by_place(self.pages, lambda page: page.choice)

    def sources_at(self, here: str) -> str:
        """Return the page a citation on `here` finds its source on, empty where that is `here`.

        A protocol page lists its own sources and `references` lists the whole run's, so a
        citation on either resolves where it stands; `Protocol.audit` and `Project.audit` judge
        whether every source cited was named. Every other page the run shares carries no sources
        list, and its citations reach the run's.

        Examples
        --------
        >>> folder = Folder((Page("Build the blocks", "01-build-the-blocks.html"),))
        >>> folder.sources_at("01-build-the-blocks.html")
        ''
        >>> folder.sources_at(folder.reagents)
        'references.html'
        """
        own = here == self.references or any(page.href == here for page in self.pages)
        return "" if own else self.references


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


@dataclass(frozen=True, slots=True)
class ProtocolFiles:
    """The pair one protocol wrote on its own: the data, and the page rendered from it.

    Parameters
    ----------
    data
        The bench protocol as JSON, which ``protocol render`` turns back into the page.
    page
        The interactive bench protocol, as one self-contained HTML page.
    """

    data: Path
    page: Path

    @property
    def paths(self) -> tuple[Path, ...]:
        """Both, in the order they were written."""
        return (self.data, self.page)


@dataclass(frozen=True, slots=True)
class RunFiles:
    """What one run wrote, whichever shape its chain length called for.

    Parameters
    ----------
    data
        The run as JSON, which ``protocol render`` turns back into the pages.
    pages
        What was written beside it: the one page, where the run is one protocol; the index, a
        page per protocol and the two shared pages, where it is a chain. The way in leads them,
        which `page` reads, so a folder's pages are not in the order they were written --
        `ProjectFiles.paths` is the one that reports that.
    """

    data: Path
    pages: tuple[Path, ...]

    @property
    def page(self) -> Path:
        """The way in: the only page where the run is one protocol, the index where it is many."""
        return self.pages[0]

    @property
    def paths(self) -> tuple[Path, ...]:
        """Every file written: the data, then the pages with the way in leading them."""
        return (self.data, *self.pages)


def page_key(value: Protocol | Project) -> str:
    """Return what a page remembers its check marks under: the key it carries, or its content.

    Two pages of one folder never share it, and every `file://` page in a browser shares one
    store, so the key is what keeps one page's marks off another. A key already minted is kept,
    which is what carries the bench's ticks across an agent's edit; everything else is digested,
    so a re-planned run starts clean.

    Examples
    --------
    >>> page_key(Protocol("Demo")) == page_key(Protocol("Demo"))
    True
    >>> page_key(Protocol("Demo", key="whichever"))
    'whichever'
    """
    return value.key or hashlib.sha256(repr(value).encode()).hexdigest()[:16]


def minted[T: (Protocol, Project)](value: T) -> T:
    """Return `value` carrying the key its page remembers the bench's marks under, and no other.

    A key already on it is left alone, and a project mints one for each of its protocols too. A
    pipeline mints before it writes the JSON, so the file an agent edits names its own store:
    `docs/adr/0002-editable-protocols.md` draws the two paths.

    Examples
    --------
    >>> minted(Protocol("Demo")).key == page_key(Protocol("Demo"))
    True
    """
    if isinstance(value, Project):
        value = replace(value, protocols=tuple(minted(one) for one in value.protocols))
    return value if value.key else replace(value, key=page_key(value))


def page_name(place: int, title: str) -> str:
    """Return the file one protocol's page is written to: its place in the chain, then its title.

    Examples
    --------
    >>> page_name(2, "LR reaction")
    '02-lr-reaction.html'
    """
    name = slug(title)
    return f"{place:02d}-{name}.html" if name else f"{place:02d}.html"


def _step_keys(steps: Sequence[Step]) -> tuple[str, ...]:
    """Return what each step is addressed by on its page, no two steps alike.

    A page must render, so a key another step has taken — or a title that slugs to nothing —
    takes the step's number and keeps taking one until the key is the page's own, rather than
    being refused.
    """
    taken: set[str] = set()
    keys = []
    for n, step in enumerate(steps, 1):
        key = step.key or str(n)
        while key in taken:
            key = f"{key}-{n}"
        taken.add(key)
        keys.append(key)
    return tuple(keys)


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
    keys = _step_keys(protocol.steps)
    # A section labels the steps under it, so it is written once, where it changes.
    sections = [""] + [step.section for step in protocol.steps]
    body = "".join(
        [
            _header(protocol, keys, toc=folder is None, place=place),
            _background(protocol.background, protocol.files),
            _materials(protocol.materials, protocol.equipment, paths=protocol.files),
            _oligos(protocol),
            _plates(protocol),
            _bill(protocol.bill),
            *(
                _step(
                    n,
                    step,
                    keys[n - 1],
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
    return _page(protocol.title, page_key(protocol), body, folder, here, _within(protocol, keys))


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
    step never stops to explain a decision. The badges carry the run's own verdicts beside the
    chain's, computed here rather than read from the data, so an edited file is judged as edited
    (`docs/adr/0002-editable-protocols.md`).
    """
    sources = folder.sources_at(folder.index)
    holes = _run_holes(project, folder)
    said = _linked(project.summary, project.files)
    summary = f'<p class="summary">{said}</p>\n' if project.summary else ""
    jumps = [(f"topic-{slug(topic.title)}", topic.title) for topic in project.background]
    parts = [
        f'<header class="intro">\n<h1>{escape(project.title)}</h1>\n{summary}'
        f"{_checks((*project.checks, *project.audit()))}"
        f"{_hole_count(holes, 'this run')}</header>\n",
        _background(project.background, project.files),
    ]
    for anchor, label, block in (
        ("flow", "How the run fits together", _flow(project, folder)),
        ("protocols", "Protocols", _protocols(folder)),
        ("schedule", "Schedule", _schedule(project, folder, sources)),
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
    sources = folder.sources_at(folder.reagents)
    blocks = (
        ("materials", "Materials", _merged_materials(project, sources)),
        ("kit", "Equipment and plasticware", _kit(project)),
        ("bill", "Bill", _bill(project.bill, sources)),
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
        page_key(project) + slug(here),
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


def _background(topics: Sequence[Topic], paths: Sequence[str] = ()) -> str:
    """Return what the reader is told before the first step, one topic to a block.

    A run says it on its index and a page written alone says it on itself, so one renderer
    serves both.
    """
    return "".join(
        f'<section class="block topic" id="topic-{escape(slug(topic.title))}">\n'
        f"<h2>{escape(topic.title)}</h2>\n"
        + "".join(f"<p>{_linked(line, paths)}</p>" for line in topic.body)
        + "\n</section>\n"
        for topic in topics
    )


def _protocols(folder: Folder) -> str:
    """Return the run a place at a time, each page carrying its store's key and its steps' keys.

    `protocol.js` reads that store on this page and counts the marks those step keys stand for,
    so how far the bench got comes from the pages themselves and is never stored twice. A place
    holding several ways of one job is one entry, with the ways nested under it.
    """
    if not folder.pages:
        return ""
    rows = "".join(_protocols_place(group) for group in folder.places)
    return (
        '<section class="block chain-list" id="protocols">\n<h2>Protocols</h2>\n'
        f'<ol class="chain-pages">{rows}</ol>\n</section>\n'
    )


def _protocols_place(group: tuple[Page, ...]) -> str:
    """One place of the run as the index lists it: one page, or the ways of one job under it."""
    ways = "".join(
        f'<li data-page-key="{escape(page.key)}" data-steps="{escape(" ".join(page.step_keys))}">'
        f'<a href="{escape(page.href)}">{escape(page.title)}</a> '
        f'<span class="page-progress muted">{_count(len(page.step_keys), "step")}</span></li>'
        for page in group
    )
    if not group[0].choice:
        return ways
    return (
        f'<li class="chain-choice"><span class="choice-job">{escape(_opening(group[0].choice))}'
        f'</span> <span class="muted">do {_one_of(len(group))}</span>'
        f'<ul class="chain-ways">{ways}</ul></li>'
    )


def _opening(job: str) -> str:
    """Return a job as it opens a line: it is written to sit inside a sentence, so lowercase."""
    return job[:1].upper() + job[1:]


#: How many ways a page spells out rather than printing as a figure. A count a reader never
#: counts out reads as a word, and every page of one run says it the same way.
SPELLED = ("", "", "two", "three", "four", "five")


def _spelled(ways: int) -> str:
    """Return how many ways a page says there are."""
    return SPELLED[ways] if ways < len(SPELLED) else str(ways)


def _one_of(ways: int) -> str:
    """Return how a page tells the bench to do one way of a job and not the others."""
    return "this one" if ways < 2 else f"one of these {_spelled(ways)}"


def _flow(project: Project, folder: Folder) -> str:
    """Return the chain as boxes and the names handed between them, in HTML and never a drawing.

    Generated from what each protocol declares it consumes and produces: a box is one place of
    the run, and above it stand the names it is handed, each saying where it came from. A place
    offering several ways of one job is one box holding them all, and only a name every way
    makes stands below it. Drawn in text, so a browser's find reaches it, it is selectable and
    it prints with no font embedded.
    """
    if not folder.pages:
        return ""
    came_from = {item.name: "the bench already holds it" for item in project.inputs}
    left: dict[str, Item] = {}
    boxes = []
    paired = tuple(zip(folder.pages, project.protocols, strict=True))
    for group in by_place(paired, lambda pair: pair[1].choice):
        wanted: dict[str, Item] = {}
        made: dict[str, Item] = {}
        for _, protocol in group:
            for item in protocol.consumes:
                wanted.setdefault(item.name, item)
            for item in protocol.produces:
                made.setdefault(item.name, item)
        needs = "".join(_flow_item(item, came_from.get(name, "")) for name, item in wanted.items())
        for name in wanted:
            left.pop(name, None)
        hand = f'<ul class="flow-hand">{needs}</ul>' if needs else ""
        # The step count stands in the Protocols list, where `protocol.js` keeps it up to date;
        # printed here too it would be the same number twice, one of them stale.
        boxes.append(f"<li>{hand}{_flow_box(group)}</li>")
        first, job = group[0][0], group[0][1].choice
        source = f"from whichever way to {job} you did" if job else f"from {first.title}"
        # Only a name every way makes stands below the box: the bench did one of them.
        every = set.intersection(*({one.name for one in p.produces} for _, p in group))
        for name, item in made.items():
            if name in every:
                came_from[name] = source
                left[name] = item
    ends = "".join(_flow_item(item, "the run ends holding it") for item in left.values())
    return (
        '<section class="block flow" id="flow">\n<h2>How the run fits together</h2>\n'
        f'<ol class="flow-chain">{"".join(boxes)}</ol>'
        + (f'<ul class="flow-hand flow-end">{ends}</ul>' if ends else "")
        + "\n</section>\n"
    )


def _flow_box(group: tuple[tuple[Page, Protocol], ...]) -> str:
    """One place of the run as a box: a protocol, or the ways of one job to pick between."""
    job = group[0][1].choice
    if not job:
        page = group[0][0]
        return (
            f'<a class="flow-box" href="{escape(page.href)}">'
            f'<span class="flow-title">{escape(page.title)}</span></a>'
        )
    ways = "".join(
        f'<li><a href="{escape(page.href)}">{escape(page.title)}</a></li>' for page, _ in group
    )
    return (
        '<div class="flow-box flow-choice">'
        f'<span class="flow-title">{escape(_opening(job))}</span> '
        f'<span class="muted">do {_one_of(len(group))}</span>'
        f'<ul class="flow-ways">{ways}</ul></div>'
    )


def _flow_item(item: Item, came_from: str) -> str:
    """One name handed between two protocols, and where it came from, or that nothing hands it."""
    said = came_from or "nothing in the run hands this over"
    dangling = "" if came_from else " is-dangling"
    return (
        f'<li class="flow-hand-item{dangling}">{_item_named(item)} '
        f'<span class="flow-from">{escape(said)}</span></li>'
    )


def _item_named(item: Item) -> str:
    """One item as every page names it: what it is called, what it is, and what it has to meet."""
    wanted = escape(", ".join(item.spec))
    spec = f' <span class="muted">· {wanted}</span>' if item.spec else ""
    return (
        f'<span class="flow-item">{escape(item.name)}</span> '
        f'<span class="muted">{escape(item.what)}</span>{spec}'
    )


#: The schedule's columns after the protocol's own name. The waiting goes in a row of its own,
#: because how long a vendor takes is whatever they say and no column is wide enough for it.
SCHEDULE_COLUMNS = ("Steps", "Held", "Hands-on", "Unattended", "Holding nothing")

#: Which of those columns is a duration rather than a count, so a figure knows how to read.
SCHEDULE_DURATIONS = (1, 2, 3)


def _schedule(project: Project, folder: Folder, sources: str = "") -> str:
    """Return the run's time plan as a table, with a hole wherever nothing states a number.

    A table and not a bar chart: most of a run's calendar is time nobody attends and nobody can
    date, and a bar draws an unknown wait as a length, which is a claim. Held is summed from the
    timers and thermocycler programs the steps already hold; hands-on is stated or it is a hole,
    and the unattended share is the difference wherever a step states both. A column nothing
    states a number in is left out, and named once under the table.

    A row is a place of the run and not a page, so the table counts what the rest of the page
    counts. The ways of one job stand under their place, which spans each column across them.
    """
    if not folder.pages:
        return ""
    times = [_time_of(protocol) for protocol in project.protocols]
    totals = _Time()
    for time in times:
        totals = totals.and_(time)
    shown = tuple(i for i, column in enumerate(totals.columns()) if column is not None)
    places = by_place(tuple(zip(folder.pages, times, strict=True)), lambda pair: pair[0].choice)
    rows = "".join(_schedule_place(group, shown, sources) for group in places)
    head = "<th>Protocol</th>" + "".join(
        f'<th class="num">{escape(SCHEDULE_COLUMNS[i])}</th>' for i in shown
    )
    # Only a duration column can go: the steps and the steps holding nothing are always counted.
    missing = [SCHEDULE_COLUMNS[i].lower() for i in range(len(SCHEDULE_COLUMNS)) if i not in shown]
    left_out = (
        f" No protocol here states {_or(missing)} time, so "
        f"{'that column is' if len(missing) == 1 else 'those columns are'} left out."
        if missing
        else ""
    )
    spans = (
        " Where one job offers several ways, the row above them gives each column's range: you "
        "do one way, and that way's own row says what it takes."
        if any(len(group) > 1 for group in places)
        else ""
    )
    total, part = _total(places, totals, shown)
    return (
        '<section class="block schedule" id="schedule">\n<h2>Schedule</h2>\n'
        '<div class="scroll"><table class="schedule"><thead><tr>'
        f"{head}</tr></thead><tbody>{rows}</tbody>"
        f"<tfoot><tr><th>Total</th>{total}</tr></tfoot></table></div>\n"
        "<p>A blank here is a number nobody has stated, never a zero. Most of this run is time "
        "nobody attends and nobody can date, so it is written down rather than drawn as a "
        f"length.{left_out}{spans}{part}</p>\n</section>\n"
    )


#: One place of the schedule: its pages, each with what that page's protocol is known to take.
type _Place = Sequence[tuple[Page, "_Time"]]


def _schedule_place(group: _Place, shown: tuple[int, ...], sources: str) -> str:
    """Return one place of the run: a protocol's row, or the ways of one job under their span."""
    rows = "".join(
        f"<tr{' class="schedule-way"' if len(group) > 1 else ''}>"
        f'<td><a href="{escape(page.href)}">{escape(page.title)}</a></td>'
        f'{time.cells(shown)}</tr><tr class="wait-row"><td colspan="{len(shown) + 1}">'
        f"{_waiting(time.waits, sources)}</td></tr>"
        for page, time in group
    )
    if len(group) == 1:
        return rows
    job = escape(_opening(group[0][0].choice))
    return (
        f'<tr class="schedule-choice"><td><span class="choice-job">{job}</span></td>'
        f"{_across(group, shown)}</tr>{rows}"
    )


def _across(group: _Place, shown: tuple[int, ...]) -> str:
    """Return a place's own cells, each column spanned across its ways, low to high.

    Both ends are a way's own measurement, so the row invents nothing; what each way states in
    full stands in its own row under this one.
    """
    hole = f'<span class="hole-none">{NO_NUMBER}</span>'
    cells = []
    for i in shown:
        said = [value for _, time in group if (value := time.values()[i]) is not None]
        figure = (
            f"{_spanned(_bare(i, min(said)), _bare(i, max(said)))}"
            f"{_over(len(said), len(group), 'way')}"
            if said
            else hole
        )
        cells.append(f'<td class="num">{figure}</td>')
    return "".join(cells)


def _bare(column: int, value: float) -> str:
    """One column's own number, in the words that column is read in."""
    return _duration(value) if column in SCHEDULE_DURATIONS else f"{round(value)}"


def _spanned(low: str, high: str) -> str:
    """Return both ends of a span, or the one figure where they read alike."""
    return low if low == high else f"{low} to {high}"


def _total(places: Sequence[_Place], totals: "_Time", shown: tuple[int, ...]) -> tuple[str, str]:
    """Return the footer's cells, and the sentence to add under the table where one is partial.

    A total of a column some protocol states nothing in is not what the run takes, and a reader
    plans a week around it. So it carries how many places it covers, beside the number, in the
    count the rest of the page prints; a place is covered where every way of it states the
    column. Where a job offers ways, how many steps the run has is a span too, so the figures
    stand alone and each protocol's own row keeps the steps behind them.

    Where no place is covered in full the count reads zero beside a measured figure, and that
    stands: the hole mark would deny a number a way did state, and a bare figure would claim
    the whole run.
    """
    spanning = any(len(group) > 1 for group in places)
    cells = []
    partial = False
    for i in shown:
        low, high, covered = 0.0, 0.0, 0
        for group in places:
            said = [value for _, time in group if (value := time.values()[i]) is not None]
            covered += len(said) == len(group)
            low, high = low + min(said, default=0.0), high + max(said, default=0.0)
        partial = partial or covered < len(places)
        figure = (
            _spanned(_bare(i, low), _bare(i, high)) if spanning else (totals.columns()[i] or "")
        )
        cells.append(f'<td class="num">{figure}{_over(covered, len(places), "protocol")}</td>')
    part = " A total that covers only part of the run says so beside it." if partial else ""
    return "".join(cells), part


def _over(said: int, of: int, noun: str) -> str:
    """How much of a whole a figure was computed over, or nothing where it covers all of it."""
    return f' <span class="muted">over {said} of {_count(of, noun)}</span>' if said < of else ""


def _or(words: Sequence[str]) -> str:
    """Return a list of words as a sentence reads one: ``a``, ``a or b``, ``a, b or c``."""
    if len(words) < 2:
        return "".join(words)
    return f"{', '.join(words[:-1])} or {words[-1]}"


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

    def cells(self, shown: tuple[int, ...]) -> str:
        """Return this row's `shown` columns, each a value or the mark for a missing one."""
        hole = f'<span class="hole-none">{NO_NUMBER}</span>'
        columns = self.columns()
        return "".join(f'<td class="num">{columns[i] or hole}</td>' for i in shown)

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

    def values(self) -> tuple[float | None, ...]:
        """Return each column's bare number, or ``None`` where nothing stated one.

        What `columns` prints, before it is worded: a span and a total are arithmetic, so they
        read the numbers rather than the sentences built from them.
        """
        return (
            float(self.steps),
            self.held if self.steps > self.blank else None,
            self.hands_on if self.hands_said else None,
            self.unattended if self.unattended_said else None,
            float(self.blank),
        )


def _taken(seconds: float, said: int, steps: int = 0) -> str | None:
    """One duration: the time, over how many steps it was stated, or ``None`` where none did.

    `steps` is left out where the row already says elsewhere how many steps stated nothing.
    """
    if not said:
        return None
    return f"{_duration(seconds)}{_over(said, steps, 'step')}"


def _waiting(waits: tuple[Wait, ...], sources: str = "") -> str:
    """Return the wait row: what is waited on and for how long, or that nothing recorded any."""
    if not waits:
        return '<span class="muted">No waiting recorded.</span>'
    return _waits(waits, sources)


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
    """Each hole by id, with where on which page it stands and what that page is called.

    A protocol page prints the hole under its own id; the bill prints it as the row whose money
    it stands in, so that is what the run's index points a reader at.
    """
    found: dict[str, tuple[Hole, tuple[str, str]]] = {}
    for row in project.bill.rows if project.bill else ():
        if row.hole:
            where = (f"{folder.reagents}#bill", "Reagents and equipment")
            found.setdefault(row.hole.id, (row.hole, where))
    for page, protocol in zip(folder.pages, project.protocols, strict=True):
        for hole in protocol.all_holes:
            found.setdefault(hole.id, (hole, (f"{page.href}#hole-{hole.id}", page.title)))
    return found


def _run_hole_list(project: Project, folder: Folder, holes: tuple[Hole, ...]) -> str:
    """Every hole the run carries, each linking to the page it stands on.

    Holes waiting on one thing on one page are gathered, because nine of them saying it nine
    times is one fact told at nine times the length. None is dropped: the gathered ones keep
    their ids on the summary and their own lines behind the disclosure.
    """
    if not holes:
        return ""
    items = "".join(
        _hole_item(group, where) for group, where in _gathered(holes, _holes_found(project, folder))
    )
    return (
        '<section class="block holes" id="holes">\n<h2>Holes</h2>\n'
        f"<p>{_count(len(holes), 'number')} this run would otherwise have to invent. "
        f"{HOLES_INTRO}</p>\n<ul>{items}</ul>\n</section>\n"
    )


def _gathered(
    holes: tuple[Hole, ...], found: Mapping[str, tuple[Hole, tuple[str, str]]]
) -> list[tuple[tuple[Hole, ...], tuple[str, str]]]:
    """Return a group per thing the holes wait on, with the page it is waited on, first said first.

    The same kind filled by the same thing on the same page is one statement; anything else
    differs where it matters and stands on its own. A group holds its holes in the order given.
    """
    gathered: dict[tuple[str, str, str], tuple[list[Hole], tuple[str, str]]] = {}
    for hole in holes:
        where = found[hole.id][1]
        key = (hole.kind, hole.filled_by, where[0])
        gathered.setdefault(key, ([], where))[0].append(hole)
    return [(tuple(group), where) for group, where in gathered.values()]


def _hole_item(group: tuple[Hole, ...], where: tuple[str, str]) -> str:
    """One hole on its own line, or several behind a disclosure saying what they all wait on."""
    page = f'<a href="{escape(where[0])}">{escape(where[1])}</a>'
    if len(group) == 1:
        return _hole(group[0], f' <span class="hole-page">Stands on {page}.</span>')
    ids = ", ".join(escape(hole.id) for hole in group)
    items = "".join(f'<li class="hole" id="hole-{escape(h.id)}">{_missing(h)}</li>' for h in group)
    return (
        '<li class="hole"><details class="hole-group">\n<summary>'
        f'<span class="hole-id">{ids}</span> <span class="hole-none">{NO_NUMBER}</span> — '
        f"{_count(len(group), 'number')}, each waiting on the same thing. "
        f'<span class="hole-kind">{escape(HOLE_KINDS[group[0].kind])}</span>'
        f'{_filled_by(group[0])} <span class="hole-page">Each stands on {page}.</span>'
        f'</summary>\n<ul class="holes-here">{items}</ul>\n</details></li>'
    )


def _merged_materials(project: Project, sources: str = "") -> str:
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
    return _materials(tuple(found.values()), (), used, project.files, note=False, sources=sources)


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
    """Every document a number was read from, once each, with what cites it.

    The run's own come first, each cited by the bill where a row of it cites one;
    `docs/adr/0018-a-project-chains-protocols.md` says why a run names any source at all.
    """
    billed = project.bill.cited if project.bill else frozenset()
    found: dict[str, Source] = dict(project.sources)
    citers: dict[str, list[str]] = {
        key: [CITED_BY_BILL] for key in project.sources if key in billed
    }
    for protocol in project.protocols:
        for key, source in protocol.sources.items():
            found.setdefault(key, source)
            citers.setdefault(key, []).append(protocol.title)
    return found, {key: ", ".join(dict.fromkeys(names)) for key, names in citers.items()}


def write_project_files(project: Project, directory: str | os.PathLike[str]) -> ProjectFiles:
    """Write `project` into `directory` as one data file and the pages rendered from it.

    The directory is made when it is not there, and every page is rendered from the data as
    written, so the two cannot disagree. Every page is computed before any is written, so each
    knows every other's title, address, step keys and key. The key goes into the data, so an
    agent editing the file and rendering again hands the bench back its own ticks. The same
    project writes the same bytes.
    """
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    data = write_project(minted(project), out / PROJECT_DATA_FILE)
    written = read_project(data)
    folder = Folder.of(written)
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


def write_protocol_files(protocol: Protocol, directory: str | os.PathLike[str]) -> ProtocolFiles:
    """Write `protocol` into `directory` as `PROTOCOL_DATA_FILE` and `PROTOCOL_FILE`.

    The directory is made when it is not there, and the page is rendered from the data as
    written, so the two cannot disagree. The data carries the key its page remembers the bench's
    check marks under, so an agent editing it and rendering again keeps the ticks already made.
    The same protocol writes the same bytes.
    """
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    data = write_protocol(minted(protocol), out / PROTOCOL_DATA_FILE)
    return ProtocolFiles(data, write_html(read_protocol(data), out / PROTOCOL_FILE))


def write_run_files(project: Project, directory: str | os.PathLike[str]) -> RunFiles:
    """Write `project` into `directory`, in the shape its chain length calls for.

    One protocol writes one self-contained page, flat in `directory`; two or more write the
    folder of linked pages, in `PROTOCOL_DIR` below it. The shape follows the chain and never
    the pipeline, so no caller chooses it:
    `docs/adr/0018-a-project-chains-protocols.md` says why.
    """
    if len(project.protocols) == 1:
        alone = write_protocol_files(_alone(project), directory)
        return RunFiles(alone.data, (alone.page,))
    folder = write_project_files(project, Path(directory) / PROTOCOL_DIR)
    return RunFiles(
        folder.data, (folder.index, *folder.protocols, folder.reagents, folder.references)
    )


def _alone(project: Project) -> Protocol:
    """Return the one protocol of a run of one, carrying what the run says over it.

    A page written alone has no index to hold the run's own title, rationale, verdicts, bill or
    sources, so each travels onto the protocol rather than being dropped for having been said a
    level up. What the run was handed is already the protocol's `consumes`.
    """
    one = project.protocols[0]
    return replace(
        one,
        title=project.title or one.title,
        summary=project.summary or one.summary,
        background=project.background or one.background,
        checks=(*one.checks, *project.checks),
        files=one.files or project.files,
        sources={**project.sources, **one.sources},
        bill=one.bill or project.bill,
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
        f'</head>\n<body data-protocol="{escape(key)}">\n{frame}'
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
    """Return the left column: every protocol of the run, so one page switches to another.

    The ways of one job are one entry with the ways under it, as the index lists them, so the
    column numbers the places of the run and not the pages.
    """
    if not folder.pages:
        return '<div class="column chain"></div>'
    items = "".join(_chain_place(group, here) for group in folder.places)
    return (
        '<nav class="column chain" aria-label="Protocols">'
        f'<h2 class="column-title">Protocols</h2><ol>{items}</ol></nav>\n'
    )


def _chain_place(group: tuple[Page, ...], here: str) -> str:
    """One place of the run in the left column, the page being read marked."""
    ways = "".join(
        f"<li{' class="is-here"' if page.href == here else ''}>"
        f'<a href="{escape(page.href)}">{escape(page.title)}</a></li>'
        for page in group
    )
    if not group[0].choice:
        return ways
    return (
        f'<li class="chain-choice"><span class="choice-job">{escape(_opening(group[0].choice))}'
        f'</span><ul class="chain-ways">{ways}</ul></li>'
    )


def _step_links(steps: Sequence[Step], keys: Sequence[str], first: int = 1) -> str:
    """Every step as a link to where it stands on its own page, numbered from `first`.

    The number is written out rather than counted by the list, because a collapsed section's
    steps are not laid out and a counter would renumber the sections after it.
    """
    return "".join(
        f'<li><a href="#step-{key}"><span class="step-mark">{n}</span> '
        f"<span>{escape(step.title)}</span></a></li>"
        for n, (key, step) in enumerate(zip(keys, steps, strict=True), first)
    )


def _step_nav(steps: Sequence[Step], keys: Sequence[str]) -> str:
    """Every step as a link, under the section it belongs to wherever the steps name one.

    One group per run of steps sharing a `Step.section`, each a `<details>` carrying its own
    steps' keys: `protocol.js` counts the marks under it and opens the one the bench is in, so
    three like rounds read as three rounds and not as one flat list. Steps naming no section
    stay the one list they were.
    """
    if not any(step.section for step in steps):
        return f'<ul class="step-nav step-list">{_step_links(steps, keys)}</ul>'
    groups, at, opened = [], 0, False
    for section, run in groupby(steps, key=lambda step: step.section):
        size = len(tuple(run))
        under = _step_links(steps[at : at + size], keys[at : at + size], at + 1)
        if section:
            groups.append(
                f'<li class="step-group"><details data-steps='
                f'"{escape(" ".join(keys[at : at + size]))}"{"" if opened else " open"}>'
                f"<summary>{escape(section)} "
                f'<span class="section-progress muted">0 of {size} done</span></summary>'
                f'<ul class="step-list">{under}</ul></details></li>'
            )
            opened = True
        else:
            groups.append(under)
        at += size
    return f'<ul class="step-nav">{"".join(groups)}</ul>'


def _within(protocol: Protocol, keys: Sequence[str]) -> str:
    """Return the right column: every step of this page, so a reader jumps within it."""
    if not protocol.steps:
        return ""
    return (
        '<nav class="column within" aria-label="Steps">'
        f'<h2 class="column-title">This protocol</h2>'
        f"{_step_nav(protocol.steps, keys)}</nav>\n"
    )


def _place(folder: Folder, here: str) -> str:
    """Where this page stands in the run, named in words, because it prints with the page.

    The ways of one job take one place, so both are protocol 4 of 6; the line then names the
    sibling, says to do one and not both, and links what to weigh where the overview carries it.
    """
    groups = folder.places
    at = next((n for n, group in enumerate(groups) if any(p.href == here for p in group)), -1)
    if at < 0:
        return ""
    total = len(groups)
    job = groups[at][0].choice
    before = _place_before(groups[at - 1]) if at else "The run starts here."
    after = _place_after(groups[at + 1]) if at + 1 < total else "The run ends here."
    return (
        f'<p class="neighbours">Protocol {at + 1} of {total}{_sibling(groups[at], here)}. '
        f"{_guidance(folder, job)}{before} {after}</p>\n"
    )


def _sibling(group: tuple[Page, ...], here: str) -> str:
    """Return what this page adds to its place, where it is one of several ways of one job."""
    if not group[0].choice:
        return ""
    others = _or([_page_link(page) for page in group if page.href != here])
    alone = "not both" if len(group) == 2 else "not the others"
    return (
        f", and one of {_spelled(len(group))} ways to {escape(group[0].choice)} — "
        f"do this one or {others}, {alone}"
    )


def _guidance(folder: Folder, job: str) -> str:
    """Return the link to what to weigh before picking a way, or nothing where nothing says."""
    if not job or job not in folder.explained:
        return ""
    where = f"{escape(folder.index)}#topic-{escape(slug(job))}"
    return f'<a href="{where}">How to choose</a> is on the overview. '


def _place_before(group: tuple[Page, ...]) -> str:
    """Return what the bench did before this place, whichever way it did it."""
    if group[0].choice:
        job = escape(group[0].choice)
        return f"Comes after whichever of the {_spelled(len(group))} ways to {job} you did."
    return f"Comes after {_page_link(group[0])}."


def _place_after(group: tuple[Page, ...]) -> str:
    """Return what the bench does next, naming every way where the next place is a choice."""
    if group[0].choice:
        job = escape(group[0].choice)
        ways = _or([_page_link(page) for page in group])
        return f"Next is one of the {_spelled(len(group))} ways to {job}: {ways}."
    return f"Next is {_page_link(group[0])}."


def _page_link(page: Page) -> str:
    """One neighbouring page as a link, so the name a reader is given is the way there."""
    return f'<a href="{escape(page.href)}">{escape(page.title)}</a>'


def _asset(name: str) -> str:
    return files("mbio.protocol").joinpath(name).read_text(encoding="utf-8")


def _fonts(body: str) -> str:
    """Each face a drawing in `body` was measured in, written once, or nothing where it draws none.

    A drawing pins every line of text to the width the package measured it at, so a page drawing
    one as text carries those faces or the browser stretches its own to fit. Every figure holding
    a drawing is marked `drawing`, whatever it draws, so the faces follow it onto any page.
    """
    drawn = 'class="drawing' in body
    return "".join(font_face(font) for font in (SANS, BOLD, MONO)) if drawn else ""


def _duration(seconds: float) -> str:
    """One duration in words, to the second under an hour and to the minute from an hour up."""
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


def _bullets(items: Iterable[str], paths: Sequence[str] = ()) -> str:
    lines = "".join(f"<li>{_linked(item, paths)}</li>" for item in items)
    return f"<ul>{lines}</ul>" if lines else ""


def _linked(text: str, paths: Sequence[str] = ()) -> str:
    """Escape `text`, linking every file of `paths` it names by that file's own name.

    Protocol prose says a file the way the bench says it — `pool.tsv`, never a path — so the
    page links the name where it is read rather than repeating a path beside it. A name the
    text does not say links nothing, so a run declares every file it writes once.
    """
    named = {name: path for path in paths if (name := Path(path).name) and name not in {".", ".."}}
    if not named:
        return escape(text)
    # Longest first, so a name spelled inside a longer one never wins the alternation.
    pattern = "|".join(re.escape(name) for name in sorted(named, key=len, reverse=True))
    out: list[str] = []
    at = 0
    for match in re.finditer(rf"(?<![\w.-])({pattern})(?![\w-]|\.\w)", text):
        out.append(escape(text[at : match.start()]))
        href = escape(named[match[0]], quote=True)
        out.append(f'<a href="{href}">{escape(match[0])}</a>')
        at = match.end()
    out.append(escape(text[at:]))
    return "".join(out)


def _copy(text: str, label: str = "Copy") -> str:
    return f'<button type="button" class="copy" data-copy="{escape(text)}">{escape(label)}</button>'


def _header(protocol: Protocol, keys: Sequence[str], *, toc: bool = True, place: str = "") -> str:
    parts = [f'<header class="intro">\n<h1>{escape(protocol.title)}</h1>\n{place}']
    if protocol.summary:
        parts.append(f'<p class="summary">{_linked(protocol.summary, protocol.files)}</p>\n')
    if protocol.overview:
        facts = "".join(
            f"<div><dt>{escape(k)}</dt><dd>{_linked(v, protocol.files)}</dd></div>"
            for k, v in protocol.overview.items()
        )
        parts.append(f'<dl class="overview">{facts}</dl>\n')
    parts.append(_handover(protocol))
    if protocol.highlights:
        lines = "".join(f"<p>{_linked(one, protocol.files)}</p>" for one in protocol.highlights)
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
            parts.append(
                f'<nav class="toc" aria-label="Steps">{_step_nav(protocol.steps, keys)}</nav>\n'
            )
    parts.append("</header>\n")
    return "".join(parts)


def _handover(protocol: Protocol) -> str:
    """Return what the bench holds before this protocol and what it is left with, or nothing.

    A protocol read on its own states them for its reader; inside a run the index draws the same
    names as a chain, and nothing is written where the protocol declares neither. A rule either
    of them carries stands under both lists.
    """
    handed = (("Have in hand", protocol.consumes), ("Leaves you with", protocol.produces))
    blocks = [
        f"<div><h3>{heading}</h3>"
        f"<ul>{''.join(f'<li>{_item_named(item)}</li>' for item in items)}</ul></div>"
        for heading, items in handed
        if items
    ]
    # A rule travels with the item, so a page handed a thing it may break states the rule even
    # where no step of it spells that thing's name.
    carried = [(item.name, rule) for _, items in handed for item in items for rule in item.rules]
    return f'<div class="handover">{"".join(blocks)}</div>\n{_rules(carried)}' if blocks else ""


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
    paths: Sequence[str] = (),
    *,
    note: bool = True,
    sources: str = "",
) -> str:
    """Everything that is not an oligo, and the hardware as one light line under it.

    `used` names, row by row, which protocols of a run take each material, for the page a whole
    run shares. A protocol's own list leaves it empty. That page is what is ordered rather than
    what is laid out, so it drops the bench note and `note` is how. `paths` is what the page
    links a file name to, so a note saying which sheet a reagent is ordered from opens it.
    `sources` is the page a citation on this one resolves against, from `Folder.sources_at`.

    A material's rules and its cautions stand under the table, so they reach a reader of this
    list and not only the steps that pipette the tube. Two tubes carrying one sentence state it
    once.
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
            + "".join(f"<td>{_linked(get(material), paths)}</td>" for _, get in shown)
            + (f"<td>{_cite(material.citation, sources)}</td>" if cited else "")
            + (f"<td>{escape(used[i])}</td>" if used else "")
            + "</tr>"
            for i, material in enumerate(materials)
        )
        table = (
            f'<div class="scroll"><table><thead><tr>{head}</tr></thead>'
            f"<tbody>{rows}</tbody></table></div>"
        )
    carried = [(m.name, rule) for m in materials for rule in m.rules]
    said = dict.fromkeys(c for m in materials for c in m.cautions)
    cautions = _cautions([Caution(one) for one in said], paths)
    line = ""
    if equipment:
        line = (
            f'<p class="equipment"><strong>Equipment:</strong> {escape(", ".join(equipment))}</p>'
        )
    return (
        '<section class="block materials" id="materials">\n<h2>Materials</h2>\n'
        f"{table}{line}{_rules(carried, sources)}{cautions}\n</section>\n"
    )


def _oligos(protocol: Protocol) -> str:
    """Render the order sheet: one row each, or what the rows share where there are many.

    A sheet long enough that nobody reads it row by row says what its rows have in common
    instead: from `OLIGO_SUMMARY` up the section states the count, the lengths, the melting
    temperatures and what each purpose covers, and the sheet goes under a closed toggle. Rows
    are in the order a plate seats them, so a reader walking the plate walks the sheet; the
    file the sheet is ordered from keeps the order its plan wrote it in.
    """
    if not protocol.oligos:
        return ""
    seats = _seats(protocol.plates)
    unseated = _Seat(len(protocol.plates), 0, "", "")
    oligos = tuple(sorted(protocol.oligos, key=lambda one: seats.get(one.name, unseated)))
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
    sheet = (
        f'<div class="scroll"><table><thead><tr>{head}</tr></thead>'
        f"<tbody>{''.join(rows)}</tbody></table></div>"
    )
    shown = sheet
    if len(oligos) >= OLIGO_SUMMARY:
        shown = (
            f"{_oligo_summary(oligos, seats)}"
            f'<details class="listing"><summary>All {len(oligos)} rows</summary>{sheet}</details>'
        )
    # Which rows warn, and the copy of every sequence, stay outside the toggle: the longer the
    # sheet the more they are wanted, and summarising is not a reason to hide either.
    return (
        '<section class="block oligos" id="oligos">\n<h2>Oligos</h2>\n'
        f"{_order_sheet(protocol.order_sheet)}{shown}{_oligo_checks(oligos)}{copy_all}"
        "\n</section>\n"
    )


class _Seat(NamedTuple):
    """Where one name sits. Ordered by the plate's place in the protocol, then reading order.

    `place` counts the wells of its own plate in reading order, so two seats a plate apart
    subtract to how many wells lie between them.
    """

    order: int
    place: int
    plate: str
    well: str


def _seats(plates: tuple[Plate, ...]) -> dict[str, _Seat]:
    """Return where each seated name sits, by the first plate that seats it."""
    found: dict[str, _Seat] = {}
    for index, plate in enumerate(plates):
        for well, held in plate.seating.items():
            at = well_at(well)
            if at is not None:
                place = at[0] * plate.columns + at[1]
                found.setdefault(held, _Seat(index, place, plate.name, well))
    return found


def _order_sheet(path: str) -> str:
    """Return the file the sheet is ordered from: a page is not what a supplier is sent."""
    if not path:
        return ""
    return (
        f'<p class="order-sheet">Order sheet: '
        f'<a href="{escape(path, quote=True)}">{escape(path)}</a></p>'
    )


def _oligo_summary(oligos: tuple[Oligo, ...], seats: Mapping[str, _Seat]) -> str:
    """Return what the rows share, and one row a purpose: where it sits and how much of it."""
    lengths = [len(oligo.sequence) for oligo in oligos]
    melting = [oligo.tm_c for oligo in oligos if oligo.tm_c is not None]
    facts = " · ".join(
        text
        for text in (
            _count(len(oligos), "oligo"),
            f"{_span(min(lengths), max(lengths))} bases",
            f"Tm {_span(min(melting), max(melting), 1)} °C" if melting else "",
        )
        if text
    )
    groups: dict[str, list[Oligo]] = {}
    places: dict[str, list[int]] = {}
    for place, oligo in enumerate(oligos):
        groups.setdefault(oligo.purpose, []).append(oligo)
        places.setdefault(oligo.purpose, []).append(place)
    rows = "".join(
        f"<tr><td>{escape(purpose)}</td>"
        f"<td>{escape(_names(group, places[purpose]))}</td>"
        f'<td class="num">{len(group)}</td>'
        f"<td>{escape(_where([seats[one.name] for one in group if one.name in seats]))}</td></tr>"
        for purpose, group in groups.items()
    )
    return (
        f'<p class="muted">{escape(facts)}</p>'
        '<div class="scroll"><table><thead><tr><th>For</th><th>Names</th>'
        f'<th class="num">Oligos</th><th>Wells</th></tr></thead><tbody>{rows}</tbody>'
        "</table></div>"
    )


def _span(low: float, high: float, places: int = 0) -> str:
    """Return a range, or the one value where both ends are it."""
    first, last = f"{low:.{places}f}", f"{high:.{places}f}"
    return first if first == last else f"{first} to {last}"


def _names(group: Sequence[Oligo], places: Sequence[int]) -> str:
    """Return a group's names, as a range only where it holds every name between its ends.

    `places` is where each of the group's rows sits in the sheet. A group the sheet interleaves
    with another is a pick and not a run, so it says so rather than naming two ends to read
    between.
    """
    if len(group) == 1:
        return group[0].name
    ends = f"{group[0].name} to {group[-1].name}"
    whole = places[-1] - places[0] + 1 == len(group)
    return ends if whole else f"some of {ends}"


def _where(seats: Sequence[_Seat]) -> str:
    """Return the wells a group occupies, as a range only where it fills one.

    A group seated among others fills no range, so it gives its count inside the two wells it
    reaches: someone counting that many wells off the first would otherwise stop short.
    """
    if not seats:
        return ""
    first, last = seats[0], seats[-1]
    if first == last:
        return f"{first.plate} {first.well}"
    if first.plate != last.plate:
        return f"{first.plate} {first.well} to {last.plate} {last.well}"
    if last.place - first.place + 1 == len(seats):
        return f"{first.plate} {first.well} to {last.well}"
    return f"{first.plate}, {len(seats)} wells from {first.well} to {last.well}"


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


def _cite(citation: Citation | None, sources: str = "") -> str:
    """One citation, as the page shows it beside the number it carries.

    `sources` is the page the sources list stands on, as `Folder.sources_at` gives it, and is
    empty where that is the page being rendered.
    """
    if citation is None:
        return ""
    where = f" {citation.locator}" if citation.locator else ""
    return (
        f'<a class="cite" href="{escape(sources)}#source-{escape(slug(citation.source))}">'
        f"{escape(citation.source + where)}</a>"
    )


def _after(citation: Citation | None, sources: str = "") -> str:
    """`_cite`, set off from the text before it; nothing when uncited."""
    return f" {_cite(citation, sources)}" if citation else ""


def _rules(rules: Iterable[tuple[str, Rule]], sources: str = "") -> str:
    """Every rule the things in this step carry, computed from the carrier, never stored.

    A rule hangs on the material, the plate or the handed item it governs, so it shows
    wherever that is and no edit to a step's prose can drop it.
    """
    items = "".join(
        f'<li class="rule is-{rule.kind}"><strong>{escape(carrier)}: '
        f"{escape('never' if rule.kind == 'forbids' else 'always')} "
        f"{escape(rule.subject)}</strong> {escape(rule.detail)}{_after(rule.citation, sources)}</li>"
        for carrier, rule in rules
    )
    return f'<ul class="rules" aria-label="Rules">{items}</ul>\n' if items else ""


def _cautions(cautions: Iterable[Caution], paths: Sequence[str] = ()) -> str:
    """Every caution, as the one paragraph both the steps and the reagents page show it in.

    One site renders it, so a sentence cannot read two ways on two pages, and a cited one
    anchors into the Sources section where a troubleshooting row does.
    """
    return "".join(
        f'<p class="caution"><strong>Caution:</strong> {_linked(one.text, paths)}'
        f"{_after(one.citation)}</p>\n"
        for one in cautions
    )


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

    Four statements, each ended: what is missing, what it waits on, what would fill it and where
    it stands. Run together they are re-parsed halfway through.

    `found` is markup naming where the hole stands, for a page that is not the one holding it.
    """
    return (
        f'<li class="hole" id="hole-{escape(hole.id)}">{_missing(hole)} '
        f'<span class="hole-kind">{escape(HOLE_KINDS[hole.kind])}</span>'
        f"{_filled_by(hole)}{found}</li>"
    )


def _missing(hole: Hole) -> str:
    """Return the id, the mark a number can never be read from, and what is not known."""
    where = f"{escape(hole.where)}: " if hole.where else ""
    return (
        f'<span class="hole-id">{escape(hole.id)}</span> '
        f'<span class="hole-none">{NO_NUMBER}</span> — {where}{_ended(escape(hole.missing))}'
    )


def _filled_by(hole: Hole) -> str:
    """Return what would close the hole, or nothing where it names none."""
    return f" <em>Filled by {escape(hole.filled_by)}.</em>" if hole.filled_by else ""


def _ended(text: str) -> str:
    """`text` with a full stop, so one statement cannot run into the next."""
    return text if text.endswith((".", "!", "?")) else f"{text}."


def _holes(protocol: Protocol) -> str:
    """Every hole the protocol carries, collected under its stable id."""
    holes = protocol.all_holes
    if not holes:
        return ""
    items = "".join(_hole(hole) for hole in holes)
    return (
        '<section class="block holes" id="holes">\n<h2>Holes</h2>\n'
        f"<p>{_count(len(holes), 'number')} this protocol would otherwise have to invent. "
        f"{HOLES_INTRO}</p>\n<ul>{items}</ul>\n</section>\n"
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
        # A plate the protocol makes carries its handling rules here, where a reagent carries
        # them in the materials table.
        + _rules((one.name, rule) for rule in one.rules)
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


def _transfer(transfer: Transfer, plates: tuple[Plate, ...]) -> str:
    """Return a transfer drawn as the pattern it repeats, or as a table where it repeats none.

    A stamp is the two plates with every well it touches filled, since the pattern is what the
    bench follows and the moves spell out one thing 96 times; those go under a closed toggle.
    A transfer that is no stamp, that names a plate `plates` does not declare, or that moves
    through a well the declared format has not got, is the table it was: nothing says what its
    wells look like, and a drawing missing the wells it cannot place would say it wrongly.
    """
    stamp = transfer.stamp
    declared = {plate.name: plate for plate in plates}
    source = declared.get(transfer.plates[0])
    destination = declared.get(transfer.plates[-1])
    taken = {move.source.well for move in transfer.moves}
    filled = {move.destination.well for move in transfer.moves}
    if stamp is None or source is None or destination is None:
        return _transfer_table(transfer)
    if not (taken <= set(source.well_names) and filled <= set(destination.well_names)):
        return _transfer_table(transfer)
    drawn = "".join(
        _stamped(f"{word} {one.name}", one, wells, transfer.title)
        for word, one, wells in (("From", source, taken), ("Into", destination, filled))
    )
    meta = _transfer_meta(transfer, stamp)
    return (
        f'<figure class="drawing transfer rows"><figcaption>{escape(transfer.title)} '
        f'<span class="muted">{escape(meta)}</span>{_after(transfer.citation)}</figcaption>'
        f'{drawn}<details class="listing"><summary>{_count(len(transfer.moves), "move")}'
        f"</summary>{_moves(transfer)}</details></figure>\n"
    )


def _transfer_table(transfer: Transfer) -> str:
    """Return a transfer as the table it has always been: where each thing goes, a row each."""
    return (
        f'<figure class="transfer"><figcaption>{escape(transfer.title)} '
        f'<span class="muted">{escape(_transfer_meta(transfer))}</span>'
        f"{_after(transfer.citation)}</figcaption>{_moves(transfer)}</figure>\n"
    )


def _transfer_meta(transfer: Transfer, stamp: Stamp | None = None) -> str:
    """Return the caption beside the title: who moves how much, and the pattern if any."""
    volume = f"{number(transfer.moves[0].volume_ul)} µL each" if stamp else ""
    return " · ".join(
        text
        for text in (
            transfer.instrument,
            _count(len(transfer.moves), "move") if stamp else f"{len(transfer.moves)} wells",
            volume or " → ".join(transfer.plates),
            stamp.words if stamp else "",
            transfer.note,
        )
        if text
    )


def _stamped(name: str, plate: Plate, wells: set[str], holds: str) -> str:
    """One plate of a stamp, every well the transfer touches filled and the rest left empty."""
    drawn = draw_plate(
        name, plate.rows, plate.columns, plate.row_labels, seating=dict.fromkeys(wells, holds)
    )
    return f'<div class="row">{drawn.element()}</div>'


def _moves(transfer: Transfer) -> str:
    """Every move as a row: from, to, and how much."""
    rows = "".join(
        f"<tr><td>{escape(move.source.plate)} {escape(move.source.well)}</td>"
        f"<td>{escape(move.destination.plate)} {escape(move.destination.well)}</td>"
        f'<td class="num">{number(move.volume_ul)}</td></tr>'
        for move in transfer.moves
    )
    return (
        '<div class="scroll"><table><thead><tr><th>From</th><th>To</th>'
        f'<th class="num">µL</th></tr></thead><tbody>{rows}</tbody></table></div>'
    )


def _bill(bill: Bill | None, sources: str = "") -> str:
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
            else escape(row.charge) + _after(row.citation, sources)
        )
        quantity = (
            f"{number(row.quantity)} {escape(row.unit)}"
            if row.quantity is not None
            else escape(row.unit) or f'<span class="hole-none">{NO_NUMBER}</span>'
        )
        cells = [
            f"<td>{escape(row.item)}</td>",
            f"<td>{escape(row.key)}</td>",
            f'<td class="num">{quantity}</td>',
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
        "<p>Every item this run consumes has a row. A quantity comes from the design and is "
        "here whatever is loaded; one the design does not size carries what the protocol "
        "states. A charge comes only from a price record; where none prices a row, the money "
        "is a hole and no figure is estimated. Price steers no part of the design.</p>\n"
        f'<div class="scroll"><table><thead><tr>{head}</tr></thead>'
        f"<tbody>{''.join(rows)}</tbody>{total}</table></div>\n</section>\n"
    )


def _sources(sources: Mapping[str, Source], cited: Mapping[str, str] | None = None) -> str:
    """Every document a number was read from, so a citation resolves on the page itself.

    `cited` names, key by key, what cites each on the page a run shares. A key it leaves out is
    listed with no citer, which is what a source nothing on the run cites has.
    """
    if not sources:
        return ""
    items = "".join(
        f'<li id="source-{escape(slug(key))}"><strong>{escape(key)}</strong> '
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
        + (
            f' <span class="cited-by">cited by {escape(cited[key])}</span>'
            if cited and key in cited
            else ""
        )
        + "</li>"
        for key, source in sources.items()
    )
    return (
        '<section class="block sources" id="sources">\n<h2>Sources</h2>\n'
        f"<ul>{items}</ul>\n</section>\n"
    )


def _step(n: int, step: Step, key: str, protocol: Protocol, base: Path, section: str = "") -> str:
    """One step, addressed by its own key, so a reworded title keeps the bench's tick.

    Everything inside it is marked by that anchor, a dot and what it is. A key is a slug and
    holds no dot, so no mark here can spell another step's anchor or another of its marks.
    """
    anchor = f"step-{key}"
    parts = [f'<p class="step-section">{escape(section)}</p>\n'] if section else []
    parts += [
        f'<section class="step" id="{anchor}">\n<h2 class="step-title"><label>'
        f'<input type="checkbox" class="done" data-key="{anchor}">'
        f'<span class="step-n">{n}</span><span>{escape(step.title)}</span></label></h2>\n'
    ]
    parts.append(_rules(protocol.rules_for(step)))
    parts.append(_cautions(protocol.cautions_for(step), protocol.files))
    if step.instructions:
        items = "".join(
            f'<li><label><input type="checkbox" data-key="{anchor}.{i}">'
            f"<span>{_linked(text, protocol.files)}</span></label></li>"
            for i, text in enumerate(step.instructions, 1)
        )
        parts.append(f'<ol class="instructions">{items}</ol>\n')
    parts += [_figure(f, base, f"step {n} {step.title!r}") for f in step.figures]
    parts += [_table(f"{anchor}.table.{i}", t) for i, t in enumerate(step.tables, 1)]
    parts += [_program(p) for p in step.programs]
    parts += [_transfer(t, protocol.plates) for t in step.transfers]
    if step.holes:
        items = "".join(_hole(hole) for hole in step.holes)
        parts.append(f'<ul class="holes-here" aria-label="Holes">{items}</ul>\n')
    if step.timers:
        timers = "".join(_timer(f"{anchor}.timer.{i}", t) for i, t in enumerate(step.timers, 1))
        parts.append(f'<div class="timers">{timers}</div>\n')
    parts.append(_waits(step.waits))
    if step.expected or step.gels:
        gels = "".join(_gel(g) for g in step.gels)
        parts.append(
            '<div class="expected"><h3>Expected result</h3>'
            f"{_bullets(step.expected, protocol.files)}{gels}</div>\n"
        )
    if step.troubleshooting:
        entries = "".join(
            f"<dt>{escape(t.problem)}</dt>"
            f"<dd>{_linked(t.solution, protocol.files)}{_after(t.citation)}</dd>"
            for t in step.troubleshooting
        )
        parts.append(f'<div class="trouble"><h3>Troubleshooting</h3><dl>{entries}</dl></div>\n')
    if step.notes:
        items = "".join(
            f"<li>{_linked(note.text, protocol.files)}{_after(note.citation)}</li>"
            for note in step.noted
        )
        parts.append(f'<div class="notes"><h3>Notes</h3><ul>{items}</ul></div>\n')
    parts.append("</section>\n")
    return "".join(parts)


def _waits(waits: tuple[Wait, ...], sources: str = "") -> str:
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
        + f"{_after(wait.citation, sources)}</li>"
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


def _temperature(step: Incubation, cycles: int | None) -> str:
    """One incubation's temperature cell: where it starts, and where `cycles` of it end.

    A stepping incubation's end is read off `last_c`, never stored; where the count is a hole,
    the step a cycle stands in for an end nothing bounds.
    """
    start = number(step.temperature_c)
    if step.delta_c is None:
        return f"{start} °C"
    last = step.temperature_c if cycles is None else step.last_c(cycles)
    end = "" if last == step.temperature_c else f" → {number(last)}"
    return f'{start}{end} °C<br><span class="muted">{number(step.delta_c)} °C a cycle</span>'


def _program(program: ThermocyclerProgram) -> str:
    meta = []
    if program.lid_temperature_c is not None:
        meta.append(f"lid {number(program.lid_temperature_c)} °C")
    steps = [step for stage in program.stages for step in stage.incubations]
    if program.duration_seconds is not None:
        # A ramp is the block changing temperature, so a program held at one has none to add.
        held = len({step.temperature_c for step in steps}) == 1 and not any(
            step.delta_c for step in steps
        )
        meta.append(_duration(program.duration_seconds) + ("" if held else " plus ramps"))
    title = escape(program.title or "Thermocycler program")
    caption = f' <span class="muted">· {escape(" · ".join(meta))}</span>' if meta else ""
    cited = {step.citation for step in steps}
    # One source behind every incubation is the program's, so it is cited once above the table.
    shared = cited.pop() if len(cited) == 1 else None
    bodies = []
    for stage in program.stages:
        # `×` marks the stage that repeats, so a count of 10 is never read as a tenth cycle.
        count = (
            f'<span class="hole-none">{NO_NUMBER}</span>'
            if stage.cycles is None
            else f"×{stage.cycles}"
            if stage.cycles > 1
            else str(stage.cycles)
        ) + _after(stage.citation)
        rows = []
        for i, step in enumerate(stage.incubations):
            time = "∞" if step.seconds is None else _duration(step.seconds)
            cycles = (
                f'<td class="num" rowspan="{len(stage.incubations)}">{count}</td>' if i == 0 else ""
            )
            rows.append(
                f"<tr><td>{escape(step.label)}{'' if shared else _after(step.citation)}</td>"
                f'<td class="num">{_temperature(step, stage.cycles)}</td>'
                f'<td class="num">{time}</td>{cycles}</tr>'
            )
        bodies.append(f'<tbody class="stage">{"".join(rows)}</tbody>')
    return (
        f'<figure class="program"><figcaption>{title}{caption}{_after(shared)}</figcaption>'
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

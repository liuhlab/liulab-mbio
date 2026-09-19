"""An HTML page parsed into a tree a test can search, as a browser would read the markup.

Each inline SVG that is well-formed XML, as `svg.document` writes one, is read by the XML parser,
which is quicker; it gives the same tree, its names lowercased as HTML's are.
"""

import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from html.parser import HTMLParser
from xml.parsers import expat

VOID = {"area", "base", "br", "col", "embed", "hr", "img", "input", "link", "meta", "source", "wbr"}


@dataclass
class Node:
    tag: str
    attrs: dict[str, str]
    children: list["Node | str"] = field(default_factory=list)

    def iter(self) -> Iterator["Node"]:
        # Each node before its children, without a generator per level.
        stack = [iter(self.children)]
        while stack:
            for child in stack[-1]:
                if isinstance(child, Node):
                    yield child
                    stack.append(iter(child.children))
                    break
            else:
                stack.pop()

    def find_all(
        self, tag: str | tuple[str, ...] = (), cls: str = "", **attrs: str
    ) -> list["Node"]:
        tags = (tag,) if isinstance(tag, str) else tag
        return [
            n
            for n in self.iter()
            if (not tags or n.tag in tags)
            and (not cls or cls in n.attrs.get("class", "").split())
            and all(n.attrs.get(k.replace("_", "-")) == v for k, v in attrs.items())
        ]

    @property
    def text(self) -> str:
        return "".join(c if isinstance(c, str) else c.text for c in self.children).strip()


class _Builder(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.root = Node("#document", {})
        self.stack = [self.root]

    def handle_starttag(self, tag, attrs) -> None:
        node = Node(tag, {k: v or "" for k, v in attrs})
        self.stack[-1].children.append(node)
        if tag not in VOID:
            self.stack.append(node)

    def handle_startendtag(self, tag, attrs) -> None:
        self.stack[-1].children.append(Node(tag, {k: v or "" for k, v in attrs}))

    def handle_endtag(self, tag) -> None:
        for i in range(len(self.stack) - 1, 0, -1):
            if self.stack[i].tag == tag:
                del self.stack[i:]
                return

    def handle_data(self, data) -> None:
        self.stack[-1].children.append(data)


_SVG = re.compile(r"<svg\b")


def _xml(stack: list[Node], image: str) -> bool:
    """Add an SVG read as XML under the node at the top of `stack`, or nothing if it is not XML."""
    parent, depth = stack[-1], len(stack)
    had = len(parent.children)
    xml = expat.ParserCreate()
    xml.buffer_text = True

    def start(tag: str, attrs: dict[str, str]) -> None:
        node = Node(tag.lower(), {key.lower(): value for key, value in attrs.items()})
        stack[-1].children.append(node)
        stack.append(node)

    def data(text: str) -> None:
        children = stack[-1].children
        if children and isinstance(children[-1], str):
            children[-1] += text
        else:
            children.append(text)

    xml.StartElementHandler = start
    xml.EndElementHandler = lambda _: stack.pop()
    xml.CharacterDataHandler = data
    try:
        xml.Parse(image, True)
    except expat.ExpatError:
        del stack[depth:]
        del parent.children[had:]
        return False
    return True


def parse(html: str) -> Node:
    builder = _Builder()
    at = 0
    while (found := _SVG.search(html, at)) and (end := html.find("</svg>", found.start())) >= 0:
        builder.feed(html[at : found.start()])
        at = end + len("</svg>")
        image = html[found.start() : at]
        if builder.rawdata or not _xml(builder.stack, image):
            builder.feed(image)
    builder.feed(html[at:])
    return builder.root

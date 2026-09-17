"""An HTML page parsed into a tree a test can search, as a browser would read the markup."""

from collections.abc import Iterator
from dataclasses import dataclass, field
from html.parser import HTMLParser

VOID = {"area", "base", "br", "col", "embed", "hr", "img", "input", "link", "meta", "source", "wbr"}


@dataclass
class Node:
    tag: str
    attrs: dict[str, str]
    children: list["Node | str"] = field(default_factory=list)

    def iter(self) -> Iterator["Node"]:
        for child in self.children:
            if isinstance(child, Node):
                yield child
                yield from child.iter()

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


def parse(html: str) -> Node:
    builder = _Builder()
    builder.feed(html)
    return builder.root

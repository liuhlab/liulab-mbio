"""The command line page lists exactly the verbs the two typer apps mount.

A verb is read off `docs/reference/cli.md` from every line of a fenced ``bash`` block starting
``pixi run``: drop those two words, keep the tokens up to the first one beginning with ``-``,
and drop the capitalised placeholders among them. The page is written so that rule finds every
verb it documents, so a verb added or renamed in the app fails here until the page follows.
"""

from pathlib import Path

import typer

from mbio.cli import app as mbio_app
from synbio.cli import app as synbio_app

PAGE = Path(__file__).resolve().parents[1] / "docs" / "reference" / "cli.md"


def _mounted(app: typer.Typer, prefix: tuple[str, ...]) -> set[str]:
    """Every verb `app` mounts, read off typer itself rather than off anything it renders."""
    found = {
        " ".join((*prefix, one.name or _spelled(one.callback))) for one in app.registered_commands
    }
    for group in app.registered_groups:
        inner = group.typer_instance
        if inner is not None:
            found |= _mounted(inner, (*prefix, group.name or inner.info.name or ""))
    return found


def _spelled(callback: object) -> str:
    """What typer calls a command whose name it takes from its callback."""
    return getattr(callback, "__name__", "").replace("_", "-")


def _listed(text: str) -> set[str]:
    """Every verb the page shows, by the rule this module's docstring states."""
    found: set[str] = set()
    shell = False
    for line in text.splitlines():
        if line.startswith("```"):
            shell = line.strip() == "```bash"
        elif shell and line.split()[:2] == ["pixi", "run"]:
            words = line.split()[2:]
            said = words[: next((at for at, one in enumerate(words) if one.startswith("-")), None)]
            found.add(" ".join(one for one in said if not one.isupper()))
    return found


def test_the_command_line_page_lists_every_verb_the_apps_mount() -> None:
    mounted = _mounted(mbio_app, ("mbio",)) | _mounted(synbio_app, ("synbio",))
    listed = _listed(PAGE.read_text(encoding="utf-8"))
    assert listed == mounted, (
        f"on the page and not in the app: {sorted(listed - mounted)}; "
        f"in the app and not on the page: {sorted(mounted - listed)}"
    )

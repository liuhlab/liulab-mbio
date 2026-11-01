"""What the verification tests share: what the page hands molviz to draw."""

from typing import Any

import molviz
import pytest

#: What each `molviz.draw` call was handed: the record, then the options.
type Handed = list[tuple[dict[str, Any], dict[str, Any]]]


@pytest.fixture
def handed(monkeypatch: pytest.MonkeyPatch) -> Handed:
    """Record what each `molviz.draw` call is handed, and draw it all the same."""
    calls: Handed = []
    draw = molviz.draw

    def spy(record: dict[str, Any], **options: Any) -> molviz.Drawing:
        calls.append((record, options))
        return draw(record, **options)

    monkeypatch.setattr(molviz, "draw", spy)
    return calls

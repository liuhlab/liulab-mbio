"""A chain of protocols read back as one, so a test can assert over a whole run at once.

A project splits a run across pages, and most of what a test asks -- does any step say this,
is this reagent bought, does the bill carry this row -- is about the run and not about which
page it landed on. `whole` puts the chain back together in the order it runs.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from liulab_mbio.protocol.model import Project, Protocol


def whole(project: Project) -> Protocol:
    """Return every protocol of `project` as one, each list in the order the chain runs.

    A reagent two protocols both buy is kept once, where it is first bought, as a reader
    working through the pages meets it. The run's own bill comes with the run's own sources, so
    a row citing the price record still resolves here.
    """
    from liulab_mbio.protocol.model import Protocol

    return Protocol(
        project.title,
        summary=project.summary,
        overview={
            label: value for one in project.protocols for label, value in one.overview.items()
        },
        highlights=tuple(line for topic in project.background for line in topic.body),
        checks=project.checks,
        consumes=project.inputs,
        materials=_once(one for p in project.protocols for one in p.materials),
        oligos=_once(one for p in project.protocols for one in p.oligos),
        equipment=_once(one for p in project.protocols for one in p.equipment),
        vessels=_once(one for p in project.protocols for one in p.vessels),
        plates=_once(one for p in project.protocols for one in p.plates),
        steps=tuple(one for p in project.protocols for one in p.steps),
        references=_once(one for p in project.protocols for one in p.references),
        sources=dict(project.sources)
        | {key: value for p in project.protocols for key, value in p.sources.items()},
        holes=_once(one for p in project.protocols for one in p.holes),
        bill=project.bill,
    )


def _once[T](found: Iterable[T]) -> tuple[T, ...]:
    """Return what `found` yields, the first of any repeat kept and the rest dropped."""
    return tuple(dict.fromkeys(found))

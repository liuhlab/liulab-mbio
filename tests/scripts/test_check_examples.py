"""What the example check refuses to call fresh, and that its table earns that answer.

Running the generators is the `docs` CI job's work, not the gate's. These are the rule's own
side: a generator reading another's committed output cannot be judged until this run has shown
that output fresh, and every such input is declared and written before it is read.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from types import ModuleType

REPO = Path(__file__).resolve().parents[2]


def _script() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "check_examples", REPO / "scripts/check_examples.py"
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # Registered before it runs: `dataclasses` resolves annotations through `sys.modules`.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_a_chained_generator_waits_for_its_input() -> None:
    """The maps cannot be judged until the plan that writes the record they draw has been."""
    script = _script()
    maps = next(one for one in script.GENERATORS if one.what == "the pUC19-GFP maps")
    record = "docs/examples/pUC19-GFP/product.dna"

    assert script.unproven(maps, frozenset()) == [record]
    assert script.unproven(maps, frozenset({record})) == []


def test_an_input_nothing_generates_is_never_waited_for() -> None:
    """A committed input is as fresh as it will ever be, so reading one holds nothing up."""
    script = _script()
    plan = next(one for one in script.GENERATORS if one.what == "the AP-1 library plan")
    written = set(script.MADE_HERE)

    assert "docs/examples/ap1-library/parts.fasta" in plan.reads
    assert "docs/examples/ap1-library/parts.fasta" not in written
    assert script.unproven(plan, frozenset(written)) == []


def test_every_generated_input_is_written_before_it_is_read() -> None:
    """Order is what earns the judgement, so no generator reads what a later one writes."""
    script = _script()
    fresh: set[str] = set()
    for generator in script.GENERATORS:
        assert script.unproven(generator, frozenset(fresh)) == [], generator.what
        fresh.update(script.owned(generator))

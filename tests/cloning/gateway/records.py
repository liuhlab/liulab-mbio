"""The worked-example Gateway records, loaded from the script that writes them.

`scripts/build_gateway_records.py` builds the donor and the destination the published example
is planned against. The tests read the same builders from there, so a change to either record
moves the tests and the committed example together instead of one of them alone.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[3] / "scripts/build_gateway_records.py"
_spec = importlib.util.spec_from_file_location("build_gateway_records", _SCRIPT)
assert _spec is not None
assert _spec.loader is not None
_module = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _module
_spec.loader.exec_module(_module)

ATTP = _module.ATTP
att_site = _module.att_site
attb_insert = _module.attb_insert
destination_vector = _module.destination_vector
donor_vector = _module.donor_vector
entry_clone = _module.entry_clone

__all__ = [
    "ATTP",
    "att_site",
    "attb_insert",
    "destination_vector",
    "donor_vector",
    "entry_clone",
]

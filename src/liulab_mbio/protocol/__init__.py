"""Bench protocols: a small model, loaded from JSON and rendered to one self-contained HTML file.

JSON keys are the field names of the classes below, lists stand for tuples, and only the
fields without a default are required::

    {"title": str, "summary": str, "overview": {label: value},
     "materials": [{"name", "sequence", "source", "storage", "note"}],
     "steps": [{"title", "instructions": [str], "cautions": [str], "notes": [str],
         "tables": [{"title", "reactions", "overage",
             "components": [{"name", "volume_ul", "stock", "final", "master_mix"}]}],
         "programs": [{"title", "lid_temperature_c",
             "stages": [{"cycles", "incubations": [{"label", "temperature_c", "seconds"}]}]}],
         "timers": [{"label", "seconds"}],
         "gels": [{"title", "ladder": {"name", "bands_bp"}, "lanes": [{"label", "bands_bp"}]}],
         "expected": [str], "troubleshooting": [{"problem", "solution"}]}],
     "references": [{"text", "url"}]}

An incubation's ``"seconds": null`` holds indefinitely.
"""

from liulab_mbio.protocol.model import (
    Component,
    Gel,
    Incubation,
    Ladder,
    Lane,
    Material,
    Protocol,
    ReactionTable,
    Reference,
    Stage,
    Step,
    ThermocyclerProgram,
    Timer,
    Troubleshooting,
    read_protocol,
)
from liulab_mbio.protocol.render import render_html, write_html

__all__ = [
    "Component",
    "Gel",
    "Incubation",
    "Ladder",
    "Lane",
    "Material",
    "Protocol",
    "ReactionTable",
    "Reference",
    "Stage",
    "Step",
    "ThermocyclerProgram",
    "Timer",
    "Troubleshooting",
    "read_protocol",
    "render_html",
    "write_html",
]

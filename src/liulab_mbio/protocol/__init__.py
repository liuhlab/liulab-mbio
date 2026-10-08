"""Bench protocols: a small model, read from and written to JSON, and rendered to one HTML file.

JSON keys are the field names of the classes below, lists stand for tuples, and only the
fields without a default are required; `write_protocol` writes every field::

    {"title": str, "summary": str, "overview": {label: short value},
     "highlights": [sentence], "checks": [{"name", "status", "detail"}],
     "consumes": [{"name", "what", "spec": [str], "storage"}],
     "produces": [{"name", "what", "spec": [str], "storage"}],
     "materials": [{"name", "supplier", "catalog", "storage", "amount", "note",
         "contains": [str], "citation": {"source", "locator"},
         "rules": [{"kind", "subject", "detail", "when", "citation"}]}],
     "oligos": [{"name", "sequence", "purpose", "tm_c", "stock", "note", "status",
         "checks": [{"name", "status", "detail"}]}],
     "equipment": [str],
     "vessels": [{"name", "kind", "catalog", "holds", "note"}],
     "plates": [{"name", "wells", "catalog", "holds", "seating": {well: name}, "note"}],
     "steps": [{"title", "section", "instructions": [str], "cautions": [str], "notes": [str],
         "tables": [{"title", "reactions", "overage",
             "components": [{"name", "volume_ul", "stock", "final", "master_mix"}]}],
         "programs": [{"title", "lid_temperature_c",
             "stages": [{"cycles", "citation",
                 "incubations": [{"label", "temperature_c", "seconds"}]}]}],
         "timers": [{"label", "seconds"}],
         "transfers": [{"title", "instrument", "note", "citation",
             "moves": [{"volume_ul", "source": {"plate", "well"},
                 "destination": {"plate", "well"}}]}],
         "holes": [{"id", "missing", "kind", "where", "filled_by", "issue"}],
         "gels": [{"title", "ladder": {"name", "bands_bp"}, "lanes": [{"label", "bands_bp"}]}],
         "expected": [str], "troubleshooting": [{"problem", "solution"}]}],
     "references": [{"text", "url"}],
     "sources": {key: {"document", "edition", "url", "read_as", "date"}},
     "holes": [{"id", "missing", "kind", "where", "filled_by", "issue"}],
     "bill": {"title", "currency", "total", "record",
         "rows": [{"item", "quantity", "unit", "key", "charge", "headroom", "citation",
             "hole"}]}}

An incubation's ``"seconds": null`` holds indefinitely, and a stage's ``"cycles": null`` leaves
the count blank where nothing sources it. A check's ``"status"`` is ``"pass"``,
``"warn"`` or ``"fail"``; an oligo's may also be absent, which says nothing judged that row. An
``"overview"`` value is a card: a few words, never a sentence.

A protocol declares what it consumes and what it produces, and nothing else about its place in a
run. A project chains protocols by those names, `write_project` writing one file of them::

    {"title": str, "summary": str,
     "inputs": [{"name", "what", "spec": [str], "storage"}],
     "protocols": [a protocol, as above],
     "checks": [{"name", "status", "detail"}], "bill": a bill, as above}

A number's provenance is its row's ``"citation"``, whose ``"source"`` keys ``"sources"``. A
number nobody published is a ``"hole"``: the field it belongs to stays empty and the hole stands
beside it, so a loader never reads a union and a reader never sees a guess. A rule hangs on the
material it belongs to, so it follows the material into every step that uses it.
"""

from liulab_mbio.protocol.model import (
    FORMATS,
    OVERVIEW_CHARS,
    Bill,
    BillRow,
    Check,
    Citation,
    Component,
    Gel,
    Hole,
    Incubation,
    Item,
    Ladder,
    Lane,
    Material,
    Move,
    Oligo,
    Plate,
    Project,
    Protocol,
    ReactionTable,
    Reference,
    Rule,
    Source,
    Stage,
    Step,
    ThermocyclerProgram,
    Timer,
    Transfer,
    Troubleshooting,
    Vessel,
    Well,
    citing,
    number,
    read_project,
    read_protocol,
    row_label,
    write_project,
    write_protocol,
)
from liulab_mbio.protocol.render import render_html, write_html

__all__ = [
    "FORMATS",
    "OVERVIEW_CHARS",
    "Bill",
    "BillRow",
    "Check",
    "Citation",
    "Component",
    "Gel",
    "Hole",
    "Incubation",
    "Item",
    "Ladder",
    "Lane",
    "Material",
    "Move",
    "Oligo",
    "Plate",
    "Project",
    "Protocol",
    "ReactionTable",
    "Reference",
    "Rule",
    "Source",
    "Stage",
    "Step",
    "ThermocyclerProgram",
    "Timer",
    "Transfer",
    "Troubleshooting",
    "Vessel",
    "Well",
    "citing",
    "number",
    "read_project",
    "read_protocol",
    "render_html",
    "row_label",
    "write_html",
    "write_project",
    "write_protocol",
]

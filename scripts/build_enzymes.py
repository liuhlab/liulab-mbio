#!/usr/bin/env python3
"""Rebuild the shipped enzyme data from REBASE and the curated supplier table.

    pixi run python scripts/build_enzymes.py [--withrefm PATH] [--out PATH]

REBASE's `withrefm` file gives the recognition site, the cut offsets, who sells each enzyme,
and which enzymes share a specificity. Everything a supplier states about its own product —
commercial name, incubation and heat-inactivation temperature, methylation sensitivity — is in
`enzyme_properties.toml` beside this script, with a cited source per value: NEB's terms of use
forbid redistributing their pages, so those facts are read, cited and re-entered, never mirrored.

The enzymes in the data file are the ones `enzyme_properties.toml` names.

Without `--withrefm` the REBASE file is downloaded. Tests read excerpts instead and never
reach the network.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import tomllib
import urllib.request
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
PROPERTIES = Path(__file__).resolve().parent / "enzyme_properties.toml"
DATA = REPO / "src/liulab_mbio/data/enzymes.json"
WITHREFM_URL = "https://rebase.neb.com/rebase/link_withrefm"

#: Which supplier fields the curated table carries, and so which need a source named.
SUPPLIER_FIELDS = (
    "commercial_name",
    "catalog_number",
    "supplier",
    "incubation_celsius",
    "heat_inactivation_celsius",
    "heat_inactivation_minutes",
    "methylation",
)

_FIELD_RE = re.compile(r"^<(\d)>(.*)$", re.MULTILINE)
_OFFSETS_RE = re.compile(r"^([A-Z]+)\((-?\d+)/(-?\d+)\)$")


@dataclass(frozen=True)
class RebaseEntry:
    """One REBASE record, reduced to the fields the data file is built from."""

    name: str
    site: str
    cuts: tuple[int, int]
    same_specificity: tuple[str, ...]
    suppliers: str


def parse_site(notation: str) -> tuple[str, int, int]:
    """Read a REBASE recognition sequence into a site and its two cut offsets.

    REBASE marks a cut inside the site with ``^`` and one outside it with offsets in
    parentheses, counted from the 3' end of the site. Both become offsets from the first base
    of the site, as `liulab_mbio.enzymes` defines them.

    Raises
    ------
    ValueError
        If REBASE gives no cut position, or gives two pairs of them, as it does for the few
        enzymes that cut on both sides of their site.

    Examples
    --------
    >>> parse_site("GGTCTC(1/5)")
    ('GGTCTC', 7, 11)
    >>> parse_site("CCC^GGG")
    ('CCCGGG', 3, 3)
    """
    notation = notation.strip()
    if offsets := _OFFSETS_RE.match(notation):
        site = offsets.group(1)
        top, bottom = (len(site) + int(offsets.group(index)) for index in (2, 3))
        return site, top, bottom
    if "^" in notation and notation.count("^") == 1 and notation.replace("^", "").isalpha():
        site = notation.replace("^", "")
        top = notation.index("^")
        return site, top, len(site) - top
    raise ValueError(f"no single pair of cut offsets in {notation!r}")


def parse_withrefm(text: str) -> dict[str, RebaseEntry]:
    """Read REBASE's `withrefm` file, keeping the records whose cut offsets it gives."""
    records: dict[str, RebaseEntry] = {}
    for block in re.split(r"\n\s*\n", text):
        fields = dict(_FIELD_RE.findall(block))
        if "1" not in fields or "3" not in fields:
            continue
        try:
            site, top, bottom = parse_site(fields["3"])
        except ValueError:
            continue
        name = fields["1"]
        records[name] = RebaseEntry(
            name=name,
            site=site,
            cuts=(top, bottom),
            same_specificity=tuple(n for n in fields.get("2", "").split(",") if n),
            suppliers=fields.get("7", "").strip(),
        )
    return records


def isoschizomers(name: str, records: dict[str, RebaseEntry]) -> tuple[str, ...]:
    """Return the enzymes a supplier sells that read AND cut exactly as `name` does.

    REBASE's own list holds every enzyme of the same specificity, which includes neoschizomers
    — same site, different cut — and hundreds nobody sells.
    """
    entry = records[name]
    return tuple(
        other
        for other in entry.same_specificity
        if (match := records.get(other))
        and match.suppliers
        and (match.site, match.cuts) == (entry.site, entry.cuts)
    )


def fetch_withrefm() -> str:
    """Download REBASE's `withrefm` file."""
    with urllib.request.urlopen(WITHREFM_URL, timeout=120) as response:
        return response.read().decode("latin-1")


def rebase_version(text: str) -> str:
    """Return the REBASE release the file announces."""
    version = re.search(r"REBASE version (\d+)", text)
    return version.group(1) if version else "unknown"


def build(text: str, properties: Mapping[str, Any]) -> dict[str, Any]:
    """Merge REBASE records with the curated supplier table into the shipped structure."""
    records = parse_withrefm(text)
    curated: Mapping[str, Mapping[str, Any]] = properties["enzymes"]
    enzymes = []
    for name in sorted(curated):
        entry = records.get(name)
        if entry is None:
            raise SystemExit(f"REBASE has no record with cut offsets for {name}")
        supplier_values = dict(curated[name])
        unverified = tuple(supplier_values.pop("unverified", ()))
        source = supplier_values.pop("source")
        field_sources = dict(supplier_values.pop("field_sources", {}))
        top, bottom = entry.cuts
        enzymes.append(
            {
                "name": name,
                "site": entry.site,
                "top_cut": top,
                "bottom_cut": bottom,
                "overhang_length": abs(bottom - top),
                "end": _end(top, bottom),
                "type": "IIS" if any(not 0 <= c <= len(entry.site) for c in entry.cuts) else "II",
                "isoschizomers": list(isoschizomers(name, records)),
                "suppliers": entry.suppliers,
                **{field: supplier_values.get(field) for field in SUPPLIER_FIELDS},
                "unverified": list(unverified),
                "provenance": {field: field_sources.get(field, source) for field in SUPPLIER_FIELDS}
                | dict.fromkeys(("site", "top_cut", "bottom_cut", "isoschizomers"), "rebase"),
            }
        )
    return {
        "rebase_version": rebase_version(text),
        "built": datetime.now(UTC).date().isoformat(),
        "sources": properties["sources"],
        "enzymes": enzymes,
    }


def _end(top: int, bottom: int) -> str:
    if top == bottom:
        return "blunt"
    return "5'" if bottom > top else "3'"


def main() -> int:
    """Write the data file, and say where it went."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--withrefm", type=Path, help="a downloaded REBASE withrefm file")
    parser.add_argument("--out", type=Path, default=DATA)
    arguments = parser.parse_args()
    text = arguments.withrefm.read_text("latin-1") if arguments.withrefm else fetch_withrefm()
    with PROPERTIES.open("rb") as handle:
        properties = tomllib.load(handle)
    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(json.dumps(build(text, properties), indent=2) + "\n")
    print(f"wrote {arguments.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

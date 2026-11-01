"""A sequence record in molviz's record format, and the options it is drawn with, as plain data.

molviz's format is the one definition of what it draws, so nothing here mirrors it as a type.
mbio hands over its biology: each feature with a CDS's translation, each primer with the bases
it does not pair with at each binding site, and every shipped enzyme's cuts, so that molviz
finds the unique cutters among them. A cut that falls off a linear record cuts nothing there and
is left out. A file's own colour is kept, SnapGene's ``noColor`` read as none, and a
``regulatory`` feature of a terminating class takes a terminator's colour, since molviz reads no
qualifiers. `molviz_options` writes every drawing's options, each enzyme by its shipped name.
The format's reference: <https://liuhlab.github.io/liulab-molviz/reference/record/>.
"""

import re
from collections.abc import Iterable
from typing import Any

from mbio.enzymes import enzymes, get_enzyme
from mbio.sequence import BindingSite, Feature, Primer, SequenceRecord, mismatches
from mbio.sites import find_sites
from mbio.translate import codon_start, translation

#: The ``/regulatory_class`` values that make a ``regulatory`` feature a terminator.
TERMINATING_CLASSES = frozenset({"terminator", "polyA_signal_sequence"})

#: The colour molviz gives a terminator, which a terminating ``regulatory`` feature takes.
TERMINATOR = "#ffaabb"

_HEX = re.compile(r"#[0-9a-fA-F]{6}")


def molviz_record(record: SequenceRecord) -> dict[str, Any]:
    """Return `record` in molviz's record format, leaving out a list that would be empty.

    Examples
    --------
    >>> molviz_record(SequenceRecord("ACGT", name="tiny"))
    {'name': 'tiny', 'circular': False, 'bases': 'ACGT'}
    """
    out: dict[str, Any] = {
        "name": record.name,
        "circular": record.topology == "circular",
        "bases": record.sequence,
    }
    for key, values in (
        ("features", [_feature(feature, record) for feature in record.features]),
        ("primers", [_primer(primer, record) for primer in record.primers]),
        ("cuts", _cuts(record)),
    ):
        if values:
            out[key] = values
    return out


def molviz_options(
    *,
    span: tuple[int, int] | None = None,
    linear: bool = False,
    sequence_view: bool = False,
    bases_per_row: int | None = None,
    both_strands: bool = True,
    features: bool = True,
    primers: bool = True,
    cut_sites: bool = True,
    hide_types: Iterable[str] = (),
    source: bool = False,
    enzymes: str | Iterable[str] | None = None,
    highlight: str | Iterable[str] = (),
) -> dict[str, Any]:
    """Return the options molviz draws a record or a figure with, each at its default left out.

    Each keyword is `mbio.plot.draw_map`'s, `span` being a region's span. molviz ignores an
    enzyme name that cuts nothing, so each name is handed over as the shipped enzyme it names:
    a commercial or isoschizomer name answers for it.

    Raises
    ------
    KeyError
        If no shipped enzyme answers to a name in `enzymes`.

    Examples
    --------
    >>> molviz_options(span=(10, 20), cut_sites=False)
    {'span': {'start': 10, 'end': 20}, 'kinds_off': ['cuts']}
    """
    off = (("features", features), ("primers", primers), ("cuts", cut_sites))
    types_off = sorted({*hide_types, *(() if source else ("source",))})
    lit = [highlight] if isinstance(highlight, str) else list(highlight)
    options: dict[str, Any] = {}
    if span is not None:
        options["span"] = {"start": span[0], "end": span[1]}
    if linear:
        options["map"] = "line"
    if sequence_view:
        options["sequence_view"] = True
    if bases_per_row is not None:
        options["bases_per_row"] = bases_per_row
    if not both_strands:
        options["bottom_strand"] = False
    if kinds_off := [kind for kind, on in off if not on]:
        options["kinds_off"] = kinds_off
    if types_off != ["source"]:
        options["types_off"] = types_off
    if enzymes is not None:
        named = (enzymes,) if isinstance(enzymes, str) else enzymes
        options["enzymes"] = list(dict.fromkeys(get_enzyme(name).name for name in named))
    if lit:
        options["highlight"] = lit
    return options


def _feature(feature: Feature, record: SequenceRecord) -> dict[str, Any]:
    own = _color(feature.color)
    segments = []
    for segment in feature.segments:
        written: dict[str, Any] = {"start": segment.start, "end": segment.end}
        if (color := _color(segment.color)) and color != own:
            written["color"] = color
        segments.append(written)
    out: dict[str, Any] = {
        "name": feature.name,
        "type": feature.type,
        "strand": int(feature.strand),
        "segments": segments,
    }
    if own := own or _terminating(feature):
        out["color"] = own
    # Read from the bases, never the `/translation` qualifier, which a join leaves stale.
    if feature.type == "CDS" and (protein := translation(feature, record)):
        out["translation"] = protein
        if (start := codon_start(feature)) != 1:
            out["codon_start"] = start
    return out


def _terminating(feature: Feature) -> str | None:
    classes = {str(one) for one in feature.qualifiers.get("regulatory_class", ())}
    if feature.type == "regulatory" and classes & TERMINATING_CLASSES:
        return TERMINATOR
    return None


def _primer(primer: Primer, record: SequenceRecord) -> dict[str, Any]:
    return {
        "name": primer.name,
        "bases": primer.sequence,
        "sites": [_site(primer, site, record) for site in primer.binding_sites],
    }


def _site(primer: Primer, site: BindingSite, record: SequenceRecord) -> dict[str, Any]:
    out: dict[str, Any] = {"start": site.start, "end": site.end, "strand": int(site.strand)}
    if found := mismatches(primer, site, record):
        out["mismatches"] = list(found)
    return out


def _cuts(record: SequenceRecord) -> list[dict[str, Any]]:
    return [
        {
            "enzyme": site.enzyme.name,
            "start": site.start,
            "end": site.end,
            "at": site.top_cut,
            "stagger": site.enzyme.bottom_cut - site.enzyme.top_cut,
        }
        for site in find_sites(record, enzymes())
        if site.cuts
    ]


def _color(value: str | None) -> str | None:
    """Return a file's colour in lower case, or ``None`` for none, SnapGene's ``noColor`` too."""
    return value.lower() if value is not None and _HEX.fullmatch(value) else None

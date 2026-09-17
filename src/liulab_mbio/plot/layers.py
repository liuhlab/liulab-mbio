"""What a map draws from a record: its items, with their names, colours and hover details.

A view takes these items and never the record. A map draws a record's features, each primer at
its binding sites, and the cut sites of the shipped unique cutters or of the enzymes named.

A feature draws in its file's colour, and a segment in its own where that differs. One the file
gives no colour takes Paul Tol's light scheme by the group its type falls in, and pale grey for a
type in none. Primers are purple and enzyme names black.

Positions a person reads, in labels and hover details, are 1-based and inclusive. A cut site is
numbered as SnapGene numbers one: by the base after which its enzymes cut the top strand.
"""

import re
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal

from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.enzymes import enzymes as shipped
from liulab_mbio.sequence import BindingSite, Feature, Primer, SequenceRecord, Strand
from liulab_mbio.sites import find_sites

#: What an item is: a feature, a primer at one binding site, or the enzymes cutting at one position.
type Kind = Literal["feature", "primer", "cut_site"]

_HEX = re.compile(r"#[0-9a-fA-F]{6}")

# The colours below are from Paul Tol's `tol_colors.py`, shipped under its licence:
#
# Copyright (c) 2022, Paul Tol
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of
#    conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of
#    conditions and the following disclaimer in the documentation and/or other materials
#    provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to
#    endorse or promote products derived from this software without specific prior written
#    permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
# AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

#: Each group of feature types, its colour in Paul Tol's light scheme, and the keys it holds:
#: INSDC's current and retired keys, and SnapGene's `LTR`. A type in no group takes `OTHER`.
GROUPS: Mapping[str, tuple[str, frozenset[str]]] = {
    "coding": (
        "#77AADD",
        frozenset({
            "CDS", "gene", "mat_peptide", "sig_peptide", "transit_peptide", "propeptide",
            "proprotein", "exon", "C_region", "D_segment", "J_segment", "N_region", "S_region",
            "V_region", "V_segment",
        }),
    ),
    "promoters and enhancers": (
        "#EE8866",
        frozenset({
            "promoter", "enhancer", "RBS", "-10_signal", "-35_signal", "CAAT_signal",
            "TATA_signal", "GC_signal", "attenuator", "misc_signal", "5'UTR", "regulatory",
        }),
    ),
    "terminators and polyA": (
        "#FFAABB",
        frozenset({"terminator", "polyA_signal", "polyA_site", "3'UTR"}),
    ),
    "replication origins": ("#EEDD88", frozenset({"rep_origin", "oriT"})),
    "binding sites": ("#99DDFF", frozenset({"protein_bind", "misc_binding", "primer_bind"})),
    "RNA": (
        "#44BB99",
        frozenset({
            "ncRNA", "misc_RNA", "mRNA", "tRNA", "rRNA", "tmRNA", "precursor_RNA",
            "prim_transcript",
        }),
    ),
    "repeats and recombination": (
        "#BBCC33",
        frozenset({
            "LTR", "repeat_region", "mobile_element", "misc_recomb", "stem_loop", "D-loop",
        }),
    ),
}  # fmt: skip

#: Pale grey, for a type in no group.
OTHER = "#DDDDDD"

#: Tol bright purple, for primers and for a `primer_bind` feature the file gives no colour.
PRIMER = "#AA3377"

#: The colour of an enzyme's name.
ENZYME = "#000000"

#: What joins the names of the enzymes cutting at one position.
SEPARATOR = " - "

#: The `/regulatory_class` values that put a `regulatory` feature among the terminators.
TERMINATING_CLASSES = frozenset({"terminator", "polyA_signal_sequence"})


@dataclass(frozen=True, slots=True)
class Span:
    """One stretch of an item, 0-based and half-open as `Segment` is, and the colour it fills."""

    start: int
    end: int
    color: str


@dataclass(frozen=True, slots=True)
class Cutter:
    """An enzyme a cut site names, and how many times it cuts the whole record."""

    name: str
    cuts: int

    @property
    def unique(self) -> bool:
        """Whether it cuts the record once, which is what makes its name bold."""
        return self.cuts == 1


@dataclass(frozen=True, slots=True)
class Item:
    """One thing a map draws.

    Parameters
    ----------
    kind
        What it is.
    name
        As the file writes it; for a cut site, its enzymes' names joined as its label joins them.
    type
        The feature type; ``primer`` or ``cut site`` for the other kinds.
    strand
        The strand it reads along.
    spans
        In top-strand order, as the feature's segments are. A primer's is its binding site. A cut
        site's is empty, at the boundary where its enzymes cut the top strand.
    label
        The text a map labels it with.
    hover
        What hovering over it shows: its name, type, span and length, as a person reads them.
    cutters
        A cut site's enzymes, in the order its label names them.
    """

    kind: Kind
    name: str
    type: str
    strand: Strand
    spans: tuple[Span, ...]
    label: str
    hover: Mapping[str, str] = field(default_factory=dict, hash=False)
    cutters: tuple[Cutter, ...] = ()

    @property
    def color(self) -> str:
        """The colour of its label: its first span's."""
        return self.spans[0].color

    @property
    def runs(self) -> tuple[tuple[str, bool], ...]:
        """Its label in stretches of text, each marked bold or not: a unique cutter's name is."""
        if not self.cutters:
            return ((self.label, False),)
        runs: list[tuple[str, bool]] = []
        for cutter in self.cutters:
            if runs:
                runs.append((SEPARATOR, False))
            runs.append((cutter.name, cutter.unique))
        runs.append((self.label.removeprefix(self.name), False))
        return tuple(runs)


def items(
    record: SequenceRecord,
    *,
    features: bool = True,
    primers: bool = True,
    cut_sites: bool = True,
    enzymes: str | Iterable[str] | None = None,
    hide_types: Iterable[str] = (),
    source: bool = False,
) -> tuple[Item, ...]:
    """Return what a map of `record` draws: its features, its primers, then its cut sites.

    Parameters
    ----------
    record
        What is drawn.
    features, primers, cut_sites
        Whether each is drawn.
    enzymes
        The names of the enzymes whose every cut site is drawn; the shipped unique cutters when
        ``None``. How many times an enzyme cuts is counted over the whole record.
    hide_types
        The feature types left off.
    source
        Whether a `source` feature is drawn.

    Raises
    ------
    KeyError
        If no shipped enzyme answers to a name in `enzymes`, as `sites.find_sites` raises.
    """
    length = len(record)
    chosen = None if enzymes is None else _chosen(enzymes)
    left_off = set(hide_types) if source else {*hide_types, "source"}
    drawn = []
    if features:
        drawn.extend(
            _feature(feature, length) for feature in record.features if feature.type not in left_off
        )
    if primers:
        drawn.extend(
            _primer(primer, site, length)
            for primer in record.primers
            for site in primer.binding_sites
        )
    if cut_sites:
        drawn.extend(_cut_sites(record, chosen))
    return tuple(drawn)


def merge_cuts(cuts: Iterable[tuple[str, int]], length: int) -> tuple[Item, ...]:
    """Return a cut site for each position cut, naming every enzyme that cuts there.

    Parameters
    ----------
    cuts
        Each cut as an enzyme's name and the 0-based boundary where it cuts the top strand, one
        for every site, so a name given twice cuts twice.
    length
        The record's.

    Returns
    -------
    tuple[Item, ...]
        In order of position, each naming its enzymes alphabetically, as SnapGene merges them.

    Examples
    --------
    >>> [site.label for site in merge_cuts([("SacI", 406), ("BanII", 406), ("EcoRI", 396)], 2686)]
    ['EcoRI (396)', 'BanII - SacI (406)']
    """
    listed = list(cuts)
    counts = Counter(name for name, _ in listed)
    at: dict[int, set[str]] = {}
    for name, position in listed:
        at.setdefault(position % length, set()).add(name)
    merged = []
    for position, names in sorted(at.items()):
        cutters = tuple(
            Cutter(one, counts[one]) for one in sorted(names, key=lambda one: (one.casefold(), one))
        )
        name = SEPARATOR.join(cutter.name for cutter in cutters)
        shown = str((position - 1) % length + 1)
        merged.append(
            Item(
                "cut_site",
                name,
                "cut site",
                Strand.NONE,
                (Span(position, position, ENZYME),),
                f"{name} ({shown})",
                {"name": name, "type": "cut site", "span": shown},
                cutters,
            )
        )
    return tuple(merged)


def default_color(feature: Feature) -> str:
    """Return the colour a feature with none of its own takes, by the group its type falls in.

    A `primer_bind` feature takes the primers' purple instead. A `regulatory` feature is grouped
    by its `/regulatory_class`: a terminator's or a polyA signal's among the terminators, any other
    among the promoters and enhancers.

    Examples
    --------
    >>> from liulab_mbio.sequence import Segment
    >>> default_color(Feature("ori", "rep_origin", (Segment(0, 10),)))
    '#eedd88'
    """
    if feature.type == "primer_bind":
        return PRIMER.lower()
    if feature.type == "regulatory" and TERMINATING_CLASSES & set(
        map(str, feature.qualifiers.get("regulatory_class", ()))
    ):
        return GROUPS["terminators and polyA"][0].lower()
    for color, types in GROUPS.values():
        if feature.type in types:
            return color.lower()
    return OTHER.lower()


def unwrapped(spans: Sequence[Span], length: int) -> tuple[tuple[int, int], ...]:
    """Return each span's start and end, counted round a circular record from the first span.

    A span that starts before the one before it ends lies past the origin, so it moves on by the
    record's length, as SnapGene reads a feature's segments.
    """
    counted: list[tuple[int, int]] = []
    for span in spans:
        start = span.start
        if counted and start < counted[-1][1]:
            start += length
        counted.append((start, start + span.end - span.start))
    return tuple(counted)


def span_text(start: int, end: int, length: int) -> str:
    """Return a span as a person reads it, 1-based and inclusive, across the origin as needed.

    Examples
    --------
    >>> span_text(2683, 2689, 2686)
    '2684 .. 3'
    """
    return f"{start % length + 1} .. {(end - 1) % length + 1}"


def text_color(fill: str) -> str:
    """Return black or white, whichever contrasts more with `fill`, black when they tie.

    Examples
    --------
    >>> text_color("#ffffff"), text_color("#993366")
    ('#000000', '#ffffff')
    """
    luminance = _luminance(fill)
    return "#000000" if (luminance + 0.05) / 0.05 >= 1.05 / (luminance + 0.05) else "#ffffff"


def outline_color(fill: str) -> str:
    """Return the colour an arrow filled with `fill` is outlined in: the fill, darkened.

    Darkened this far, even a white arrow's outline contrasts with a white page.
    """
    red, green, blue = _channels(fill)
    return "#{:02x}{:02x}{:02x}".format(*(round(0.45 * channel) for channel in (red, green, blue)))


def _feature(feature: Feature, length: int) -> Item:
    own = _given(feature.color)
    default = default_color(feature)
    spans = tuple(
        Span(segment.start, segment.end, _given(segment.color) or own or default)
        for segment in feature.segments
    )
    runs: list[list[int]] = []
    for start, end in unwrapped(spans, length):
        if runs and start == runs[-1][1]:
            runs[-1][1] = end
        else:
            runs.append([start, end])
    bases = sum(segment.end - segment.start for segment in feature.segments)
    hover = {
        "name": feature.name,
        "type": feature.type,
        "span": ", ".join(span_text(start, end, length) for start, end in runs),
        "length": f"{bases} bp",
    }
    return Item("feature", feature.name, feature.type, feature.strand, spans, feature.name, hover)


def _primer(primer: Primer, site: BindingSite, length: int) -> Item:
    span = span_text(site.start, site.end, length)
    hover = {
        "name": primer.name,
        "type": "primer",
        "span": span,
        "length": f"{site.end - site.start} bp",
    }
    return Item(
        "primer",
        primer.name,
        "primer",
        site.strand,
        (Span(site.start, site.end, PRIMER.lower()),),
        f"{primer.name} ({span})",
        hover,
    )


def _chosen(names: str | Iterable[str]) -> tuple[Enzyme, ...]:
    """Return the shipped enzymes `names` answer to, each once however many names it has."""
    named = (names,) if isinstance(names, str) else names
    return tuple({enzyme.name: enzyme for enzyme in map(get_enzyme, named)}.values())


def _cut_sites(record: SequenceRecord, chosen: tuple[Enzyme, ...] | None) -> tuple[Item, ...]:
    """Return the cut sites of `chosen`, or of the shipped unique cutters when that is ``None``."""
    # A site whose cut falls off a linear record cuts nothing there.
    found = [
        site for site in find_sites(record, shipped() if chosen is None else chosen) if site.cuts
    ]
    counts = Counter(site.enzyme.name for site in found)
    return merge_cuts(
        (
            (site.enzyme.name, site.top_cut)
            for site in found
            if chosen is not None or counts[site.enzyme.name] == 1
        ),
        len(record),
    )


def _given(color: str | None) -> str | None:
    """Return a file's colour lower-cased, or `None` for none, SnapGene's `noColor` included."""
    return color.lower() if color is not None and _HEX.fullmatch(color) else None


def _channels(color: str) -> tuple[int, int, int]:
    return int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)


def _luminance(color: str) -> float:
    """WCAG 2's relative luminance."""
    linear = [
        value / 12.92 if value <= 0.04045 else ((value + 0.055) / 1.055) ** 2.4
        for value in (channel / 255 for channel in _channels(color))
    ]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]

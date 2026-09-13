"""Restriction enzymes: the shipped records and how their cut offsets are written.

A cut offset is a 0-based boundary counted from the first base of the recognition site as
written 5' to 3', like a slice point: offset ``k`` severs a strand between site bases ``k - 1``
and ``k``. `Enzyme.top_cut` is the cut in the strand carrying the site, `Enzyme.bottom_cut` the
cut in the complementary strand, written in that same frame — the bottom strand is severed
between the partners of bases ``k - 1`` and ``k``.

An offset outside ``0..len(site)`` is a cut outside the site, which is what makes an enzyme
Type IIS. REBASE's ``GGTCTC(1/5)`` for BsaI is ``top_cut=7, bottom_cut=11``.
"""

import json
from collections.abc import Mapping
from dataclasses import KW_ONLY, dataclass, field
from functools import cache
from importlib.resources import files
from typing import Literal

from liulab_mbio.sequence import IUPAC_DNA, Segment, SequenceRecord, Strand

#: Which strand the surviving single strand of a cut belongs to, or neither.
type EndType = Literal["5'", "3'", "blunt"]

#: Type IIS cuts outside its recognition site; Type II cuts within it.
type EnzymeType = Literal["II", "IIS"]


@dataclass(frozen=True, slots=True)
class Enzyme:
    """One restriction enzyme and where it cuts.

    Parameters
    ----------
    name
        The REBASE name, such as ``"BsaI"``.
    site
        The recognition site, IUPAC, 5' to 3'. Stored upper-case.
    top_cut, bottom_cut
        Cut offsets, as described in the module docstring.
    isoschizomers
        Enzymes a supplier sells that read and cut exactly this site. Neoschizomers — the same
        site cut elsewhere, as XmaI cuts SmaI's — are not among them.
    commercial_name, catalog_number, supplier
        The product the other supplier fields describe, such as ``"BsaI-HFv2"``, ``"R3733"``.
    incubation_celsius
        The supplier's digestion temperature, which is not always a protocol's temperature.
    heat_inactivation_celsius, heat_inactivation_minutes
        ``None`` when the supplier says heat does not inactivate the enzyme, and also when it
        says nothing — which is what `unverified` then names.
    methylation
        The supplier's Dam, Dcm and CpG sensitivity, keyed ``"dam"``, ``"dcm"``, ``"cpg"``.
    unverified
        Fields no source stated. Never a guess.

    Raises
    ------
    ValueError
        If `site` is empty or holds a letter outside `IUPAC_DNA`.
    """

    name: str
    site: str
    _: KW_ONLY
    top_cut: int
    bottom_cut: int
    isoschizomers: tuple[str, ...] = ()
    commercial_name: str | None = None
    catalog_number: str | None = None
    supplier: str | None = None
    incubation_celsius: int | None = None
    heat_inactivation_celsius: int | None = None
    heat_inactivation_minutes: int | None = None
    methylation: Mapping[str, str] = field(default_factory=dict, hash=False)
    unverified: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Upper-case the site and refuse one that is not IUPAC DNA."""
        site = self.site.upper()
        if not site:
            raise ValueError(f"enzyme {self.name!r} has no recognition site")
        if bad := set(site) - IUPAC_DNA:
            raise ValueError(f"{self.name} site is not IUPAC DNA: {''.join(sorted(bad))}")
        object.__setattr__(self, "site", site)

    @property
    def type(self) -> EnzymeType:
        """``"IIS"`` when either cut falls outside the site, ``"II"`` otherwise."""
        return "IIS" if any(not 0 <= cut <= len(self.site) for cut in self.cut_offsets) else "II"

    @property
    def cut_offsets(self) -> tuple[int, int]:
        """The two cut offsets, top strand first."""
        return self.top_cut, self.bottom_cut

    @property
    def overhang_length(self) -> int:
        """How many bases the two cuts leave single-stranded."""
        return abs(self.bottom_cut - self.top_cut)

    @property
    def end(self) -> EndType:
        """The end left behind: ``"5'"``, ``"3'"`` or ``"blunt"``.

        Examples
        --------
        >>> Enzyme("KpnI", "GGTACC", top_cut=5, bottom_cut=1).end
        "3'"
        """
        if self.bottom_cut == self.top_cut:
            return "blunt"
        return "5'" if self.bottom_cut > self.top_cut else "3'"

    def cut_positions(self, start: int, strand: Strand = Strand.FORWARD) -> tuple[int, int]:
        """Where this enzyme cuts a record holding its site at `start`.

        Parameters
        ----------
        start
            The 0-based index where the site begins on the top strand. A reverse-strand site is
            the reverse complement of `site` there.
        strand
            Which strand carries the site.

        Returns
        -------
        tuple[int, int]
            The top-strand cut and the bottom-strand cut, in record coordinates and in that
            order whatever the strand. Either may fall outside the record, which is what a
            circular record's caller reduces modulo its length.

        Examples
        --------
        >>> Enzyme("BsaI", "GGTCTC", top_cut=7, bottom_cut=11).cut_positions(4)
        (11, 15)
        """
        if strand == Strand.REVERSE:
            end = start + len(self.site)
            return end - self.bottom_cut, end - self.top_cut
        return start + self.top_cut, start + self.bottom_cut

    def overhang(self, record: SequenceRecord, start: int, strand: Strand = Strand.FORWARD) -> str:
        """Return the top-strand bases this enzyme leaves single-stranded, 5' to 3'.

        Empty for a blunt cutter. On a circular record the overhang reads across the origin.

        Raises
        ------
        ValueError
            If either cut falls off the end of a linear record.
        """
        top, bottom = self.cut_positions(start, strand)
        low, high = min(top, bottom), max(top, bottom)
        if low == high:
            return ""
        length = len(record)
        if record.topology == "circular":
            low, high = low % length, low % length + (high - low)
        elif not (low >= 0 and high <= length):
            raise ValueError(f"{self.name} cut at {low}-{high} falls off a linear record")
        return record.extract(Segment(low, high))


@cache
def _shipped() -> tuple[tuple[Enzyme, ...], Mapping[str, tuple[Enzyme, ...]]]:
    data = json.loads((files("liulab_mbio") / "data" / "enzymes.json").read_text())
    records = tuple(
        Enzyme(
            entry["name"],
            entry["site"],
            top_cut=entry["top_cut"],
            bottom_cut=entry["bottom_cut"],
            isoschizomers=tuple(entry["isoschizomers"]),
            commercial_name=entry["commercial_name"],
            catalog_number=entry["catalog_number"],
            supplier=entry["supplier"],
            incubation_celsius=entry["incubation_celsius"],
            heat_inactivation_celsius=entry["heat_inactivation_celsius"],
            heat_inactivation_minutes=entry["heat_inactivation_minutes"],
            methylation=dict(entry["methylation"]),
            unverified=tuple(entry["unverified"]),
        )
        for entry in data["enzymes"]
    )
    index: dict[str, tuple[Enzyme, ...]] = {}
    for record in records:
        for shared in record.isoschizomers:
            index[shared.casefold()] = (*index.get(shared.casefold(), ()), record)
    # Written second, and never accumulated: an enzyme shipped under a name owns that name,
    # whatever else lists it as an isoschizomer.
    for record in records:
        for owned in (record.commercial_name, record.name):
            if owned is not None:
                index[owned.casefold()] = (record,)
    return records, index


def enzymes() -> tuple[Enzyme, ...]:
    """Return every shipped enzyme, ordered by name."""
    return _shipped()[0]


def get_enzyme(name: str) -> Enzyme:
    """Return one enzyme by its name, its commercial name, or an isoschizomer's name.

    Raises
    ------
    KeyError
        If no shipped enzyme answers to `name`, or if `name` is an isoschizomer of more than
        one of them, which would leave the supplier properties ambiguous.

    Examples
    --------
    >>> get_enzyme("BsaI-HFv2").site
    'GGTCTC'
    """
    found = _shipped()[1].get(name.casefold())
    if not found:
        raise KeyError(f"no shipped enzyme is called {name!r}")
    enzyme, *shared = found
    if shared:
        names = ", ".join(sorted(record.name for record in found))
        raise KeyError(f"{name!r} is an isoschizomer of more than one shipped enzyme: {names}")
    return enzyme

"""Build each part's synthesis sequence, and the sheet a vendor is given.

One part is one synthesised block: the 5' external stuffer, the coding sequence, the internal
stuffer for its position, the barcode naming it, and the 3' external stuffer. The external enzyme
releases the part from the block on the scheme's two overhangs, and the internal enzyme excises
the stuffer the part carries, which is what the next round opens it on.

The overhangs written into a block are the standard's rather than the scheme's own, because a
scheme's stuffers carry whatever overhangs it was written with and `design_standard` chooses the
set the part lists cost least. Nothing else in a stuffer moves.

A block's coding bases stop where the overhang either side of it spells the rest: an overhang
finishes the upstream part's last codon and spells whole codons of its own, so a part's own bases
are its amino acids less the ones charged to a junction. `liulab_mbio.library.standard` decides
which those are, and this writes them.
"""

from collections.abc import Mapping, Sequence
from dataclasses import KW_ONLY, dataclass
from typing import NoReturn

from liulab_mbio.barcodes import (
    GC_BAND,
    MAX_HOMOPOLYMER,
    MIN_DISTANCE,
    SEED,
    BarcodeRules,
    design_barcodes,
)
from liulab_mbio.codons import CodonUsage, codon_usage
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.library.scheme import Scheme
from liulab_mbio.library.standard import End, PartList, Standard, Terminus, junction_residues
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand
from liulab_mbio.sites import CutSite, Domestication, domesticate, find_sites
from liulab_mbio.translate import SiteNotRemovableError, reverse_translate

#: The columns of the synthesis order sheet.
SHEET_COLUMNS = ("name", "sequence", "length", "position", "barcode")

#: The columns of the amino-acid change table, wild type beside synthesised.
CHANGE_COLUMNS = ("part", "position", "end", "wild_type", "synthesised")


@dataclass(frozen=True, slots=True)
class Part:
    """One synthesised block, as it is ordered and as it reads.

    Parameters
    ----------
    name
        The member of the part list, as it is named and ordered.
    position
        The position this part fills.
    index
        Which position that is, counting from zero, so a round knows when this part goes in.
    sequence
        The whole block, 5' to 3', as a vendor synthesises it.
    barcode
        The bases naming this part in the finished barcode block.
    protein
        The amino acids the product reads here, the junctions included: what the part now is,
        rather than what it was.
    coding
        Where the part's own coding bases lie in `sequence`. The bases the junctions spell lie
        outside it, in the stuffers either side.
    changes
        Each synonymous codon change that took a forbidden site out of the coding bases.
    """

    name: str
    position: str
    _: KW_ONLY
    index: int
    sequence: str
    barcode: str
    protein: str
    coding: Segment
    changes: tuple[Domestication, ...] = ()

    @property
    def length(self) -> int:
        """How many bases are ordered."""
        return len(self.sequence)

    @property
    def coding_sequence(self) -> str:
        """The part's own coding bases, short of a whole codon where a junction takes the rest."""
        return self.sequence[self.coding.start : self.coding.end]


def barcode_phase(scheme: Scheme) -> int:
    """How many bases into a codon the finished barcode block begins.

    The terminal position's internal stuffer is never excised, and everything ahead of it is whole
    codons, so the block opens that stuffer's length into a codon. Every barcode of the block sits
    at this phase, a barcode and its scar being whole codons together.
    """
    return len(scheme.internal_stuffer(-1)) % 3


def barcode_rules(
    scheme: Scheme,
    *,
    distance: int = MIN_DISTANCE,
    max_homopolymer: int | None = MAX_HOMOPOLYMER,
    gc_band: tuple[float, float] | None = GC_BAND,
) -> BarcodeRules:
    """Return the rules this scheme's barcodes hold to, its phase read off its own geometry.

    The dials are `liulab_mbio.barcodes`' own; the length, the scar, the phase and the forbidden
    enzymes are the scheme's and are not for a caller to restate.
    """
    return BarcodeRules(
        scheme.barcode_length,
        scar=scheme.cloning_scar,
        phase=barcode_phase(scheme),
        forbidden=_enzymes(scheme),
        distance=distance,
        max_homopolymer=max_homopolymer,
        gc_band=gc_band,
    )


def design_parts(
    scheme: Scheme,
    part_lists: Sequence[PartList],
    standard: Standard,
    *,
    host: str,
    rules: BarcodeRules | None = None,
    seed: int = SEED,
) -> tuple[Part, ...]:
    """Build every part's synthesis sequence, in position order and then list order.

    Each part is written into the scheme's stuffers with the standard's overhangs, given a barcode
    of its own, and checked to carry the scheme's sites and no others. The same inputs return the
    same bases: the barcodes are drawn from `seed`, and nothing else here is random.

    Parameters
    ----------
    scheme
        The architecture the parts are built into.
    part_lists
        One per position, in the scheme's order: each member's name and the protein it codes for.
    standard
        The overhang standard the build works to, and what it charges each part. Design it with
        `liulab_mbio.library.standard.design_standard` over these same part lists.
    host
        The name of the codon usage table the coding bases are written for.
    rules
        What every barcode holds to. `barcode_rules` for this scheme by default.
    seed
        The seed the barcodes are drawn with, offset by the position, so one seed settles a build.

    Raises
    ------
    ValueError
        If the part lists do not match the scheme, the standard was not designed for them, a
        junction leaves a part no coding bases of its own, or a block spells a site the scheme
        does not expect. `liulab_mbio.translate.SiteNotRemovableError` where that site lies in a
        part's coding bases, naming the part and the site.
    KeyError
        If no shipped codon usage table is called `host`, or the scheme names an enzyme this
        package does not ship.
    liulab_mbio.barcodes.SpaceExhaustedError
        If a part list is larger than the barcodes its rules allow.
    """
    _check_lists(scheme, part_lists)
    _check_standard(scheme, standard)
    design = _Design(
        scheme,
        standard,
        table=codon_usage(host),
        host=host,
        enzymes=_enzymes(scheme),
        charged={(one.position, one.part, one.end): one for one in standard.termini},
    )
    held = rules if rules is not None else barcode_rules(scheme)
    made: list[Part] = []
    for index, parts in enumerate(part_lists):
        codes = design_barcodes(len(parts), held, seed=seed + index)
        made.extend(
            _built(design, index, name, protein, barcode)
            for (name, protein), barcode in zip(parts.items(), codes, strict=True)
        )
    return tuple(made)


def synthesis_sheet(parts: Sequence[Part]) -> str:
    """Return these parts as a tab-separated sheet, one row each, in the order given.

    The columns are `SHEET_COLUMNS`: the name to order it under, the whole block 5' to 3', its
    length, the position it fills and its barcode. The barcode stands on the same row, so the
    sheet ordered from is what decodes the sequencing afterwards.
    """
    rows = ["\t".join(SHEET_COLUMNS)]
    rows.extend(
        "\t".join((part.name, part.sequence, str(part.length), part.position, part.barcode))
        for part in parts
    )
    return "\n".join(rows) + "\n"


def change_table(standard: Standard) -> str:
    """Return every amino acid the overhang standard moved, wild type beside synthesised.

    The columns are `CHANGE_COLUMNS`. A standard that changed nothing returns the header alone.
    """
    rows = ["\t".join(CHANGE_COLUMNS)]
    rows.extend(
        "\t".join((one.part, one.position, one.end, one.wild_type, one.synthesised))
        for one in standard.changes
    )
    return "\n".join(rows) + "\n"


@dataclass(frozen=True, slots=True)
class _Design:
    """What every part of one build is written from."""

    scheme: Scheme
    standard: Standard
    _: KW_ONLY
    table: CodonUsage
    host: str
    enzymes: tuple[Enzyme, ...]
    charged: Mapping[tuple[str, str, End], Terminus]


def _enzymes(scheme: Scheme) -> tuple[Enzyme, ...]:
    """Every enzyme a block is allowed to spell a site for, and nowhere but its stuffers."""
    return (scheme.internal, scheme.external, *scheme.blunt)


def _check_lists(scheme: Scheme, part_lists: Sequence[PartList]) -> None:
    """Refuse part lists that are not one per position, or a position with nothing to fill it."""
    names = [position.name for position in scheme.positions]
    if len(part_lists) != len(names):
        raise ValueError(
            f"this scheme has {len(names)} position(s) — {', '.join(names)} — and "
            f"{len(part_lists)} part list(s) were given"
        )
    for name, parts in zip(names, part_lists, strict=True):
        if not parts:
            raise ValueError(f"the part list for position {name!r} holds no part")


def _check_standard(scheme: Scheme, standard: Standard) -> None:
    """Refuse a standard that was not designed for this scheme."""
    if len(standard.entry_overhangs) != scheme.position_count:
        raise ValueError(
            f"this scheme has {scheme.position_count} position(s) and the standard holds "
            f"{len(standard.entry_overhangs)} entry overhang(s)"
        )
    if standard.scar_overhang != scheme.cloning_scar:
        raise ValueError(
            f"the standard's cloning scar {standard.scar_overhang!r} is not this scheme's "
            f"{scheme.cloning_scar!r}: every part's 3' end leaves the scheme's own"
        )


def _built(design: _Design, index: int, name: str, protein: str, barcode: str) -> Part:
    """Write one part's block, take any forbidden site out of its coding bases, and check it."""
    scheme, standard = design.scheme, design.standard
    position = scheme.positions[index]
    entry = standard.entry_overhangs[index]
    following = standard.entry_overhangs[(index + 1) % scheme.position_count]
    donated = junction_residues(len(entry))[0]
    head = design.charged.get((position.name, name, "5'"))
    tail = design.charged.get((position.name, name, "3'"))
    _check_charged(design, index, name, position.name, tail, donated)
    protein = protein.upper()
    synthesised = _synthesised(protein, head, tail)
    coding = _coding(design, name, position.name, synthesised, head, tail, following, donated)
    stuffer = position.internal_stuffer_prefix[: -len(following)] + following
    regions = (
        position.external_stuffer_5[: -len(entry)] + entry,
        coding,
        stuffer + scheme.internal_stuffer_core,
        barcode,
        standard.scar_overhang + position.external_stuffer_3[len(standard.scar_overhang) :],
    )
    at = _offsets(regions)
    span = Segment(at[1], at[2])
    record = _record(name, "".join(regions), span, 3 * _residues(synthesised, head, tail))
    cleaned, report = domesticate(record, design.enzymes, usage=design.table)
    _check_sites(design, cleaned, name, position.name, regions, at)
    return Part(
        name,
        position.name,
        index=index,
        sequence=cleaned.sequence,
        barcode=barcode,
        protein=synthesised,
        coding=span,
        changes=report.changes,
    )


def _offsets(regions: Sequence[str]) -> tuple[int, ...]:
    """Where each region begins in the block, and where the last one ends."""
    at = [0]
    for bases in regions:
        at.append(at[-1] + len(bases))
    return tuple(at)


def _residues(synthesised: str, head: Terminus | None, tail: Terminus | None) -> int:
    """How many of a part's amino acids its own coding bases spell whole codons for."""
    first = len(head.wild_type) if head is not None else 0
    last = len(synthesised) - (len(tail.wild_type) if tail is not None else 0)
    return max(0, last - first)


def _synthesised(protein: str, head: Terminus | None, tail: Terminus | None) -> str:
    """Return the protein as the standard spells it, its charged ends replaced."""
    first = len(head.wild_type) if head is not None else 0
    last = len(protein) - (len(tail.wild_type) if tail is not None else 0)
    lead = head.synthesised if head is not None else ""
    trail = tail.synthesised if tail is not None else ""
    return lead + protein[first:last] + trail


def _coding(
    design: _Design,
    name: str,
    position: str,
    synthesised: str,
    head: Terminus | None,
    tail: Terminus | None,
    following: str,
    donated: int,
) -> str:
    """Return the coding bases a block carries: whole codons, then what the next junction takes.

    The amino acids a junction spells live in the stuffers either side, so they are not written
    here. The last one the part still owns keeps only the bases the overhang does not supply.
    """
    whole = _residues(synthesised, head, tail)
    first = len(head.wild_type) if head is not None else 0
    bases = reverse_translate(synthesised[first : first + whole], host=design.host) if whole else ""
    if tail is not None and donated:
        codon = _junction_codon(
            design, name, position, synthesised[first + whole], following[: (3 - donated) % 3]
        )
        bases += codon[:donated]
    if not bases:
        raise ValueError(
            f"part {name!r} at position {position!r} keeps no coding bases of its own: the "
            f"junctions either side spell all {len(synthesised)} of its amino acids"
        )
    return bases


def _junction_codon(design: _Design, name: str, position: str, residue: str, suffix: str) -> str:
    """Return the codon this host spells `residue` with most often, of those ending in `suffix`."""
    every = design.table.synonymous(reverse_translate(residue, host=design.host))
    found = next((codon for codon in every if codon.endswith(suffix)), None)
    if found is None:
        raise ValueError(
            f"part {name!r} at position {position!r} ends on {residue}, which no codon spells "
            f"ending in {suffix!r}: design the standard over these same part lists"
        )
    return found


def _check_charged(
    design: _Design,
    index: int,
    name: str,
    position: str,
    tail: Terminus | None,
    donated: int,
) -> None:
    """Refuse a standard that charges this part no residue where its junction must take one."""
    terminal = index == design.scheme.position_count - 1
    if donated and not terminal and tail is None:
        raise ValueError(
            f"the standard charges part {name!r} at position {position!r} nothing at its 3' end, "
            "where this scheme's junction takes a residue: design it over these same part lists"
        )


def _record(name: str, block: str, span: Segment, whole: int) -> SequenceRecord:
    """Return the block, the codons a synonymous change may move annotated as a coding sequence.

    The bases a junction takes are left outside it, so domestication cannot move an overhang.
    """
    features = (
        (Feature(name, "CDS", (Segment(span.start, span.start + whole),), strand=Strand.FORWARD),)
        if whole
        else ()
    )
    return SequenceRecord(block, name=name, features=features)


def _check_sites(
    design: _Design,
    record: SequenceRecord,
    name: str,
    position: str,
    regions: Sequence[str],
    at: Sequence[int],
) -> None:
    """Refuse a block spelling a site anywhere but where the scheme's stuffers put one.

    The stuffers are searched on their own and their hits moved to where the block puts them, so
    a site the coding bases, the barcode or a junction between regions spells is the difference.
    """
    expected = {
        (site.enzyme.name, at[which] + site.start, int(site.strand))
        for which in (0, 2, 4)
        for site in find_sites(SequenceRecord(regions[which]), design.enzymes)
    }
    found = {
        (site.enzyme.name, site.start, int(site.strand)): site
        for site in find_sites(record, design.enzymes)
    }
    if strange := sorted(set(found) - expected):
        _refuse(name, position, found[strange[0]], at)
    if missing := sorted(expected - set(found)):
        enzyme, start, _ = missing[0]
        raise ValueError(
            f"part {name!r} at position {position!r} lost the {enzyme} site its stuffer puts at "
            f"{start}: the block cannot be cut as the scheme says"
        )


def _refuse(name: str, position: str, site: CutSite, at: Sequence[int]) -> NoReturn:
    """Name the part, the site and the region it fell in, and say which kind of refusal it is."""
    coding = at[1] <= site.start and site.end <= at[2]
    where = (
        "its coding sequence"
        if coding
        else "its barcode"
        if at[3] <= site.start and site.end <= at[4]
        else "a junction between the regions of its block"
    )
    strand = "forward" if site.strand == Strand.FORWARD else "reverse"
    message = (
        f"part {name!r} at position {position!r}: {site.enzyme.name} reads a site at "
        f"{site.start} on the {strand} strand, in {where}, where the scheme expects none"
    )
    raise SiteNotRemovableError(message) if coding else ValueError(message)

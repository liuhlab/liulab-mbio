"""iGGA: the one method this package plans, as constants, checked once when they are imported.

`IGGA` is the method. Its overhangs, enzymes and stuffers are the DNA of molecules already on the
shelf — the four ccdB cassettes, the DMX vector, the part carrier — so turning one of them would
not change the tube, which is what makes it a constant rather than an input.
`docs/adr/0010-method-in-code.md` says why, and `liulab_synbio.igga.project` holds what a
project chooses instead.

Nothing here states the overhang a design works from. A part's entry overhang is read off the 5'
external stuffer by cutting it with the external enzyme, and the cloning scar off the 3' one, so
the overhangs are the ones the ordered DNA will leave. `INTERFACE_OVERHANGS` names the four the
method's own junctions carry, and the cargo pair is checked against the DNA.
"""

from collections.abc import Mapping
from dataclasses import KW_ONLY, dataclass
from types import MappingProxyType
from typing import Literal, NoReturn

from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.sites import find_sites

#: What the method checks when it is imported, and what a project checks when it is read. A
#: refusal leads with the name of the invariant it broke.
type Invariant = Literal[
    "internal-stuffer-prefix",
    "internal-stuffer-cuts",
    "external-stuffer-5",
    "external-stuffer-3",
    "barcode-frame",
    "terminal-frame",
    "enzyme-regions",
    "interface-overhangs",
]


def refuse(invariant: Invariant, detail: str) -> NoReturn:
    """Refuse, naming the invariant first so a caller reports which one broke."""
    raise ValueError(f"{invariant}: {detail}")


def _cuts(bases: str, enzyme: Enzyme) -> bool:
    """Whether `enzyme` reads a site in `bases`."""
    return bool(find_sites(SequenceRecord(bases), enzyme))


def _overhangs(bases: str, enzyme: Enzyme) -> set[str]:
    """Return every overhang `enzyme` leaves inside `bases`, as the top strand spells it.

    A cut falling off the end leaves nothing to read and is not one of these, and neither is the
    blunt end of a chopper.
    """
    return {site.overhang for site in find_sites(SequenceRecord(bases), enzyme) if site.overhang}


def _one_overhang(bases: str, enzyme: Enzyme, invariant: Invariant, where: str) -> str:
    """Return the one overhang `enzyme` leaves inside `bases`, refusing unless there is one.

    Raises
    ------
    ValueError
        Naming `invariant`, if the stuffer carries no such cut or more than one.
    """
    found = _overhangs(bases, enzyme)
    if not found:
        refuse(
            invariant,
            f"{where} is {bases!r}, in which {enzyme.name} leaves no overhang: a stuffer carries "
            "that enzyme's site and the bases its cut leaves single-stranded",
        )
    if len(found) > 1:
        refuse(
            invariant,
            f"{where} is {bases!r}, in which {enzyme.name} leaves {len(found)} different "
            f"overhangs, {', '.join(sorted(found))}: a part enters on one",
        )
    return found.pop()


@dataclass(frozen=True, slots=True)
class Scheme:
    """The architecture of one method, checked on construction.

    Every position carries the same stuffers, so what a position is is a name, which a project
    supplies.

    Parameters
    ----------
    name
        What the method is called, so a report can say what a design assumed.
    internal_enzyme
        The name of the enzyme that opens the library built so far by excising its internal
        stuffer. Its sites are the ones the product keeps, which is what lets the next round open
        it.
    external_enzyme
        The name of the enzyme that releases a part from its synthesised block. Its sites leave on
        the discarded external stuffers.
    blunt_enzymes
        The names of the blunt choppers that cut the pieces a round throws away, so that neither
        can ligate back. Any number.
    internal_stuffer_prefix
        The bases an internal stuffer opens with, before the shared core. It ends with the entry
        overhang, which is how a part carries its place.
    internal_stuffer_core
        The bases every internal stuffer shares, following the prefix.
    external_stuffer_5
        The 5' flank of the synthesised block, 5' to 3' as it is ordered, up to and including the
        bases the external enzyme's cut leaves single-stranded. Those bases are the entry
        overhang, which the method reads rather than states.
    external_stuffer_3
        The 3' flank, from the bases that cut leaves single-stranded — the cloning scar — to the
        end of the block.
    cloning_scar
        The bases every part's 3' end leaves at its junction.
    source
        Where these values came from, so provenance travels with the method.
    reserved
        The names of further enzymes a block must be free of, beyond the ones above. A step
        outside the rounds — seating a part in its carrier, or the last transfer into a working
        vector — cuts the cargo too, and its enzyme has no role here to be named by. A block
        spells none of these anywhere, its stuffers included. A project adds to this list and
        never replaces it.

    Raises
    ------
    ValueError
        If an invariant fails, which names the invariant: one of `Invariant`.
    KeyError
        If a named enzyme is not one this package ships.
    """

    name: str
    _: KW_ONLY
    internal_enzyme: str
    external_enzyme: str
    blunt_enzymes: tuple[str, ...]
    internal_stuffer_prefix: str
    internal_stuffer_core: str
    external_stuffer_5: str
    external_stuffer_3: str
    cloning_scar: str
    source: str = ""
    reserved: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Normalise the sequences, then check every invariant in turn."""
        object.__setattr__(self, "blunt_enzymes", tuple(self.blunt_enzymes))
        object.__setattr__(self, "reserved", tuple(self.reserved))
        for field_name in (
            "internal_stuffer_prefix",
            "internal_stuffer_core",
            "external_stuffer_5",
            "external_stuffer_3",
            "cloning_scar",
        ):
            object.__setattr__(self, field_name, getattr(self, field_name).upper())
        self._check_stuffers()
        self._check_terminal_frame()
        self._check_enzyme_regions()

    @property
    def internal(self) -> Enzyme:
        """The enzyme that opens the library built so far."""
        return get_enzyme(self.internal_enzyme)

    @property
    def external(self) -> Enzyme:
        """The enzyme that releases a part from its block."""
        return get_enzyme(self.external_enzyme)

    @property
    def blunt(self) -> tuple[Enzyme, ...]:
        """The blunt choppers, in the order this method names them."""
        return tuple(get_enzyme(name) for name in self.blunt_enzymes)

    @property
    def blunt_for_the_destination(self) -> tuple[Enzyme, ...]:
        """The choppers acting where the library is opened: those that shred its own stuffer.

        Read off the DNA, like every other overhang here: a chopper is in this tube where it
        reads the internal stuffer core, which is the piece that tube throws away.
        """
        return tuple(one for one in self.blunt if _cuts(self.internal_stuffer_core, one))

    @property
    def reserved_enzymes(self) -> tuple[Enzyme, ...]:
        """The enzymes a step outside the rounds reserves, in the order this method names them."""
        return tuple(get_enzyme(name) for name in self.reserved)

    @property
    def internal_stuffer(self) -> str:
        """The internal stuffer a part carries: its prefix and the shared core."""
        return self.internal_stuffer_prefix + self.internal_stuffer_core

    @property
    def entry_overhang(self) -> str:
        """The overhang admitting a part, read off the 5' external stuffer."""
        return _one_overhang(
            self.external_stuffer_5, self.external, "external-stuffer-5", "the 5' external stuffer"
        )

    @property
    def scar_overhang(self) -> str:
        """The overhang the 3' external stuffer yields, shared by every position."""
        return _one_overhang(
            self.external_stuffer_3, self.external, "external-stuffer-3", "the 3' external stuffer"
        )

    def barcode_block_length(self, barcode_length: int, positions: int) -> int:
        """How long a finished barcode block is: one barcode a position, joined by the scar."""
        return positions * barcode_length + (positions - 1) * len(self.cloning_scar)

    def retained_length(self, barcode_length: int, positions: int) -> int:
        """Return what the product keeps past its last part: terminal stuffer and barcode block.

        The frame of everything downstream rides on this, which is why the internal stuffer prefix
        is longer than an overhang where it has to be.
        """
        return len(self.internal_stuffer) + self.barcode_block_length(barcode_length, positions)

    def _check_stuffers(self) -> None:
        """Check a part enters on one overhang, leaves the stated scar, and that two cuts open it.

        Reading either overhang is itself the check that one external stuffer carries one cut.
        """
        entry, scar = self.entry_overhang, self.scar_overhang
        if scar != self.cloning_scar:
            refuse(
                "external-stuffer-3",
                f"the 3' external stuffer yields {scar!r}, which is not this method's cloning "
                f"scar {self.cloning_scar!r}",
            )
        if entry == scar:
            refuse(
                "internal-stuffer-cuts",
                f"a part enters on {entry!r}, which is also this method's cloning scar, so "
                f"{self.internal.name} leaves that overhang at both ends of the internal "
                "stuffer: excising it leaves two ends that anneal to each other, and a part goes "
                "in either way round",
            )
        if not self.internal_stuffer_prefix.endswith(entry):
            refuse(
                "internal-stuffer-prefix",
                f"the internal stuffer prefix is {self.internal_stuffer_prefix!r}, which does "
                f"not end with {entry!r}, the overhang a part enters on",
            )
        left = _overhangs(self.internal_stuffer, self.internal)
        if entry not in left:
            spells = ", ".join(sorted(left)) or "no overhang at all"
            refuse(
                "internal-stuffer-prefix",
                f"{self.internal.name} leaves {spells} in the internal stuffer, and not the "
                f"{entry!r} the next part enters on",
            )
        if scar not in left:
            spells = ", ".join(sorted(left)) or "no overhang at all"
            refuse(
                "internal-stuffer-cuts",
                f"{self.internal.name} leaves {spells} in the internal stuffer, and not the "
                f"cloning scar {scar!r}: a stuffer carries both of the cuts that open it",
            )

    def _check_terminal_frame(self) -> None:
        """Check what the product keeps past its last part is a whole number of codons.

        A barcode and the scar joining it to the last are whole codons together, which the
        project checks, so the barcode block is whole codons however many positions there are and
        whatever a barcode is long. What is left to check is the stuffer, less the scar the block
        does not repeat at its 5' end.
        """
        over = (len(self.internal_stuffer) - len(self.cloning_scar)) % 3
        if over:
            refuse(
                "terminal-frame",
                f"the internal stuffer is {len(self.internal_stuffer)} bases and the cloning scar "
                f"{len(self.cloning_scar)}, leaving the product {over} base(s) past a whole "
                "number of codons where it reads through what it keeps",
            )

    def _check_enzyme_regions(self) -> None:
        """Check no enzyme reads into a region another owns, and that every chopper chops."""
        internal, external = self.internal, self.external
        for end, bases in (("5'", self.external_stuffer_5), ("3'", self.external_stuffer_3)):
            if _cuts(bases, internal):
                refuse(
                    "enzyme-regions",
                    f"{internal.name} reads a site in the {end} external stuffer, which "
                    f"{external.name} owns: opening the library would cut the donor block as well",
                )
        if _cuts(self.internal_stuffer, external):
            refuse(
                "enzyme-regions",
                f"{external.name} reads a site in the internal stuffer, which {internal.name} "
                "owns: releasing the part would cut it apart",
            )
        discarded = (self.internal_stuffer_core, self.external_stuffer_5, self.external_stuffer_3)
        for chopper in self.blunt:
            if not any(_cuts(bases, chopper) for bases in discarded):
                refuse(
                    "enzyme-regions",
                    f"the blunt enzyme {chopper.name} reads no site in the internal stuffer core "
                    "or either external stuffer, so it chops none of the pieces a round discards",
                )
        kept = (self.internal_stuffer, self.external_stuffer_5, self.external_stuffer_3)
        for held in self.reserved_enzymes:
            if any(_cuts(bases, held) for bases in kept):
                refuse(
                    "enzyme-regions",
                    f"the reserved enzyme {held.name} reads a site in a stuffer of this method, "
                    "and a block spells a reserved site nowhere: the step that reserved it would "
                    "cut every part",
                )


#: The junction overhangs the method's shelved DNA carries, each pair as the two ends of one
#: molecule. Four overhangs over four pairs, every one of them one base closing a codon plus a
#: codon: ``T|ATG``, ``A|GGA``, ``T|TCC``, ``C|TAA``. The cargo pair is the one a library build
#: works on, and `IGGA` reads it off its own stuffers.
INTERFACE_OVERHANGS: Mapping[str, tuple[str, str]] = MappingProxyType(
    {
        "working cassette": ("TATG", "CTAA"),
        "N-terminal part": ("TATG", "AGGA"),
        "cargo": ("AGGA", "TTCC"),
        "C-terminal part": ("TTCC", "CTAA"),
    }
)

#: The method. Its invariants are checked here, once, when this module is imported.
IGGA = Scheme(
    "iGGA",
    internal_enzyme="BbsI",
    external_enzyme="BsaI",
    blunt_enzymes=("SrfI", "PmeI"),
    internal_stuffer_prefix="AGGA",
    internal_stuffer_core="AAGTCTTCAGCCCGGGCAGAAGACAATTCC",
    external_stuffer_5="GTTTAAACACATTCAGCGGGTCTCAAGGA",
    external_stuffer_3="TTCCTGAGACCCGCTGAATGTGTTTAAAC",
    cloning_scar="TTCC",
    reserved=("BsmBI",),
    source=(
        "docs/synthesis-and-assembly.md, the iGGA cargo and pipeline sections. The 34-base "
        "internal stuffer is the method page's own; the external stuffers lay BsaI and PmeI at "
        "the offsets their own cut offsets ask for, so every overhang is read off the DNA rather "
        "than stated."
    ),
)


#: The enzyme that cuts a cargo fragment out of its oligo, templated on the oligo rather than
#: carried by a primer: all of a gene's fragments share one primer pair, so a primer-borne site
#: would give every fragment of that gene the same overhang. It is one the method already
#: reserves, which is what keeps a designed block free of it.
SYNTHESIS_ENZYME = "BsmBI"

#: How the orthogonal primer set divides across its three roles, in the order the set is allotted
#: in. 96 inner primers makes one batch exactly one plate of PCR2, which is what fixes that
#: number; the other two index the batch together and so multiply.
ORTHOGONAL_SPLIT: tuple[tuple[str, int], ...] = (
    ("P2 gene reverse", 96),
    ("P1 batch forward", 35),
    ("P3 batch outer", 34),
)

#: The share of designs assembled from that many fragments with a perfect clone among four
#: colonies, as Lund et al. 2024 measured it. The only measured curve the fragment count is read
#: against; a count absent from it was not measured and nothing interpolates one.
LUND_SUCCESS: Mapping[int, float] = MappingProxyType(
    {2: 1.0, 3: 0.938, 5: 0.846, 8: 0.667, 12: 0.400, 16: 0.0}
)


def _check_interface() -> None:
    """Check the method reads the cargo pair off its own DNA.

    The pairs are what the ccdB cassettes and the DMX vector carry, and the stuffers are what a
    part is ordered with. They are the same molecules, so they agree or one of them is wrong.
    """
    carried = (IGGA.entry_overhang, IGGA.scar_overhang)
    if carried != INTERFACE_OVERHANGS["cargo"]:
        refuse(
            "interface-overhangs",
            f"the stuffers enter a part on {carried[0]} and scar it with {carried[1]}, where the "
            f"method's cargo junction is {INTERFACE_OVERHANGS['cargo']}",
        )


_check_interface()

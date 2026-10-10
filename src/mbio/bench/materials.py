"""What a material brings with it: the parameters that belong to the thing, not to the step.

A number attaches to the thing that changes it. An electroporation program belongs to the
cells: NEB #C3020 pulses at 2.0 kV, 200 Ω and 25 µF, Endura at 1800 V, 600 Ω and 10 µF, and
they disagree on every setting. So the program is keyed by the cells' catalogue number and
reached by a function, as `inactivation.heat_inactivation` reaches an enzyme's own. Swap the
cells and the number changes while the step does not.

A rule the material carries travels with it the same way. T7 DNA ligase refuses blunt ends
only because the PEG stays low, so the two rules below hang on the ligase and on its buffer
and reach every step that uses either. A caution rides the same key.

Every number here is sourced in ``docs/research/bench-numbers.md``.
"""

from collections.abc import Mapping
from dataclasses import KW_ONLY, dataclass
from types import MappingProxyType

from mbio.protocol.model import Citation, Material, Rule, Source

#: Every document the parameters below were read from, keyed as a `Citation` names it. A
#: protocol copies the ones it cites into its own `sources`.
SOURCES: Mapping[str, Source] = MappingProxyType(
    {
        "M0318": Source(
            "New England Biolabs #M0318 T7 DNA Ligase product page",
            url="https://www.neb.com/en/products/m0318-t7-dna-ligase",
            read_as="r.jina.ai",
            date="2026-10-06",
        ),
        "C3020": Source(
            "New England Biolabs #C3020 NEB 10-beta electrocompetent cells, electroporation "
            "protocol",
            url="https://www.neb.com/en/products/c3020-neb-10-beta-electrocompetent-e-coli",
            read_as="r.jina.ai",
            date="2026-10-06",
        ),
        "MA133": Source(
            "Endura Competent Cells manual, Lucigen / Biosearch Technologies",
            edition="MA133 26Feb2018",
            read_as="plain curl",
            date="2026-10-06",
        ),
    }
)

#: The buffer T7 DNA ligase ships with, and what it brings into the tube. 1x is 7.5% PEG 6000,
#: well under the 20% that would force blunt ligation. M0318 product page, reaction conditions.
STICKTOGETHER_PEG_PERCENT = 7.5

#: The PEG 6000 concentration at which T7 DNA ligase gains measurable blunt-end activity.
#: M0318 product page.
BLUNT_FORCING_PEG_PERCENT = 20.0


@dataclass(frozen=True, slots=True)
class Electroporation:
    """The pulse one strain of electrocompetent cells is given, and the cuvette it is given in.

    Parameters
    ----------
    cells
        What the catalogue number is.
    volts, ohms, microfarads
        The program. These belong to the cells; no two strains here agree on them.
    cuvette_mm
        The gap of the cuvette, chilled.
    cells_ul
        One transformation's worth.
    time_constant_ms
        What the supplier says to expect, low to high.
    recovery
        The medium and the hold the supplier gives, in its own words.
    citation
        Where it was read.
    """

    cells: str
    volts: float
    ohms: float
    microfarads: float
    _: KW_ONLY
    cuvette_mm: float
    cells_ul: float
    time_constant_ms: tuple[float, float]
    recovery: str
    citation: Citation


#: Each strain's own program, keyed by catalogue number.
ELECTROPORATION: Mapping[str, Electroporation] = MappingProxyType(
    {
        "C3020": Electroporation(
            "NEB 10-beta electrocompetent E. coli",
            2000.0,
            200.0,
            25.0,
            cuvette_mm=1.0,
            cells_ul=25.0,
            time_constant_ms=(4.8, 5.1),
            recovery="975 µL of 37 °C NEB 10-beta/Stable Outgrowth Medium immediately, then "
            "1 hour at 37 °C and 250 rpm",
            citation=Citation("C3020", "electroporation protocol"),
        ),
        "60242": Electroporation(
            "Endura electrocompetent cells",
            1800.0,
            600.0,
            10.0,
            # The manual writes the gap in centimetres; this field is millimetres.
            cuvette_mm=1.0,
            cells_ul=25.0,
            time_constant_ms=(3.5, 4.5),
            recovery="975 µL of the supplied Recovery Medium within 10 seconds of the pulse, "
            "then 1 hour at 37 °C and 250 rpm",
            citation=Citation("MA133", "p. 4-5"),
        ),
    }
)

#: What each material brings into the tube besides itself, keyed by catalogue number.
CONTAINS: Mapping[str, tuple[str, ...]] = MappingProxyType({"B0535": ("PEG 6000",)})

#: What must not, or must, happen where a material is used, keyed by catalogue number. Both T7
#: ligase rules are here rather than in a step's prose, because prose an agent edits can be
#: deleted and the method's escape removal rests on the first of them.
RULES: Mapping[str, tuple[Rule, ...]] = MappingProxyType(
    {
        "M0318": (
            Rule(
                "forbids",
                "PEG",
                "Never add PEG to a T7 ligase reaction. Blunt-end ligation is forced above "
                f"{BLUNT_FORCING_PEG_PERCENT:g}% (w/v) PEG 6000, and StickTogether buffer "
                f"already carries {STICKTOGETHER_PEG_PERCENT:g}%; the refusal that removes "
                "escapees holds only while nothing adds more.",
                citation=Citation("M0318", "reaction conditions"),
            ),
            Rule(
                "forbids",
                "heat inactivation",
                "Do not heat-inactivate T7 DNA ligase in PEG: transformation will be "
                "inhibited. Go from the ligation straight to the clean-up.",
                when="PEG",
                citation=Citation("M0318", "note 4"),
            ),
        )
    }
)

#: What to watch out for wherever a polymerase is pipetted, and wherever cells are pulsed.
POLYMERASE_ON_ICE = "Keep the polymerase on ice."
CUVETTE_ON_ICE = "Keep the cells and the cuvette on ice; a warm cuvette arcs."

#: What to watch out for where a material is used, keyed by catalogue number. ADR 0020 puts it
#: here: the same sentence written into each step that needs it is as many places to fix, and
#: an agent editing the protocol JSON can delete any of them.
CAUTIONS: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {
        "M0491": (POLYMERASE_ON_ICE,),
        "M0267": (POLYMERASE_ON_ICE,),
        "60242": (CUVETTE_ON_ICE,),
        "C3020": (CUVETTE_ON_ICE,),
    }
)

#: The spin-column kits a run may name, keyed by catalogue number: the product's own name and
#: who sells it. Each is one a note under `docs/research/` already reads a number off, so
#: nothing here is a product this package met for the first time. A kit the table does not hold
#: is named as the caller wrote it and carries no number, since nothing invents one.
KITS: Mapping[str, tuple[str, str]] = MappingProxyType(
    {
        # docs/research/restriction-ligation.md, from the two Monarch manuals.
        "T1120": ("Monarch Spin DNA Gel Extraction Kit", "New England Biolabs"),
        "T1130": ("Monarch Spin PCR & DNA Cleanup Kit", "New England Biolabs"),
        # docs/research/bench-numbers.md, from Qian's supplementary protocol.
        "D4003": ("DNA Clean & Concentrator-5", "Zymo Research"),
        "D4007": ("Zymoclean Gel DNA Recovery Kit", "Zymo Research"),
    }
)

#: The kit a step names where the run names none. It is the one this package already reads its
#: column recovery from, so the page's numbers and its materials row name the same product. A
#: lab that uses another one names it per run; no lab's habit is written in here.
DEFAULT_CLEANUP_KIT = "T1130"

#: Every catalogue number these parameters are keyed by. A number is keyed without its pack
#: size, because a pack size changes nothing about the thing in the tube.
_KEYED = ELECTROPORATION.keys() | RULES.keys() | CONTAINS.keys() | CAUTIONS.keys() | KITS.keys()


def _key(catalog: str) -> str:
    """Return the catalogue number a parameter is keyed by, without its pack size."""
    stem = catalog.strip().lstrip("#").upper().split("-")[0]
    return stem if stem in _KEYED else stem.rstrip("SLIH")


def electroporation(catalog: str) -> Electroporation | None:
    """Return the electroporation program these cells are pulsed with, or ``None`` for unknown.

    Examples
    --------
    >>> one = electroporation("#C3020")
    >>> one.volts, one.ohms, one.microfarads
    (2000.0, 200.0, 25.0)
    >>> electroporation("60242-2").volts
    1800.0
    """
    return ELECTROPORATION.get(_key(catalog))


def rules(catalog: str) -> tuple[Rule, ...]:
    """Return what the material with this catalogue number forbids or requires.

    Examples
    --------
    >>> [rule.subject for rule in rules("M0318L")]
    ['PEG', 'heat inactivation']
    """
    return RULES.get(_key(catalog), ())


def contains(catalog: str) -> tuple[str, ...]:
    """Return what this material brings into the tube besides itself."""
    return CONTAINS.get(_key(catalog), ())


def cautions(catalog: str) -> tuple[str, ...]:
    """Return what to watch out for wherever the material with this number is used.

    Examples
    --------
    >>> cautions("M0491S")
    ('Keep the polymerase on ice.',)
    """
    return CAUTIONS.get(_key(catalog), ())


def kit(named: str = "", *, note: str = "") -> Material:
    """Return the kit a run named, as a material carrying who sells it and its number.

    Parameters
    ----------
    named
        A catalogue number `KITS` holds, or the product's own name. Empty names
        `DEFAULT_CLEANUP_KIT`.
    note
        What it is there for.

    Examples
    --------
    >>> kit().name, kit("D4003").supplier
    ('Monarch Spin PCR & DNA Cleanup Kit', 'Zymo Research')
    >>> named = kit("Wizard SV Gel and PCR Clean-Up System")
    >>> named.name, named.supplier, named.catalog
    ('Wizard SV Gel and PCR Clean-Up System', '', '')
    """
    text = (named or DEFAULT_CLEANUP_KIT).strip()
    found = KITS.get(_key(text))
    if found is None:
        return material(text, note=note)
    sold, supplier = found
    return material(sold, supplier=supplier, catalog=text.lstrip("#").upper(), note=note)


def material(
    name: str,
    *,
    supplier: str = "",
    catalog: str = "",
    storage: str = "",
    amount: str = "",
    note: str = "",
    citation: Citation | None = None,
) -> Material:
    """Return a material carrying its own parameters: what it brings, rules and cautions.

    The caller names the thing and what travels with it is looked up, so a material built here
    cannot reach a protocol without them. One built around this carries none of them.

    Examples
    --------
    >>> ligase = material("T7 DNA Ligase", supplier="NEB", catalog="#M0318L")
    >>> len(ligase.rules), material("StickTogether", catalog="#B0535S").contains
    (2, ('PEG 6000',))
    >>> material("Q5 DNA Polymerase", catalog="M0491").cautions
    ('Keep the polymerase on ice.',)
    """
    return Material(
        name,
        supplier=supplier,
        catalog=catalog,
        storage=storage,
        amount=amount,
        note=note,
        contains=contains(catalog),
        rules=rules(catalog),
        cautions=cautions(catalog),
        citation=citation,
    )

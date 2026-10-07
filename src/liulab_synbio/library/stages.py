"""What the iGGA round asks of the bench model: its plates, its materials' rules and its holes.

Every number here is quoted from ``docs/research/bench-numbers.md``, which names the document
each came from. Nothing is invented: where the method needs a number nobody published, a `Hole`
stands in its place and says which ticket it is routed to.

This is the pooled experiment's side of the bench. The validating experiment's plates, its
barcode kit, its depth floors and its per-well pass rule are `liulab_synbio.dmx`'s. A round
makes a pool nobody picks from, so neither side borrows the other's gate.

The plates and the rules are `liulab_mbio.bench`'s. This module only says what each stage holds
and which materials carry which rule — a different method would call the same functions with
its own numbers.
"""

from liulab_mbio.bench import materials
from liulab_mbio.protocol.model import Citation, Hole, Material, Source, Vessel

#: The documents a citation in this protocol resolves against: the ones a material brings with
#: it. The round's own numbers are the paper's and travel as references, not as cited rows.
SOURCES: dict[str, Source] = dict(materials.SOURCES)

#: The ligase this method's round runs on, and the buffer it runs in. Both carry their own
#: rules, so neither can be used in a protocol that does not show them: never add PEG to a T7
#: ligase reaction, and do not heat-inactivate it in PEG.
LIGASE = materials.material(
    "T7 DNA Ligase",
    supplier="NEB",
    catalog="#M0318L",
    storage="-20 °C",
    note="3,000,000 U/mL. Refuses blunt ends, which is what removes escapees.",
    citation=Citation("M0318", "reaction conditions"),
)
LIGASE_BUFFER = materials.material(
    "StickTogether DNA Ligase Buffer",
    supplier="NEB",
    catalog="#B0535S",
    storage="-20 °C",
    note="Supplied 2x. 1x is 7.5% PEG 6000, well under the 20% that would force blunt ligation.",
    citation=Citation("M0318", "reaction conditions"),
)

#: What this method cannot write completely. Each is a number nobody published, routed to the
#: ticket that would decide it; the ids are the research note's own, so a reference still
#: resolves.
HOLES: tuple[Hole, ...] = (
    Hole(
        "H22",
        "no selection antibiotic for a round",
        "undecided",
        where="each round, plating",
        filled_by="the vector's own marker",
        issue="liuhlab/liulab-mbio#264",
    ),
    Hole(
        "H23",
        "no published document describes the split digest; the one source is the paper itself",
        "unpublished",
        where="the split digest, units per reaction",
        filled_by="nothing; it is the method's own",
        issue="liuhlab/liulab-mbio#264",
    ),
    Hole(
        "H24",
        "every mass in the final assembly is unsourced, the SPRI ratio with them",
        "undecided",
        where="final assembly, DNA in",
        filled_by="a pilot, or carrying the round's own numbers across and saying so",
        issue="liuhlab/liulab-mbio#264",
    ),
)


def titre_plates(number: int) -> tuple[Vessel, Vessel]:
    """Return the two plates one round is bounded by: a dilution, and the no-donor control.

    The control is the cut destination carried through the ligation with no donor added, so
    what grows on it is the parental destination the design rests on not surviving. Subtracting
    it from the dilution is what the round is judged on. Both grow during the outgrowth, so
    they cost no day. The source plates nothing at any round; ``D18`` of
    ``docs/research/synthesis-and-assembly-departures.md`` says why this method does.

    Raises
    ------
    ValueError
        If the round does not count from one.

    Examples
    --------
    >>> dilution, control = titre_plates(2)
    >>> dilution.name, control.name
    ('round 2 titre', 'round 2 no-donor control')
    """
    if number < 1:
        raise ValueError(f"rounds count from one, got {number}")
    selection = "the destination vector's own marker, which this plan does not name"
    return (
        Vessel(
            f"round {number} titre",
            kind="LB agar plate",
            holds="a measured dilution of the recovery",
            note=f"grown during the outgrowth; {selection}",
        ),
        Vessel(
            f"round {number} no-donor control",
            kind="LB agar plate",
            holds="the same digest taken through the ligation with no donor added",
            note=f"its colonies are subtracted from the titre; {selection}",
        ),
    )


def method_materials() -> tuple[Material, ...]:
    """Return the materials whose own parameters and rules this method depends on."""
    return (LIGASE, LIGASE_BUFFER)

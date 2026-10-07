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

#: The documents a citation in this protocol could resolve against: the ones a material brings
#: with it. The round's own numbers are the paper's and travel as references, not as cited rows.
#: Copied whole rather than picked over, because which of them a run cites depends on the steps
#: it builds; `liulab_mbio.protocol.citing` drops the rest before the protocol is returned.
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

#: What no source sets for the final assembly: no published document runs the step, so the DNA
#: into it and the ratio it is cleaned up at are the method's own. Carried by the steps it bites
#: in, and listed among the protocol's holes.
FINAL_MASSES = Hole(
    "H24",
    "every mass in the final assembly is unsourced, the SPRI ratio with them",
    "undecided",
    where="final assembly, DNA in",
    filled_by="a pilot, or carrying the round's own numbers across and saying so",
    issue="liuhlab/liulab-mbio#264",
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
    FINAL_MASSES,
)


#: What the pool route cannot write, and would have to invent a number to fill. A block is
#: split for synthesis and nothing says what puts it back together: the method assembles cargo
#: into the DMX vector, and a block already carries the flanks that vector would supply, so the
#: pieces have no destination to close into and no reaction sized against one.
POOL_HOLES: tuple[Hole, ...] = (
    Hole(
        "H25",
        "nothing names what the assembled pieces close into, or what selects the closure",
        "undecided",
        where="block assembly, the destination",
        filled_by="splitting the cargo rather than the whole block, so the pieces close into "
        "the DMX backbone the parent's own BsmBI digest leaves",
        issue="liuhlab/liulab-mbio#373",
    ),
    Hole(
        "H26",
        "no source sets the DNA in, the enzyme, the ligase or the cycling for a block assembly",
        "undecided",
        where="block assembly, the reaction and its program",
        filled_by="NEB's own Golden Gate table, once the assembly has a destination to be "
        "sized against",
        issue="liuhlab/liulab-mbio#373",
    ),
)


#: What neither read of the finished library can write. The method says what each one spans and
#: that it is a long read, and names no primer against it: #343 is open on that, and the two
#: reads hit the same wall. Carried by the linkage read, where the amplicon is longest.
READ_PRIMERS = Hole(
    "H27",
    "no primer pair, amplicon length or read depth for either library read",
    "undecided",
    where="the linkage and representation reads",
    filled_by="whichever way #343 is decided: the plan designs both pairs against the finished "
    "record, or the step says what each must span and leaves the pair to the reader",
    issue="liuhlab/liulab-mbio#343",
)

#: What neither read can be judged by. The source's own figures are what it reached, not a mark
#: it set, and what this library must reach follows from the screen downstream, which the method
#: cannot see. Carried by the representation read, which is the one repeated.
READ_PASS_MARK = Hole(
    "H28",
    "no pass mark for either library read: what share of the library must be seen, how even the "
    "counts must be, or what linkage fidelity passes",
    "undecided",
    where="the linkage and representation reads, what carries the library forward",
    filled_by="the screen downstream, which is what sets the representation the library has to "
    "reach",
)

#: What a build that names no working vector cannot say. The working vector is the user's own
#: stock, chosen per application, and the cargo enzyme, the assembly's cycling and the final
#: vector's length all follow from it. Carried by the step that picks it.
WORKING_VECTOR = Hole(
    "H31",
    "no working vector is named, so nothing names the cargo enzyme, its cycling or the final "
    "vector's length",
    "lab",
    where="final assembly, the vector the library moves into",
    filled_by="the project naming a working vector, which is a stock the lab holds and an "
    "application chooses",
)

#: What a library whose backbone does not present the cargo cannot say. The cargo is released by
#: the enzyme the method releases a part with, and the sites that free it belong to the vector the
#: rounds ran in, not to the cargo. Carried by the release step.
CARGO_RELEASE = Hole(
    "H32",
    "the backbone this library was built in presents no cut that frees its cargo, so no release "
    "digest can be written",
    "undecided",
    where="final assembly, releasing the cargo",
    filled_by="a destination vector carrying the releasing sites outboard of the cargo, as the "
    "method's own DMX vector does",
    issue="liuhlab/liulab-mbio#225",
)

#: What Twist's banded cycle count does not reach. It is stated against KAPA HiFi HotStart or
#: Twist's own TrueAmp mix, and the same guide's FAQ answers that another polymerase may give
#: worse uniformity at the same count. This method runs Q5, which nobody has published a count
#: for. Carried by PCR1, which prints the band's fewest.
#: ``docs/research/oligo-pool-pcr-cycles.md`` section 2.
PCR1_POLYMERASE = Hole(
    "H29",
    "no cycle count published for Q5 on a Twist oligo pool",
    "unpublished",
    where="PCR1, the polymerase the cycle count is stated against",
    filled_by="a count Twist states against Q5, or a pilot on this pool",
    issue="liuhlab/liulab-mbio#225",
)

#: What no source covers at all. Twist's table amplifies the pool as it arrives; PCR2's template
#: is PCR1's product, and the published subpool counts are first amplifications off a pool. So
#: PCR2 prints the stopping rule and no number.
#: ``docs/research/oligo-pool-pcr-cycles.md`` sections 4 and 6.
PCR2_CYCLES = Hole(
    "H30",
    "no cycle count for PCR2, which pulls one block out of an already-amplified batch",
    "unpublished",
    where="PCR2, the cycle count",
    filled_by="a pilot titrated against the heteroduplex hump on capillary electrophoresis, or a "
    "real-time run stopped before the curve plateaus",
    issue="liuhlab/liulab-mbio#225",
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

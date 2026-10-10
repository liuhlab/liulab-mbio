"""What the iGGA round asks of the bench model: its plates, its materials' rules and its holes.

Every number here is quoted from ``docs/research/bench-numbers.md``, which names the document
each came from. Nothing is invented: where the method needs a number nobody published, a `Hole`
stands in its place.

This is the pooled experiment's side of the bench. The validating experiment's plates, its
barcode kit, its depth floors and its per-well pass rule are `synbio.dmx`'s. A round
makes a pool nobody picks from, so neither side borrows the other's gate.

The plates and the rules are `mbio.bench`'s. This module only says what each stage holds
and which materials carry which rule — a different method would call the same functions with
its own numbers.
"""

from mbio.bench import materials
from mbio.bench.phenotype import SELECTION, SELECTION_PLATE, selection_marker
from mbio.protocol.model import Citation, Hole, Material, Source, Vessel
from mbio.sequence import SequenceRecord
from synbio.igga.bench import DIGEST_SOURCE, DIGEST_SOURCE_KEY

#: The documents a citation in this protocol could resolve against: the ones a material brings
#: with it, and the specifications a digest's unit count is read from. The round's own numbers
#: are the paper's, and the rows resting on them cite it.
#: Copied whole rather than picked over, because which of them a run cites depends on the steps
#: it builds; `mbio.protocol.citing` drops the rest before the protocol is returned.
SOURCES: dict[str, Source] = dict(materials.SOURCES) | {DIGEST_SOURCE_KEY: DIGEST_SOURCE}

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


def selection_for(record: SequenceRecord) -> str:
    """Return what to select this vector's transformants on, or ``""`` where nothing names it.

    The marker is a fact in the record, so the drug follows from it rather than from the
    document the rest of the stage was quoted from: departure D11 rebuilds this method's DMX
    vector AmpR to KanR, and the paper's carbenicillin does not carry over with it.
    """
    marker = selection_marker(record)
    if marker is None:
        return ""
    drug = SELECTION.get(marker.name.lower(), "")
    return SELECTION_PLATE.get(drug, drug)


#: What a round has no drug to plate on, raised only where `selection_for` cannot name one.
ROUND_SELECTION = Hole(
    "H22",
    "no selection antibiotic for a round",
    "undecided",
    where="each round, plating",
    filled_by="the vector's own marker",
)


#: What no source sizes in the final assembly. The paper runs its own last transfer by the
#: methods it gives for a round, which sizes the release and the clean-up around this step. It
#: does not size the step itself, where the cargo arrives unpurified in the release volume and
#: no published reaction holds it. The cargo's own mass is not part of the hole: the release
#: frees it, so it computes.
FINAL_MASSES = Hole(
    "H24",
    "no mass for the one-pot assembly: the cargo arrives unpurified in the whole release, and "
    "the working vector's mass against it is not sourced",
    "undecided",
    where="final assembly, the one-pot reaction",
    filled_by="a pilot",
)


#: What this method cannot write completely. The ids are the research note's own, so a reference
#: still resolves. `FINAL_MASSES` is not here: a build may size that reaction, so the step that
#: would print the mass is the one place that knows whether the hole stands.
HOLES: tuple[Hole, ...] = (
    Hole(
        "H23",
        "nobody has published two enzymes run in one tube, the second added to the first's "
        "reaction with no heat kill before the beads; the one source is the paper itself",
        "unpublished",
        where="the split digest, the procedure",
        filled_by="nothing; it is the method's own",
    ),
)


def holes_for(record: SequenceRecord) -> tuple[Hole, ...]:
    """Return the holes a run against this vector carries.

    `ROUND_SELECTION` joins them only where the record names no marker to plate on: a vector
    that names one answers the hole rather than carrying it.
    """
    return HOLES if selection_for(record) else (ROUND_SELECTION, *HOLES)


#: What the linkage read alone cannot be judged by. The representation read is held to Joung's
#: pooled-library bar, in `mbio.bench.coverage`; nothing published says what share of
#: reads must carry a barcode that still names its part, and the source's own figure is what one
#: library reached rather than a mark it set.
READ_PASS_MARK = Hole(
    "H28",
    "no pass mark for the linkage read: what share of reads must carry a barcode that still "
    "names its part",
    "undecided",
    where="the linkage read",
    filled_by="this run, or a source that sets a mark rather than reporting what it reached",
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
    filled_by="this run naming a working vector, which is a stock the lab holds and an "
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

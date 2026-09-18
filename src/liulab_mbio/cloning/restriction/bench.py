"""What a restriction and ligation cloning takes at the bench: its digest and its ligation.

Every number below is NEB's, through ``docs/research/restriction-ligation.md``, and each says
which section it was read from. Functions return `liulab_mbio.protocol` values, so a protocol
prints them unchanged. What any cloning pipeline shares -- DNA amounts, gels, validation and the
heat inactivation off an enzyme's own record -- is `liulab_mbio.bench`.

Three numbers the note states as a range rather than a value, so the pick is named as a pick:
the vector into a ligation, the insert-to-vector ratio, and what a gel extraction recovers. The
ligation table NEB prints is 50 ng of a 4 kb vector against 37.5 ng of a 1 kb insert, which is
an example of the rule and not the rule; the rule is picomoles, and that is what is computed.
"""

from collections.abc import Sequence

from liulab_mbio.bench.amounts import Amount, dna_amount, to_pmol
from liulab_mbio.bench.reactions import fits, reaction_table
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.protocol.model import (
    Component,
    Incubation,
    ReactionTable,
    Reference,
    Stage,
    ThermocyclerProgram,
)

# --------------------------------------------------------------------------------------
# The digest
# --------------------------------------------------------------------------------------

#: NEB's "typical" digestion: DNA in micrograms, the volume it is cut in and the time it is cut
#: for. One unit cuts 1 µg in 50 µL in 60 minutes, and the table asks for a 5-10 fold
#: overdigestion. §2.
DIGEST_MICROGRAMS = 1.0
DIGEST_VOLUME_UL = 50.0
DIGEST_SECONDS = 3600

#: Each enzyme in that table: units, and the volume they generally come in, µL. §2.
DIGEST_UNITS = 10
ENZYME_UL = 1.0

#: The 10X buffer in it, µL. §2.
BUFFER_UL = 5.0

#: What a digest's buffer line is called where nothing sourced says the enzymes share one.
UNNAMED_BUFFER = "restriction enzyme buffer"

#: Where to look up a pair nothing here can judge. §3, §4.
BUFFER_FINDER = "NEBcloner or NEB's Double Digest Finder"

#: What the enzymes' own volume may not pass, as a fraction of the reaction. More than this is
#: over 5% glycerol, which is one of the conditions that contributes to star activity. §2, §3.
MAX_ENZYME_FRACTION = 0.1

#: What the DNA solution may not pass, as a fraction of the reaction: a spin-column eluate
#: carries salt, and salt inhibits the digest. Technical Guide, §2.
MAX_DNA_FRACTION = 0.25

#: The conditions NEB names as contributing to star activity, each with its own countermeasure.
#: Printed, never predicted: NEB's own caveat is that their weight varies from enzyme to enzyme.
#: §2.
STAR_ACTIVITY: tuple[str, ...] = (
    "Keep the enzymes under 10% of the reaction, which keeps glycerol under 5%.",
    "Use the buffer the supplier supplies each enzyme in wherever you can.",
    "Use the fewest units and the shortest incubation that digest completely.",
    "Keep the DNA free of organic solvents, and use Mg²⁺ rather than another metal.",
)


def digest_amount(
    dna: tuple[str, int],
    *,
    micrograms: float = DIGEST_MICROGRAMS,
    concentration_ng_ul: float | None = None,
) -> Amount:
    """Return what one digest takes, as picomoles of the DNA's own length.

    `dna` is the name and the length in base pairs of the plasmid being cut. Without
    `concentration_ng_ul` the volume is a placeholder, so nothing can be shown not to fit.

    Raises
    ------
    ValueError
        If the length, the concentration or `micrograms` is not positive.

    Examples
    --------
    >>> digest_amount(("pUC19", 2686)).nanograms
    1000.0
    """
    if micrograms <= 0:
        raise ValueError(f"digest of {dna[0]!r}: micrograms must be positive, got {micrograms}")
    name, length_bp = dna
    nanograms = 1000.0 * micrograms
    return dna_amount(
        name,
        length_bp,
        pmol=to_pmol(nanograms, length_bp),
        concentration_ng_ul=concentration_ng_ul,
    )


def shared_buffer(enzymes: Sequence[Enzyme]) -> str | None:
    """Return the buffer every one of these enzymes is supplied in, or ``None`` for no verdict.

    Two products supplied in the same buffer digest together in it, which is NEB's own rule.
    Two supplied in different ones carry **no verdict** rather than a pass: judging that pair
    takes the activity percentages, which no shippable source states for the enzymes here. The
    comparison is per product and not per enzyme name -- BamHI and BamHI-HF differ. §3, §4.

    Examples
    --------
    >>> from liulab_mbio.enzymes import get_enzyme
    >>> shared_buffer([get_enzyme("EcoRI"), get_enzyme("BamHI")])
    'rCutSmart Buffer'
    """
    supplied = {enzyme.supplied_buffer for enzyme in enzymes}
    return supplied.pop() if len(supplied) == 1 and None not in supplied else None


def digest_reaction(
    enzymes: Sequence[Enzyme],
    amount: Amount,
    *,
    volume_ul: float = DIGEST_VOLUME_UL,
    title: str = "",
    reactions: int = 1,
) -> ReactionTable:
    """Return the digest that cuts one plasmid with these enzymes in one tube.

    The buffer line names the buffer the enzymes share, and says nothing more than "restriction
    enzyme buffer" where `shared_buffer` has no verdict.

    Raises
    ------
    ValueError
        If no enzyme is given, if the enzymes take more than `MAX_ENZYME_FRACTION` of the
        reaction, or if the DNA does not fit it.
    """
    if not enzymes:
        raise ValueError("a digest needs at least one enzyme")
    taken = ENZYME_UL * len(enzymes)
    if taken > MAX_ENZYME_FRACTION * volume_ul:
        raise ValueError(
            f"{len(enzymes)} enzymes take {taken:g} µL of a {volume_ul:g} µL reaction, over the "
            f"{MAX_ENZYME_FRACTION:.0%} that keeps glycerol under 5%; scale the reaction up"
        )
    components = [
        Component(
            f"10X {shared_buffer(enzymes) or UNNAMED_BUFFER}",
            BUFFER_UL,
            stock="10X",
            final="1X",
        ),
        *(
            Component(enzyme.supplier_label, ENZYME_UL, final=f"{DIGEST_UNITS} units")
            for enzyme in enzymes
        ),
    ]
    return reaction_table(
        (amount,), components, volume_ul=volume_ul, title=title, reactions=reactions
    )


# --------------------------------------------------------------------------------------
# The ligation
# --------------------------------------------------------------------------------------

#: The ligation NEB prints for T4 DNA Ligase (NEB #M0202): the reaction, its 10X buffer and the
#: ligase, µL. §7.
LIGATION_VOLUME_UL = 20.0
LIGATION_BUFFER_UL = 2.0
LIGASE_UL = 1.0

#: T4 DNA Ligase as NEB sells it, units per µL. §7; the concentrate is five times this.
LIGASE_UNITS_UL = 400.0

#: Vector into one ligation, femtomoles. NEB recommends 20-30 fmol (§7); the low end is this
#: package's pick inside that range.
VECTOR_FMOL = 20.0

#: Insert to vector molar ratio. NEB calls 1:1 to 1:10 optimal for a single insertion (§7); 3 is
#: this package's pick inside that range, and the ratio NEB's own worked example uses.
INSERT_RATIO = 3.0

#: What the vector and insert together should reach, ng/µL: below the floor a fragment closes on
#: itself instead of joining its partner. §7, printed rather than computed -- a plan does not
#: know the concentrations the bench will measure.
LIGATION_NG_UL = (1.0, 10.0)

#: The incubations NEB's own table gives at room temperature, seconds: cohesive ends, and blunt
#: ends or a single-base overhang. The twelvefold difference is the only cost a supplier states
#: for a blunt ligation. §7, §8.
ROOM_CELSIUS = 25.0
COHESIVE_SECONDS = 600
BLUNT_SECONDS = 7200

#: The other half of the same table: either kind, overnight at this temperature. §7.
OVERNIGHT_CELSIUS = 16.0

#: How the ligase is stopped before the transformation. §7.
LIGASE_KILL_CELSIUS = 65.0
LIGASE_KILL_SECONDS = 600

#: What of the finished ligation goes into 50 µL of competent cells, µL. §7.
TRANSFORM_UL = (1.0, 5.0)


def ligation_amounts(
    backbone: tuple[str, int],
    insert: tuple[str, int],
    *,
    ratio: float = INSERT_RATIO,
    femtomoles: float = VECTOR_FMOL,
    volume_ul: float = LIGATION_VOLUME_UL,
    backbone_ng_ul: float | None = None,
    insert_ng_ul: float | None = None,
) -> tuple[Amount, Amount]:
    """Return what one ligation takes, the backbone first.

    `backbone` and `insert` are each a name and a length in base pairs. The backbone gets
    `femtomoles` of itself and the insert `ratio` times as many, weighed at its own length: a
    short insert weighs less for the same number of molecules, so matching the two masses would
    not match their molecules.

    Raises
    ------
    ValueError
        If a length, a concentration or `ratio` is not positive, or if the DNA does not fit
        `volume_ul` beside the buffer and the ligase.

    Examples
    --------
    >>> opened, released = ligation_amounts(("pUC19 backbone", 2665), ("GFP", 723))
    >>> opened.nanograms, released.nanograms
    (32.83, 26.72)
    """
    if ratio <= 0:
        raise ValueError(f"molar ratio must be positive, got {ratio}")
    if femtomoles <= 0:
        raise ValueError(f"vector femtomoles must be positive, got {femtomoles}")
    pmol = femtomoles / 1000.0
    opened = dna_amount(*backbone, pmol=pmol, concentration_ng_ul=backbone_ng_ul)
    released = dna_amount(*insert, pmol=pmol * ratio, concentration_ng_ul=insert_ng_ul)
    fits(
        (opened, released),
        volume_ul=volume_ul,
        taken_ul=LIGATION_BUFFER_UL + LIGASE_UL,
        what=f"ligation of {insert[0]} into {backbone[0]}",
    )
    return opened, released


def ligation_reaction(
    amounts: Sequence[Amount], *, volume_ul: float = LIGATION_VOLUME_UL, reactions: int = 1
) -> ReactionTable:
    """Return the ligation reaction for these fragments, the backbone first.

    Raises
    ------
    ValueError
        If fewer than two fragments are given, or if the DNA does not fit the reaction.
    """
    if len(amounts) < 2:
        raise ValueError("a ligation joins at least two fragments")
    components = [
        Component("T4 DNA Ligase Buffer", LIGATION_BUFFER_UL, stock="10X", final="1X"),
        Component(
            "T4 DNA Ligase",
            LIGASE_UL,
            stock=f"{LIGASE_UNITS_UL:g} U/µL",
            final=f"{LIGASE_UNITS_UL * LIGASE_UL:g} units",
        ),
    ]
    return reaction_table(
        amounts,
        components,
        volume_ul=volume_ul,
        title="Ligation, T4 DNA Ligase (NEB #M0202)",
        reactions=reactions,
    )


def ligation_program(*, blunt: bool = False) -> ThermocyclerProgram:
    """Return the ligation at room temperature and the heat step that stops it.

    `blunt` picks NEB's longer incubation, which its table gives to a blunt end or a single-base
    overhang. Either kind may run overnight at `OVERNIGHT_CELSIUS` instead.
    """
    return ThermocyclerProgram(
        (
            Stage(
                (Incubation("Ligate", ROOM_CELSIUS, BLUNT_SECONDS if blunt else COHESIVE_SECONDS),)
            ),
            Stage((Incubation("Heat inactivation", LIGASE_KILL_CELSIUS, LIGASE_KILL_SECONDS),)),
        ),
        title="Ligation",
    )


# --------------------------------------------------------------------------------------
# Getting from the digest to the ligation
# --------------------------------------------------------------------------------------

#: What a gel extraction recovers, as a fraction: below `LARGE_FRAGMENT_BP`, and at or above it.
#: Monarch Spin DNA Gel Extraction Kit (NEB #T1120), §9.
GEL_RECOVERY = (0.70, 0.90)
LARGE_GEL_RECOVERY = (0.50, 0.70)
LARGE_FRAGMENT_BP = 10000


def gel_recovery(length_bp: int) -> tuple[float, float]:
    """Return what a gel extraction recovers of a fragment this long, as a fraction.

    Raises
    ------
    ValueError
        If the length is not positive.

    Examples
    --------
    >>> gel_recovery(723)
    (0.7, 0.9)
    """
    if length_bp <= 0:
        raise ValueError(f"a fragment is at least one base pair, got {length_bp}")
    return LARGE_GEL_RECOVERY if length_bp >= LARGE_FRAGMENT_BP else GEL_RECOVERY


#: NEB's four transformation controls and what each should give relative to the others. No
#: supplier states an absolute colony count for this method, so what a plan promises is a ratio.
#: §12.
CONTROLS: tuple[str, ...] = (
    "Uncut vector, 100 pg to 1 ng: the cells are viable and the antibiotic is right.",
    "Cut vector, no ligase: under 1% of the colonies the uncut vector gave.",
    "Vector-only ligation: the same as the cut-vector control, its ends being unable to rejoin.",
    "A transformation efficiency under 10⁴ cfu/µg means the cells, not the ligation.",
)


#: Where the numbers above come from, ready for a protocol's reference list.
REFERENCES: tuple[Reference, ...] = (
    Reference(
        "NEB, Optimizing Restriction Endonuclease Reactions, for the digest and the star-activity "
        "conditions",
        url="https://www.neb.com/en-us/tools-and-resources/usage-guidelines/optimizing-restriction-endonuclease-reactions",
    ),
    Reference(
        "NEB, Double Digests, for cutting with two enzymes in one buffer",
        url="https://www.neb.com/en-us/tools-and-resources/usage-guidelines/double-digests",
    ),
    Reference(
        "NEB, Ligation Protocol with T4 DNA Ligase (NEB #M0202), for the ligation, its ratios "
        "and both incubations",
        url="https://www.neb.com/en-us/protocols/dna-ligation-with-t4-dna-ligase-m0202",
    ),
    Reference(
        "NEB, Monarch Spin DNA Gel Extraction Kit instruction manual, NEB #T1120, version "
        "2.0 10.25, for what a gel extraction recovers",
        url="https://www.neb.com/-/media/nebus/files/manuals/manualt1120.pdf",
    ),
    Reference(
        "NEB, Heat Inactivation, for stopping a digest and for the enzymes heat does not stop",
        url="https://www.neb.com/en-us/tools-and-resources/usage-guidelines/heat-inactivation",
    ),
    Reference(
        "NEB, Troubleshooting Guide for Cloning, for the four transformation controls",
        url="https://www.neb.com/en-us/tools-and-resources/troubleshooting-guides/troubleshooting-guide-for-cloning",
    ),
)

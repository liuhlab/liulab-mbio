"""The library-assembly bench protocol: its reactions, its programs, and the step order.

The protocol covers every round as one experiment. A round is two digests, a clean-up, a
ligation, a second clean-up, an electroporation, a growth and a prep, and the rounds run in the
order the scheme fills its positions.

Every number comes from `liulab_mbio.library.bench` and `liulab_mbio.library.coverage`, each
sourced in ``docs/research/protein-library-assembly.md``, or from the design the plan computed.
What the method leaves unpublished -- the ligase's units, the buffer's strength, the
electroporation settings -- is one line saying so rather than a number invented here.

The shared builders in `liulab_mbio.bench.steps` are shaped for a PCR, a gel and a heat-shock
transformation, and this method runs none of the three, so only `listed` is reused.
"""

from collections.abc import Sequence
from dataclasses import KW_ONLY, dataclass

from liulab_mbio import checks as judged
from liulab_mbio.bench.amounts import REFERENCES as AMOUNT_REFERENCES
from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.steps import listed
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.library.bench import (
    DIGEST_CELSIUS,
    DIGEST_SECONDS,
    DIGEST_VOLUME_UL,
    ENZYME_UL,
    GROWTH_CELSIUS,
    LIGATION_SECONDS,
    LIGATION_VOLUME_UL,
    MOLAR_RATIO,
    OUTGROWTH_SECONDS,
    RECOVERY_SECONDS,
    SPRI_AFTER_DIGEST,
    SPRI_AFTER_LIGATION,
    TRANSFORMATION_NG,
)
from liulab_mbio.library.bench import REFERENCES as BENCH_REFERENCES
from liulab_mbio.library.coverage import REFERENCES as COVERAGE_REFERENCES
from liulab_mbio.library.coverage import RoundCoverage
from liulab_mbio.library.parts import Part
from liulab_mbio.library.rounds import Round
from liulab_mbio.library.scheme import Scheme
from liulab_mbio.library.standard import PartList, Standard
from liulab_mbio.library.vector import Destination
from liulab_mbio.protocol import (
    OVERVIEW_CHARS,
    Check,
    Component,
    Incubation,
    Material,
    Protocol,
    ReactionTable,
    Reference,
    Stage,
    Step,
    ThermocyclerProgram,
    Timer,
    Troubleshooting,
)
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.sites import find_sites

#: The buffer both digests run in, and the ligase and buffer the ligation runs in. The method
#: names all three and publishes neither the ligase's units nor either buffer's strength.
CUTSMART = "CutSmart Buffer"
LIGASE = "T7 DNA Ligase"
LIGASE_BUFFER = "StickTogether DNA Ligase Buffer"

#: What each clean-up is done with. The method gives the two volume ratios and not the product.
SPRI_BEADS = "SPRI paramagnetic beads"

#: The electrocompetent strain the method names, and who sells it.
STRAIN = "Endura ElectroCompetent Cells"
STRAIN_SUPPLIER = "Lucigen"
STRAIN_CATALOG = "60242-2"

#: The hardware a round needs, which no reagent table covers.
EQUIPMENT: tuple[str, ...] = (
    f"Incubator or heat block at {DIGEST_CELSIUS:g} °C",
    "Magnetic rack for the bead clean-ups",
    "Electroporator and cuvettes",
    f"Shaking incubator at {GROWTH_CELSIUS:g} °C",
    "Spectrophotometer or fluorometer",
)


@dataclass(frozen=True, slots=True)
class RoundBench:
    """What one round takes at the bench, computed from that round's own lengths.

    Parameters
    ----------
    number
        Which round it is, counting from one.
    position
        The position it fills, which names its part list.
    destination_digest, donor_digest
        What goes into each of the two digests: the library built so far, and the part list.
    ligation
        What the ligation takes, the opened destination first and the released part list second.
    transformation
        The most of the purified ligation one electroporation takes.
    coverage
        What this round has to cover, and what the colonies asked for leave out.
    """

    number: int
    position: str
    _: KW_ONLY
    destination_digest: Amount
    donor_digest: Amount
    ligation: tuple[Amount, Amount]
    transformation: Amount
    coverage: RoundCoverage


def choppers(scheme: Scheme) -> tuple[tuple[Enzyme, ...], tuple[Enzyme, ...]]:
    """Return the blunt enzymes each digest carries, the internal digest's first.

    A chopper belongs to the digest whose discarded piece it cuts: the internal stuffer core the
    internal digest excises, or the external stuffers the external digest leaves behind. One that
    reads a site in both is named by both, and the scheme has already refused one that reads a
    site in neither.
    """
    core = SequenceRecord(scheme.internal_stuffer_core)
    flanks = [
        SequenceRecord(bases)
        for position in scheme.positions
        for bases in (position.external_stuffer_5, position.external_stuffer_3)
    ]
    inside = tuple(one for one in scheme.blunt if find_sites(core, one))
    outside = tuple(one for one in scheme.blunt if any(find_sites(flank, one) for flank in flanks))
    return inside, outside


def digest_reaction(
    dna: Amount,
    enzymes: Sequence[Enzyme],
    *,
    volume_ul: float = DIGEST_VOLUME_UL,
    reactions: int = 1,
) -> ReactionTable:
    """Return one digest: the DNA, each enzyme in turn, and the buffer and water that fill it.

    The buffer is one line with the water because the method names `CUTSMART` and not the
    strength it is supplied at, so its own volume is the supplier's to set.

    Raises
    ------
    ValueError
        If the DNA and the enzymes do not fit `volume_ul`.
    """
    components = [
        Component(
            dna.name,
            dna.volume_ul,
            final=f"{dna.pmol:g} pmol ({dna.nanograms:g} ng)",
            master_mix=False,
        ),
        *(Component(one.supplier_label, ENZYME_UL) for one in enzymes),
    ]
    used = sum(component.volume_ul for component in components)
    if used >= volume_ul:
        raise ValueError(
            f"the DNA and enzymes take {used:g} µL of a {volume_ul:g} µL digest; "
            "concentrate the DNA or scale the digest up"
        )
    components.append(
        Component(
            f"{CUTSMART} and nuclease-free water",
            round(volume_ul - used, 2),
            final=f"to {volume_ul:g} µL",
        )
    )
    return ReactionTable(
        tuple(components),
        title=f"Digest with {listed([one.name for one in enzymes])}",
        reactions=reactions,
    )


def digest_program(enzymes: Sequence[Enzyme]) -> ThermocyclerProgram:
    """Return the two hours a digest runs: the first enzyme alone, then the rest added to it."""
    first, rest = enzymes[0], tuple(enzymes[1:])
    stages = [Stage((Incubation(first.name, DIGEST_CELSIUS, DIGEST_SECONDS),))]
    if rest:
        added = listed([one.name for one in rest])
        stages.append(
            Stage((Incubation(f"{first.name} and {added}", DIGEST_CELSIUS, DIGEST_SECONDS),))
        )
    return ThermocyclerProgram(
        tuple(stages), title=f"Digest with {listed([one.name for one in enzymes])}"
    )


def ligation_reaction(
    amounts: tuple[Amount, Amount],
    *,
    volume_ul: float = LIGATION_VOLUME_UL,
    reactions: int = 1,
) -> ReactionTable:
    """Return one round's ligation: the two digest fragments, then the ligase, buffer and water.

    The last line carries all three together: the method publishes neither the ligase's units nor
    its volume, so the split between them is the supplier's to set and is not invented here.

    Raises
    ------
    ValueError
        If the two fragments do not fit `volume_ul`.
    """
    components = [
        Component(
            one.name,
            one.volume_ul,
            final=f"{one.pmol:g} pmol ({one.nanograms:g} ng)",
            master_mix=False,
        )
        for one in amounts
    ]
    used = sum(component.volume_ul for component in components)
    if used >= volume_ul:
        raise ValueError(
            f"the two fragments take {used:g} µL of a {volume_ul:g} µL ligation; "
            "concentrate them or scale the ligation up"
        )
    components.append(
        Component(
            f"{LIGASE} in {LIGASE_BUFFER}, and nuclease-free water",
            round(volume_ul - used, 2),
            final=f"to {volume_ul:g} µL",
        )
    )
    return ReactionTable(tuple(components), title="Ligation", reactions=reactions)


def growth_program() -> ThermocyclerProgram:
    """Return the recovery and the outgrowth, both at `GROWTH_CELSIUS` and not at 37 °C.

    The outgrowth carries the shorter of the two times the method gives; the step says the range.
    """
    return ThermocyclerProgram(
        (
            Stage((Incubation("Recovery, shaking", GROWTH_CELSIUS, RECOVERY_SECONDS),)),
            Stage((Incubation("Outgrowth", GROWTH_CELSIUS, OUTGROWTH_SECONDS[0]),)),
        ),
        title=f"Recovery and outgrowth at {GROWTH_CELSIUS:g} °C",
    )


def protocol(
    *,
    scheme: Scheme,
    vector: SequenceRecord,
    destination: Destination,
    part_lists: Sequence[PartList],
    standard: Standard,
    parts: Sequence[Part],
    rounds: Sequence[Round],
    bench: Sequence[RoundBench],
    constructs: int,
    checks: Sequence[judged.Check],
    host: str,
    sheet: str,
    barcodes: str,
) -> Protocol:
    """Return the bench protocol for one planned library, ready to render.

    Each argument is the `liulab_mbio.library.plan.LibraryPlan` field or property of that name;
    `sheet` and `barcodes` are what the plan calls the two files a step points at. The steps run
    in the order someone does them: order the blocks, pool each part list, then every round in
    turn, and finally read the barcode block back.
    """
    inside, outside = choppers(scheme)
    return Protocol(
        f"Library assembly: {len(part_lists)} part lists into {vector.name or 'the vector'}",
        summary=(
            f"Join {len(parts)} synthesised parts into {constructs} distinct constructs in "
            f"{len(rounds)} rounds. Each round opens the library with {scheme.internal.name}, "
            f"releases one part list with {scheme.external.name}, ligates the two, and "
            "transforms, grows and preps the result for the round after it."
        ),
        overview=_overview(scheme, standard, rounds, bench, constructs, part_lists, host),
        highlights=_highlights(scheme, destination, standard, rounds, bench, constructs, barcodes),
        checks=_checks(checks),
        materials=_materials(scheme, part_lists, vector, sheet),
        equipment=EQUIPMENT,
        steps=_steps(scheme, parts, part_lists, rounds, bench, inside, outside, sheet, barcodes),
        references=_references(scheme),
    )


def _card(value: str, short: str) -> str:
    """Return the fact where a card holds it, and the shorter form where it does not."""
    return value if len(value) <= OVERVIEW_CHARS else short


def _overview(
    scheme: Scheme,
    standard: Standard,
    rounds: Sequence[Round],
    bench: Sequence[RoundBench],
    constructs: int,
    part_lists: Sequence[PartList],
    host: str,
) -> dict[str, str]:
    """Return the facts to check before starting, each short enough to be a card."""
    last = bench[-1]
    product = rounds[-1].product
    sizes = ", ".join(
        f"{position.name} {len(parts)}"
        for position, parts in zip(scheme.positions, part_lists, strict=True)
    )
    return {
        "Scheme": _card(scheme.name, f"{scheme.position_count} positions"),
        "Part lists": _card(sizes, f"{len(part_lists)} lists"),
        "Parts": f"{sum(len(one) for one in part_lists)} synthesised blocks",
        "Constructs": f"{constructs:,} distinct",
        "Rounds": f"{len(rounds)}, one a part list",
        "Opened with": scheme.internal.supplier_label,
        "Released with": scheme.external.supplier_label,
        "Blunt enzymes": listed([one.name for one in scheme.blunt]),
        "Entry overhangs": _card(listed(standard.entry_overhangs), "one a position"),
        "Cloning scar": standard.scar_overhang,
        "Barcode": f"{scheme.barcode_length} bp, block {scheme.barcode_block_length} bp",
        "Codon usage": host,
        "Product": _card(f"{product.name}, {len(product)} bp", f"{len(product)} bp"),
        "Coverage": f"{last.coverage.coverage:g}x, {last.coverage.colonies:,} colonies at the end",
        "Amino acids changed": f"{standard.cost} over {len(standard.changes)} part end(s)",
    }


def _highlights(
    scheme: Scheme,
    destination: Destination,
    standard: Standard,
    rounds: Sequence[Round],
    bench: Sequence[RoundBench],
    constructs: int,
    barcodes: str,
) -> tuple[str, ...]:
    """Return what the facts mean, a sentence each."""
    last = bench[-1]
    said = [
        f"One round appends one part list to every member of the library at once, so "
        f"{len(rounds)} rounds make {constructs} constructs out of "
        f"{sum(one.coverage.part_list_size for one in bench)} synthesised parts.",
        f"The product keeps {scheme.internal.name}'s sites, which is what lets the next round "
        f"open it, and loses {scheme.external.name}'s with the external stuffers.",
    ]
    if destination.edit is None:
        said.append(
            f"Your vector already carried an internal stuffer, so every part enters position "
            f"{scheme.positions[0].name} on {standard.entry_overhangs[0]}, which is the overhang "
            "that stuffer spells: DNA that exists cannot be re-chosen."
        )
    else:
        said.append(
            f"Your vector carried no internal stuffer, so one was put in at "
            f"{destination.stuffer.start} carrying {standard.entry_overhangs[0]}, the overhang "
            "this design chose. Read the edit before you order the vector."
        )
    if standard.cost:
        said.append(
            f"The overhang standard changes {standard.cost} amino acid(s) across "
            f"{len(standard.changes)} part end(s); changes.tsv has wild type beside synthesised, "
            "and that is what you are about to pay for."
        )
    else:
        said.append(
            "The overhang standard changes no amino acid: every part list already spells its "
            "junctions."
        )
    said.append(
        f"The last round needs {last.coverage.colonies:,} colonies for "
        f"{last.coverage.coverage:g}x coverage of its {last.coverage.products:,} products. A "
        "round short of that loses members no later round can put back."
    )
    said.append(
        f"The finished barcode block is {scheme.barcode_block_length} bp and reads in the "
        f"reverse of the order the rounds ran; {barcodes} says which barcode names which part."
    )
    return tuple(said)


def _checks(checks: Sequence[judged.Check]) -> tuple[Check, ...]:
    """Return the plan's verdicts, one badge each, so a warning is seen and not read."""
    return tuple(
        Check(check.name, check.status, detail=check.detail)
        for check in checks
        if check.status is not None
    )


def _enzyme_material(enzyme: Enzyme, note: str) -> Material:
    """Return one enzyme as a material, its own record carrying the catalogue number."""
    return Material(
        enzyme.commercial_name or enzyme.name,
        supplier=enzyme.supplier or "",
        catalog=enzyme.catalog_number or "",
        storage="-20 °C",
        amount=f"{ENZYME_UL:g} µL per digest",
        note=note,
    )


def _materials(
    scheme: Scheme,
    part_lists: Sequence[PartList],
    vector: SequenceRecord,
    sheet: str,
) -> tuple[Material, ...]:
    """Every reagent and consumable the protocol asks for. The parts are the order sheet."""
    inside, outside = choppers(scheme)
    made = [
        Material(
            f"{position.name} part list",
            storage="-20 °C",
            amount=f"{len(one)} members, pooled",
            note=f"synthesised blocks, ordered from {sheet}",
        )
        for position, one in zip(scheme.positions, part_lists, strict=True)
    ]
    made.append(
        Material(
            f"{vector.name or 'destination'} vector",
            storage="-20 °C",
            note="opened by the internal digest in round 1",
        )
    )
    made.append(_enzyme_material(scheme.internal, "opens the library, excising its stuffer"))
    made.append(_enzyme_material(scheme.external, "releases a part from its synthesised block"))
    for one in scheme.blunt:
        jobs = [
            what
            for what, named in (("the excised stuffer", inside), ("the donor backbone", outside))
            if one in named
        ]
        made.append(_enzyme_material(one, f"cuts {listed(jobs)}, so it cannot ligate back"))
    made.append(Material(CUTSMART, storage="-20 °C", note="both digests run in it"))
    made.append(
        Material(
            f"{LIGASE} and {LIGASE_BUFFER}",
            storage="-20 °C",
            note="the method names neither the units nor the volume; follow the supplier",
        )
    )
    made.append(
        Material(
            SPRI_BEADS,
            amount=f"{SPRI_AFTER_DIGEST:g} volumes after a digest, "
            f"{SPRI_AFTER_LIGATION:g} after a ligation",
            note="the two ratios differ; both elute in water",
        )
    )
    made.append(
        Material(
            STRAIN,
            supplier=STRAIN_SUPPLIER,
            catalog=STRAIN_CATALOG,
            storage="-80 °C",
            amount="one aliquot per round",
        )
    )
    made.append(Material("Recovery medium", amount="one outgrowth per round"))
    made.append(
        Material(
            "Selective broth and plates",
            note="the destination vector's own antibiotic, which this plan does not name",
        )
    )
    made.append(Material("Plasmid prep kit", amount="one prep per round"))
    made.append(Material("Electroporation cuvettes", amount=f"{len(part_lists)}, one per round"))
    return tuple(made)


def _references(scheme: Scheme) -> tuple[Reference, ...]:
    """Where the numbers come from, and where the scheme itself came from."""
    items = [*BENCH_REFERENCES, *COVERAGE_REFERENCES, *AMOUNT_REFERENCES]
    if scheme.source:
        items.append(Reference(f"The scheme this build was given: {scheme.source}"))
    return tuple(items)


def _steps(
    scheme: Scheme,
    parts: Sequence[Part],
    part_lists: Sequence[PartList],
    rounds: Sequence[Round],
    bench: Sequence[RoundBench],
    inside: Sequence[Enzyme],
    outside: Sequence[Enzyme],
    sheet: str,
    barcodes: str,
) -> tuple[Step, ...]:
    """Return every step in the order it happens, the rounds one after another."""
    made = [_order_step(parts, sheet), _pool_step(scheme, part_lists)]
    for one, row in zip(rounds, bench, strict=True):
        made.extend(_round_steps(scheme, one, row, inside, outside, len(rounds)))
    made.append(_confirm_step(scheme, rounds, barcodes))
    return tuple(made)


def _order_step(parts: Sequence[Part], sheet: str) -> Step:
    """Order every block, which is what the whole design comes down to."""
    lengths = [part.length for part in parts]
    return Step(
        "Order the synthesised parts",
        instructions=(
            f"Order every row of {sheet} as a double-stranded synthesised block.",
            "Ask for sequence-verified material: a block carries the scheme's sites at fixed "
            "offsets, and a base out of place stops it being cut where the design says.",
        ),
        expected=(
            f"{len(parts)} blocks, {min(lengths)} to {max(lengths)} bp.",
            "Each block reads: 5' external stuffer, coding bases, internal stuffer, barcode, "
            "3' external stuffer.",
        ),
        notes=(
            f"{sheet} carries each block's barcode on its own row, so the sheet you order from "
            "is also what decodes the sequencing afterwards.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The vendor cannot synthesise a block",
                "It is usually a repeat or a GC run in the coding bases. Plan again with another "
                "codon usage table; the stuffers and the overhangs are fixed by the scheme.",
            ),
        ),
    )


def _pool_step(scheme: Scheme, part_lists: Sequence[PartList]) -> Step:
    """Pool each part list, which is what a round joins in one tube."""
    return Step(
        "Pool each part list",
        instructions=(
            "Resuspend every block and measure each concentration.",
            "Pool the members of each part list in equal molar amounts, one tube per position.",
        ),
        expected=tuple(
            f"{position.name}: one tube holding {len(one)} member(s)."
            for position, one in zip(scheme.positions, part_lists, strict=True)
        ),
        notes=(
            "Library coverage is counted on equally represented members, so an uneven pool loses "
            "members that no later round can put back.",
        ),
        troubleshooting=(
            Troubleshooting(
                "One member is much more dilute than the rest",
                "Pool to the lowest member rather than to the mean; a member short here is short "
                "in every round after it.",
            ),
        ),
    )


def _round_steps(
    scheme: Scheme,
    one: Round,
    row: RoundBench,
    inside: Sequence[Enzyme],
    outside: Sequence[Enzyme],
    total: int,
) -> list[Step]:
    """Return the eight steps of one round, in the order they happen."""
    number = row.number
    opened, released = row.ligation
    internal = (scheme.internal, *inside)
    external = (scheme.external, *outside)
    return [
        _open_step(scheme, one, row, internal, opened),
        _release_step(scheme, row, external, released),
        _digest_cleanup_step(row, opened, released),
        _ligation_step(row, opened, released),
        _ligation_cleanup_step(row),
        _electroporation_step(row),
        _growth_step(row),
        _prep_step(scheme, one, row, number == total),
    ]


def _digest_instructions(enzymes: Sequence[Enzyme]) -> tuple[str, ...]:
    """Return the two hours a digest runs, the second enzyme going in after the first."""
    first, rest = enzymes[0], tuple(enzymes[1:])
    said = [f"Set the digest below up with {first.name} alone, and start the program."]
    if rest:
        said.append(
            f"After the first hour add {listed([one.name for one in rest])} and run the second."
        )
    return tuple(said)


def _open_step(
    scheme: Scheme, one: Round, row: RoundBench, enzymes: Sequence[Enzyme], opened: Amount
) -> Step:
    """Open the library built so far, which is what the round's part enters."""
    chopped = listed([enzyme.name for enzyme in enzymes[1:]]) or "nothing else"
    return Step(
        f"Round {row.number}: open {row.destination_digest.name} with {scheme.internal.name}",
        instructions=_digest_instructions(enzymes),
        tables=(digest_reaction(row.destination_digest, enzymes),),
        programs=(digest_program(enzymes),),
        expected=(
            f"{row.destination_digest.name} opens on {one.entry_overhang} and "
            f"{one.scar_overhang}, leaving {opened.length_bp} bp of backbone.",
            f"The {one.excised.length} bp internal stuffer comes out, and {chopped} cuts it so "
            "it cannot go back in.",
        ),
        notes=(
            f"{scheme.internal.name} is the enzyme the product keeps sites for. That is the "
            "design working, not a site left behind.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Mostly empty library later",
                "The destination was not opened to completion. Run a little of the digest on a "
                "gel before ligating the next round.",
            ),
        ),
    )


def _release_step(
    scheme: Scheme, row: RoundBench, enzymes: Sequence[Enzyme], released: Amount
) -> Step:
    """Cut the round's whole part list out of its blocks, in one tube."""
    chopped = listed([enzyme.name for enzyme in enzymes[1:]]) or "nothing else"
    return Step(
        f"Round {row.number}: release the {row.position} part list with {scheme.external.name}",
        instructions=_digest_instructions(enzymes),
        tables=(digest_reaction(row.donor_digest, enzymes),),
        programs=(digest_program(enzymes),),
        expected=(
            f"Every member of the {row.position} part list is released, about "
            f"{released.length_bp} bp, on the same two overhangs the destination now offers.",
            f"{chopped} cuts the two external stuffers left behind, so neither can ligate back.",
        ),
        notes=(
            "One tube takes the whole part list: every member carries the same stuffers and "
            "differs only in its coding bases and its barcode.",
            "The amounts are weighed at the pool's mean block length, "
            f"{row.donor_digest.length_bp} bp.",
        ),
        troubleshooting=(
            Troubleshooting(
                "A member is missing from the library later",
                "Its block was under-represented in the pool, or its own digest was incomplete. "
                "Check the pool before repeating the round.",
            ),
        ),
    )


def _digest_cleanup_step(row: RoundBench, opened: Amount, released: Amount) -> Step:
    """Take the enzymes and the shredded pieces away, and measure what is left."""
    return Step(
        f"Round {row.number}: clean both digests up",
        instructions=(
            f"Add {SPRI_AFTER_DIGEST:g} volumes of {SPRI_BEADS} to each digest and elute in water.",
            "Measure both concentrations; the ligation table asks for picomoles, not nanograms.",
        ),
        expected=(
            f"{opened.name}: {opened.pmol:g} pmol is {opened.nanograms:g} ng at "
            f"{opened.length_bp} bp.",
            f"{released.name}: {released.pmol:g} pmol is {released.nanograms:g} ng at "
            f"{released.length_bp} bp.",
        ),
        notes=(
            f"{SPRI_AFTER_DIGEST:g} volumes here and {SPRI_AFTER_LIGATION:g} after the ligation. "
            "The two ratios differ; they are not one number used twice.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Too dilute for the ligation",
                "Elute in a smaller volume, or scale the ligation up.",
            ),
        ),
    )


def _ligation_step(row: RoundBench, opened: Amount, released: Amount) -> Step:
    """Join the opened library and the released part list, at a molar ratio and not a mass one."""
    return Step(
        f"Round {row.number}: ligate the {row.position} part list into the library",
        instructions=(
            "Pipette the two DNAs into the tube first, then the ligase, its buffer and water.",
            f"Hold at room temperature for {LIGATION_SECONDS // 60} minutes.",
        ),
        tables=(ligation_reaction(row.ligation),),
        timers=(Timer("Ligation", LIGATION_SECONDS),),
        expected=(
            f"Nothing visible. {released.pmol:g} pmol of the released part list against "
            f"{opened.pmol:g} pmol of the opened library is the {MOLAR_RATIO:g}:1 molar ratio.",
        ),
        notes=(
            "Picomoles, not nanograms: the shorter fragment weighs less at the same ratio.",
            f"The method names {LIGASE} and {LIGASE_BUFFER} and publishes neither the units nor "
            "the volume, so the last line of the table is the supplier's own.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Few colonies after this round",
                "Check that both digests went to completion; an uncut end cannot ligate.",
            ),
        ),
    )


def _ligation_cleanup_step(row: RoundBench) -> Step:
    """Desalt the ligation, which is what stops the cuvette arcing."""
    return Step(
        f"Round {row.number}: clean the ligation up",
        instructions=(
            f"Add {SPRI_AFTER_LIGATION:g} volume of {SPRI_BEADS} and elute in water.",
            "Measure the concentration.",
        ),
        expected=(
            f"Enough for one electroporation: {row.transformation.nanograms:g} ng is "
            f"{row.transformation.pmol:g} pmol at {row.transformation.length_bp} bp.",
        ),
        notes=(
            f"One volume here, not the {SPRI_AFTER_DIGEST:g} the digests took. Eluting in water "
            "rather than a buffer is what keeps the salt out of the cuvette.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The eluate is still salty",
                "Repeat the clean-up; salt carried into the cuvette arcs the pulse.",
            ),
        ),
    )


def _electroporation_step(row: RoundBench) -> Step:
    """Get the whole ligation into cells, which is where the library's size is won or lost."""
    return Step(
        f"Round {row.number}: electroporate into {STRAIN}",
        instructions=(
            f"Thaw one aliquot of {STRAIN} on ice.",
            f"Add at most {TRANSFORMATION_NG:g} ng of the purified ligation and mix without "
            "making bubbles.",
            "Pulse with the settings the cell supplier gives for your cuvette.",
        ),
        cautions=("Keep the cells and the cuvette on ice; a warm cuvette arcs.",),
        expected=(
            "A pulse with no arc, and a time constant in the range the cell supplier's manual "
            "gives.",
        ),
        notes=(
            "The method names the instrument and not the voltage, the capacitance, the "
            "resistance or the cuvette gap, so the settings are the cell supplier's.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The cuvette arcs",
                "The DNA carries salt. Clean it up again, elute in water, and use a fresh "
                "aliquot of cells.",
            ),
        ),
    )


def _growth_step(row: RoundBench) -> Step:
    """Recover and grow, both at 30 °C, and count what the round actually got."""
    coverage = row.coverage
    return Step(
        f"Round {row.number}: recover and grow at {GROWTH_CELSIUS:g} °C",
        instructions=(
            "Add recovery medium straight away and shake for the first hour.",
            "Plate a measured dilution on selection to count the transformants, and grow the "
            "rest in selective broth.",
        ),
        programs=(growth_program(),),
        expected=(
            f"At least {coverage.colonies:,} colonies, scaled up from the dilution: "
            f"{coverage.coverage:g}x over the {coverage.products:,} distinct products this round "
            "can make.",
            f"At that count the chance a named product is missing is "
            f"{coverage.absent_probability:.3g}.",
        ),
        notes=(
            f"Both steps run at {GROWTH_CELSIUS:g} °C and not at 37 °C. That is a library "
            "precaution rather than an oversight.",
            f"The program carries the shorter outgrowth; {OUTGROWTH_SECONDS[0] // 3600} to "
            f"{OUTGROWTH_SECONDS[1] // 3600} hours is the range the method gives.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Fewer colonies than the count above",
                "The round has lost library members and no later round can put them back. "
                "Electroporate more of the ligation, or run the round again.",
            ),
        ),
    )


def _prep_step(scheme: Scheme, one: Round, row: RoundBench, last: bool) -> Step:
    """Prep the round's plasmid, which is either the next destination or the finished library."""
    where = (
        "This is the finished library."
        if last
        else f"This is the destination round {row.number + 1} opens."
    )
    return Step(
        f"Round {row.number}: prep the library",
        instructions=(
            "Harvest the whole culture rather than a single colony.",
            "Prep the plasmid and measure the concentration.",
        ),
        expected=(
            f"One pool of plasmid, {len(one.product)} bp per molecule. {where}",
            f"{scheme.internal.name} still cuts it and {scheme.external.name} no longer does.",
        ),
        notes=(
            "Harvesting the whole culture is what keeps the library a library; picking colonies "
            "here throws away everything not picked.",
        ),
        troubleshooting=(
            Troubleshooting(
                "A low yield",
                "The culture was grown too briefly, or the antibiotic was wrong for this "
                "backbone. Check the vector's own marker.",
            ),
        ),
    )


def _confirm_step(scheme: Scheme, rounds: Sequence[Round], barcodes: str) -> Step:
    """Read the barcode block back, which is what links a construct to its parts."""
    final = rounds[-1]
    order = listed([one.position for one in reversed(rounds)])
    return Step(
        "Confirm the library",
        instructions=(
            f"Sequence across the {scheme.barcode_block_length} bp barcode block of the finished "
            "library, at "
            f"{final.block.start}-{final.block.end} of the representative construct.",
            f"Decode each read against {barcodes}.",
        ),
        expected=(
            f"The block reads {order}, each barcode separated from the last by the cloning scar "
            f"{scheme.cloning_scar}.",
            f"Every construct carries one member of each part list, and the {len(rounds)} "
            "barcodes say which.",
        ),
        notes=(
            "The block reads in the reverse of the order the rounds ran: each round inserted its "
            "barcode ahead of the ones already there.",
            "This plan designs no sequencing primers.",
        ),
        troubleshooting=(
            Troubleshooting(
                "A read carries fewer barcodes than there are rounds",
                "A round did not go in. Check that round's prep against its own record before "
                "blaming the sequencing.",
            ),
        ),
    )

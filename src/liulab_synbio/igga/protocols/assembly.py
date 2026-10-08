"""Protocol 05: join the part lists round by round, then read the finished library back.

A round is two digests, a clean-up, a ligation, a second clean-up, an electroporation, a growth
and a prep, and the rounds run in the order the scheme fills its positions. The two reads at
the end are the library's own: linkage once, representation at every bottleneck after it.

Every number comes from `liulab_synbio.igga.bench` and `liulab_mbio.bench.coverage`, each
sourced in ``docs/research/protein-library-assembly.md``, or from the design the plan computed.
What the method leaves unpublished -- the ligase's units, the buffer's strength -- is one line
saying so rather than a number invented here.
"""

import math
from collections.abc import Mapping, Sequence

from liulab_mbio.barcodes import deletion_ambiguity
from liulab_mbio.bench.amounts import Amount, to_nanograms
from liulab_mbio.bench.steps import card, listed
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.protocol.figures import SOURCE as FIGURE_SOURCE
from liulab_mbio.protocol.figures import SOURCE_KEY as FIGURE_SOURCE_KEY
from liulab_mbio.protocol.model import (
    Figure,
    Hole,
    Reference,
    Source,
    Step,
    Timer,
    Troubleshooting,
    number,
)
from liulab_mbio.protocol.model import Item as Handed
from liulab_synbio.igga import stages
from liulab_synbio.igga.bench import (
    DIGEST_VOLUME_UL,
    ENZYME_UL,
    GROWTH_CELSIUS,
    LIGASE,
    LIGASE_BUFFER,
    LIGATION_SECONDS,
    MOLAR_RATIO,
    OUTGROWTH_SECONDS,
    SPRI_AFTER_DIGEST,
    SPRI_AFTER_LIGATION,
    SPRI_BEADS,
    STRAIN,
    STRAIN_CATALOG,
    TRANSFORMATION_NG,
    RoundBench,
    digest_program,
    digest_reaction,
    growth_program,
    ligation_reaction,
    pool_floor_ng_ul,
    pulse,
)
from liulab_synbio.igga.figures import assembly_rows
from liulab_synbio.igga.method import Scheme
from liulab_synbio.igga.parts import Part
from liulab_synbio.igga.protocols.protocol import Protocol, figured, labelled
from liulab_synbio.igga.protocols.run import (
    ROUND_EQUIPMENT,
    Run,
    as_platform,
    marks_sentence,
    with_pair,
)
from liulab_synbio.igga.reads import ReadPair
from liulab_synbio.igga.rounds import Round

#: What the page is headed and what the chain names it by.
ASSEMBLY = "Library assembly in rounds"


class Assembly(Protocol):
    """Pool each part list, run every round, and read linkage and representation back."""

    round_reagents = True

    def title(self, run: Run) -> str:
        """Return the page's own heading."""
        return ASSEMBLY

    def summary(self, run: Run) -> str:
        """Return what the rounds come to, which is the whole library."""
        return (
            f"Join {len(run.part_lists)} part lists into {run.constructs} distinct constructs in "
            f"{len(run.rounds)} rounds, each round opening the library with "
            f"{run.scheme.internal.name} and ligating one part list into it, then read the "
            "finished library back."
        )

    def steps(self, run: Run) -> tuple[Step, ...]:
        """Return the pool, every round's eight steps, and the two reads that close the protocol."""
        made = [labelled(_pool_step(run), "Pool the part lists")]
        vector = _vector_record(run)
        for place, (one, row) in enumerate(zip(run.rounds, run.bench, strict=True), 1):
            rows = assembly_rows(run.rounds, lit=place, vector=vector, at=run.records_at)
            made += [labelled(step, f"Round {place}") for step in _round_steps(run, one, row, rows)]
        reads = run.reads
        made += [
            labelled(
                _linkage_step(run, None if reads is None else reads.linkage),
                "Read the library back",
            ),
            labelled(
                _representation_step(run, None if reads is None else reads.representation),
                "Read the library back",
            ),
        ]
        return tuple(made)

    def consumes(self, run: Run) -> tuple[Handed, ...]:
        """Return the cargo, as the protocol before it left it, and the vector round 1 opens.

        A run ordering its blocks whole archives nothing, so its cargo is the vendor's tube.
        """
        if run.validation is not None:
            cargo = (run.picked, run.calls)
        else:
            cargo = (run.archive,) if run.pool else (run.ordered,)
        return (*cargo, *run.blocks[:1])

    def produces(self, run: Run) -> tuple[Handed, ...]:
        """Return the pooled library after the last round, as a plasmid prep."""
        return (run.prep,)

    def equipment(self, run: Run) -> tuple[str, ...]:
        """Return the hardware a round needs, which no reagent table covers."""
        return ROUND_EQUIPMENT

    def references(self, run: Run) -> tuple[Reference, ...]:
        """Where a round's numbers come from, and where the scheme itself came from."""
        return run.round_references

    def sources(self, run: Run) -> dict[str, Source]:
        """Return what the method's own materials and the step figures are cited to."""
        return dict(stages.SOURCES) | {FIGURE_SOURCE_KEY: FIGURE_SOURCE}

    def holes(self, run: Run) -> tuple[Hole, ...]:
        """Return what the destination's own record leaves unanswered."""
        return stages.holes_for(run.vector)

    def overview(self, run: Run) -> dict[str, str]:
        """Return the facts to check before starting, each short enough to be a card."""
        last = run.bench[-1]
        product = run.rounds[-1].product
        blocks = sum(len(one) for one in run.part_lists)
        sizes = ", ".join(
            f"{position} {len(parts)}"
            for position, parts in zip(run.positions, run.part_lists, strict=True)
        )
        scheme, standard, validation = run.scheme, run.standard, run.validation
        return {
            "Method": card(scheme.name, f"{len(run.positions)} positions"),
            "Part lists": card(sizes, f"{len(run.part_lists)} lists"),
            "Parts": (
                f"{blocks} blocks, assembled from {run.pool.pool.count} oligos"
                if run.pool
                else f"{blocks} synthesised blocks"
            ),
            "Constructs": f"{run.constructs:,} distinct",
            "Rounds": f"{len(run.rounds)}, one a part list",
            "Opened with": scheme.internal.supplier_label,
            "Released with": scheme.external.supplier_label,
            "Blunt enzymes": listed([one.name for one in scheme.blunt]),
            "Entry overhangs": card(listed(standard.entry_overhangs), "one a position"),
            "Cloning scar": standard.scar_overhang,
            "Barcode": (
                f"{run.barcode_length} bp, block "
                f"{scheme.barcode_block_length(run.barcode_length, len(run.positions))} bp"
            ),
            "Codon usage": run.host,
            "Product": card(f"{product.name}, {len(product)} bp", f"{len(product)} bp"),
            "Completeness": card(
                f"{last.coverage.completeness:g} chance nothing is missing, "
                f"{last.coverage.colonies:,} colonies at the end",
                f"{last.coverage.coverage:.0f}x",
            ),
            "Amino acids changed": (f"{standard.cost} over {len(standard.changes)} part end(s)"),
            "Designs read back": (
                card(
                    f"{len(validation.designs)}, every one"
                    if validation.floor == 0
                    else f"{len(validation.designs)}, from {validation.floor} fragment(s)",
                    f"by {validation.route.name}",
                )
                if validation
                else "none; the library stays polyclonal"
            ),
        }


def _pool_step(run: Run) -> Step:
    """Pool each part list, which is what a round joins in one tube.

    Every number here is the round's own donor digest read backwards. That digest takes
    `DIGEST_NG` of the pool in the volume its two enzymes leave it, so the pool's total, its
    concentration floor and each member's share are fixed downstream rather than chosen.
    """
    bench, parts, pool = run.bench, run.parts, run.pool
    floor = _stated_floor()
    left = DIGEST_VOLUME_UL - 2 * ENZYME_UL
    members = {position: _at(parts, position) for position in {row.position for row in bench}}
    first = (
        "Clean up each assembly and measure each concentration."
        if pool
        else "Spin each tube down, resuspend in TE pH 8.0 or 10 mM Tris-HCl pH 8.0, 50 °C for "
        "15-20 min, and measure each concentration."
    )
    return Step(
        "Pool each part list",
        key="pool-part-lists",
        instructions=(
            first,
            "Pool the members of each part list in equal picomoles, one tube a position: each "
            "member gives the position's total divided by its member count, which is unequal "
            "masses because the members differ in length.",
            f"Bring each pool to at least {floor:g} ng/µL. The round's {DIGEST_VOLUME_UL:g} µL "
            f"digest leaves {left:g} µL for the DNA after its two {ENZYME_UL:g} µL enzymes, and "
            "the whole total has to arrive in it.",
        ),
        expected=tuple(
            f"{row.position}: one tube, {len(members[row.position])} member(s), at least "
            f"{row.donor_digest.nanograms:,.0f} ng at {floor:g} ng/µL or above, "
            f"{number(row.donor_digest.pmol / len(members[row.position]))} pmol of each member."
            for row in bench
            if members[row.position]
        ),
        notes=(
            "Library coverage is counted on equally represented members, so an uneven pool loses "
            "members that no later round can put back.",
            "Equal picomoles are unequal masses: a short member weighs less than a long one for "
            "the same number of molecules, and weighing them equally would not pool them equally.",
            *_pool_masses(bench, members),
        ),
        troubleshooting=(
            Troubleshooting(
                "One member is much more dilute than the rest",
                "Pool to the lowest member rather than to the mean; a member short here is short "
                "in every round after it.",
            ),
            Troubleshooting(
                f"The pool is below {floor:g} ng/µL",
                f"Concentrate it, by SPRI at the round's own {SPRI_AFTER_DIGEST:g}x ratio, or "
                "scale the digest up so the same mass arrives in a larger volume.",
            ),
        ),
    )


def _at(parts: Sequence[Part], position: str) -> tuple[Part, ...]:
    """Return the parts filling one position, in the order they were designed."""
    return tuple(one for one in parts if one.position == position)


def _stated_floor() -> float:
    """Return the pool's floor as a page states it, rounded up so it is never below the real one."""
    return math.ceil(pool_floor_ng_ul() * 10) / 10


def _pool_masses(
    bench: Sequence[RoundBench], members: Mapping[str, Sequence[Part]]
) -> tuple[str, ...]:
    """Return what a position's equal picomoles weigh, shortest member to longest."""
    said = []
    for row in bench:
        each = members[row.position]
        if not each:
            continue
        pmol = row.donor_digest.pmol / len(each)
        shortest, longest = min(each, key=_len_of), max(each, key=_len_of)
        said.append(
            f"{row.position}: {number(pmol)} pmol a member is "
            f"{number(to_nanograms(pmol, shortest.length))} ng of its shortest at "
            f"{shortest.length:,} bp and {number(to_nanograms(pmol, longest.length))} ng of its "
            f"longest at {longest.length:,} bp."
        )
    return tuple(said)


def _len_of(one: Part) -> int:
    """How many bases the part is ordered as."""
    return one.length


def _vector_record(run: Run) -> str:
    """Return what the plan calls the vector round 1 opens, or nothing where it writes none.

    The first block vector is that vector: each is this build's own destination respelt for the
    position whose blocks it holds, and position one's enters on the overhang it already spells.
    The file is returned as the run names it, for a figure to say where the records sit.
    """
    if not run.block_vectors:
        return ""
    _, file = run.block_vectors[0]
    return file


def _round_steps(run: Run, one: Round, row: RoundBench, rows: Figure | None) -> list[Step]:
    """Return the eight steps of one round, in the order they happen.

    `rows` is drawn on the ligation, the one step of the eight that changes a molecule.
    """
    scheme = run.scheme
    opened, released = row.ligation
    internal = (scheme.internal, *run.inside)
    external = (scheme.external, *run.outside)
    return [
        _open_step(scheme, one, row, internal, opened),
        _release_step(scheme, row, external, released),
        _digest_cleanup_step(row, opened, released),
        figured(_ligation_step(row, opened, released), rows),
        _ligation_cleanup_step(row),
        _electroporation_step(row),
        _growth_step(row, run.selection),
        _prep_step(scheme, one, row, row.number == len(run.rounds)),
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
        key=f"round-{row.number}-open",
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
        key=f"round-{row.number}-release",
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
        key=f"round-{row.number}-digest-cleanup",
        instructions=(
            f"Add {SPRI_AFTER_DIGEST:g} volumes of {SPRI_BEADS} to each digest and elute in water.",
            "Measure both concentrations; the ligation table asks for picomoles, not nanograms.",
        ),
        expected=(
            f"{opened.name}: {number(opened.pmol)} pmol is {opened.nanograms:g} ng at "
            f"{opened.length_bp} bp.",
            f"{released.name}: {number(released.pmol)} pmol is {released.nanograms:g} ng at "
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
        key=f"round-{row.number}-ligate",
        instructions=(
            "Pipette the two DNAs into the tube first, then the ligase, its buffer and water.",
            f"Hold at room temperature for {LIGATION_SECONDS // 60} minutes.",
        ),
        tables=(ligation_reaction(row.ligation),),
        timers=(Timer("Ligation", LIGATION_SECONDS),),
        expected=(
            f"Nothing visible. {number(released.pmol)} pmol of the released part list against "
            f"{number(opened.pmol)} pmol of the opened library is the "
            f"{MOLAR_RATIO:g}:1 molar ratio.",
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
        key=f"round-{row.number}-ligation-cleanup",
        instructions=(
            f"Add {SPRI_AFTER_LIGATION:g} volume of {SPRI_BEADS} and elute in water.",
            "Measure the concentration.",
        ),
        expected=(
            f"Enough for one electroporation: {row.transformation.nanograms:g} ng is "
            f"{number(row.transformation.pmol)} pmol at {row.transformation.length_bp} bp.",
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
    shot = pulse()
    low, high = shot.time_constant_ms
    return Step(
        f"Round {row.number}: electroporate into {STRAIN}",
        key=f"round-{row.number}-electroporate",
        instructions=(
            f"Thaw one aliquot of {STRAIN} on ice, {shot.cells_ul:g} µL a pulse.",
            f"Add at most {TRANSFORMATION_NG:g} ng of the purified ligation and mix without "
            "making bubbles.",
            f"Pulse at {shot.volts:g} V, {shot.ohms:g} Ω and {shot.microfarads:g} µF in a "
            f"{shot.cuvette_mm:g} mm cuvette.",
        ),
        expected=(f"A pulse with no arc, and a time constant of {low:g} to {high:g} ms.",),
        notes=(
            "The method names the instrument and none of the settings, so the program is the "
            f"cells' own: it is keyed by {STRAIN_CATALOG} and changes when the cells do. The "
            "materials table says where it was read.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The cuvette arcs",
                "The DNA carries salt. Clean it up again, elute in water, and use a fresh "
                "aliquot of cells.",
            ),
        ),
    )


def _growth_step(row: RoundBench, selection: str) -> Step:
    """Recover and grow, both at 30 °C, and bound the round on net colonies.

    Two plates, both growing during the outgrowth: a measured dilution of the recovery, and the
    same cut destination carried through the ligation with no donor added. What grows on the
    control is parental destination that survived, so the round is judged on the difference.
    """
    coverage = row.coverage
    dilution, control = stages.titre_plates(row.number)
    return Step(
        f"Round {row.number}: recover and grow at {GROWTH_CELSIUS:g} °C",
        key=f"round-{row.number}-grow",
        instructions=(
            "Add recovery medium straight away and shake for the first hour.",
            f"Plate a measured dilution of the recovery on {selection or 'selection'} as "
            f"{dilution.name}, and grow the rest in selective broth.",
            f"Plate the no-donor ligation from the same digest as {control.name}, at the same "
            "dilution.",
        ),
        programs=(growth_program(),),
        expected=(
            f"At least {coverage.colonies:,} net colonies, scaled up from the dilution: the "
            f"floor for the {coverage.completeness:g} chance this design asked for that none of "
            f"its {coverage.products:,} distinct products is missing, equally represented.",
            f"At that count the chance a named product is missing is "
            f"{number(coverage.absent_probability)}.",
            f"{control.name} should be near empty beside it; its colonies come off the count.",
        ),
        notes=(
            f"Both steps run at {GROWTH_CELSIUS:g} °C and not at 37 °C. That is a library "
            "precaution rather than an oversight.",
            f"The program carries the shorter outgrowth; {OUTGROWTH_SECONDS[0] // 3600} to "
            f"{OUTGROWTH_SECONDS[1] // 3600} hours is the range the method gives.",
            f"The {coverage.completeness:g} was chosen for this run, and follows from the "
            "representation the screen downstream asks for. No source sets it, and nothing here "
            f"defaults it; it works out at {coverage.coverage:.0f}x this round's products.",
            "The control measures the chain the design rests on — two cuts, a blunt chopper, a "
            "ligase that refuses blunt ends and the 2x clean-up — rather than assuming it.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Fewer net colonies than the count above",
                "The round has lost library members and no later round can put them back. "
                "Electroporate more of the ligation, or run the round again.",
            ),
            Troubleshooting(
                "The no-donor control is not empty",
                "The destination was not opened to completion, or the blunt chopper missed. "
                "What grows is the round before's library, one position short, and it reaches "
                "the linkage read as truncated members.",
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
        key=f"round-{row.number}-prep",
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


def _linkage_step(run: Run, pair: ReadPair | None) -> Step:
    """Read the whole cargo back, which is what says the barcode block still names its parts.

    Read once. A barcode that names the wrong member is the one fault no later round repairs,
    so this is where the library is carried forward or a round is sent back.

    A build stating `linkage_fidelity` sets the mark nobody published, which is what H28 says
    closes it; a build stating none carries the hole.
    """
    scheme, rounds, positions = run.scheme, run.rounds, run.positions
    final = rounds[-1]
    ambiguous = max(
        deletion_ambiguity([one.barcode for one in run.parts if one.index == index])
        for index in range(len(positions))
    )
    order = listed([one.position for one in reversed(rounds)])
    fidelity = run.linkage_fidelity
    barcodes = run.barcodes
    return Step(
        "Read linkage",
        key="read-linkage",
        instructions=(
            f"Amplify the whole cargo out of the finished library, from the vector before the "
            f"first {rounds[0].entry_overhang} to the vector past the final "
            f"{rounds[0].scar_overhang}, so one read carries a member's parts and its barcode "
            f"block together{with_pair(pair, run.read_sheet)}.",
            f"Sequence the amplicon{as_platform(pair)}.",
            f"Decode each read's block against {barcodes}, then read the coding bases beside it "
            "against the member that block names.",
        ),
        expected=(
            "A table from barcode combination to cargo: each read decodes to one member of each "
            "part list, and the coding bases it carries are that member's.",
            f"The {scheme.barcode_block_length(run.barcode_length, len(positions))} bp block reads "
            f"{order}, each barcode separated from the last by the cloning scar "
            f"{scheme.cloning_scar}, at {final.block.start}-{final.block.end} of the "
            "representative construct.",
        ),
        notes=(
            "The block reads in the reverse of the order the rounds ran: each round inserted its "
            "barcode ahead of the ones already there.",
            f"{ambiguous:.1%} of the single-base deletions a barcode can carry leave a read "
            "another barcode of the same part list could leave, which no read can be assigned "
            "through.",
            (
                f"This build passes the linkage read at {fidelity:.1%} of reads carrying a "
                "barcode that still names its part. That is this run's own mark, not a "
                "published one."
                if fidelity is not None
                else "Takacsi-Nagy's Figures 1D and 1E read about 95% of their reads carrying a "
                "valid barcode at every position, and nearly 90% of the library correctly "
                "linked. That is what one source reached, not a mark this library is held to."
            ),
        ),
        troubleshooting=(
            Troubleshooting(
                "A read carries fewer barcodes than there are rounds",
                "A round did not go in. Check that round's prep against its own record before "
                "blaming the sequencing.",
            ),
            Troubleshooting(
                "A combination decodes to a member the coding bases are not",
                f"That round joined a block to the wrong barcode, so {barcodes} no longer names "
                "what the library holds. Rebuild that round rather than carrying the table "
                "forward; no later round repairs it.",
            ),
        ),
        holes=() if fidelity is not None else (stages.READ_PASS_MARK,),
    )


def _representation_step(run: Run, pair: ReadPair | None) -> Step:
    """Count which combinations the library holds and how evenly, over the barcode block alone.

    This is the read a bottleneck repeats, so it spans the block and nothing else: the forward
    anchor is the internal stuffer every member keeps, which is a method constant.
    """
    scheme, rounds = run.scheme, run.rounds
    constructs, barcodes, marks = run.constructs, run.barcodes, run.marks
    block = scheme.barcode_block_length(run.barcode_length, len(run.positions))
    return Step(
        "Read representation",
        key="read-representation",
        instructions=(
            f"Amplify across the {block} bp barcode block alone, forward from the "
            f"{len(scheme.internal_stuffer)} bp internal stuffer every member keeps and back "
            f"from the vector past the final {rounds[0].scar_overhang}"
            f"{with_pair(pair, run.read_sheet)}.",
            f"Sequence the amplicon{as_platform(pair)}, decode each read against {barcodes}, "
            "and count the reads each barcode combination gets.",
        ),
        expected=(
            f"A count for each of up to {constructs:,} distinct combinations: how many of them "
            "are seen at all, and how evenly they are read.",
            "A combination with no reads is a member the library has lost.",
            marks_sentence(constructs, marks, "of the combinations"),
        ),
        notes=(
            "Read representation again after every later bottleneck — the final assembly, and "
            "anything downstream that resamples the library. Linkage is read once; this one is "
            "read at each.",
            "The short amplicon is what makes repeating it cheap: linkage spans the whole cargo, "
            "this spans the block.",
            "Where you amplify the block to read it, carry any sample index on a primer "
            "rather than ligating it on, and keep the barcodes away from where a primer "
            "anneals: both cost more read counts than what a barcode spells does.",
        ),
        troubleshooting=(
            Troubleshooting(
                "A combination is missing",
                "It was lost at a round or at a bottleneck since, and no later step puts it "
                "back. Read the rounds' own titre plates before re-reading this one.",
            ),
            Troubleshooting(
                "The counts are heavily skewed",
                "Members differ in length and a bottleneck can favour the short ones. Imkeller's "
                "Table 2 prices the skew in screen coverage: a 90th/10th ratio of 2.5 wants "
                "200-fold, 5 wants 300-fold and 10 wants 400-fold, so a skewed library costs "
                "cells downstream rather than failing here.",
            ),
        ),
    )

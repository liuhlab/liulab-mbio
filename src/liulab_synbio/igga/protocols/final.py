"""Protocol 06: release the finished cargo and move it into the vector the application asks for.

One tube, staged: the cargo is freed and the enzymes that freed it are killed before the working
vector, the cargo enzyme and the ligase go in. A build naming no working vector still runs these
steps, because they are a stage of the method, and what the vector would have fixed is a hole.
"""

from collections.abc import Sequence

from liulab_mbio.bench.coverage import (
    RepresentationMarks,
    absent_probability,
    colonies_for_completeness,
)
from liulab_mbio.bench.goldengate import assembly_program
from liulab_mbio.bench.inactivation import heat_inactivations
from liulab_mbio.bench.steps import listed
from liulab_mbio.cloning.plan import PRODUCT_FILE
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.protocol.figures import SOURCE as FIGURE_SOURCE
from liulab_mbio.protocol.figures import SOURCE_KEY as FIGURE_SOURCE_KEY
from liulab_mbio.protocol.figures import ligation_figure
from liulab_mbio.protocol.model import (
    Figure,
    Incubation,
    Reference,
    Source,
    Stage,
    Step,
    ThermocyclerProgram,
    Troubleshooting,
    number,
)
from liulab_mbio.protocol.model import Item as Handed
from liulab_mbio.sequence import Segment, SequenceRecord
from liulab_mbio.sites import find_sites
from liulab_synbio.igga import stages
from liulab_synbio.igga.bench import (
    DIGEST_CELSIUS,
    DIGEST_SECONDS,
    GROWTH_CELSIUS,
    LIGASE,
    LIGASE_BUFFER,
    SPRI_AFTER_LIGATION,
    SPRI_BEADS,
    STRAIN,
    digest_amount,
    digest_reaction,
    final_assembly_amounts,
    final_assembly_reaction,
    growth_program,
    pulse,
    ratio,
)
from liulab_synbio.igga.method import Scheme
from liulab_synbio.igga.project import FinalAssembly
from liulab_synbio.igga.protocols.protocol import Protocol, figured, labelled
from liulab_synbio.igga.protocols.run import (
    CUVETTES,
    FINAL_SELECTIVE,
    PREP_KIT,
    WORKING_ITEM,
    Run,
    marks_sentence,
    round_equipment,
    with_pair,
)
from liulab_synbio.igga.reads import ReadPair
from liulab_synbio.igga.vector import Working, released_cargo

#: What the page is headed and what the chain names it by.
FINAL = "Final cargo ligation"


class FinalLigation(Protocol):
    """Free the cargo from the library backbone and close it into the working vector."""

    def title(self, run: Run) -> str:
        """Return the page's own heading."""
        return FINAL

    def summary(self, run: Run) -> str:
        """Return what the move comes to."""
        return (
            f"Release the finished cargo from {run.vector.name or 'the library'} and move it "
            "into the vector the application asks for."
        )

    def steps(self, run: Run) -> tuple[Step, ...]:
        """Return the five steps that move the finished library into the working vector."""
        scheme, working = run.scheme, run.working
        product = run.rounds[-1].product
        span = released_cargo(product, scheme)
        freeing = (scheme.external, *_shredders(scheme, product, span))
        reads = run.reads
        return (
            labelled(
                _pick_working_step(scheme, working, run.working_file),
                "Choose the working vector",
            ),
            labelled(_free_step(scheme, product, span, freeing), "Move the library across"),
            labelled(
                figured(
                    _assemble_step(scheme, product, span, working, run.final_assembly),
                    _junction_figure(scheme, product, span, run.records_at),
                ),
                "Move the library across",
            ),
            labelled(
                _growth_step(run.constructs, run.bench[-1].coverage.completeness, working),
                "Move the library across",
            ),
            labelled(
                _representation_step(
                    run.barcodes,
                    working,
                    None if reads is None else reads.final_representation,
                    run.constructs,
                    run.marks,
                    run.read_sheet,
                ),
                "Read the library back",
            ),
        )

    def consumes(self, run: Run) -> tuple[Handed, ...]:
        """Return the last round's prep, and the backbone the library ends in where one was named."""
        if run.working:
            return (run.prep, Handed(WORKING_ITEM, "the backbone the library ends in"))
        return (run.prep,)

    def produces(self, run: Run) -> tuple[Handed, ...]:
        """Return the finished library, pooled and ready for the screen."""
        return (run.library,)

    def shares(self, run: Run) -> tuple[str, ...]:
        """Return the rounds' reagents this move takes: the working vector, not the destination."""
        return (run.working_reagent, FINAL_SELECTIVE, PREP_KIT, CUVETTES)

    def equipment(self, run: Run) -> tuple[str, ...]:
        """Return the hardware this protocol needs, which is a round's."""
        return round_equipment(run)

    def references(self, run: Run) -> tuple[Reference, ...]:
        """Where the numbers come from, and where the scheme itself came from."""
        return run.round_references

    def sources(self, run: Run) -> dict[str, Source]:
        """Return what the method's own materials and the junction figure are cited to."""
        return dict(stages.SOURCES) | {FIGURE_SOURCE_KEY: FIGURE_SOURCE}


def _shredders(scheme: Scheme, product: SequenceRecord, span: Segment | None) -> tuple[Enzyme, ...]:
    """Return the blunt enzymes cutting the backbone the cargo leaves, read off the record."""
    if span is None:
        return ()
    return tuple(
        one
        for one in scheme.blunt
        if any(not product.covers(span, site.span) for site in find_sites(product, one))
    )


def _junction_figure(
    scheme: Scheme, product: SequenceRecord, span: Segment | None, at: str
) -> Figure | None:
    """Return the end the freed library ligates to the working vector on, at base level."""
    if span is None:
        return None
    entry = len(scheme.entry_overhang)
    return ligation_figure(
        product,
        path=at + PRODUCT_FILE,
        junction=(span.start, span.start + entry),
        enzymes=(scheme.external.name,),
        caption=(
            f"The end the library is freed on, which the working vector takes: "
            f"{scheme.external.name} leaves {scheme.entry_overhang} here."
        ),
    )


def _pick_working_step(scheme: Scheme, working: Working | None, file: str = "") -> Step:
    """Pick the vector the library moves into, which is what fixes the cargo enzyme."""
    entry, scar = scheme.entry_overhang, scheme.scar_overhang
    if working is None:
        return Step(
            "Pick the working vector",
            key="pick-working-vector",
            instructions=(
                "Take one tube of the working vector stock for the application this library is "
                "built for.",
                f"Confirm its own cargo enzyme opens it on {entry} and {scar}, giving up the "
                "ccdB cassette the library displaces.",
            ),
            expected=(
                "One opened backbone and the ccdB cassette beside it, and nothing else cut.",
            ),
            notes=(
                "The cargo enzyme is chosen against every molecule in this pot, so it has to be "
                "settled before the blocks are designed: the pipeline reads it off the working "
                "vector before it designs any block, so no block spells it.",
            ),
            holes=(stages.WORKING_VECTOR,),
        )
    record, cargo = working.record, working.enzyme
    stuffer = working.destination.stuffer
    held = f"Take one tube of {record.name or 'the working vector'} stock"
    made = (
        (f"{held}, which is the backbone the build named with the ccdB cassette already in it.",)
        if working.destination.edit is None
        else (
            f"{held}: the backbone the build named with this plan's ccdB cassette put in, "
            f"{len(record)} bp in all.",
            f"Have it made to {file or 'the record this plan wrote'}, which is the record "
            "every length below is read off. The backbone alone does not open.",
        )
    )
    return Step(
        f"Pick the working vector and confirm {cargo.name} opens it",
        key="pick-working-vector",
        instructions=(
            *made,
            f"Digest a little of it with {cargo.supplier_label} and run it on a gel.",
        ),
        expected=(
            f"{len(record)} bp opens on {entry} and {scar}, giving up its "
            f"{stuffer.end - stuffer.start} bp ccdB cassette.",
            f"Two bands and no more: {cargo.name} reads this vector nowhere else.",
        ),
        notes=(
            f"{working.cargo.check.detail}.",
            "ccdB is what makes the assembly self-selecting: a vector that took no cargo keeps "
            "the cassette and kills its host.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The stock will not open",
                "The cassette is not where the record says, or the enzyme has lost activity. "
                "Sequence the stock before assembling into it.",
            ),
        ),
    )


def _free_step(
    scheme: Scheme,
    product: SequenceRecord,
    span: Segment | None,
    enzymes: Sequence[Enzyme],
) -> Step:
    """Free the cargo from the backbone the rounds ran in, then kill what freed it."""
    if span is None:
        return Step(
            "Release the cargo from the library backbone",
            key="release-cargo",
            instructions=(
                f"Digest the finished library with {scheme.external.name} and the blunt enzyme "
                "that shreds the backbone it leaves.",
                "Heat-kill both. Nothing is purified: the working vector goes into this tube.",
            ),
            expected=(
                f"The whole cargo free, on {scheme.entry_overhang} and {scheme.scar_overhang}.",
            ),
            notes=(
                "The sites that free the cargo belong to the vector the rounds ran in, not to "
                "the cargo, so a backbone without them cannot release it.",
            ),
            holes=(stages.CARGO_RELEASE,),
        )
    named = listed([one.name for one in enzymes])
    return Step(
        f"Release the cargo with {named}",
        key="release-cargo",
        instructions=(
            f"Digest the finished library with {named} at {DIGEST_CELSIUS:g} °C.",
            "Heat-kill, then leave the tube alone: nothing is purified between the two stages.",
        ),
        tables=(
            digest_reaction(
                digest_amount((product.name or "the finished library", len(product))), enzymes
            ),
        ),
        programs=(
            ThermocyclerProgram(
                (Stage((Incubation(named, DIGEST_CELSIUS, DIGEST_SECONDS),)),),
                title=f"Digest with {named}",
            ),
            *heat_inactivations(enzymes),
        ),
        expected=(
            f"The whole {span.end - span.start} bp cargo free, on {scheme.entry_overhang} and "
            f"{scheme.scar_overhang}, out of {len(product)} bp of library.",
            f"{listed([one.name for one in enzymes[1:]]) or 'Nothing else'} cuts the backbone it "
            "came out of, so that backbone cannot close again.",
        ),
        notes=(
            "One tube, two stages. Killing the releasing enzymes before the working vector goes "
            "in is what keeps them off its backbone.",
        ),
        troubleshooting=(
            Troubleshooting(
                "The library backbone survives the assembly",
                "The blunt enzyme missed, or the heat kill was short. Both show up as colonies "
                "carrying the round's own vector rather than the working one.",
            ),
        ),
    )


def _assemble_step(
    scheme: Scheme,
    product: SequenceRecord,
    span: Segment | None,
    working: Working | None,
    sized: FinalAssembly | None,
) -> Step:
    """Join the freed cargo to the opened working vector, in the tube the release left.

    `sized` is what this build measured the vector at. Without one no published reaction sizes
    it, so the amounts are H24 rather than a figure. The cargo is never stated: the release
    fixes it, and the table computes the ratio the two meet at.
    """
    vector = "the working vector" if sized is None else f"{sized.vector_ng:g} ng of working vector"
    measured_note = (
        ()
        if sized is None
        else (
            f"{sized.vector_ng:g} ng of vector is what this run measured, not a published "
            "figure. The cargo is not measured out: the release tube goes in whole.",
        )
    )
    if working is None:
        return Step(
            "Assemble the cargo into the working vector",
            key="assemble-into-working-vector",
            instructions=(
                f"Add {vector}, its cargo enzyme and the ligase to the release tube, and "
                "run the enzyme's own Golden Gate cycling.",
            ),
            expected=(
                "One circular final vector a member, the ccdB cassette displaced by the cargo.",
            ),
            notes=(
                "This is the one reaction where the working vector meets material the rounds "
                "made; the rounds all finish first.",
                *measured_note,
            ),
            holes=(stages.WORKING_VECTOR, *(() if sized is not None else (stages.FINAL_MASSES,))),
        )
    cargo = working.enzyme
    tables = ()
    met = ""
    if sized is not None and span is not None:
        amounts = final_assembly_amounts(
            (f"{product.name or 'the library'} cargo, in the release", span.end - span.start),
            (
                f"{working.record.name or 'the working vector'}, opened",
                len(working.record) - _cassette_length(working),
            ),
            library_bp=len(product),
            vector_ng=sized.vector_ng,
        )
        tables = (final_assembly_reaction(amounts, cargo),)
        met = (
            f"The release delivers the cargo at {number(ratio(*amounts))}:1 over the vector, "
            "which is what the digest frees and not a ratio anyone sets."
        )
    joined = (
        f"about {len(working.record) - _cassette_length(working) + span.end - span.start} bp, "
        if span is not None
        else ""
    )
    return Step(
        f"Assemble the cargo into {working.record.name or 'the working vector'} with {cargo.name}",
        key="assemble-into-working-vector",
        instructions=(
            f"Add {vector}, {cargo.supplier_label} and {LIGASE} in "
            f"{LIGASE_BUFFER} to the release tube.",
            "Run the cycling below without purifying anything first.",
        ),
        tables=tables,
        programs=(assembly_program(cargo, fragments=2, library=True),),
        expected=(
            f"One circular final vector a member, {joined}joined on "
            f"{scheme.entry_overhang} and {scheme.scar_overhang}.",
            "The ccdB cassette is displaced, so a vector that took no cargo kills its host.",
        ),
        notes=(
            f"{len(product)} bp of library goes in and the cargo alone comes out: the backbone "
            "the rounds ran in is shredded and stays behind.",
            "The cycling is NEB's longer single-insert program, which it gives for library "
            "preparation rather than for cloning one gene.",
            *measured_note,
            *((met,) if met else ()),
        ),
        troubleshooting=(
            Troubleshooting(
                "Colonies that are still ccdB-positive",
                "The working vector was not opened to completion. Check the pick step's gel "
                "before repeating.",
            ),
        ),
        holes=() if sized is not None else (stages.FINAL_MASSES,),
    )


def _cassette_length(working: Working) -> int:
    """How many bases the ccdB cassette takes out of the working vector."""
    stuffer = working.destination.stuffer
    return stuffer.end - stuffer.start


def _growth_step(constructs: int, completeness: float, working: Working | None) -> Step:
    """Clean the assembly up and get all of it into cells, which is the library's last bottleneck."""
    shot = pulse()
    colonies = colonies_for_completeness(constructs, completeness)
    return Step(
        f"Clean the assembly up and electroporate into {STRAIN}",
        key="electroporate-and-grow",
        instructions=(
            f"Add {SPRI_BEADS} at {SPRI_AFTER_LIGATION:g}x the volume and elute in water, as "
            "every round did after its ligation.",
            f"Pulse at {shot.volts:g} V, {shot.ohms:g} Ω and {shot.microfarads:g} µF in a "
            f"{shot.cuvette_mm:g} mm cuvette.",
            f"Recover and grow at {GROWTH_CELSIUS:g} °C, as every round did.",
        ),
        programs=(growth_program(),),
        expected=(
            f"At least {colonies:,} net colonies: the floor for the {completeness:g} chance "
            f"this design asked for that none of its {constructs:,} distinct members is missing, "
            "equally represented.",
            f"At that count the chance a named member is missing is "
            f"{number(absent_probability(constructs, colonies))}.",
            "Near-empty plates from a no-cargo control beside it; what grows there is working "
            "vector that kept its ccdB cassette.",
        ),
        notes=(
            "This is a bottleneck like a round's, and the library can only lose members here. "
            "Electroporate all of the assembly rather than a measured part of it.",
            "The settings are the cells' own and are keyed by their catalogue number, not set "
            "by this method.",
        )
        + (
            ()
            if working is None
            else (
                f"ccdB does the selecting here, so the plates carry "
                f"{working.record.name or 'the working vector'}'s own marker.",
            )
        ),
        troubleshooting=(
            Troubleshooting(
                "Fewer net colonies than the count above",
                "The library has lost members in the transfer. Nothing downstream puts them "
                "back; repeat the assembly from more of the released cargo.",
            ),
        ),
    )


def _representation_step(
    barcodes: str,
    working: Working | None,
    pair: ReadPair | None,
    constructs: int,
    marks: RepresentationMarks,
    read_sheet: str,
) -> Step:
    """Read the library again on the other side of the move, which is the only way to size the loss."""
    where = working.record.name if working is not None else "the working vector"
    return Step(
        "Read representation in the final vector",
        key="read-final-representation",
        instructions=(
            f"Amplify across the barcode block again{with_pair(pair, read_sheet)}.",
            f"Sequence, decode each read against {barcodes}, and compare the counts with the "
            "read taken in the library backbone.",
        ),
        expected=(
            "The same combinations, at a similar evenness. A combination seen before the move "
            f"and not after it was lost in the transfer into {where or 'the working vector'}.",
            marks_sentence(constructs, marks, "of what survived the move"),
        ),
        notes=(
            "Linkage is read once, in the library backbone; a barcode still names the same part "
            "after the move, because the move carries the whole cargo in one piece.",
            "The forward anchor is the same retained internal stuffer, which travels with the "
            "cargo; the reverse anchor moved with the vector, so this pair is designed against "
            "the final record rather than reused from the read before the move.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Markedly fewer combinations than before the move",
                "The transfer was the bottleneck, not the rounds. Compare the colony count with "
                "the step above before rebuilding anything.",
            ),
        ),
    )

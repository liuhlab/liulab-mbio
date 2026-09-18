"""The restriction and ligation bench protocol: its own steps, its own notes, and their order.

The steps any bench shares are `liulab_mbio.bench.steps`. This module runs them around the two
digests, the gel purification and the ligation, and adds what only this method has to say. Every
number is computed by `liulab_mbio.cloning.restriction.plan` or by the modules it calls, and
`liulab_mbio.cloning.restriction.bench` is where each one's source is written down.

Two things this method says and Golden Gate does not. Its junction is not scarless: the two ends
came from a recognition site and ligating them puts that site back, so the protocol names what
each junction spells. And no supplier states a colony count for it, so the plate is described by
NEB's four controls and their ratios rather than by a number.
"""

from collections.abc import Mapping, Sequence

from liulab_mbio import checks as judged
from liulab_mbio.bench import REFERENCES as BENCH_REFERENCES
from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.gels import agarose_percent, choose_ladder
from liulab_mbio.bench.inactivation import heat_inactivation
from liulab_mbio.bench.oligos import oligo_row
from liulab_mbio.bench.pcr import COLONY_PCR_MASTER_MIX, colony_pcr_master_mix_component
from liulab_mbio.bench.phenotype import Phenotype
from liulab_mbio.bench.steps import (
    CELLS_UL,
    COLONY_PCR_TITLE,
    HEAT_SHOCK_CELSIUS,
    IPTG_UM,
    OUTGROWTH_CELSIUS,
    OUTGROWTH_UL,
    PLATE_REFERENCE,
    SEQUENCING_TITLE,
    XGAL_UG_ML,
    badges,
    card,
    catalogued,
    colony_pcr_step,
    enzyme_material,
    listed,
    phenotype_sentences,
    quantify_step,
    sequencing_step,
    transform_step,
)
from liulab_mbio.bench.validation import ColonyCheck, SangerRead
from liulab_mbio.cloning.restriction.bench import (
    BLUNT_SECONDS,
    BUFFER_FINDER,
    COHESIVE_SECONDS,
    CONTROLS,
    DIGEST_SECONDS,
    LIGATION_NG_UL,
    MAX_DNA_FRACTION,
    OVERNIGHT_CELSIUS,
    REFERENCES,
    ROOM_CELSIUS,
    STAR_ACTIVITY,
    TRANSFORM_UL,
    digest_amount,
    digest_reaction,
    gel_recovery,
    ligation_program,
    ligation_reaction,
    shared_buffer,
)
from liulab_mbio.cloning.restriction.digest import Diagnostic, Piece
from liulab_mbio.cloning.restriction.ligation import Junction, Ligation
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.primers.evaluation import PrimerReport
from liulab_mbio.primers.thresholds import PrimerRole, Thresholds
from liulab_mbio.protocol.model import (
    Gel,
    Ladder,
    Lane,
    Material,
    Protocol,
    Reference,
    Step,
    Timer,
    Troubleshooting,
)
from liulab_mbio.sequence import SequenceRecord

#: Who sells the products this protocol names. Every catalogue number it prints comes out of an
#: enzyme record or a product name a supplier wrote; none is written here.
SUPPLIER = "New England Biolabs"

#: The strain a protocol names unless the caller picks one. Blue/white screening needs a host
#: that supplies the rest of the lacZ fragment the vector carries, which this one does.
DEFAULT_HOST = "NEB 5-alpha Competent E. coli (C2987)"

#: The hardware a run needs, which no reagent table covers.
EQUIPMENT: tuple[str, ...] = (
    "Thermocycler with a heated lid",
    "Agarose gel rig, power supply and a gel imager",
    "Scalpel or gel cutting tips, and a clean cutting surface",
    "Microcentrifuge",
    "Spectrophotometer or fluorometer",
    f"Heat block or water bath at {HEAT_SHOCK_CELSIUS:g} °C",
    f"Shaking incubator and a plate incubator at {OUTGROWTH_CELSIUS:g} °C",
)


def protocol(
    *,
    vector: SequenceRecord,
    source: SequenceRecord,
    enzymes: Sequence[Enzyme],
    vector_pieces: Sequence[Piece],
    source_pieces: Sequence[Piece],
    ligation: Ligation,
    digests: Sequence[Amount],
    amounts: Sequence[Amount],
    diagnostic: Diagnostic,
    colony: ColonyCheck,
    reads: Sequence[SangerRead],
    read_reports: Sequence[PrimerReport],
    phenotype: Phenotype,
    checks: Sequence[judged.Check],
    host: str,
    thresholds: Mapping[PrimerRole, Thresholds],
) -> Protocol:
    """Return the bench protocol for one planned cloning, ready to render.

    Each argument is the `liulab_mbio.cloning.restriction.plan.Plan` field or property of that
    name. The steps run in the order someone does them: the two digests, the gel that separates
    the pieces and the extraction that recovers them, quantification, the ligation, the
    transformation and plating, the colony PCR and the sequencing.
    """
    backbone, insert = ligation.pieces
    named = listed([enzyme.name for enzyme in enzymes])
    return Protocol(
        f"Restriction and ligation: {insert.name} into {vector.name}",
        summary=(
            f"Cut {vector.name} and {source.name} with {named}, gel-purify the "
            f"{backbone.length} bp backbone and the {insert.length} bp insert, ligate them, and "
            "confirm the clone by colony PCR and sequencing."
        ),
        overview=_overview(vector, source, enzymes, ligation, phenotype),
        highlights=_highlights(ligation, phenotype),
        checks=badges(checks),
        materials=_materials(
            vector=vector,
            source=source,
            enzymes=enzymes,
            colony=colony,
            host=host,
            phenotype=phenotype,
            sizes=_sizes(vector_pieces, source_pieces),
        ),
        oligos=(
            *(
                oligo_row(report, purpose=COLONY_PCR_TITLE, thresholds=thresholds["colony PCR"])
                for report in colony.reports
            ),
            *(
                oligo_row(report, purpose=SEQUENCING_TITLE, thresholds=thresholds["sequencing"])
                for report in read_reports
            ),
        ),
        equipment=EQUIPMENT,
        steps=_steps(
            vector=vector,
            source=source,
            enzymes=enzymes,
            vector_pieces=vector_pieces,
            source_pieces=source_pieces,
            ligation=ligation,
            digests=digests,
            amounts=amounts,
            diagnostic=diagnostic,
            colony=colony,
            reads=reads,
            phenotype=phenotype,
            host=host,
        ),
        references=_references(phenotype),
    )


def _overview(
    vector: SequenceRecord,
    source: SequenceRecord,
    enzymes: Sequence[Enzyme],
    ligation: Ligation,
    phenotype: Phenotype,
) -> dict[str, str]:
    """Return the facts to check before starting, each short enough to be a card."""
    backbone, insert = ligation.pieces
    facts = {
        "Vector": f"{vector.name}, {len(vector)} bp",
        "Insert": f"{insert.name}, {insert.length} bp from {source.name}",
        "Enzymes": card(
            listed([enzyme.supplier_label for enzyme in enzymes]), f"{len(enzymes)} enzymes"
        ),
        "Backbone": f"{backbone.length} bp",
        "Junctions": card(
            listed([one.label for one in ligation.junctions]),
            f"{len(ligation.junctions)} junctions",
        ),
        "Product": f"{ligation.product.name}, {len(ligation.product)} bp",
    }
    selection = _selection(phenotype)
    if selection:
        facts["Selection"] = selection
    return facts


def _selection(phenotype: Phenotype) -> str:
    """Return what to select transformants on, or nothing where the vector annotates no marker."""
    if phenotype.antibiotic:
        return phenotype.antibiotic
    return f"{phenotype.marker.name} marker" if phenotype.marker is not None else ""


def _highlights(ligation: Ligation, phenotype: Phenotype) -> tuple[str, ...]:
    """Return what the facts mean, a sentence each: the junctions first, then the phenotype."""
    backbone, insert = ligation.pieces
    return (
        f"One ligation joins two fragments: {backbone.name} ({backbone.length} bp) and "
        f"{insert.name} ({insert.length} bp).",
        f"This method's junction is not scarless: {_spelled(ligation.junctions)}",
        *phenotype_sentences(phenotype, [insert.name]),
    )


def _spelled(junctions: Sequence[Junction]) -> str:
    """Say what each junction spells, and that the product gains those bases."""
    said = listed(
        [
            f"{one.spells} at {one.start} ({one.enzyme})"
            if one.enzyme
            else f"{one.label} at {one.start}, which spells no enzyme's site"
            for one in junctions
        ]
    )
    return (
        "the two ends came from a recognition site and ligating them puts it back, so the "
        f"product reads {said}."
    )


def _sizes(*digests: Sequence[Piece]) -> tuple[int, ...]:
    """Every band the two digests give, which is what a ladder and a gel are chosen for."""
    return tuple(piece.length for pieces in digests for piece in pieces)


def _materials(
    *,
    vector: SequenceRecord,
    source: SequenceRecord,
    enzymes: Sequence[Enzyme],
    colony: ColonyCheck,
    host: str,
    phenotype: Phenotype,
    sizes: tuple[int, ...],
) -> tuple[Material, ...]:
    """Every reagent and consumable the protocol asks for. The oligos are the order sheet."""
    ladders = dict.fromkeys((choose_ladder(sizes).name, colony.ladder.name))
    return (
        Material(f"{vector.name} plasmid", storage="-20 °C", note="cut to make the backbone"),
        Material(f"{source.name} plasmid", storage="-20 °C", note="cut to release the insert"),
        *(
            enzyme_material(enzyme, note=f"cuts both plasmids once, leaving {_left(enzyme)}")
            for enzyme in enzymes
        ),
        _buffer_material(enzymes),
        Material(
            "T4 DNA Ligase and its 10X reaction buffer",
            supplier=SUPPLIER,
            catalog="M0202",
            storage="-20 °C",
            note="thaw a fresh buffer aliquot: its ATP goes over freeze-thaws",
        ),
        Material("Agarose, 1X TAE or TBE, and a DNA stain"),
        *(catalogued(name, supplier=SUPPLIER) for name in ladders),
        Material("Gel extraction spin columns", supplier=SUPPLIER, catalog="T1120"),
        Material(
            host,
            supplier=SUPPLIER if host == DEFAULT_HOST else "",
            storage="-80 °C",
            amount=f"{CELLS_UL:g} µL per transformation",
        ),
        Material(
            "SOC or NEB 10-beta/Stable Outgrowth Medium",
            amount=f"{OUTGROWTH_UL:g} µL per transformation",
        ),
        Material(_plate(phenotype), amount="one plate per transformation"),
        catalogued(
            COLONY_PCR_MASTER_MIX,
            supplier=SUPPLIER,
            storage="-20 °C",
            amount=f"{colony_pcr_master_mix_component().volume_ul:g} µL per reaction",
        ),
    )


def _left(enzyme: Enzyme) -> str:
    """Say what kind of end an enzyme leaves, for a material's note."""
    return "a blunt end" if enzyme.end == "blunt" else f"a {enzyme.end} overhang"


def _buffer_material(enzymes: Sequence[Enzyme]) -> Material:
    """Return the digest's buffer as a material, named only where the records name one."""
    buffer = shared_buffer(enzymes)
    if buffer is None:
        return Material(
            "A buffer both enzymes work in",
            supplier=SUPPLIER,
            storage="-20 °C",
            note=f"the shipped records do not agree on one; look the pair up in {BUFFER_FINDER}",
        )
    return Material(
        f"{buffer}, 10X",
        supplier=SUPPLIER,
        storage="-20 °C",
        note="both enzymes are supplied in it, so they digest together in it",
    )


def _plate(phenotype: Phenotype) -> str:
    """Return what to pour the selection plates with."""
    antibiotic = phenotype.antibiotic or "the vector's own antibiotic"
    if phenotype.blue_white:
        return f"LB agar plates with {antibiotic}, {XGAL_UG_ML} µg/mL X-gal and {IPTG_UM} µM IPTG"
    return f"LB agar plates with {antibiotic}"


def _steps(
    *,
    vector: SequenceRecord,
    source: SequenceRecord,
    enzymes: Sequence[Enzyme],
    vector_pieces: Sequence[Piece],
    source_pieces: Sequence[Piece],
    ligation: Ligation,
    digests: Sequence[Amount],
    amounts: Sequence[Amount],
    diagnostic: Diagnostic,
    colony: ColonyCheck,
    reads: Sequence[SangerRead],
    phenotype: Phenotype,
    host: str,
) -> tuple[Step, ...]:
    """Return the steps in the order they happen, the shared ones carrying this method's notes."""
    backbone, insert = ligation.pieces
    return (
        _digest_step(vector, enzymes, vector_pieces, digests[0], keeping=backbone),
        _digest_step(source, enzymes, source_pieces, digests[1], keeping=insert),
        _purify_step(vector, source, vector_pieces, source_pieces, keeping=(backbone, insert)),
        quantify_step(amounts),
        _ligation_step(ligation, amounts),
        transform_step(
            host,
            phenotype,
            inserts=[insert.name],
            colonies="No supplier states a colony count for this method, so run the controls "
            "below and read the plate against them rather than against a number.",
            expected=CONTROLS,
        ),
        colony_pcr_step(colony, junctions=len(ligation.junctions)),
        _diagnostic_step(diagnostic, product=ligation.product),
        sequencing_step(
            reads, junctions=[one.label for one in ligation.junctions], inserts=[insert.name]
        ),
    )


def _diagnostic_step(diagnostic: Diagnostic, *, product: SequenceRecord) -> Step:
    """Cut a miniprep with the cloning pair, which drops the insert back out of a correct clone."""
    clone, empty = diagnostic.names
    named = listed([enzyme.name for enzyme in diagnostic.enzymes])
    ladder = choose_ladder(diagnostic.bands)
    percent = agarose_percent(diagnostic.bands)
    amount = digest_amount((clone, len(product)))
    return Step(
        f"Check a miniprep by digesting it with {named}",
        instructions=(
            "Miniprep two or three of the colonies the PCR called correct.",
            f"Mix the reaction below, {amount.nanograms:g} ng of miniprep first.",
            f"Incubate at {_celsius(diagnostic.enzymes)} for {DIGEST_SECONDS // 60} minutes.",
            f"Run the whole digest on a {percent:g}% agarose gel beside the ladder.",
        ),
        tables=(digest_reaction(diagnostic.enzymes, amount, title="Diagnostic digest"),),
        timers=(Timer("Diagnostic digest", DIGEST_SECONDS),),
        gels=(
            Gel(
                ladder,
                (Lane(clone, diagnostic.clone), Lane(empty, diagnostic.empty)),
                title="Diagnostic digest",
            ),
        ),
        expected=(
            f"A correct clone: {_bands(diagnostic.clone)}.",
            f"A colony carrying {empty} instead: {_bands(diagnostic.empty)}.",
            *_off_the_gel(diagnostic.bands, ladder),
        ),
        notes=(
            "The junctions put both recognition sites back, so the pair that made the clone is "
            "what cuts the insert out of it again.",
            "Sequence only a miniprep that gives the clone's bands.",
        ),
        troubleshooting=(
            Troubleshooting(
                "Every miniprep gives the vector's bands",
                "The background is uncut or religated vector; run the controls under the "
                "transformation to find which.",
            ),
            Troubleshooting(
                "One band, at the plasmid's full length",
                "Only one site cut. Check the enzymes' methylation sensitivity against the "
                "strain the miniprep was grown in.",
            ),
        ),
    )


def _digest_step(
    record: SequenceRecord,
    enzymes: Sequence[Enzyme],
    pieces: Sequence[Piece],
    amount: Amount,
    *,
    keeping: Piece,
) -> Step:
    """Cut one plasmid with both enzymes in one tube."""
    named = listed([enzyme.name for enzyme in enzymes])
    room = f"{MAX_DNA_FRACTION:.0%}"
    return Step(
        f"Digest {record.name} with {named}",
        instructions=(
            f"Mix the reaction below, {amount.nanograms:g} ng of {record.name} first.",
            f"Incubate at {_celsius(enzymes)} for {DIGEST_SECONDS // 60} minutes.",
        ),
        tables=(digest_reaction(enzymes, amount, title=f"{record.name} digest"),),
        timers=(Timer(f"{record.name} digest", DIGEST_SECONDS),),
        expected=tuple(
            f"{piece.name}: {piece.length} bp"
            + (
                f", which is what goes on ({piece.left_enzyme.name} and "
                f"{piece.right_enzyme.name} ends)."
                if piece is keeping
                else "."
            )
            for piece in pieces
        ),
        notes=(
            _buffer_note(enzymes),
            f"Keep the DNA solution under {room} of the reaction; a column eluate carries salt, "
            "and salt leaves the digest incomplete.",
            *STAR_ACTIVITY,
        ),
        troubleshooting=(
            Troubleshooting(
                "An uncut band remains",
                "Add more units or incubate longer, and check the enzymes' methylation "
                "sensitivity against the strain the plasmid was grown in.",
            ),
            Troubleshooting(
                "Extra bands",
                "Star activity: use fewer units, a shorter incubation and the supplied buffer.",
            ),
        ),
    )


def _buffer_note(enzymes: Sequence[Enzyme]) -> str:
    """Say which buffer the digest runs in, or that nothing sourced here can say."""
    buffer = shared_buffer(enzymes)
    if buffer is None:
        return (
            "Nothing sourced here says these enzymes share a buffer, so this carries no "
            f"verdict rather than a pass: look the pair up in {BUFFER_FINDER} before putting "
            "both in one tube."
        )
    return (
        f"Both enzymes are supplied in {buffer}, which is NEB's own rule for digesting two of "
        "them together."
    )


def _celsius(enzymes: Sequence[Enzyme]) -> str:
    """Say what to incubate a digest at, or name the temperatures where the records disagree."""
    wanted = sorted(
        {enzyme.incubation_celsius for enzyme in enzymes if enzyme.incubation_celsius is not None}
    )
    if len(wanted) == 1:
        return f"{wanted[0]} °C"
    named = listed([f"{one} °C" for one in wanted]) or "a temperature their records do not state"
    return f"{named} -- these enzymes want different ones, which NEB answers with two digests"


def _purify_step(
    vector: SequenceRecord,
    source: SequenceRecord,
    vector_pieces: Sequence[Piece],
    source_pieces: Sequence[Piece],
    *,
    keeping: tuple[Piece, Piece],
) -> Step:
    """Run both digests out, cut the two bands that go on, and recover them."""
    backbone, insert = keeping
    sizes = _sizes(vector_pieces, source_pieces)
    ladder = choose_ladder(sizes)
    percent = agarose_percent(sizes)
    low, high = gel_recovery(min(backbone.length, insert.length))
    return Step(
        "Separate the digests on a gel and recover the two fragments",
        instructions=(
            f"Pour a {percent:g}% agarose gel and load each whole digest beside the ladder.",
            "Run until the bands below are apart.",
            f"Cut out the {backbone.length} bp {backbone.name} band and the "
            f"{insert.length} bp {insert.name} band.",
            "Recover each slice on a spin column and elute in the smallest volume the kit allows.",
        ),
        cautions=(
            "Keep the gel off the ultraviolet box for as long as you can; use blue light where "
            "there is one.",
        ),
        gels=(
            Gel(
                ladder,
                (
                    Lane(vector.name, tuple(piece.length for piece in vector_pieces)),
                    Lane(source.name, tuple(piece.length for piece in source_pieces)),
                ),
                title="Digests",
            ),
        ),
        expected=(
            f"{vector.name}: {_bands(_sizes(vector_pieces))}.",
            f"{source.name}: {_bands(_sizes(source_pieces))}.",
            f"Each slice recovers {low:.0%} to {high:.0%} of the DNA that was in it.",
            *_off_the_gel(sizes, ladder),
        ),
        notes=(
            "A gel is what separates the two pieces of each digest from one another, so a piece "
            "that is not wanted cannot religate into the backbone.",
            "It also takes the enzymes away, so nothing has to be heat inactivated before the "
            f"ligation.{_heat(vector_pieces, source_pieces)}",
        ),
        troubleshooting=(
            Troubleshooting(
                "Two bands did not separate",
                "Run the gel further, or pour it at a percentage that resolves that size range.",
            ),
            Troubleshooting(
                "Little DNA comes off the column",
                "Elute twice through the same column, and keep the agarose percentage low.",
            ),
        ),
    )


def _off_the_gel(sizes: tuple[int, ...], ladder: Ladder) -> tuple[str, ...]:
    """Say which bands run off the bottom, where any is smaller than the ladder measures."""
    smallest = min(ladder.bands_bp)
    lost = sorted({size for size in sizes if size < smallest})
    if not lost:
        return ()
    named = listed([f"{size} bp" for size in lost])
    return (f"{named} runs past the bottom of this ladder.",)


def _bands(sizes: Sequence[int]) -> str:
    """Name the bands one lane gives, largest first as a gel reads them."""
    return listed([f"{size} bp" for size in sorted(sizes, reverse=True)])


def _heat(*digests: Sequence[Piece]) -> str:
    """Name the enzymes heat cannot stop, where any of them is used, or say nothing."""
    enzymes = {
        enzyme.name: enzyme
        for pieces in digests
        for piece in pieces
        for enzyme in (piece.left_enzyme, piece.right_enzyme)
    }
    unstoppable = [
        name for name, enzyme in sorted(enzymes.items()) if heat_inactivation(enzyme) is None
    ]
    if not unstoppable:
        return ""
    return (
        f" That matters here: the supplier states no heat inactivation for "
        f"{listed(unstoppable)}, so a gel or a column is the only way to stop "
        f"{'them' if len(unstoppable) > 1 else 'it'}."
    )


def _ligation_step(ligation: Ligation, amounts: Sequence[Amount]) -> Step:
    """Set the ligation up, run it, and stop the ligase."""
    low, high = TRANSFORM_UL
    floor, ceiling = LIGATION_NG_UL
    seconds = BLUNT_SECONDS if ligation.blunt else COHESIVE_SECONDS
    return Step(
        "Ligate the insert into the backbone",
        instructions=(
            "Mix the reaction below, the DNA first and the ligase last.",
            f"Hold at {ROOM_CELSIUS:g} °C for {seconds // 60} minutes, or overnight at "
            f"{OVERNIGHT_CELSIUS:g} °C.",
            "Then run the heat step below, and chill on ice.",
        ),
        tables=(ligation_reaction(amounts),),
        programs=(ligation_program(blunt=ligation.blunt),),
        expected=(
            *(
                f"The product carries {one.label} at {one.start}, joining {one.before} to "
                f"{one.after}."
                for one in ligation.junctions
            ),
            f"{low:g} to {high:g} µL of this goes into the cells; the rest keeps at -20 °C.",
        ),
        notes=(
            "Picomoles, not nanograms: the table asks for a molar ratio, and the shorter "
            "fragment weighs less at the same ratio.",
            f"Keep the two fragments together at {floor:g} to {ceiling:g} ng/µL. Below that a "
            "fragment closes on itself instead of joining its partner.",
            "Junction positions are 0-based, on the product.",
        ),
        troubleshooting=(
            Troubleshooting(
                "No colonies later",
                "At least one fragment has to carry a 5' phosphate; vary the ratio, and use a "
                "fresh buffer aliquot, because its ATP goes over freeze-thaws.",
            ),
            Troubleshooting(
                "Empty vector on the plate",
                "The backbone band carried some uncut plasmid; run the gel further and cut the "
                "band clean.",
            ),
        ),
    )


def _references(phenotype: Phenotype) -> tuple[Reference, ...]:
    """Where the numbers come from."""
    items = [*REFERENCES, *BENCH_REFERENCES]
    if phenotype.blue_white:
        items.append(PLATE_REFERENCE)
    return tuple(items)

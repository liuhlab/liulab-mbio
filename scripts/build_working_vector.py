#!/usr/bin/env python3
"""Build the demo's working vector, and the bench protocol that makes it.

    pixi run python scripts/build_working_vector.py [--from PATH] [--out DIR]

The AP-1 example ends by moving its library into a working vector. The one it names is
pLVX-TetOne-Puro-GFP, Addgene 171123, with every BsaI and BsmBI site taken out, so the demo
shows the domestication instead of shipping a backbone someone cleaned by hand. The record it
writes is what the library plan then reads, and the protocol beside it is its own sitting: the
library run opens with this record already in hand and has no domestication step.

Three of the eight sites lie in a coding sequence and `liulab_mbio.sites.domesticate` takes
those out by itself. The other five lie in none, so this script states the base it changes at
each one. `docs/research/working-vector-plvx-tetone.md` section 7 prices all eight, and
`docs/research/domestication-methods.md` section 3.2 is why the route is to order the plasmid
whole: that note's own worked case leaves the two long terminal repeat sites standing, and once
#485 gave them a precedent and they came into the edit set, no oligo-based route reaches them.
Every site is counted again on the result, which is what holds the arithmetic honest.

The parent ships as GenBank in `tests/data/`, converted from Addgene's own render: a file a
shipped record is built from is an input, not a reference document. The script never reaches
the network.

A DNA sequence carries no licence. The parent is deposited and the domesticated record is the
lab's own molecule, which is why it ships.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path

from liulab_mbio.bench.steps import enzyme_material
from liulab_mbio.edits import replace
from liulab_mbio.enzymes import Enzyme, enzymes, get_enzyme
from liulab_mbio.io import read_record
from liulab_mbio.protocol.model import (
    Check,
    Hole,
    Item,
    Protocol,
    Reference,
    Step,
    read_protocol,
    write_protocol,
)
from liulab_mbio.protocol.render import minted, write_html
from liulab_mbio.sequence import Feature, Segment, SequenceRecord
from liulab_mbio.sites import DomesticationReport, domesticate, find_sites
from liulab_mbio.translate import translate

REPO = Path(__file__).resolve().parents[1]

#: pLVX-TetOne-Puro-GFP, Addgene 171123, read from the depositor's own render. 9,895 bp circular.
PARENT = REPO / "tests/data/pLVX-TetOne-Puro-GFP.gb"

#: Where the record and its protocol go: the example the library plan is run against.
OUT = REPO / "docs/examples/ap1-library"

#: What this writes, as the example's page and `scripts/check_examples.py` name them.
RECORD_FILE = "working-vector.gb"
PROTOCOL_DATA_FILE = "working-vector-domestication.json"
PROTOCOL_FILE = "working-vector-domestication.html"

#: The two enzymes the record is cleared of: the one a round releases a part with, and the one
#: the final transfer into this vector uses.
CLEARED = ("BsaI", "BsmBI")

NAME = "pLVX-TetOne-dom"


@dataclass(frozen=True)
class Edit:
    """One base this build changes where the site lies in no coding sequence.

    `liulab_mbio.sites.domesticate` hands these back untouched, because taking a site out of a
    non-coding element changes what the record spells and the package leaves that to its caller.
    These five are the demo's choice.
    """

    position: int
    before: str
    after: str
    enzyme: str
    element: str


#: The five bases this build changes by hand, in the order they lie on the plasmid. Section 7 of
#: `docs/research/working-vector-plvx-tetone.md` prices every one of them. The two long terminal
#: repeat copies carry the same change at the same offset, which is the one condition the
#: published precedent rides on.
EDITS = (
    Edit(457, "C", "T", "BsaI", "5' long terminal repeat, transcript position +4"),
    Edit(3725, "G", "A", "BsaI", "hPGK promoter"),
    Edit(3881, "C", "G", "BsmBI", "hPGK promoter"),
    Edit(6474, "G", "A", "BsaI", "the linker joining the WPRE to the HIV-1 remnant"),
    Edit(7115, "C", "T", "BsaI", "3' long terminal repeat, transcript position +4"),
)

#: What the two holes this record leaves are called. Neither is a number a bench page may
#: invent: a precedent is not a measurement.
LTR_HOLE = Hole(
    "H33",
    "no titre for a lentiviral vector whose transcript start has been changed",
    "undecided",
    where="the change at transcript position +4 in both long terminal repeats",
    filled_by="one packaging run titring this vector beside the parent it was built from",
)
PROMOTER_HOLE = Hole(
    "H34",
    "no expression figure for either changed base of the hPGK promoter",
    "undecided",
    where="the two changes in the promoter driving the transactivator",
    filled_by="a reporter assay reading this promoter beside the parent's",
)

REFERENCES = (
    Reference(
        "Addgene 171123, pLVX-TetOne-Puro-GFP: the parent plasmid and its verified sequence.",
        url="https://www.addgene.org/171123/",
    ),
    Reference(
        "Haellman, V. et al. (2021) Mammalian synthetic biology platform VAMSyB. "
        "Metabolic Engineering 66, 41-50. The vector carrying the same change in both repeats.",
        url="https://doi.org/10.1016/j.ymben.2021.04.003",
    ),
    Reference(
        "Peterman, N. et al. (2025) Nucleic Acids Research 53, gkaf528. The hPGK promoter of "
        "Addgene 239691 carries the same change at the same base.",
        url="https://doi.org/10.1093/nar/gkaf528",
    ),
    Reference(
        "New England Biolabs, Protocol for a single restriction enzyme digest: 10 units an "
        "hour for 1 microgram of DNA.",
        url="https://www.neb.com/en-us/protocols/0001/01/01/optimizing-restriction-enzyme-reactions",
    ),
    Reference(
        "Twist Bioscience, gene synthesis sequence acceptance criteria: a direct repeat longer "
        "than 200 bases is a high-complexity sequence.",
        url="https://www.twistbioscience.com/faq/gene-synthesis",
    ),
)


def rebuild(parent: SequenceRecord) -> tuple[SequenceRecord, DomesticationReport]:
    """Return the parent with no BsaI and no BsmBI site left, and what the package changed itself.

    The synonymous changes go first, because they are the ones a codon table decides rather than
    a person. Each stated base is then checked against the record before it is changed, and the
    sites are counted again at the end: the edit arithmetic is never trusted.

    Raises
    ------
    ValueError
        If a stated position lies outside the record, reads another base, or leaves a site
        standing.
    """
    counted = {name: len(find_sites(parent, get_enzyme(name))) for name in CLEARED}
    record, report = domesticate(parent, CLEARED)
    for edit in EDITS:
        record = _changed(record, edit)
    left = find_sites(record, [get_enzyme(name) for name in CLEARED])
    if left:
        where = ", ".join(f"{one.enzyme.name} at {one.start}" for one in left)
        raise ValueError(f"the rebuilt record still reads {where}")
    return _annotated(record, report, counted), report


def _changed(record: SequenceRecord, edit: Edit) -> SequenceRecord:
    """Change one stated base, refusing a position this record does not hold as stated.

    The span is checked against the record's own length rather than through
    `SequenceRecord.fits`, which admits an end past the length because that is how
    `docs/adr/0001-coordinates.md` spells a span across the origin.
    """
    if not 0 <= edit.position < len(record):
        raise ValueError(
            f"{edit.enzyme} at {edit.position}: the record is {len(record)} bases long"
        )
    was = record.sequence[edit.position]
    if was != edit.before:
        raise ValueError(
            f"position {edit.position} reads {was!r} and this build states {edit.before!r}: "
            "the parent is not the one these bases were measured on"
        )
    changed, _ = replace(record, edit.position, edit.position + 1, edit.after)
    return changed


def _annotated(
    record: SequenceRecord, report: DomesticationReport, counted: dict[str, int]
) -> SequenceRecord:
    """Name the record, and mark every base the build changed, whoever changed it."""
    coding = tuple((one.position, one.site.enzyme.name) for one in report.changes for one in (one,))
    stated = tuple((edit.position, edit.enzyme) for edit in EDITS)
    marks = tuple(
        Feature(
            f"{enzyme} site removed",
            "misc_difference",
            (Segment(position, position + 1),),
        )
        for position, enzyme in sorted((*coding, *stated))
    )
    return SequenceRecord(
        record.sequence,
        topology=record.topology,
        name=NAME,
        features=(*record.features, *marks),
        primers=record.primers,
        notes={"Description": _described(counted, len(report.changes))},
    )


def _described(counted: dict[str, int], synonymous: int) -> str:
    """Say what the record is, counting its sites off the parent rather than stating them."""
    taken = ", ".join(f"{count} {name}" for name, count in counted.items())
    return (
        f"pLVX-TetOne-Puro-GFP (Addgene 171123) with its {taken} sites taken out: "
        f"{synonymous} by a synonymous codon, {len(EDITS)} by a stated base. Built by "
        "scripts/build_working_vector.py"
    )


def unchanged_proteins(
    parent: SequenceRecord, record: SequenceRecord, report: DomesticationReport
) -> tuple[str, ...]:
    """Return the coding sequences a synonymous change moved a codon in, each still reading alike.

    Every one is translated on both records and compared, so what the page says about the
    protein is a measurement of the file it ships beside.

    Raises
    ------
    ValueError
        If any of them no longer reads as the parent reads it.
    """
    named = []
    for change in report.changes:
        feature = change.feature
        found = next((one for one in record.features if one.name == feature.name), None)
        if found is None:
            raise ValueError(f"the rebuilt record annotates no {feature.name!r}")
        if translate(record.extract(found)) != translate(parent.extract(feature)):
            raise ValueError(f"{feature.name} no longer reads as the parent reads it")
        named.append(feature.name)
    return tuple(dict.fromkeys(named))


def single_cutter(record: SequenceRecord) -> Enzyme:
    """Return a shipped enzyme cutting `record` exactly once, for the digest's own control.

    A lane that has to run is what tells an enzyme that did not cut from an enzyme that was not
    added, so the control is read off the record rather than named by hand.

    Raises
    ------
    ValueError
        If no shipped enzyme cuts this record exactly once.
    """
    for enzyme in sorted(enzymes(), key=lambda one: one.name):
        if len(find_sites(record, enzyme)) == 1:
            return enzyme
    raise ValueError("no shipped enzyme cuts this record once, so the digest has no control")


def protocol(
    parent: SequenceRecord, record: SequenceRecord, report: DomesticationReport
) -> Protocol:
    """Return the bench protocol that turns `parent` into `record`.

    Every number it states is counted off the two records here, so the page and the file beside
    it cannot disagree.
    """
    counted = {name: len(find_sites(parent, get_enzyme(name))) for name in CLEARED}
    still = {
        name: len(find_sites(record, get_enzyme(name)))
        for name in ("BbsI", "SapI")
        if find_sites(record, get_enzyme(name))
    }
    control = single_cutter(record)
    proteins = unchanged_proteins(parent, record, report)
    length = f"{len(record):,}"
    return Protocol(
        "Domesticating the working vector",
        summary=(
            f"The library's working vector is pLVX-TetOne-Puro-GFP, which reads "
            f"{counted['BsaI']} BsaI and {counted['BsmBI']} BsmBI sites. This orders the same "
            f"plasmid with all {counted['BsaI'] + counted['BsmBI']} taken out and confirms what "
            "arrives. It is a sitting of its own: the library run opens with the record it "
            "leaves behind."
        ),
        overview={
            "Parent": "pLVX-TetOne-Puro-GFP, Addgene 171123",
            "Plasmid": f"{length} bp circular",
            "Sites taken out": f"{counted['BsaI']} BsaI, {counted['BsmBI']} BsmBI",
            "Bases changed": f"{len(report.changes) + len(EDITS)}",
            "Route": "the whole plasmid, ordered already changed",
            "Protein changed": "none",
        },
        highlights=_highlights(report, still, proteins),
        checks=_checks(record, report, still, proteins),
        consumes=(
            Item(
                "pLVX-TetOne-Puro-GFP",
                "The parent plasmid, ordered from Addgene as deposit 171123",
                spec=(f"circular, {length} bp",),
                storage="-20 °C",
            ),
        ),
        produces=(
            Item(
                "domesticated working vector",
                "The same plasmid with no BsaI and no BsmBI site, which the library run opens with",
                spec=(f"circular, {length} bp", "no BsaI or BsmBI site"),
                storage="-20 °C",
            ),
        ),
        materials=(
            enzyme_material(
                get_enzyme("BsaI"), amount="10 units a digest", note="reads nothing in the result"
            ),
            enzyme_material(
                get_enzyme("BsmBI"),
                amount="10 units a digest",
                note="reads nothing in the result",
            ),
            enzyme_material(control, amount="10 units a digest", note="the control lane, one cut"),
        ),
        equipment=(
            "Heat block or thermocycler holding 37 °C and 55 °C",
            "Agarose gel rig and a transilluminator",
        ),
        files=(RECORD_FILE,),
        steps=_steps(record, control, length),
        references=REFERENCES,
        holes=(LTR_HOLE, PROMOTER_HOLE),
    )


def _highlights(
    report: DomesticationReport, still: dict[str, int], proteins: tuple[str, ...]
) -> tuple[str, ...]:
    """Return what a reader has to read rather than scan: the route, and what it leaves."""
    left = ", ".join(f"{name} {count}" for name, count in still.items())
    return (
        "Two sites sit at the same offset in the two long terminal repeats, and 610 bases "
        "around each one are identical between the copies. No primer, assembly overlap or "
        "Type IIS overhang is unique to one copy, so mutagenesis and assembly both fail at "
        "those two.",
        "Leaving those two standing would make this a cheaper job. They are changed here "
        "because a published vector carries the change, and changing them is what puts the "
        "whole plasmid on the order form rather than a fragment of it.",
        f"{len(report.changes)} of the sites sit in a coding sequence, where a synonymous codon "
        f"takes the site away: {' and '.join(proteins)} read exactly as the parent reads them.",
        f"The other {len(EDITS)} cost a stated base each. A published vector already carries the "
        "change made in both repeats, and another carries the change made at the BsaI site of "
        "the hPGK promoter. Neither group measured what the change cost, so neither number is "
        "stated here.",
        "The BsmBI site in the hPGK promoter has no precedent anywhere. The base chosen neither "
        "makes nor breaks a CpG, which is the most this run can say for it.",
        f"Only BsaI and BsmBI were taken out. The result still reads {left}, so an assembly "
        "whose own enzyme is one of those meets those sites unchanged.",
    )


def _checks(
    record: SequenceRecord,
    report: DomesticationReport,
    still: dict[str, int],
    proteins: tuple[str, ...],
) -> tuple[Check, ...]:
    """Return the three verdicts on the record, each counted off the file this writes."""
    left = ", ".join(f"{name} reads {count} sites" for name, count in still.items())
    return (
        Check(
            "sites",
            "pass",
            f"no BsaI and no BsmBI site is left in {len(record):,} bases",
        ),
        Check(
            "protein",
            "pass",
            f"{' and '.join(proteins)} translate as the parent translates them, across "
            f"{len(report.changes)} changed codons",
        ),
        Check("other enzymes", "warn", f"{left}; neither was in this job"),
    )


def _steps(record: SequenceRecord, control: Enzyme, length: str) -> tuple[Step, ...]:
    """Return the four steps: order it, grow it, show nothing cuts it, read both repeats."""
    return (
        Step(
            "Order the domesticated plasmid",
            key="order",
            instructions=(
                f"Send {RECORD_FILE} to a supplier that builds whole plasmids, as one clonal "
                f"construct of {length} bp.",
                "Declare the two 634 bp direct repeats on the order form.",
                "Ask for the plasmid verified over its whole sequence.",
            ),
            cautions=(
                "A supplier whose clonal genes stop at 7 kb cannot build this plasmid; ask for "
                "a long-construct quote before ordering.",
            ),
            notes=(
                "A direct repeat longer than 200 bases is a high-complexity sequence, which is "
                "priced and scheduled apart from the rest, so declaring it up front saves a "
                "rejected order (Twist Bioscience, gene synthesis sequence acceptance "
                "criteria).",
                "The change at transcript position +4 goes into both repeats because a vector "
                "whose two repeats differ there changes what reverse transcription copies "
                "(Haellman et al. 2021, which made the same change in both).",
            ),
            expected=("The supplier accepts the sequence and quotes a build.",),
        ),
        Step(
            "Transform and pick a colony",
            key="transform",
            instructions=(
                "Transform 1 µL of the plasmid into a recombination-deficient cloning strain.",
                "Plate on LB with 100 µg/mL carbenicillin and grow overnight at 30 °C.",
                "Pick one colony into 5 mL of LB with carbenicillin, grow overnight at 30 °C, "
                "and miniprep it.",
            ),
            cautions=(
                "Use a recombination-deficient strain and a 5 mL culture; a larger one is not "
                "worth the deletions it returns.",
            ),
            expected=(
                "Colonies by the next morning at the latest, and a miniprep of a few micrograms.",
            ),
        ),
        Step(
            "Show that neither enzyme cuts",
            key="confirm-digest",
            instructions=(
                "Digest 1 µg of the miniprep with 10 units of BsaI-HFv2 in rCutSmart Buffer for "
                "1 hour at 37 °C.",
                "Digest a second 1 µg with 10 units of BsmBI-v2 in NEBuffer r3.1 for 1 hour at "
                "55 °C.",
                f"Digest a third 1 µg with 10 units of {control.commercial_name or control.name}.",
                "Run all three beside 1 µg of undigested miniprep on a 0.8% agarose gel.",
            ),
            expected=(
                "The BsaI and BsmBI lanes look like the undigested lane.",
                f"The control lane is one band at about {length} bp.",
            ),
            notes=(
                "An enzyme that did not cut and an enzyme nobody added look the same on a gel, "
                "which is what the control lane tells apart (New England Biolabs, single "
                "enzyme digest protocol).",
            ),
        ),
        Step(
            "Read both repeats",
            key="confirm-sequence",
            instructions=(
                "Send 500 ng of the miniprep for whole-plasmid sequencing, and ask for the "
                "whole plasmid assembled rather than reads from primers.",
                f"Compare the assembled consensus with {RECORD_FILE} base by base.",
                "Confirm that both long terminal repeats carry their change, one at a time.",
            ),
            expected=(
                f"The consensus matches {RECORD_FILE} at every base, both repeats included.",
            ),
        ),
    )


def _check_written(record: SequenceRecord, path: Path) -> None:
    """Refuse to leave behind a file that reads back as a different record.

    Raises
    ------
    ValueError
        Naming the first thing that came back differently.
    """
    again = read_record(path)
    if (again.sequence, again.topology) != (record.sequence, record.topology):
        raise ValueError(f"{path} does not read back as the record the build made")
    if _marks(again) != _marks(record):
        raise ValueError(f"{path} does not read back every base the build changed")


def _marks(record: SequenceRecord) -> list[tuple[str, int]]:
    """Every base the build marked as changed, which a written file has to carry back."""
    return [
        (one.name, int(one.segments[0].start))
        for one in record.features
        if one.name.endswith("site removed")
    ]


def write(record: SequenceRecord, made: Protocol, into: Path) -> list[Path]:
    """Write the record, the protocol and the page it renders to, and say where each went."""
    into.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location(
        "build_dmx_vector", Path(__file__).parent / "build_dmx_vector.py"
    )
    assert spec is not None
    assert spec.loader is not None
    sibling = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = sibling
    spec.loader.exec_module(sibling)
    sibling.write_genbank(record, into / RECORD_FILE)
    _check_written(record, into / RECORD_FILE)
    data = write_protocol(minted(made), into / PROTOCOL_DATA_FILE)
    page = write_html(read_protocol(data), into / PROTOCOL_FILE)
    return [into / RECORD_FILE, data, page]


def main(argv: list[str] | None = None) -> int:
    """Build the record and its protocol, printing what the result measures."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from", dest="source", type=Path, default=PARENT, help="the parent")
    parser.add_argument("--out", type=Path, default=OUT, help="the directory they go in")
    args = parser.parse_args(argv)
    if not args.source.exists():
        print(
            f"{args.source} is not here. The parent ships in tests/data/; "
            "docs/research/working-vector-plvx-tetone.md says where it came from.",
            file=sys.stderr,
        )
        return 1
    parent = read_record(args.source)
    record, report = rebuild(parent)
    made = protocol(parent, record, report)
    for path in write(record, made, args.out):
        print(f"{path}")
    print(f"  {record.name}, {len(record)} bp {record.topology}")
    for one in report.changes:
        print(
            f"  {one.site.enzyme.name} at {one.site.start}: codon {one.codon_index} "
            f"{one.old_codon} -> {one.new_codon} ({one.amino_acid})"
        )
    for edit in EDITS:
        print(f"  {edit.enzyme} at {edit.position}: {edit.before} -> {edit.after}, {edit.element}")
    for name in ("BsaI", "BsmBI", "BbsI", "SapI", "PaqCI"):
        print(f"  {name:6} {len(find_sites(record, get_enzyme(name)))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

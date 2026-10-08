"""The method's gate, on a design that is known good and on seven it is known to break.

The good half is `docs/examples/ap1-library`, read as committed files rather than planned again:
the gate judges finished records, so what it reads here is what an agent would hand it. Each bad
half breaks one rule of that same design, and every break is one that was actually measured --
the BsmBI sites that reached designed coding regions before #220, a fragment out of frame, a
barcode spelling a stop where the construct reads it, two barcodes inside the distance rule, and
a block one base short of the length convention it was built to.

The rounds run in the DMX vector, so the gate judges the destination as one where a primer that
reads a well binds it. The AP-1 destination is a minimal stand-in no such primer reads, so the
DMX vector is also built here from the method's own stuffers, which is where its cassette comes
from.
"""

import csv
from collections import defaultdict
from pathlib import Path

import pytest

from mbio.enzymes import Enzyme, get_enzyme
from mbio.io import read_record
from mbio.ligase import LigaseProfile
from mbio.sequence import Segment, SequenceRecord, reverse_complement
from mbio.sites import CutSite, released
from synbio.igga.gate import (
    WELL_PRIMERS,
    Verdict,
    check_dmx_vector,
    check_library,
    check_reaction,
    final_assembly_reactions,
    library_reactions,
)
from synbio.igga.method import IGGA
from synbio.igga.project import read_build
from synbio.igga.vector import Cassette

DEMO = Path(__file__).parents[3] / "docs" / "examples" / "ap1-library"

#: Where a break is cut into a block: inside the coding bases, past the 29-base 5' external
#: stuffer and well short of the internal stuffer every block ends with.
INSIDE_CODING = 100

#: BsmBI's recognition site, which is the one no block may spell anywhere.
BSMBI = "CGTCTC"

#: Bases spelling no site of the method's enzymes, laid round the DMX vector's cassette.
FILLER = "ACGATCGTTA" * 20


@pytest.fixture(scope="module")
def build():
    return read_build(DEMO / "project.json")


@pytest.fixture(scope="module")
def blocks():
    found = defaultdict(list)
    with (DEMO / "parts.tsv").open(encoding="utf-8") as sheet:
        for row in csv.DictReader(sheet, delimiter="\t"):
            found[row["position"]].append(SequenceRecord(row["sequence"], name=row["name"]))
    return dict(found)


@pytest.fixture(scope="module")
def barcodes():
    found = defaultdict(list)
    with (DEMO / "barcodes.tsv").open(encoding="utf-8") as table:
        for row in csv.DictReader(table, delimiter="\t"):
            found[row["position"]].append(row["barcode"])
    return dict(found)


@pytest.fixture(scope="module")
def destination():
    return read_record(DEMO / "vector.gb")


@pytest.fixture(scope="module")
def products():
    return [read_record(DEMO / name) for name in ("round-1.dna", "round-2.dna", "product.dna")]


@pytest.fixture(scope="module")
def judge(build, destination, blocks, barcodes, products):
    """Judge the AP-1 design, with whatever part of it a test swapped out."""

    def run(**changed) -> Verdict:
        given = {
            "destination": destination,
            "blocks": blocks,
            "barcodes": barcodes,
            "products": products,
            **changed,
        }
        return check_library(build, **given)

    return run


@pytest.fixture(scope="module")
def good(judge) -> Verdict:
    return judge()


def _swap(blocks, position, index, sequence):
    """Return the blocks with one of them replaced."""
    changed = {name: list(found) for name, found in blocks.items()}
    changed[position][index] = SequenceRecord(sequence, name=f"broken {position}")
    return changed


def _codes(barcodes, position, index, code):
    """Return the barcode table with one barcode replaced."""
    changed = {name: list(found) for name, found in barcodes.items()}
    changed[position][index] = code
    return changed


def _dmx_vector(broken=""):
    """A DMX-shaped vector, with `broken` laid inside the forward primer's footprint.

    The cassette is the method's own external stuffers, which carry a blunt site outboard of
    each releasing site, and a primer that reads a well binds either side of it. The cargo
    carries an internal stuffer, so the blunt site a part brings is there to be judged too.
    """
    forward, reverse = WELL_PRIMERS
    cargo = FILLER[:100] + IGGA.internal_stuffer + FILLER[:100]
    return SequenceRecord(
        FILLER
        + forward[:10]
        + broken
        + forward[10:]
        + FILLER[:20]
        + IGGA.external_stuffer_5
        + cargo
        + IGGA.external_stuffer_3
        + FILLER[:20]
        + reverse_complement(reverse)
        + FILLER,
        topology="circular",
        name="DMX-test",
    )


# The known-good half.


def test_the_gate_passes_the_ap1_design_it_was_given(good):
    assert good.status == "pass"
    assert not good.failures


def test_the_gate_judges_every_molecule_of_the_design(good, build):
    named = {one.name for one in good.judgements}
    assert named == {
        "destination opens",
        "donor releases",
        "ends distinguishable",
        "ligation fidelity",
        "cargo sites",
        "cargo frame",
        "well primers",
        "blunt sites",
        "barcode spacing",
        "barcode reading",
        "barcode block",
        "product opens",
        "terminal block",
        "terminal stop",
    }
    cargo = [one for one in good.judgements if one.name == "cargo sites"]
    assert len(cargo) == 72


def test_a_fidelity_nothing_judges_carries_no_verdict(good):
    assert good["ligation fidelity"].status is None
    assert good["ligation fidelity"].check.value > 0.9


# A ligase matrix the user holds. The rates are this file's own, chosen either side of
# `STRONG_LIGATION`; the sheet name is the one Bilotti 2022's workbook carries, and none of
# that workbook is here.
AP1_OVERHANGS = ("AGGA", "AGAT", "GCAT", "TTCC")
SHEET = "File S7. T7 PEG"


def _t7_peg(**rates: float):
    """A sheet of a workbook the user holds, joining each AP-1 overhang at the rate asked for.

    A rate is per 100,000 events, which the filler row makes exact.
    """
    joined = {one: rates.get(one, 300.0) for one in AP1_OVERHANGS}
    counts = {one: {reverse_complement(one): round(rate * 10)} for one, rate in joined.items()}
    counts["AAAA"] = {"TTTT": 1_000_000 - sum(row[one] for row in counts.values() for one in row)}
    return LigaseProfile(
        Path("File S1_NAR.xlsx"),
        conditions=SHEET,
        overhang_length=4,
        observations=1_000_000,
        counts=counts,
    )


def test_without_a_ligase_matrix_the_ligation_is_judged_as_it_always_was(good):
    assert "ligation on-target rate" not in {one.name for one in good.judgements}


def test_a_ligase_matrix_adds_a_check_and_moves_the_fidelity_score_not_at_all(good, judge):
    judged = judge(profile=_t7_peg())

    assert judged["ligation fidelity"].check == good["ligation fidelity"].check
    assert judged["ligation on-target rate"].status == "pass"
    assert judged.status == "pass"


def test_an_overhang_the_ligase_joins_rarely_warns_and_names_the_sheet(judge):
    judged = judge(profile=_t7_peg(AGGA=60.8))

    one = judged["ligation on-target rate"]
    assert one.status == "warn"
    assert "AGGA at 60.8" in one.check.detail
    assert SHEET in one.check.detail
    assert judged.status == "warn"
    assert not judged.failures


def test_the_chain_is_two_digests_and_a_ligation_a_round(build, destination, blocks, products):
    made = library_reactions(build, destination=destination, blocks=blocks, products=products)
    assert [one.kind for one in made] == ["digest", "digest", "ligation"] * 3
    assert made[0]["destination"].one is destination
    assert made[3]["destination"].one is products[0]
    assert len(made[1]["donor"]) == 24


def test_the_gate_hands_back_the_findings_behind_its_checks(good):
    assert good.findings == ()


# The known-bad half: one rule at a time.


def test_a_bsmbi_site_in_a_designed_coding_region_is_caught(judge, blocks):
    one = str(blocks["N"][0].sequence)
    broken = one[:INSIDE_CODING] + BSMBI + one[INSIDE_CODING + len(BSMBI) :]
    verdict = judge(blocks=_swap(blocks, "N", 0, broken))
    assert verdict.status == "fail"
    failed = verdict["cargo sites"]
    assert failed.status == "fail"
    assert "carries 1 site it must not, of BsmBI" in failed.check.detail
    assert {one.enzyme.name for one in failed.findings} == {"BsmBI"}


def test_a_fragment_out_of_frame_is_caught(judge, blocks):
    one = str(blocks["N"][0].sequence)
    verdict = judge(blocks=_swap(blocks, "N", 0, one[:INSIDE_CODING] + one[INSIDE_CODING + 1 :]))
    assert verdict.status == "fail"
    failed = verdict["cargo frame"]
    assert failed.status == "fail"
    assert failed.check.value == 2
    assert "a cargo another part follows" in failed.check.detail


def test_a_capping_block_one_base_short_is_caught(judge, blocks):
    one = str(blocks["C"][0].sequence)
    verdict = judge(blocks=_swap(blocks, "C", 0, one[:INSIDE_CODING] + one[INSIDE_CODING + 1 :]))
    assert verdict.status == "fail"
    failed = next(
        judged
        for judged in verdict.failures
        if judged.name == "cargo frame" and judged.where.startswith("the C")
    )
    assert failed.check.value == 0
    assert "the block that ends the chain" in failed.check.detail


def test_a_stop_codon_at_the_barcodes_frame_offset_is_caught(judge, barcodes):
    verdict = judge(barcodes=_codes(barcodes, "N", 0, "TATAATTGTCT"))
    assert verdict.status == "fail"
    failed = verdict["barcode reading"]
    assert failed.status == "fail"
    assert "a stop, in the frame the construct reads" in failed.check.detail


def test_two_barcodes_inside_the_distance_rule_are_caught(judge, barcodes):
    near = barcodes["N"][0][:-1] + ("A" if barcodes["N"][0][-1] != "A" else "C")
    verdict = judge(barcodes=_codes(barcodes, "N", 1, near))
    assert verdict.status == "fail"
    failed = verdict["barcode spacing"]
    assert failed.status == "fail"
    assert failed.check.value == 1
    assert "under the 3 two parts of one list need" in failed.check.detail


def test_a_product_that_does_not_carry_the_barcodes_the_table_names_is_caught(
    judge, build, products
):
    bases = str(products[-1].sequence)
    at = bases.find(build.scheme.internal_stuffer_core) + len(build.scheme.internal_stuffer_core)
    broken = bases[:at] + "AAATTTGGGCC" + bases[at + build.barcode.length :]
    verdict = judge(
        products=[*products[:-1], SequenceRecord(broken, topology="circular", name="broken")]
    )
    assert verdict.status == "fail"
    failed = verdict["barcode block"]
    assert failed.status == "fail"
    assert "names no part of position C" in failed.check.detail


def test_a_destination_the_round_cannot_open_cleanly_is_caught(judge, destination):
    bases = str(destination.sequence)
    broken = bases[:100] + "GAAGAC" + bases[106:]
    verdict = judge(destination=SequenceRecord(broken, topology="circular", name="broken"))
    assert verdict.status == "fail"
    failed = verdict["destination opens"]
    assert failed.status == "fail"
    assert (
        "BbsI leaves 1 of 1 destination of round 1 opening uncut in the 2 places"
        in failed.check.detail
    )


# The DMX vector, which the validating experiment amplifies a well of.


def test_a_blunt_site_between_a_releasing_site_and_the_primer_is_where_it_belongs(build):
    """The two the cassette carries are judged; the one the cargo's stuffer carries is not."""
    vector = _dmx_vector()
    bases = str(vector.sequence)
    forward = WELL_PRIMERS[0]
    clear = f"{bases.index(forward) + len(forward) + 1} .. {bases.index('GGTCTC')}"
    (judged,) = check_dmx_vector(vector, build=build)
    assert judged.status == "pass"
    assert "GCCCGGGC" in bases
    assert judged.check.value == bases.count("GTTTAAAC")
    assert clear in judged.check.detail


def test_a_blunt_site_inside_a_primers_footprint_is_caught(build):
    vector = _dmx_vector(broken="GTTTAAAC")
    bases = str(vector.sequence)
    reads = WELL_PRIMERS[0][10:]
    clear = bases.index(reads) + len(reads) + 1
    (judged,) = check_dmx_vector(vector, build=build)
    assert judged.status == "fail"
    assert f"PmeI at {bases.index('GTTTAAAC') + 1} .." in judged.check.detail
    assert f"outside {clear} .." in judged.check.detail
    (finding,) = judged.findings
    assert isinstance(finding, CutSite)
    assert finding.enzyme.name == "PmeI"


def test_a_vector_a_primer_no_longer_reads_says_so_rather_than_passing(build):
    bases = str(_dmx_vector().sequence).replace(WELL_PRIMERS[0], FILLER[: len(WELL_PRIMERS[0])])
    (judged,) = check_dmx_vector(
        SequenceRecord(bases, topology="circular", name="unread"), build=build
    )
    assert judged.status == "fail"
    assert "binds it in 0 places" in judged.check.detail


def test_a_cargo_carrying_an_annealing_region_is_caught(judge, blocks):
    """A second priming site in a well leaves that well's read uncallable."""
    one = str(blocks["N"][0].sequence)
    region = WELL_PRIMERS[0]
    broken = one[:INSIDE_CODING] + region + one[INSIDE_CODING + len(region) :]
    verdict = judge(blocks=_swap(blocks, "N", 0, broken))
    assert verdict.status == "fail"
    failed = verdict["well primers"]
    assert failed.status == "fail"
    assert "a second priming site" in failed.check.detail
    (finding,) = failed.findings
    assert isinstance(finding, Segment)


def test_a_destination_no_well_primer_reads_is_not_judged_as_one(judge):
    """A destination carrying a cassette and nothing else: no primer reads it, so nothing judges it."""
    stand_in = SequenceRecord(
        FILLER + IGGA.internal_stuffer + FILLER, topology="circular", name="stand-in"
    )

    with pytest.raises(KeyError):
        judge(destination=stand_in)["blunt sites"]


def test_a_dmx_destination_is_judged_where_a_build_accepts_one(judge):
    assert judge(destination=_dmx_vector())["blunt sites"].status == "pass"


def test_a_dmx_destination_a_blunt_site_breaks_fails_the_build(judge):
    verdict = judge(destination=_dmx_vector(broken="GTTTAAAC"))
    assert verdict.status == "fail"
    assert verdict["blunt sites"].status == "fail"


#: The enzyme that admits cargo to a working vector here, as #256 chose for pLVX.
CARGO = "PaqCI"


def _carrying(cassette: Cassette) -> SequenceRecord:
    """A circular plasmid of filler giving `cassette` up to its own enzyme."""
    return SequenceRecord(cassette.bases + FILLER, topology="circular", name="stand-in")


def _cassette(enzyme: Enzyme, payload: str, *, left: str = "", right: str = "") -> Cassette:
    """A piece `enzyme` frees, cutting out onto `left` at one end and `right` at the other."""
    reach = ("TA" * enzyme.top_cut)[: enzyme.top_cut - len(enzyme.site)]
    bases = (
        (left or IGGA.entry_overhang)
        + reach
        + reverse_complement(enzyme.site)
        + payload
        + enzyme.site
        + reach
        + (right or IGGA.scar_overhang)
    )
    return Cassette(bases, enzyme, (), "the digest this stand-in is cut in")


@pytest.fixture(scope="module")
def final():
    """The two molecules of a final assembly: a library to free cargo from, and a working vector."""
    library = _carrying(_cassette(IGGA.external, "ACGTTGCA" * 12))
    working = _carrying(_cassette(get_enzyme(CARGO), "ACGTTGCA" * 20))
    return library, working


def test_the_final_assembly_is_two_digests_and_a_ligation(build, final):
    library, working = final
    made = final_assembly_reactions(
        build, library=library, working=working, cargo=get_enzyme(CARGO)
    )
    judged = [one for reaction in made for one in check_reaction(reaction, build=build)]

    assert [one.kind for one in made] == ["digest", "digest", "ligation"]
    assert made[-1]["destination"].cutter == get_enzyme(CARGO)
    assert made[-1]["donor"].cutter == IGGA.external
    assert [one.status for one in judged] == ["pass", "pass", "pass", None]
    assert "AGGA, TTCC" in judged[2].check.detail


def test_a_working_vector_opening_on_the_wrong_ends_fails_the_ligation(build, final):
    """The digest passes it: two cuts is two cuts. Only the ends it leaves say it is wrong."""
    library, _ = final
    askew = _carrying(_cassette(get_enzyme(CARGO), "ACGTTGCA" * 20, left="AAAA", right="TTTT"))
    made = final_assembly_reactions(build, library=library, working=askew, cargo=get_enzyme(CARGO))
    judged = [one for reaction in made for one in check_reaction(reaction, build=build)]

    assert [one.status for one in judged[:2]] == ["pass", "pass"]
    assert judged[2].status == "fail"
    assert "AAAA, AGGA, TTCC, TTTT" in judged[2].check.detail


def test_a_build_names_both_the_working_vector_and_its_enzyme_or_neither(judge, final):
    _, working = final
    with pytest.raises(ValueError, match="name both, or neither"):
        judge(working=working)


def test_a_donor_held_in_a_circular_backbone_is_judged_on_its_cargo(judge, blocks):
    """Two pieces come off a circular donor, and the cloning scar says which one is the part."""
    one = blocks["N"][0]
    entry, scar = IGGA.entry_overhang, IGGA.scar_overhang
    cargo = one.sequence[len(IGGA.external_stuffer_5) - len(entry) :]
    cargo = cargo[: len(cargo) - len(IGGA.external_stuffer_3) + len(scar)]
    circular = SequenceRecord(
        IGGA.external_stuffer_5[: -len(entry)]
        + cargo
        + IGGA.external_stuffer_3[len(scar) :]
        + FILLER,
        topology="circular",
        name="N in a DMX backbone",
    )

    assert len(released(circular, IGGA.external)) == 2
    assert judge(blocks=_swap_record(blocks, "N", 0, circular)).status == "pass"


def _swap_record(blocks, position, index, record):
    """Return the blocks with one of them replaced by a record as it stands."""
    changed = {name: list(found) for name, found in blocks.items()}
    changed[position][index] = record
    return changed

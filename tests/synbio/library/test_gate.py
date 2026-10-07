"""The method's gate, on a design that is known good and on seven it is known to break.

The good half is `docs/examples/ap1-library`, read as committed files rather than planned again:
the gate judges finished records, so what it reads here is what an agent would hand it. Each bad
half breaks one rule of that same design, and every break is one that was actually measured --
the BsmBI sites that reached designed coding regions before #220, a fragment out of frame, a
barcode spelling a stop where the construct reads it, two barcodes inside the distance rule, and
a block one base short of the length convention it was built to.
"""

import csv
from collections import defaultdict
from pathlib import Path

import pytest

from liulab_mbio.io import read_record
from liulab_mbio.sequence import SequenceRecord
from liulab_synbio.library.gate import Verdict, check_library, library_reactions
from liulab_synbio.library.project import read_project

DEMO = Path(__file__).parents[3] / "docs" / "examples" / "ap1-library"

#: Where a break is cut into a block: inside the coding bases, past the 29-base 5' external
#: stuffer and well short of the internal stuffer every block ends with.
INSIDE_CODING = 100

#: BsmBI's recognition site, which is the one no block may spell anywhere.
BSMBI = "CGTCTC"


@pytest.fixture(scope="module")
def project():
    return read_project(DEMO / "project.json")


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
def judge(project, destination, blocks, barcodes, products):
    """Judge the AP-1 design, with whatever part of it a test swapped out."""

    def run(**changed) -> Verdict:
        given = {
            "destination": destination,
            "blocks": blocks,
            "barcodes": barcodes,
            "products": products,
            **changed,
        }
        return check_library(project, **given)

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


# The known-good half.


def test_the_gate_passes_the_ap1_design_it_was_given(good):
    assert good.status == "pass"
    assert not good.failures


def test_the_gate_judges_every_molecule_of_the_design(good, project):
    named = {one.name for one in good.judgements}
    assert named == {
        "destination opens",
        "donor releases",
        "ends distinguishable",
        "ligation fidelity",
        "cargo sites",
        "cargo frame",
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


def test_the_chain_is_two_digests_and_a_ligation_a_round(project, destination, blocks, products):
    made = library_reactions(project, destination=destination, blocks=blocks, products=products)
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
    assert "carries 1 site(s) it must not, of BsmBI" in failed.check.detail
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
    judge, project, products
):
    bases = str(products[-1].sequence)
    at = bases.find(project.scheme.internal_stuffer_core) + len(
        project.scheme.internal_stuffer_core
    )
    broken = bases[:at] + "AAATTTGGGCC" + bases[at + project.barcode.length :]
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
    assert "are not cut by BbsI in the 2 places this method cuts them" in failed.check.detail

"""The verdicts a Gateway plan carries beyond one reaction's: the insert, the host, the markers
and the frame."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from liulab_mbio.bench.phenotype import read_phenotype
from liulab_mbio.cloning.gateway import plan_gateway
from liulab_mbio.cloning.gateway.bench import DEFAULT_HOST, PROPAGATION_HOST
from liulab_mbio.cloning.gateway.checks import plan_checks
from liulab_mbio.cloning.gateway.design import C_TERMINAL_FRAME, N_TERMINAL_FRAME
from liulab_mbio.cloning.gateway.recombination import PlannedReaction, recombine
from liulab_mbio.cloning.plan import status
from liulab_mbio.sequence import reverse_complement

from .records import att_site, destination_vector, donor_vector, entry_clone

if TYPE_CHECKING:
    from collections.abc import Sequence

    from liulab_mbio.checks import Check
    from liulab_mbio.cloning.gateway import Plan
    from liulab_mbio.cloning.gateway.design import Fusion
    from liulab_mbio.sequence import SequenceRecord

#: What a fusion primer adds at each end, read on the top strand: the forward tail's bases as
#: `liulab_mbio.cloning.gateway.design` writes them, and the reverse tail's turned over, because
#: that tail is written on the other strand. The frame judged here is the one those primers make.
FRAME_BASES = (N_TERMINAL_FRAME, reverse_complement(C_TERMINAL_FRAME))


def _check(checks: Sequence[Check], name: str) -> Check:
    """Return the check of that name, from a plan's checks or from the plan-level ones."""
    return next(check for check in checks if check.name == name)


def _judged(
    carrier: SequenceRecord,
    acceptor: SequenceRecord,
    *,
    host: str = DEFAULT_HOST,
    fusion: Fusion = "none",
) -> tuple[Check, ...]:
    """Return the plan-level checks for one LR reaction over these two records.

    Every check here reads the reaction: its product, its junctions, what moved and what the
    product says about itself. None reads an oligo, so none is designed. `status` over the
    result is the worst of them, as a plan's own is over all of its.
    """
    made = recombine(carrier, acceptor, reaction="LR")
    first, second = made.junctions
    lr = PlannedReaction(
        made,
        (),
        read_phenotype(
            made.product,
            (first.start, second.end),
            vector=made.backbone.record,
            span=made.cassette,
        ),
    )
    return plan_checks(lr=lr, bp=None, host=host, fusion=fusion)


def _fusion_checks(
    gfp: SequenceRecord, added: tuple[str, str], fusion: Fusion
) -> tuple[Check, ...]:
    """Judge a fusion from an insert whose stop codon is gone, into a tagged destination vector."""
    return _judged(
        entry_clone(gfp.sequence[:-3], added=added),
        destination_vector(tag="6xHis"),
        fusion=fusion,
    )


@pytest.mark.parametrize("orientation", ["forward", "reverse"])
def test_an_att_site_inside_the_insert_is_found_either_way_round(
    gfp: SequenceRecord, destination: SequenceRecord, orientation: str
) -> None:
    site = att_site("attB1")
    carried = site if orientation == "forward" else reverse_complement(site)
    checks = _judged(entry_clone(gfp.sequence + carried), destination)

    check = _check(checks, "insert att sites")
    assert (check.status, check.value) == ("fail", 1)
    assert "attB1" in check.detail
    assert status(checks) == "fail"


def test_the_ccdb_rule_judges_the_f_prime_strain_and_names_what_grows_a_vector(
    gateway_plan: Plan, entry: SequenceRecord, destination: SequenceRecord
) -> None:
    refused = plan_gateway(entry, destination, host="One Shot TOP10F'")

    assert _check(gateway_plan.checks, "ccdB host").status == "pass"
    assert _check(refused.checks, "ccdB host").status == "fail"
    assert "ccdA" in _check(refused.checks, "ccdB host").detail
    assert refused.status == "fail"


def test_the_ccdb_resistant_host_rule_is_guidance_and_carries_no_verdict(
    entry: SequenceRecord, destination: SequenceRecord
) -> None:
    plan = plan_gateway(entry, destination, host="One Shot ccdB Survival 2 T1R")
    vector_host = _check(plan.checks, "ccdB vector host")
    selecting = _check(plan.checks, "ccdB host")

    assert vector_host.status is None
    assert PROPAGATION_HOST in vector_host.detail
    assert selecting.status is None
    assert "no manual forbids" in selecting.detail
    assert plan.status == "pass"


def test_a_selection_that_cannot_tell_the_two_clones_apart_is_reported(
    gateway_plan: Plan, gfp: SequenceRecord, destination: SequenceRecord
) -> None:
    same = plan_gateway(entry_clone(gfp.sequence, marker="AmpR"), destination)
    apart_marker = _check(gateway_plan.checks, "LR markers")
    same_marker = _check(same.checks, "LR markers")

    assert (apart_marker.status, apart_marker.value) == ("pass", 2)
    assert apart_marker.detail.startswith("ampicillin or carbenicillin")
    assert "kanamycin on the entry clone" in apart_marker.detail
    assert same_marker.status == "warn"
    assert "unreacted entry clone grows on it too" in same_marker.detail
    assert same.status == "warn"


@pytest.mark.parametrize(
    ("marker", "plate"),
    [
        ("KanR", "LB agar plates with kanamycin"),
        ("ZeoR", "Low Salt LB agar plates with Zeocin"),
        ("PuroR", "LB agar plates with the antibiotic PuroR selects"),
    ],
)
def test_the_bp_plate_names_the_donor_marker_and_translates_only_what_is_sourced(
    insert: SequenceRecord, destination: SequenceRecord, marker: str, plate: str
) -> None:
    plan = plan_gateway(insert, destination, donor=donor_vector(marker=marker))

    assert plate in {material.name for material in plan.protocol().materials}


def test_the_frame_is_judged_at_each_end_a_fusion_reads_through(gfp: SequenceRecord) -> None:
    checks = _fusion_checks(gfp, FRAME_BASES, "both")
    n_end = _check(checks, "N-terminal frame")
    c_end = _check(checks, "C-terminal frame")

    assert (n_end.status, n_end.value) == ("pass", 2)
    assert "TSLYKKAGS" in n_end.detail
    assert "in frame with 6xHis upstream" in n_end.detail
    assert (c_end.status, c_end.value) == ("pass", 1)
    assert "YPAFLYKVV" in c_end.detail
    assert status(checks) == "pass"


def test_no_frame_is_judged_where_no_fusion_reads_through(gateway_plan: Plan) -> None:
    assert {"N-terminal frame", "C-terminal frame"}.isdisjoint(
        check.name for check in gateway_plan.checks
    )


def test_two_added_bases_that_spell_a_stop_break_the_n_terminal_frame(
    gfp: SequenceRecord,
) -> None:
    check = _check(_fusion_checks(gfp, ("AA", FRAME_BASES[1]), "N-terminal"), "N-terminal frame")

    assert check.status == "fail"
    assert check.detail.endswith("spelling TSLYKKAG*, a stop")


def test_a_gene_keeping_its_stop_codon_breaks_the_c_terminal_frame(
    gfp: SequenceRecord,
) -> None:
    checks = _judged(
        entry_clone(gfp.sequence, added=FRAME_BASES),
        destination_vector(tag="6xHis"),
        fusion="C-terminal",
    )
    check = _check(checks, "C-terminal frame")

    assert check.status == "fail"
    assert check.detail == "GFP ends in a stop codon, which a C-terminal fusion has to lose"


def test_every_check_reaches_the_page_and_the_unjudged_one_shows_no_verdict(
    gateway_plan: Plan,
) -> None:
    badges = gateway_plan.protocol().checks

    assert [badge.name for badge in badges] == [check.name for check in gateway_plan.checks]
    assert [badge.name for badge in badges if badge.status is None] == ["ccdB vector host"]

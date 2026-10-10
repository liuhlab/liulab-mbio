"""Reading a Sanger trace into the calls it trusts."""

from pathlib import Path

import pytest

from mbio.sequence import Feature, Segment, SequenceRecord
from mbio.verification.judge import verify
from mbio.verification.result import SequencingResult
from mbio.verification.trace import read_trace


@pytest.fixture(scope="module")
def clean(data_dir: Path) -> SequencingResult:
    """A 3730xl read with qualities: mCherry, a linker, then EGFP."""
    return read_trace(data_dir / "3730.ab1")


def test_a_clean_trace_trusts_its_mott_span_and_calls_no_mixed_base_inside_it(clean):
    """A clean trace's trusted span is abi-trim's, and holds only definite calls."""
    assert clean.name == "3730"
    assert clean.quality is not None
    assert len(clean.quality) == len(clean.bases) == 1165
    assert clean.trusted == (14, 1089)
    assert set(clean.bases[14:1089]) <= set("ACGT")


def test_a_second_peak_makes_a_call_mixed_and_the_basecaller_s_own_codes_are_kept(clean):
    """The basecaller wrote G at 6 and K at 8; a T peak at 6 makes it K too."""
    assert clean.bases[6] == "K"
    assert clean.bases[8] == "K"


def test_a_trace_with_no_quality_trusts_nothing(data_dir: Path):
    """Every quality in a 310 trace is zero, so no span of it is trusted."""
    result = read_trace(data_dir / "310.ab1")
    assert result.trusted == (0, 0)
    assert result.quality is not None
    assert not any(result.quality)


def test_a_trace_verifies_against_its_own_calls_until_a_substitution_is_planted(clean):
    """A record built from the read's trusted calls reads true, and one changed base fails it."""
    start, end = clean.trusted
    calls = clean.bases[start:end]
    insert = (Feature("insert", "misc_feature", (Segment(100, 900),)),)
    assert verify(SequenceRecord(calls), (clean,), insert).verified
    planted = calls[:500] + next(one for one in "ACGT" if one != calls[500]) + calls[501:]
    found = verify(SequenceRecord(planted), (clean,), insert)
    assert [(one.kind, one.start, one.end) for one in found.disagreements] == [
        ("substitution", 500, 501)
    ]
    assert found.checks[0].status == "fail"


def test_a_trace_with_no_quality_says_so_through_verify(data_dir: Path):
    """A 310 trace trusts nothing, so it carries no verdict, and the detail says why."""
    result = read_trace(data_dir / "310.ab1")
    found = verify(SequenceRecord("GATTACA"), (result,), ())
    assert found.result_checks[0].status is None
    assert found.result_checks[0].detail == "the trace carries no quality values"

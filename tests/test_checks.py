from liulab_mbio.checks import Check, worst, worst_of


def test_the_worst_verdict_is_the_result() -> None:
    assert worst(("pass", "fail", "warn")) == "fail"
    assert worst(("warn", "pass")) == "warn"
    assert worst(("pass", "pass")) == "pass"


def test_a_check_with_no_verdict_never_makes_the_result_worse() -> None:
    checks = (Check("tm", "warn", 62.3), Check("end_stability", None, -4.0))
    assert worst(check.status for check in checks) == "warn"
    assert worst_of(checks) == "warn"
    assert worst(("pass", None)) == "pass"
    # Nothing judged at all is not a failure.
    assert worst((None, None)) == "pass"
    assert worst(()) == "pass"

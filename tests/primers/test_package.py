from liulab_mbio import primers

#: Every public name `liulab_mbio.primers` defined while it was one module.
PUBLIC_NAMES = (
    "ONETAQ",
    "PHUSION",
    "Q5",
    "TAQ",
    "TARGET_TM",
    "THRESHOLDS",
    "Band",
    "PairReport",
    "Polymerase",
    "PrimerReport",
    "PrimingSite",
    "Reading",
    "Thresholds",
    "amplicon_sizes",
    "design_pair",
    "design_primer",
    "evaluate_pair",
    "evaluate_primer",
    "find_binding_sites",
    "find_priming_sites",
    "melting_temperature",
    "reading",
)


def test_every_public_name_still_imports_from_the_package() -> None:
    assert [name for name in PUBLIC_NAMES if not hasattr(primers, name)] == []

"""A number attaches to the thing that changes it, and a rule travels with its material."""

import pytest

from mbio.bench import materials


def test_the_electroporation_program_belongs_to_the_cells_and_not_to_the_step() -> None:
    neb, endura = materials.electroporation("#C3020"), materials.electroporation("60242-2")
    assert neb is not None
    assert endura is not None
    settings = (
        (neb.volts, neb.ohms, neb.microfarads),
        (endura.volts, endura.ohms, endura.microfarads),
    )
    assert settings == ((2000.0, 200.0, 25.0), (1800.0, 600.0, 10.0))
    # They disagree on every setting, which is the point: swapping the cells changes the number
    # while the step stays as it was.
    assert not any(one == other for one, other in zip(*settings, strict=True))
    # Both manuals give the same gap, NEB's in millimetres and Endura's as 0.1 cm.
    assert (neb.cuvette_mm, endura.cuvette_mm) == (1.0, 1.0)


def test_cells_nobody_shipped_a_program_for_have_none_rather_than_a_guess() -> None:
    assert materials.electroporation("#C3040H") is None


@pytest.mark.parametrize("catalog", ["M0318", "#M0318S", "M0318L"])
def test_a_pack_size_changes_nothing_about_the_thing_in_the_tube(catalog: str) -> None:
    assert [rule.subject for rule in materials.rules(catalog)] == ["PEG", "heat inactivation"]


def test_the_ligase_carries_both_peg_rules_wherever_it_is_used() -> None:
    ligase = materials.material("T7 DNA Ligase", catalog="#M0318L")
    forbidding, inactivating = ligase.rules
    assert forbidding.kind == "forbids"
    assert "Never add PEG" in forbidding.detail
    assert inactivating.when == "PEG"
    assert inactivating.citation is not None
    # The PEG is the buffer's, not the enzyme's, so the condition can tell them apart.
    assert ligase.contains == ()
    assert materials.material("buffer", catalog="#B0535S").contains == ("PEG 6000",)


def test_a_material_nothing_rules_carries_no_rule() -> None:
    assert materials.material("Water").rules == ()


@pytest.mark.parametrize(
    ("catalog", "carried"),
    [
        ("M0491", materials.POLYMERASE_ON_ICE),
        ("#M0267S", materials.POLYMERASE_ON_ICE),
        ("60242-2", materials.CUVETTE_ON_ICE),
        ("#C3020", materials.CUVETTE_ON_ICE),
    ],
)
def test_a_caution_is_one_sentence_the_catalogue_number_carries(catalog: str, carried: str) -> None:
    assert materials.material("the tube", catalog=catalog).cautions == (carried,)


def test_a_material_nothing_cautions_carries_no_caution() -> None:
    assert materials.material("Water").cautions == ()


def test_a_kit_named_by_catalogue_number_brings_its_product_name_and_supplier() -> None:
    named = materials.kit("#T1130S")
    assert (named.name, named.supplier, named.catalog) == (
        "Monarch Spin PCR & DNA Cleanup Kit",
        "New England Biolabs",
        "T1130S",
    )


def test_a_kit_nobody_here_lists_is_named_as_given_and_invents_no_number() -> None:
    named = materials.kit("Wizard SV Gel and PCR Clean-Up System")
    assert (named.name, named.supplier, named.catalog) == (
        "Wizard SV Gel and PCR Clean-Up System",
        "",
        "",
    )


def test_the_default_kit_is_the_one_the_packages_own_numbers_were_read_from() -> None:
    assert materials.kit().catalog == materials.DEFAULT_CLEANUP_KIT == "T1130"

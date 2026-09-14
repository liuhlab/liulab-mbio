"""The shipped data file, read through the loader."""

import pytest

from liulab_mbio.enzymes import Enzyme, enzymes, get_enzyme

# The PoC set: the pUC19 multiple cloning site, the Type IIS enzymes used for Golden Gate, and
# the two blunt eight-base cutters the iterative library scheme needs.
EXPECTED = {
    "AarI", "BamHI", "BbsI", "BpiI", "BsaI", "BsmBI", "BspQI", "BtgZI", "EcoRI", "Esp3I",
    "HindIII", "KpnI", "NcoI", "NdeI", "NotI", "PaqCI", "PmeI", "PstI", "SacI", "SalI",
    "SapI", "SbfI", "SmaI", "SphI", "SrfI", "XbaI", "XhoI", "XmaI",
}  # fmt: skip


def test_the_data_file_ships_the_whole_poc_set() -> None:
    assert {enzyme.name for enzyme in enzymes()} == EXPECTED


def test_bsai_carries_every_field_the_bench_needs() -> None:
    bsai = get_enzyme("BsaI")
    assert (bsai.site, bsai.cut_offsets) == ("GGTCTC", (7, 11))  # REBASE GGTCTC(1/5)
    assert (bsai.type, bsai.end, bsai.overhang_length) == ("IIS", "5'", 4)
    assert (bsai.commercial_name, bsai.catalog_number) == ("BsaI-HFv2", "R3733")
    assert bsai.supplier == "New England Biolabs"
    assert (bsai.incubation_celsius, bsai.heat_inactivation_celsius) == (37, 80)
    assert bsai.methylation["dam"] == "not sensitive"
    assert "Eco31I" in bsai.isoschizomers


def test_sapi_leaves_the_three_base_overhang_golden_gate_designs_around() -> None:
    sapi = get_enzyme("SapI")
    assert (sapi.site, sapi.cut_offsets) == ("GCTCTTC", (8, 11))  # REBASE GCTCTTC(1/4)
    assert (sapi.end, sapi.overhang_length) == ("5'", 3)


def test_the_two_blunt_eight_base_cutters_cut_inside_their_site() -> None:
    srfi, pmei = get_enzyme("SrfI"), get_enzyme("PmeI")
    assert (srfi.site, srfi.cut_offsets) == ("GCCCGGGC", (4, 4))  # REBASE GCCC^GGGC
    assert (pmei.site, pmei.cut_offsets) == ("GTTTAAAC", (4, 4))  # REBASE GTTT^AAAC
    assert (srfi.catalog_number, pmei.catalog_number) == ("R0629", "R0560")
    for enzyme in (srfi, pmei):
        assert (enzyme.type, enzyme.end, enzyme.overhang_length) == ("II", "blunt", 0)
        assert enzyme.unverified == ()


def test_a_neoschizomer_is_not_listed_as_an_isoschizomer() -> None:
    # SmaI and XmaI share CCCGGG and cut it differently.
    assert "XmaI" not in get_enzyme("SmaI").isoschizomers
    assert get_enzyme("SmaI").end == "blunt"
    assert get_enzyme("XmaI").end == "5'"


def test_an_enzyme_is_found_by_commercial_name() -> None:
    assert get_enzyme("BsaI-HFv2") is get_enzyme("BsaI")
    assert get_enzyme("BsmBI-v2").incubation_celsius == 55


def test_an_enzyme_is_found_by_an_isoschizomer_nobody_in_the_set_is_named_after() -> None:
    assert get_enzyme("Eco31I") is get_enzyme("BsaI")


def test_a_name_in_the_set_beats_the_same_name_as_an_isoschizomer() -> None:
    # BpiI is BbsI's isoschizomer and ships as its own record, with its own supplier.
    assert get_enzyme("BpiI").supplier == "Thermo Fisher Scientific"
    assert get_enzyme("BbsI").supplier == "New England Biolabs"


def test_looking_up_is_not_case_sensitive() -> None:
    assert get_enzyme("bsai") is get_enzyme("BsaI")


def test_an_isoschizomer_of_two_shipped_enzymes_is_refused_rather_than_guessed() -> None:
    # BstV2I reads GAAGAC like both BbsI and BpiI, which differ in supplier and properties.
    with pytest.raises(KeyError, match="BbsI"):
        get_enzyme("BstV2I")


def test_an_unknown_name_is_refused() -> None:
    with pytest.raises(KeyError, match="NoSuchI"):
        get_enzyme("NoSuchI")


def test_an_enzyme_with_no_stated_heat_inactivation_says_so() -> None:
    # NEB's page says KpnI-HF is not inactivated by heat; Thermo states none for AarI.
    assert get_enzyme("KpnI").heat_inactivation_celsius is None
    assert "heat_inactivation_celsius" not in get_enzyme("KpnI").unverified
    assert "heat_inactivation_celsius" in get_enzyme("AarI").unverified


def test_every_record_is_self_consistent() -> None:
    for enzyme in enzymes():
        assert isinstance(enzyme, Enzyme)
        assert set(enzyme.methylation) == {"dam", "dcm", "cpg"}
        assert enzyme.overhang_length == abs(enzyme.bottom_cut - enzyme.top_cut)
        for isoschizomer in enzyme.isoschizomers:
            assert isoschizomer != enzyme.name

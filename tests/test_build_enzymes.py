"""The regeneration script's parsers, on excerpts rather than on the network."""

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO = Path(__file__).resolve().parents[1]


def _script() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "build_enzymes", REPO / "scripts/build_enzymes.py"
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # Registered before it runs: `dataclasses` resolves annotations through `sys.modules`.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


build = _script()

# Three records as REBASE writes them in withrefm, cut down to the fields the script reads.
WITHREFM = """
REBASE version 609                                              withrefm.609

<1>BsaI
<2>Eco31I,Bso31I,Uba1216I
<3>GGTCTC(1/5)
<4>-4(6)
<5>Bacillus stearothermophilus 6-55
<6>Z. Chen
<7>NSV
<8>Zhu, Z., Xu, S.-Y., Unpublished observations.

<1>SmaI
<2>Cfr9I,XmaI
<3>CCC^GGG
<4>2(4)
<5>Serratia marcescens Sb
<6>ATCC 25419
<7>BIJNQRSVX
<8>Endow, S.A., Roberts, R.J., (1977) J. Mol. Biol., vol. 112, pp. 521-529.

<1>XmaI
<2>Cfr9I,SmaI
<3>C^CCGGG
<4>2(4)
<5>Xanthomonas malvacearum
<6>ATCC 9924
<7>INRV
<8>Endow, S.A., Roberts, R.J., (1977) J. Mol. Biol., vol. 112, pp. 521-529.

<1>Uba1216I
<2>BsaI
<3>GGTCTC(1/5)
<4>
<5>Unidentified bacterium
<6>I. Kobayashi
<7>
<8>Unpublished observations.

<1>Eco31I
<2>BsaI,Bso31I,Uba1216I
<3>GGTCTC(1/5)
<4>-4(6)
<5>Escherichia coli RFL31
<6>A. Janulaitis
<7>BV
<8>Kriukiene, E., Lubiene, J., Lagunavicius, A., Lubys, A., Unpublished observations.
"""


def test_a_type_iis_notation_puts_both_cuts_outside_the_site() -> None:
    assert build.parse_site("GGTCTC(1/5)") == ("GGTCTC", 7, 11)
    assert build.parse_site("GCTCTTC(1/4)") == ("GCTCTTC", 8, 11)


def test_a_caret_cuts_both_strands_symmetrically() -> None:
    assert build.parse_site("CCC^GGG") == ("CCCGGG", 3, 3)
    assert build.parse_site("C^CCGGG") == ("CCCGGG", 1, 5)
    assert build.parse_site("GGTAC^C") == ("GGTACC", 5, 1)


def test_a_negative_offset_counts_back_from_the_end_of_the_site() -> None:
    assert build.parse_site("AACTGG(-5/-1)") == ("AACTGG", 1, 5)


@pytest.mark.parametrize("notation", ["?", "GAATTC", "(8/13)GACNNNNNNTGG(12/7)"])
def test_a_site_whose_cuts_rebase_does_not_give_is_refused(notation: str) -> None:
    with pytest.raises(ValueError, match="cut"):
        build.parse_site(notation)


def test_reading_withrefm_keeps_the_fields_the_data_file_needs() -> None:
    records = build.parse_withrefm(WITHREFM)
    assert records["BsaI"].site == "GGTCTC"
    assert records["BsaI"].cuts == (7, 11)
    assert records["BsaI"].suppliers == "NSV"
    assert records["BsaI"].same_specificity == ("Eco31I", "Bso31I", "Uba1216I")


def test_an_isoschizomer_is_sold_by_someone_and_cuts_the_same_way() -> None:
    records = build.parse_withrefm(WITHREFM)
    # Uba1216I cuts as BsaI does but nobody sells it; Eco31I is sold by two suppliers.
    assert build.isoschizomers("BsaI", records) == ("Eco31I",)
    # XmaI carries SmaI's site but cuts it elsewhere, so it is no isoschizomer of SmaI.
    assert build.isoschizomers("SmaI", records) == ()

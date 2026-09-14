import pytest

from liulab_mbio.enzymes import Enzyme
from liulab_mbio.sequence import SequenceRecord, Strand
from liulab_mbio.sites import has_site
from liulab_mbio.translate import (
    SiteNotRemovableError,
    optimize_coding_sequence,
    optimize_protein,
    reverse_translate,
    translate,
)

# The only table on this branch. The host is a parameter, so another table flows through it.
HOST = "e-coli-k12"

# M G L K T *, holding one BsaI site. The leucine codon E. coli prefers, CTG, takes the site
# away and spells a BbsI site that was not there before.
TRAP = "ATGGGTCTCAAGACCTAA"

# M L L S G *, written with codons E. coli uses rarely and holding no site anyone forbids here.
CODED = "ATGTTACTATCCGGATAA"

# Two tryptophan codons, which have no synonym between them.
TRPI = Enzyme("TrpI", "TGGTGG", top_cut=3, bottom_cut=3)

# Every glycine codon begins GG, so no synonymous change escapes this pattern.
BLOCKI = Enzyme("BlockI", "GGNNNN", top_cut=3, bottom_cut=3)


def test_a_protein_is_written_with_the_codons_the_host_counts_most() -> None:
    assert reverse_translate("MW*", host=HOST) == "ATGTGGTAA"
    # E. coli spells leucine CTG far more often than the other five codons.
    assert reverse_translate("L", host=HOST) == "CTG"
    assert translate("ATGTGGTAA") == "MW*"


def test_the_protein_round_trips_through_the_coding_sequence() -> None:
    protein = "MSGRWKLQNHGVDE*"
    gene = optimize_protein(protein, host=HOST, forbidden=["NotI", "PstI", "NdeI"])
    assert translate(gene.dna) == protein
    assert gene.protein == protein
    assert len(gene.dna) == 3 * len(protein)
    assert not has_site(SequenceRecord(gene.dna), ["NotI", "PstI", "NdeI"])


def test_a_forbidden_site_the_codons_spell_is_taken_away() -> None:
    plain = reverse_translate("MSGRW*", host=HOST)
    assert has_site(SequenceRecord(plain), "NotI")

    gene = optimize_protein("MSGRW*", host=HOST, forbidden=["NotI"], name="part A")
    assert not has_site(SequenceRecord(gene.dna), "NotI")
    assert gene.protein == "MSGRW*"
    assert (gene.name, gene.host, gene.forbidden, gene.recoded) == (
        "part A",
        HOST,
        ("NotI",),
        True,
    )
    (change,) = gene.changes
    assert (change.old_codon, change.new_codon, change.amino_acid) == ("CGC", "CGT", "R")
    assert change.site.enzyme.name == "NotI"
    assert gene.protein[change.codon_index] == change.amino_acid


def test_taking_one_site_away_does_not_spell_another() -> None:
    loose = optimize_coding_sequence(TRAP, host=HOST, forbidden=["BsaI"])
    (only,) = loose.changes
    assert (only.old_codon, only.new_codon) == ("CTC", "CTG")
    assert has_site(SequenceRecord(loose.dna), "BbsI")

    gene = optimize_coding_sequence(TRAP, host=HOST, forbidden=["BsaI", "BbsI"])
    (change,) = gene.changes
    assert (change.old_codon, change.new_codon, change.amino_acid) == ("GGT", "GGC", "G")
    assert not has_site(SequenceRecord(gene.dna), ["BsaI", "BbsI"])
    assert gene.protein == translate(TRAP)


def test_a_site_on_the_reverse_strand_is_found_and_taken_away() -> None:
    gene = optimize_coding_sequence("ATGGAGACCTAA", host=HOST, forbidden=["BsaI"])
    (change,) = gene.changes
    assert change.site.strand == Strand.REVERSE
    assert (change.old_codon, change.new_codon, change.amino_acid) == ("GAG", "GAA", "E")
    assert not has_site(SequenceRecord(gene.dna), "BsaI")


def test_a_coded_sequence_is_checked_and_not_written_again() -> None:
    gene = optimize_coding_sequence(CODED, host=HOST, forbidden=["BsaI", "BbsI", "NotI"])
    assert (gene.dna, gene.changes, gene.recoded) == (CODED, (), False)
    assert gene.protein == "MLLSG*"
    # Written from the protein instead, the same host would have chosen other codons.
    assert reverse_translate(gene.protein, host=HOST) != CODED


def test_a_site_whose_amino_acids_have_one_codon_each_is_refused_by_name() -> None:
    with pytest.raises(SiteNotRemovableError, match="part A: the TrpI site at 3-9 on the forward"):
        optimize_protein("MWW*", host=HOST, forbidden=[TRPI], name="part A")


def test_a_site_every_synonymous_change_keeps_is_refused_too() -> None:
    with pytest.raises(SiteNotRemovableError, match="every synonymous change there spells another"):
        optimize_coding_sequence(TRAP, host=HOST, forbidden=["BsaI", BLOCKI])


def test_the_same_request_writes_the_same_sequence() -> None:
    forbidden = ["NotI", "BsaI"]
    assert optimize_protein("MSGRW*", host=HOST, forbidden=forbidden) == optimize_protein(
        "MSGRW*", host=HOST, forbidden=forbidden
    )


def test_a_letter_that_is_not_an_amino_acid_is_refused() -> None:
    with pytest.raises(ValueError, match="BJ"):
        reverse_translate("MBJW", host=HOST)


def test_a_length_that_is_not_whole_codons_is_refused() -> None:
    with pytest.raises(ValueError, match="not a whole number of codons"):
        optimize_coding_sequence("ATGGT", host=HOST)


def test_a_base_that_is_not_definite_is_refused() -> None:
    with pytest.raises(ValueError, match="not one of ACGT"):
        optimize_coding_sequence("ATGNNNTAA", host=HOST)


def test_a_stop_before_the_end_is_refused() -> None:
    with pytest.raises(ValueError, match="truncate the protein"):
        optimize_coding_sequence("ATGTAAGGCTAA", host=HOST)
    with pytest.raises(ValueError, match="truncate the protein"):
        reverse_translate("M*G*", host=HOST)


def test_a_host_the_package_does_not_ship_is_refused() -> None:
    with pytest.raises(KeyError, match="martian"):
        optimize_protein("MW*", host="martian")


def test_an_enzyme_the_package_does_not_ship_is_refused() -> None:
    with pytest.raises(KeyError, match="NotAnEnzyme"):
        optimize_protein("MW*", host=HOST, forbidden=["NotAnEnzyme"])

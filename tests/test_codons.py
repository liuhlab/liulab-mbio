import pytest

from liulab_mbio.codons import codon_tables, codon_usage


def test_the_shipped_table_counts_the_whole_e_coli_k12_genome() -> None:
    usage = codon_usage()
    assert (usage.name, usage.accession) == ("e-coli-k12", "U00096.3")
    # Every complete coding sequence of MG1655, counted twice by hand before it was shipped.
    assert (usage.cds_count, usage.codon_count) == (4317, 1342016)
    assert usage.codon_count == sum(usage.counts.values())


def test_the_mammal_tables_count_one_coding_sequence_for_each_protein_coding_gene() -> None:
    human, mouse = codon_usage("human"), codon_usage("mouse")
    # GENCODE's canonical set, less the coding sequences it marks as running off an end.
    assert (human.cds_count, human.codon_count) == (19597, 11327553)
    assert (mouse.cds_count, mouse.codon_count) == (21479, 11894511)
    assert (human.accession, mouse.accession) == ("GCF_000001405.40", "GCF_000001635.27")


def test_all_sixty_four_codons_are_used_somewhere_in_every_genome_shipped() -> None:
    # A zero cell makes a fraction undefined, which is what sinks the small published tables.
    for name in codon_tables():
        counts = codon_usage(name).counts
        assert len(counts) == 64
        assert min(counts.values()) > 0


def test_which_synonymous_codon_is_favoured_depends_on_the_host() -> None:
    # What domestication turns on: E. coli spells aspartate GAT, and a mammal spells it GAC.
    assert codon_usage().synonymous("GAC")[0] == "GAT"
    assert codon_usage("human").synonymous("GAC")[0] == "GAC"
    assert codon_usage("mouse").synonymous("GAC")[0] == "GAC"


def test_a_codon_reads_back_as_the_amino_acid_it_spells() -> None:
    usage = codon_usage()
    assert (usage.amino_acid("ATG"), usage.amino_acid("TGG")) == ("M", "W")
    assert usage.amino_acid("TAA") == "*"


def test_the_fraction_of_a_codon_is_its_share_of_its_own_amino_acid() -> None:
    usage = codon_usage()
    assert usage.fraction("GAT") + usage.fraction("GAC") == pytest.approx(1.0)
    # E. coli spells aspartate GAT more often than GAC.
    assert usage.fraction("GAT") > usage.fraction("GAC")
    assert usage.fraction("ATG") == 1.0


def test_per_thousand_shares_out_every_codon_counted() -> None:
    usage = codon_usage()
    assert sum(usage.per_thousand(codon) for codon in usage.counts) == pytest.approx(1000.0)


def test_synonymous_codons_come_back_with_the_most_used_first() -> None:
    usage = codon_usage()
    assert usage.synonymous("GAC") == ("GAT", "GAC")
    assert usage.synonymous("GGG")[0] == "GGC"
    assert set(usage.synonymous("TAA")) == {"TAA", "TAG", "TGA"}


def test_the_tables_shipped_are_listed_by_name() -> None:
    assert "e-coli-k12" in codon_tables()


def test_an_unknown_table_is_refused() -> None:
    with pytest.raises(KeyError, match="martian"):
        codon_usage("martian")


def test_a_codon_that_is_not_one_is_refused() -> None:
    with pytest.raises(KeyError, match="AT"):
        codon_usage().fraction("AT")

"""The codon usage build script's counting, on an excerpt rather than on the network."""

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

REPO = Path(__file__).resolve().parents[2]


def _script() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "build_codon_usage", REPO / "scripts/build_codon_usage.py"
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # Registered before it runs: `dataclasses` resolves annotations through `sys.modules`.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


build = _script()

# Two complete coding sequences, one truncated mid-codon, one with an ambiguous base, written
# the way NCBI's `fasta_cds_na` writes them.
CDS = """>lcl|U00096.3_cds_AAC73112.1_1 [gene=thrA]
ATGAAAGTGTAA
>lcl|U00096.3_cds_AAC73113.1_2 [gene=thrB]
ATGAAAAAATGA
>lcl|U00096.3_cds_AAC73114.1_3 [partial=true]
ATGAA
>lcl|U00096.3_cds_AAC73115.1_4 [gene=unknown]
ATGNNNTAA
"""


def test_every_complete_coding_sequence_is_counted() -> None:
    counts, used = build.count_codons(CDS)
    assert used == 2
    assert counts["ATG"] == 2
    assert counts["AAA"] == 3
    assert counts["GTG"] == 1
    assert (counts["TAA"], counts["TGA"]) == (1, 1)


def test_a_sequence_that_is_not_whole_codons_is_skipped_rather_than_counted_in_part() -> None:
    assert build.count_codons(">short\nATGAA\n") == ({}, 0)


def test_a_sequence_carrying_an_ambiguous_base_is_skipped() -> None:
    assert build.count_codons(">ambiguous\nATGNNNTAA\n") == ({}, 0)


def test_a_sequence_wrapped_over_several_lines_is_read_as_one() -> None:
    counts, used = build.count_codons(">wrapped\nATGAAA\nGTGTAA\n")
    assert (used, counts["GTG"]) == (1, 1)


def test_every_codon_gets_a_cell_even_where_the_genome_never_spells_it() -> None:
    table = build.build({host.name: CDS for host in build.HOSTS})["tables"][0]
    assert len(table["counts"]) == 64
    assert table["counts"]["CGG"] == 0
    assert (table["cds_count"], table["codon_count"]) == (2, 8)


# A GENCODE excerpt: a canonical two-exon transcript on the minus strand with its stop codon on
# its own row, another isoform of the same gene, a canonical CDS running off its start, and a
# mitochondrial one.
CANONICAL = ' tag "Ensembl_canonical";'


def _attributes(transcript: str, extra: str = "") -> str:
    return f'gene_id "G1"; transcript_id "{transcript}"; transcript_type "protein_coding";{extra}'


GTF = "##description: excerpt\n" + "\n".join(
    "\t".join(row)
    for row in (
        ("chr1", "HAVANA", "transcript", "8", "26", ".", "-", ".", _attributes("T1", CANONICAL)),
        ("chr1", "HAVANA", "CDS", "21", "26", ".", "-", "0", _attributes("T1", CANONICAL)),
        ("chr1", "HAVANA", "CDS", "11", "16", ".", "-", "0", _attributes("T1", CANONICAL)),
        ("chr1", "HAVANA", "stop_codon", "8", "10", ".", "-", "0", _attributes("T1", CANONICAL)),
        ("chr1", "HAVANA", "CDS", "11", "26", ".", "-", "0", _attributes("T2")),
        (
            "chr1",
            "HAVANA",
            "CDS",
            "100",
            "105",
            ".",
            "+",
            "0",
            _attributes("T3", CANONICAL + ' tag "cds_start_NF";'),
        ),
        ("chrM", "ENSEMBL", "CDS", "1", "6", ".", "+", "0", _attributes("T4", CANONICAL)),
    )
)

# The plus strand under T1: the stop codon, the second exon, an intron, then the first exon.
CHR1 = "CCCCCCC" + "TTA" + "GGGCAC" + "GGGG" + "TTTCAT" + "CCCC"


def _fetch(chrom: str, start: int, end: int, strand: str) -> str:
    bases = CHR1[start:end].lower()
    return bases[::-1].translate(str.maketrans("acgt", "tgca")) if strand == "-" else bases


def test_only_whole_canonical_nuclear_coding_sequences_are_read_from_a_gtf() -> None:
    assert build.canonical_cds(GTF.splitlines()) == {
        "T1": ("chr1", "-", [(20, 26), (10, 16), (7, 10)])
    }


def test_a_minus_strand_transcript_is_spliced_five_prime_to_three_prime() -> None:
    _, _, spans = build.canonical_cds(GTF.splitlines())["T1"]
    assert build.spliced(_fetch, "chr1", "-", spans) == "ATGAAAGTGCCCTAA"


def test_the_built_table_names_the_genome_it_counted() -> None:
    built = build.build({host.name: CDS for host in build.HOSTS})
    assert built["tables"][0]["accession"] == "U00096.3"
    assert "no restrictions" in built["source"]["terms"]

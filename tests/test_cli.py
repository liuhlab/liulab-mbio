"""The root command line: the version and codon-optimise verbs, and that both sub-apps mount."""

import re

from typer.testing import CliRunner

from liulab_mbio import __version__
from liulab_mbio.cli import app


def test_version_comes_from_the_installed_metadata() -> None:
    # hatch-vcs derives it from git, so the value moves. That it exists is the claim.
    assert __version__


def test_the_version_verb_prints_it() -> None:
    result = CliRunner().invoke(app, ["version"])
    assert result.exit_code == 0
    assert result.output.strip() == __version__


def test_a_bare_invocation_prints_help_rather_than_nothing() -> None:
    result = CliRunner().invoke(app, [])
    assert result.exit_code != 0
    assert result.output.strip()


def test_both_pipelines_are_mounted() -> None:
    for verb in ("goldengate", "protocol"):
        result = CliRunner().invoke(app, [verb])
        # `no_args_is_help` on the sub-app, so it prints its own help and exits non-zero.
        assert result.exit_code != 0
        assert result.output.strip()


def test_the_codon_optimize_verb_writes_a_protein_as_dna() -> None:
    result = CliRunner().invoke(
        app, ["codon-optimize", "MW*", "--kind", "protein", "--host", "e-coli-k12"]
    )
    assert result.exit_code == 0, result.output
    lines = re.sub(r"\x1b\[[0-9;]*m", "", result.output).splitlines()
    assert lines[0].startswith("coding sequence: 9 bp, 3 aa, host e-coli-k12, 0 codon change(s)")
    assert lines[-1] == "ATGTGGTAA"


def test_the_codon_optimize_verb_refuses_a_kind_that_is_neither() -> None:
    result = CliRunner().invoke(
        app, ["codon-optimize", "MW*", "--kind", "rna", "--host", "e-coli-k12"]
    )
    assert result.exit_code == 1
    assert "--kind is 'protein' or 'dna'" in result.output

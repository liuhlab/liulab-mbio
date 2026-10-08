"""The synbio package: that it imports, carries one version with mbio, and its root app runs."""

import re

from typer.testing import CliRunner

import liulab_mbio
import liulab_synbio
from liulab_synbio.cli import app


def test_both_packages_report_the_one_distribution_version() -> None:
    # hatch-vcs derives it from git, so the value moves. That it exists and is shared is
    # the claim: one distribution ships both import packages.
    assert liulab_synbio.__version__
    assert liulab_synbio.__version__ == liulab_mbio.__version__


def test_the_version_verb_prints_the_installed_version() -> None:
    result = CliRunner().invoke(app, ["version"])
    assert result.exit_code == 0
    assert result.output.strip() == liulab_synbio.__version__


def test_a_bare_invocation_prints_help_rather_than_nothing() -> None:
    result = CliRunner().invoke(app, [])
    assert result.exit_code != 0
    assert result.output.strip()


def test_the_help_lists_the_verbs_it_mounts() -> None:
    result = CliRunner().invoke(app, ["--help"])
    assert result.exit_code == 0
    # Typer prints the help through rich, which styles it; the names are what is asserted.
    listed = re.sub(r"\x1b\[[0-9;]*m", "", result.output)
    assert "version" in listed
    assert "igga" in listed
    assert "dmx" in listed


def test_the_igga_pipeline_is_mounted() -> None:
    result = CliRunner().invoke(app, ["igga"])
    # `no_args_is_help` on the sub-app, so it prints its own help and exits non-zero.
    assert result.exit_code != 0
    assert result.output.strip()

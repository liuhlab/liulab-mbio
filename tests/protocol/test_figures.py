"""The figure library chooses a `Figure`'s fields from a design, and never draws one."""

from pathlib import Path

import pytest

from mbio.plot.sequence_view import LIMIT
from mbio.protocol import Figure, Protocol, Step, render_html
from mbio.protocol.figures import (
    LIGATION_CITATION,
    SOURCE,
    SOURCE_KEY,
    ligation_figure,
    tail_figure,
)
from mbio.sequence import BindingSite, Primer, SequenceRecord, Strand

from ..html import parse

PLASMID = SequenceRecord("ACGT" * 25, name="p", topology="circular")


def _figure(**over: object) -> Figure:
    fields: dict[str, object] = {
        "path": "product.dna",
        "junction": (40, 44),
        "enzymes": ["BsaI"],
        "caption": "The entry junction",
        "citation": LIGATION_CITATION,
    }
    return ligation_figure(PLASMID, **(fields | over))  # type: ignore[arg-type]


def test_a_ligation_figure_draws_the_junction_in_context_with_both_strands() -> None:
    """The base-level figure: a span, the bases both ways, the cuts through them and the frame."""
    figure = _figure(context=10)

    assert figure.records == ("product.dna",)
    assert figure.span == (30, 54)
    assert figure.sequence_view
    assert figure.linear
    assert figure.enzymes == ("BsaI",)


def test_a_junction_near_the_origin_of_a_plasmid_runs_across_it() -> None:
    """A span across the origin ends past the length, which is the one coordinate rule."""
    assert _figure(junction=(2, 6), context=10).span == (92, 116)


def test_a_junction_near_the_end_of_a_linear_record_stops_at_it() -> None:
    """Nothing is drawn off a linear molecule, so the context is what is there."""
    linear = SequenceRecord("ACGT" * 25, name="p")
    figure = ligation_figure(
        linear, path="p.dna", junction=(94, 98), enzymes=["BsaI"], caption="The end"
    )
    assert figure.span == (70, 100)


def test_a_ligation_figure_cites_the_note_its_equivalent_was_read_in() -> None:
    """A figure's provenance is a note's: one `Citation` keying one `Source`."""
    one = Protocol("t", steps=(Step("Join", figures=(_figure(),)),), sources={SOURCE_KEY: SOURCE})
    assert one.cited == frozenset({SOURCE_KEY})
    assert next(c for c in one.audit() if c.name == "sources").status == "pass"


def test_a_figure_computed_from_a_design_cites_nothing() -> None:
    """The note covers one method's own ligation, so no other caller is made to cite it."""
    assert _figure(citation=None).citation is None
    assert (
        ligation_figure(
            PLASMID, path="product.dna", junction=(40, 44), enzymes=[], caption="A join"
        ).citation
        is None
    )


def test_context_below_zero_is_refused() -> None:
    with pytest.raises(ValueError, match="0 or more bases"):
        _figure(context=-1)


def test_a_span_longer_than_a_sequence_view_holds_is_refused_where_it_is_chosen() -> None:
    """A junction on a chromosome would draw every base of it, which no page can show."""
    long = SequenceRecord("A" * (LIMIT + 1), name="chromosome")
    with pytest.raises(ValueError, match="at most"):
        ligation_figure(
            long, path="c.dna", junction=(0, LIMIT + 1), enzymes=["BsaI"], caption="All of it"
        )


def test_the_figure_a_library_chose_renders_like_any_other(data_dir: Path) -> None:
    """A library function's only product is a spec, so the page draws it with everything else."""
    from mbio.io import read_record

    record = read_record(data_dir / "pUC19.dna")
    figure = ligation_figure(
        record,
        path="pUC19.dna",
        junction=(410, 414),
        enzymes=["EcoRI"],
        caption="The MCS, cut",
        context=20,
    )
    one = Protocol("Cut", steps=(Step("Cut it", figures=(figure,)),), sources={SOURCE_KEY: SOURCE})

    page = parse(render_html(one, base=data_dir))

    [drawn] = page.find_all("figure", cls="map")
    assert drawn.find_all("svg")


def test_a_tail_figure_draws_the_end_its_primer_makes_and_lights_the_cut() -> None:
    """A reverse primer's end is the record's end, and the primer and the enzyme are lit."""
    forward = Primer("f", "GGTCTCAACGT", binding_sites=(BindingSite(7, 11, Strand.FORWARD),))
    reverse = Primer("r", "GGTCTCATTGC", binding_sites=(BindingSite(29, 33, Strand.REVERSE),))
    amplicon = SequenceRecord("A" * 40, name="a", primers=(forward, reverse))

    def drawn(primer: str) -> Figure:
        return tail_figure(
            amplicon, path="a.dna", primer=primer, enzymes=["BsaI"], caption="An end", context=3
        )

    assert drawn("f").span == (0, 14)
    assert drawn("r").span == (26, 40)
    assert drawn("r").highlight == ("r", "BsaI")
    assert drawn("r").sequence_view
    with pytest.raises(ValueError, match="no primer"):
        drawn("gone")
    # A primer of the same name carried over from a template anneals inside; the one making the
    # end is the one drawn.
    carried = Primer("f", "ACGT", binding_sites=(BindingSite(20, 24, Strand.FORWARD),))
    amplicon = SequenceRecord("A" * 40, name="a", primers=(carried, forward, reverse))
    assert drawn("f").span == (0, 14)

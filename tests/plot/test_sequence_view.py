"""The sequence view laid out, items in and shapes out: rows holding the stretch in order, bars
and names under the bases, translations centred on their codons, primers beside the strand they
spell with their tails and mismatches, cuts through both strands, and labels apart from each other
and clear of the bases, on crowded and seeded random records, whole and in regions."""

import math
import random
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from itertools import combinations, pairwise

import pytest

from liulab_mbio.plot import layers, sequence_view, svg
from liulab_mbio.plot.fonts import BOLD, SANS
from liulab_mbio.plot.labels import Box, Point
from liulab_mbio.sequence import (
    BindingSite,
    Feature,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)

LENGTH = 1200
#: How far two shapes may reach into each other and still count as touching.
TOUCH = 1e-6


@dataclass(frozen=True)
class Given:
    """A record's items and bases, and how the sequence view draws them."""

    items: tuple[layers.Item, ...]
    bases: str
    circular: bool = True
    span: tuple[int, int] | None = None
    bases_per_row: int = 60
    both_strands: bool = True


def _laid(given: Given) -> sequence_view.SequenceView:
    return sequence_view.layout(
        given.items,
        bases=given.bases,
        circular=given.circular,
        span=given.span,
        bases_per_row=given.bases_per_row,
        both_strands=given.both_strands,
    )


def _record(
    *features: Feature, bases: str = "", circular: bool = True, primers: tuple[Primer, ...] = ()
) -> SequenceRecord:
    return SequenceRecord(
        bases or "ACGT" * (LENGTH // 4),
        topology="circular" if circular else "linear",
        features=features,
        primers=primers,
    )


def _items(record: SequenceRecord) -> tuple[layers.Item, ...]:
    return layers.items(record, cut_sites=False, translations=True)


def _feature(name: str, start: int, end: int, strand: Strand = Strand.FORWARD) -> Feature:
    return Feature(name, "misc_feature", (Segment(start, end),), strand=strand)


def _bases(seed: int = 1) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(LENGTH))


def _primer(
    name: str,
    bases: str,
    start: int,
    end: int,
    strand: Strand = Strand.FORWARD,
    *,
    tail: str = "",
    changed: tuple[int, ...] = (),
) -> Primer:
    """A primer binding `bases` from `start` to `end`, a `tail` 5' of it, and its base changed at
    each offset from its 3' end in `changed`."""
    site = "".join(bases[at % len(bases)] for at in range(start, end))
    letters = list(tail + (reverse_complement(site) if strand == Strand.REVERSE else site))
    for offset in changed:
        letters[-1 - offset] = {"A": "C", "C": "G", "G": "T", "T": "A"}.get(
            letters[-1 - offset], "A"
        )
    return Primer(name, "".join(letters), binding_sites=(BindingSite(start, end, strand),))


def _crowded() -> Given:
    """Names crowding a few rows, over features packed in tracks, CDS joined across the origin,
    primers overlapping with tails and mismatches, and cut sites a base or two apart."""
    rng = random.Random(1)
    bases = "".join(rng.choice("ACGT") for _ in range(LENGTH))
    features = [_feature(f"crowded name {i:02d}", 300 + 2 * i, 303 + 2 * i) for i in range(30)]
    features += [
        _feature(f"origin {i}", (1190 + 2 * i) % LENGTH, (1190 + 2 * i) % LENGTH + 3)
        for i in range(12)
    ]
    features += [_feature(f"long {i}", 40 * i, 40 * i + 500, Strand.REVERSE) for i in range(6)]
    features += [
        Feature("both ways", "misc_feature", (Segment(1100, 1300),), strand=Strand.BOTH),
        Feature(
            "joined across",
            "CDS",
            (Segment(1150, 1180), Segment(1190, 1220), Segment(20, 61)),
            strand=Strand.FORWARD,
            qualifiers={"codon_start": (2,)},
        ),
        Feature(
            "reverse coding",
            "CDS",
            (Segment(290, 350), Segment(400, 700)),
            strand=Strand.REVERSE,
        ),
        _feature("a name far longer than its feature", 800, 801),
        _feature("fits inside easily", 850, 1000),
    ]
    primers = [
        _primer(
            f"primer {i}",
            bases,
            290 + 3 * i,
            310 + 3 * i,
            Strand.REVERSE if i % 3 == 2 else Strand.FORWARD,
            tail="GGTCTCA" * (i % 3),
            changed=tuple(range(0, 20, 4 + i)),
        )
        for i in range(10)
    ]
    primers += [
        _primer("across", bases, 1185, 1215, tail="A" * 30, changed=(2,)),
        _primer("back across", bases, 1195, 1205, Strand.REVERSE, tail="C" * 12),
    ]
    names = ["Cutter", "Other", "Third", "Fourth", "Fifth"]
    cuts = [(names[i % 5], 300 + i + i // 3) for i in range(20)]
    cuts += [("Cutter", 0), ("Other", 1199), ("Third", 1199), ("Long name for a cutter", 310)]
    staggers = {"Cutter": 4, "Other": -4, "Third": 0, "Fourth": 2, "Fifth": -20}
    return Given(
        (
            *_items(_record(*features, bases=bases, primers=tuple(primers))),
            *layers.merge_cuts(cuts, LENGTH, staggers=staggers),
        ),
        bases,
    )


def _random(seed: int, circular: bool) -> tuple[tuple[layers.Item, ...], str]:
    """Features with random names, strands, types and joins, primers with random tails and
    mismatches, and cut sites with random staggers, half of each clustered, on random bases."""
    rng = random.Random(seed)
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789 -()"
    bases = "".join(rng.choice("ACGTACGTN") for _ in range(LENGTH))
    centres = [rng.randrange(LENGTH) for _ in range(3)]

    def span(size: int) -> tuple[int, int]:
        start = rng.randrange(LENGTH)
        if rng.random() < 0.5:
            start = (rng.choice(centres) + rng.randint(-20, 20)) % LENGTH
        if not circular:
            start = min(start, LENGTH - size)
        return start, start + size

    def name() -> str:
        return "".join(rng.choice(letters) for _ in range(rng.randint(1, 30)))

    features = []
    for i in range(rng.randint(20, 80)):
        start, end = span(rng.choice([rng.randint(1, 6), rng.randint(6, 30), rng.randint(30, 600)]))
        segments = (Segment(start, end),)
        if end - start > 60 and rng.random() < 0.3:
            last = (end - 20) % LENGTH
            segments = (Segment(start, start + 20), Segment(last, last + 20))
        features.append(
            Feature(
                f"{i}{name()}",
                rng.choice(["CDS", "misc_feature"]),
                segments,
                strand=rng.choice(list(Strand)),
                qualifiers={"codon_start": (rng.randint(1, 3),)},
            )
        )
    primers = []
    for i in range(rng.randint(0, 20)):
        start, end = span(rng.randint(8, 40))
        tail = "".join(rng.choice("ACGT") for _ in range(rng.choice([0, rng.randint(1, 50)])))
        changed = tuple(rng.sample(range(end - start), rng.randint(0, 3)))
        strand = rng.choice([Strand.FORWARD, Strand.REVERSE])
        primers.append(
            _primer(f"{i}{name()}", bases, start, end, strand, tail=tail, changed=changed)
        )
    cutters = [name() for _ in range(6)]
    cuts = [(rng.choice(cutters), span(0)[0]) for _ in range(rng.randint(0, 30))]
    staggers = {cutter: rng.randint(-5, 5) for cutter in cutters}
    record = _record(*features, bases=bases, circular=circular, primers=tuple(primers))
    return (*_items(record), *layers.merge_cuts(cuts, LENGTH, staggers=staggers)), bases


def _random_given(seed: int) -> Given:
    """Random items drawn whole or in a region, across the origin too, in rows of any size."""
    rng = random.Random(-seed)
    circular = seed % 3 != 1
    items, bases = _random(seed, circular)
    span = None
    if seed % 3 == 2:
        start = rng.randrange(LENGTH)
        span = (start, start + rng.randint(1, LENGTH))
    return Given(
        items,
        bases,
        circular,
        span,
        bases_per_row=(45, 100, 1, 60, 13, 250)[seed],
        both_strands=rng.random() < 0.7,
    )


#: One record per stretch the view draws, per strands it shows, and per crowd no other record runs
#: the layout through: a circular record whole, one across its origin, one region inside it of a
#: base a row, a linear record, and a crowd on one strand in rows of seven.
RECORDS: dict[str, Callable[[], Given]] = {
    "crowded": _crowded,
    "crowded, across the origin": lambda: Given(
        _crowded().items, _crowded().bases, span=(1100, 1400)
    ),
    "crowded, one strand, 7 a row": lambda: Given(
        _crowded().items, _crowded().bases, bases_per_row=7, both_strands=False
    ),
    **{f"random {seed}": lambda seed=seed: _random_given(seed) for seed in (1, 2)},
}


@pytest.fixture(scope="module", params=list(RECORDS))
def given(request: pytest.FixtureRequest) -> Given:
    return RECORDS[request.param]()


@pytest.fixture(scope="module")
def laid_out(given: Given) -> sequence_view.SequenceView:
    return _laid(given)


def _overlap(one: Box, other: Box) -> bool:
    return (
        one.x < other.x + other.width - TOUCH
        and other.x < one.x + one.width - TOUCH
        and one.y < other.y + other.height - TOUCH
        and other.y < one.y + one.height - TOUCH
    )


def _crosses(leader: tuple[Point, Point], box: Box) -> bool:
    """Whether a vertical leader passes through the box's inside."""
    (x, y0), (_, y1) = leader
    low, high = sorted((y0, y1))
    return (
        box.x + TOUCH < x < box.x + box.width - TOUCH
        and low < box.y + box.height - TOUCH
        and box.y + TOUCH < high
    )


def _boxes(letters: svg.Letters, *, words: bool = True) -> list[tuple[str, Box]]:
    """Each word of a line of letters, or the whole line, and the box it takes: from its first
    letter's start to its last's end, as tall as the face."""
    font, scale = letters.font, letters.size / letters.font.units_per_em
    drawn, places = font.letters(letters.text, letters.size), letters.places
    assert len(drawn) == len(places)
    assert {place.rotate for place in places} == {0.0}
    assert len({place.y for place in places}) == 1
    top = places[0].y - font.ascender * scale
    height = (font.ascender - font.descender) * scale
    found, first = [], 0
    for word in letters.text.split(" ") if words else [letters.text]:
        last = first + len(word) - 1 if words else len(places) - 1
        width = places[last].x + drawn[last].advance - places[first].x
        found.append((word, Box(places[first].x, top, width, height)))
        first = last + 2
    return found


def _name_box(name: sequence_view.Name) -> Box:
    [(_, box)] = _boxes(name.letters, words=False)
    return box


def _residues(row: sequence_view.Row) -> list[tuple[str, Box, str]]:
    """Each amino acid a row draws: its three letters, its box and its colour."""
    return [
        (word, box, line.fill)
        for translation in row.translations
        for line in (translation.letters, translation.stops)
        if line
        for word, box in _boxes(line)
    ]


def _band(bar: sequence_view.Bar) -> Box:
    return Box(bar.start, bar.body.y, bar.end - bar.start, bar.body.height)


def test_rows_hold_the_stretch_in_order_as_many_bases_each_as_asked(
    given: Given, laid_out: sequence_view.SequenceView
) -> None:
    start, end = given.span or (0, len(given.bases))
    rows = laid_out.rows
    assert rows[0].start == start
    assert rows[-1].end == end
    assert all(one.end == other.start for one, other in pairwise(rows))
    assert all(row.end - row.start == given.bases_per_row for row in rows[:-1])
    assert 0 < rows[-1].end - rows[-1].start <= given.bases_per_row
    assert all(row.bases.width == (row.end - row.start) * sequence_view.CELL for row in rows)


def test_rows_stack_down_the_view_none_reaching_into_the_next(
    laid_out: sequence_view.SequenceView,
) -> None:
    extents = [row.extent for row in laid_out.rows]
    assert all(one.y + one.height < other.y for one, other in pairwise(extents))
    canvas = laid_out.extent
    for box in extents:
        assert canvas.x <= box.x
        assert box.x + box.width <= canvas.x + canvas.width
        assert canvas.y <= box.y
        assert box.y + box.height <= canvas.y + canvas.height


def test_no_two_labels_names_or_amino_acids_overlap(laid_out: sequence_view.SequenceView) -> None:
    for row in laid_out.rows:
        boxes = [
            *(label.box for label in row.labels),
            *(_name_box(name) for name in row.names),
            *(box for _, box, _ in _residues(row)),
        ]
        assert not [(a, b) for a, b in combinations(boxes, 2) if _overlap(a, b)]


def test_every_label_lies_above_the_bases_and_every_bar_name_and_amino_acid_under_them(
    laid_out: sequence_view.SequenceView,
) -> None:
    for row in laid_out.rows:
        top, bottom = row.bases.y, row.bases.y + row.bases.height
        assert all(label.box.y + label.box.height <= top for label in row.labels)
        assert all(_band(bar).y >= bottom for bar in row.bars)
        assert all(_name_box(name).y >= bottom for name in row.names)
        assert all(box.y >= bottom for _, box, _ in _residues(row))


def test_no_amino_acid_or_name_underneath_meets_a_bar(laid_out: sequence_view.SequenceView) -> None:
    for row in laid_out.rows:
        boxes = [box for _, box, _ in _residues(row)]
        boxes += [_name_box(name) for name in row.names if not name.inside]
        assert not [(box, bar) for box in boxes for bar in row.bars if _overlap(box, _band(bar))]


def test_no_leader_crosses_a_label_other_than_its_own_and_each_runs_up_from_the_bases(
    laid_out: sequence_view.SequenceView,
) -> None:
    for row in laid_out.rows:
        for label in row.labels:
            (x0, y0), (x1, y1) = label.leader
            assert x0 == x1
            assert y0 == pytest.approx(row.bases.y)
            assert y1 == pytest.approx(label.box.y + label.box.height)
            assert label.box.x - TOUCH <= x1 <= label.box.x + label.box.width + TOUCH
            crossed = [
                other.item.label
                for other in row.labels
                if other is not label and _crosses(label.leader, other.box)
            ]
            assert not crossed


def test_each_name_lies_inside_its_bar_less_its_points_or_underneath_it_within_the_row(
    given: Given, laid_out: sequence_view.SequenceView
) -> None:
    width = given.bases_per_row * sequence_view.CELL
    for row in laid_out.rows:
        for name in row.names:
            box, body = _name_box(name), name.bar.body
            if name.inside:
                assert body.x - TOUCH <= box.x
                assert box.x + box.width <= body.x + body.width + TOUCH
                assert body.y - TOUCH <= box.y
                assert box.y + box.height <= body.y + body.height + TOUCH
                continue
            assert box.y >= body.y + body.height
            if box.width <= width:
                assert box.x >= -TOUCH
                assert box.x + box.width <= width + TOUCH


def _strand(row: sequence_view.Row, which: str) -> Box | None:
    """The box a row's top or bottom strand takes, if the row draws it."""
    for shape in row.shapes:
        if isinstance(shape, svg.Group) and shape.classes == ("strand", which):
            [letters] = shape.shapes
            assert isinstance(letters, svg.Letters)
            [(_, box)] = _boxes(letters, words=False)
            return box
    return None


def _line_box(points: Iterable[Point]) -> Box:
    xs, ys = zip(*points, strict=True)
    return Box(min(xs), min(ys), max(xs) - min(xs), max(ys) - min(ys))


def _within(inner: Box, outer: Box) -> bool:
    return (
        outer.x - TOUCH <= inner.x
        and inner.x + inner.width <= outer.x + outer.width + TOUCH
        and outer.y - TOUCH <= inner.y
        and inner.y + inner.height <= outer.y + outer.height + TOUCH
    )


def test_everything_a_row_draws_lies_within_its_extent(
    laid_out: sequence_view.SequenceView,
) -> None:
    for row in laid_out.rows:
        drawn = [
            row.bases,
            *(box for label in row.labels for box in (label.box, _line_box(label.leader))),
            *(_band(bar) for bar in row.bars),
            *(_name_box(name) for name in row.names),
            *(box for _, box, _ in _residues(row)),
            *(arrow.bounds for arrow in row.arrows),
            *(_line_box(tail.points) for tail in row.tails),
            *(mark.box for mark in row.mismatches),
            *(_line_box(line) for cut in row.cuts for line in cut.lines),
        ]
        assert [box for box in drawn if not _within(box, row.extent)] == []


def test_each_primer_lies_beside_the_strand_it_spells_apart_from_the_rest_its_marks_on_it(
    laid_out: sequence_view.SequenceView,
) -> None:
    for row in laid_out.rows:
        top = _strand(row, "top")
        assert top is not None
        lowest = _strand(row, "bottom") or top
        parts = [(arrow.item, arrow.bounds) for arrow in row.arrows]
        parts += [(tail.item, _line_box(tail.points)) for tail in row.tails]
        for item, box in parts:
            assert _within(box, row.bases)
            if item.strand == Strand.REVERSE:
                assert box.y >= lowest.y + lowest.height
            else:
                assert box.y + box.height <= top.y
        assert not [
            (one, other)
            for (item, one), (another, other) in combinations(parts, 2)
            if item is not another and _overlap(one, other)
        ]
        for mark in row.mismatches:
            assert [
                arrow
                for arrow in row.arrows
                if arrow.item is mark.item and _within(mark.box, arrow.bounds)
            ]


def test_each_cut_is_drawn_among_the_bases(laid_out: sequence_view.SequenceView) -> None:
    for row in laid_out.rows:
        for cut in row.cuts:
            assert cut.lines
            assert all(_within(_line_box(line), row.bases) for line in cut.lines)


def test_bars_lie_on_the_bases_they_draw() -> None:
    items = _items(
        _record(_feature("one", 0, 1), _feature("two", 55, 65, Strand.REVERSE), circular=False)
    )
    laid = _laid(Given(items, "ACGT" * (LENGTH // 4), circular=False))
    first, second = laid.rows[:2]
    cell = sequence_view.CELL
    assert [(bar.item.name, bar.start, bar.end) for bar in first.bars] == [
        ("one", 0, cell),
        ("two", 55 * cell, 60 * cell),
    ]
    assert [(bar.item.name, bar.start, bar.end) for bar in second.bars] == [("two", 0, 5 * cell)]


def test_a_bar_is_pointed_where_its_feature_ends_on_its_strand_and_flat_where_a_row_cuts_it() -> (
    None
):
    features = (
        _feature("forward", 50, 70),
        _feature("reverse", 50, 70, Strand.REVERSE),
        _feature("both", 50, 70, Strand.BOTH),
        _feature("none", 50, 70, Strand.NONE),
    )
    laid = _laid(Given(_items(_record(*features)), "ACGT" * (LENGTH // 4)))
    points = {
        (bar.item.name, index): (bar.point_start, bar.point_end)
        for index, row in enumerate(laid.rows[:2])
        for bar in row.bars
    }
    assert points == {
        ("forward", 0): (False, False),
        ("forward", 1): (False, True),
        ("reverse", 0): (True, False),
        ("reverse", 1): (False, False),
        ("both", 0): (True, False),
        ("both", 1): (False, True),
        ("none", 0): (False, False),
        ("none", 1): (False, False),
    }


def test_features_pack_in_tracks_by_earliest_start_in_each_row() -> None:
    features = (_feature("long", 0, 50), _feature("short", 10, 20), _feature("apart", 52, 58))
    laid = _laid(Given(_items(_record(*features)), "ACGT" * (LENGTH // 4)))
    middles = {bar.item.name: bar.middle for bar in laid.rows[0].bars}
    assert middles["long"] == middles["apart"] < middles["short"]


def test_features_running_on_into_the_next_row_keep_their_tracks_in_order() -> None:
    features = (
        _feature("third", 110, 130),
        _feature("second", 100, 160),
        _feature("first", 0, 200),
    )
    laid = _laid(Given(_items(_record(*features)), "ACGT" * (LENGTH // 4)))
    for row in laid.rows[1:3]:
        order = sorted({(bar.middle, bar.item.name) for bar in row.bars})
        assert [name for _, name in order] == ["first", "second", "third"]


def test_a_name_goes_inside_when_it_fits_underneath_when_its_track_is_clear_else_above() -> None:
    features = (
        _feature("inside", 0, 40),
        _feature("underneath", 100, 101),
        _feature("crowded out", 150, 151),
        _feature("its neighbour", 153, 154),
    )
    laid = _laid(Given(_items(_record(*features)), "ACGT" * (LENGTH // 4)))
    names = {name.bar.item.name: name.inside for row in laid.rows for name in row.names}
    assert names == {"inside": True, "underneath": False}
    labels = [label.item.name for row in laid.rows for label in row.labels]
    assert sorted(labels) == ["crowded out", "its neighbour"]


def _translation(given: Given) -> list[tuple[int, str, float, str]]:
    """Each amino acid drawn: its row, its three letters, the base it is centred on and its colour."""
    return [
        (index, word, row.start + (box.x + box.width / 2) / sequence_view.CELL - 0.5, fill)
        for index, row in enumerate(_laid(given).rows)
        for word, box, fill in _residues(row)
    ]


def test_a_translation_reads_across_joined_segments_from_its_codon_start_a_stop_in_red() -> None:
    bases = "AATGGCAAAAGTAAC"
    cds = Feature(
        "cds",
        "CDS",
        (Segment(0, 6), Segment(10, 14)),
        strand=Strand.FORWARD,
        qualifiers={"codon_start": (2,)},
    )
    given = Given(_items(_record(cds, bases=bases, circular=False)), bases, circular=False)
    # ATG, then GC across the gap to G, then TAA; each centred on its middle base.
    assert _translation(given) == [
        (0, "Met", pytest.approx(2), "#252525"),
        (0, "Ala", pytest.approx(5), "#252525"),
        (0, "Ter", pytest.approx(12), sequence_view.STOP),
    ]


def test_a_reverse_translation_reads_the_bottom_strand_across_rows_and_the_origin() -> None:
    # Read on the bottom strand from its last base, the CDS spells ATG GCG TAA across the origin.
    bases = "GCCATCCCCTTAC"
    cds = Feature("cds", "CDS", (Segment(9, 18),), strand=Strand.REVERSE)
    given = Given(_items(_record(cds, bases=bases)), bases, bases_per_row=5)
    assert _translation(given) == [
        (0, "Ala", pytest.approx(0), "#252525"),
        (0, "Met", pytest.approx(3), "#252525"),
        (2, "Ter", pytest.approx(10), sequence_view.STOP),
    ]


def test_a_primer_is_an_arrow_at_its_binding_site_beside_the_strand_it_spells() -> None:
    bases = _bases()
    primers = (
        _primer("forward", bases, 30, 52),
        _primer("reverse", bases, 70, 95, Strand.REVERSE),
    )
    items = _items(_record(bases=bases, circular=False, primers=primers))
    first, second = _laid(Given(items, bases, circular=False)).rows[:2]
    cell = sequence_view.CELL
    arrows = [
        (
            index,
            arrow.item.name,
            arrow.start / cell,
            arrow.end / cell,
            arrow.head_start,
            arrow.head_end,
        )
        for index, row in enumerate((first, second))
        for arrow in row.arrows
    ]
    # Each points to its 3' end: the forward primer right, the reverse one left.
    assert arrows == [(0, "forward", 30, 52, False, True), (1, "reverse", 10, 35, True, False)]
    top, bottom = _strand(first, "top"), _strand(second, "bottom")
    assert top is not None
    assert bottom is not None
    assert first.arrows[0].bounds.y + first.arrows[0].bounds.height <= top.y
    assert second.arrows[0].bounds.y >= bottom.y + bottom.height
    assert [[label.item.label for label in row.labels] for row in (first, second)] == [
        ["forward (31 .. 52)"],
        ["reverse (71 .. 95)"],
    ]


def _tail(row: sequence_view.Row, name: str) -> list[tuple[float, float]]:
    """A primer's tail in a row: each point along its line, as a base boundary and a height."""
    [tail] = [tail for tail in row.tails if tail.item.name == name]
    return [(point.x / sequence_view.CELL, point.y) for point in tail.points]


def test_a_tail_runs_on_from_the_arrows_back_a_base_to_a_cell_bent_away_from_the_bases() -> None:
    bases = _bases()
    primers = (
        _primer("forward", bases, 30, 52, tail="GGTCTCAAA"),
        _primer("reverse", bases, 70, 95, Strand.REVERSE, tail="GAATTC"),
        _primer("over the edge", bases, 62, 80, tail="CCCCC"),
    )
    items = _items(_record(bases=bases, circular=False, primers=primers))
    rows = _laid(Given(items, bases, circular=False)).rows[:2]
    middles = {
        (index, arrow.item.name): arrow.middle
        for index, row in enumerate(rows)
        for arrow in row.arrows
    }
    # Up from the forward arrow's back over the 9 bases before its binding site, down from the
    # reverse arrow's over the 6 after.
    forward, reverse = _tail(rows[0], "forward"), _tail(rows[1], "reverse")
    assert [x for x, _ in forward] == [30, 29, 21]
    assert forward[0][1] == middles[0, "forward"] > forward[1][1] == forward[2][1]
    assert [x for x, _ in reverse] == [35, 36, 41]
    assert reverse[0][1] == middles[1, "reverse"] < reverse[1][1] == reverse[2][1]
    # A row's edge cuts the tail flat, and it bends where it meets its arrow.
    before, after = _tail(rows[0], "over the edge"), _tail(rows[1], "over the edge")
    assert [x for x, _ in before] == [60, 57]
    assert before[0][1] == before[1][1]
    assert [x for x, _ in after] == [2, 1, 0]
    assert after[0][1] == middles[1, "over the edge"] > after[1][1]
    # A region's edge cuts it off.
    region = _laid(Given(items, bases, circular=False, span=(25, 120)))
    assert [x for x, _ in _tail(region.rows[0], "forward")] == [5, 4, 0]


def test_each_base_a_primer_does_not_pair_with_is_marked_on_its_arrow() -> None:
    bases = _bases()
    primers = (
        _primer("forward", bases, 30, 52, tail="GGTCTC", changed=(0, 5)),
        _primer("reverse", bases, 70, 95, Strand.REVERSE, tail="GAATTC", changed=(0, 3)),
        _primer("across", bases, 1190, 1210, changed=(2, 15)),
    )
    laid = _laid(Given(_items(_record(bases=bases, primers=primers)), bases))
    cell = sequence_view.CELL
    marked = set()
    for row in laid.rows:
        for mark in row.mismatches:
            base = int(mark.box.x // cell)
            assert base * cell < mark.box.x
            assert mark.box.x + mark.box.width < (base + 1) * cell
            marked.add((mark.item.name, (row.start + base) % LENGTH))
    # Counted from the 3' end: back from a forward site's end, on from a reverse site's start.
    assert marked == {
        ("forward", 51),
        ("forward", 46),
        ("reverse", 70),
        ("reverse", 73),
        ("across", 7),
        ("across", 1194),
    }


def _cuts(row: sequence_view.Row) -> dict[str, list[tuple[str, float] | tuple[str, float, float]]]:
    """Each cut site's lines in a row, by its name: through the top or the bottom strand at a base
    boundary, or along the rail from one boundary to another."""
    top, bottom = _strand(row, "top"), _strand(row, "bottom")
    assert top is not None
    cell = sequence_view.CELL
    found: dict[str, list[tuple[str, float] | tuple[str, float, float]]] = {}
    for cut in row.cuts:
        for (x1, y1), (x2, y2) in cut.lines:
            lines = found.setdefault(cut.item.name, [])
            if x1 != x2:
                assert y1 == y2
                assert top.y + top.height < y1 < (bottom.y if bottom else math.inf)
                lines.append(("rail", x1 / cell, x2 / cell))
                continue
            low, high = sorted((y1, y2))
            [strand] = [
                which
                for which, box in (("top", top), ("bottom", bottom))
                if box and low <= box.y and box.y + box.height <= high
            ]
            lines.append((strand, x1 / cell))
    return found


def test_a_cut_is_drawn_through_the_top_strand_along_the_rail_and_through_the_bottom_strand() -> (
    None
):
    cuts = layers.merge_cuts(
        [
            ("Five", 100),
            ("Three", 140),
            ("Blunt", 150),
            ("Over", 178),
            ("Edge", 240),
            ("KpnI", 270),
            ("XmaI", 270),
        ],
        LENGTH,
        staggers={"Five": 4, "Three": -4, "Over": 4, "Edge": -4, "KpnI": -4, "XmaI": 4},
    )
    bases = _bases()
    rows = _laid(Given(cuts, bases, circular=False)).rows[1:5]
    assert [_cuts(row) for row in rows] == [
        {"Five": [("top", 40), ("rail", 40, 44), ("bottom", 44)]},
        {
            "Three": [("top", 20), ("rail", 16, 20), ("bottom", 16)],
            "Blunt": [("top", 30), ("bottom", 30)],
            "Over": [("top", 58), ("rail", 58, 60)],
        },
        # The overhang runs on into the next row; one at a row's edge stays with its overhang.
        {
            "Over": [("rail", 0, 2), ("bottom", 2)],
            "Edge": [("top", 60), ("rail", 56, 60), ("bottom", 56)],
        },
        {
            "KpnI - XmaI": [
                ("top", 30),
                ("rail", 26, 30),
                ("bottom", 26),
                ("rail", 30, 34),
                ("bottom", 34),
            ]
        },
    ]
    assert [[label.item.name for label in row.labels] for row in rows] == [
        ["Five"],
        ["Three", "Blunt", "Over"],
        ["Edge"],
        ["KpnI - XmaI"],
    ]
    one = _laid(Given(cuts, bases, circular=False, both_strands=False)).rows[1:4]
    assert [_cuts(row) for row in one] == [
        {"Five": [("top", 40)]},
        {"Three": [("top", 20)], "Blunt": [("top", 30)], "Over": [("top", 58)]},
        {"Edge": [("top", 60)]},
    ]
    # Across the origin of a circular record.
    wraps = layers.merge_cuts([("Wraps", 1198)], LENGTH, staggers={"Wraps": 4})
    circle = _laid(Given(wraps, bases)).rows
    assert (_cuts(circle[-1]), _cuts(circle[0])) == (
        {"Wraps": [("top", 58), ("rail", 58, 60)]},
        {"Wraps": [("rail", 0, 2), ("bottom", 2)]},
    )


def test_enzymes_cutting_at_one_position_are_one_label_a_name_a_line_bold_where_unique() -> None:
    cuts = layers.merge_cuts(
        [("SacI", 406), ("BanII", 406), ("Ecl136II", 406), ("BanII", 900), ("EcoRI", 396)], LENGTH
    )
    row = _laid(Given(cuts, _bases())).rows[6]
    lines = {
        group.data["name"]: [
            (shape.text, shape.font, shape.y)
            for shape in group.shapes
            if isinstance(shape, svg.Text)
        ]
        for group in row.shapes
        if isinstance(group, svg.Group) and group.classes == ("cut_site", "label")
    }
    assert {name: [(text, font) for text, font, _ in line] for name, line in lines.items()} == {
        "EcoRI": [("EcoRI", BOLD)],
        "BanII - Ecl136II - SacI": [("BanII", SANS), ("Ecl136II", BOLD), ("SacI", BOLD)],
    }
    boxes = {label.item.name: label.box for label in row.labels}
    stacked = boxes["BanII - Ecl136II - SacI"]
    line = (SANS.ascender - SANS.descender) / SANS.units_per_em * sequence_view.LABEL_SIZE
    assert stacked.height == pytest.approx(boxes["EcoRI"].height + 2 * line)
    baselines = [y for _, _, y in lines["BanII - Ecl136II - SacI"]]
    assert baselines == sorted(baselines)
    assert all(stacked.y < y < stacked.y + stacked.height for y in baselines)

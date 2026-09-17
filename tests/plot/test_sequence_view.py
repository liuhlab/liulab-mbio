"""The sequence view laid out, items in and shapes out: rows holding the stretch in order, bars
and names under the bases, translations centred on their codons, and labels apart from each other
and clear of the bases, on crowded and seeded random records, whole and in regions."""

import random
from collections.abc import Callable
from dataclasses import dataclass
from itertools import combinations, pairwise

import pytest

from liulab_mbio.plot import layers, sequence_view, svg
from liulab_mbio.plot.labels import Box, Point
from liulab_mbio.sequence import BindingSite, Feature, Primer, Segment, SequenceRecord, Strand

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


def _crowded() -> Given:
    """Names crowding a few rows, over features packed in tracks, CDS joined across the origin,
    and primers and cut sites the view leaves for now."""
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
    primers = (Primer("primer", "ACGT", binding_sites=(BindingSite(300, 320, Strand.FORWARD),)),)
    cuts = layers.merge_cuts([("Cutter", 310), ("Other", 1199)], LENGTH)
    return Given((*_items(_record(*features, bases=bases, primers=primers)), *cuts), bases)


def _random(seed: int, circular: bool) -> tuple[tuple[layers.Item, ...], str]:
    """Features with random names, strands, types and joins, half clustered, on random bases."""
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
    return _items(_record(*features, bases=bases, circular=circular)), bases


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


RECORDS: dict[str, Callable[[], Given]] = {
    "crowded": _crowded,
    "crowded, across the origin": lambda: Given(
        _crowded().items, _crowded().bases, span=(1100, 1400)
    ),
    "crowded, one strand, 7 a row": lambda: Given(
        _crowded().items, _crowded().bases, bases_per_row=7, both_strands=False
    ),
    "crowded, 150 a row": lambda: Given(_crowded().items, _crowded().bases, bases_per_row=150),
    **{f"random {seed}": lambda seed=seed: _random_given(seed) for seed in range(6)},
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

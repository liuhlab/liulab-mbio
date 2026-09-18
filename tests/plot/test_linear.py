"""The linear map laid out, items in and shapes out: arrows cut cleanly where the stretch drawn cuts
them, names inside or underneath, and labels apart from each other and clear of the drawing, hiding
past the cap, on crowded and seeded random records, whole and in regions, and along longer lines."""

import random
from collections.abc import Callable
from dataclasses import dataclass, replace
from itertools import combinations, pairwise

import pytest

from liulab_mbio.plot import layers, linear, svg
from liulab_mbio.plot.labels import Box, Point
from liulab_mbio.sequence import BindingSite, Feature, Primer, Segment, SequenceRecord, Strand

from . import crowds

LENGTH = 3000
#: How far two shapes may reach into each other and still count as touching.
TOUCH = 1e-6


@dataclass(frozen=True)
class Given:
    """Items to lay out, and the stretch of the record they are drawn in."""

    items: tuple[layers.Item, ...]
    circular: bool = True
    span: tuple[int, int] | None = None
    length: int = LENGTH
    width: float = linear.WIDTH


def _feature(name: str, start: int, end: int, strand: Strand = Strand.FORWARD) -> Feature:
    return Feature(name, "misc_feature", (Segment(start, end),), strand=strand)


def _primer(name: str, start: int, end: int, strand: Strand = Strand.FORWARD) -> Primer:
    return Primer(name, "ACGT", binding_sites=(BindingSite(start, end, strand),))


def _items(
    *features: Feature, primers: tuple[Primer, ...] = (), circular: bool = True
) -> tuple[layers.Item, ...]:
    record = SequenceRecord(
        "A" * LENGTH,
        topology="circular" if circular else "linear",
        features=features,
        primers=primers,
    )
    return layers.items(record, cut_sites=False)


def _laid(given: Given) -> linear.LinearMap:
    return linear.layout(
        given.items,
        name="map",
        length=given.length,
        circular=given.circular,
        span=given.span,
        width=given.width,
    )


def _crowded() -> tuple[layers.Item, ...]:
    """Labels crowding one stretch and both ends of the line, over features packed in rows.

    Features cross the origin and each other, primers overlap one another, and enzymes cut at the
    same positions and a base apart.
    """
    features = [_feature(f"crowded name {i:02d}", 380 + 6 * i, 390 + 6 * i) for i in range(30)]
    features += [
        _feature(f"origin {i}", (2990 + 4 * i) % LENGTH, (2990 + 4 * i) % LENGTH + 6)
        for i in range(10)
    ]
    features += [_feature(f"long {i}", 100 * i, 100 * i + 1200, Strand.REVERSE) for i in range(6)]
    features += [
        _feature("across", 2800, 3300, Strand.BOTH),
        Feature(
            "joined across",
            "CDS",
            (Segment(2900, 2950), Segment(2980, 3020), Segment(60, 90)),
            strand=Strand.FORWARD,
        ),
        _feature("a name far longer than its feature", 1500, 1510),
        _feature("fits inside easily", 1600, 2200),
    ]
    primers = [_primer(f"primer {i}", 370 + 5 * i, 390 + 5 * i) for i in range(10)]
    primers += [_primer(f"across {i}", 2990 - i, 3010 - i, Strand.REVERSE) for i in range(3)]
    cuts = [(name, 385 + i) for i in range(30) for name in (f"Cutter{i}", f"Also{i}")[: i % 2 + 1]]
    cuts += [("AtTheOrigin", 0), ("BeforeIt", 2999), ("InTheMiddle", 1500), ("InTheMiddle", 400)]
    record = SequenceRecord(
        "A" * LENGTH,
        topology="circular",
        features=tuple(features),
        primers=tuple(primers),
    )
    return (*layers.items(record), *layers.merge_cuts(cuts, LENGTH))


def _random(seed: int, circular: bool) -> tuple[layers.Item, ...]:
    """Features, primers and cut sites with random names, half clustered; across the origin too on
    a circular record."""
    rng = random.Random(seed)
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789 -()"
    centres = [rng.randrange(LENGTH) for _ in range(3)]

    def span(size: int) -> tuple[int, int]:
        start = rng.randrange(LENGTH)
        if rng.random() < 0.5:
            start = (rng.choice(centres) + rng.randint(-60, 60)) % LENGTH
        if not circular:
            start = min(start, LENGTH - size)
        return start, start + size

    def name() -> str:
        return "".join(rng.choice(letters) for _ in range(rng.randint(1, 30)))

    features = []
    for i in range(rng.randint(10, 60)):
        start, end = span(rng.choice([rng.randint(1, 30), rng.randint(30, 900)]))
        segments = (Segment(start, end),)
        if end - start > 60 and rng.random() < 0.3:
            last = (end - 20) % LENGTH
            segments = (Segment(start, start + 20), Segment(last, last + 20))
        strand = rng.choice(list(Strand))
        features.append(Feature(f"{i}{name()}", "misc_feature", segments, strand=strand))
    primers = []
    for i in range(rng.randint(0, 10)):
        strand = rng.choice([Strand.FORWARD, Strand.REVERSE])
        primers.append(_primer(f"{i}{name()}", *span(rng.randint(15, 40)), strand))
    cuts = [(f"{i}{name()}", span(0)[0] or 1) for i in range(rng.randint(0, 30))]
    record = SequenceRecord(
        "A" * LENGTH,
        topology="circular" if circular else "linear",
        features=tuple(features),
        primers=tuple(primers),
    )
    return (*layers.items(record), *layers.merge_cuts(cuts, LENGTH))


def _random_region(seed: int) -> Given:
    """Random items drawn in a random stretch, which may run across the origin."""
    rng = random.Random(-seed)
    start = rng.randrange(LENGTH)
    return Given(_random(seed, circular=True), span=(start, start + rng.randint(40, LENGTH)))


#: The records of `crowds`, which `RECORDS` lays out along longer lines too.
CROWDS: dict[str, Callable[[], Given]] = {
    "pUC19 and its unique 6+ cutters": lambda: Given(crowds.puc19_crowded(), length=2686),
    "fifty EcoRI sites": lambda: Given(
        layers.items(crowds.ecori_crowd(), enzymes=["EcoRI", "HindIII"])
    ),
}

RECORDS: dict[str, Callable[[], Given]] = {
    "crowded": lambda: Given(_crowded()),
    "crowded, across the origin": lambda: Given(_crowded(), span=(2700, 3500)),
    "crowded, a region": lambda: Given(_crowded(), span=(360, 700)),
    "pUC19 and its unique 6+ cutters": CROWDS["pUC19 and its unique 6+ cutters"],
    **{f"random {seed}": lambda seed=seed: Given(_random(seed, True)) for seed in range(4)},
    **{
        f"random linear {seed}": lambda seed=seed: Given(_random(seed, False), circular=False)
        for seed in range(4, 8)
    },
    **{f"random region {seed}": lambda seed=seed: _random_region(seed) for seed in range(8, 12)},
    **{
        f"{name}, {times} times as long": lambda crowd=crowd, times=times: replace(
            crowd(), width=times * linear.WIDTH
        )
        for name, crowd in CROWDS.items()
        for times in (2, 4, 8)
    },
}


@pytest.fixture(scope="module", params=list(RECORDS))
def given(request: pytest.FixtureRequest) -> Given:
    return RECORDS[request.param]()


@pytest.fixture(scope="module")
def laid_out(given: Given) -> linear.LinearMap:
    return _laid(given)


def _overlap(one: Box, other: Box) -> bool:
    return (
        one.x < other.x + other.width - TOUCH
        and other.x < one.x + one.width - TOUCH
        and one.y < other.y + other.height - TOUCH
        and other.y < one.y + one.height - TOUCH
    )


def _crosses(leader: tuple[Point, Point], box: Box) -> bool:
    """Whether the segment passes through the box's inside, by clipping it to the box."""
    (x0, y0), (x1, y1) = leader
    low, high = 0.0, 1.0
    for p, q in (
        (-(x1 - x0), x0 - (box.x + TOUCH)),
        (x1 - x0, box.x + box.width - TOUCH - x0),
        (-(y1 - y0), y0 - (box.y + TOUCH)),
        (y1 - y0, box.y + box.height - TOUCH - y0),
    ):
        if p == 0:
            if q < 0:
                return False
        elif p < 0:
            low = max(low, q / p)
        else:
            high = min(high, q / p)
    return low < high


def _box(name: linear.Name) -> Box:
    """The box a name's letters take: from its first letter to its last, as tall as the face."""
    shape = name.letters
    font, scale = shape.font, shape.size / shape.font.units_per_em
    letters = font.letters(shape.text, shape.size)
    places = shape.places
    assert len(places) == len(letters)
    assert {place.rotate for place in places} == {0.0}
    assert len({place.y for place in places}) == 1
    left, baseline = places[0].x, places[0].y
    assert [place.x - left for place in places] == pytest.approx([one.x for one in letters])
    top = baseline - font.ascender * scale
    width = letters[-1].x + letters[-1].advance
    return Box(left, top, width, (font.ascender - font.descender) * scale)


def _band(arrow: linear.Arrow) -> Box:
    return Box(arrow.start, arrow.middle - arrow.width / 2, arrow.end - arrow.start, arrow.width)


def test_no_two_labels_or_names_overlap(laid_out: linear.LinearMap) -> None:
    boxes = [label.box for label in laid_out.labels] + [_box(name) for name in laid_out.names]
    assert not [(a, b) for a, b in combinations(boxes, 2) if _overlap(a, b)]


def test_every_label_lies_clear_of_the_line_and_the_features(laid_out: linear.LinearMap) -> None:
    assert laid_out.top < 0
    assert all(label.box.y + label.box.height <= laid_out.top for label in laid_out.labels)
    assert all(_band(arrow).y > laid_out.top for arrow in laid_out.arrows)
    assert all(_box(name).y > laid_out.top for name in laid_out.names)


def test_no_leader_crosses_a_label_other_than_its_own(laid_out: linear.LinearMap) -> None:
    crossings = [
        (label.item.label, other.item.label)
        for label in laid_out.labels
        for other in laid_out.labels
        if other is not label and _crosses(label.leader, other.box)
    ]
    assert not crossings


def test_each_leader_runs_straight_up_from_the_line_to_its_label(
    given: Given, laid_out: linear.LinearMap
) -> None:
    for label in laid_out.labels:
        (x0, y0), (x1, y1) = label.leader
        assert x0 == x1
        assert 0 <= x0 <= given.width
        assert laid_out.top <= y0 < 0
        assert y1 == pytest.approx(label.box.y + label.box.height)
        assert label.box.x - TOUCH <= x1 <= label.box.x + label.box.width + TOUCH


def test_each_name_lies_inside_its_arrow_less_its_heads_or_underneath_clear_of_every_arrow(
    laid_out: linear.LinearMap,
) -> None:
    for name in laid_out.names:
        box, body = _box(name), name.arrow.body
        if name.inside:
            assert body.x - TOUCH <= box.x
            assert box.x + box.width <= body.x + body.width + TOUCH
            assert body.y - TOUCH <= box.y
            assert box.y + box.height <= body.y + body.height + TOUCH
        else:
            assert box.y > body.y + body.height
            assert not [a.item.name for a in laid_out.arrows if _overlap(box, _band(a))]


def test_every_label_lies_within_its_rise_above_the_drawing(laid_out: linear.LinearMap) -> None:
    assert all(label.box.y >= laid_out.top - linear.RISE for label in laid_out.labels)


def test_labels_hide_cut_sites_first_then_primers_then_features(
    laid_out: linear.LinearMap,
) -> None:
    order = []
    for item in laid_out.hidden:
        cuts = min((cutter.cuts for cutter in item.cutters), default=0)
        order.append((["cut_site", "primer", "feature"].index(item.kind), -cuts, -len(item.label)))
    assert order == sorted(order)


def test_a_notice_says_what_hid_at_the_bottom_right_clear_of_everything(
    laid_out: linear.LinearMap,
) -> None:
    notices = [
        shape
        for shape in laid_out.shapes
        if isinstance(shape, svg.Group) and "notice" in shape.classes
    ]
    if not laid_out.hidden:
        assert not notices
        return
    [text] = [one for group in notices for one in group.shapes]
    assert isinstance(text, svg.Text)
    assert text.text == layers.notice(laid_out.hidden)
    scale = text.size / text.font.units_per_em
    top, right = text.y - text.font.ascender * scale, text.x + text.font.width(text.text, text.size)
    drawn = [_band(arrow) for arrow in laid_out.arrows] + [_box(name) for name in laid_out.names]
    assert top >= max(one.y + one.height for one in drawn)
    boxes = [*drawn, *(label.box for label in laid_out.labels)]
    assert right >= max(one.x + one.width for one in boxes) - TOUCH
    extent = laid_out.extent
    assert right <= extent.x + extent.width + TOUCH
    assert text.y - text.font.descender * scale <= extent.y + extent.height + TOUCH


def test_puc19_with_its_99_unique_6_cutters_hides_no_label() -> None:
    laid = _laid(Given(crowds.puc19_crowded(), length=2686))
    assert not laid.hidden
    assert len(laid.labels) == 42


def test_only_a_crowd_past_its_rise_hides_and_its_cut_sites_first() -> None:
    # Cut sites a base apart crowd a primer and two boxed names; far along the line, an enzyme that
    # cuts three times, earlier in the order than any of them, stands alone.
    cuts = [(f"Crowd{'x' * (i % 5)}{i}", 100 + i) for i in range(60)]
    cuts += [("Often", 2500), ("Often", 2600), ("Often", 2700)]
    items = (
        *_items(
            _feature("crowded out", 110, 120),
            _feature("its neighbour", 140, 150),
            primers=(_primer("in the crowd", 120, 140),),
        ),
        *layers.merge_cuts(cuts, LENGTH),
    )
    laid = _laid(Given(items))
    assert laid.hidden
    assert {item.name for item in laid.hidden} < {name for name, _ in cuts[:60]}
    lengths = [len(item.label) for item in laid.hidden]
    assert lengths == sorted(lengths, reverse=True)
    shown = {label.item.name for label in laid.labels}
    assert {"Often", "crowded out", "its neighbour", "in the crowd"} <= shown


def test_every_arrow_lies_along_the_line(given: Given, laid_out: linear.LinearMap) -> None:
    assert all(0 <= arrow.start < arrow.end <= given.width for arrow in laid_out.arrows)


def test_a_label_hidden_for_lack_of_room_shows_along_a_line_8_times_as_long() -> None:
    items = _crowded()
    hidden = _laid(Given(items)).hidden
    shown = {id(label.item) for label in _laid(Given(items, width=8 * linear.WIDTH)).labels}
    assert [item for item in hidden if id(item) in shown]


def test_a_longer_line_numbers_more_positions_as_often_as_its_points_per_base_allow() -> None:
    def numbers(length: int, times: int) -> list[str]:
        laid = linear.layout((), name="map", length=length, width=times * linear.WIDTH)
        [scale] = [
            shape
            for shape in laid.shapes
            if isinstance(shape, svg.Group) and "scale" in shape.classes
        ]
        return [shape.text for shape in scale.shapes if isinstance(shape, svg.Text)]

    counts = [len(numbers(LENGTH, times)) for times in (1, 2, 4, 8)]
    assert all(fewer < more for fewer, more in pairwise(counts))
    # As many points to a base, as often numbered.
    short = numbers(LENGTH // 8, 1)
    assert numbers(LENGTH, 8)[: len(short)] == short


def test_features_pack_in_rows_by_earliest_start() -> None:
    laid = _laid(
        Given(
            _items(
                _feature("long", 100, 900),
                _feature("short", 150, 200),
                _feature("apart", 1000, 1100),
            )
        )
    )
    middles = {arrow.item.name: arrow.middle for arrow in laid.arrows}
    assert 0 < middles["long"] == middles["apart"] < middles["short"]


def test_a_feature_running_past_either_end_of_a_region_is_cut_cleanly_there() -> None:
    items = _items(
        _feature("past both", 900, 2100, Strand.BOTH),
        _feature("past the end", 1500, 2100, Strand.REVERSE),
        _feature("past the start", 900, 1500, Strand.REVERSE),
        _feature("outside", 2500, 2600),
        circular=False,
    )
    laid = _laid(Given(items, circular=False, span=(1000, 2000)))
    arrows = {arrow.item.name: arrow for arrow in laid.arrows}
    assert sorted(arrows) == ["past both", "past the end", "past the start"]
    half = linear.WIDTH / 2
    assert {name: (arrow.start, arrow.end) for name, arrow in arrows.items()} == {
        "past both": pytest.approx((0, linear.WIDTH)),
        "past the end": pytest.approx((half, linear.WIDTH)),
        "past the start": pytest.approx((0, half)),
    }
    heads = {name: (arrow.head_start, arrow.head_end) for name, arrow in arrows.items()}
    assert heads == {
        "past both": (False, False),
        "past the end": (True, False),
        "past the start": (False, False),
    }


def test_a_feature_across_the_origin_of_a_circular_record_opened_lies_cut_at_both_ends() -> None:
    joined = Feature(
        "joined", "CDS", (Segment(2700, 2800), Segment(2900, 3100)), strand=Strand.FORWARD
    )
    laid = _laid(Given(_items(joined)))
    scale = linear.WIDTH / LENGTH
    pieces = sorted((a.start, a.end, a.head_start, a.head_end) for a in laid.arrows)
    assert pieces == [
        pytest.approx((0, 100 * scale, False, True)),
        pytest.approx((2700 * scale, 2800 * scale, False, False)),
        pytest.approx((2900 * scale, linear.WIDTH, False, False)),
    ]
    # The two segments at the end are one bar, joined across their gap, in one row.
    assert len({arrow.middle for arrow in laid.arrows if arrow.start > 0}) == 1


def test_a_region_across_the_origin_draws_what_lies_either_side_of_it() -> None:
    items = _items(
        _feature("before", 2900, 2950), _feature("after", 20, 80), _feature("no", 500, 600)
    )
    laid = _laid(Given(items, span=(2800, 3200)))
    scale = linear.WIDTH / 400
    assert {arrow.item.name: (arrow.start, arrow.end) for arrow in laid.arrows} == {
        "before": pytest.approx((100 * scale, 150 * scale)),
        "after": pytest.approx((220 * scale, 280 * scale)),
    }


def test_a_name_goes_inside_when_it_fits_underneath_when_its_row_is_clear_and_else_above() -> None:
    items = _items(
        _feature("inside", 100, 700),
        _feature("underneath", 1500, 1510),
        _feature("crowded out", 2500, 2510),
        _feature("its neighbour", 2530, 2540),
    )
    laid = _laid(Given(items))
    placed = {name.arrow.item.name: name.inside for name in laid.names}
    assert placed == {"inside": True, "underneath": False}
    assert sorted(label.item.name for label in laid.labels) == ["crowded out", "its neighbour"]
    [inside] = [name for name in laid.names if name.inside]
    assert inside.letters.fill == layers.text_color(inside.arrow.span.color)


def test_close_labels_rise_in_a_staircase_and_one_beside_them_stays_low() -> None:
    cuts = layers.merge_cuts(
        [("Beside", 700), *((f"Close{i}", 1500 + 3 * i) for i in range(5)), ("Far", 2500)], LENGTH
    )
    laid = _laid(Given(cuts))
    rows = {label.item.name: label.box.y for label in laid.labels}
    close = [rows[f"Close{i}"] for i in range(5)]
    assert close == sorted(close)
    assert len(set(close)) == 5
    assert rows["Beside"] == rows["Far"] == close[-1]


def test_a_cut_at_either_end_of_the_stretch_is_labelled_there_once() -> None:
    cuts = layers.merge_cuts([("AtTheOrigin", 0), ("AtTheStart", 1000), ("AtTheEnd", 2000)], LENGTH)
    whole = {label.item.name: label.leader[0].x for label in _laid(Given(cuts)).labels}
    assert whole == pytest.approx({"AtTheOrigin": linear.WIDTH, "AtTheStart": 240, "AtTheEnd": 480})
    region = _laid(Given(cuts, span=(1000, 2000)))
    assert {label.item.name: label.leader[0].x for label in region.labels} == pytest.approx(
        {"AtTheStart": 0, "AtTheEnd": linear.WIDTH}
    )


def test_primers_lie_above_the_line_overlapping_ones_further_up() -> None:
    primers = (_primer("one", 100, 130), _primer("over", 120, 150, Strand.REVERSE))
    laid = _laid(Given(_items(primers=primers)))
    arrows = {arrow.item.name: arrow for arrow in laid.arrows}
    assert arrows["over"].middle < arrows["one"].middle < 0
    assert (arrows["one"].head_start, arrows["one"].head_end) == (False, True)
    assert (arrows["over"].head_start, arrows["over"].head_end) == (True, False)
    assert laid.top <= arrows["over"].middle - arrows["over"].width / 2
    labels = {label.item.name: label.leader[0] for label in laid.labels}
    assert labels["one"].y < arrows["one"].middle

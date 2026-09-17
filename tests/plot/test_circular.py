"""The circular map laid out, items in and shapes out: arrows as the record has them, and labels
apart from each other and from the drawing, on crowded and seeded random records."""

import math
import random
from collections import Counter
from itertools import combinations

import pytest

from liulab_mbio.plot import circular, layers
from liulab_mbio.plot.labels import Box, Point
from liulab_mbio.sequence import BindingSite, Feature, Primer, Segment, SequenceRecord, Strand

LENGTH = 3000
#: How far two shapes may reach into each other and still count as touching.
TOUCH = 1e-6


def _layout(record: SequenceRecord) -> circular.CircularMap:
    return circular.layout(layers.items(record), name=record.name, length=len(record))


def _feature(name: str, start: int, end: int, strand: Strand = Strand.FORWARD) -> Feature:
    return Feature(name, "misc_feature", (Segment(start, end),), strand=strand)


def _primer(name: str, start: int, end: int, strand: Strand = Strand.FORWARD) -> Primer:
    return Primer(name, "ACGT", binding_sites=(BindingSite(start, end, strand),))


def _crowded() -> tuple[layers.Item, ...]:
    """Labels crowding one stretch, the origin, the bottom and the sides, over features on rings.

    Primers overlap one another, and enzymes cut at the same positions.
    """
    features = [_feature(f"crowded name {i:02d}", 380 + 6 * i, 390 + 6 * i) for i in range(40)]
    features += [
        _feature(f"origin {i}", (2990 + 4 * i) % LENGTH, (2990 + 4 * i) % LENGTH + 6)
        for i in range(12)
    ]
    features += [
        _feature(f"a much longer name at the bottom {i}", 1480 + 3 * i, 1490 + 3 * i)
        for i in range(12)
    ]
    features += [_feature(f"side {i}", 740 + 2 * i, 760 + 2 * i) for i in range(8)]
    features += [_feature(f"side {i}", 2240 + 2 * i, 2260 + 2 * i) for i in range(8)]
    features += [_feature(f"long {i}", 100 * i, 100 * i + 1200, Strand.REVERSE) for i in range(6)]
    features.append(_feature("dead on top", 2994, 3006))
    features.append(_feature("dead on the bottom", 1494, 1506))
    primers = [_primer(f"primer {i}", 370 + 5 * i, 390 + 5 * i) for i in range(10)]
    primers += [_primer(f"across {i}", 2990 - i, 3010 - i, Strand.REVERSE) for i in range(3)]
    cuts = [
        (name, 385 + 3 * i) for i in range(30) for name in (f"Cutter{i}", f"Also{i}")[: i % 2 + 1]
    ]
    cuts += [("AtTheOrigin", 0), ("BeforeIt", 2999), ("AtTheBottom", 1500), ("AtTheBottom", 400)]
    record = SequenceRecord(
        "A" * LENGTH,
        topology="circular",
        name="crowded",
        features=tuple(features),
        primers=tuple(primers),
    )
    return (*layers.items(record), *layers.merge_cuts(cuts, LENGTH))


def _random(seed: int) -> tuple[layers.Item, ...]:
    """Features, primers and cut sites with random names, half clustered, some across the origin."""
    rng = random.Random(seed)
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789 -()"
    centres = [rng.randrange(LENGTH) for _ in range(3)]

    def position() -> int:
        if rng.random() < 0.5:
            return (rng.choice(centres) + rng.randint(-60, 60)) % LENGTH
        return rng.randrange(LENGTH)

    def name() -> str:
        return "".join(rng.choice(letters) for _ in range(rng.randint(1, 30)))

    features = []
    for i in range(rng.randint(10, 100)):
        start = position()
        size = rng.choice([rng.randint(1, 30), rng.randint(30, 900)])
        features.append(_feature(f"{i}{name()}", start, start + size, rng.choice(list(Strand))))
    primers = []
    for i in range(rng.randint(0, 10)):
        start = position()
        strand = rng.choice([Strand.FORWARD, Strand.REVERSE])
        primers.append(_primer(f"{i}{name()}", start, start + rng.randint(15, 40), strand))
    cuts = [(f"{i}{name()}", position()) for i in range(rng.randint(0, 30))]
    record = SequenceRecord(
        "A" * LENGTH,
        topology="circular",
        name=f"random {seed}",
        features=tuple(features),
        primers=tuple(primers),
    )
    return (*layers.items(record), *layers.merge_cuts(cuts, LENGTH))


RECORDS = {
    "crowded": _crowded,
    **{f"random {seed}": lambda seed=seed: _random(seed) for seed in range(12)},
}


@pytest.fixture(scope="module", params=list(RECORDS))
def given(request: pytest.FixtureRequest) -> tuple[layers.Item, ...]:
    return RECORDS[request.param]()


@pytest.fixture(scope="module")
def laid_out(given: tuple[layers.Item, ...]) -> circular.CircularMap:
    return circular.layout(given, name="map", length=LENGTH)


def _overlap(one: Box, other: Box) -> bool:
    return (
        one.x < other.x + other.width - TOUCH
        and other.x < one.x + one.width - TOUCH
        and one.y < other.y + other.height - TOUCH
        and other.y < one.y + one.height - TOUCH
    )


def _nearest_to_centre(box: Box) -> float:
    dx = max(box.x, 0.0, -(box.x + box.width))
    dy = max(box.y, 0.0, -(box.y + box.height))
    return math.hypot(dx, dy)


def _closest_to_centre(leader: tuple[Point, Point]) -> float:
    (x0, y0), (x1, y1) = leader
    dx, dy = x1 - x0, y1 - y0
    t = min(1.0, max(0.0, -(x0 * dx + y0 * dy) / (dx * dx + dy * dy or 1.0)))
    return math.hypot(x0 + t * dx, y0 + t * dy)


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


def test_every_item_is_labelled_once(
    given: tuple[layers.Item, ...], laid_out: circular.CircularMap
) -> None:
    assert Counter(id(label.item) for label in laid_out.labels) == Counter(map(id, given))


def test_no_two_labels_overlap(laid_out: circular.CircularMap) -> None:
    boxes = [label.box for label in laid_out.labels]
    assert not [(a, b) for a, b in combinations(boxes, 2) if _overlap(a, b)]


def test_no_label_overlaps_the_drawing(laid_out: circular.CircularMap) -> None:
    assert all(_nearest_to_centre(label.box) >= laid_out.radius for label in laid_out.labels)


def test_no_leader_crosses_a_label_other_than_its_own(laid_out: circular.CircularMap) -> None:
    crossings = [
        (label.item.label, other.item.label)
        for label in laid_out.labels
        for other in laid_out.labels
        if other is not label and _crosses(label.leader, other.box)
    ]
    assert not crossings


def test_no_leader_enters_the_drawing(laid_out: circular.CircularMap) -> None:
    assert all(
        _closest_to_centre(label.leader) >= laid_out.radius - TOUCH for label in laid_out.labels
    )


def test_labels_and_the_drawing_lie_on_the_canvas(laid_out: circular.CircularMap) -> None:
    extent = laid_out.extent
    reach = Box(-laid_out.radius, -laid_out.radius, 2 * laid_out.radius, 2 * laid_out.radius)
    for box in [reach, *(label.box for label in laid_out.labels)]:
        assert extent.x <= box.x
        assert box.x + box.width <= extent.x + extent.width
        assert extent.y <= box.y
        assert box.y + box.height <= extent.y + extent.height


def test_each_segment_is_an_arrow_with_a_head_where_its_strand_points() -> None:
    joined = (Segment(100, 160), Segment(200, 260), Segment(300, 360))
    record = SequenceRecord(
        "A" * LENGTH,
        topology="circular",
        features=(
            Feature("forward", "CDS", joined, strand=Strand.FORWARD),
            Feature("reverse", "CDS", joined, strand=Strand.REVERSE),
            Feature("both", "CDS", joined, strand=Strand.BOTH),
            Feature("none", "CDS", joined, strand=Strand.NONE),
        ),
    )
    heads = {}
    for arrow in _layout(record).arrows:
        heads.setdefault(arrow.item.name, []).append((arrow.head_start, arrow.head_end))
    assert heads == {
        "forward": [(False, False), (False, False), (False, True)],
        "reverse": [(True, False), (False, False), (False, False)],
        "both": [(True, False), (False, False), (False, True)],
        "none": [(False, False), (False, False), (False, False)],
    }


def test_a_segment_across_the_origin_is_one_arrow() -> None:
    record = SequenceRecord(
        "A" * 2686,
        topology="circular",
        features=(Feature("across", "misc_feature", (Segment(2600, 2750),)),),
    )
    [arrow] = _layout(record).arrows
    assert arrow.start < math.tau < arrow.end
    assert arrow.end - arrow.start == pytest.approx(math.tau * 150 / 2686)


def test_features_that_overlap_go_on_rings_further_in_the_longest_outermost() -> None:
    record = SequenceRecord(
        "A" * LENGTH,
        topology="circular",
        features=(
            _feature("short", 150, 200),
            _feature("long", 100, 900),
            _feature("apart", 1000, 1100),
        ),
    )
    radii = {arrow.item.name: arrow.radius for arrow in _layout(record).arrows}
    assert radii["long"] == radii["apart"] > radii["short"]


def test_a_primer_is_an_arrow_outside_the_backbone_at_each_binding_site() -> None:
    sites = (BindingSite(2950, 3010, Strand.FORWARD), BindingSite(100, 130, Strand.REVERSE))
    record = SequenceRecord(
        "A" * LENGTH,
        topology="circular",
        primers=(Primer("both", "ACGT", binding_sites=sites), _primer("over", 110, 140)),
    )
    laid = _layout(record)
    arrows = {(arrow.item.name, arrow.span.start): arrow for arrow in laid.arrows}
    across, reverse, over = arrows["both", 2950], arrows["both", 100], arrows["over", 110]
    assert across.start < math.tau < across.end
    assert (across.head_start, across.head_end) == (False, True)
    assert (reverse.head_start, reverse.head_end) == (True, False)
    assert circular.RADIUS < across.radius == reverse.radius < over.radius
    assert over.radius + over.width / 2 < laid.radius


def test_a_cut_sites_leader_points_at_its_cut() -> None:
    items = layers.merge_cuts([("EcoRI", 750), ("BanII", 2250), ("SacI", 2250)], LENGTH)
    laid = circular.layout(items, name="cuts", length=LENGTH)
    angles = {
        label.item.label: math.atan2(label.leader[0].x, -label.leader[0].y) % math.tau
        for label in laid.labels
    }
    assert angles == {
        "EcoRI (750)": pytest.approx(math.pi / 2),
        "BanII - SacI (2250)": pytest.approx(3 * math.pi / 2),
    }

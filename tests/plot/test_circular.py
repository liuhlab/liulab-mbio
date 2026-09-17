"""The circular map laid out, items in and shapes out: arrows as the record has them, names on the
arrows they fit, and labels apart from each other and from the drawing, hiding past the cap, on
crowded and seeded random records."""

import math
import random
from collections import Counter
from dataclasses import dataclass
from itertools import combinations

import pytest

from liulab_mbio.enzymes import enzymes
from liulab_mbio.plot import circular, layers, svg
from liulab_mbio.plot.fonts import SANS
from liulab_mbio.plot.labels import Box, Point
from liulab_mbio.sequence import BindingSite, Feature, Primer, Segment, SequenceRecord, Strand
from liulab_mbio.sites import find_sites

from . import crowds

LENGTH = 3000
#: How far two shapes may reach into each other and still count as touching.
TOUCH = 1e-6


@dataclass(frozen=True)
class Given:
    """Items to lay out, and the length of the record they lie on."""

    items: tuple[layers.Item, ...]
    length: int = LENGTH


def _layout(record: SequenceRecord) -> circular.CircularMap:
    return circular.layout(layers.items(record), name=record.name, length=len(record))


def _feature(name: str, start: int, end: int, strand: Strand = Strand.FORWARD) -> Feature:
    return Feature(name, "misc_feature", (Segment(start, end),), strand=strand)


def _primer(name: str, start: int, end: int, strand: Strand = Strand.FORWARD) -> Primer:
    return Primer(name, "ACGT", binding_sites=(BindingSite(start, end, strand),))


def _one(name: str, start: int, end: int, strand: Strand = Strand.FORWARD) -> circular.CircularMap:
    """One feature laid out, which may cross the origin."""
    feature = _feature(name, start % LENGTH, start % LENGTH + end - start, strand)
    record = SequenceRecord("A" * LENGTH, topology="circular", features=(feature,))
    return circular.layout(layers.items(record, cut_sites=False), name="one", length=LENGTH)


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


def _near_fit() -> tuple[layers.Item, ...]:
    """Names all round the circle on arrows about as long as they are, five rings deep."""
    rng = random.Random(1)
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789 -()αΔ"
    features = []
    for i in range(12):
        for ring in range(5):
            size = 240 - 10 * ring - rng.randint(0, 9)
            name = "".join(rng.choice(letters) for _ in range(rng.randint(1, 10)))
            start = (LENGTH * i // 12 - size // 2) % LENGTH
            features.append(_feature(name, start, start + size, rng.choice(list(Strand))))
    record = SequenceRecord("A" * LENGTH, topology="circular", features=tuple(features))
    return layers.items(record)


RECORDS = {
    "crowded": lambda: Given(_crowded()),
    "near fit": lambda: Given(_near_fit()),
    "pUC19 and its unique 6+ cutters": lambda: Given(crowds.puc19_crowded(), 2686),
    **{f"random {seed}": lambda seed=seed: Given(_random(seed)) for seed in range(12)},
}


@pytest.fixture(scope="module", params=list(RECORDS))
def given(request: pytest.FixtureRequest) -> Given:
    return RECORDS[request.param]()


@pytest.fixture(scope="module")
def laid_out(given: Given) -> circular.CircularMap:
    return circular.layout(given.items, name="map", length=given.length)


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


def _letter_boxes(name: circular.Name) -> list[list[Point]]:
    """Each letter's box as drawn, corner by corner: its advance along its turned baseline, and the
    face's height across it."""
    shape = name.letters
    font, scale = shape.font, shape.size / shape.font.units_per_em
    boxes = []
    for letter, (x, y, rotate) in zip(
        font.letters(shape.text, shape.size), shape.places, strict=True
    ):
        turn = math.radians(rotate)
        along, up = (math.cos(turn), math.sin(turn)), (math.sin(turn), -math.cos(turn))
        boxes.append(
            [
                Point(x + a * along[0] + h * up[0], y + a * along[1] + h * up[1])
                for a, h in (
                    (0.0, font.descender * scale),
                    (letter.advance, font.descender * scale),
                    (letter.advance, font.ascender * scale),
                    (0.0, font.ascender * scale),
                )
            ]
        )
    return boxes


def test_every_item_is_named_on_its_arrow_labelled_or_hidden_once(
    given: Given, laid_out: circular.CircularMap
) -> None:
    named = [id(name.arrow.item) for name in laid_out.names]
    labelled = [id(label.item) for label in laid_out.labels]
    hidden = [id(item) for item in laid_out.hidden]
    assert Counter(named + labelled + hidden) == Counter(map(id, given.items))


def test_every_label_lies_within_reach_of_the_centre(laid_out: circular.CircularMap) -> None:
    assert all(
        label.box.y >= -circular.REACH and label.box.y + label.box.height <= circular.REACH
        for label in laid_out.labels
    )


def _hides_before(item: layers.Item) -> tuple[int, int, int]:
    """Cut sites first, cutting most often first; then primers; then features; longest first."""
    cuts = min((cutter.cuts for cutter in item.cutters), default=0)
    return ["cut_site", "primer", "feature"].index(item.kind), -cuts, -len(item.label)


def test_labels_hide_cut_sites_first_then_primers_then_features(
    laid_out: circular.CircularMap,
) -> None:
    order = [_hides_before(item) for item in laid_out.hidden]
    assert order == sorted(order)


def _notice(shapes: tuple[svg.Shape, ...]) -> list[tuple[str, Box]]:
    """Each notice drawn, as its words and the box its text takes."""
    found = []
    for shape in shapes:
        if isinstance(shape, svg.Group) and "notice" in shape.classes:
            [text] = shape.shapes
            assert isinstance(text, svg.Text)
            scale = text.size / text.font.units_per_em
            top = text.y - text.font.ascender * scale
            height = (text.font.ascender - text.font.descender) * scale
            found.append(
                (text.text, Box(text.x, top, text.font.width(text.text, text.size), height))
            )
    return found


def test_a_notice_says_what_hid_at_the_bottom_right_clear_of_everything(
    laid_out: circular.CircularMap,
) -> None:
    notices = _notice(laid_out.shapes)
    if not laid_out.hidden:
        assert not notices
        return
    [(words, box)] = notices
    assert words == layers.notice(laid_out.hidden)
    boxes = [label.box for label in laid_out.labels]
    assert box.y >= max([laid_out.radius, *(one.y + one.height for one in boxes)])
    assert (
        box.x + box.width >= max([laid_out.radius, *(one.x + one.width for one in boxes)]) - TOUCH
    )
    extent = laid_out.extent
    assert extent.x <= box.x
    assert box.x + box.width <= extent.x + extent.width + TOUCH
    assert box.y + box.height <= extent.y + extent.height + TOUCH


def test_puc19_with_its_99_unique_6_cutters_hides_no_label() -> None:
    items = crowds.puc19_crowded()
    laid = circular.layout(items, name="pUC19", length=2686)
    assert sum(item.kind == "cut_site" for item in items) == 38
    assert not laid.hidden
    assert len(laid.labels) == 44


def test_the_pinned_cutters_cut_where_the_shipped_enzymes_do() -> None:
    pinned = dict(crowds.unique_6_cutters())
    shipped = {site.enzyme.name: site.top_cut for site in find_sites(crowds.puc19(), enzymes())}
    assert len(pinned) == 99
    both = pinned.keys() & shipped.keys()
    assert len(both) == 16
    assert {name: pinned[name] for name in both} == {name: shipped[name] for name in both}


def test_only_a_crowd_past_the_reach_hides_and_its_cut_sites_first() -> None:
    # Near the top, cut sites crowd a primer and two boxed names; lower down the same column, an
    # enzyme that cuts three times, earlier in the order than any of them, stands alone.
    cuts = [(f"Crowd{'x' * (i % 5)}{i}", 60 + 2 * i) for i in range(30)]
    cuts += [("Often", 1250), ("Often", 1260), ("Often", 1270)]
    record = SequenceRecord(
        "A" * LENGTH,
        topology="circular",
        features=(_feature("tiny", 70, 72), _feature("small", 100, 102)),
        primers=(_primer("in the crowd", 80, 100),),
    )
    items = (*layers.items(record), *layers.merge_cuts(cuts, LENGTH))
    laid = circular.layout(items, name="crowd", length=LENGTH)
    assert laid.hidden
    assert {item.name for item in laid.hidden} < {name for name, _ in cuts[:30]}
    lengths = [len(item.label) for item in laid.hidden]
    assert lengths == sorted(lengths, reverse=True)
    shown = {label.item.name for label in laid.labels}
    assert {"Often", "tiny", "small", "in the crowd"} <= shown


def test_each_name_on_an_arrow_lies_inside_its_arrow_less_its_heads(
    laid_out: circular.CircularMap,
) -> None:
    outside = []
    for name in laid_out.names:
        arrow = name.arrow
        first = arrow.start + arrow.head * arrow.head_start
        last = arrow.end - arrow.head * arrow.head_end
        for corners in _letter_boxes(name):
            edges = list(zip(corners, corners[1:] + corners[:1], strict=True))
            if (
                max(math.hypot(x, y) for x, y in corners) > arrow.radius + arrow.width / 2 + TOUCH
                or min(map(_closest_to_centre, edges)) < arrow.radius - arrow.width / 2 - TOUCH
                or any(
                    (math.atan2(x, -y) - first + TOUCH) % math.tau > last - first + 2 * TOUCH
                    for x, y in corners
                )
            ):
                outside.append(name.letters.text)
    assert not outside


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


def test_a_name_fits_on_an_arrow_on_any_ring_and_one_that_does_not_is_boxed() -> None:
    laid_out = circular.layout(_near_fit(), name="near fit", length=LENGTH)
    assert {name.arrow.radius for name in laid_out.names} == {a.radius for a in laid_out.arrows}
    assert laid_out.labels


def test_a_name_goes_on_its_arrow_once_the_band_less_its_head_holds_it() -> None:
    width = SANS.width("ori", circular.LABEL_SIZE)
    on_arrow = []
    for size in range(40, 140):
        [arrow] = (laid_out := _one("ori", 1000, 1000 + size)).arrows
        room = arrow.radius * (arrow.end - arrow.start - arrow.head)
        # No size threshold: a name is boxed only while its band is short of the name and an em.
        if laid_out.names:
            assert width <= room
        else:
            assert room < width + circular.LABEL_SIZE
        on_arrow.append(bool(laid_out.names))
    assert on_arrow == sorted(on_arrow)
    assert any(on_arrow)


@pytest.mark.parametrize("degrees", [0, 60, 120, 180, 240, 300])
def test_a_name_on_an_arrow_curves_along_its_band_and_reads_upright(degrees: int) -> None:
    middle = LENGTH * degrees // 360
    [name] = _one("ΔlacZ AVAW", middle - 150, middle + 150, Strand.NONE).names
    shape, arrow = name.letters, name.arrow
    letters = SANS.letters(shape.text, shape.size)
    rise = (SANS.ascender + SANS.descender) / 2 / SANS.units_per_em * shape.size
    middles, alongs = [], []
    for letter, (x, y, rotate) in zip(letters, shape.places, strict=True):
        # Upright: each letter's top points up the page.
        assert abs(rotate) < 90
        turn = math.radians(rotate)
        along, up = (math.cos(turn), math.sin(turn)), (math.sin(turn), -math.cos(turn))
        middles.append(
            Point(
                x + letter.advance / 2 * along[0] + rise * up[0],
                y + letter.advance / 2 * along[1] + rise * up[1],
            )
        )
        alongs.append(along)
    # Each letter's middle lies on the band's middle, turned square to it.
    for (x, y), along in zip(middles, alongs, strict=True):
        assert math.hypot(x, y) == pytest.approx(arrow.radius)
        assert x * along[0] + y * along[1] == pytest.approx(0, abs=1e-6)
    # Letters follow on in reading order, as far apart round the band as the font tables set them.
    for i in range(1, len(letters)):
        (x0, y0), (x1, y1) = middles[i - 1], middles[i]
        assert (x1 - x0) * alongs[i - 1][0] + (y1 - y0) * alongs[i - 1][1] > 0
        apart = arrow.radius * math.acos(max(-1.0, min(1.0, (x0 * x1 + y0 * y1) / arrow.radius**2)))
        spaced = (
            letters[i].x + letters[i].advance / 2 - letters[i - 1].x - letters[i - 1].advance / 2
        )
        assert apart == pytest.approx(spaced)

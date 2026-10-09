import pytest

from mbio.enzymes import get_enzyme
from mbio.reaction import Pool, Reaction
from mbio.sequence import SequenceRecord


@pytest.fixture
def one() -> SequenceRecord:
    return SequenceRecord("AAAAGGTCTCGTTTTCCCC", name="one")


@pytest.fixture
def other() -> SequenceRecord:
    return SequenceRecord("CCCCGGTCTCGAAAATTTT", name="other")


def test_a_pool_holds_one_species_or_many(one: SequenceRecord, other: SequenceRecord) -> None:
    single = Pool("destination", [one])
    assert single.one is one
    assert len(Pool("donor", [one, other])) == 2
    assert list(Pool("donor", [one, other])) == [one, other]


def test_a_pool_of_several_refuses_to_answer_as_one(
    one: SequenceRecord, other: SequenceRecord
) -> None:
    pool = Pool("donor", [one, other])
    with pytest.raises(ValueError, match="holds 2 species"):
        _ = pool.one


def test_a_role_nothing_fills_is_refused() -> None:
    with pytest.raises(ValueError, match="holds no molecule"):
        Pool("destination", [])


def test_a_reaction_binds_each_role_to_its_pool(one: SequenceRecord, other: SequenceRecord) -> None:
    made = Reaction(
        "one pot",
        "one-pot",
        pools=(Pool("destination", [one]), Pool("donor", [other])),
        enzymes=("BsaI",),
    )
    assert made.roles == ("destination", "donor")
    assert made["donor"].one is other
    assert made.records == (one, other)
    assert "insert" not in made
    assert made.acting == (get_enzyme("BsaI"),)


def test_a_reaction_reads_an_enzyme_by_name_or_by_record(one: SequenceRecord) -> None:
    named = Reaction("cut", "digest", pools=(Pool("donor", [one]),), enzymes=("BsaI",))
    held = Reaction("cut", "digest", pools=(Pool("donor", [one]),), enzymes=(get_enzyme("BsaI"),))
    assert named.acting == held.acting


def test_a_ligation_a_ligase_alone_drives_names_no_enzyme(one: SequenceRecord) -> None:
    assert Reaction("join", "ligation", pools=(Pool("insert", [one]),)).acting == ()


def test_an_array_of_vessels_is_one_reaction(one: SequenceRecord) -> None:
    seating = Reaction("seating", "one-pot", pools=(Pool("carrier", [one]),), vessels=72)
    assert seating.vessels == 72


@pytest.mark.parametrize(
    ("changed", "says"),
    [
        ({"pools": ()}, "holds no molecule"),
        ({"vessels": 0}, "at least one tube"),
        ({"enzymes": ("NotAnEnzyme",)}, "NotAnEnzyme"),
    ],
)
def test_a_tube_nothing_could_be_is_refused(
    one: SequenceRecord, changed: dict[str, object], says: str
) -> None:
    given: dict[str, object] = {"pools": (Pool("donor", [one]),), **changed}
    with pytest.raises((ValueError, KeyError), match=says):
        Reaction("bad", "digest", **given)  # type: ignore[arg-type]


def test_a_role_bound_twice_is_refused(one: SequenceRecord, other: SequenceRecord) -> None:
    with pytest.raises(ValueError, match="binds a role twice"):
        Reaction(
            "bad", "digest", pools=(Pool("donor", [one]), Pool("donor", [other])), enzymes=("BsaI",)
        )


def test_a_role_this_reaction_does_not_bind_is_a_key_error(one: SequenceRecord) -> None:
    with pytest.raises(KeyError):
        Reaction("cut", "digest", pools=(Pool("donor", [one]),))["destination"]


def test_a_pool_names_the_enzyme_that_made_its_ends():
    """A ligation joins ends it did not cut, so the pool says which enzyme left them."""
    one = SequenceRecord("ACGT")

    assert Pool("donor", [one]).cutter is None
    assert Pool("donor", [one], enzyme="BsaI").cutter == get_enzyme("BsaI")
    assert Pool("donor", [one], enzyme=get_enzyme("BsaI")).cutter == get_enzyme("BsaI")
    with pytest.raises(KeyError):
        Pool("donor", [one], enzyme="NotAnEnzyme")

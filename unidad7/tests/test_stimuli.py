"""!
@file test_stimuli.py
@brief Unit tests for stimuli.py.
"""

from collections import Counter

import pytest

from rs_code import RSCode
from stimuli import Case, StimulusGenerator


@pytest.fixture(scope="module", params=["qr", "small"])
def code(request):
    return RSCode.qr_2m() if request.param == "qr" else RSCode.tp3_reduced()


@pytest.fixture
def gen(code):
    return StimulusGenerator(code, seed=31)


def check_case(code, case):
    """!
    @brief Every field of a case must agree with the others.
    """
    assert isinstance(case, Case)
    assert code.encode(case.message) == case.codeword
    assert code.add_errors(case.codeword, case.errors) == case.received
    assert code.syndromes(case.received) == case.syndromes
    assert all(v != 0 for v in case.errors.values())
    assert all(0 <= p < code.n for p in case.errors)
    assert case.n_errors == len(case.errors)
    assert case.correctable == (case.n_errors <= code.t)


class TestMessages:
    """!
    @brief Message patterns.
    """

    def test_patterns(self, code, gen):
        top = code.field.order - 1
        assert gen.message("zeros") == [0] * code.k
        assert gen.message("ones") == [top] * code.k
        assert gen.message("ramp") == [i % code.field.order for i in range(code.k)]
        msg = gen.message()
        assert len(msg) == code.k and all(0 <= s <= top for s in msg)

    def test_unknown_pattern(self, gen):
        with pytest.raises(ValueError):
            gen.message("checkerboard")


class TestCases:
    """!
    @brief Each generator returns consistent cases of the requested kind.
    """

    def test_clean(self, code, gen):
        for pattern in ("random", "zeros", "ones", "ramp"):
            case = gen.clean(pattern)
            check_case(code, case)
            assert case.kind == f"clean_{pattern}"
            assert case.n_errors == 0 and not any(case.syndromes)

    def test_explicit_case(self, code, gen):
        case = gen.case({0: 1, 3: 0}, message=[1] * code.k)
        check_case(code, case)
        assert case.errors == {0: 1} and case.kind == "custom"
        with pytest.raises(ValueError):
            gen.case({code.n: 1})

    def test_random_errors(self, code, gen):
        for count in range(code.n + 1):
            case = gen.random_errors(count)
            check_case(code, case)
            assert case.n_errors == count and case.kind == f"errors_{count}"

    def test_restricted_positions_and_values(self, code, gen):
        case = gen.random_errors(2, positions=[4, 5, 6], values=[3], kind="mine")
        check_case(code, case)
        assert set(case.errors) <= {4, 5, 6} and set(case.errors.values()) == {3}
        assert case.kind == "mine"

    def test_data_and_parity(self, code, gen):
        for count in range(1, code.t + 1):
            data = gen.data_errors(count)
            parity = gen.parity_errors(count)
            check_case(code, data)
            check_case(code, parity)
            assert all(p >= code.nsym for p in data.errors)
            assert all(p < code.nsym for p in parity.errors)

    def test_extreme_values(self, code, gen):
        case = gen.extreme_values(code.t)
        check_case(code, case)
        assert set(case.errors.values()) <= {1, code.field.order - 1}

    def test_edges(self, code, gen):
        case = gen.edge_errors()
        check_case(code, case)
        assert set(case.errors) == {0, code.n - 1}

    def test_burst(self, code, gen):
        case = gen.burst(code.t, start=2)
        check_case(code, case)
        assert sorted(case.errors) == list(range(2, 2 + code.t))
        for _ in range(20):
            case = gen.burst(3)
            positions = sorted(case.errors)
            assert positions == list(range(positions[0], positions[0] + 3))
        assert gen.burst(code.n).n_errors == code.n

    def test_random_syndromes(self, code, gen):
        syn = gen.random_syndromes()
        assert len(syn) == 2 * code.t
        assert all(0 <= s < code.field.order for s in syn)

    @pytest.mark.parametrize("call", [
        lambda g: g.random_errors(-1),
        lambda g: g.random_errors(g.code.n + 1),
        lambda g: g.random_errors(3, positions=[0, 1]),
        lambda g: g.random_errors(1, values=[0, 1]),
        lambda g: g.random_errors(1, values=[]),
        lambda g: g.burst(0),
        lambda g: g.burst(g.code.n + 1),
        lambda g: g.burst(3, start=-1),
        lambda g: g.burst(3, start=g.code.n - 2),
    ])
    def test_rejects_bad_arguments(self, gen, call):
        with pytest.raises(ValueError):
            call(gen)


class TestCollections:
    """!
    @brief suite() and stream().
    """

    def test_suite_covers_every_category(self, code, gen):
        cases = gen.suite(per_count=2)
        for case in cases:
            check_case(code, case)
        counts = Counter(c.n_errors for c in cases)
        assert all(counts[e] > 0 for e in range(code.t + 1))
        assert any(not c.correctable for c in cases)
        assert counts[code.n] > 0
        kinds = {c.kind for c in cases}
        assert {"clean_zeros", "clean_ones", "clean_ramp", "clean_random", "edges",
                f"burst_{code.t}", f"beyond_{code.t + 1}"} <= kinds

    def test_suite_is_reproducible(self, code):
        a = StimulusGenerator(code, seed=5).suite()
        b = StimulusGenerator(code, seed=5).suite()
        c = StimulusGenerator(code, seed=6).suite()
        assert a == b
        assert a != c

    def test_stream_default_weights(self, code, gen):
        cases = list(gen.stream(300))
        assert len(cases) == 300
        counts = Counter(c.n_errors for c in cases)
        assert set(counts) <= set(range(code.t + 2)) | {2 * code.t}
        assert counts[0] > counts[1]

    def test_stream_custom_weights(self, code, gen):
        cases = list(gen.stream(50, {2: 1, 5: 0}))
        assert all(c.n_errors == 2 for c in cases)

    @pytest.mark.parametrize("weights", [{-1: 1}, {100: 1}, {1: 0}])
    def test_stream_rejects_bad_weights(self, gen, weights):
        with pytest.raises(ValueError):
            list(gen.stream(5, weights))

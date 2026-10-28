"""!
@file test_sc_bind.py
@brief SystemC called from Python: one elaboration, several sc_start() calls.
"""

import sc_bind


def test_systemc_runs_in_steps():
    first = sc_bind.run(10)
    assert first > 0
    assert sc_bind.run(5) == first + 5
    assert sc_bind.run(0) == first + 5
    assert sc_bind.now_ns() == 150


def test_release():
    assert sc_bind.release().startswith("3.0")

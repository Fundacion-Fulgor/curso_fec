"""!
@file test_run_demo.py
@brief The RTL must keep reproducing Table 4.3 of Zhang, so the demo does not break.
"""

import shutil

import pytest

from run_demo import run


def test_demo_matches_the_book(tmp_path):
    if shutil.which("verilator") is None:
        pytest.skip("verilator was not found")
    assert run(tmp_path)

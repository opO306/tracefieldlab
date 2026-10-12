import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ff_prescan import prescan, refine


def test_prescan_n1_values():
    r = prescan([1, 0.8, 0.8], 0.0009114239, 0.3019805958, 1.8286248401)
    assert 0.02 < r["min_dist"] < 0.05 and r["vmax_half"] < 1e-6 and r["x2_minus_x3"] < -0.3


def test_refine_reduces_residual():
    r = refine([1, 0.8, 0.8], 0.0009114239, 0.3019805958, 1.8286248401 / 2)
    assert r["resid"] < 1e-11
    assert abs(r["x"] - 0.0009114239) < 1e-8 and abs(r["tau"] - 0.91431242) < 1e-6

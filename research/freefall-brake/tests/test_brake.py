import os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ff_common import dyad, DI
from ff_brake import rest_velocity_check, brake_problem

dyad.set_prec(128)


def _P(x2, x3):
    pos = [DI(F(-1, 2)), DI(0), DI(x2), DI(0), DI(x3), DI(F(1, 3))]
    return pos + [DI(0)] * 6


def test_nonsingular_conservation_matrix_ok():
    ok, detail = rest_velocity_check([1, 1, 1], (6, 7, 8), _P(F(1, 2), F(1, 5)))
    assert ok, detail


def test_singular_when_x2_equals_x3():
    ok, _ = rest_velocity_check([1, 1, 1], (6, 7, 8), _P(F(1, 5), F(1, 5)))
    assert not ok


def test_box_straddling_singularity_fails():
    P = _P(F(1, 2), F(1, 5))
    P[2] = DI(F(1, 10), F(3, 10))            # x2 box contains x3 = 1/5
    ok, _ = rest_velocity_check([1, 1, 1], (6, 7, 8), P)
    assert not ok


def test_brake_problem_shape():
    p = brake_problem("t", [1, F(4, 5), F(4, 5)])
    assert p.param_comps == [4, 5] and p.eq_comps == [6, 7, 8] and p.tol_control
    assert p.base[:4] == [F(-1, 2), 0, F(1, 2), 0] and p.base[6:] == [0] * 6

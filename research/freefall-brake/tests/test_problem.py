import math, os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pytest
from mpmath import mp
from ff_common import dyad, ts
from ff_problem import Problem, run_problem

dyad.set_prec(128)
mp.dps = 60
R = F(1, 2 ** 40)


def contains_real(pair, x, eps=F(1, 10 ** 50)):
    lo, hi = F(pair[0]), F(pair[1])
    q = F(mp.nstr(x, 58))
    return lo <= q - eps and q + eps <= hi


def harmonic():
    return Problem("C1a_harmonic", ts.make_harmonic(), [F(1), F(0)], [], [1], [F(0)], F(1, 8))


def kepler(base_x=F(3, 2)):
    return Problem("C1b_kepler", ts.make_kepler(1), [base_x, F(0), F(0), F(0)],
                   [3], [0, 1], [F(-1, 2), F(0)], F(1, 64))


def test_c1a_harmonic_pass_contains_pi():
    res = run_problem(harmonic(), [F(math.pi)], [R])
    assert res["verdict"] == "PASS", res["reason"]
    assert contains_real(res["K"][0], mp.pi)


def test_c1a_wrong_centre_refused():
    res = run_problem(harmonic(), [F(3)], [R])
    assert res["verdict"] == "REFUSED"


def test_c1b_kepler_pass_contains_exact():
    res = run_problem(kepler(), [F(1 / math.sqrt(3)), F(math.pi)], [R, R])
    assert res["verdict"] == "PASS", res["reason"]
    assert contains_real(res["K"][0], 1 / mp.sqrt(3))
    assert contains_real(res["K"][1], mp.pi)


def test_collision_start_is_refused_not_crash():
    res = run_problem(kepler(base_x=F(0)), [F(1, 2), F(1)], [R, R])
    assert res["verdict"] == "REFUSED" and "singular" in res["reason"]


def test_zero_radius_rejected():
    with pytest.raises(ValueError):
        run_problem(harmonic(), [F(math.pi)], [F(0)])


def test_non_dyadic_base_is_handled():
    p = Problem("h13", ts.make_harmonic(), [F(1, 3), F(0)], [], [1], [F(0)], F(1, 8))
    res = run_problem(p, [F(math.pi)], [R])
    assert res["verdict"] == "PASS", res["reason"]
    assert contains_real(res["K"][0], mp.pi)


def test_tau_box_at_zero_is_refused_not_crash():
    res = run_problem(harmonic(), [F(0)], [R])
    assert res["verdict"] == "REFUSED"
